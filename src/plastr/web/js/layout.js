/**
 * layout.js — layout engine + the main-thread controller that drives it.
 *
 * Two exports:
 *
 *   LayoutEngine     pure maths, no DOM and no worker APIs, so it can run
 *                    either inside `layout-worker.js` (the normal case) or on
 *                    the main thread as a fallback when workers are blocked.
 *
 *   LayoutController owns the module Worker, converts UI intent ("rearrange the
 *                    selection, circular") into engine commands, and copies
 *                    positions back into the GraphModel as ticks arrive.
 *
 * Force model
 * -----------
 * Particles feel:
 *   * Barnes-Hut repulsion (O(n log n)) from every other particle,
 *   * spring attraction along links, between the two *correct* end particles,
 *   * a soft bending spring between particles i and i+2 of the same segment,
 *     which stops polylines curling into knots,
 *   * gravity towards the centroid, keeping components from drifting away,
 * and after integration the consecutive-particle distances inside each segment
 * are restored by a few position-based constraint passes, which is what makes
 * the intra-segment springs behave as *stiff* springs at any time step.
 */

import { END_START, END_END } from './graph.js';

const MAX_TREE_DEPTH = 30;

/* ====================================================================== */
/*  Barnes-Hut quadtree over a flat particle array                        */
/* ====================================================================== */

class QuadTree {
  constructor() {
    this.cap = 0;
    this.n = 0;
    this._grow(4096);
    this.stack = new Int32Array(8192);
  }

  _grow(cap) {
    const old = this.cap;
    this.cap = cap;
    const cx = new Float32Array(cap); const cy = new Float32Array(cap);
    const mass = new Float32Array(cap); const size = new Float32Array(cap);
    const bx = new Float32Array(cap); const by = new Float32Array(cap);
    const body = new Int32Array(cap); const child = new Int32Array(cap * 4);
    if (old) {
      cx.set(this.cx); cy.set(this.cy); mass.set(this.mass); size.set(this.size);
      bx.set(this.bx); by.set(this.by); body.set(this.body); child.set(this.child);
    }
    this.cx = cx; this.cy = cy; this.mass = mass; this.size = size;
    this.bx = bx; this.by = by; this.body = body; this.child = child;
  }

  _alloc(bx, by, size) {
    if (this.n >= this.cap) this._grow(this.cap * 2);
    const i = this.n++;
    this.mass[i] = 0; this.cx[i] = 0; this.cy[i] = 0; this.body[i] = -1;
    this.bx[i] = bx; this.by[i] = by; this.size[i] = size;
    const c = i * 4;
    this.child[c] = -1; this.child[c + 1] = -1; this.child[c + 2] = -1; this.child[c + 3] = -1;
    return i;
  }

  _childIndex(node, x, y) {
    const h = this.size[node] * 0.5;
    const mx = this.bx[node] + h;
    const my = this.by[node] + h;
    const q = (x >= mx ? 1 : 0) + (y >= my ? 2 : 0);
    let c = this.child[node * 4 + q];
    if (c < 0) {
      c = this._alloc((q & 1) ? mx : this.bx[node], (q & 2) ? my : this.by[node], h);
      this.child[node * 4 + q] = c; // re-read: _alloc may have reallocated
    }
    return c;
  }

  /** Build the tree over particles [0, count). */
  build(px, py, count) {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (let i = 0; i < count; i++) {
      const x = px[i], y = py[i];
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    if (!Number.isFinite(minX)) { minX = minY = -1; maxX = maxY = 1; }
    const span = Math.max(maxX - minX, maxY - minY, 1) * 1.02 + 2;
    this.n = 0;
    this._alloc(minX - 1, minY - 1, span);
    for (let i = 0; i < count; i++) this._insert(i, px[i], py[i], px, py);
  }

  _insert(i, x, y, px, py) {
    let node = 0, depth = 0;
    for (;;) {
      const m = this.mass[node];
      if (m === 0) { // empty node: becomes a leaf holding this body
        this.mass[node] = 1; this.cx[node] = x; this.cy[node] = y; this.body[node] = i;
        return;
      }
      if (this.body[node] >= 0) { // occupied leaf: push the sitting body down
        const j = this.body[node];
        this.body[node] = -1;
        if (depth < MAX_TREE_DEPTH) {
          const c = this._childIndex(node, px[j], py[j]);
          this.mass[c] = 1; this.cx[c] = px[j]; this.cy[c] = py[j]; this.body[c] = j;
        }
      }
      // accumulate into this (now internal) node
      const inv = 1 / (m + 1);
      this.cx[node] = (this.cx[node] * m + x) * inv;
      this.cy[node] = (this.cy[node] * m + y) * inv;
      this.mass[node] = m + 1;
      if (depth >= MAX_TREE_DEPTH) return; // coincident points: stop here
      node = this._childIndex(node, x, y);
      depth++;
    }
  }

  /**
   * Repulsive force on particle `i` at (x, y). Writes into out[0], out[1].
   * `rep` is the Fruchterman-Reingold constant k^2; magnitude is k^2/d.
   */
  force(i, x, y, theta2, rep, out) {
    const stack = this.stack;
    let sp = 0;
    stack[sp++] = 0;
    let fx = 0, fy = 0;
    const limit = stack.length - 5;
    while (sp > 0) {
      const node = stack[--sp];
      const m = this.mass[node];
      if (m === 0) continue;
      const b = this.body[node];
      if (b === i) continue;
      let dx = x - this.cx[node];
      let dy = y - this.cy[node];
      let d2 = dx * dx + dy * dy;
      const s = this.size[node];
      if (b >= 0 || s * s < theta2 * d2 || sp >= limit) {
        if (d2 < 0.25) { // coincident: nudge apart deterministically enough
          dx = (i % 7) - 3 + 0.5;
          dy = (i % 5) - 2 + 0.5;
          d2 = dx * dx + dy * dy + 0.5;
        }
        const f = rep * m / d2;
        fx += f * dx; fy += f * dy;
      } else {
        const c = node * 4;
        const c0 = this.child[c], c1 = this.child[c + 1];
        const c2 = this.child[c + 2], c3 = this.child[c + 3];
        if (c0 >= 0) stack[sp++] = c0;
        if (c1 >= 0) stack[sp++] = c1;
        if (c2 >= 0) stack[sp++] = c2;
        if (c3 >= 0) stack[sp++] = c3;
      }
    }
    out[0] = fx; out[1] = fy;
  }
}

/**
 * Shelf-pack boxes into rows, largest first, targeting an aspect ratio.
 *
 * Writes `x`/`y` (top-left) into each box and returns the packed extent. Used
 * wherever whole components have to be placed: a uniform grid sized by the
 * largest component gives a two-segment plasmid the same cell as a
 * three-hundred-segment chromosome, which on a calibrated scale is thousands of
 * empty world units between neighbours.
 */
function shelfPack(boxes, aspect) {
  if (!boxes.length) return { w: 0, h: 0 };
  boxes.sort((a, b) => (b.w * b.h) - (a.w * a.h));
  let area = 0;
  for (const b of boxes) area += b.w * b.h;
  let rowWidth = Math.sqrt(area * (aspect || 1.4));
  for (const b of boxes) if (b.w > rowWidth) rowWidth = b.w;

  let cx = 0, cy = 0, rowH = 0, maxX = 0;
  for (const b of boxes) {
    if (cx > 0 && cx + b.w > rowWidth) { cx = 0; cy += rowH; rowH = 0; }
    b.x = cx;
    b.y = cy;
    cx += b.w;
    if (cx > maxX) maxX = cx;
    if (b.h > rowH) rowH = b.h;
  }
  return { w: maxX, h: cy + rowH };
}

/* ====================================================================== */
/*  Layout engine                                                          */
/* ====================================================================== */

/** Segment count the `repulsion` parameter is calibrated against. */
export const REFERENCE_SEGMENTS = 10;

/**
 * How much to damp repulsion for a graph of `n` segments.
 *
 * Repulsion accumulates over every particle while a particle's link forces do
 * not, so one value cannot serve every graph: undamped, the setting that lays
 * out a 10-segment graph turned a 514-segment one into a hairball with links
 * fifteen times longer than the segments they join.
 *
 * This is keyed to the **segment** count rather than the particle count on
 * purpose. Particles per segment is a rendering choice -- changing it should
 * not silently change the physics, which is exactly what happened when the
 * polyline resolution was reduced and repulsion jumped seventeen-fold.
 */
export function repulsionScaleFor(n) {
  if (!n || n <= REFERENCE_SEGMENTS) return 1;
  // Square root, not 1.5: measured across a 10-segment and a 125-segment
  // graph, a steeper law over-damped the large one or blew up the small one.
  return Math.max(1e-3, Math.pow(REFERENCE_SEGMENTS / n, 0.5));
}

export const DEFAULT_PARAMS = Object.freeze({
  // Repulsion accumulates over every particle, so a graph with many segments
  // naturally claims more room. These values were tuned so a small bacterial
  // graph settles at roughly a few times its longest segment rather than
  // flinging its components kilometres apart.
  // Calibrated by measuring the median link end-gap over real graphs: Bandage
  // gives a link an ideal length of a quarter of one polyline step, and this is
  // the value that reproduces that on both a 10-segment and a 285-segment
  // graph. The previous 0.05 predated the per-graph length calibration and left
  // small graphs with gaps twice the length of the contigs they joined.
  repulsion: 0.01,
  linkStrength: 1.80,
  // Rest length of a link, **as a fraction of the mean intra-segment spacing**.
  // Bandage gives real edges an ideal length of 5.0 against a node segment
  // length of 20.0, so a link is a quarter of one polyline step and joined
  // contigs sit end to end. An absolute value cannot work: once the drawn scale
  // is calibrated per graph, a fixed 24 units is longer than an entire contig
  // on a large assembly, and the drawing becomes dots joined by long leaders.
  linkRestFrac: 0.25,
  // Keeps a long contig reading as a smooth sweep rather than a squiggle.
  bendStrength: 0.90,
  // How firmly a closed molecule is held to a circle. Strong enough to survive
  // repulsion and the link springs, weak enough that dragging a contig still
  // deforms the ring under the pointer rather than fighting it.
  ringStrength: 0.35,
  // Off. Straightening the joins between contigs looks like the way to make a
  // chain flow, and it is a trap: the term has no length limit, so it pulls
  // whole chains straight and the drawing stretches into a line. Measured on
  // the demo graph, even 0.05 took the extent from 952 to 3777 units and the
  // aspect ratio from 1.0 to 5.1. Smoothness across joins is a property of the
  // drawing, not of the layout, and is handled in the renderer instead.
  jointStrength: 0.35,
  // Turns gentler than this (degrees) are left exactly as the layout put them.
  jointRelaxAbove: 55,
  // Curvature stiffness, applied over triples of consecutive vertices. This is
  // what lets a polyline carry many vertices without buckling, so long contigs
  // can bend into sweeping curves instead of being held straight by starving
  // them of vertices.
  curvature: 0.55,
  // Pull towards the component's own centre, not the whole drawing's. A global
  // centroid packs every component into one disc, which is why unrelated
  // plasmids used to be threaded through the chromosome.
  gravity: 0.02,
  damping: 0.82,
  theta: 0.75,
  maxIter: 600,
  // Width-to-height ratio aimed for when packing components into rows.
  packAspect: 1.4,
  constraintPasses: 8,
  animFrames: 42,
  // The three constants below are world units, and world units now mean
  // something fixed: the calibration puts the mean node at 40 and one polyline
  // step at 20 whatever the assembly's size (GraphModel.calibrate). They were
  // chosen against the old scale, where a bacterial contig came out ~5 units
  // long and a 70-unit rank gap was fourteen contigs of empty space.
  // `componentPad` is Bandage's own `componentSeparation` (settings.cpp:37);
  // the other two are one polyline step and a little under two of them, which
  // leaves a rearranged graph a few times looser than the force layout instead
  // of an order of magnitude looser.
  componentPad: 50,
  rankGap: 20,
  rowGap: 36,
});

export const LAYOUT_MODES = ['force', 'linear', 'circular', 'grid'];

export class LayoutEngine {
  constructor() {
    this.ready = false;
    this.tree = new QuadTree();
    this._acc = new Float64Array(2);
    this.mode = 'force';
    this.iter = 0;
    this.alpha = 0;
    this.running = false;
  }

  /** Install the structural arrays. Called once per graph revision. */
  init(d) {
    this.nParticles = d.nParticles | 0;
    this.nSegments = d.nSegments | 0;
    this.px = d.px instanceof Float32Array ? d.px : Float32Array.from(d.px || []);
    this.py = d.py instanceof Float32Array ? d.py : Float32Array.from(d.py || []);
    this.vx = new Float32Array(this.nParticles);
    this.vy = new Float32Array(this.nParticles);
    this.fx = new Float32Array(this.nParticles);
    this.fy = new Float32Array(this.nParticles);
    this.mobile = new Uint8Array(this.nParticles);
    this._preX = new Float32Array(this.nParticles);
    this._preY = new Float32Array(this.nParticles);
    this.segP0 = d.segP0; this.segK = d.segK;
    this.segRest = d.segRest; this.segLen = d.segLen; this.segComp = d.segComp;
    // A closed contig's polyline is a loop: its last vertex neighbours its
    // first. Without that the ring is held by nothing and the curvature term
    // straightens it back into a strand.
    this.segClosed = d.segClosed || new Uint8Array(this.nSegments);
    // Per segment: is it part of a closed molecule, and what radius should
    // that ring have? Folded into per-component arrays by _buildComponentMap,
    // which owns the dense component numbering.
    this.segRing = d.segRing || null;
    this.segRingRadius = d.segRingRadius || null;
    this.linkFrom = d.linkFrom; this.linkTo = d.linkTo;
    this.linkFromEnd = d.linkFromEnd; this.linkToEnd = d.linkToEnd;
    this.nLinks = this.linkFrom ? this.linkFrom.length : 0;
    // Mean particle spacing sets the natural length scale of the whole model.
    let sum = 0;
    for (let i = 0; i < this.nSegments; i++) sum += this.segRest[i];
    this.unit = this.nSegments ? Math.max(4, sum / this.nSegments) : 25;
    this.repulsionScale = repulsionScaleFor(this.nSegments);
    this._buildComponentMap();
    this.ready = this.nParticles > 0;
    this.running = false;
    return this.ready;
  }

  /** Overwrite positions (after the user drags nodes around). */
  setPositions(px, py) {
    if (!this.ready) return;
    const n = Math.min(this.nParticles, px.length, py.length);
    for (let i = 0; i < n; i++) { this.px[i] = px[i]; this.py[i] = py[i]; }
    this.vx.fill(0); this.vy.fill(0);
  }

  particleOf(seg, whichEnd) {
    return whichEnd === END_END ? this.segP0[seg] + this.segK[seg] - 1 : this.segP0[seg];
  }

  /** The vertex one step inward from an end; the end itself if there is none. */
  innerOf(seg, whichEnd) {
    const k = this.segK[seg];
    if (k < 2) return this.segP0[seg];
    return whichEnd === END_END ? this.segP0[seg] + k - 2 : this.segP0[seg] + 1;
  }

  /**
   * Dense particle -> component index, plus the scratch accumulators the
   * per-component gravity pass needs. Component ids in `segComp` are whatever
   * the server assigned, so they are compacted to 0..n-1 here once rather than
   * looked up every iteration.
   */
  _buildComponentMap() {
    const dense = new Map();
    this._particleComp = new Int32Array(this.nParticles);
    for (let s = 0; s < this.nSegments; s++) {
      const raw = this.segComp ? this.segComp[s] : 0;
      let c = dense.get(raw);
      if (c === undefined) { c = dense.size; dense.set(raw, c); }
      const p0 = this.segP0[s], k = this.segK[s];
      for (let i = 0; i < k; i++) this._particleComp[p0 + i] = c;
    }
    this._nComp = Math.max(1, dense.size);
    this._gcx = new Float64Array(this._nComp);
    this._gcy = new Float64Array(this._nComp);
    this._gcn = new Int32Array(this._nComp);

    // Fold the per-segment ring flags onto the dense component numbering.
    this.compRing = new Uint8Array(this._nComp);
    this.compRadius = new Float64Array(this._nComp);
    if (this.segRing && this.segRingRadius) {
      for (let s = 0; s < this.nSegments; s++) {
        if (!this.segRing[s]) continue;
        const c = this._particleComp[this.segP0[s]];
        this.compRing[c] = 1;
        this.compRadius[c] = this.segRingRadius[s];
      }
    }
  }

  _componentCount() { return this._nComp || 1; }

  _compOfParticle(i) { return this._particleComp ? this._particleComp[i] : 0; }

  /**
   * Pack connected components into tidy rows instead of letting repulsion push
   * them apart forever.
   *
   * Repulsion between components has nothing to balance it once gravity is
   * per-component, so unrelated plasmids drift away without bound and the
   * interesting component ends up a speck. Bandage sidesteps this by laying
   * each component out separately and packing them with a fixed separation
   * (`minDistCC`); this does the same, rigidly translating each component so
   * its shape -- which the force model owns -- is never disturbed.
   *
   * Largest component first, then shelf-packing into rows whose width targets
   * the view's aspect ratio, which reproduces Bandage's familiar "big component
   * on top, small ones in rows underneath" arrangement.
   */
  _packComponents() {
    const nComp = this._componentCount();
    if (nComp < 2) return;
    const { px, py, nParticles } = this;

    const minX = new Float64Array(nComp).fill(Infinity);
    const minY = new Float64Array(nComp).fill(Infinity);
    const maxX = new Float64Array(nComp).fill(-Infinity);
    const maxY = new Float64Array(nComp).fill(-Infinity);
    for (let i = 0; i < nParticles; i++) {
      const c = this._compOfParticle(i);
      const x = px[i], y = py[i];
      if (x < minX[c]) minX[c] = x;
      if (x > maxX[c]) maxX[c] = x;
      if (y < minY[c]) minY[c] = y;
      if (y > maxY[c]) maxY[c] = y;
    }

    const pad = Math.max(this.unit * 2.5, this.params.componentPad || 0);
    const boxes = [];
    for (let c = 0; c < nComp; c++) {
      if (!Number.isFinite(minX[c])) continue;
      boxes.push({ c, w: maxX[c] - minX[c] + pad, h: maxY[c] - minY[c] + pad });
    }
    if (boxes.length < 2) return;
    shelfPack(boxes, this.params.packAspect);

    // Translate each component from where it sits to its slot.
    const dx = new Float64Array(nComp);
    const dy = new Float64Array(nComp);
    for (const b of boxes) {
      dx[b.c] = b.x - minX[b.c] + pad / 2;
      dy[b.c] = b.y - minY[b.c] + pad / 2;
    }
    for (let i = 0; i < nParticles; i++) {
      const c = this._compOfParticle(i);
      px[i] += dx[c];
      py[i] += dy[c];
    }
  }

  /**
   * Prepare a run.
   * @param {string} mode one of LAYOUT_MODES
   * @param {Int32Array|null} subset mobile segment indices (null = everything)
   * @param {object} params  DEFAULT_PARAMS overrides, plus an optional
   *   `pinned` array of particle indices to hold still for the whole run
   */
  begin(mode, subset, params) {
    if (!this.ready) return false;
    this.params = { ...DEFAULT_PARAMS, ...(params || {}) };
    this.mode = LAYOUT_MODES.includes(mode) ? mode : 'force';
    this.iter = 0;
    this.settledFor = 0;

    // Mobility mask
    this.mobile.fill(0);
    let segList;
    if (subset && subset.length) {
      segList = Array.from(subset).filter((s) => s >= 0 && s < this.nSegments);
    } else {
      segList = new Array(this.nSegments);
      for (let i = 0; i < this.nSegments; i++) segList[i] = i;
    }
    if (!segList.length) return false;
    for (const s of segList) {
      const p0 = this.segP0[s], k = this.segK[s];
      for (let i = 0; i < k; i++) this.mobile[p0 + i] = 1;
    }
    // Pinned particles are held still even though their segment is mobile. This
    // is what lets a dragged vertex stay under the pointer while the constraint
    // solver relaxes the rest of its polyline -- and the graph -- around it.
    const pinned = this.params.pinned;
    // Component packing translates whole components rigidly, immobile
    // particles included, so it would drag a pinned vertex out from under the
    // pointer. A run with pins leaves the components where they are.
    this.packing = !(pinned && pinned.length);
    if (!this.packing) {
      for (let i = 0; i < pinned.length; i++) {
        // `| 0` alone would turn a caller's undefined into particle 0 and
        // silently freeze the first vertex of the graph.
        const p = Number(pinned[i]);
        if (Number.isInteger(p) && p >= 0 && p < this.nParticles) this.mobile[p] = 0;
      }
    }
    this.segList = segList;

    if (this.mode === 'force') {
      this.vx.fill(0); this.vy.fill(0);
      this.alpha = 1;
      this.maxIter = Math.max(20, this.params.maxIter | 0);
    } else {
      const targets = this._computeTargets(this.mode, segList);
      if (!targets) return false;
      this.tx = targets.x; this.ty = targets.y;
      this.fromX = Float32Array.from(this.px);
      this.fromY = Float32Array.from(this.py);
      this.frames = Math.max(1, this.params.animFrames | 0);
      this.alpha = 1;
    }
    this.running = true;
    return true;
  }

  stop() { this.running = false; }

  /** Advance one iteration/frame. Returns {done, iter, alpha}. */
  step() {
    if (!this.running || !this.ready) return { done: true, iter: this.iter, alpha: 0 };
    if (this.mode === 'force') return this._stepForce();
    return this._stepAnimate();
  }

  /* ------------------------------------------------------------- force */

  _stepForce() {
    const {
      px, py, vx, vy, fx, fy, mobile, nParticles,
      segP0, segK, segRest, nSegments,
    } = this;
    const P = this.params;
    const unit = this.unit;
    const alpha = this.alpha;

    fx.fill(0); fy.fill(0);

    /* --- 1. Barnes-Hut repulsion --------------------------------------- */
    // Repulsion accumulates over every particle, while a particle's link forces
    // do not grow with the graph. Left unnormalised, the same `repulsion` value
    // that lays out a 10-node graph nicely blows a 500-node graph into a
    // hairball with links 15x longer than the segments they join. Normalising
    // by particle count makes the parameter mean the same thing at any scale;
    // the exponent and reference count were fitted by measuring mean link
    // length (in units of mean segment length) across real graphs.
    const rep = P.repulsion * (1.8 * unit) * (1.8 * unit) * this.repulsionScale;
    const theta2 = P.theta * P.theta;
    this.tree.build(px, py, nParticles);
    const acc = this._acc;
    for (let i = 0; i < nParticles; i++) {
      if (!mobile[i]) continue;
      this.tree.force(i, px[i], py[i], theta2, rep, acc);
      fx[i] += acc[0]; fy[i] += acc[1];
    }

    /* --- 2. link springs ----------------------------------------------- */
    const ks = P.linkStrength;
    const rest = unit * (P.linkRestFrac !== undefined ? P.linkRestFrac : 0.25);
    for (let l = 0; l < this.nLinks; l++) {
      // A link from a segment to itself is a circular contig closing up. Its
      // two particles are already held a fixed distance apart by the polyline's
      // own distance constraints, so a spring pulling them to `rest` fights the
      // constraint solver forever: the pair buzzes at full speed, `moved` never
      // falls, and the layout burns its whole iteration budget without ever
      // settling. Bandage drops these from the layout graph for the same
      // reason; the renderer draws the loop instead.
      if (this.linkFrom[l] === this.linkTo[l]) continue;
      const a = this.particleOf(this.linkFrom[l], this.linkFromEnd[l]);
      const b = this.particleOf(this.linkTo[l], this.linkToEnd[l]);
      if (a === b) continue;
      const dx = px[b] - px[a];
      const dy = py[b] - py[a];
      const d = Math.sqrt(dx * dx + dy * dy) || 1e-3;
      const f = ks * (d - rest) / d;
      const ux = dx * f, uy = dy * f;
      if (mobile[a]) { fx[a] += ux; fy[a] += uy; }
      if (mobile[b]) { fx[b] -= ux; fy[b] -= uy; }
    }

    /* --- 3. bending springs (keep polylines straight-ish) --------------- */
    const kb = P.bendStrength;
    if (kb > 0) {
      for (let s = 0; s < nSegments; s++) {
        const p0 = segP0[s], k = segK[s];
        if (k < 3) continue;
        const r2 = segRest[s] * 2;
        for (let i = 0; i + 2 < k; i++) {
          const a = p0 + i, b = p0 + i + 2;
          const dx = px[b] - px[a];
          const dy = py[b] - py[a];
          const d = Math.sqrt(dx * dx + dy * dy) || 1e-3;
          const f = kb * (d - r2) / d;
          const ux = dx * f, uy = dy * f;
          if (mobile[a]) { fx[a] += ux; fy[a] += uy; }
          if (mobile[b]) { fx[b] -= ux; fy[b] -= uy; }
        }
      }
    }

    /* --- 3b. curvature stiffness --------------------------------------- */
    // Penalise the turn at each interior vertex by pulling it towards the
    // midpoint of its neighbours. Unlike the i<->i+2 distance spring this acts
    // on the *shape* rather than the span, so a chain stays smooth however many
    // vertices it has -- which is what lets long contigs sweep instead of
    // buckling into a zigzag. Bandage gets the same effect from FMMM's
    // multilevel solve; this is the cheap local equivalent.
    // Pull the middle of a triple towards the midpoint of its neighbours, and
    // push the neighbours back by half each so the term cannot translate the
    // chain as a whole -- only straighten it.
    const smoothTriple = (a, b, c, k) => {
      const mx = (px[a] + px[c]) * 0.5;
      const my = (py[a] + py[c]) * 0.5;
      const ux = (mx - px[b]) * k;
      const uy = (my - py[b]) * k;
      if (mobile[b]) { fx[b] += ux * 2; fy[b] += uy * 2; }
      if (mobile[a]) { fx[a] -= ux; fy[a] -= uy; }
      if (mobile[c]) { fx[c] -= ux; fy[c] -= uy; }
    };

    const kc = P.curvature || 0;
    if (kc > 0) {
      for (let s = 0; s < nSegments; s++) {
        const p0 = segP0[s], k = segK[s];
        if (k < 3) continue;
        // A loop has a turn at every vertex, including the two either side of
        // the seam; skipping those leaves a visible corner where the ring
        // closes.
        const closed = this.segClosed[s];
        const first = closed ? 0 : 1;
        const last = closed ? k - 1 : k - 2;
        for (let i = first; i <= last; i++) {
          smoothTriple(p0 + (i - 1 + k) % k, p0 + i, p0 + (i + 1) % k, kc);
        }
      }
    }

    /* --- 3b2. straighten the joins between contigs ---------------------- */
    // The curvature term above only smooths a contig's own vertices, and a
    // typical contig has three of them, so it is nearly a straight bar. What
    // the eye actually follows is the chain of bars, and nothing wanted
    // consecutive bars to line up: joins turned through 62 degrees on average
    // and the worst tenth doubled back, which is what made the drawing look
    // shaky rather than fluid.
    //
    // Applying the same curvature over the triples that span a link -- the
    // vertex before the join, the join, and the vertex after -- lets a run of
    // contigs settle into one sweeping curve. Where several contigs meet, the
    // terms compete and balance, which is the fan a branch point should have.
    const kj = P.jointStrength || 0;
    if (kj > 0) {
      // Only the kinks. Smoothing every join pulls whole chains straight, and a
      // straight chain takes more room than a curled one, so the drawing
      // stretches: at full strength on every join the demo graph went from 952
      // units across to 6522 and from square to a 10:1 streak. Leaving turns
      // gentler than `jointRelaxAbove` alone keeps the layout's own shape and
      // spends the force only where the eye actually catches a corner.
      const cosLimit = Math.cos((P.jointRelaxAbove || 50) * Math.PI / 180);
      const relax = (a, b, c) => {
        const ux = px[b] - px[a], uy = py[b] - py[a];
        const vx2 = px[c] - px[b], vy2 = py[c] - py[b];
        const lu = Math.hypot(ux, uy), lv = Math.hypot(vx2, vy2);
        if (lu < 1e-6 || lv < 1e-6) return;
        const cos = (ux * vx2 + uy * vy2) / (lu * lv);
        if (cos >= cosLimit) return; // already gentle enough to leave alone
        // Ramp in over the sharper half of the range so there is no step at the
        // threshold, which would make contigs flicker between two shapes.
        smoothTriple(a, b, c, kj * Math.min(1, (cosLimit - cos) / 1.2));
      };
      for (let l = 0; l < this.nLinks; l++) {
        if (this.linkFrom[l] === this.linkTo[l]) continue;
        const a = this.particleOf(this.linkFrom[l], this.linkFromEnd[l]);
        const b = this.particleOf(this.linkTo[l], this.linkToEnd[l]);
        const ai = this.innerOf(this.linkFrom[l], this.linkFromEnd[l]);
        const bi = this.innerOf(this.linkTo[l], this.linkToEnd[l]);
        if (ai === a || bi === b) continue;
        relax(ai, a, b);
        relax(a, b, bi);
      }
    }

    /* --- 3c. hold closed molecules open as rings ------------------------ */
    // A circular molecule has no intrinsic 2D shape -- every drawing of the
    // loop is equally true -- so the layout may as well pick the one that says
    // "this is circular" at a glance. Left to itself a ring of contigs is
    // neutrally stable: bending it costs nothing, so it crumples into a blob
    // that hides the one fact worth showing. This pulls each vertex towards the
    // radius the loop's own length implies, which restores the shape without
    // pinning anything: links, repulsion and dragging all still act.
    const kr = P.ringStrength || 0;
    if (kr > 0 && this.compRing && this.compRadius) {
      const nComp = this._componentCount();
      const cxs = this._gcx, cys = this._gcy, cns = this._gcn;
      cxs.fill(0); cys.fill(0); cns.fill(0);
      for (let i = 0; i < nParticles; i++) {
        const c = this._compOfParticle(i);
        cxs[c] += px[i]; cys[c] += py[i]; cns[c]++;
      }
      for (let c = 0; c < nComp; c++) if (cns[c]) { cxs[c] /= cns[c]; cys[c] /= cns[c]; }
      for (let i = 0; i < nParticles; i++) {
        if (!mobile[i]) continue;
        const c = this._compOfParticle(i);
        if (!this.compRing[c]) continue;
        const target = this.compRadius[c];
        if (!(target > 0)) continue;
        const dx = px[i] - cxs[c];
        const dy = py[i] - cys[c];
        const r = Math.hypot(dx, dy);
        if (r < 1e-6) continue;
        const f = kr * (target - r) / r;
        fx[i] += dx * f;
        fy[i] += dy * f;
      }
    }

    /* --- 4. gravity towards each component's own centroid --------------- */
    const g = P.gravity;
    if (g > 0) {
      const nComp = this._componentCount();
      const cxs = this._gcx, cys = this._gcy, cns = this._gcn;
      cxs.fill(0); cys.fill(0); cns.fill(0);
      for (let i = 0; i < nParticles; i++) {
        if (!mobile[i]) continue;
        const c = this._compOfParticle(i);
        cxs[c] += px[i]; cys[c] += py[i]; cns[c]++;
      }
      for (let c = 0; c < nComp; c++) {
        if (cns[c]) { cxs[c] /= cns[c]; cys[c] /= cns[c]; }
      }
      for (let i = 0; i < nParticles; i++) {
        if (!mobile[i]) continue;
        const c = this._compOfParticle(i);
        if (!cns[c]) continue;
        fx[i] += (cxs[c] - px[i]) * g;
        fy[i] += (cys[c] - py[i]) * g;
      }
    }

    /* --- 5. integrate --------------------------------------------------- */
    const damp = P.damping;
    const maxStep = unit * 0.6 * alpha + 0.35;
    const maxForce = unit * 20;
    let moved = 0;
    for (let i = 0; i < nParticles; i++) {
      if (!mobile[i]) continue;
      let ax = fx[i], ay = fy[i];
      if (ax > maxForce) ax = maxForce; else if (ax < -maxForce) ax = -maxForce;
      if (ay > maxForce) ay = maxForce; else if (ay < -maxForce) ay = -maxForce;
      let nvx = (vx[i] + ax) * damp;
      let nvy = (vy[i] + ay) * damp;
      if (nvx > maxStep) nvx = maxStep; else if (nvx < -maxStep) nvx = -maxStep;
      if (nvy > maxStep) nvy = maxStep; else if (nvy < -maxStep) nvy = -maxStep;
      vx[i] = nvx; vy[i] = nvy;
      px[i] += nvx; py[i] += nvy;
      const m = Math.abs(nvx) + Math.abs(nvy);
      if (m > moved) moved = m;
    }

    /* --- 6. stiff intra-segment bonds via position constraints ---------- */
    // The projection moves particles without touching their velocities, so the
    // displacement it undoes comes straight back next step as kinetic energy.
    // Recording the pre-projection positions lets the correction be folded back
    // into the velocities, which is what stops long chains from shivering.
    const passes = Math.max(0, P.constraintPasses | 0);
    const preX = passes > 0 ? this._preX : null;
    const preY = passes > 0 ? this._preY : null;
    if (preX) {
      for (let i = 0; i < nParticles; i++) { preX[i] = px[i]; preY[i] = py[i]; }
    }
    for (let pass = 0; pass < passes; pass++) {
      for (let s = 0; s < nSegments; s++) {
        const p0 = segP0[s], k = segK[s];
        if (k < 2) continue;
        const r = segRest[s];
        // One more bond on a closed contig, joining its last vertex back to its
        // first. That single constraint is what holds a circular molecule shut:
        // its self-link is deliberately kept out of the spring set, so without
        // this the ring has nothing to close it and unrolls into a strand.
        const bonds = this.segClosed[s] ? k : k - 1;
        for (let i = 0; i < bonds; i++) {
          const a = p0 + i, b = p0 + ((i + 1) % k);
          const ma = mobile[a], mb = mobile[b];
          if (!ma && !mb) continue;
          let dx = px[b] - px[a];
          let dy = py[b] - py[a];
          let d = Math.sqrt(dx * dx + dy * dy);
          if (d < 1e-4) { dx = r; dy = 0; d = r; }
          const corr = (d - r) / d;
          if (ma && mb) {
            const hx = dx * corr * 0.5, hy = dy * corr * 0.5;
            px[a] += hx; py[a] += hy;
            px[b] -= hx; py[b] -= hy;
          } else if (ma) {
            px[a] += dx * corr; py[a] += dy * corr;
          } else {
            px[b] -= dx * corr; py[b] -= dy * corr;
          }
        }
      }
    }

    if (preX) {
      // Fold the projection's displacement into the velocities rather than
      // discarding it, and re-measure `moved` on the positions the drawing will
      // actually use, so the settle test reflects the final state.
      moved = 0;
      for (let i = 0; i < nParticles; i++) {
        if (!mobile[i]) continue;
        const cx = px[i] - preX[i];
        const cy = py[i] - preY[i];
        vx[i] += cx;
        vy[i] += cy;
        const m = Math.abs(vx[i]) + Math.abs(vy[i]);
        if (m > moved) moved = m;
      }
    }

    /* --- 6b. keep components packed ------------------------------------ */
    // Every few steps rather than every step: packing is a rigid translation so
    // it never disturbs a component's shape, but doing it constantly makes the
    // components visibly twitch while the sim is still moving them.
    if (this.packing && this._nComp > 1 && (this.iter % 4) === 0) this._packComponents();

    /* --- 7. cool ------------------------------------------------------- */
    this.iter++;
    const t = this.iter / this.maxIter;
    this.alpha = t >= 1 ? 0 : Math.pow(1 - t, 1.4);
    // Settling is judged relative to the drawing's own scale. An absolute
    // threshold means a graph calibrated to large world units never looks
    // settled and always burns the full iteration budget.
    if (moved < unit * 0.0025) this.settledFor++; else this.settledFor = 0;
    const done = this.iter >= this.maxIter || this.settledFor > 24;
    if (done) {
      if (this.packing && this._nComp > 1) this._packComponents();
      this.running = false;
    }
    return { done, iter: this.iter, alpha: this.alpha };
  }

  /* --------------------------------------------------- animated modes */

  _stepAnimate() {
    this.iter++;
    const t = Math.min(1, this.iter / this.frames);
    // smoothstep for a calm ease-in-out
    const e = t * t * (3 - 2 * t);
    const { px, py, tx, ty, fromX, fromY, mobile, nParticles } = this;
    for (let i = 0; i < nParticles; i++) {
      if (!mobile[i]) continue;
      px[i] = fromX[i] + (tx[i] - fromX[i]) * e;
      py[i] = fromY[i] + (ty[i] - fromY[i]) * e;
    }
    this.alpha = 1 - e;
    const done = t >= 1;
    if (done) this.running = false;
    return { done, iter: this.iter, alpha: this.alpha };
  }

  /* ------------------------------------------------ deterministic modes */

  /** Adjacency restricted to `segList` (indices into segList, not segments). */
  _inducedAdjacency(segList) {
    const pos = new Int32Array(this.nSegments).fill(-1);
    segList.forEach((s, i) => { pos[s] = i; });
    const adj = segList.map(() => []);
    for (let l = 0; l < this.nLinks; l++) {
      const a = pos[this.linkFrom[l]];
      const b = pos[this.linkTo[l]];
      if (a < 0 || b < 0) continue;
      adj[a].push({ o: b, selfEnd: this.linkFromEnd[l], otherEnd: this.linkToEnd[l] });
      if (a !== b) adj[b].push({ o: a, selfEnd: this.linkToEnd[l], otherEnd: this.linkFromEnd[l] });
    }
    return adj;
  }

  /** Connected components of the induced subgraph, as arrays of local indices. */
  _inducedComponents(segList, adj) {
    const seen = new Uint8Array(segList.length);
    const comps = [];
    const queue = new Int32Array(segList.length);
    for (let s = 0; s < segList.length; s++) {
      if (seen[s]) continue;
      let head = 0, tail = 0;
      queue[tail++] = s; seen[s] = 1;
      const comp = [];
      while (head < tail) {
        const u = queue[head++];
        comp.push(u);
        for (const e of adj[u]) {
          if (!seen[e.o]) { seen[e.o] = 1; queue[tail++] = e.o; }
        }
      }
      comps.push(comp);
    }
    // Biggest first so the eye starts at the important structure.
    comps.sort((a, b) => {
      let la = 0, lb = 0;
      for (const i of a) la += this.segLen[segList[i]];
      for (const i of b) lb += this.segLen[segList[i]];
      return lb - la;
    });
    return comps;
  }

  /**
   * BFS a component recording, for each segment, its rank offset from the root
   * and whether it must be drawn reversed so links meet head-to-tail.
   */
  _traverse(comp, adj, segList) {
    // Root: prefer a dead end (fewest connections), tie-break on length.
    let root = comp[0];
    for (const u of comp) {
      const du = adj[u].length, dr = adj[root].length;
      if (du < dr || (du === dr && this.segLen[segList[u]] > this.segLen[segList[root]])) root = u;
    }
    const rank = new Map();
    const flip = new Map();
    const order = [];
    const seen = new Set([root]);
    rank.set(root, 0); flip.set(root, false);
    const queue = [root];
    while (queue.length) {
      const u = queue.shift();
      order.push(u);
      for (const e of adj[u]) {
        if (seen.has(e.o)) continue;
        seen.add(e.o);
        const uFlip = flip.get(u);
        // Which drawn side of u does this link leave from?
        const leavesLeft = uFlip ? (e.selfEnd === END_END) : (e.selfEnd === END_START);
        const placedRight = !leavesLeft;
        rank.set(e.o, rank.get(u) + (placedRight ? 1 : -1));
        // Orient the neighbour so its attached end faces back towards u.
        flip.set(e.o, placedRight ? (e.otherEnd === END_END) : (e.otherEnd === END_START));
        queue.push(e.o);
      }
    }
    // Isolated members (shouldn't happen inside one component, but be safe).
    for (const u of comp) {
      if (!rank.has(u)) { rank.set(u, 0); flip.set(u, false); order.push(u); }
    }
    return { rank, flip, order };
  }

  /** Write a straight horizontal polyline for a segment into the target arrays. */
  _emitLine(seg, x, y, flip, tx, ty) {
    const p0 = this.segP0[seg], k = this.segK[seg], L = this.segLen[seg];
    const step = k > 1 ? L / (k - 1) : 0;
    for (let i = 0; i < k; i++) {
      tx[p0 + i] = flip ? (x + L - step * i) : (x + step * i);
      ty[p0 + i] = y;
    }
  }

  _computeTargets(mode, segList) {
    const tx = Float32Array.from(this.px);
    const ty = Float32Array.from(this.py);
    const adj = this._inducedAdjacency(segList);
    const comps = this._inducedComponents(segList, adj);
    if (!comps.length) return null;

    if (mode === 'linear') this._layoutLinear(segList, adj, comps, tx, ty);
    else if (mode === 'circular') this._layoutCircular(segList, adj, comps, tx, ty);
    else this._layoutGrid(segList, comps, tx, ty);

    this._recentre(segList, tx, ty);
    return { x: tx, y: ty };
  }

  /** Keep the rearranged block where it already was on screen. */
  _recentre(segList, tx, ty) {
    let ox = 0, oy = 0, nx = 0, ny = 0, n = 0;
    for (const s of segList) {
      const p0 = this.segP0[s], k = this.segK[s];
      for (let i = 0; i < k; i++) {
        ox += this.px[p0 + i]; oy += this.py[p0 + i];
        nx += tx[p0 + i]; ny += ty[p0 + i];
        n++;
      }
    }
    if (!n) return;
    const dx = (ox - nx) / n, dy = (oy - ny) / n;
    if (!dx && !dy) return;
    for (const s of segList) {
      const p0 = this.segP0[s], k = this.segK[s];
      for (let i = 0; i < k; i++) { tx[p0 + i] += dx; ty[p0 + i] += dy; }
    }
  }

  _layoutLinear(segList, adj, comps, tx, ty) {
    const P = this.params;
    let yCursor = 0;
    let prevHalf = null;
    for (const comp of comps) {
      const { rank, flip } = this._traverse(comp, adj, segList);
      // Normalise ranks to start at 0.
      let minRank = Infinity, maxRank = -Infinity;
      for (const u of comp) {
        const r = rank.get(u);
        if (r < minRank) minRank = r;
        if (r > maxRank) maxRank = r;
      }
      const nRanks = maxRank - minRank + 1;
      const rows = Array.from({ length: nRanks }, () => []);
      for (const u of comp) rows[rank.get(u) - minRank].push(u);
      // Column x positions from the widest segment in each rank.
      const colX = new Float64Array(nRanks);
      const colW = new Float64Array(nRanks);
      let x = 0;
      for (let r = 0; r < nRanks; r++) {
        let w = 0;
        for (const u of rows[r]) w = Math.max(w, this.segLen[segList[u]]);
        colW[r] = w;
        colX[r] = x;
        x += w + P.rankGap;
      }
      // Measure the row stack before placing anything: the gap between two
      // components is half of this one plus half of the *next* one, and
      // advancing by this one's height twice let a tall fan of short contigs
      // land straight on top of the thin chain above it.
      let height = 0;
      for (let r = 0; r < nRanks; r++) height = Math.max(height, rows[r].length * P.rowGap);
      yCursor += (prevHalf === null ? 0 : prevHalf + P.componentPad) + height / 2;
      for (let r = 0; r < nRanks; r++) {
        const list = rows[r];
        list.sort((a, b) => this.segLen[segList[b]] - this.segLen[segList[a]]);
        for (let i = 0; i < list.length; i++) {
          const u = list[i];
          const seg = segList[u];
          const L = this.segLen[seg];
          const px0 = colX[r] + (colW[r] - L) / 2;
          const py0 = yCursor + (i - (list.length - 1) / 2) * P.rowGap;
          this._emitLine(seg, px0, py0, flip.get(u), tx, ty);
        }
      }
      prevHalf = height / 2;
    }
  }

  _layoutCircular(segList, adj, comps, tx, ty) {
    const P = this.params;
    // The arc between consecutive nodes is a link, so it gets a link's ideal
    // length. An absolute floor here used to keep a whole graph's circles
    // apart by more than the nodes on them were long.
    const gap = this.unit * (P.linkRestFrac !== undefined ? P.linkRestFrac : 0.25);
    // Lay each component on its own circle, then pack the circles.
    const circles = comps.map((comp) => {
      const { order, flip } = this._traverse(comp, adj, segList);
      let circumference = 0;
      for (const u of order) circumference += this.segLen[segList[u]] + gap;
      // One short node still needs an arc it can be read along.
      const R = Math.max(this.unit, circumference / (2 * Math.PI));
      const d = R * 2 + P.componentPad;
      return { order, flip, R, w: d, h: d };
    });
    shelfPack(circles, P.packAspect);
    circles.forEach((c) => {
      const ox = c.x + c.w / 2;
      const oy = c.y + c.h / 2;
      let angle = 0;
      const R = c.R;
      for (const u of c.order) {
        const seg = segList[u];
        const p0 = this.segP0[seg], k = this.segK[seg], L = this.segLen[seg];
        const da = L / R;
        const isFlipped = c.flip.get(u);
        for (let i = 0; i < k; i++) {
          const t = k > 1 ? i / (k - 1) : 0;
          const a = angle + (isFlipped ? (1 - t) : t) * da;
          tx[p0 + i] = ox + Math.cos(a) * R;
          ty[p0 + i] = oy + Math.sin(a) * R;
        }
        angle += da + gap / R;
      }
    });
  }

  _layoutGrid(segList, comps, tx, ty) {
    const P = this.params;
    // Each component keeps its current shape; only its position changes.
    const boxes = comps.map((comp) => {
      let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
      for (const u of comp) {
        const seg = segList[u];
        const p0 = this.segP0[seg], k = this.segK[seg];
        for (let i = 0; i < k; i++) {
          const x = this.px[p0 + i], y = this.py[p0 + i];
          if (x < minX) minX = x;
          if (x > maxX) maxX = x;
          if (y < minY) minY = y;
          if (y > maxY) maxY = y;
        }
      }
      if (!Number.isFinite(minX)) { minX = minY = 0; maxX = maxY = 1; }
      return { comp, minX, minY, maxX, maxY, w: maxX - minX, h: maxY - minY };
    });
    const pad = P.componentPad;
    for (const b of boxes) { b.w += pad; b.h += pad; }
    shelfPack(boxes, P.packAspect);
    boxes.forEach((b) => {
      const dx = b.x + b.w / 2 - (b.minX + b.maxX) / 2;
      const dy = b.y + b.h / 2 - (b.minY + b.maxY) / 2;
      for (const u of b.comp) {
        const seg = segList[u];
        const p0 = this.segP0[seg], k = this.segK[seg];
        for (let i2 = 0; i2 < k; i2++) {
          tx[p0 + i2] = this.px[p0 + i2] + dx;
          ty[p0 + i2] = this.py[p0 + i2] + dy;
        }
      }
    });
  }
}

/* ====================================================================== */
/*  Main-thread controller                                                 */
/* ====================================================================== */

/**
 * Drives a LayoutEngine that normally lives in `js/layout-worker.js`
 * (`new Worker(url, {type: 'module'})`). If the browser refuses to create the
 * worker — e.g. the page was opened over file:// — the same engine is run in
 * small time-sliced chunks on the main thread instead, so the feature never
 * simply disappears.
 */
export class LayoutController {
  constructor(graph, handlers = {}) {
    this.graph = graph;
    this.onTick = handlers.onTick || (() => {});
    this.onDone = handlers.onDone || (() => {});
    this.onError = handlers.onError || (() => {});
    this.onState = handlers.onState || (() => {});
    this.worker = null;
    this.fallback = null;
    this.dirty = true;
    this.running = false;
    this.lastMode = 'force';
    this._fallbackTimer = 0;
  }

  /** Mark the structure as stale; the next start() re-uploads it. */
  markDirty() { this.dirty = true; }

  get usingWorker() { return !!this.worker; }

  _ensure() {
    if (this.worker || this.fallback) return;
    try {
      // Resolved against this module's URL, i.e. `js/layout-worker.js`.
      const url = new URL('./layout-worker.js', import.meta.url);
      this.worker = new Worker(url, { type: 'module' });
      this.worker.onmessage = (ev) => this._handle(ev.data);
      this.worker.onerror = (ev) => {
        // A broken worker must not kill the layout feature.
        const msg = (ev && ev.message) || 'layout worker failed';
        this.worker = null;
        this.fallback = new LayoutEngine();
        this.dirty = true;
        this.onError(msg + ' — falling back to main-thread layout');
      };
    } catch (e) {
      this.worker = null;
      this.fallback = new LayoutEngine();
      this.dirty = true;
      this.onError('Web Worker unavailable (' + ((e && e.message) || e) + '); layout runs on the main thread');
    }
  }

  _post(msg, transfer) {
    if (this.worker) this.worker.postMessage(msg, transfer || []);
  }

  _handle(msg) {
    if (!msg || !msg.type) return;
    // Discard anything belonging to a run we have already superseded. Pressing
    // Rearrange while a layout was running used to let the old run's 'done'
    // land after the new one had started, so the UI showed "settled" and hid
    // Stop while the worker was still iterating.
    if (msg.runId !== undefined && this._runId !== undefined && msg.runId !== this._runId) return;
    if (msg.type === 'tick' || msg.type === 'done') {
      const g = this.graph;
      if (msg.px && msg.py && msg.px.length === g.px.length) {
        g.px.set(msg.px);
        g.py.set(msg.py);
        g.updateBounds();
      }
      if (msg.type === 'tick') {
        this.onTick(msg);
      } else {
        this.running = false;
        this.onState(false);
        this.onDone(msg);
      }
    } else if (msg.type === 'error') {
      this.running = false;
      this.onState(false);
      this.onError(msg.message || 'layout error');
    }
  }

  /**
   * Start a layout run.
   * @param {object} o
   * @param {string} o.mode  force|linear|circular|grid
   * @param {Int32Array|number[]|null} o.subset mobile segment indices
   * @param {object} o.params
   */
  start(o = {}) {
    const g = this.graph;
    if (!g || g.isEmpty) { this.onError('Nothing to lay out — no graph is loaded'); return false; }
    this._ensure();
    this.stop();
    this._runId = (this._runId || 0) + 1;

    if (this.dirty) {
      const data = g.toLayoutArrays();
      if (this.worker) {
        this._post({ cmd: 'init', data }, [
          data.px.buffer, data.py.buffer, data.segP0.buffer, data.segK.buffer,
          data.segRest.buffer, data.segLen.buffer, data.segComp.buffer,
          data.linkFrom.buffer, data.linkTo.buffer,
          data.linkFromEnd.buffer, data.linkToEnd.buffer,
        ]);
      } else if (this.fallback) {
        this.fallback.init(data);
      }
      this.dirty = false;
    } else {
      // Positions may have changed on the main thread (drag, session restore).
      const px = Float32Array.from(g.px);
      const py = Float32Array.from(g.py);
      if (this.worker) this._post({ cmd: 'positions', px, py }, [px.buffer, py.buffer]);
      else if (this.fallback) this.fallback.setPositions(px, py);
    }

    const subset = o.subset && o.subset.length
      ? (o.subset instanceof Int32Array ? o.subset : Int32Array.from(o.subset))
      : null;
    const mode = o.mode || 'force';
    this.lastMode = mode;
    this.running = true;
    this.onState(true);

    if (this.worker) {
      const msg = { cmd: 'run', mode, subset, params: o.params || {}, runId: this._runId };
      this._post(msg, subset ? [subset.buffer] : []);
    } else if (this.fallback) {
      const ok = this.fallback.begin(mode, subset, o.params || {});
      if (!ok) {
        this.running = false;
        this.onState(false);
        this.onError('Layout could not start (empty selection?)');
        return false;
      }
      this._runFallback();
    }
    return true;
  }

  _runFallback() {
    const engine = this.fallback;
    const tick = () => {
      if (!this.running || !engine) return;
      const t0 = (typeof performance !== 'undefined' ? performance.now() : Date.now());
      let res = { done: false, iter: engine.iter, alpha: engine.alpha };
      do {
        res = engine.step();
      } while (!res.done && ((typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0) < 10);
      this._handle({
        type: res.done ? 'done' : 'tick',
        px: engine.px, py: engine.py,
        iter: res.iter, alpha: res.alpha, mode: engine.mode,
      });
      if (!res.done && this.running) {
        this._fallbackTimer = setTimeout(tick, 0);
      }
    };
    this._fallbackTimer = setTimeout(tick, 0);
  }

  stop() {
    if (!this.running) return;
    this.running = false;
    if (this._fallbackTimer) { clearTimeout(this._fallbackTimer); this._fallbackTimer = 0; }
    if (this.fallback) this.fallback.stop();
    this._post({ cmd: 'stop', runId: this._runId });
    this.onState(false);
  }

  /** Push main-thread position edits (node dragging) into the engine. */
  syncPositions() {
    const g = this.graph;
    if (!g || g.isEmpty) return;
    if (this.worker) {
      const px = Float32Array.from(g.px);
      const py = Float32Array.from(g.py);
      this._post({ cmd: 'positions', px, py }, [px.buffer, py.buffer]);
    } else if (this.fallback && this.fallback.ready) {
      this.fallback.setPositions(g.px, g.py);
    }
  }

  destroy() {
    this.stop();
    if (this.worker) { this.worker.terminate(); this.worker = null; }
    this.fallback = null;
  }
}

export default LayoutController;
