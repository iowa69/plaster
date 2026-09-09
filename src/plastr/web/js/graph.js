/**
 * graph.js — the client-side assembly graph model.
 *
 * A segment is *not* a dot. Following Bandage's convention (and the design
 * spec), every segment is drawn as a polyline made of `k` layout particles with
 * a distinguished *start* end (5' of the `+` strand) and *end* end. Links
 * attach to one of those two ends depending on the orientation signs, per the
 * table in docs/API.md:
 *
 *   from_orient | to_orient | joins
 *   ------------+-----------+------------------
 *        +      |     +     | A.end   -> B.start
 *        +      |     -     | A.end   -> B.end
 *        -      |     +     | A.start -> B.start
 *        -      |     -     | A.start -> B.end
 *
 * i.e. `from_orient === '+'` means the link *leaves* the END particle, and
 * `to_orient === '+'` means it *arrives* at the START particle.
 *
 * Positions live in two flat Float32Arrays (px/py) indexed by particle so they
 * can be handed to the layout worker without any per-object marshalling.
 */

/* ------------------------------------------------------------------ utils */

export function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

/** End identifiers used everywhere in the codebase. */
export const END_START = 0;
export const END_END = 1;

/**
 * Contiguity classes, in Bandage's precedence order: the *lower* value always
 * wins, because a segment reached by several routes keeps the strongest verdict
 * it earned (`DeBruijnNode::upgradeContiguityStatus`, debruijnnode.cpp:312).
 */
export const CONTIGUITY = Object.freeze({
  STARTING: 0,
  CONTIGUOUS: 1,
  MAYBE_CONTIGUOUS: 2,
  NOT_CONTIGUOUS: 3,
});

/** User-facing names, indexed by the values above. */
export const CONTIGUITY_LABELS = Object.freeze([
  'Starting', 'Contiguous', 'Maybe contiguous', 'Not contiguous',
]);

// Bandage keeps five internal levels (globals.h:38). The two contiguous ones
// are painted the same green (settings.cpp:101-102), so they collapse into
// CONTIGUITY.CONTIGUOUS on the way out; the search has to keep them apart
// because the either-strand test is only run for oriented nodes the
// strand-specific test did not already claim.
const ST_STARTING = 0;
const ST_CONTIG_STRAND = 1;
const ST_CONTIG_EITHER = 2;
const ST_MAYBE = 3;
const ST_NOT = 4;

/** Geometry defaults; the UI can override `lengthScale` (node size slider). */
export const DEFAULT_GEOM = Object.freeze({
  minLen: 5,           // shortest polyline in world units; short nodes are stubs
  maxLen: 40000,       // a sanity ceiling, not a working limit
  // World units per megabase. This is only the *fallback* used before a graph
  // has been calibrated (and when `autoLength` is off); `GraphModel.calibrate`
  // replaces it per graph. See `calibrateScale` for why a fixed value cannot
  // work: it makes the drawn size of a contig depend on nothing but its length,
  // so a 5 kb contig in a 5 Mb assembly lands on `minLen` and every contig in
  // the graph becomes the same stub.
  unitsPerMegabase: 1000,
  // Bandage's calibration targets: the *mean* node should come out this many
  // world units long, and the whole drawing at least this long in total.
  // (AssemblyGraph::determineGraphInfo, assemblygraph.cpp:408-420.)
  meanNodeLength: 40,
  minTotalGraphLength: 500,
  autoLength: true,    // false pins the scale to `unitsPerMegabase`
  lengthScale: 1,      // user multiplier from the "node length" slider
  // One polyline vertex per `particleSpacing` world units, with no practical
  // ceiling: this is Bandage's `nodeSegmentLength` (20.0). A long contig has to
  // be free to acquire many vertices, because that is the only way it can bend
  // into the sweeping curves that make a Bandage picture readable. Buckling is
  // held off by the layout's curvature term, not by starving the polyline.
  particleSpacing: 20,
  minParticles: 2,
  maxParticles: 400,
});

/**
 * Fewest vertices a contig that closes on itself is drawn with.
 *
 * A closed molecule has to read as a ring, and a ring needs enough vertices to
 * look like one: at two it is a line whatever the layout does, at six a
 * hexagon. Sixteen is round to the eye at any zoom a whole plasmid is viewed
 * at, and costs nothing -- these are single contigs, not the whole graph.
 */
export const RING_MIN_PARTICLES = 16;

/**
 * Vertices to spread around a closed molecule assembled in several contigs, so
 * the ring reads as a circle rather than as a polygon with a corner per contig.
 */
export const RING_VERTEX_BUDGET = 48;

/**
 * Node-drag deformation, Bandage's constants.
 *
 * `strength` and `falloffPower` reproduce `GraphicsItemNode::shiftPoints`
 * (graphicsitemnode.cpp:577-585) with `g_settings->dragStrength` (settings.cpp:66):
 * the grabbed vertex takes the whole movement and its neighbours a share that
 * decays as 2^(-d^1.8 / strength), so a short stub moves rigidly while a long
 * ribbon bends. `passes` is ours: Bandage lets the polyline stretch, but our
 * drawn length carries the contig's size, so the spacing is restored afterwards.
 */
export const DEFAULT_DRAG = Object.freeze({
  strength: 100,
  falloffPower: 1.8,
  passes: 6,
  anchorEnds: false,
});

/**
 * World units per base for a whole graph, Bandage's way.
 *
 * Bandage calibrates once per graph so the *mean* node is ~40 units long
 * whatever the assembly's size, which is what keeps the same picture legible
 * for a 40 kb phage and a 5 Mb chromosome. A fixed units-per-megabase cannot:
 * at 1000 u/Mb a 5 kb contig is 5 units long, which is the `minLen` floor, so
 * on any real bacterial assembly every contig collapses to an identical stub
 * and the drawing stops carrying length at all.
 *
 * @param {number} nSegments  segments that will be drawn
 * @param {number} totalBases total bases across those segments
 */
export function calibrateScale(nSegments, totalBases, geom = DEFAULT_GEOM) {
  const n = Number(nSegments) || 0;
  const total = Number(totalBases) || 0;
  const target = Math.max(n * (geom.meanNodeLength || 40), geom.minTotalGraphLength || 500);
  const megabases = total / 1e6;
  const perMegabase = megabases > 0 ? target / megabases : 10000;
  return perMegabase / 1e6;
}

/**
 * Draw length of a segment, **linear in sequence length**.
 *
 * Bandage draws a node's length in proportion to its bases, which is what makes
 * a long contig read as a long sweeping ribbon and a collapsed repeat as a
 * stub. A square-root mapping flattens that: a 32 kb contig and a 2.4 kb repeat
 * end up nearly the same size on screen and the picture stops carrying the
 * information people open a graph viewer to see.
 *
 * `perBase` comes from the graph (see `GraphModel.calibrate`) so the assembly
 * fills a consistent area whatever its size.
 */
export function drawLengthFor(bp, geom = DEFAULT_GEOM, perBase = null) {
  const L = Number(bp);
  const safe = Number.isFinite(L) && L > 0 ? L : 1;
  const scale = perBase === null ? perBaseScaleFor(geom) : perBase;
  return clamp(safe * scale * (geom.lengthScale || 1), geom.minLen, geom.maxLen);
}

/** World units per base. Absolute, so it does not depend on the graph. */
export function perBaseScaleFor(geom = DEFAULT_GEOM) {
  return (geom.unitsPerMegabase || 1000) / 1e6;
}

/** Particle count for a polyline, so long ribbons bend smoothly. */
export function particleCountFor(drawLen, geom = DEFAULT_GEOM) {
  const k = Math.round(drawLen / (geom.particleSpacing || 25));
  return clamp(k, geom.minParticles, geom.maxParticles);
}

/**
 * Which end of the *source* segment a link leaves from.
 * '+' -> END particle, '-' -> START particle.
 */
export function fromEndOf(fromOrient) {
  return fromOrient === '-' ? END_START : END_END;
}

/**
 * Which end of the *target* segment a link arrives at.
 * '+' -> START particle, '-' -> END particle.
 */
export function toEndOf(toOrient) {
  return toOrient === '-' ? END_END : END_START;
}

/**
 * Key identifying one directed use of a link: which end of what joins what.
 * A number rather than a string because it is built once per link and looked
 * up once per path step; the multiplier bounds it at ~1M segments, which is far
 * beyond anything the drawing accepts.
 */
export function linkKey(a, aEnd, b, bEnd) {
  return (a * 2 + aEnd) * 2147483 + (b * 2 + bEnd);
}

/** Particle index for a given end of a segment record. */
export function endParticle(seg, whichEnd) {
  return whichEnd === END_END ? seg.p0 + seg.k - 1 : seg.p0;
}

/**
 * The vertex one step *inward* from an end.
 *
 * An edge leaves a node along the node's own terminal direction, so the drawing
 * needs both the endpoint and its neighbour to know which way that is. For a
 * two-vertex segment the neighbour is the opposite end, which still gives the
 * right direction.
 */
export function innerParticle(seg, whichEnd) {
  if (seg.k < 2) return endParticle(seg, whichEnd);
  return whichEnd === END_END ? seg.p0 + seg.k - 2 : seg.p0 + 1;
}

/**
 * Restore the spacing of a polyline by Gauss-Seidel projection.
 *
 * The same maths the layout engine runs after every integration step, but
 * standalone so a drag can relax a chain on the main thread without waking the
 * worker. `fixed[i]` holds vertex `i` still; a chain with a fixed vertex
 * relaxes *towards* it, which is what makes a dragged node bend rather than
 * stretch.
 *
 * @param {Float32Array} px
 * @param {Float32Array} py
 * @param {number} p0    first particle of the chain
 * @param {number} k     particle count
 * @param {number} rest  wanted distance between consecutive particles
 * @param {Uint8Array|null} fixed  per-vertex "do not move" flags, or null
 * @param {number} passes
 */
export function relaxChain(px, py, p0, k, rest, fixed, passes) {
  if (k < 2 || !(rest > 0)) return;
  for (let pass = 0; pass < passes; pass++) {
    for (let i = 0; i + 1 < k; i++) {
      const a = p0 + i, b = a + 1;
      const fa = fixed ? fixed[i] : 0;
      const fb = fixed ? fixed[i + 1] : 0;
      if (fa && fb) continue;
      let dx = px[b] - px[a];
      let dy = py[b] - py[a];
      let d = Math.sqrt(dx * dx + dy * dy);
      if (d < 1e-4) { dx = rest; dy = 0; d = rest; }
      const corr = (d - rest) / d;
      if (!fa && !fb) {
        const hx = dx * corr * 0.5, hy = dy * corr * 0.5;
        px[a] += hx; py[a] += hy;
        px[b] -= hx; py[b] -= hy;
      } else if (fb) {
        px[a] += dx * corr; py[a] += dy * corr;
      } else {
        px[b] -= dx * corr; py[b] -= dy * corr;
      }
    }
  }
}

/** Format helpers shared by several panels. */
export function fmtBp(n) {
  const v = Number(n);
  if (!Number.isFinite(v)) return '—';
  if (Math.abs(v) >= 1e9) return (v / 1e9).toFixed(2) + ' Gb';
  if (Math.abs(v) >= 1e6) return (v / 1e6).toFixed(2) + ' Mb';
  if (Math.abs(v) >= 1e3) return (v / 1e3).toFixed(1) + ' kb';
  return v.toFixed(0) + ' bp';
}
export function fmtInt(n) {
  const v = Number(n);
  return Number.isFinite(v) ? v.toLocaleString('en-US', { maximumFractionDigits: 0 }) : '—';
}
export function fmtNum(n, dp = 2) {
  const v = Number(n);
  return Number.isFinite(v) ? v.toFixed(dp) : '—';
}

/* ------------------------------------------------------------- GraphModel */

export class GraphModel {
  constructor() {
    this.geom = { ...DEFAULT_GEOM };
    this.clear();
  }

  clear() {
    /** @type {Array<Object>} segment records */
    this.segments = [];
    /** @type {Map<string, Object>} */
    this.byName = new Map();
    /** @type {Array<Object>} resolved links */
    this.links = [];
    /**
     * Paths resolved against the drawn graph; see `_resolvePaths` for the
     * record shape. Raw payload order is preserved.
     * @type {Array<Object>}
     */
    this.paths = [];
    this.pathsByName = new Map();
    /** segment index -> path records that walk through it */
    this.pathsBySegment = new Map();
    this.droppedPathMembers = 0;
    /** link lookup keyed by both attachment points, for path resolution */
    this.linkAt = new Map();
    /** null, or the summary of the last contiguity search */
    this.contiguity = null;
    /** @type {Array<string>} reference sequence names seen in ref_hits */
    this.references = [];

    this.perBase = null;
    this.nParticles = 0;
    this.px = new Float32Array(0);
    this.py = new Float32Array(0);
    /** particle -> segment index */
    this.particleSeg = new Int32Array(0);
    /** 4 floats per segment: minX, minY, maxX, maxY */
    this.bbox = new Float32Array(0);

    /** @type {Array<{id:(number|string), segs:number[], length:number}>} */
    this.components = [];
    this.componentIndex = new Map(); // component id -> index into this.components

    /** adjacency: segIdx -> array of {seg, link, selfEnd, otherEnd} */
    this.adj = [];

    this.truncated = false;
    this.shown = 0;
    this.total = 0;

    this.stats = {
      depth: { min: 0, max: 0, has: false },
      gc: { min: 0, max: 1, has: false },
      length: { min: 1, max: 1 },
      totalLength: 0,
      hasSequences: false,
    };

    /** name -> best search hit (set by the search panel) */
    this.searchHits = new Map();
    this.droppedLinks = 0;
    this.revision = 0;
  }

  get isEmpty() { return this.segments.length === 0; }

  /**
   * Fix the bp -> world-unit scale for this graph.
   *
   * Called once per payload, before segments are built, so every drawn length
   * in the model shares one scale. `lengthScale` (the node-length slider) is
   * applied on top in `drawLengthFor`, so re-scaling never needs recalibration.
   */
  calibrate(nSegments, totalBases) {
    this.perBase = this.geom.autoLength === false
      ? perBaseScaleFor(this.geom)
      : calibrateScale(nSegments, totalBases, this.geom);
    return this.perBase;
  }

  /**
   * Re-map every segment's drawn length; used by the node-length slider.
   *
   * Particle counts are re-derived too: a segment that grows tenfold needs the
   * vertices to bend with, and one that shrinks should give them back. Callers
   * must re-seed positions after this, which `revision++` signals.
   */
  rescale(lengthScale) {
    this.geom.lengthScale = Number(lengthScale) || 1;

    // Keep the picture. The slider changes how long a contig is drawn, not
    // where it sits, so each polyline is resampled onto its new vertex count
    // along its own current shape and then stretched about its own centre.
    const oldPx = this.px;
    const oldPy = this.py;
    const old = this.segments.map(seg => ({ p0: seg.p0, k: seg.k, drawLen: seg.drawLen }));

    let p = 0;
    for (const seg of this.segments) {
      seg.drawLen = drawLengthFor(seg.length, this.geom, this.perBase);
      seg.k = seg.closed
        ? Math.max(particleCountFor(seg.drawLen, this.geom), RING_MIN_PARTICLES)
        : particleCountFor(seg.drawLen, this.geom);
      seg.p0 = p;
      p += seg.k;
    }
    this._reallocParticles(p);

    for (let s = 0; s < this.segments.length; s++) {
      this._resamplePolyline(this.segments[s], old[s], oldPx, oldPy);
    }

    for (const l of this.links) {
      const a = this.segments[l.a];
      const b = this.segments[l.b];
      l.pa = endParticle(a, l.aEnd);
      l.pb = endParticle(b, l.bEnd);
      l.paIn = innerParticle(a, l.aEnd);
      l.pbIn = innerParticle(b, l.bEnd);
    }
    this.revision++;
  }

  /**
   * Lay `seg`'s new vertices along the shape it had under `old`.
   *
   * Vertices are placed at equal arc-length fractions of the old polyline, then
   * the whole run is scaled about its midpoint by the ratio of the new drawn
   * length to the old one. A segment with no usable old geometry is left at the
   * origin for the layout to seed.
   */
  _resamplePolyline(seg, old, oldPx, oldPy) {
    const k = seg.k;
    if (!old || old.k < 1 || old.p0 + old.k > oldPx.length) return;

    // Arc length along the old polyline.
    const cum = new Float64Array(old.k);
    for (let i = 1; i < old.k; i++) {
      const dx = oldPx[old.p0 + i] - oldPx[old.p0 + i - 1];
      const dy = oldPy[old.p0 + i] - oldPy[old.p0 + i - 1];
      cum[i] = cum[i - 1] + Math.hypot(dx, dy);
    }
    const total = cum[old.k - 1];

    let cx = 0;
    let cy = 0;
    for (let i = 0; i < k; i++) {
      const t = k > 1 ? (i / (k - 1)) * total : 0;
      let x;
      let y;
      if (total <= 0) {
        x = oldPx[old.p0];
        y = oldPy[old.p0];
      } else {
        let j = 1;
        while (j < old.k - 1 && cum[j] < t) j++;
        const span = cum[j] - cum[j - 1];
        const f = span > 0 ? (t - cum[j - 1]) / span : 0;
        x = oldPx[old.p0 + j - 1] + (oldPx[old.p0 + j] - oldPx[old.p0 + j - 1]) * f;
        y = oldPy[old.p0 + j - 1] + (oldPy[old.p0 + j] - oldPy[old.p0 + j - 1]) * f;
      }
      this.px[seg.p0 + i] = x;
      this.py[seg.p0 + i] = y;
      cx += x;
      cy += y;
    }

    // Stretch about the centre so the node grows in place rather than dragging
    // one end across the drawing.
    const ratio = old.drawLen > 0 ? seg.drawLen / old.drawLen : 1;
    if (!(ratio > 0) || Math.abs(ratio - 1) < 1e-6 || k < 2) return;
    cx /= k;
    cy /= k;
    for (let i = 0; i < k; i++) {
      this.px[seg.p0 + i] = cx + (this.px[seg.p0 + i] - cx) * ratio;
      this.py[seg.p0 + i] = cy + (this.py[seg.p0 + i] - cy) * ratio;
    }
  }

  /** Resize the particle arrays and rebuild the particle -> segment map. */
  _reallocParticles(p) {
    this.nParticles = p;
    this.px = new Float32Array(p);
    this.py = new Float32Array(p);
    this.particleSeg = new Int32Array(p);
    for (const seg of this.segments) {
      for (let i = 0; i < seg.k; i++) this.particleSeg[seg.p0 + i] = seg.idx;
    }
  }

  /**
   * Rebuild the model from a `GET /api/graph` payload.
   * `keepPositions` re-uses the previous coordinates of segments with the same
   * name, so a re-fetch after an operation does not scramble the picture.
   */
  setData(payload, { keepPositions = true } = {}) {
    const prev = keepPositions ? this._snapshotPositions() : null;
    const geom = this.geom;

    const rawSegs = (payload && Array.isArray(payload.segments)) ? payload.segments : [];
    const rawLinks = (payload && Array.isArray(payload.links)) ? payload.links : [];

    // Calibrate before any segment is sized, so one scale covers the graph.
    // Duplicate names are skipped below; counting them here would shrink the
    // scale slightly, which is harmless, but the totals are cheap to get right.
    {
      const seen = new Set();
      let n = 0;
      let totalBases = 0;
      for (let i = 0; i < rawSegs.length; i++) {
        const s = rawSegs[i] || {};
        const name = String(s.name !== undefined && s.name !== null ? s.name : 'seg_' + i);
        if (seen.has(name)) continue;
        seen.add(name);
        n++;
        const L = Number(s.length);
        if (Number.isFinite(L) && L > 0) totalBases += L;
      }
      this.calibrate(n, totalBases);
    }

    this.segments = [];
    this.byName = new Map();
    let p = 0;

    // Which contigs close on themselves. Known before the segments are built
    // because a closed contig needs enough vertices to be drawn as a circle at
    // all -- two vertices can only ever be a straight line, however the layout
    // is seeded.
    const closedNames = new Set();
    for (let i = 0; i < rawLinks.length; i++) {
      const l = rawLinks[i] || {};
      if (l.from !== undefined && String(l.from) === String(l.to)) closedNames.add(String(l.from));
    }

    for (let i = 0; i < rawSegs.length; i++) {
      const s = rawSegs[i] || {};
      const name = String(s.name !== undefined && s.name !== null ? s.name : 'seg_' + i);
      if (this.byName.has(name)) continue; // duplicate names would break lookups
      const length = Number.isFinite(Number(s.length)) ? Number(s.length) : 0;
      const drawLen = drawLengthFor(length, geom, this.perBase);
      const closed = closedNames.has(name);
      const k = closed
        ? Math.max(particleCountFor(drawLen, geom), RING_MIN_PARTICLES)
        : particleCountFor(drawLen, geom);
      const depth = (s.depth === null || s.depth === undefined) ? null : Number(s.depth);
      const gc = (s.gc === null || s.gc === undefined) ? null : Number(s.gc);
      const rec = {
        idx: this.segments.length,
        name,
        length,
        depth: Number.isFinite(depth) ? depth : null,
        gc: Number.isFinite(gc) ? gc : null,
        component: (s.component === null || s.component === undefined) ? 0 : s.component,
        degStart: Number(s.deg_start) || 0,
        degEnd: Number(s.deg_end) || 0,
        circular: !!s.circular,
        closed,
        refHits: Array.isArray(s.ref_hits) ? s.ref_hits : [],
        // Filled in by `determineContiguity`; null while no search is active.
        contiguity: null,
        drawLen,
        k,
        p0: p,
      };
      rec.bestHit = bestRefHit(rec.refHits);
      p += k;
      this.byName.set(name, rec);
      this.segments.push(rec);
    }

    this.nParticles = p;
    this.px = new Float32Array(p);
    this.py = new Float32Array(p);
    this.particleSeg = new Int32Array(p);
    this.bbox = new Float32Array(this.segments.length * 4);
    for (const seg of this.segments) {
      for (let i = 0; i < seg.k; i++) this.particleSeg[seg.p0 + i] = seg.idx;
    }

    /* ---- links: resolve names to particles, drop dangling ones ---- */
    this.links = [];
    this.linkAt = new Map();
    this.droppedLinks = 0;
    for (let i = 0; i < rawLinks.length; i++) {
      const l = rawLinks[i] || {};
      const a = this.byName.get(String(l.from));
      const b = this.byName.get(String(l.to));
      if (!a || !b) { this.droppedLinks++; continue; } // truncated graph
      const fEnd = fromEndOf(l.from_orient);
      const tEnd = toEndOf(l.to_orient);
      this.links.push({
        idx: this.links.length,
        from: a.name, to: b.name,
        fromOrient: l.from_orient === '-' ? '-' : '+',
        toOrient: l.to_orient === '-' ? '-' : '+',
        overlap: Number(l.overlap) || 0,
        a: a.idx, b: b.idx,
        aEnd: fEnd, bEnd: tEnd,
        pa: endParticle(a, fEnd),
        pb: endParticle(b, tEnd),
        // One vertex inward from each attachment point. The renderer extends
        // these through the endpoints to get the Bézier control points, so an
        // edge continues the direction the node was already travelling.
        paIn: innerParticle(a, fEnd),
        pbIn: innerParticle(b, tEnd),
        // A link from a segment to itself is the graph's way of saying the
        // contig closes into a circle. Drawn as a straight chord it would run
        // through the node body and vanish under it, so the renderer bows it
        // out sideways and the layout leaves it out of the spring set.
        selfLoop: a.idx === b.idx,
      });
      // Both directions: a path step names the join, not which way it was
      // written, so either attachment order has to find the same link.
      const li = this.links.length - 1;
      this.linkAt.set(linkKey(a.idx, fEnd, b.idx, tEnd), li);
      this.linkAt.set(linkKey(b.idx, tEnd, a.idx, fEnd), li);
    }

    this._resolvePaths(payload && payload.paths);
    this.contiguity = null; // a new payload invalidates any search
    this.truncated = !!(payload && payload.truncated);
    this.shown = Number(payload && payload.shown) || this.segments.length;
    this.total = Number(payload && payload.total) || this.segments.length;

    // References: prefer the server-provided list, fall back to scanning hits.
    const refSet = new Set();
    if (payload && Array.isArray(payload.references)) {
      for (const r of payload.references) if (r) refSet.add(String(r));
    }
    for (const seg of this.segments) {
      for (const h of seg.refHits) if (h && h.ref) refSet.add(String(h.ref));
    }
    this.references = [...refSet].sort();

    this._buildAdjacency();
    this._buildComponents();
    this._refineRingResolution();
    this._computeStats();
    this._restorePositions(prev);
    this.updateBounds();
    this.revision++;
    return this;
  }

  /* ------------------------------------------------------------ internals */

  /**
   * Give contigs on a closed molecule enough vertices to follow its curve.
   *
   * A ring is drawn as the contigs around it, so with two vertices each they
   * are chords and the molecule comes out as a polygon: seven contigs make a
   * heptagon, which reads as a shape rather than as a circle. Enough vertices
   * between them and the ring is round.
   *
   * Only whole closed molecules qualify, so this costs nothing on an ordinary
   * graph, and it runs after components are known -- which is why it is a
   * second pass rather than part of building the segments.
   */
  _refineRingResolution() {
    let changed = false;
    for (const comp of this.components) {
      const walk = this._walkComponent(comp);
      if (!walk.ring || walk.selfLooped) continue;
      // Spread the ring's vertex budget over its contigs by length, so a long
      // contig curves and a stub is not given vertices it cannot use.
      let total = 0;
      for (const si of comp.segs) total += this.segments[si].drawLen;
      if (!(total > 0)) continue;
      for (const si of comp.segs) {
        const seg = this.segments[si];
        const share = Math.round(RING_VERTEX_BUDGET * (seg.drawLen / total));
        const want = Math.max(seg.k, Math.min(RING_MIN_PARTICLES, Math.max(3, share)));
        if (want !== seg.k) { seg.k = want; changed = true; }
      }
    }
    if (!changed) return;

    let p = 0;
    for (const seg of this.segments) { seg.p0 = p; p += seg.k; }
    this._reallocParticles(p);
    for (const l of this.links) {
      const a = this.segments[l.a];
      const b = this.segments[l.b];
      l.pa = endParticle(a, l.aEnd);
      l.pb = endParticle(b, l.bEnd);
      l.paIn = innerParticle(a, l.aEnd);
      l.pbIn = innerParticle(b, l.bEnd);
    }
  }

  _snapshotPositions() {
    if (!this.segments.length) return null;
    const m = new Map();
    for (const seg of this.segments) {
      const xs = new Float32Array(seg.k);
      const ys = new Float32Array(seg.k);
      for (let i = 0; i < seg.k; i++) { xs[i] = this.px[seg.p0 + i]; ys[i] = this.py[seg.p0 + i]; }
      m.set(seg.name, { xs, ys });
    }
    return m;
  }

  _restorePositions(prev) {
    // Seed everything first so nothing ends up at (0,0).
    this.seedPositions();
    if (!prev || !prev.size) return;
    // Re-use old coordinates where we can; new segments keep their seed but are
    // nudged next to a neighbour that already has a position, when possible.
    const restored = new Uint8Array(this.segments.length);
    for (const seg of this.segments) {
      const old = prev.get(seg.name);
      if (!old) continue;
      if (old.xs.length === seg.k) {
        for (let i = 0; i < seg.k; i++) { this.px[seg.p0 + i] = old.xs[i]; this.py[seg.p0 + i] = old.ys[i]; }
      } else {
        // Different particle count (size scale changed): resample the polyline.
        const n = old.xs.length;
        for (let i = 0; i < seg.k; i++) {
          const t = seg.k === 1 ? 0 : (i / (seg.k - 1)) * (n - 1);
          const j = Math.min(n - 1, Math.floor(t));
          const f = t - j;
          const j2 = Math.min(n - 1, j + 1);
          this.px[seg.p0 + i] = old.xs[j] * (1 - f) + old.xs[j2] * f;
          this.py[seg.p0 + i] = old.ys[j] * (1 - f) + old.ys[j2] * f;
        }
      }
      restored[seg.idx] = 1;
    }
    for (const seg of this.segments) {
      if (restored[seg.idx]) continue;
      const nb = (this.adj[seg.idx] || []).find((e) => restored[e.seg]);
      if (!nb) continue;
      const other = this.segments[nb.seg];
      const ox = this.px[endParticle(other, nb.otherEnd)];
      const oy = this.py[endParticle(other, nb.otherEnd)];
      const ang = Math.random() * Math.PI * 2;
      this._placeSegment(seg, ox + Math.cos(ang) * 24, oy + Math.sin(ang) * 24, ang);
    }
  }

  /**
   * Resolve `P`/`W` line walks against the segments that are actually drawn.
   *
   * A truncated view can leave a path pointing at segments nobody drew, so
   * missing members are dropped and counted rather than silently skewing the
   * walk: a path that lost a member is no longer a contiguous walk, and the
   * join across the hole is flagged the same way a scaffold gap is.
   *
   * Each record is:
   *   {idx, name, segs: Int32Array, orients: string, links: Int32Array,
   *    gapAfter: Uint8Array, segSet: Set, linkSet: Set,
   *    dropped, total, broken, bases}
   * where `links[j]` joins `segs[j]` to `segs[j+1]` (-1 when the join is a
   * scaffold gap, crosses a dropped member, or has no link in the graph) and
   * `gapAfter[j]` marks the joins that are not graph edges at all.
   */
  _resolvePaths(rawPaths) {
    this.paths = [];
    this.pathsByName = new Map();
    this.pathsBySegment = new Map();
    this.droppedPathMembers = 0;
    if (!Array.isArray(rawPaths)) return;

    for (const raw of rawPaths) {
      if (!raw) continue;
      const rawSteps = Array.isArray(raw.steps) ? raw.steps : [];
      const gapSet = new Set((Array.isArray(raw.gaps) ? raw.gaps : []).map(Number));
      const segs = [];
      const orients = [];
      const gapAfter = [];
      let dropped = 0;
      let bases = 0;
      // True while the walk has been interrupted since the last kept member --
      // by a scaffold gap in the file, or by a member we are not drawing.
      let broken = false;
      for (let i = 0; i < rawSteps.length; i++) {
        const step = String(rawSteps[i]);
        const last = step.charAt(step.length - 1);
        const signed = last === '+' || last === '-';
        const seg = this.byName.get(signed ? step.slice(0, -1) : step);
        if (!seg) { dropped++; broken = true; continue; }
        if (segs.length) gapAfter[segs.length - 1] = broken ? 1 : 0;
        segs.push(seg.idx);
        orients.push(last === '-' ? '-' : '+');
        bases += seg.length;
        broken = gapSet.has(i);
      }
      if (segs.length) gapAfter[segs.length - 1] = 0;

      const links = new Int32Array(Math.max(0, segs.length - 1)).fill(-1);
      let brokenJoins = 0;
      for (let j = 0; j + 1 < segs.length; j++) {
        if (gapAfter[j]) { brokenJoins++; continue; }
        const key = linkKey(segs[j], fromEndOf(orients[j]),
                            segs[j + 1], toEndOf(orients[j + 1]));
        const li = this.linkAt.get(key);
        if (li === undefined) { brokenJoins++; continue; }
        links[j] = li;
      }

      const rec = {
        idx: this.paths.length,
        name: String(raw.name === undefined || raw.name === null
          ? 'path_' + this.paths.length : raw.name),
        segs: Int32Array.from(segs),
        orients: orients.join(''),
        links,
        gapAfter: Uint8Array.from(gapAfter),
        segSet: new Set(segs),
        linkSet: new Set([...links].filter((l) => l >= 0)),
        dropped,
        total: rawSteps.length,
        broken: brokenJoins,
        bases,
      };
      this.droppedPathMembers += dropped;
      this.paths.push(rec);
      this.pathsByName.set(rec.name, rec);
      for (const si of rec.segSet) {
        let list = this.pathsBySegment.get(si);
        if (!list) { list = []; this.pathsBySegment.set(si, list); }
        list.push(rec);
      }
    }
  }

  /** Resolved path by name, or null. */
  pathByName(name) { return this.pathsByName.get(String(name)) || null; }

  /** Every resolved path that walks through a segment. */
  pathsForSegment(segIdx) { return this.pathsBySegment.get(segIdx) || []; }

  _buildAdjacency() {
    this.adj = this.segments.map(() => []);
    for (const l of this.links) {
      this.adj[l.a].push({ seg: l.b, link: l.idx, selfEnd: l.aEnd, otherEnd: l.bEnd });
      if (l.a !== l.b) {
        this.adj[l.b].push({ seg: l.a, link: l.idx, selfEnd: l.bEnd, otherEnd: l.aEnd });
      }
    }
  }

  /**
   * Group segments into components. The server supplies a `component` field;
   * when it is missing or constant we fall back to a union over the links so
   * component colouring and per-component layout still work.
   */
  _buildComponents() {
    const n = this.segments.length;
    const groups = new Map();
    // Only fall back to union-find when the server genuinely did not supply the
    // field. Treating "every id happens to be 0" as missing renumbered the ids
    // of a single-component view to local indices, so filtering to a component
    // and then filtering again asked the server for a component number it had
    // never issued, and the canvas went blank with no error.
    let usable = false;
    for (const seg of this.segments) {
      if (seg.component !== null && seg.component !== undefined) { usable = true; break; }
    }
    if (!usable && this.links.length) {
      // Union-find over links.
      const parent = new Int32Array(n);
      for (let i = 0; i < n; i++) parent[i] = i;
      const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
      for (const l of this.links) {
        const ra = find(l.a), rb = find(l.b);
        if (ra !== rb) parent[ra] = rb;
      }
      for (const seg of this.segments) seg.component = find(seg.idx);
    }
    for (const seg of this.segments) {
      const key = seg.component;
      let g = groups.get(key);
      if (!g) { g = { id: key, segs: [], length: 0 }; groups.set(key, g); }
      g.segs.push(seg.idx);
      g.length += seg.length;
    }
    this.components = [...groups.values()].sort((a, b) => b.length - a.length || b.segs.length - a.segs.length);
    this.componentIndex = new Map();
    this.components.forEach((c, i) => {
      this.componentIndex.set(c.id, i);
      for (const s of c.segs) this.segments[s].compIndex = i;
    });
  }

  _computeStats() {
    let dMin = Infinity, dMax = -Infinity, hasD = false;
    let gMin = Infinity, gMax = -Infinity, hasG = false;
    let lMin = Infinity, lMax = 1, total = 0;
    for (const s of this.segments) {
      if (s.depth !== null && Number.isFinite(s.depth)) {
        hasD = true;
        if (s.depth < dMin) dMin = s.depth;
        if (s.depth > dMax) dMax = s.depth;
      }
      if (s.gc !== null && Number.isFinite(s.gc)) {
        hasG = true;
        if (s.gc < gMin) gMin = s.gc;
        if (s.gc > gMax) gMax = s.gc;
      }
      if (s.length < lMin) lMin = s.length;
      if (s.length > lMax) lMax = s.length;
      total += s.length;
    }
    if (!hasD) { dMin = 0; dMax = 1; }
    if (!hasG) { gMin = 0; gMax = 1; }
    if (!Number.isFinite(lMin)) lMin = 1;
    if (dMax <= dMin) dMax = dMin + 1;
    if (gMax <= gMin) { gMin = Math.max(0, gMin - 0.01); gMax = gMin + 0.02; }
    this.stats = {
      depth: { min: dMin, max: dMax, has: hasD },
      gc: { min: gMin, max: gMax, has: hasG },
      length: { min: Math.max(1, lMin), max: Math.max(1, lMax) },
      totalLength: total,
      hasSequences: this.stats ? this.stats.hasSequences : false,
    };
  }

  /* ------------------------------------------------------------ geometry */

  /** Lay one segment out as a straight polyline centred on (cx, cy). */
  _placeSegment(seg, cx, cy, angle) {
    const half = seg.drawLen / 2;
    const dx = Math.cos(angle), dy = Math.sin(angle);
    const x0 = cx - dx * half, y0 = cy - dy * half;
    const step = seg.k > 1 ? seg.drawLen / (seg.k - 1) : 0;
    for (let i = 0; i < seg.k; i++) {
      this.px[seg.p0 + i] = x0 + dx * step * i;
      this.py[seg.p0 + i] = y0 + dy * step * i;
    }
  }

  /**
   * Walk a component, following its links, and say whether it closes.
   *
   * The order is a depth-first walk from a dead end where there is one, so a
   * contig's neighbours in the walk are its neighbours in the graph. That is
   * what lets the seed lay the component out along its own structure instead of
   * scattering it and asking the force model to discover the structure again.
   *
   * `cyclic` counts links rather than inspecting the walk: a connected
   * component with at least as many links as contigs contains a cycle, and a
   * tree has exactly one fewer. It is the cheap, exact test.
   */
  _walkComponent(comp) {
    const segs = comp.segs;
    const inComp = new Set(segs);
    const links = new Set();
    for (const si of segs) {
      for (const e of this.adj[si] || []) {
        if (inComp.has(e.seg)) links.add(e.link);
      }
    }

    // Start at a dead end if the component has one, otherwise at the longest
    // contig, so a linear component is walked end to end rather than from its
    // middle.
    let root = segs[0];
    let bestDeg = Infinity;
    for (const si of segs) {
      const deg = (this.adj[si] || []).length;
      const better = deg < bestDeg
        || (deg === bestDeg && this.segments[si].drawLen > this.segments[root].drawLen);
      if (better) { root = si; bestDeg = deg; }
    }

    const order = [];
    const seen = new Set([root]);
    const stack = [root];
    while (stack.length) {
      const u = stack.pop();
      order.push(u);
      // Longest neighbour last, so it is popped first and the walk follows the
      // backbone rather than wandering off down a short bubble arm.
      const next = (this.adj[u] || [])
        .filter((e) => inComp.has(e.seg) && !seen.has(e.seg))
        .sort((a, b) => this.segments[a.seg].drawLen - this.segments[b.seg].drawLen);
      for (const e of next) { seen.add(e.seg); stack.push(e.seg); }
    }
    for (const si of segs) if (!seen.has(si)) order.push(si);

    // Independent cycles through the component: 1 is a plain loop, 0 a tree,
    // and a large number means a tangle whose shape is not a ring at all.
    const excess = links.size - segs.length + 1;

    // A ring seed is only right when the component really is a closed molecule:
    // a single contig that links to itself, or contigs joined nose to tail with
    // nothing else attached. Force a ring on a tangle and every extra link
    // becomes a chord pulling the circle shut, which the layout then has to
    // undo -- measurably worse than not seeding a ring at all.
    const selfLooped = segs.length === 1
      && (this.adj[segs[0]] || []).some((e) => e.seg === segs[0]);
    const closedWalk = excess === 1
      && segs.every((si) => (this.adj[si] || []).filter((e) => inComp.has(e.seg)).length === 2);

    return { order, ring: selfLooped || closedWalk, selfLooped, excess };
  }

  /**
   * Give every particle a starting position that already has the component's
   * shape: a component containing a cycle is laid on a ring, one without on a
   * line, each contig following the curve rather than sitting as a chord.
   *
   * Seeding this way is what makes a circular replicon look circular. A ring is
   * a stable equilibrium of the force model, so once the contigs start on one
   * they stay on it; started from a scatter, the same component settles into a
   * knot that happens to be connected, and nothing about the picture says the
   * molecule closes. It also converges much faster, because the layout is
   * refining a shape rather than discovering it.
   */
  seedPositions() {
    const comps = this.components;
    if (!comps.length) return;
    const cols = Math.max(1, Math.ceil(Math.sqrt(comps.length)));

    // The gap between consecutive contigs on the curve. Links rest at a
    // fraction of a polyline step, so leaving room for one here means the seed
    // already satisfies them and the force model has nothing to undo.
    const gap = (this.geom.particleSpacing || 20) * 0.25;

    // Everything here is in drawn units, and the bp -> unit scale is calibrated
    // per graph, so a contig can be five units long or five thousand. The grid
    // cell has to come from the shape each component will actually be seeded
    // into, or components start overlapping and the force model spends its
    // budget pulling apart a tangle it was handed.
    const plans = comps.map((comp) => {
      const walk = this._walkComponent(comp);
      let perimeter = 0;
      for (const si of walk.order) perimeter += this.segments[si].drawLen + gap;
      // A closed molecule keeps no gap after its last contig: the walk has to
      // meet its own start, and a spare gap would leave the ring visibly open.
      if (walk.ring) perimeter -= walk.selfLooped ? gap : 0;
      // A ring of this perimeter is that wide across; a line is its own length.
      const extent = walk.ring ? perimeter / Math.PI : perimeter;
      return { ...walk, perimeter, extent };
    });

    let cell = 0;
    for (const p of plans) cell = Math.max(cell, p.extent * 1.25 + 40);

    plans.forEach((plan, ci) => {
      if (!(plan.perimeter > 0)) return;
      const gx = (ci % cols) * cell;
      const gy = Math.floor(ci / cols) * cell;

      if (plan.ring) {
        // Radius chosen so the contigs laid end to end exactly close the ring.
        const R = plan.perimeter / (2 * Math.PI);
        let s = 0;
        for (const si of plan.order) {
          const seg = this.segments[si];
          this._placeSegmentOnRing(seg, gx, gy, R, s, plan.perimeter);
          s += seg.drawLen + gap;
        }
      } else {
        let s = -plan.perimeter / 2;
        for (const si of plan.order) {
          const seg = this.segments[si];
          this._placeSegment(seg, gx + s + seg.drawLen / 2, gy, 0);
          s += seg.drawLen + gap;
        }
      }
    });
    this.updateBounds();
  }

  /**
   * Lay one contig along an arc of a ring, vertex by vertex.
   *
   * Placing the contig as a straight chord would leave a polygon whose corners
   * the force model then has to round off; following the arc means the ring is
   * already smooth, which is the point of seeding it as a ring at all.
   */
  _placeSegmentOnRing(seg, cx, cy, R, sStart, perimeter) {
    // A closed contig's own last vertex is the neighbour of its first, so its
    // vertices tile k gaps around the ring rather than k - 1 along a strand.
    const spans = seg.closed ? seg.k : Math.max(1, seg.k - 1);
    const step = seg.drawLen / spans;
    for (let i = 0; i < seg.k; i++) {
      const a = ((sStart + step * i) / perimeter) * Math.PI * 2;
      this.px[seg.p0 + i] = cx + Math.cos(a) * R;
      this.py[seg.p0 + i] = cy + Math.sin(a) * R;
    }
  }

  /** Recompute per-segment bounding boxes (used for culling and hit tests). */
  updateBounds() {
    for (const seg of this.segments) this._updateSegBounds(seg);
  }

  _updateSegBounds(seg) {
    const { px, py, bbox } = this;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (let i = 0; i < seg.k; i++) {
      const x = px[seg.p0 + i], y = py[seg.p0 + i];
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
    if (!Number.isFinite(minX)) { minX = maxX = minY = maxY = 0; }
    const o = seg.idx * 4;
    bbox[o] = minX; bbox[o + 1] = minY; bbox[o + 2] = maxX; bbox[o + 3] = maxY;
  }

  /** World bounding box of some (or all) segments: [minX, minY, maxX, maxY]. */
  boundsOf(indices) {
    const list = indices && indices.length ? indices : this.segments.map((s) => s.idx);
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const i of list) {
      const o = i * 4;
      if (this.bbox[o] < minX) minX = this.bbox[o];
      if (this.bbox[o + 1] < minY) minY = this.bbox[o + 1];
      if (this.bbox[o + 2] > maxX) maxX = this.bbox[o + 2];
      if (this.bbox[o + 3] > maxY) maxY = this.bbox[o + 3];
    }
    if (!Number.isFinite(minX)) return [-100, -100, 100, 100];
    if (maxX - minX < 1) { minX -= 50; maxX += 50; }
    if (maxY - minY < 1) { minY -= 50; maxY += 50; }
    return [minX, minY, maxX, maxY];
  }

  /** Translate whole segments (used by node dragging). */
  translateSegments(indices, dx, dy) {
    for (const i of indices) {
      const seg = this.segments[i];
      if (!seg) continue;
      for (let j = 0; j < seg.k; j++) { this.px[seg.p0 + j] += dx; this.py[seg.p0 + j] += dy; }
      const o = i * 4;
      this.bbox[o] += dx; this.bbox[o + 1] += dy; this.bbox[o + 2] += dx; this.bbox[o + 3] += dy;
    }
  }

  /* --------------------------------------------------------- contiguity */

  /**
   * Bandage's doubled graph, as two adjacency lists over *oriented* nodes.
   *
   * Bandage stores each segment twice -- once per strand -- and every link
   * twice, as an edge and its reverse complement. Contiguity is defined on that
   * doubled graph, so the search cannot run on our bidirected adjacency
   * directly. An oriented node is `segment * 2 + reversed`, which makes the
   * reverse complement `node ^ 1`; a link joining (a, aEnd) to (b, bEnd) is an
   * edge leaving the node that exits through `aEnd` and entering the node that
   * enters through `bEnd`.
   */
  _orientedAdjacency() {
    const n = this.segments.length * 2;
    const out = Array.from({ length: n }, () => new Set());
    const inn = Array.from({ length: n }, () => new Set());
    for (const l of this.links) {
      const u = l.a * 2 + (l.aEnd === END_END ? 0 : 1);
      const v = l.b * 2 + (l.bEnd === END_START ? 0 : 1);
      out[u].add(v); inn[v].add(u);
      out[v ^ 1].add(u ^ 1); inn[u ^ 1].add(v ^ 1);
    }
    return {
      out: out.map((s) => [...s]),
      inn: inn.map((s) => [...s]),
    };
  }

  /** Accept names, indices or segment records for a starting set. */
  _resolveSegmentList(list) {
    const out = [];
    const seen = new Set();
    for (const item of (Array.isArray(list) ? list : [list])) {
      let seg = null;
      if (item === null || item === undefined) continue;
      if (typeof item === 'object') seg = this.segments[item.idx] || null;
      else if (typeof item === 'number') seg = this.segments[item] || null;
      else seg = this.byName.get(String(item)) || null;
      if (seg && !seen.has(seg.idx)) { seen.add(seg.idx); out.push(seg.idx); }
    }
    return out;
  }

  /**
   * Bandage's "determine contiguity", segment by segment.
   *
   * Faithful to `DeBruijnNode::determineContiguity` (debruijnnode.cpp:163) and
   * the two traversals it leans on in debruijnedge.cpp. For each edge leaving a
   * starting node -- in *both* directions, since Bandage walks a node's
   * incoming edges backwards as well as its outgoing ones forwards -- it traces
   * every path up to `steps` nodes long. A node on *any* of those paths is only
   * MAYBE_CONTIGUOUS; a node on *all* of them is CONTIGUOUS, because then the
   * walk cannot avoid it. Every node so touched then gets the converse test:
   * if all of *its* onward paths lead back to the start, it is CONTIGUOUS too.
   *
   * Bandage does the intersection twice, once ignoring strand and once folding
   * each path's reverse complements in, which catches the case where two routes
   * pass through the same segment on opposite strands. We keep that, then
   * collapse to the four classes the drawing distinguishes by taking each
   * segment's better strand -- exactly what `GraphicsItemNode::setColour` does
   * for a node drawn without an arrow (graphicsitemnode.cpp:381-387).
   *
   * @param {Array} starts    segment names, indices or records
   * @param {object} [opts]
   * @param {number} [opts.steps=15]  Bandage's `contiguitySearchSteps`
   * @param {number} [opts.budget]    node visits before the search gives up
   * @returns {object|null} summary, also left on `this.contiguity`
   */
  determineContiguity(starts, opts = {}) {
    const startIdx = this._resolveSegmentList(starts);
    if (!startIdx.length) { this.clearContiguity(); return null; }

    const steps = Math.max(1, Math.round(Number(opts.steps) || 15));
    // The path trace is exponential in the branching factor, which is fine for
    // the tangles Bandage was written for and ruinous for a repeat-rich one. A
    // visit budget keeps the browser responsive; the result says when it bit.
    const budget = Math.max(1000, Math.round(Number(opts.budget) || 2e6));
    const { out, inn } = this._orientedAdjacency();
    const status = new Uint8Array(this.segments.length * 2).fill(ST_NOT);
    let spent = 0;
    let truncated = false;

    const upgrade = (nd, s) => { if (s < status[nd]) status[nd] = s; };
    const timesIn = (path, nd) => {
      let c = 0;
      for (let i = 0; i < path.length; i++) if (path[i] === nd) c++;
      return c;
    };

    // DeBruijnEdge::tracePaths (debruijnedge.cpp:143). `node` is the node the
    // edge being followed leads to; the starting node is never in a path.
    const trace = (node, forward, left, soFar, startOn, allPaths) => {
      if (++spent > budget) { truncated = true; return; }
      const path = soFar.concat(node);
      if (--left === 0) { allPaths.push(path); return; }
      const nexts = forward ? out[node] : inn[node];
      if (!nexts.length) { allPaths.push(path); return; }
      for (const nn of nexts) {
        // Back at the start: the walk has closed a loop, so the path is done.
        if (nn === startOn) { allPaths.push(path); continue; }
        // Twice already means we are going round a cycle, not exploring.
        if (timesIn(path, nn) < 2) trace(nn, forward, left, path, startOn, allPaths);
        if (truncated) return;
      }
    };

    // DeBruijnNode::getNodesCommonToAllPaths (debruijnnode.cpp:241).
    const commonToAll = (allPaths, withRC) => {
      if (!allPaths.length) return [];
      let common = new Set(allPaths[0]);
      for (let i = 1; i < allPaths.length && common.size; i++) {
        const other = new Set();
        for (const nd of allPaths[i]) { other.add(nd); if (withRC) other.add(nd ^ 1); }
        const next = new Set();
        for (const nd of common) if (other.has(nd)) next.add(nd);
        common = next;
      }
      return common;
    };

    // DeBruijnEdge::leadsOnlyToNode (debruijnedge.cpp:228).
    const leadsOnlyTo = (node, forward, left, target, soFar, withRC) => {
      if (++spent > budget) { truncated = true; return false; }
      const path = soFar.concat(node);
      // Back where this check started: the walk could be circular DNA that
      // never reaches the target, so it does not count as leading there.
      if (node === path[0]) return false;
      if (node === target) return true;
      if (withRC && (node ^ 1) === target) return true;
      if (--left === 0) return false;
      const nexts = forward ? out[node] : inn[node];
      if (!nexts.length) return false;
      for (const nn of nexts) {
        if (timesIn(path, nn) < 2
            && !leadsOnlyTo(nn, forward, left, target, path, withRC)) return false;
      }
      return true;
    };

    // DeBruijnNode::doesPathLeadOnlyToNode (debruijnnode.cpp:294).
    const leadsBack = (node, target, withRC) => {
      for (const nn of out[node]) {
        if (leadsOnlyTo(nn, true, steps, target, [node], withRC)) return true;
      }
      for (const nn of inn[node]) {
        if (leadsOnlyTo(nn, false, steps, target, [node], withRC)) return true;
      }
      return false;
    };

    for (const si of startIdx) {
      const startOn = si * 2;
      upgrade(startOn, ST_STARTING);
      const checked = new Set();
      const edges = out[startOn].map((nd) => [nd, true])
        .concat(inn[startOn].map((nd) => [nd, false]));

      for (const [first, forward] of edges) {
        const allPaths = [];
        const wasTruncated = truncated;
        trace(first, forward, steps, [], startOn, allPaths);
        for (const path of allPaths) {
          for (const nd of path) { upgrade(nd, ST_MAYBE); checked.add(nd); }
        }
        // An incomplete path set makes the intersection too generous: a path
        // that was never traced cannot vote a node out. Rather than claim a
        // contiguity the graph has not been shown to have, leave these nodes
        // on the MAYBE_CONTIGUOUS they already earned.
        if (truncated && !wasTruncated) continue;
        for (const nd of commonToAll(allPaths, false)) upgrade(nd, ST_CONTIG_STRAND);
        for (const nd of commonToAll(allPaths, true)) {
          upgrade(nd, ST_CONTIG_EITHER);
          upgrade(nd ^ 1, ST_CONTIG_EITHER);
        }
      }

      for (const nd of checked) {
        if (truncated) break;
        const was = status[nd];
        if (was !== ST_CONTIG_STRAND && leadsBack(nd, startOn, false)) {
          upgrade(nd, ST_CONTIG_STRAND);
        }
        if (was !== ST_CONTIG_STRAND && was !== ST_CONTIG_EITHER
            && leadsBack(nd, startOn, true)) {
          upgrade(nd, ST_CONTIG_EITHER);
          upgrade(nd ^ 1, ST_CONTIG_EITHER);
        }
      }
    }

    const counts = [0, 0, 0, 0];
    for (const seg of this.segments) {
      const s = Math.min(status[seg.idx * 2], status[seg.idx * 2 + 1]);
      seg.contiguity = s === ST_STARTING ? CONTIGUITY.STARTING
        : (s === ST_CONTIG_STRAND || s === ST_CONTIG_EITHER) ? CONTIGUITY.CONTIGUOUS
          : s === ST_MAYBE ? CONTIGUITY.MAYBE_CONTIGUOUS : CONTIGUITY.NOT_CONTIGUOUS;
      counts[seg.contiguity]++;
    }

    this.contiguity = {
      starts: startIdx.map((i) => this.segments[i].name),
      steps,
      truncated,
      visits: spent,
      counts: {
        starting: counts[CONTIGUITY.STARTING],
        contiguous: counts[CONTIGUITY.CONTIGUOUS],
        maybe: counts[CONTIGUITY.MAYBE_CONTIGUOUS],
        not: counts[CONTIGUITY.NOT_CONTIGUOUS],
      },
    };
    return this.contiguity;
  }

  /** Forget a contiguity search; every segment goes back to `contiguity: null`. */
  clearContiguity() {
    for (const seg of this.segments) seg.contiguity = null;
    this.contiguity = null;
  }

  /* ------------------------------------------------------- node dragging */

  /**
   * Index (within the segment) of the polyline vertex nearest a world point.
   * Bandage's `m_grabIndex` (graphicsitemnode.cpp:471).
   */
  nearestVertex(segIdx, wx, wy) {
    const seg = this.segments[segIdx];
    if (!seg) return -1;
    let best = 0;
    let bestD = Infinity;
    for (let i = 0; i < seg.k; i++) {
      const dx = this.px[seg.p0 + i] - wx;
      const dy = this.py[seg.p0 + i] - wy;
      const d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  }

  /**
   * Bend a segment by dragging one of its vertices.
   *
   * The grabbed vertex takes the whole displacement and the rest of the chain
   * a share that falls away with index distance (Bandage's `shiftPoints`), then
   * the polyline's own spacing is projected back with the grabbed vertex
   * pinned, so the ribbon keeps the drawn length that carries the contig's
   * size. The two end vertices move like any other, which is what keeps the
   * neighbouring segments' edges attached: the renderer draws them from the end
   * particles, so they follow wherever the ends land. Pass `anchorEnds` to hold
   * the ends still instead and bend only the middle.
   *
   * @param {number} segIdx
   * @param {number} vertex  index within the segment, 0..k-1
   * @param {number} dx
   * @param {number} dy
   * @param {object} [opts]  overrides for DEFAULT_DRAG
   */
  deformSegment(segIdx, vertex, dx, dy, opts = {}) {
    const seg = this.segments[segIdx];
    if (!seg || !(dx || dy)) return;
    const v = clamp(vertex | 0, 0, seg.k - 1);
    const o = { ...DEFAULT_DRAG, ...opts };
    const strength = o.strength > 0 ? o.strength : DEFAULT_DRAG.strength;
    const anchored = !!o.anchorEnds && seg.k > 2;

    for (let i = 0; i < seg.k; i++) {
      if (i !== v && anchored && (i === 0 || i === seg.k - 1)) continue;
      const w = i === v ? 1
        : Math.pow(2, -Math.pow(Math.abs(i - v), o.falloffPower) / strength);
      this.px[seg.p0 + i] += dx * w;
      this.py[seg.p0 + i] += dy * w;
    }

    if (seg.k > 1 && o.passes > 0) {
      const fixed = new Uint8Array(seg.k);
      if (anchored) { fixed[0] = 1; fixed[seg.k - 1] = 1; }
      fixed[v] = 1;
      relaxChain(this.px, this.py, seg.p0, seg.k,
                 seg.drawLen / (seg.k - 1), fixed, o.passes | 0);
    }
    this._updateSegBounds(seg);
  }

  /**
   * A stateful handle for one drag gesture, for pointer code that only has a
   * cursor position to give:
   *
   *   const drag = graph.beginVertexDrag(si, wx, wy);
   *   drag.moveTo(wx, wy);   // on every pointermove
   *   drag.end();            // on pointerup
   *
   * `drag.particle` is the absolute particle index of the grabbed vertex, which
   * is what `LayoutController.start({params: {pinned: [...]}})` wants if the
   * force layout is to keep relaxing the rest of the graph around it.
   */
  beginVertexDrag(segIdx, wx, wy, opts = {}) {
    const seg = this.segments[segIdx];
    if (!seg) return null;
    const vertex = this.nearestVertex(segIdx, wx, wy);
    const self = this;
    return {
      seg: segIdx,
      vertex,
      particle: seg.p0 + vertex,
      moveBy(dx, dy) { self.deformSegment(segIdx, vertex, dx, dy, opts); },
      moveTo(x, y) {
        self.deformSegment(segIdx, vertex,
          x - self.px[seg.p0 + vertex], y - self.py[seg.p0 + vertex], opts);
      },
      end() { self._updateSegBounds(seg); },
    };
  }

  /* ------------------------------------------------------------- queries */

  segmentByName(name) { return this.byName.get(String(name)) || null; }

  /** Indices of every segment in a component, addressed by component *id*. */
  segmentsInComponent(compId) {
    const i = this.componentIndex.get(compId);
    if (i === undefined) return [];
    return this.components[i].segs.slice();
  }

  degree(seg) { return (seg.degStart || 0) + (seg.degEnd || 0); }

  /** Human-readable neighbour list for the selection panel. */
  neighbourNames(seg) {
    const out = [];
    for (const e of (this.adj[seg.idx] || [])) out.push(this.segments[e.seg].name);
    return [...new Set(out)];
  }

  /** Apply search hits (`POST /api/search`), keeping the best per segment. */
  setSearchHits(hits) {
    this.searchHits = new Map();
    if (!Array.isArray(hits)) return;
    for (const h of hits) {
      if (!h || h.segment === undefined || h.segment === null) continue;
      const key = String(h.segment);
      const cur = this.searchHits.get(key);
      const score = Number(h.bitscore);
      const curScore = cur ? Number(cur.bitscore) : -Infinity;
      if (!cur || (Number.isFinite(score) && score > curScore) ||
          (!Number.isFinite(curScore) && Number(h.identity) > Number(cur.identity))) {
        this.searchHits.set(key, h);
      }
    }
  }

  clearSearchHits() { this.searchHits = new Map(); }

  /**
   * Pack the structural data the layout worker needs. Everything is a typed
   * array so it can be transferred rather than cloned.
   */
  /**
   * Which components are closed molecules, and how big a ring each should be.
   *
   * A circular molecule carries no information about its own 2D shape: every
   * way of drawing the loop is equally faithful, so the layout is free to pick
   * the one a reader can actually recognise, and that is a circle. The engine
   * uses these to hold ring components open; a chain of contigs joined nose to
   * tail is otherwise neutrally stable and crumples for free.
   */
  ringComponents() {
    // Keyed by segment, not by component. The engine builds its own dense
    // component numbering from `segComp`, which is in segment order and need
    // not agree with the order of `this.components`; indexing these by the
    // latter silently applies one component's ring to another.
    const segRing = new Uint8Array(this.segments.length);
    const segRingRadius = new Float32Array(this.segments.length);
    const gap = (this.geom.particleSpacing || 20) * 0.25;
    for (const comp of this.components) {
      const walk = this._walkComponent(comp);
      if (!walk.ring) continue;
      let perimeter = 0;
      for (const si of comp.segs) perimeter += this.segments[si].drawLen + gap;
      if (walk.selfLooped) perimeter -= gap;
      if (!(perimeter > 0)) continue;
      const radius = perimeter / (2 * Math.PI);
      for (const si of comp.segs) { segRing[si] = 1; segRingRadius[si] = radius; }
    }
    return { segRing, segRingRadius };
  }

  toLayoutArrays() {
    const n = this.segments.length;
    const segP0 = new Int32Array(n);
    const segK = new Int32Array(n);
    const segRest = new Float32Array(n);
    const segLen = new Float32Array(n);
    const segComp = new Int32Array(n);
    // A contig that links to itself is a closed molecule. The layout has to
    // know, because its polyline is a loop rather than a strand: the last
    // vertex is a neighbour of the first, and nothing else would hold the ring
    // shut once the self-link is left out of the spring set.
    const segClosed = new Uint8Array(n);
    for (const s of this.segments) {
      segP0[s.idx] = s.p0;
      segK[s.idx] = s.k;
      // A closed polyline has k gaps between its k vertices, not k - 1.
      segRest[s.idx] = s.closed
        ? s.drawLen / Math.max(1, s.k)
        : (s.k > 1 ? s.drawLen / (s.k - 1) : s.drawLen);
      segLen[s.idx] = s.drawLen;
      segComp[s.idx] = s.compIndex === undefined ? 0 : s.compIndex;
      segClosed[s.idx] = s.closed ? 1 : 0;
    }
    const m = this.links.length;
    const linkFrom = new Int32Array(m);
    const linkTo = new Int32Array(m);
    const linkFromEnd = new Uint8Array(m);
    const linkToEnd = new Uint8Array(m);
    for (let i = 0; i < m; i++) {
      const l = this.links[i];
      linkFrom[i] = l.a; linkTo[i] = l.b;
      linkFromEnd[i] = l.aEnd; linkToEnd[i] = l.bEnd;
    }
    return {
      nParticles: this.nParticles,
      nSegments: n,
      px: Float32Array.from(this.px),
      py: Float32Array.from(this.py),
      segP0, segK, segRest, segLen, segComp, segClosed,
      linkFrom, linkTo, linkFromEnd, linkToEnd,
      ...this.ringComponents(),
    };
  }

  /** Serialisable view state for `POST /api/session`. */
  exportLayout() {
    const positions = {};
    for (const seg of this.segments) {
      const xs = new Array(seg.k), ys = new Array(seg.k);
      for (let i = 0; i < seg.k; i++) {
        xs[i] = Math.round(this.px[seg.p0 + i] * 10) / 10;
        ys[i] = Math.round(this.py[seg.p0 + i] * 10) / 10;
      }
      positions[seg.name] = { x: xs, y: ys };
    }
    return { version: 1, positions };
  }

  /** Restore what `exportLayout()` produced. Unknown segments are ignored. */
  importLayout(layout) {
    if (!layout || !layout.positions) return 0;
    let n = 0;
    for (const [name, p] of Object.entries(layout.positions)) {
      const seg = this.byName.get(name);
      if (!seg || !p || !Array.isArray(p.x) || !Array.isArray(p.y)) continue;
      const count = Math.min(seg.k, p.x.length, p.y.length);
      for (let i = 0; i < count; i++) {
        this.px[seg.p0 + i] = Number(p.x[i]) || 0;
        this.py[seg.p0 + i] = Number(p.y[i]) || 0;
      }
      // If the stored polyline was shorter, extend along the last direction.
      for (let i = count; i < seg.k; i++) {
        this.px[seg.p0 + i] = this.px[seg.p0 + count - 1] + (i - count + 1);
        this.py[seg.p0 + i] = this.py[seg.p0 + count - 1];
      }
      n++;
    }
    this.updateBounds();
    return n;
  }
}

/** Pick the most convincing reference hit for a segment (longest block wins). */
export function bestRefHit(hits) {
  if (!Array.isArray(hits) || !hits.length) return null;
  let best = null, bestScore = -Infinity;
  for (const h of hits) {
    if (!h) continue;
    const len = Number(h.block_len);
    const span = Number.isFinite(len) ? len
      : Math.abs((Number(h.q_en) || 0) - (Number(h.q_st) || 0));
    const score = span * (h.is_primary ? 1.5 : 1) * (Number(h.identity) || 1);
    if (score > bestScore) { bestScore = score; best = h; }
  }
  return best;
}

/** Reference sequence lengths, harvested from `r_len` on the alignment hits. */
export function referenceLengths(graph) {
  const out = new Map();
  for (const seg of graph.segments) {
    for (const h of seg.refHits) {
      if (!h || !h.ref) continue;
      const len = Number(h.r_len);
      if (Number.isFinite(len) && len > 0) {
        const cur = out.get(String(h.ref)) || 0;
        if (len > cur) out.set(String(h.ref), len);
      }
    }
  }
  return out;
}

export default GraphModel;
