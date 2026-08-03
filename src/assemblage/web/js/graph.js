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

/** Geometry defaults; the UI can override `lengthScale` (node size slider). */
export const DEFAULT_GEOM = Object.freeze({
  minLen: 26,          // shortest polyline in world units
  maxLen: 900,         // longest polyline in world units
  lengthScale: 2.2,    // multiplies sqrt(bp)
  particleSpacing: 25, // world units of polyline per particle
  minParticles: 2,
  maxParticles: 12,
});

/**
 * Draw length of a segment: `clamp(minLen, maxLen, scale * sqrt(length))`.
 * Segments with an unknown/zero length still get the minimum length so they
 * remain clickable.
 */
export function drawLengthFor(bp, geom = DEFAULT_GEOM) {
  const L = Number(bp);
  const safe = Number.isFinite(L) && L > 0 ? L : 1;
  return clamp(geom.lengthScale * Math.sqrt(safe), geom.minLen, geom.maxLen);
}

/** Particle count for a polyline: `clamp(2, 12, round(drawLen / 25))`. */
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

/** Particle index for a given end of a segment record. */
export function endParticle(seg, whichEnd) {
  return whichEnd === END_END ? seg.p0 + seg.k - 1 : seg.p0;
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
    /** @type {Array<Object>} */
    this.paths = [];
    /** @type {Array<string>} reference sequence names seen in ref_hits */
    this.references = [];

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
   * Rebuild the model from a `GET /api/graph` payload.
   * `keepPositions` re-uses the previous coordinates of segments with the same
   * name, so a re-fetch after an operation does not scramble the picture.
   */
  setData(payload, { keepPositions = true } = {}) {
    const prev = keepPositions ? this._snapshotPositions() : null;
    const geom = this.geom;

    const rawSegs = (payload && Array.isArray(payload.segments)) ? payload.segments : [];
    const rawLinks = (payload && Array.isArray(payload.links)) ? payload.links : [];

    this.segments = [];
    this.byName = new Map();
    let p = 0;

    for (let i = 0; i < rawSegs.length; i++) {
      const s = rawSegs[i] || {};
      const name = String(s.name !== undefined && s.name !== null ? s.name : 'seg_' + i);
      if (this.byName.has(name)) continue; // duplicate names would break lookups
      const length = Number.isFinite(Number(s.length)) ? Number(s.length) : 0;
      const drawLen = drawLengthFor(length, geom);
      const k = particleCountFor(drawLen, geom);
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
        refHits: Array.isArray(s.ref_hits) ? s.ref_hits : [],
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
      });
    }

    this.paths = (payload && Array.isArray(payload.paths)) ? payload.paths : [];
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
    this._computeStats();
    this._restorePositions(prev);
    this.updateBounds();
    this.revision++;
    return this;
  }

  /* ------------------------------------------------------------ internals */

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
    let usable = false;
    for (const seg of this.segments) {
      if (seg.component !== 0 && seg.component !== null && seg.component !== undefined) usable = true;
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
   * Give every particle a sensible starting position: components are dropped
   * into a coarse grid, segments spiral outwards inside their component. A good
   * seed means the force layout converges in a few hundred iterations instead
   * of a few thousand.
   */
  seedPositions() {
    const comps = this.components;
    if (!comps.length) return;
    const cols = Math.max(1, Math.ceil(Math.sqrt(comps.length)));
    // Cell size follows the biggest component so nothing starts on top of
    // anything else.
    let cell = 0;
    for (const c of comps) {
      let span = 0;
      for (const si of c.segs) span += this.segments[si].drawLen;
      cell = Math.max(cell, Math.sqrt(Math.max(1, span)) * 6 + 120);
    }
    comps.forEach((comp, ci) => {
      const gx = (ci % cols) * cell;
      const gy = Math.floor(ci / cols) * cell;
      const n = comp.segs.length;
      // Golden-angle spiral: even coverage, no clumping.
      const golden = Math.PI * (3 - Math.sqrt(5));
      const spread = Math.sqrt(Math.max(1, n)) * 26 + 30;
      comp.segs.forEach((si, k) => {
        const seg = this.segments[si];
        const r = spread * Math.sqrt((k + 0.5) / n);
        const a = k * golden;
        this._placeSegment(seg, gx + Math.cos(a) * r, gy + Math.sin(a) * r, a + Math.PI / 2);
      });
    });
    this.updateBounds();
  }

  /** Recompute per-segment bounding boxes (used for culling and hit tests). */
  updateBounds() {
    const { px, py, bbox } = this;
    for (const seg of this.segments) {
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
  toLayoutArrays() {
    const n = this.segments.length;
    const segP0 = new Int32Array(n);
    const segK = new Int32Array(n);
    const segRest = new Float32Array(n);
    const segLen = new Float32Array(n);
    const segComp = new Int32Array(n);
    for (const s of this.segments) {
      segP0[s.idx] = s.p0;
      segK[s.idx] = s.k;
      segRest[s.idx] = s.k > 1 ? s.drawLen / (s.k - 1) : s.drawLen;
      segLen[s.idx] = s.drawLen;
      segComp[s.idx] = s.compIndex === undefined ? 0 : s.compIndex;
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
      segP0, segK, segRest, segLen, segComp,
      linkFrom, linkTo, linkFromEnd, linkToEnd,
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
