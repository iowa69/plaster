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
      seg.k = particleCountFor(seg.drawLen, this.geom);
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

    for (let i = 0; i < rawSegs.length; i++) {
      const s = rawSegs[i] || {};
      const name = String(s.name !== undefined && s.name !== null ? s.name : 'seg_' + i);
      if (this.byName.has(name)) continue; // duplicate names would break lookups
      const length = Number.isFinite(Number(s.length)) ? Number(s.length) : 0;
      const drawLen = drawLengthFor(length, geom, this.perBase);
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
   * Give every particle a sensible starting position: components are dropped
   * into a coarse grid, segments spiral outwards inside their component. A good
   * seed means the force layout converges in a few hundred iterations instead
   * of a few thousand.
   */
  seedPositions() {
    const comps = this.components;
    if (!comps.length) return;
    const cols = Math.max(1, Math.ceil(Math.sqrt(comps.length)));

    // Everything here is in drawn units, so it has to be derived from the drawn
    // sizes rather than from fixed constants: the bp -> unit scale is
    // calibrated per graph, so a segment can be five units long or five
    // thousand. Seeding at a fixed radius piles every contig on top of its
    // neighbours and leaves the force model to untangle a knot it did not need
    // to be given.
    const spreadOf = (comp) => {
      const n = Math.max(1, comp.segs.length);
      let total = 0;
      for (const si of comp.segs) total += this.segments[si].drawLen;
      const mean = total / n;
      return Math.sqrt(n) * Math.max(mean, 1) * 0.75 + mean;
    };

    // Cell size follows the biggest component so nothing starts on top of
    // anything else.
    let cell = 0;
    for (const c of comps) cell = Math.max(cell, spreadOf(c) * 2.4 + 40);

    comps.forEach((comp, ci) => {
      const gx = (ci % cols) * cell;
      const gy = Math.floor(ci / cols) * cell;
      const n = comp.segs.length;
      // Golden-angle spiral: even coverage, no clumping.
      const golden = Math.PI * (3 - Math.sqrt(5));
      const spread = spreadOf(comp);
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
