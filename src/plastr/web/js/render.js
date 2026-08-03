/**
 * render.js — Canvas 2D renderer, colour mapping, hit testing, interaction,
 * and PNG / SVG export.
 *
 * Design notes
 * ------------
 * * Every segment is a polyline through its layout particles. Its screen width
 *   comes from read depth (Bandage-style) and is expressed in world units so it
 *   scales with zoom.
 * * Only segments whose bounding box intersects the viewport are drawn, and
 *   they are batched by (colour, quantised width) so a 10k-segment graph costs a
 *   few dozen `stroke()` calls rather than 10k.
 * * `paint()` is target-agnostic: the live canvas, the PNG export canvas and the
 *   SVG serialiser all consume the same scene description, so what you export is
 *   exactly what you see.
 */

import { fmtBp } from './graph.js';

/* ====================================================================== */
/*  Colour ramps                                                           */
/* ====================================================================== */

const VIRIDIS = ['#440154', '#482878', '#3e4989', '#31688e', '#26828e',
  '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'];
const MAGMA = ['#000c24', '#1c1044', '#4f127b', '#812581', '#b5367a',
  '#e55964', '#fb8761', '#fec287', '#fcfdbf'];
const COOLWARM = ['#3b4cc0', '#6f91f2', '#a9c6fd', '#cfd4d8', '#f5c0a6',
  '#ee8468', '#b40426'];
const HEAT = ['#ffe08a', '#ffb14e', '#fa8775', '#ea5f94', '#cd34b5'];

export const CATEGORICAL = [
  '#4f9cf9', '#f2994a', '#27ae60', '#eb5757', '#bb6bd9', '#f2c94c',
  '#56ccf2', '#6fcf97', '#ff8fab', '#9b9bff', '#00c2a8', '#c98b3a',
  '#7ec8e3', '#d68fd6', '#8ecf6f', '#ff6b6b', '#a0d8b3', '#e8a5c4',
];

/** Number of quantised entries used for a continuous ramp. */
const RAMP_STEPS = 32;
/** Number of quantised line-width buckets. */
const WIDTH_BUCKETS = 12;

function hexToRgb(h) {
  const s = h.replace('#', '');
  const v = s.length === 3
    ? [s[0] + s[0], s[1] + s[1], s[2] + s[2]]
    : [s.slice(0, 2), s.slice(2, 4), s.slice(4, 6)];
  return [parseInt(v[0], 16) || 0, parseInt(v[1], 16) || 0, parseInt(v[2], 16) || 0];
}

function rgbStr(c) { return 'rgb(' + (c[0] | 0) + ',' + (c[1] | 0) + ',' + (c[2] | 0) + ')'; }

/** Sample a list of hex stops at t in [0, 1]. */
export function sampleRamp(stops, t) {
  const n = stops.length;
  if (n === 0) return '#888888';
  if (n === 1) return stops[0];
  const x = Math.max(0, Math.min(1, Number.isFinite(t) ? t : 0)) * (n - 1);
  const i = Math.min(n - 2, Math.floor(x));
  const f = x - i;
  const a = hexToRgb(stops[i]);
  const b = hexToRgb(stops[i + 1]);
  return rgbStr([a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f]);
}

function buildRamp(stops, steps = RAMP_STEPS) {
  const out = new Array(steps);
  for (let i = 0; i < steps; i++) out[i] = sampleRamp(stops, steps === 1 ? 0 : i / (steps - 1));
  return out;
}

/** Deterministic pseudo-random colour for an arbitrary key. */
export function hashColour(key, sat = 62, light = 58) {
  const s = String(key);
  let h = 2166136261;
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  const hue = Math.abs(h) % 360;
  return `hsl(${hue} ${sat}% ${light}%)`;
}

/** Read the theme tokens the renderer needs out of the stylesheet. */
export function readTheme() {
  const fallback = {
    canvasBg: '#0b0e13', grid: '#161c25', node: '#7f8ea3', link: '#5c6b82',
    sel: '#ffd166', text: '#dbe3ee', dim: '#8b97a8', faint: '#5f6b7c', accent: '#4f9cf9',
    outline: 'rgba(28,36,48,0.85)',
  };
  try {
    const cs = getComputedStyle(document.documentElement);
    const g = (n, fb) => ((cs.getPropertyValue(n) || '').trim() || fb);
    return {
      canvasBg: g('--canvas-bg', fallback.canvasBg),
      grid: g('--canvas-grid', fallback.grid),
      node: g('--node-default', fallback.node),
      link: g('--link-strong', fallback.link),
      sel: g('--sel', fallback.sel),
      text: g('--text', fallback.text),
      dim: g('--text-dim', fallback.dim),
      faint: g('--text-faint', fallback.faint),
      accent: g('--accent', fallback.accent),
      outline: g('--node-outline', fallback.outline),
    };
  } catch {
    return fallback;
  }
}

/* ====================================================================== */
/*  Colour mapper                                                          */
/* ====================================================================== */

export class ColourMapper {
  constructor() {
    this.mode = 'uniform';
    this.palette = ['#7f8ea3'];
    this.segPal = new Int32Array(0);
    this.legend = { type: 'none' };
  }

  setMode(mode) { this.mode = mode || 'uniform'; }

  /**
   * Recompute `palette` + `segPal` (an index per segment) and the legend
   * description that the UI renders.
   */
  update(graph, theme) {
    const n = graph.segments.length;
    this.segPal = new Int32Array(n);
    const grey = theme.faint;

    switch (this.mode) {
      case 'depth': return this._continuous(graph, theme, {
        ramp: VIRIDIS,
        label: 'depth (x)',
        log: true,
        value: (s) => (s.depth === null ? null : s.depth),
        lo: graph.stats.depth.min,
        hi: graph.stats.depth.max,
        available: graph.stats.depth.has,
        missingMsg: 'no depth in this graph',
        format: (v) => (v >= 100 ? v.toFixed(0) : v.toFixed(1)) + '×',
      });

      case 'gc': return this._continuous(graph, theme, {
        ramp: COOLWARM,
        label: 'GC fraction',
        log: false,
        value: (s) => (s.gc === null ? null : s.gc),
        lo: graph.stats.gc.min,
        hi: graph.stats.gc.max,
        available: graph.stats.gc.has,
        missingMsg: 'no GC content in this graph',
        format: (v) => (v * 100).toFixed(1) + '%',
      });

      case 'length': return this._continuous(graph, theme, {
        ramp: MAGMA,
        label: 'segment length',
        log: true,
        value: (s) => s.length,
        lo: graph.stats.length.min,
        hi: graph.stats.length.max,
        available: n > 0,
        missingMsg: 'no segments',
        format: (v) => fmtBp(v),
      });

      case 'component': {
        const comps = graph.components;
        this.palette = comps.map((_, i) => CATEGORICAL[i % CATEGORICAL.length]);
        if (!this.palette.length) this.palette = [theme.node];
        for (const s of graph.segments) this.segPal[s.idx] = (s.compIndex || 0) % this.palette.length;
        this.legend = {
          type: 'cat',
          label: 'connected component',
          items: comps.slice(0, 14).map((c, i) => ({
            colour: this.palette[i % this.palette.length],
            label: `#${i + 1} · ${c.segs.length} seg · ${fmtBp(c.length)}`,
          })),
          more: Math.max(0, comps.length - 14),
        };
        return this;
      }

      case 'random': {
        // One colour per node, seeded from its name so it is stable across
        // reloads. This is Bandage's default and it is the most legible way to
        // follow an individual contig through a tangle.
        this.palette = graph.segments.map((s) => hashColour('seg:' + s.name, 58, 62));
        if (!this.palette.length) this.palette = [theme.node];
        for (const s of graph.segments) this.segPal[s.idx] = s.idx;
        this.legend = {
          type: 'none',
          label: 'random per segment',
          note: 'each contig gets its own colour, stable across reloads',
        };
        return this;
      }

      case 'randomComponent': {
        const comps = graph.components;
        this.palette = comps.map((c) => hashColour('comp:' + c.id));
        if (!this.palette.length) this.palette = [theme.node];
        for (const s of graph.segments) this.segPal[s.idx] = (s.compIndex || 0) % this.palette.length;
        this.legend = {
          type: 'cat',
          label: 'random per component',
          items: comps.slice(0, 10).map((c, i) => ({
            colour: this.palette[i % this.palette.length],
            label: `#${i + 1} · ${c.segs.length} seg`,
          })),
          more: Math.max(0, comps.length - 10),
        };
        return this;
      }

      case 'reference': {
        const refs = graph.references;
        this.palette = [grey].concat(refs.map((r, i) => CATEGORICAL[i % CATEGORICAL.length]));
        const idx = new Map();
        refs.forEach((r, i) => idx.set(r, i + 1));
        let unaligned = 0;
        for (const s of graph.segments) {
          const hit = s.bestHit;
          const p = hit && hit.ref !== undefined ? idx.get(String(hit.ref)) : undefined;
          this.segPal[s.idx] = p === undefined ? 0 : p;
          if (p === undefined) unaligned++;
        }
        this.legend = {
          type: 'cat',
          label: refs.length ? 'reference sequence' : 'reference (none loaded)',
          items: [{ colour: grey, label: `unaligned (${unaligned})` }].concat(
            refs.slice(0, 13).map((r, i) => ({ colour: this.palette[i + 1], label: r })),
          ),
          more: Math.max(0, refs.length - 13),
          note: refs.length ? '' : 'Insert a reference to colour by chromosome.',
        };
        return this;
      }

      case 'search': {
        const ramp = buildRamp(HEAT, 8);
        this.palette = [grey].concat(ramp);
        const hits = graph.searchHits;
        let nHit = 0;
        for (const s of graph.segments) {
          const h = hits.get(s.name);
          if (!h) { this.segPal[s.idx] = 0; continue; }
          nHit++;
          const id = Number(h.identity);
          const t = Number.isFinite(id) ? Math.max(0, Math.min(1, (id - 0.7) / 0.3)) : 1;
          this.segPal[s.idx] = 1 + Math.round(t * (ramp.length - 1));
        }
        this.legend = {
          type: 'scale',
          label: `search hits (${nHit})`,
          ramp,
          ticks: ['70%', '85%', '100%'],
          note: nHit ? '' : 'Run a sequence search to highlight hits.',
        };
        return this;
      }

      case 'uniform':
      default:
        this.palette = [theme.node];
        this.segPal.fill(0);
        this.legend = { type: 'none', label: 'uniform' };
        return this;
    }
  }

  _continuous(graph, theme, cfg) {
    const ramp = buildRamp(cfg.ramp);
    this.palette = [theme.faint].concat(ramp);
    if (!cfg.available) {
      this.segPal.fill(0);
      this.legend = { type: 'none', label: cfg.label, note: cfg.missingMsg };
      return this;
    }
    let lo = Number(cfg.lo), hi = Number(cfg.hi);
    if (!Number.isFinite(lo)) lo = 0;
    if (!Number.isFinite(hi) || hi <= lo) hi = lo + 1;
    const tf = cfg.log
      ? (v) => (Math.log(Math.max(v, 1e-6)) - Math.log(Math.max(lo, 1e-6))) /
               (Math.log(Math.max(hi, 1e-6)) - Math.log(Math.max(lo, 1e-6)) || 1)
      : (v) => (v - lo) / (hi - lo);
    for (const s of graph.segments) {
      const v = cfg.value(s);
      if (v === null || !Number.isFinite(v)) { this.segPal[s.idx] = 0; continue; }
      const t = Math.max(0, Math.min(1, tf(v)));
      this.segPal[s.idx] = 1 + Math.round(t * (ramp.length - 1));
    }
    const mid = cfg.log ? Math.sqrt(Math.max(lo, 1e-6) * hi) : (lo + hi) / 2;
    this.legend = {
      type: 'scale',
      label: cfg.label + (cfg.log ? ' (log)' : ''),
      ramp,
      ticks: [cfg.format(lo), cfg.format(mid), cfg.format(hi)],
      missingColour: theme.faint,
    };
    return this;
  }
}

/* ====================================================================== */
/*  small geometry helpers                                                 */
/* ====================================================================== */

function distPtSeg(px, py, x1, y1, x2, y2) {
  const dx = x2 - x1, dy = y2 - y1;
  const l2 = dx * dx + dy * dy;
  let t = l2 ? ((px - x1) * dx + (py - y1) * dy) / l2 : 0;
  if (t < 0) t = 0; else if (t > 1) t = 1;
  const cx = x1 + t * dx, cy = y1 + t * dy;
  return Math.hypot(px - cx, py - cy);
}

/** Liang-Barsky: does the line segment touch the axis-aligned rectangle? */
function segHitsRect(x1, y1, x2, y2, rx0, ry0, rx1, ry1) {
  if ((x1 >= rx0 && x1 <= rx1 && y1 >= ry0 && y1 <= ry1) ||
      (x2 >= rx0 && x2 <= rx1 && y2 >= ry0 && y2 <= ry1)) return true;
  const dx = x2 - x1, dy = y2 - y1;
  const p = [-dx, dx, -dy, dy];
  const q = [x1 - rx0, rx1 - x1, y1 - ry0, ry1 - y1];
  let t0 = 0, t1 = 1;
  for (let i = 0; i < 4; i++) {
    if (p[i] === 0) { if (q[i] < 0) return false; continue; }
    const r = q[i] / p[i];
    if (p[i] < 0) { if (r > t1) return false; if (r > t0) t0 = r; }
    else { if (r < t0) return false; if (r < t1) t1 = r; }
  }
  return true;
}

/* ====================================================================== */
/*  Renderer                                                               */
/* ====================================================================== */

export class Renderer {
  constructor(canvas, graph) {
    this.canvas = canvas;
    this.graph = graph;
    this.ctx = canvas.getContext('2d', { alpha: false });
    this.theme = readTheme();
    this.colour = new ColourMapper();

    /** world point at the centre of the viewport + pixels-per-world-unit */
    this.view = { cx: 0, cy: 0, scale: 1 };

    this.opts = {
      showLinks: true,
      showArrows: false,
      showLabels: false,
      showGrid: false,
      depthWidth: true,
      widthScale: 1,
      // World units, constant like Bandage's node width. Paired with the
      // absolute length scale this makes a short contig a stub and a long one
      // a ribbon, at the same drawn size in any graph.
      baseWidth: 5.5,
    };

    this.selected = new Set();   // segment indices
    this.hoverIdx = -1;
    this.boxRect = null;         // screen-space rectangle while box selecting

    this.segWidth = new Float32Array(0);
    this.segBucket = new Int32Array(0);
    this.bucketWidth = new Float32Array(WIDTH_BUCKETS);

    this._buckets = [];
    this._usedKeys = [];
    this._visible = [];
    this._visibleLinks = [];
    this._dpr = 1;
    this.width = 1;
    this.height = 1;

    this.handlers = {};
    this._frame = 0;
    this._needsDraw = false;
    this._spaceDown = false;
    this._drag = null;
    this._pointerId = null;

    this.resize();
  }

  /* ------------------------------------------------------------ sizing */

  resize() {
    const rect = this.canvas.getBoundingClientRect();
    const w = Math.max(1, Math.round(rect.width));
    const h = Math.max(1, Math.round(rect.height));
    const dpr = Math.min(3, window.devicePixelRatio || 1);
    this.width = w; this.height = h; this._dpr = dpr;
    this.canvas.width = Math.round(w * dpr);
    this.canvas.height = Math.round(h * dpr);
    this.requestDraw();
  }

  refreshTheme() {
    this.theme = readTheme();
    this.updateStyle();
    this.requestDraw();
  }

  /* ------------------------------------------------------- coordinates */

  worldToScreenX(x) { return (x - this.view.cx) * this.view.scale + this.width / 2; }
  worldToScreenY(y) { return (y - this.view.cy) * this.view.scale + this.height / 2; }
  screenToWorld(sx, sy) {
    return [
      (sx - this.width / 2) / this.view.scale + this.view.cx,
      (sy - this.height / 2) / this.view.scale + this.view.cy,
    ];
  }

  /** Zoom keeping the world point under (sx, sy) fixed. */
  zoomAt(sx, sy, factor) {
    const [wx, wy] = this.screenToWorld(sx, sy);
    const s = Math.max(1e-4, Math.min(80, this.view.scale * factor));
    this.view.scale = s;
    this.view.cx = wx - (sx - this.width / 2) / s;
    this.view.cy = wy - (sy - this.height / 2) / s;
    this._emit('viewchange');
    this.requestDraw();
  }

  panBy(dxScreen, dyScreen) {
    this.view.cx -= dxScreen / this.view.scale;
    this.view.cy -= dyScreen / this.view.scale;
    this._emit('viewchange');
    this.requestDraw();
  }

  /** Fit the given segment indices (or the whole graph) into the viewport. */
  fitToView(indices, padding = 0.08) {
    const g = this.graph;
    if (g.isEmpty) { this.view = { cx: 0, cy: 0, scale: 1 }; this.requestDraw(); return; }
    const [minX, minY, maxX, maxY] = g.boundsOf(indices && indices.length ? indices : null);
    const w = Math.max(1, maxX - minX);
    const h = Math.max(1, maxY - minY);
    const s = Math.min(this.width / (w * (1 + padding * 2)), this.height / (h * (1 + padding * 2)));
    this.view.scale = Math.max(1e-4, Math.min(20, s));
    this.view.cx = (minX + maxX) / 2;
    this.view.cy = (minY + maxY) / 2;
    this._emit('viewchange');
    this.requestDraw();
  }

  centreOn(segIdx) {
    const g = this.graph;
    const seg = g.segments[segIdx];
    if (!seg) return;
    const o = segIdx * 4;
    this.view.cx = (g.bbox[o] + g.bbox[o + 2]) / 2;
    this.view.cy = (g.bbox[o + 1] + g.bbox[o + 3]) / 2;
    this._emit('viewchange');
    this.requestDraw();
  }

  /* ------------------------------------------------------------- style */

  /** Recompute per-segment colour index and line width. Call after any change. */
  updateStyle() {
    const g = this.graph;
    this.colour.update(g, this.theme);
    const n = g.segments.length;
    this.segWidth = new Float32Array(n);
    this.segBucket = new Int32Array(n);
    if (!n) return;

    // Median depth gives a stable reference point for the width scale, so a
    // handful of very deep repeat nodes don't flatten everything else.
    let median = 1;
    if (this.opts.depthWidth && g.stats.depth.has) {
      const ds = [];
      for (const s of g.segments) if (s.depth !== null && Number.isFinite(s.depth) && s.depth > 0) ds.push(s.depth);
      if (ds.length) { ds.sort((a, b) => a - b); median = ds[ds.length >> 1] || 1; }
    }
    const base = this.opts.baseWidth * this.opts.widthScale;
    let minW = Infinity, maxW = -Infinity;
    for (const s of g.segments) {
      let f = 1;
      if (this.opts.depthWidth && g.stats.depth.has && s.depth !== null && Number.isFinite(s.depth)) {
        f = Math.sqrt(Math.max(s.depth, 1e-3) / median);
        f = Math.max(0.35, Math.min(3.2, f));
      }
      const w = base * f;
      this.segWidth[s.idx] = w;
      if (w < minW) minW = w;
      if (w > maxW) maxW = w;
    }
    if (!Number.isFinite(minW)) { minW = base; maxW = base; }
    if (maxW <= minW) maxW = minW + 1e-6;
    for (let b = 0; b < WIDTH_BUCKETS; b++) {
      this.bucketWidth[b] = minW + (maxW - minW) * (b / (WIDTH_BUCKETS - 1));
    }
    const span = maxW - minW;
    for (let i = 0; i < n; i++) {
      const t = span > 0 ? (this.segWidth[i] - minW) / span : 0;
      this.segBucket[i] = Math.max(0, Math.min(WIDTH_BUCKETS - 1, Math.round(t * (WIDTH_BUCKETS - 1))));
    }
  }

  setOption(key, value) {
    this.opts[key] = value;
    if (key === 'depthWidth' || key === 'widthScale' || key === 'baseWidth') this.updateStyle();
    this.requestDraw();
  }

  setColourMode(mode) {
    this.colour.setMode(mode);
    this.updateStyle();
    this.requestDraw();
    return this.colour.legend;
  }

  /* --------------------------------------------------------- selection */

  setSelection(indices) {
    this.selected = new Set(indices);
    this.requestDraw();
    this._emit('selection');
  }

  selectNames(names) {
    const g = this.graph;
    const set = new Set();
    for (const n of names || []) {
      const s = g.segmentByName(n);
      if (s) set.add(s.idx);
    }
    this.selected = set;
    this.requestDraw();
    this._emit('selection');
  }

  selectedNames() {
    const out = [];
    for (const i of this.selected) {
      const s = this.graph.segments[i];
      if (s) out.push(s.name);
    }
    return out;
  }

  clearSelection() {
    if (!this.selected.size) return;
    this.selected = new Set();
    this.requestDraw();
    this._emit('selection');
  }

  selectAllVisible() {
    this.selected = new Set(this._visible.length ? this._visible : this.graph.segments.map((s) => s.idx));
    this.requestDraw();
    this._emit('selection');
  }

  /* ------------------------------------------------------- scene / draw */

  /** Compute the visible segment and link lists for the current viewport. */
  _scene(W, H) {
    const g = this.graph;
    const s = this.view.scale;
    const halfW = W / (2 * s), halfH = H / (2 * s);
    // Margin covers half the widest line so nothing pops at the edge.
    const margin = 40 / s;
    const x0 = this.view.cx - halfW - margin;
    const x1 = this.view.cx + halfW + margin;
    const y0 = this.view.cy - halfH - margin;
    const y1 = this.view.cy + halfH + margin;

    const vis = this._visible;
    vis.length = 0;
    const bb = g.bbox;
    for (let i = 0; i < g.segments.length; i++) {
      const o = i * 4;
      if (bb[o + 2] < x0 || bb[o] > x1 || bb[o + 3] < y0 || bb[o + 1] > y1) continue;
      vis.push(i);
    }
    const visSet = vis.length === g.segments.length ? null : new Set(vis);

    const vlinks = this._visibleLinks;
    vlinks.length = 0;
    if (this.opts.showLinks) {
      for (let i = 0; i < g.links.length; i++) {
        const l = g.links[i];
        if (visSet && !visSet.has(l.a) && !visSet.has(l.b)) continue;
        vlinks.push(i);
      }
    }
    return { vis, vlinks, rect: [x0, y0, x1, y1] };
  }

  /**
   * Append a segment's polyline to the current path. Written without closures
   * because this runs once per visible segment per frame.
   */
  _segPath(ctx, si, smooth) {
    const g = this.graph;
    const seg = g.segments[si];
    const p0 = seg.p0, k = seg.k;
    const s = this.view.scale;
    const ox = this.width / 2 - this.view.cx * s;
    const oy = this.height / 2 - this.view.cy * s;
    const px = g.px, py = g.py;
    if (k < 2) {
      const x = px[p0] * s + ox, y = py[p0] * s + oy;
      ctx.moveTo(x, y);
      ctx.lineTo(x + 0.01, y);
      return;
    }
    if (!smooth || k === 2) {
      ctx.moveTo(px[p0] * s + ox, py[p0] * s + oy);
      for (let i = 1; i < k; i++) ctx.lineTo(px[p0 + i] * s + ox, py[p0 + i] * s + oy);
      return;
    }
    // Quadratic spline through particle midpoints: smooth without overshoot.
    ctx.moveTo(px[p0] * s + ox, py[p0] * s + oy);
    for (let i = 1; i < k - 1; i++) {
      const cxp = px[p0 + i] * s + ox, cyp = py[p0 + i] * s + oy;
      const nxp = px[p0 + i + 1] * s + ox, nyp = py[p0 + i + 1] * s + oy;
      ctx.quadraticCurveTo(cxp, cyp, (cxp + nxp) / 2, (cyp + nyp) / 2);
    }
    ctx.lineTo(px[p0 + k - 1] * s + ox, py[p0 + k - 1] * s + oy);
  }

  /** Screen-space polyline points for a segment (used by the SVG writer). */
  _segPoints(si) {
    const g = this.graph;
    const seg = g.segments[si];
    const pts = [];
    for (let i = 0; i < seg.k; i++) {
      pts.push([this.worldToScreenX(g.px[seg.p0 + i]), this.worldToScreenY(g.py[seg.p0 + i])]);
    }
    return pts;
  }

  requestDraw() {
    if (this._needsDraw) return;
    this._needsDraw = true;
    this._frame = requestAnimationFrame(() => {
      this._needsDraw = false;
      this.draw();
    });
  }

  draw() {
    this.paint(this.ctx, this.width, this.height, this._dpr, { interactive: true });
  }

  /**
   * Render the whole scene into any 2D context.
   * @param {CanvasRenderingContext2D} ctx
   * @param {number} W logical width in CSS pixels
   * @param {number} H logical height in CSS pixels
   * @param {number} ratio device pixel ratio / export multiplier
   */
  paint(ctx, W, H, ratio, o = {}) {
    const g = this.graph;
    const th = this.theme;
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.fillStyle = th.canvasBg;
    ctx.fillRect(0, 0, W, H);
    if (g.isEmpty) return;

    const scene = this._scene(W, H);
    const scale = this.view.scale;

    if (this.opts.showGrid) this._paintGrid(ctx, W, H);

    /* ---- links ---- */
    if (this.opts.showLinks && scene.vlinks.length) {
      ctx.strokeStyle = th.link;
      ctx.globalAlpha = 0.55;
      ctx.lineWidth = Math.max(0.5, Math.min(3, 1.1 * Math.sqrt(scale)));
      ctx.beginPath();
      for (const li of scene.vlinks) {
        const l = g.links[li];
        const ax = this.worldToScreenX(g.px[l.pa]);
        const ay = this.worldToScreenY(g.py[l.pa]);
        const bx = this.worldToScreenX(g.px[l.pb]);
        const by = this.worldToScreenY(g.py[l.pb]);
        ctx.moveTo(ax, ay);
        ctx.lineTo(bx, by);
      }
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    /* ---- selection halo (drawn under the segments) ---- */
    if (this.selected.size) {
      ctx.strokeStyle = th.sel;
      ctx.globalAlpha = 0.5;
      ctx.lineCap = 'round';
      ctx.lineJoin = 'round';
      for (const si of scene.vis) {
        if (!this.selected.has(si)) continue;
        const lw = Math.max(3, Math.min(90, this.segWidth[si] * scale + 6));
        ctx.lineWidth = lw;
        ctx.beginPath();
        this._segPath(ctx, si, this.graph.segments[si].drawLen * scale > 30);
        ctx.stroke();
      }
      ctx.globalAlpha = 1;
    }

    /* ---- segments, batched by (colour, width bucket) ---- */
    const pal = this.colour.palette;
    const nPal = pal.length;
    const keys = this._usedKeys;
    keys.length = 0;
    const buckets = this._buckets;
    const needed = nPal * WIDTH_BUCKETS;
    while (buckets.length < needed) buckets.push([]);
    for (const si of scene.vis) {
      const key = (this.colour.segPal[si] % nPal) * WIDTH_BUCKETS + this.segBucket[si];
      const arr = buckets[key];
      if (arr.length === 0) keys.push(key);
      arr.push(si);
    }
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';

    // Two passes per width bucket: a dark outline, then the coloured body on
    // top. That is what gives a Bandage node its ribbon look and keeps
    // neighbouring contigs of similar colour distinguishable.
    const outline = o.outline !== false && this.opts.outlineNodes !== false;
    for (const pass of outline ? [0, 1] : [1]) {
      for (const key of keys) {
        const arr = buckets[key];
        const colIdx = Math.floor(key / WIDTH_BUCKETS);
        const wIdx = key % WIDTH_BUCKETS;
        // Keep a floor in *screen* pixels: zoomed out to fit a whole assembly
        // a world-space width collapses to a hairline and the ribbons that
        // carry the length information stop reading as ribbons at all.
        const body = Math.max(2.4, Math.min(90, this.bucketWidth[wIdx] * scale));
        if (pass === 0) {
          // Skip the outline when the body is too thin for it to read.
          if (body < 2.2) continue;
          ctx.strokeStyle = th.outline || 'rgba(20,26,34,0.85)';
          ctx.lineWidth = body + Math.min(2.6, Math.max(1, body * 0.28));
        } else {
          ctx.strokeStyle = pal[colIdx] || th.node;
          ctx.lineWidth = body;
        }
        ctx.beginPath();
        for (const si of arr) {
          this._segPath(ctx, si, g.segments[si].drawLen * scale > 30);
        }
        ctx.stroke();
      }
    }
    for (const key of keys) buckets[key].length = 0;

    /* ---- hover highlight ---- */
    if (o.interactive && this.hoverIdx >= 0 && this.hoverIdx < g.segments.length) {
      ctx.strokeStyle = th.accent;
      ctx.lineWidth = Math.max(1.2, Math.min(90, this.segWidth[this.hoverIdx] * scale + 1.5));
      ctx.globalAlpha = 0.9;
      ctx.beginPath();
      this._segPath(ctx, this.hoverIdx, true);
      ctx.stroke();
      ctx.globalAlpha = 1;
    }

    /* ---- direction arrows ---- */
    if (this.opts.showArrows && scene.vis.length <= 4000 && scale > 0.08) {
      ctx.fillStyle = th.dim;
      for (const si of scene.vis) {
        const seg = g.segments[si];
        if (seg.drawLen * scale < 22) continue;
        const pEnd = seg.p0 + seg.k - 1;
        const pPrev = seg.p0 + Math.max(0, seg.k - 2);
        const x2 = this.worldToScreenX(g.px[pEnd]);
        const y2 = this.worldToScreenY(g.py[pEnd]);
        const x1 = this.worldToScreenX(g.px[pPrev]);
        const y1 = this.worldToScreenY(g.py[pPrev]);
        let dx = x2 - x1, dy = y2 - y1;
        const d = Math.hypot(dx, dy) || 1;
        dx /= d; dy /= d;
        const size = Math.max(3.5, Math.min(11, this.segWidth[si] * scale * 0.9 + 2.5));
        ctx.beginPath();
        ctx.moveTo(x2 + dx * size, y2 + dy * size);
        ctx.lineTo(x2 - dy * size * 0.55, y2 + dx * size * 0.55);
        ctx.lineTo(x2 + dy * size * 0.55, y2 - dx * size * 0.55);
        ctx.closePath();
        ctx.fill();
      }
    }

    /* ---- labels ---- */
    if (this.opts.showLabels && scene.vis.length <= 600) {
      ctx.fillStyle = th.text;
      ctx.strokeStyle = th.canvasBg;
      ctx.lineWidth = 3;
      ctx.font = '11px system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      for (const si of scene.vis) {
        const seg = g.segments[si];
        if (seg.drawLen * scale < 46) continue;
        const mid = seg.p0 + (seg.k >> 1);
        const x = this.worldToScreenX(g.px[mid]);
        const y = this.worldToScreenY(g.py[mid]) - Math.max(8, this.segWidth[si] * scale * 0.6 + 7);
        const label = seg.name;
        ctx.strokeText(label, x, y);
        ctx.fillText(label, x, y);
      }
    }

    /* ---- box-select rubber band ---- */
    if (o.interactive && this.boxRect) {
      const r = this.boxRect;
      ctx.strokeStyle = th.accent;
      ctx.fillStyle = th.accent;
      ctx.globalAlpha = 0.12;
      ctx.fillRect(r.x, r.y, r.w, r.h);
      ctx.globalAlpha = 0.9;
      ctx.lineWidth = 1;
      ctx.setLineDash([4, 3]);
      ctx.strokeRect(r.x + 0.5, r.y + 0.5, r.w, r.h);
      ctx.setLineDash([]);
      ctx.globalAlpha = 1;
    }
  }

  _paintGrid(ctx, W, H) {
    const s = this.view.scale;
    // Pick a world spacing that lands between 60 and 600 screen pixels.
    let spacing = Math.pow(10, Math.ceil(Math.log10(80 / s)));
    if (spacing * s > 600) spacing /= 2;
    if (spacing * s < 40) spacing *= 2;
    const [wx0, wy0] = this.screenToWorld(0, 0);
    const [wx1, wy1] = this.screenToWorld(W, H);
    const startX = Math.floor(wx0 / spacing) * spacing;
    const startY = Math.floor(wy0 / spacing) * spacing;
    const cols = Math.ceil((wx1 - wx0) / spacing) + 1;
    const rows = Math.ceil((wy1 - wy0) / spacing) + 1;
    if (cols > 400 || rows > 400 || !Number.isFinite(cols) || !Number.isFinite(rows)) return;
    ctx.strokeStyle = this.theme.grid;
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let i = 0; i <= cols; i++) {
      const x = Math.round(this.worldToScreenX(startX + i * spacing)) + 0.5;
      ctx.moveTo(x, 0); ctx.lineTo(x, H);
    }
    for (let j = 0; j <= rows; j++) {
      const y = Math.round(this.worldToScreenY(startY + j * spacing)) + 0.5;
      ctx.moveTo(0, y); ctx.lineTo(W, y);
    }
    ctx.stroke();
  }

  /* --------------------------------------------------------- hit tests */

  /** Segment index under a screen point, or -1. */
  hitTest(sx, sy, tolPx = 7) {
    const g = this.graph;
    if (g.isEmpty) return -1;
    const [wx, wy] = this.screenToWorld(sx, sy);
    const tol = tolPx / this.view.scale;
    const candidates = this._visible.length ? this._visible : g.segments.map((s) => s.idx);
    let best = -1, bestD = Infinity;
    for (const si of candidates) {
      const halfW = this.segWidth[si] / 2;
      const reach = tol + halfW;
      const o = si * 4;
      if (wx < g.bbox[o] - reach || wx > g.bbox[o + 2] + reach ||
          wy < g.bbox[o + 1] - reach || wy > g.bbox[o + 3] + reach) continue;
      const seg = g.segments[si];
      let d = Infinity;
      if (seg.k === 1) {
        d = Math.hypot(wx - g.px[seg.p0], wy - g.py[seg.p0]);
      } else {
        for (let i = 0; i + 1 < seg.k; i++) {
          const a = seg.p0 + i;
          const dd = distPtSeg(wx, wy, g.px[a], g.py[a], g.px[a + 1], g.py[a + 1]);
          if (dd < d) d = dd;
        }
      }
      d -= halfW;
      if (d < tol && d < bestD) { bestD = d; best = si; }
    }
    return best;
  }

  /** Segment indices intersecting a screen-space rectangle. */
  boxSelect(rect) {
    const g = this.graph;
    const [x0, y0] = this.screenToWorld(rect.x, rect.y);
    const [x1, y1] = this.screenToWorld(rect.x + rect.w, rect.y + rect.h);
    const rx0 = Math.min(x0, x1), rx1 = Math.max(x0, x1);
    const ry0 = Math.min(y0, y1), ry1 = Math.max(y0, y1);
    const out = [];
    for (let si = 0; si < g.segments.length; si++) {
      const o = si * 4;
      if (g.bbox[o + 2] < rx0 || g.bbox[o] > rx1 || g.bbox[o + 3] < ry0 || g.bbox[o + 1] > ry1) continue;
      const seg = g.segments[si];
      let hit = false;
      if (seg.k === 1) {
        const x = g.px[seg.p0], y = g.py[seg.p0];
        hit = x >= rx0 && x <= rx1 && y >= ry0 && y <= ry1;
      } else {
        for (let i = 0; i + 1 < seg.k && !hit; i++) {
          const a = seg.p0 + i;
          hit = segHitsRect(g.px[a], g.py[a], g.px[a + 1], g.py[a + 1], rx0, ry0, rx1, ry1);
        }
      }
      if (hit) out.push(si);
    }
    return out;
  }

  /* -------------------------------------------------------- interaction */

  /**
   * Wire up pointer/wheel/keyboard interaction.
   * `handlers` may define: select(indices, additive), hover(index, ev),
   * dragEnd(indices), viewchange(), dblclick(index).
   */
  attach(handlers = {}) {
    this.handlers = handlers;
    const c = this.canvas;

    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      const r = c.getBoundingClientRect();
      const f = Math.exp(-e.deltaY * (e.deltaMode === 1 ? 0.03 : 0.0015));
      this.zoomAt(e.clientX - r.left, e.clientY - r.top, f);
    }, { passive: false });

    c.addEventListener('pointerdown', (e) => this._onDown(e));
    c.addEventListener('pointermove', (e) => this._onMove(e));
    c.addEventListener('pointerup', (e) => this._onUp(e));
    c.addEventListener('pointercancel', (e) => this._onUp(e));
    c.addEventListener('pointerleave', () => {
      if (this.hoverIdx !== -1) { this.hoverIdx = -1; this.requestDraw(); }
      this._emit('hover', -1, null);
    });
    c.addEventListener('dblclick', (e) => {
      const r = c.getBoundingClientRect();
      const si = this.hitTest(e.clientX - r.left, e.clientY - r.top);
      if (si >= 0) { this.centreOn(si); this._emit('dblclick', si); }
      else this.fitToView();
    });
    c.addEventListener('contextmenu', (e) => e.preventDefault());

    window.addEventListener('keydown', (e) => {
      if (e.code === 'Space' && !this._isTypingTarget(e.target)) {
        this._spaceDown = true;
        c.classList.add('pannable');
        e.preventDefault();
      }
    });
    window.addEventListener('keyup', (e) => {
      if (e.code === 'Space') { this._spaceDown = false; c.classList.remove('pannable'); }
    });
    window.addEventListener('blur', () => { this._spaceDown = false; c.classList.remove('pannable'); });
  }

  _isTypingTarget(t) {
    if (!t || !t.tagName) return false;
    const tag = t.tagName.toLowerCase();
    return tag === 'input' || tag === 'textarea' || tag === 'select' || t.isContentEditable;
  }

  _localPoint(e) {
    const r = this.canvas.getBoundingClientRect();
    return [e.clientX - r.left, e.clientY - r.top];
  }

  _onDown(e) {
    if (this.graph.isEmpty) return;
    const [sx, sy] = this._localPoint(e);
    const si = this.hitTest(sx, sy);
    this._pointerId = e.pointerId;
    try { this.canvas.setPointerCapture(e.pointerId); } catch { /* not fatal */ }
    this.canvas.focus({ preventScroll: true });

    const wantPan = e.button === 1 || this._spaceDown || (e.button === 0 && si < 0 && !e.shiftKey);
    const wantBox = e.button === 0 && si < 0 && e.shiftKey;

    if (wantPan) {
      this._drag = { mode: 'pan', sx, sy };
      this.canvas.classList.add('panning');
    } else if (wantBox) {
      this._drag = { mode: 'box', sx, sy, additive: true };
      this.boxRect = { x: sx, y: sy, w: 0, h: 0 };
    } else if (si >= 0) {
      if (e.shiftKey) {
        if (this.selected.has(si)) this.selected.delete(si); else this.selected.add(si);
        this._emit('selection');
      } else if (!this.selected.has(si)) {
        this.selected = new Set([si]);
        this._emit('selection');
      }
      const [wx, wy] = this.screenToWorld(sx, sy);
      this._drag = { mode: 'node', wx, wy, moved: false, si };
      this.requestDraw();
    } else {
      this._drag = { mode: 'pan', sx, sy };
      this.canvas.classList.add('panning');
    }
  }

  _onMove(e) {
    if (this.graph.isEmpty) return;
    const [sx, sy] = this._localPoint(e);
    const d = this._drag;

    if (!d) {
      const si = this.hitTest(sx, sy);
      if (si !== this.hoverIdx) {
        this.hoverIdx = si;
        this.canvas.classList.toggle('overnode', si >= 0);
        this.requestDraw();
      }
      this._emit('hover', si, e);
      return;
    }

    if (d.mode === 'pan') {
      this.panBy(sx - d.sx, sy - d.sy);
      d.sx = sx; d.sy = sy;
    } else if (d.mode === 'box') {
      this.boxRect = {
        x: Math.min(d.sx, sx), y: Math.min(d.sy, sy),
        w: Math.abs(sx - d.sx), h: Math.abs(sy - d.sy),
      };
      this.requestDraw();
    } else if (d.mode === 'node') {
      const [wx, wy] = this.screenToWorld(sx, sy);
      const dx = wx - d.wx, dy = wy - d.wy;
      if (dx || dy) {
        d.moved = true;
        this.graph.translateSegments([...this.selected], dx, dy);
        d.wx = wx; d.wy = wy;
        this.requestDraw();
      }
    }
  }

  _onUp(e) {
    const d = this._drag;
    this._drag = null;
    this.canvas.classList.remove('panning');
    if (this._pointerId !== null) {
      try { this.canvas.releasePointerCapture(this._pointerId); } catch { /* fine */ }
      this._pointerId = null;
    }
    if (!d) return;

    if (d.mode === 'box') {
      const r = this.boxRect;
      this.boxRect = null;
      if (r && (r.w > 2 || r.h > 2)) {
        const hits = this.boxSelect(r);
        if (e && e.shiftKey) for (const h of hits) this.selected.add(h);
        else this.selected = new Set(hits);
        this._emit('selection');
      }
      this.requestDraw();
    } else if (d.mode === 'node' && d.moved) {
      this.graph.updateBounds();
      this._emit('dragend', [...this.selected]);
      this.requestDraw();
    }
  }

  _emit(name, ...args) {
    const h = this.handlers && this.handlers[name];
    if (typeof h === 'function') {
      try { h(...args); } catch (err) { console.error('renderer handler', name, err); }
    }
  }

  /* -------------------------------------------------------- PNG export */

  /**
   * Render the current view into a PNG blob at `mult`x the on-screen size.
   * Resolves to null if the browser refuses to encode.
   */
  toPNG(mult = 2) {
    const W = this.width, H = this.height;
    const c = document.createElement('canvas');
    c.width = Math.max(1, Math.round(W * mult));
    c.height = Math.max(1, Math.round(H * mult));
    const ctx = c.getContext('2d');
    if (!ctx) return Promise.resolve(null);
    this.paint(ctx, W, H, mult, { interactive: false });
    return new Promise((resolve) => {
      if (c.toBlob) c.toBlob((b) => resolve(b), 'image/png');
      else resolve(null);
    });
  }

  /* -------------------------------------------------------- SVG export */

  /**
   * Serialise the current view as a standalone SVG document. Uses exactly the
   * same culling and geometry as `paint()`, so the figure matches the screen.
   */
  toSVG(o = {}) {
    const g = this.graph;
    const th = this.theme;
    const W = this.width, H = this.height;
    const scale = this.view.scale;
    const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
      .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    const num = (v) => (Math.round(v * 100) / 100);

    const parts = [];
    parts.push(`<?xml version="1.0" encoding="UTF-8"?>`);
    parts.push(`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`);
    parts.push(`<title>Plastr graph view</title>`);
    parts.push(`<rect x="0" y="0" width="${W}" height="${H}" fill="${esc(th.canvasBg)}"/>`);

    if (g.isEmpty) { parts.push('</svg>'); return parts.join('\n'); }
    const scene = this._scene(W, H);

    // links
    if (this.opts.showLinks && scene.vlinks.length) {
      const lw = Math.max(0.5, Math.min(3, 1.1 * Math.sqrt(scale)));
      const d = [];
      for (const li of scene.vlinks) {
        const l = g.links[li];
        d.push(`M${num(this.worldToScreenX(g.px[l.pa]))} ${num(this.worldToScreenY(g.py[l.pa]))}`
             + `L${num(this.worldToScreenX(g.px[l.pb]))} ${num(this.worldToScreenY(g.py[l.pb]))}`);
      }
      parts.push(`<path d="${d.join('')}" fill="none" stroke="${esc(th.link)}" stroke-width="${num(lw)}" stroke-opacity="0.55"/>`);
    }

    // selection halos
    if (this.selected.size) {
      parts.push(`<g fill="none" stroke="${esc(th.sel)}" stroke-opacity="0.5" stroke-linecap="round" stroke-linejoin="round">`);
      for (const si of scene.vis) {
        if (!this.selected.has(si)) continue;
        const lw = Math.max(3, Math.min(90, this.segWidth[si] * scale + 6));
        parts.push(`<path d="${this._svgPathData(si, num)}" stroke-width="${num(lw)}"/>`);
      }
      parts.push('</g>');
    }

    // segments
    const pal = this.colour.palette;
    const nPal = pal.length;
    parts.push('<g fill="none" stroke-linecap="round" stroke-linejoin="round">');
    for (const si of scene.vis) {
      const col = pal[this.colour.segPal[si] % nPal] || th.node;
      const lw = Math.max(0.7, Math.min(90, this.segWidth[si] * scale));
      const seg = g.segments[si];
      parts.push(`<path d="${this._svgPathData(si, num)}" stroke="${esc(col)}" stroke-width="${num(lw)}">`
        + `<title>${esc(seg.name)} · ${esc(fmtBp(seg.length))}</title></path>`);
    }
    parts.push('</g>');

    // arrows
    if (this.opts.showArrows && scene.vis.length <= 4000 && scale > 0.08) {
      parts.push(`<g fill="${esc(th.dim)}">`);
      for (const si of scene.vis) {
        const seg = g.segments[si];
        if (seg.drawLen * scale < 22) continue;
        const pEnd = seg.p0 + seg.k - 1;
        const pPrev = seg.p0 + Math.max(0, seg.k - 2);
        const x2 = this.worldToScreenX(g.px[pEnd]);
        const y2 = this.worldToScreenY(g.py[pEnd]);
        const x1 = this.worldToScreenX(g.px[pPrev]);
        const y1 = this.worldToScreenY(g.py[pPrev]);
        let dx = x2 - x1, dy = y2 - y1;
        const dd = Math.hypot(dx, dy) || 1;
        dx /= dd; dy /= dd;
        const s2 = Math.max(3.5, Math.min(11, this.segWidth[si] * scale * 0.9 + 2.5));
        parts.push(`<polygon points="${num(x2 + dx * s2)},${num(y2 + dy * s2)} `
          + `${num(x2 - dy * s2 * 0.55)},${num(y2 + dx * s2 * 0.55)} `
          + `${num(x2 + dy * s2 * 0.55)},${num(y2 - dx * s2 * 0.55)}"/>`);
      }
      parts.push('</g>');
    }

    // labels
    if (this.opts.showLabels && scene.vis.length <= 600) {
      parts.push(`<g font-family="system-ui, sans-serif" font-size="11" text-anchor="middle" fill="${esc(th.text)}">`);
      for (const si of scene.vis) {
        const seg = g.segments[si];
        if (seg.drawLen * scale < 46) continue;
        const mid = seg.p0 + (seg.k >> 1);
        const x = this.worldToScreenX(g.px[mid]);
        const y = this.worldToScreenY(g.py[mid]) - Math.max(8, this.segWidth[si] * scale * 0.6 + 7);
        parts.push(`<text x="${num(x)}" y="${num(y)}" stroke="${esc(th.canvasBg)}" stroke-width="3" `
          + `paint-order="stroke" >${esc(seg.name)}</text>`);
      }
      parts.push('</g>');
    }

    if (o.caption) {
      parts.push(`<text x="10" y="${H - 10}" font-family="system-ui, sans-serif" font-size="11" `
        + `fill="${esc(th.dim)}">${esc(o.caption)}</text>`);
    }
    parts.push('</svg>');
    return parts.join('\n');
  }

  _svgPathData(si, num) {
    const pts = this._segPoints(si);
    if (!pts.length) return '';
    if (pts.length < 3) {
      return pts.map((p, i) => (i ? 'L' : 'M') + num(p[0]) + ' ' + num(p[1])).join(' ');
    }
    const d = ['M' + num(pts[0][0]) + ' ' + num(pts[0][1])];
    for (let i = 1; i < pts.length - 1; i++) {
      const mx = (pts[i][0] + pts[i + 1][0]) / 2;
      const my = (pts[i][1] + pts[i + 1][1]) / 2;
      d.push('Q' + num(pts[i][0]) + ' ' + num(pts[i][1]) + ' ' + num(mx) + ' ' + num(my));
    }
    const last = pts[pts.length - 1];
    d.push('L' + num(last[0]) + ' ' + num(last[1]));
    return d.join(' ');
  }
}

export default Renderer;
