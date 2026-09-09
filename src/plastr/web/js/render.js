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
 * * Anything placed *on* a node — a label, an alignment sub-span, a path run —
 *   is placed by arc length the way Bandage places it, never by vertex index:
 *   the layout lets particles bunch where a contig bends.
 * * `paint()` is target-agnostic: the live canvas, the PNG export canvas and the
 *   SVG serialiser all consume the same scene description, so what you export is
 *   exactly what you see.
 */

import { fmtBp, fmtNum, perBaseScaleFor, END_START, END_END } from './graph.js';

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
/**
 * Size of the quantised hue wheel behind the 'random' colour mode.
 *
 * A palette entry per segment gives every node its own batch key, so the
 * painter's (colour, width) batching collapses to one `stroke()` per node on
 * exactly the graphs it exists to rescue. 256 hues sit 1.4 degrees apart —
 * closer than two contigs can be told apart anyway — and put a ceiling on the
 * batch count instead.
 */
const HUE_STEPS = 256;
const HUE_PALETTE = [];
for (let i = 0; i < HUE_STEPS; i++) {
  HUE_PALETTE.push(`hsl(${Math.round((i * 360) / HUE_STEPS)} 58% 62%)`);
}
/**
 * Number of quantised line-width buckets.
 *
 * Width carries read depth, so coarse buckets throw that signal away: at 12 a
 * 4x-depth contig and a mean-depth one land in the same bucket at most working
 * zooms. Bucketing exists to keep the batch count down, and only the
 * `strokeStyle` change is expensive, so this can be generous.
 */
const WIDTH_BUCKETS = 48;

/**
 * How far an edge's control point continues past the node end, in world units.
 * Bandage's `edgeLength`, which is also the ideal spring length of a link, so
 * the curve's reach matches the gap the layout leaves.
 */
const EDGE_EXTENSION_WORLD = 5;

/**
 * Catmull-Rom tension, as the fraction of the neighbour span that becomes the
 * Bezier control arm. 1/6 is the textbook value that makes the spline pass
 * through its points with a continuous tangent; a little under keeps a sharp
 * turn from overshooting into a loop.
 */
const SPLINE_TENSION = 0.155;

/**
 * Edge width in world units, so it scales with the drawing like a contig does.
 * A contig at mean depth is five units wide, so an edge is a little under a
 * third of that -- clearly a connector, never mistakable for a contig.
 */
const EDGE_WIDTH_WORLD = 1.5;

/** Widest a node body should get, in CSS pixels, when a fit zooms in. */
const MAX_NODE_BODY_PX = 44;

/**
 * How fast label text grows with zoom. Bandage draws text in scene coordinates
 * scaled by 1 / (1 + (zoom - 1) * factor) (`drawTextPathAtLocation`,
 * graphicsitemnode.cpp:263-295), so on screen it grows as
 * zoom / (1 + (zoom - 1) * factor): nearly constant, but not quite. Text that
 * scaled with the drawing is a wall of letters at any working zoom, and text
 * pinned to an exact pixel size stops relating to the node it names at all.
 */
const TEXT_ZOOM_FACTOR = 0.7;
const LABEL_BASE_PX = 11;
const LABEL_MIN_PX = 8;
const LABEL_MAX_PX = 26;

/** Rainbow parts per query (`blastRainbowPartsPerQuery`, settings.cpp:50). */
const RAINBOW_PARTS_PER_QUERY = 100;

/**
 * Samples used to approximate a link's Bezier when hit testing it. A link is a
 * few pixels wide, so a dozen chords sit well inside the tolerance the pointer
 * needs anyway and solving the cubic would buy nothing.
 */
const LINK_SAMPLES = 16;

function hexToRgb(h) {
  const s = h.replace('#', '');
  const v = s.length === 3
    ? [s[0] + s[0], s[1] + s[1], s[2] + s[2]]
    : [s.slice(0, 2), s.slice(2, 4), s.slice(4, 6)];
  return [parseInt(v[0], 16) || 0, parseInt(v[1], 16) || 0, parseInt(v[2], 16) || 0];
}

function rgbStr(c) { return 'rgb(' + (c[0] | 0) + ',' + (c[1] | 0) + ',' + (c[2] | 0) + ')'; }

/**
 * HSV to a CSS colour, for the rainbow hit parts (`QColor::setHsvF`).
 *
 * Not the same sweep as CSS `hsl()`: at full saturation hsl passes through
 * muddy mid-tones where hsv stays vivid, and a rainbow whose whole job is to be
 * read as a position cannot afford a dull stretch in the middle.
 */
function hsvColour(h, s, v) {
  const i = Math.floor(((h % 1) + 1) % 1 * 6);
  const f = ((h % 1) + 1) % 1 * 6 - i;
  const p = v * (1 - s), q = v * (1 - f * s), t = v * (1 - (1 - f) * s);
  const rgb = [[v, t, p], [q, v, p], [p, v, t], [p, q, v], [t, p, v], [v, p, q]][i % 6];
  return rgbStr([rgb[0] * 255, rgb[1] * 255, rgb[2] * 255]);
}

/** Where an alignment identity sits on the 70-100% scale the legend prints. */
function identityFraction(identity) {
  const id = Number(identity);
  return Number.isFinite(id) ? Math.max(0, Math.min(1, (id - 0.7) / 0.3)) : 1;
}

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

/** FNV-1a over a key, so a colour drawn from it survives a reload. */
function hashKey(key) {
  const s = String(key);
  let h = 2166136261;
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    h = Math.imul(h, 16777619);
  }
  return Math.abs(h);
}

/** Deterministic pseudo-random colour for an arbitrary key. */
export function hashColour(key, sat = 62, light = 58) {
  return `hsl(${hashKey(key) % 360} ${sat}% ${light}%)`;
}

/**
 * Value at a fractional index in a sorted array, interpolating between the two
 * neighbours (`AssemblyGraph::getValueUsingFractionalIndex`).
 */
function fractionalIndex(sorted, index) {
  const n = sorted.length;
  if (!n) return 0;
  if (n === 1) return sorted[0];
  const whole = Math.floor(index);
  if (whole < 0) return sorted[0];
  if (whole >= n - 1) return sorted[n - 1];
  const f = index - whole;
  return sorted[whole] * (1 - f) + sorted[whole + 1] * f;
}

/** First and third quartiles of the drawn segments' depth, or null. */
function depthQuartiles(graph) {
  const ds = [];
  for (const s of graph.segments) {
    if (s.depth !== null && Number.isFinite(s.depth)) ds.push(s.depth);
  }
  if (!ds.length) return null;
  ds.sort((a, b) => a - b);
  return [
    fractionalIndex(ds, (ds.length - 1) / 4),
    fractionalIndex(ds, ((ds.length - 1) * 3) / 4),
  ];
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
    /** Manual [low, high] depth cutoffs; null takes Bandage's auto quartiles. */
    this.depthRange = null;
    /** Reference name -> palette index, so a hit sub-span matches its node. */
    this.refIndex = new Map();
    /** The identity ramp behind 'search', kept for per-hit sub-spans. */
    this.searchRamp = [];
  }

  /** Palette colour for a reference sequence; the unaligned grey if unknown. */
  refColour(ref) {
    const i = this.refIndex.get(String(ref));
    return this.palette[i === undefined ? 0 : i] || this.palette[0];
  }

  /** Heat colour for an alignment identity, matching the 'search' legend. */
  identityColour(identity) {
    const ramp = this.searchRamp;
    if (!ramp.length) return this.palette[this.palette.length - 1];
    return ramp[Math.round(identityFraction(identity) * (ramp.length - 1))];
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
      case 'depth': {
        // Bandage's automatic cutoffs are the first and third quartiles of node
        // depth rather than the extremes (`autoDepthValue`,
        // graphicsitemnode.cpp:851-866). Against min/max a single 500x repeat
        // squeezes every ordinary node into the bottom of the ramp and the
        // colouring stops saying anything; against the quartiles the ramp
        // spends itself on the bulk of the graph and the outliers clamp.
        let [lo, hi] = this.depthRange || depthQuartiles(graph)
          || [graph.stats.depth.min, graph.stats.depth.max];
        if (!(hi > lo)) { lo = graph.stats.depth.min; hi = graph.stats.depth.max; }
        return this._continuous(graph, theme, {
          ramp: VIRIDIS,
          label: 'depth (x)',
          log: true,
          value: (s) => (s.depth === null ? null : s.depth),
          lo,
          hi,
          available: graph.stats.depth.has,
          missingMsg: 'no depth in this graph',
          format: (v) => (v >= 100 ? v.toFixed(0) : v.toFixed(1)) + '×',
          note: this.depthRange ? '' : 'auto range: 1st to 3rd quartile; outside it colours clamp',
        });
      }

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
        // follow an individual contig through a tangle. The hue is quantised
        // (see HUE_STEPS) so the painter can still batch.
        this.palette = HUE_PALETTE;
        for (const s of graph.segments) this.segPal[s.idx] = hashKey('seg:' + s.name) % HUE_STEPS;
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
        this.refIndex = idx;
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
        this.searchRamp = ramp;
        const hits = graph.searchHits;
        let nHit = 0;
        for (const s of graph.segments) {
          const h = hits.get(s.name);
          if (!h) { this.segPal[s.idx] = 0; continue; }
          nHit++;
          this.segPal[s.idx] = 1 + Math.round(identityFraction(h.identity) * (ramp.length - 1));
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
      note: cfg.note || '',
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

/**
 * A pen that records into an SVG path string instead of a canvas.
 *
 * The partial-path and arrowhead builders are geometry, not painting, and the
 * export is only trustworthy while it shares them: a second implementation of
 * "the stretch of this node between two fractions" is a second thing to drift.
 */
class SvgPen {
  constructor(num) { this.d = []; this.num = num; }
  moveTo(x, y) { this.d.push('M' + this.num(x) + ' ' + this.num(y)); }
  lineTo(x, y) { this.d.push('L' + this.num(x) + ' ' + this.num(y)); }
  closePath() { this.d.push('Z'); }
  toString() { return this.d.join(''); }
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
      // Bandage stacks up to four lines on a node and each line is its own
      // toggle (`getNodeText`, graphicsitemnode.cpp:820-840). `showLabels` is
      // the master switch over the four.
      showLabels: false,
      labelName: false,
      labelLength: false,
      labelDepth: false,
      labelCustom: false,
      labelHalo: true,
      labelTextSize: LABEL_BASE_PX,
      showGrid: false,
      // Paint each alignment as a sub-span at its own place along the node
      // rather than flat-colouring the whole node by its best hit.
      showHitSpans: true,
      rainbowHits: false,
      showScaleBar: true,
      depthWidth: true,
      widthScale: 1,
      // World units, constant like Bandage's node width (`averageNodeWidth`).
      // Paired with the absolute length scale this makes a short contig a stub
      // and a long one a ribbon, at the same drawn size in any graph.
      baseWidth: 5,
      // Shape of the depth-to-width curve: `depthPower` is how fast width
      // follows depth and `depthEffectOnWidth` how much of that reaches the
      // drawing at all. Bandage's defaults (settings.cpp:39-41).
      depthPower: 0.5,
      depthEffectOnWidth: 0.5,
    };

    this.selected = new Set();   // segment indices
    this.selectedLinks = new Set();
    this.hoverIdx = -1;
    this.hoverLink = -1;
    this._hoverPt = null;        // last cursor position, for the link readout
    this.boxRect = null;         // screen-space rectangle while box selecting

    /** Highlighted path: segment index -> [startFraction, endFraction]. */
    this.pathRuns = new Map();
    this.pathLinks = new Set();
    this.pathName = '';

    /**
     * Every search hit, grouped by segment. The model keeps only the best hit
     * per segment, which is enough to colour a node and not enough to say
     * *where* on the node the hit is; see `setSearchHits`.
     */
    this.searchHitsBySeg = new Map();
    this._queryExtent = new Map();
    this._cumBuf = new Float64Array(64);

    this.segWidth = new Float32Array(0);
    this.segBucket = new Int32Array(0);
    this.bucketWidth = new Float32Array(WIDTH_BUCKETS);
    this._arrowBuf = new Float64Array(6);

    this._buckets = [];
    this._usedKeys = [];
    this._visible = [];
    this._visibleLinks = [];
    this._unselected = [];
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

  /**
   * The most zoomed-in a fit should ever go.
   *
   * Framing one short contig would otherwise solve for whatever scale makes a
   * five-unit stub fill the window, which paints a single ribbon a hundred
   * pixels thick and tells the user nothing. Cap it where the widest node body
   * is still a ribbon rather than a slab, so framing a small selection centres
   * it at a readable size instead of swallowing the screen.
   */
  _maxUsefulScale() {
    let widest = 0;
    for (let i = 0; i < this.segWidth.length; i++) {
      if (this.segWidth[i] > widest) widest = this.segWidth[i];
    }
    if (!(widest > 0)) return 20;
    return Math.min(20, Math.max(1, MAX_NODE_BODY_PX / widest));
  }

  /** Fit the given segment indices (or the whole graph) into the viewport. */
  fitToView(indices, padding = 0.08) {
    const g = this.graph;
    if (g.isEmpty) { this.view = { cx: 0, cy: 0, scale: 1 }; this.requestDraw(); return; }
    const [minX, minY, maxX, maxY] = g.boundsOf(indices && indices.length ? indices : null);
    const w = Math.max(1, maxX - minX);
    const h = Math.max(1, maxY - minY);
    const s = Math.min(this.width / (w * (1 + padding * 2)), this.height / (h * (1 + padding * 2)));
    this.view.scale = Math.max(1e-4, Math.min(this._maxUsefulScale(), s));
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

    // Bandage's width curve, exactly: depth relative to the *mean depth of the
    // drawn nodes*, raised to `depthPower`, damped towards 1 by
    // `depthEffectOnWidth` (getNodeWidth, graphicsitemnode.cpp:889-895 and
    // 917-924). The damping is what keeps the drawing readable: at the default
    // half-power, half-effect a 100x repeat is five times a mean node's width
    // rather than a hundred times, and no clamp is needed to rescue it.
    let mean = 0;
    if (this.opts.depthWidth && g.stats.depth.has) {
      // Weighted by length, as Bandage's mean drawn depth is (getMeanDepth,
      // assemblygraph.cpp:241-263): a thousand short misassembled fragments
      // must not outvote the megabase contig they were broken off, or every
      // real node comes out over-wide.
      let sum = 0, bases = 0, plain = 0, count = 0;
      for (const s of g.segments) {
        if (s.depth === null || !Number.isFinite(s.depth)) continue;
        const len = s.length > 0 ? s.length : 0;
        sum += s.depth * len;
        bases += len;
        plain += s.depth;
        count++;
      }
      // A graph that declares no lengths still has depths worth scaling by.
      mean = bases > 0 ? sum / bases : (count ? plain / count : 0);
    }
    const base = this.opts.baseWidth * this.opts.widthScale;
    const power = this.opts.depthPower;
    const effect = this.opts.depthEffectOnWidth;
    let minW = Infinity, maxW = -Infinity;
    for (const s of g.segments) {
      // Bandage's fallback when there is no mean to divide by.
      let rel = 1;
      if (mean > 0 && s.depth !== null && Number.isFinite(s.depth)) {
        rel = Math.max(0, s.depth / mean);
      }
      const w = Math.max(0, base * ((Math.pow(rel, power) - 1) * effect + 1));
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
    if (key === 'depthWidth' || key === 'widthScale' || key === 'baseWidth'
        || key === 'depthPower' || key === 'depthEffectOnWidth') this.updateStyle();
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

  selectNames(names, add = false) {
    const g = this.graph;
    const set = add ? new Set(this.selected) : new Set();
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
    if (!this.selected.size && !this.selectedLinks.size) return;
    this.selected = new Set();
    this.selectedLinks = new Set();
    this.requestDraw();
    this._emit('selection');
  }

  /** Select links by index into `graph.links`. */
  setLinkSelection(indices) {
    this.selectedLinks = new Set(indices);
    this.requestDraw();
    this._emit('selection');
  }

  /** The selected links themselves, for a panel that wants to describe them. */
  selectedLinkList() {
    const out = [];
    for (const i of this.selectedLinks) {
      const l = this.graph.links[i];
      if (l) out.push(l);
    }
    return out;
  }

  /**
   * Highlight a resolved path as a run along its member segments.
   *
   * Deliberately tolerant about what a path is, because more than one thing in
   * the model is one: the resolved record `_resolvePaths` builds
   * (`{name, segs, links}`, already matched against the drawn graph), a GFA `P`
   * line (`{name, steps: ['3+', '7-']}`), a walk with per-segment extents
   * (`{name, segments: [{name, start, end}]}`, fractions or base coordinates),
   * a bare list of segment names, or the name of a path in the model.
   * Anything it cannot resolve to a drawn segment is skipped rather than
   * refused, so a path that runs off the end of a truncated graph still draws
   * the part that is on screen.
   */
  setHighlightPath(path) {
    this.pathRuns = new Map();
    this.pathLinks = new Set();
    this.pathName = '';
    let p = path;
    if (typeof p === 'string') {
      p = (this.graph.pathByName ? this.graph.pathByName(path) : null)
        || (this.graph.paths || []).find((q) => q && String(q.name) === path) || null;
    }
    if (!p) { this.requestDraw(); return this; }

    // A resolved path already carries segment indices and the links between
    // them, so there is nothing to look up and nothing to guess.
    if (p.segs && p.segs.length !== undefined && typeof p.segs[0] === 'number') {
      this.pathName = String(p.name || '');
      for (const si of p.segs) if (this.graph.segments[si]) this.pathRuns.set(si, [0, 1]);
      for (const li of (p.linkSet || p.links || [])) if (li >= 0) this.pathLinks.add(li);
      this.requestDraw();
      return this;
    }

    const steps = Array.isArray(p) ? p : (p.segments || p.steps || p.nodes || []);
    if (!Array.isArray(p)) this.pathName = String(p.name || '');

    const order = [];
    for (const step of steps) {
      const raw = typeof step === 'string' ? step : (step && (step.name ?? step.segment));
      if (raw === undefined || raw === null) continue;
      const name = String(raw);
      // GFA writes the orientation onto the step ('3+'); a resolved path need
      // not, and a segment may legitimately be called '3+'.
      let seg = this.graph.segmentByName(name);
      if (!seg && /[+-]$/.test(name)) seg = this.graph.segmentByName(name.slice(0, -1));
      if (!seg) continue;
      let s = 0, e = 1;
      if (step && typeof step === 'object') {
        const a = Number(step.start ?? step.startFraction);
        const b = Number(step.end ?? step.endFraction);
        if (Number.isFinite(a) && Number.isFinite(b) && b > a) {
          // Fractions and base coordinates are both natural to write; only
          // bases can exceed 1.
          const denom = (a > 1 || b > 1) ? (seg.length || 1) : 1;
          s = Math.max(0, a / denom);
          e = Math.min(1, b / denom);
        }
      }
      this.pathRuns.set(seg.idx, [s, e]);
      order.push(seg.idx);
    }

    // The joins matter as much as the members: a path drawn as disconnected
    // runs does not read as a route through the graph.
    if (order.length > 1) {
      const want = new Set();
      for (let i = 0; i + 1 < order.length; i++) {
        want.add(order[i] + ':' + order[i + 1]);
        want.add(order[i + 1] + ':' + order[i]);
      }
      for (const l of this.graph.links) {
        if (want.has(l.a + ':' + l.b)) this.pathLinks.add(l.idx);
      }
    }
    this.requestDraw();
    return this;
  }

  clearHighlightPath() { this.setHighlightPath(null); }

  /**
   * Take the full search hit list (`POST /api/search` -> `hits`) so hits can be
   * drawn as sub-spans. Falls back to the model's best-hit-per-segment map when
   * nothing calls this, which still places one span per segment.
   */
  setSearchHits(hits) {
    this.searchHitsBySeg = new Map();
    this._queryExtent = new Map();
    for (const h of Array.isArray(hits) ? hits : []) {
      if (!h || h.segment === undefined || h.segment === null) continue;
      const key = String(h.segment);
      const arr = this.searchHitsBySeg.get(key);
      if (arr) arr.push(h); else this.searchHitsBySeg.set(key, [h]);
      // A search hit carries no query length, so the query's own aligned extent
      // stands in as the rainbow's denominator: exact when the query aligns to
      // its ends, and a fixed stretch of the same rainbow when it does not.
      const q = String(h.query);
      if ((Number(h.q_en) || 0) > (this._queryExtent.get(q) || 0)) {
        this._queryExtent.set(q, Number(h.q_en) || 0);
      }
    }
    this.requestDraw();
    return this;
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
   *
   * `trimPx` stops the stroke short of the final vertex, leaving room for the
   * arrowhead wedge (`_arrowPath`) to finish the node.
   */

  /**
   * The vertex just inside each neighbouring contig, for every contig end.
   *
   * A turn between two contigs used to happen entirely in the short link
   * between them, which put a corner in a drawing that should read as one
   * flowing strand. Knowing what is on the far side of each end lets the spline
   * enter and leave a contig already aimed at its neighbours, so the turn is
   * spread along the contigs instead of concentrated at the join.
   *
   * Cached per graph revision: this walks the adjacency, and paint() runs on
   * every frame.
   */
  _neighbourVertices() {
    const g = this.graph;
    if (this._ghostRev === g.revision && this._ghostStart) {
      return { start: this._ghostStart, end: this._ghostEnd };
    }
    const n = g.segments.length;
    const start = new Int32Array(n).fill(-1);
    const end = new Int32Array(n).fill(-1);
    for (const seg of g.segments) {
      for (const e of g.adj[seg.idx] || []) {
        if (e.seg === seg.idx) continue; // a self-link has no far side
        const other = g.segments[e.seg];
        if (!other) continue;
        // One step into the neighbour, so the direction is the neighbour's, not
        // the near-zero step across the gap.
        const ghost = e.otherEnd === END_END
          ? other.p0 + Math.max(0, other.k - 2)
          : other.p0 + Math.min(other.k - 1, 1);
        const slot = e.selfEnd === END_START ? start : end;
        if (slot[seg.idx] < 0) slot[seg.idx] = ghost;
      }
    }
    this._ghostRev = g.revision;
    this._ghostStart = start;
    this._ghostEnd = end;
    return { start, end };
  }

  _segPath(ctx, si, smooth, trimPx = 0) {
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
    // A closed contig is a ring: its vertices wrap, so the spline runs all the
    // way round and the path is closed rather than stopping one vertex short
    // and leaving a notch at the seam. It has no ends, so nothing is trimmed
    // for an arrowhead either.
    if (seg.closed && k >= 3) {
      const X = (i) => px[p0 + (i % k)] * s + ox;
      const Y = (i) => py[p0 + (i % k)] * s + oy;
      ctx.moveTo((X(0) + X(1)) / 2, (Y(0) + Y(1)) / 2);
      for (let i = 1; i <= k; i++) {
        ctx.quadraticCurveTo(X(i), Y(i), (X(i) + X(i + 1)) / 2, (Y(i) + Y(i + 1)) / 2);
      }
      ctx.closePath();
      return;
    }
    let ex = px[p0 + k - 1] * s + ox, ey = py[p0 + k - 1] * s + oy;
    if (trimPx > 0) {
      const bx = px[p0 + k - 2] * s + ox, by = py[p0 + k - 2] * s + oy;
      const dx = ex - bx, dy = ey - by;
      const d = Math.hypot(dx, dy);
      // Clamped to the final leg: a node shorter than the wedge collapses to
      // nothing here and is drawn as the wedge alone, which is Bandage's
      // degenerate pure-triangle node.
      if (d > 0) {
        const t = Math.min(trimPx, d) / d;
        ex -= dx * t; ey -= dy * t;
      }
    }
    if (!smooth) {
      ctx.moveTo(px[p0] * s + ox, py[p0] * s + oy);
      for (let i = 1; i < k - 1; i++) ctx.lineTo(px[p0 + i] * s + ox, py[p0 + i] * s + oy);
      ctx.lineTo(ex, ey);
      return;
    }

    // Catmull-Rom through the contig's own vertices, with a control point taken
    // from the contig on the far side of each end. A contig has three vertices
    // on average, so a spline over its own points alone is barely a curve and
    // the drawing is really a chain of straight bars meeting at angles -- the
    // thing that makes a graph look shaky rather than fluid. Borrowing the
    // neighbours' vertices bends each contig towards what it joins, and the
    // turn is spread over two contig bodies instead of one short link.
    const ghosts = this._neighbourVertices();
    const gs = ghosts.start[si];
    const ge = ghosts.end[si];
    const X = (i) => px[p0 + i] * s + ox;
    const Y = (i) => py[p0 + i] * s + oy;

    // Point before the first vertex, and after the last.
    let bx, by;
    if (gs >= 0) { bx = px[gs] * s + ox; by = py[gs] * s + oy; }
    else { bx = 2 * X(0) - X(1); by = 2 * Y(0) - Y(1); }   // mirror: straight out
    let axx, ayy;
    if (ge >= 0 && trimPx <= 0) { axx = px[ge] * s + ox; ayy = py[ge] * s + oy; }
    else { axx = 2 * ex - X(k - 2); ayy = 2 * ey - Y(k - 2); }

    const pt = (i) => {
      if (i < 0) return [bx, by];
      if (i > k - 1) return [axx, ayy];
      if (i === k - 1) return [ex, ey];
      return [X(i), Y(i)];
    };

    ctx.moveTo(X(0), Y(0));
    for (let i = 0; i < k - 1; i++) {
      const [x0, y0] = pt(i - 1);
      const [x1, y1] = pt(i);
      const [x2, y2] = pt(i + 1);
      const [x3, y3] = pt(i + 2);
      ctx.bezierCurveTo(
        x1 + (x2 - x0) * SPLINE_TENSION, y1 + (y2 - y0) * SPLINE_TENSION,
        x2 - (x3 - x1) * SPLINE_TENSION, y2 - (y3 - y1) * SPLINE_TENSION,
        x2, y2,
      );
    }
  }

  /**
   * Append one link to the current path, Bandage-style.
   *
   * A link is a cubic Bézier whose control points continue each node's own
   * terminal direction (`GraphicsItemEdge::calculateAndSetPath`). That
   * tangential join is what makes connected contigs read as one flowing strand
   * instead of a set of bars wired together, and it is the most recognisable
   * property of a Bandage drawing. The extension is clamped to half the gap so
   * a short link cannot loop back on itself.
   *
   * A link from a segment to itself is a circular contig closing up. Drawn as a
   * chord it would run straight through the node body and be hidden by it, so
   * it is bowed out sideways instead.
   */
  _linkPath(ctx, l, scale) {
    const q = this._linkGeom(l, scale);
    if (q.skip) return;
    ctx.moveTo(q.ax, q.ay);
    if (q.straight) ctx.lineTo(q.bx, q.by);
    else ctx.bezierCurveTo(q.c1x, q.c1y, q.c2x, q.c2y, q.bx, q.by);
  }

  /**
   * Drawn width of a link. `edgeWidth` is the user's setting; the square root
   * of the zoom keeps a link from swelling into a ribbon of its own when you
   * zoom into a tangle.
   */
  _linkWidth(scale) {
    // Edges are a first-class part of the picture, not a hint that two contigs
    // are related: they are what the eye follows from one contig to the next,
    // so they scale with the view like the contigs do rather than being pinned
    // to a hairline. The old cap of 3 px meant that at any zoom past the first
    // the edges thinned away to nothing while the contigs grew.
    const w = Number(this.opts.edgeWidth) > 0 ? Number(this.opts.edgeWidth) : EDGE_WIDTH_WORLD;
    return Math.max(0.8, Math.min(14, w * scale));
  }

  /**
   * Screen-space geometry of one link. Shared by the canvas painter and the SVG
   * writer so an export cannot drift from what is on screen.
   */
  _linkGeom(l, scale) {
    const g = this.graph;
    const ax = this.worldToScreenX(g.px[l.pa]);
    const ay = this.worldToScreenY(g.py[l.pa]);
    const bx = this.worldToScreenX(g.px[l.pb]);
    const by = this.worldToScreenY(g.py[l.pb]);

    // Inward neighbours give each end its direction. Older models may not carry
    // them, in which case fall back to a straight chord.
    if (l.paIn === undefined || l.pbIn === undefined) {
      return { ax, ay, bx, by, straight: true };
    }
    const aix = this.worldToScreenX(g.px[l.paIn]);
    const aiy = this.worldToScreenY(g.py[l.paIn]);
    const bix = this.worldToScreenX(g.px[l.pbIn]);
    const biy = this.worldToScreenY(g.py[l.pbIn]);

    // Outward unit vectors: away from the inner neighbour, through the end.
    let adx = ax - aix, ady = ay - aiy;
    const adl = Math.hypot(adx, ady) || 1;
    adx /= adl; ady /= adl;
    let bdx = bx - bix, bdy = by - biy;
    const bdl = Math.hypot(bdx, bdy) || 1;
    bdx /= bdl; bdy /= bdl;

    const ext = EDGE_EXTENSION_WORLD * scale;

    // A closed contig is already drawn as a ring, so its own self-link needs no
    // bow: the loop the reader sees IS the link. Bowing it as well draws a
    // second arc across a circle that is already closed.
    if (l.selfLoop && g.segments[l.a] && g.segments[l.a].closed) {
      return { ax, ay, bx, by, straight: true, skip: true };
    }
    if (l.selfLoop) {
      // Bow the loop off the node along the end's normal so it stays visible.
      const nx = -ady, ny = adx;
      const r = Math.max(ext * 3, 10);
      return {
        ax, ay, bx, by, straight: false,
        c1x: ax + adx * r + nx * r, c1y: ay + ady * r + ny * r,
        c2x: bx + bdx * r + nx * r, c2y: by + bdy * r + ny * r,
      };
    }

    const dist = Math.hypot(bx - ax, by - ay);
    const e = Math.min(ext, dist / 2);
    return {
      ax, ay, bx, by, straight: false,
      c1x: ax + adx * e, c1y: ay + ady * e,
      c2x: bx + bdx * e, c2y: by + bdy * e,
    };
  }

  /**
   * Painted screen width of a segment's ribbon.
   *
   * This is the *bucketed* width the painter actually strokes, floor and
   * ceiling included, so hit testing, halos and the SVG export all agree with
   * the pixels on screen rather than with the unquantised ideal.
   */
  _widthPx(si) {
    return Math.max(1, Math.min(90, this.bucketWidth[this.segBucket[si]] * this.view.scale));
  }

  /**
   * Screen-space wedge that turns a node's flat end into a point: apex on the
   * final vertex, base half a width back, so the sides run at 45 degrees. That
   * is exactly the notch Bandage subtracts from the node body
   * (GraphicsItemNode::shape, graphicsitemnode.cpp:443-467) — the arrowhead is
   * the end of the node, not an ornament parked beside it.
   *
   * Fills `out` with [tipX, tipY, leftX, leftY, rightX, rightY].
   */
  _arrowPoints(si, bodyPx, out) {
    const g = this.graph;
    const seg = g.segments[si];
    if (seg.k < 2) return false;
    const pEnd = seg.p0 + seg.k - 1;
    const x2 = this.worldToScreenX(g.px[pEnd]);
    const y2 = this.worldToScreenY(g.py[pEnd]);
    let dx = x2 - this.worldToScreenX(g.px[pEnd - 1]);
    let dy = y2 - this.worldToScreenY(g.py[pEnd - 1]);
    const d = Math.hypot(dx, dy);
    if (!d) return false;
    dx /= d; dy /= d;
    const half = bodyPx / 2;
    // A node shorter than half its own width cannot give the wedge its full
    // 45 degrees, so the base falls back to the previous vertex and the node
    // becomes a plain triangle, as it does in Bandage.
    const back = Math.min(half, d);
    const bx = x2 - dx * back, by = y2 - dy * back;
    out[0] = x2; out[1] = y2;
    out[2] = bx - dy * half; out[3] = by + dx * half;
    out[4] = bx + dy * half; out[5] = by - dx * half;
    return true;
  }

  /** Append one arrowhead wedge to the current path. */
  _arrowPath(ctx, si, bodyPx) {
    const a = this._arrowBuf;
    if (!this._arrowPoints(si, bodyPx, a)) return;
    ctx.moveTo(a[0], a[1]);
    ctx.lineTo(a[2], a[3]);
    ctx.lineTo(a[4], a[5]);
    ctx.closePath();
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

  /* ----------------------------------------------------- arc length */

  /**
   * Fill `_cumBuf` with the cumulative world-space length along a segment's
   * polyline and return the total.
   *
   * Everything Bandage places on a node — its labels, its hit spans, a path
   * run — is placed by arc length, never by vertex index (`getNodePathLength`,
   * graphicsitemnode.cpp:648). The two are not the same: the layout lets
   * particles bunch where a contig bends, so the middle vertex of a curved node
   * can sit a long way from its middle.
   */
  _arc(si) {
    const g = this.graph;
    const seg = g.segments[si];
    if (this._cumBuf.length < seg.k) this._cumBuf = new Float64Array(seg.k * 2);
    const cum = this._cumBuf;
    cum[0] = 0;
    for (let i = 1; i < seg.k; i++) {
      const a = seg.p0 + i;
      cum[i] = cum[i - 1] + Math.hypot(g.px[a] - g.px[a - 1], g.py[a] - g.py[a - 1]);
    }
    return cum[seg.k - 1];
  }

  /** Screen point at arc-length fraction `f` (`findLocationOnPath`). */
  _pointAtFraction(si, f) {
    const g = this.graph;
    const seg = g.segments[si];
    if (seg.k < 2) {
      return [this.worldToScreenX(g.px[seg.p0]), this.worldToScreenY(g.py[seg.p0])];
    }
    const total = this._arc(si);
    const cum = this._cumBuf;
    const target = Math.max(0, Math.min(1, f)) * total;
    let i = 1;
    while (i < seg.k - 1 && cum[i] < target) i++;
    const d = cum[i] - cum[i - 1];
    const t = d > 0 ? Math.max(0, Math.min(1, (target - cum[i - 1]) / d)) : 0;
    const a = seg.p0 + i - 1;
    return [
      this.worldToScreenX(g.px[a] + (g.px[a + 1] - g.px[a]) * t),
      this.worldToScreenY(g.py[a] + (g.py[a + 1] - g.py[a]) * t),
    ];
  }

  /**
   * Append the stretch of a segment between two arc-length fractions to the
   * current path (`makePartialPath`, graphicsitemnode.cpp:600). Straight
   * between vertices rather than splined: a sub-span is stroked at the node's
   * own width directly over the body, so it has to follow the same line the
   * body's spline is already close to, and a second spline over a partial
   * vertex range would drift off it at the cut ends.
   */
  _partialPath(ctx, si, f0, f1) {
    if (f1 < f0) { const t = f0; f0 = f1; f1 = t; }
    const g = this.graph;
    const seg = g.segments[si];
    // A highlighted path outlives the graph it was resolved against: redrawing
    // with a narrower scope leaves it holding indices that no longer exist.
    if (!seg || seg.k < 2 || !(f1 > f0)) return false;
    const total = this._arc(si);
    const cum = this._cumBuf;
    const a = Math.max(0, f0) * total;
    const b = Math.min(1, f1) * total;
    let started = false;
    for (let i = 1; i < seg.k; i++) {
      if (cum[i] < a) continue;
      const p = seg.p0 + i - 1;
      const dx = g.px[p + 1] - g.px[p], dy = g.py[p + 1] - g.py[p];
      const d = cum[i] - cum[i - 1];
      if (!started) {
        started = true;
        const t = d > 0 ? (a - cum[i - 1]) / d : 0;
        ctx.moveTo(this.worldToScreenX(g.px[p] + dx * t), this.worldToScreenY(g.py[p] + dy * t));
      }
      if (cum[i] >= b) {
        const t = d > 0 ? Math.max(0, (b - cum[i - 1]) / d) : 1;
        ctx.lineTo(this.worldToScreenX(g.px[p] + dx * t), this.worldToScreenY(g.py[p] + dy * t));
        return true;
      }
      ctx.lineTo(this.worldToScreenX(g.px[p + 1]), this.worldToScreenY(g.py[p + 1]));
    }
    return started;
  }

  /* ------------------------------------------------------- hit spans */

  /** The palette colour a segment's body is painted in. */
  _segColour(si) {
    const pal = this.colour.palette;
    return pal[this.colour.segPal[si] % pal.length] || this.theme.node;
  }

  /**
   * The sub-spans an alignment claims of one node, in paint order, or null.
   *
   * An alignment is a property of a *stretch* of a contig, not of the contig
   * (`getBlastHitParts`, blasthit.cpp:49-101). Flat-colouring a node by its best
   * hit hides that it is half chromosome and half plasmid, which is the first
   * thing anyone opens a graph viewer to see. `lenPx` is the node's drawn
   * length, used only to stop the rainbow cutting parts thinner than a pixel.
   */
  _hitSpans(seg, lenPx) {
    if (!this.opts.showHitSpans) return null;
    const rainbow = !!this.opts.rainbowHits;
    const out = [];

    if (this.colour.mode === 'reference') {
      const hits = seg.refHits;
      if (!hits || !hits.length) return null;
      for (const h of hits) {
        const qlen = Number(h.q_len) > 0 ? Number(h.q_len) : seg.length;
        const s = Number(h.q_st) / qlen, e = Number(h.q_en) / qlen;
        if (!(qlen > 0) || !Number.isFinite(s) || !Number.isFinite(e) || e <= s) continue;
        // Bandage's rainbow runs along the sequence you searched *with*. Here
        // the node is the query and the chromosome is the subject, so the
        // rainbow runs along the reference and the colours on a node say which
        // part of the chromosome each stretch of it came from. On the minus
        // strand the node walks the reference backwards.
        const rlen = Number(h.r_len);
        const q0 = Number(h.r_st) / rlen, q1 = Number(h.r_en) / rlen;
        if (rainbow && rlen > 0 && Number.isFinite(q0) && Number.isFinite(q1)) {
          const back = Number(h.strand) === -1;
          this._rainbowParts(out, s, e, back ? q1 : q0, back ? q0 : q1, lenPx);
        } else {
          out.push({ colour: this.colour.refColour(h.ref), s, e });
        }
      }
      return out.length ? out : null;
    }

    if (this.colour.mode === 'search') {
      let hits = this.searchHitsBySeg.get(seg.name);
      if (!hits) {
        const best = this.graph.searchHits.get(seg.name);
        if (!best) return null;
        hits = [best];
      }
      const len = seg.length;
      if (!(len > 0)) return null;
      for (const h of hits) {
        const s = Number(h.s_st) / len, e = Number(h.s_en) / len;
        if (!Number.isFinite(s) || !Number.isFinite(e) || e <= s) continue;
        const qlen = this._queryExtent.get(String(h.query)) || 0;
        const q0 = Number(h.q_st) / qlen, q1 = Number(h.q_en) / qlen;
        if (rainbow && qlen > 0 && Number.isFinite(q0) && Number.isFinite(q1)) {
          const back = Number(h.strand) === -1;
          this._rainbowParts(out, s, e, back ? q1 : q0, back ? q0 : q1, lenPx);
        } else {
          out.push({ colour: this.colour.identityColour(h.identity), s, e });
        }
      }
      return out.length ? out : null;
    }
    return null;
  }

  /**
   * Cut one hit into parts coloured by position along the query
   * (blasthit.cpp:53-88). The 0.9 keeps the far end off red so the start and
   * the end of a query can never be confused for each other; the part count is
   * capped against the drawn length because parts under a pixel buy nothing and
   * cost a stroke each.
   */
  _rainbowParts(out, s, e, qFrom, qTo, lenPx) {
    let parts = Math.ceil(RAINBOW_PARTS_PER_QUERY * Math.abs(qTo - qFrom));
    const room = Math.floor((e - s) * lenPx * 2);
    if (parts > room) parts = room;
    if (parts < 1) parts = 1;
    const ns = (e - s) / parts;
    const qs = (qTo - qFrom) / parts;
    for (let i = 0; i < parts; i++) {
      const q = Math.max(0, Math.min(1, qFrom + qs * i));
      out.push({ colour: hsvColour(q * 0.9, 1, 1), s: s + ns * i, e: s + ns * (i + 1) });
    }
  }

  /** Stroke one segment's body, and its arrowhead, in a given colour. */
  _strokeSegment(ctx, si, colour, widthPx, arrows) {
    const tip = arrows && widthPx >= 2.5 ? widthPx / 2 : 0;
    ctx.strokeStyle = colour;
    ctx.lineWidth = widthPx;
    ctx.beginPath();
    this._segPath(ctx, si, this.graph.segments[si].drawLen * this.view.scale > 30, tip);
    ctx.stroke();
    if (!tip) return;
    ctx.fillStyle = colour;
    ctx.beginPath();
    this._arrowPath(ctx, si, widthPx);
    ctx.fill();
  }

  /* ---------------------------------------------------------- labels */

  /** The stack of lines to write on a node (`getNodeText`). */
  _labelLines(seg) {
    const o = this.opts;
    const lines = [];
    // A custom label is whatever the model has attached to the segment; it goes
    // first, as it does in Bandage.
    const custom = seg.label ?? seg.customLabel;
    if (o.labelCustom && custom) lines.push(String(custom));
    // `showLabels` used to mean the name on its own, and a stored session or a
    // host that knows only the master switch still means that by it.
    if (o.labelName || !(o.labelLength || o.labelDepth || (o.labelCustom && custom))) {
      lines.push(seg.name);
    }
    if (o.labelLength) lines.push(fmtBp(seg.length));
    if (o.labelDepth && seg.depth !== null) lines.push(fmtNum(seg.depth, 1) + '×');
    return lines;
  }

  /** Label size in screen pixels; see TEXT_ZOOM_FACTOR. */
  _labelPx() {
    const z = Math.max(0.05, this.view.scale);
    const base = Number(this.opts.labelTextSize) > 0 ? Number(this.opts.labelTextSize) : LABEL_BASE_PX;
    const px = (base * z) / (1 + (z - 1) * TEXT_ZOOM_FACTOR);
    return Math.max(LABEL_MIN_PX, Math.min(LABEL_MAX_PX, px));
  }

  _labelFont(px) {
    return `${px.toFixed(1)}px system-ui, -apple-system, "Segoe UI", Roboto, sans-serif`;
  }

  /* -------------------------------------------------------- scale bar */

  /**
   * Metrics for the scale bar: a round number of bases and the screen distance
   * they occupy.
   *
   * Bandage has no scale bar because until the draw length was calibrated there
   * was nothing to measure. Now a node is `perBase` world units per base, so a
   * ruler on the canvas is honest — for every node above the `minLen` floor.
   * A stub shorter than that is drawn longer than it is, which is a floor worth
   * having and the one thing this bar does not describe.
   */
  _scaleBarMetrics(W) {
    const g = this.graph;
    const perBase = (g.perBase === null || g.perBase === undefined
      ? perBaseScaleFor(g.geom) : g.perBase) * (g.geom.lengthScale || 1);
    const pxPerBase = perBase * this.view.scale;
    if (!Number.isFinite(pxPerBase) || pxPerBase <= 0) return null;
    // A 1-2-5 step, so the bar always reads as a round number of bases.
    const target = 150 / pxPerBase;
    const mag = Math.pow(10, Math.floor(Math.log10(target)));
    let bp = mag;
    for (const m of [2, 5, 10]) if (mag * m <= target) bp = mag * m;
    const len = bp * pxPerBase;
    if (!(len >= 24) || len > W * 0.6) return null;
    return { bp, len };
  }

  _paintScaleBar(ctx, W, H) {
    const m = this._scaleBarMetrics(W);
    if (!m) return;
    const th = this.theme;
    const x = 16, y = H - 20;
    ctx.strokeStyle = th.dim;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(x + 0.5, y - 4); ctx.lineTo(x + 0.5, y + 4);
    ctx.moveTo(x + 0.5, y); ctx.lineTo(x + m.len + 0.5, y);
    ctx.moveTo(x + m.len + 0.5, y - 4); ctx.lineTo(x + m.len + 0.5, y + 4);
    ctx.stroke();
    ctx.fillStyle = th.dim;
    ctx.font = '11px system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'alphabetic';
    ctx.fillText(fmtBp(m.bp), x + m.len / 2, y - 7);
  }

  /* ------------------------------------------------------- link paint */

  /**
   * Links that are part of something: a highlighted path, the selection, or
   * whatever the pointer is over. Drawn over the plain pass, because an edge
   * you have picked out is worth nothing if it stays buried in the mesh.
   */
  _paintLinkAccents(ctx, scene, scale, interactive) {
    const g = this.graph;
    const th = this.theme;
    const base = this._linkWidth(scale);
    const groups = [
      [this.pathLinks, th.accent, base * 3 + 1, 0.85],
      [this.selectedLinks, th.sel, base * 2.4 + 1, 0.95],
    ];
    for (const [set, colour, w, alpha] of groups) {
      if (!set.size) continue;
      ctx.strokeStyle = colour;
      ctx.lineWidth = w;
      ctx.globalAlpha = alpha;
      ctx.beginPath();
      for (const li of scene.vlinks) if (set.has(li)) this._linkPath(ctx, g.links[li], scale);
      ctx.stroke();
    }
    if (interactive && this.hoverLink >= 0 && this.hoverLink < g.links.length) {
      ctx.strokeStyle = th.accent;
      ctx.lineWidth = base * 2.6 + 1.5;
      ctx.globalAlpha = 1;
      ctx.beginPath();
      this._linkPath(ctx, g.links[this.hoverLink], scale);
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }

  /** One line describing a link, overlap included. */
  linkDescription(l) {
    if (!l) return '';
    return `${l.from}${l.fromOrient} → ${l.to}${l.toOrient}`
      + (l.overlap ? ` · ${fmtBp(l.overlap)} overlap` : ' · blunt join')
      + (l.selfLoop ? ' · self loop' : '');
  }

  /**
   * A readout for the hovered link, drawn on the canvas.
   *
   * Only when nothing has claimed `linkhover`: a host that puts the link in its
   * own tooltip or status strip should own the whole job, and two readouts for
   * one edge is worse than none.
   */
  _paintLinkReadout(ctx, W, H) {
    if (this.handlers && typeof this.handlers.linkhover === 'function') return;
    const l = this.graph.links[this.hoverLink];
    const pt = this._hoverPt;
    if (!l || !pt) return;
    const th = this.theme;
    const text = this.linkDescription(l);
    ctx.font = '11px system-ui, -apple-system, "Segoe UI", Roboto, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    const w = ctx.measureText(text).width + 16;
    const h = 22;
    const x = Math.max(4, Math.min(pt[0] + 14, W - w - 4));
    const y = Math.max(4, Math.min(pt[1] + 14, H - h - 4));
    ctx.fillStyle = th.canvasBg;
    ctx.globalAlpha = 0.92;
    ctx.fillRect(x, y, w, h);
    ctx.globalAlpha = 1;
    ctx.strokeStyle = th.outline || th.faint;
    ctx.lineWidth = 1;
    ctx.strokeRect(x + 0.5, y + 0.5, w - 1, h - 1);
    ctx.fillStyle = th.text;
    ctx.fillText(text, x + 8, y + h / 2);
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
    // Below this the wedge would be a pixel or two and only make the node look
    // frayed, so the arrowheads (and with them the flat cap they need) are left
    // off entirely.
    const arrows = this.opts.showArrows && scale > 0.08;

    if (this.opts.showGrid) this._paintGrid(ctx, W, H);

    /* ---- links ---- */
    if (this.opts.showLinks && scene.vlinks.length) {
      ctx.strokeStyle = th.link;
      ctx.lineCap = 'round';
      ctx.lineJoin = 'round';
      ctx.lineWidth = this._linkWidth(scale);
      ctx.beginPath();
      for (const li of scene.vlinks) {
        if (this.pathLinks.has(li) || this.selectedLinks.has(li)) continue;
        this._linkPath(ctx, g.links[li], scale);
      }
      ctx.stroke();
      this._paintLinkAccents(ctx, scene, scale, o.interactive !== false);
    }

    /* ---- highlighted path: the rim goes under the bodies ---- */
    this._paintPathRuns(ctx, scene, arrows, true);

    /* ---- segments, batched by (colour, width bucket) ---- */
    // Selected nodes are held back and drawn last, so the selection is never
    // buried under the tangle it was made in (Bandage raises a selected item to
    // the front, graphicsitemnode.cpp:180-195).
    const unsel = this._unselected;
    unsel.length = 0;
    for (const si of scene.vis) if (!this.selected.has(si)) unsel.push(si);

    const pal = this.colour.palette;
    const nPal = pal.length;
    const keys = this._usedKeys;
    keys.length = 0;
    const buckets = this._buckets;
    const needed = nPal * WIDTH_BUCKETS;
    while (buckets.length < needed) buckets.push([]);
    for (const si of unsel) {
      const key = (this.colour.segPal[si] % nPal) * WIDTH_BUCKETS + this.segBucket[si];
      const arr = buckets[key];
      if (arr.length === 0) keys.push(key);
      arr.push(si);
    }
    // Bandage strokes the node body with a flat cap and then cuts the
    // arrowhead out of it, so an arrowed node is still one solid shape.
    ctx.lineCap = 'butt';
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
        // Keep a floor in *screen* pixels so ribbons stay visible when zoomed
        // out, but a soft one: at 2.4 px every node hit the floor together at
        // ordinary zooms and the depth signal in the width vanished, while the
        // outline pass below could never trigger because `body` was never
        // allowed under its own threshold.
        const body = Math.max(1, Math.min(90, this.bucketWidth[wIdx] * scale));
        const rim = Math.min(2.6, Math.max(1, body * 0.28));
        const tip = arrows && body >= 2.5 ? body / 2 : 0;
        if (pass === 0) {
          // Skip the outline when the body is too thin for it to read.
          if (body < 2.2) continue;
          ctx.strokeStyle = th.outline || 'rgba(20,26,34,0.85)';
          ctx.lineWidth = body + rim;
        } else {
          ctx.strokeStyle = pal[colIdx] || th.node;
          ctx.lineWidth = body;
        }
        ctx.beginPath();
        for (const si of arr) {
          this._segPath(ctx, si, g.segments[si].drawLen * scale > 30, tip);
        }
        ctx.stroke();
        if (!tip) continue;
        // The wedge is a filled shape rather than a stroke, so it needs a path
        // of its own; it rides in the same batch to keep the style changes down.
        ctx.beginPath();
        for (const si of arr) this._arrowPath(ctx, si, body);
        if (pass === 0) {
          ctx.lineWidth = rim;
          ctx.stroke();
        } else {
          ctx.fillStyle = ctx.strokeStyle;
          ctx.fill();
        }
      }
    }
    for (const key of keys) buckets[key].length = 0;

    /* ---- alignment sub-spans ---- */
    this._paintHitSpans(ctx, unsel, scale, arrows);

    /* ---- selected nodes, raised to the front ---- */
    if (this.selected.size) {
      ctx.lineCap = 'butt';
      ctx.lineJoin = 'round';
      const sel = [];
      for (const si of scene.vis) if (this.selected.has(si)) sel.push(si);
      ctx.globalAlpha = 0.55;
      for (const si of sel) {
        ctx.strokeStyle = th.sel;
        ctx.lineWidth = this._widthPx(si) + 6;
        ctx.beginPath();
        this._segPath(ctx, si, g.segments[si].drawLen * scale > 30);
        ctx.stroke();
      }
      ctx.globalAlpha = 1;
      // The node keeps its own colour inside the halo: a selection that
      // recolours what it selects hides the very property you selected it for.
      for (const si of sel) this._strokeSegment(ctx, si, this._segColour(si), this._widthPx(si), arrows);
      this._paintHitSpans(ctx, sel, scale, arrows);
    }

    /* ---- highlighted path: the shading goes over them ---- */
    this._paintPathRuns(ctx, scene, arrows, false);

    /* ---- hover: brighten, never replace ---- */
    if (o.interactive && this.hoverIdx >= 0 && this.hoverIdx < g.segments.length) {
      const hw = this._widthPx(this.hoverIdx);
      this._strokeSegment(ctx, this.hoverIdx, this._segColour(this.hoverIdx), hw, arrows);
      this._paintHitSpans(ctx, [this.hoverIdx], scale, arrows);
      // A wash of white over the node's own colour, rather than the accent in
      // place of it: on a rainbow graph a recoloured node reads as a different
      // node, and the colouring the user turned on is exactly what they are
      // pointing at it to read.
      ctx.globalAlpha = 0.28;
      this._strokeSegment(ctx, this.hoverIdx, '#ffffff', hw, arrows);
      ctx.globalAlpha = 1;
    }

    /* ---- labels ---- */
    this._paintLabels(ctx, scene);

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

    if (o.interactive && this.hoverLink >= 0) this._paintLinkReadout(ctx, W, H);
    if (this.opts.showScaleBar) this._paintScaleBar(ctx, W, H);
  }

  /**
   * Paint every alignment sub-span on a list of segments, over the node bodies
   * that are already down. Bandage clips the parts to the node outline so they
   * cannot spill past an arrowhead; here the body is a stroke rather than a
   * filled outline, so a span that reaches the tip is trimmed the same way the
   * body is and the wedge is refilled in the span's own colour.
   */
  _paintHitSpans(ctx, list, scale, arrows) {
    if (!this.opts.showHitSpans) return;
    const g = this.graph;
    const mode = this.colour.mode;
    if (mode !== 'reference' && mode !== 'search') return;
    ctx.lineCap = 'butt';
    ctx.lineJoin = 'round';
    for (const si of list) {
      const seg = g.segments[si];
      const lenPx = seg.drawLen * scale;
      const parts = this._hitSpans(seg, lenPx);
      if (!parts) continue;
      const body = this._widthPx(si);
      const tip = arrows && body >= 2.5 ? body / 2 : 0;
      const cut = tip && lenPx > 0 ? 1 - Math.min(tip, lenPx) / lenPx : 1;
      ctx.lineWidth = body;
      for (const p of parts) {
        ctx.strokeStyle = p.colour;
        // A span that claims the whole node has to claim its end caps too, or
        // the body colour underneath shows as a nub past either end.
        ctx.lineCap = 'butt';
        if (p.s < cut) {
          ctx.beginPath();
          if (this._partialPath(ctx, si, p.s, Math.min(p.e, cut))) ctx.stroke();
        }
        if (tip && p.e >= cut - 1e-9) {
          ctx.fillStyle = p.colour;
          ctx.beginPath();
          this._arrowPath(ctx, si, body);
          ctx.fill();
        }
      }
    }
    ctx.lineCap = 'butt';
  }

  /**
   * The run a highlighted path claims of each member node: outlined around it
   * and shaded along it (`pathHighlightNode3`, graphicsitemnode.cpp:1090-1103).
   *
   * Two passes, because both have to leave the node itself readable: the
   * outline goes *under* the bodies so it shows as a rim rather than swallowing
   * the node, and the shading goes over them at low alpha so the node's own
   * colour — and any alignment span painted on it — still reads through. A path
   * is a route through the drawing, not a recolouring of it.
   */
  _paintPathRuns(ctx, scene, arrows, under) {
    if (!this.pathRuns.size) return;
    const g = this.graph;
    const r = scene.rect;
    ctx.lineCap = 'butt';
    ctx.lineJoin = 'round';
    ctx.strokeStyle = under ? this.theme.accent : '#ffffff';
    ctx.globalAlpha = under ? 0.95 : 0.2;
    for (const [si, run] of this.pathRuns) {
      const o = si * 4;
      if (g.bbox[o + 2] < r[0] || g.bbox[o] > r[2] || g.bbox[o + 3] < r[1] || g.bbox[o + 1] > r[3]) continue;
      const body = this._widthPx(si);
      ctx.lineWidth = under ? body + 6 : body;
      ctx.beginPath();
      if (this._partialPath(ctx, si, run[0], run[1])) ctx.stroke();
    }
    ctx.globalAlpha = 1;
    ctx.lineCap = 'butt';
  }

  /**
   * Node labels: the stack of enabled lines, centred on the node's arc-length
   * centre and counter-scaled so the text stays near a constant size while the
   * drawing zooms (graphicsitemnode.cpp:209-234).
   */
  _paintLabels(ctx, scene) {
    const g = this.graph;
    const th = this.theme;
    if (!this.opts.showLabels || scene.vis.length > 600) return;
    const px = this._labelPx();
    const halo = this.opts.labelHalo !== false;
    ctx.font = this._labelFont(px);
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.lineJoin = 'round';
    // Bandage strokes the outline at twice its thickness and fills over it, so
    // the halo is a rim around the glyph rather than a fattened glyph.
    ctx.lineWidth = Math.max(2, px * 0.34);
    ctx.strokeStyle = halo ? 'rgba(255,255,255,0.92)' : th.canvasBg;
    ctx.fillStyle = halo ? '#10141b' : th.text;
    const lh = px * 1.15;
    for (const si of scene.vis) {
      const seg = g.segments[si];
      const lines = this._labelLines(seg);
      if (!lines.length) continue;
      // Nothing readable fits on a node a few pixels long, and a graph's worth
      // of overlapping stacks is worse than no labels at all.
      if (seg.drawLen * this.view.scale < 30) continue;
      const c = this._pointAtFraction(si, 0.5);
      const top = c[1] - ((lines.length - 1) * lh) / 2;
      for (let i = 0; i < lines.length; i++) {
        ctx.strokeText(lines[i], c[0], top + i * lh);
        ctx.fillText(lines[i], c[0], top + i * lh);
      }
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
      // The painted width, not the ideal one: the ribbon on screen is the
      // bucketed width with a pixel floor, and the clickable region has to be
      // the shape the user can actually see.
      const halfW = this._widthPx(si) / (2 * this.view.scale);
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

  /**
   * Link index under a screen point, or -1.
   *
   * The Bezier is sampled rather than solved: a link is a couple of pixels
   * wide, so the chord error over a sixteenth of a curve is far inside the
   * tolerance the pointer needs anyway. Every point of a cubic lies inside its
   * control hull, so the hull rejects most links before any sampling happens.
   */
  hitTestLink(sx, sy, tolPx = 6) {
    if (!this.opts.showLinks || this.graph.isEmpty) return -1;
    const g = this.graph;
    const scale = this.view.scale;
    let best = -1, bestD = Infinity;
    for (const li of this._visibleLinks) {
      const q = this._linkGeom(g.links[li], scale);
      const x1 = q.straight ? q.ax : Math.min(q.ax, q.c1x, q.c2x, q.bx);
      const x2 = q.straight ? q.bx : Math.max(q.ax, q.c1x, q.c2x, q.bx);
      const y1 = q.straight ? q.ay : Math.min(q.ay, q.c1y, q.c2y, q.by);
      const y2 = q.straight ? q.by : Math.max(q.ay, q.c1y, q.c2y, q.by);
      if (sx < Math.min(x1, x2) - tolPx || sx > Math.max(x1, x2) + tolPx ||
          sy < Math.min(y1, y2) - tolPx || sy > Math.max(y1, y2) + tolPx) continue;
      let d;
      if (q.straight) {
        d = distPtSeg(sx, sy, q.ax, q.ay, q.bx, q.by);
      } else {
        d = Infinity;
        let px = q.ax, py = q.ay;
        for (let i = 1; i <= LINK_SAMPLES; i++) {
          const t = i / LINK_SAMPLES;
          const u = 1 - t;
          const nx = u * u * u * q.ax + 3 * u * u * t * q.c1x + 3 * u * t * t * q.c2x + t * t * t * q.bx;
          const ny = u * u * u * q.ay + 3 * u * u * t * q.c1y + 3 * u * t * t * q.c2y + t * t * t * q.by;
          const dd = distPtSeg(sx, sy, px, py, nx, ny);
          if (dd < d) d = dd;
          px = nx; py = ny;
        }
      }
      if (d < tolPx && d < bestD) { bestD = d; best = li; }
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
      if (this.hoverIdx !== -1 || this.hoverLink !== -1) {
        this.hoverIdx = -1;
        this.hoverLink = -1;
        this.requestDraw();
      }
      this._hoverPt = null;
      this._emit('hover', -1, null, null);
      this._emit('linkhover', null, null);
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

  /**
   * What the pointer is over, passed as a third argument to `hover` so a host
   * can describe a link without a handler of its own.
   */
  hoverInfo() {
    if (this.hoverIdx >= 0) {
      return { kind: 'segment', segment: this.graph.segments[this.hoverIdx], link: null };
    }
    if (this.hoverLink >= 0) {
      const link = this.graph.links[this.hoverLink];
      return { kind: 'link', segment: null, link, text: this.linkDescription(link) };
    }
    return null;
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
    const li = si < 0 ? this.hitTestLink(sx, sy) : -1;
    this._pointerId = e.pointerId;
    try { this.canvas.setPointerCapture(e.pointerId); } catch { /* not fatal */ }
    this.canvas.focus({ preventScroll: true });

    const wantPan = e.button === 1 || this._spaceDown
      || (e.button === 0 && si < 0 && li < 0 && !e.shiftKey);
    const wantBox = e.button === 0 && si < 0 && li < 0 && e.shiftKey;

    if (e.button === 0 && si < 0 && li >= 0 && !this._spaceDown) {
      if (e.shiftKey) {
        if (this.selectedLinks.has(li)) this.selectedLinks.delete(li);
        else this.selectedLinks.add(li);
      } else {
        this.selected = new Set();
        this.selectedLinks = new Set([li]);
      }
      // Still a pan drag underneath, so grabbing the canvas near an edge and
      // pulling does what it does everywhere else.
      this._drag = { mode: 'pan', sx, sy };
      this._emit('selection');
      this.requestDraw();
      return;
    }

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
      // A node wins over a link it ends on: the node is the thing the pointer
      // is plainly over, and the link's own reach starts where the node stops.
      const li = si >= 0 ? -1 : this.hitTestLink(sx, sy);
      const moved = !this._hoverPt || this._hoverPt[0] !== sx || this._hoverPt[1] !== sy;
      this._hoverPt = [sx, sy];
      if (si !== this.hoverIdx || li !== this.hoverLink || (li >= 0 && moved)) {
        this.hoverIdx = si;
        this.hoverLink = li;
        this.canvas.classList.toggle('overnode', si >= 0 || li >= 0);
        this.requestDraw();
      }
      this._emit('hover', si, e, this.hoverInfo());
      this._emit('linkhover', li >= 0 ? this.graph.links[li] : null, e);
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
      const lw = this._linkWidth(scale);
      const d = [];
      for (const li of scene.vlinks) {
        const q = this._linkGeom(g.links[li], scale);
        if (q.skip) continue;
        d.push(`M${num(q.ax)} ${num(q.ay)}` + (q.straight
          ? `L${num(q.bx)} ${num(q.by)}`
          : `C${num(q.c1x)} ${num(q.c1y)} ${num(q.c2x)} ${num(q.c2y)} ${num(q.bx)} ${num(q.by)}`));
      }
      parts.push(`<path d="${d.join('')}" fill="none" stroke="${esc(th.link)}" stroke-width="${num(lw)}" stroke-linecap="round"/>`);
    }

    // links that belong to a path or the selection are drawn over the rest
    if (this.opts.showLinks && (this.pathLinks.size || this.selectedLinks.size)) {
      const lw = this._linkWidth(scale);
      for (const [set, colour, w] of [[this.pathLinks, th.accent, lw * 3 + 1],
        [this.selectedLinks, th.sel, lw * 2.4 + 1]]) {
        const d = [];
        for (const li of scene.vlinks) {
          if (!set.has(li)) continue;
          const q = this._linkGeom(g.links[li], scale);
          d.push(`M${num(q.ax)} ${num(q.ay)}` + (q.straight
            ? `L${num(q.bx)} ${num(q.by)}`
            : `C${num(q.c1x)} ${num(q.c1y)} ${num(q.c2x)} ${num(q.c2y)} ${num(q.bx)} ${num(q.by)}`));
        }
        if (d.length) {
          parts.push(`<path d="${d.join('')}" fill="none" stroke="${esc(colour)}" stroke-width="${num(w)}"/>`);
        }
      }
    }

    // the path rim, under the bodies
    if (this.pathRuns.size) {
      parts.push(`<g fill="none" stroke="${esc(th.accent)}" stroke-opacity="0.95" stroke-linecap="butt">`);
      for (const [si, run] of this.pathRuns) {
        const pen = new SvgPen(num);
        if (this._partialPath(pen, si, run[0], run[1])) {
          parts.push(`<path d="${pen}" stroke-width="${num(this._widthPx(si) + 6)}"/>`);
        }
      }
      parts.push('</g>');
    }

    // selection halos
    if (this.selected.size) {
      parts.push(`<g fill="none" stroke="${esc(th.sel)}" stroke-opacity="0.55" stroke-linecap="round" stroke-linejoin="round">`);
      for (const si of scene.vis) {
        if (!this.selected.has(si)) continue;
        parts.push(`<path d="${this._svgPathData(si, num)}" stroke-width="${num(this._widthPx(si) + 6)}"/>`);
      }
      parts.push('</g>');
    }

    // segments, with the arrowhead wedges the canvas painter draws. Selected
    // nodes come last, in the same raise-to-front order as on screen.
    const arrows = this.opts.showArrows && scale > 0.08;
    const pal = this.colour.palette;
    const nPal = pal.length;
    const wedges = [];
    const pts = this._arrowBuf;
    const order = scene.vis.filter((si) => !this.selected.has(si))
      .concat(scene.vis.filter((si) => this.selected.has(si)));
    parts.push(`<g fill="none" stroke-linecap="${arrows ? 'butt' : 'round'}" stroke-linejoin="round">`);
    for (const si of order) {
      const col = pal[this.colour.segPal[si] % nPal] || th.node;
      // The painted width, so the export matches the screen exactly.
      const lw = this._widthPx(si);
      const tip = arrows && lw >= 2.5 ? lw / 2 : 0;
      const seg = g.segments[si];
      parts.push(`<path d="${this._svgPathData(si, num, tip)}" stroke="${esc(col)}" stroke-width="${num(lw)}">`
        + `<title>${esc(seg.name)} · ${esc(fmtBp(seg.length))}</title></path>`);
      if (tip && this._arrowPoints(si, lw, pts)) {
        wedges.push(`<polygon fill="${esc(col)}" points="${num(pts[0])},${num(pts[1])} `
          + `${num(pts[2])},${num(pts[3])} ${num(pts[4])},${num(pts[5])}"/>`);
      }
      // alignment sub-spans, over the body they belong to
      const lenPx = seg.drawLen * scale;
      const spans = this._hitSpans(seg, lenPx) || [];
      const cut = tip && lenPx > 0 ? 1 - Math.min(tip, lenPx) / lenPx : 1;
      for (const p of spans) {
        if (p.s < cut) {
          const pen = new SvgPen(num);
          if (this._partialPath(pen, si, p.s, Math.min(p.e, cut))) {
            parts.push(`<path d="${pen}" stroke="${esc(p.colour)}" stroke-width="${num(lw)}" stroke-linecap="butt"/>`);
          }
        }
        if (tip && p.e >= cut - 1e-9 && this._arrowPoints(si, lw, pts)) {
          wedges.push(`<polygon fill="${esc(p.colour)}" points="${num(pts[0])},${num(pts[1])} `
            + `${num(pts[2])},${num(pts[3])} ${num(pts[4])},${num(pts[5])}"/>`);
        }
      }
      // and the highlighted path's run along it, shaded over the body the way
      // the canvas painter shades it
      const run = this.pathRuns.get(si);
      if (run) {
        const pen = new SvgPen(num);
        if (this._partialPath(pen, si, run[0], run[1])) {
          parts.push(`<path d="${pen}" stroke="#ffffff" stroke-width="${num(lw)}" `
            + `stroke-opacity="0.2" stroke-linecap="butt"/>`);
        }
      }
    }
    parts.push('</g>');
    if (wedges.length) parts.push('<g>' + wedges.join('') + '</g>');

    // labels: the same stack, centred on the same arc-length centre
    if (this.opts.showLabels && scene.vis.length <= 600) {
      const px = this._labelPx();
      const lh = px * 1.15;
      const halo = this.opts.labelHalo !== false;
      parts.push(`<g font-family="system-ui, sans-serif" font-size="${num(px)}" text-anchor="middle" `
        + `fill="${esc(halo ? '#10141b' : th.text)}" stroke="${esc(halo ? '#ffffff' : th.canvasBg)}" `
        + `stroke-width="${num(Math.max(2, px * 0.34))}" stroke-linejoin="round" paint-order="stroke">`);
      for (const si of scene.vis) {
        const lines = this._labelLines(g.segments[si]);
        if (!lines.length || g.segments[si].drawLen * scale < 30) continue;
        const c = this._pointAtFraction(si, 0.5);
        const top = c[1] - ((lines.length - 1) * lh) / 2 + px * 0.35;
        for (let i = 0; i < lines.length; i++) {
          parts.push(`<text x="${num(c[0])}" y="${num(top + i * lh)}">${esc(lines[i])}</text>`);
        }
      }
      parts.push('</g>');
    }

    if (this.opts.showScaleBar) {
      const m = this._scaleBarMetrics(W);
      if (m) {
        const x = 16, y = H - 20;
        parts.push(`<path d="M${x} ${y - 4}V${y + 4}M${x} ${y}H${num(x + m.len)}M${num(x + m.len)} ${y - 4}V${y + 4}" `
          + `fill="none" stroke="${esc(th.dim)}" stroke-width="1"/>`);
        parts.push(`<text x="${num(x + m.len / 2)}" y="${y - 7}" font-family="system-ui, sans-serif" `
          + `font-size="11" text-anchor="middle" fill="${esc(th.dim)}">${esc(fmtBp(m.bp))}</text>`);
      }
    }

    if (o.caption) {
      parts.push(`<text x="10" y="${H - 10}" font-family="system-ui, sans-serif" font-size="11" `
        + `fill="${esc(th.dim)}">${esc(o.caption)}</text>`);
    }
    parts.push('</svg>');
    return parts.join('\n');
  }

  _svgPathData(si, num, trimPx = 0) {
    const pts = this._segPoints(si);
    if (!pts.length) return '';
    if (trimPx > 0 && pts.length > 1) {
      const a = pts[pts.length - 2], b = pts[pts.length - 1];
      const dx = b[0] - a[0], dy = b[1] - a[1];
      const d = Math.hypot(dx, dy);
      if (d > 0) {
        const t = Math.min(trimPx, d) / d;
        pts[pts.length - 1] = [b[0] - dx * t, b[1] - dy * t];
      }
    }
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
