/**
 * ui.js — panel wiring.
 *
 * Owns the `app` object that the other modules talk to, and connects every
 * control in index.html to the API. Nothing here draws: the renderer owns the
 * canvas and this module only tells it what to show.
 */

import api, { ApiError, API_BASE } from './api.js';
import GraphModel, { fmtBp, fmtInt, fmtNum, bestRefHit } from './graph.js';
import { LayoutController, DEFAULT_PARAMS } from './layout.js';
import Renderer, { readTheme } from './render.js';
import initScaffoldPanel from './scaffold.js';

const $ = (id) => document.getElementById(id);
const on = (id, event, fn) => {
  const el = $(id);
  if (el) el.addEventListener(event, fn);
  return el;
};

function esc(s) {
  return String(s ?? '').replace(/[&<>"']/g, (c) => (
    { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]
  ));
}

const num = (id, fallback = 0) => {
  const v = parseFloat($(id)?.value);
  return Number.isFinite(v) ? v : fallback;
};
const int = (id, fallback = 0) => {
  const v = parseInt($(id)?.value, 10);
  return Number.isFinite(v) ? v : fallback;
};
const checked = (id) => !!$(id)?.checked;
const radio = (name) => document.querySelector(`input[name="${name}"]:checked`)?.value;

/* ------------------------------------------------------------------ the app */

export function createApp() {
  const graph = new GraphModel();
  const canvas = $('canvas');
  const renderer = new Renderer(canvas, graph);

  const app = {
    graph,
    renderer,
    serverStatus: null,
    report: null,
    scaffoldPanel: null,
    layout: null,
    busyDepth: 0,
  };

  /* ---- toasts, busy, errors ---- */

  app.toast = (message, kind = 'info', ms = 4200) => {
    const box = $('toasts');
    if (!box) return;
    const el = document.createElement('div');
    el.className = `toast ${kind}`;
    el.innerHTML = `<span>${esc(message)}</span><button class="toast-x" aria-label="dismiss">×</button>`;
    el.querySelector('.toast-x').addEventListener('click', () => el.remove());
    box.appendChild(el);
    if (ms) setTimeout(() => el.remove(), ms);
  };

  app.showError = (err) => {
    // The server puts a human-readable sentence in `error`; show it verbatim.
    const message = err instanceof ApiError
      ? (err.payload?.error || err.message)
      : (err?.message || String(err));
    app.toast(message, 'error', 9000);
    setStatus(message);
    console.error(err);
  };

  app.setBusy = (busy, label = 'Working…') => {
    app.busyDepth = Math.max(0, app.busyDepth + (busy ? 1 : -1));
    const on_ = app.busyDepth > 0;
    const busyEl = $('busy');
    if (busyEl) busyEl.hidden = !on_;
    const l = $('busy-label');
    if (l && on_) l.textContent = label;
    document.body.classList.toggle('is-busy', on_);
  };

  const setStatus = (msg) => { const el = $('status-msg'); if (el) el.textContent = msg; };
  app.setStatus = setStatus;
  app.status = setStatus; // scaffold.js reports progress through this

  /** Run an async action with the busy indicator and uniform error handling. */
  const guard = async (label, fn) => {
    app.setBusy(true, label);
    try {
      return await fn();
    } catch (err) {
      app.showError(err);
      return null;
    } finally {
      app.setBusy(false);
    }
  };
  app.guard = guard;

  /* ---- status / graph loading ---- */

  app.refreshStatus = async () => {
    const status = await api.status();
    app.serverStatus = status;
    $('file-name').textContent = status.source_path
      ? status.source_path.split('/').pop()
      : 'no assembly loaded';
    $('file-name').title = status.source_path || '';
    const g = status.graph;
    $('graph-summary').textContent = g
      ? `${fmtInt(g.segments)} segments · ${fmtInt(g.links)} links · ${fmtBp(g.total_length)}`
      : '';
    $('undo-depth').textContent = status.undo_depth ? `(${status.undo_depth})` : '';
    $('op-undo').disabled = !status.undo_depth;
    $('ref-drop').disabled = !status.reference;
    if (status.reference) {
      $('ref-path').value = status.reference.path || $('ref-path').value;
    }
    document.body.classList.toggle('has-reference', !!status.reference);
    document.body.classList.toggle('has-graph', !!status.loaded);
    return status;
  };

  /* ---- graph scope ---- */

  /**
   * Query parameters for GET /api/graph, from the Graph scope panel.
   *
   * `scope` and its arguments are Bandage's graph scope: the *server* decides
   * which segments a scope names, because only it holds the sequences and the
   * whole graph. A server that does not know these parameters ignores them and
   * returns the entire graph, so `min_length` and `max_nodes` always go too.
   */
  function scopeParams() {
    const scope = $('scope-mode')?.value || 'entire';
    const params = {
      scope,
      min_length: int('lod-minlen', 0),
      max_nodes: int('lod-maxnodes', 15000),
    };
    if (scope === 'around') {
      params.nodes = ($('scope-nodes')?.value || '').split(/[,\s]+/).filter(Boolean).join(',');
      params.distance = int('scope-distance', 0);
      params.match = radio('scope-match') || 'exact';
    } else if (scope === 'depth') {
      if ($('scope-depth-min')?.value !== '') params.depth_min = num('scope-depth-min', 0);
      if ($('scope-depth-max')?.value !== '') params.depth_max = num('scope-depth-max', 0);
    } else if (scope === 'component') {
      const comp = $('lod-component')?.value;
      if (comp !== '' && comp !== undefined && comp !== null) params.component = parseInt(comp, 10);
    }
    return params;
  }

  /**
   * GET /api/graph with whatever parameters the scope panel produced.
   *
   * `api.graph()` forwards only the three it was written for; rather than wait
   * for it to grow the other five this goes at the endpoint directly, and
   * raises the same ApiError the rest of the app already handles.
   */
  async function fetchGraph(params) {
    const query = new URLSearchParams();
    for (const [k, v] of Object.entries(params)) {
      if (v === null || v === undefined || v === '') continue;
      query.append(k, String(v));
    }
    const res = await fetch(`${API_BASE}/graph?${query}`, { headers: { Accept: 'application/json' } });
    const text = await res.text();
    let data = null;
    try { data = JSON.parse(text); } catch { /* an error page, not JSON */ }
    if (!res.ok) {
      throw new ApiError(data?.error || `HTTP ${res.status}`, data?.detail || text.slice(0, 400),
        res.status, '/graph');
    }
    return data || {};
  }

  /**
   * Say what the server drew and what it left out, and why.
   *
   * `shown`/`total`/`truncated` come from the payload. A server with graph
   * scopes may also send `dropped`: a list of `{reason, count}` naming what each
   * filter removed, which is printed verbatim. Links are counted here instead,
   * because whether a link resolved is only known once the client has matched
   * its ends against the segments it actually received.
   */
  function renderScopeReport(payload) {
    const box = $('scope-report');
    if (!box) return;
    if (!payload) { box.innerHTML = ''; return; }

    const shown = Number(payload.shown ?? graph.segments.length) || 0;
    const total = Number(payload.total ?? shown) || 0;
    const hidden = Math.max(0, total - shown);

    const notes = [];
    const dropped = Array.isArray(payload.dropped) ? payload.dropped : [];
    for (const d of dropped) {
      notes.push(`${fmtInt(d.count)} ${esc(d.reason || 'segment(s) dropped')}`);
    }
    // A reason is only guessed here for a server too old to send one: printed
    // beside the server's own account it would state a second, different rule.
    if (!dropped.length && hidden) {
      notes.push(payload.truncated
        ? `${fmtInt(hidden)} segment(s) past the ${fmtInt(int('lod-maxnodes', 15000))}-node cap.`
        : `${fmtInt(hidden)} segment(s) fall outside this scope.`);
    }
    if (payload.truncated) notes.push('Raise “max nodes drawn” to see more.');
    if (graph.droppedLinks) {
      notes.push(`${fmtInt(graph.droppedLinks)} link(s) hidden because the segment at one end is not drawn.`);
    }

    box.innerHTML =
      `<div class="kv"><span>Drawn</span><b>${fmtInt(shown)} / ${fmtInt(total)} segments</b></div>` +
      `<div class="kv"><span>Links</span><b>${fmtInt(graph.links.length)}</b></div>` +
      (notes.length
        ? `<ul class="drop-list">${notes.map((n) => `<li>${n}</li>`).join('')}</ul>`
        : '<p class="result-note">Nothing was dropped.</p>');
  }

  /** Fetch the graph and hand it to the renderer. */
  app.reloadGraph = async ({ relayout = false, keepPositions = true } = {}) => {
    const payload = await fetchGraph(scopeParams());
    graph.setData(payload, { keepPositions });
    if (!keepPositions) graph.seedPositions();
    graph.updateBounds();
    // The renderer caches a colour index and a line width per segment; those
    // arrays are sized to the old graph until this is called.
    renderer.updateStyle();
    app.layout?.markDirty();

    fillComponentSelects();
    renderer.setColourMode($('colour-mode')?.value || 'random');
    updateColourLegend();
    updateCounts();
    const hint = $('empty-hint');
    if (hint) hint.hidden = !graph.isEmpty;

    if (payload.truncated) {
      app.toast(
        `Showing the ${fmtInt(payload.shown)} longest of ${fmtInt(payload.total)} segments. ` +
        'Raise "max nodes" or filter by length to see more.',
        'warn', 7000,
      );
    }
    renderScopeReport(payload);
    // Fit now so something is on screen while the layout settles, and fit
    // again when it finishes, since the graph moves a long way in between.
    renderer.fitToView(null);
    if (relayout) {
      app._fitAfterLayout = true;
      startLayout({ mode: radio('rr-mode') || 'force', scope: 'all' });
    }
    return payload;
  };

  app.loadFile = async (path, format) => guard('Loading assembly…', async () => {
    await api.load(path, format || null);
    await app.refreshStatus();
    await app.reloadGraph({ relayout: true, keepPositions: false });
    app.scaffoldPanel?.clear();
    clearReferenceResults();
    setStatus(`Loaded ${path}`);
    app.toast('Assembly loaded', 'ok');
  });

  /* ---- counts / legend ---- */

  function updateCounts() {
    const el = $('status-counts');
    if (el) {
      el.textContent = graph.isEmpty
        ? '—'
        : `${fmtInt(graph.segments.length)} shown / ${fmtInt(graph.total)} total`;
    }
    const sel = $('status-sel');
    if (sel) {
      const n = renderer.selectedNames().length;
      sel.textContent = n ? `${fmtInt(n)} selected` : 'nothing selected';
    }
  }
  app.updateCounts = updateCounts;

  /** Render whatever legend the renderer's colour mapper produced. */
  function updateColourLegend() {
    const box = $('colour-legend');
    if (!box) return;
    const legend = renderer.colour?.legend;
    if (!legend) { box.innerHTML = ''; return; }

    const caption = legend.label ? `<p class="muted small">${esc(legend.label)}</p>` : '';
    const note = legend.note ? `<p class="muted small">${esc(legend.note)}</p>` : '';

    if (legend.type === 'scale') {
      const stops = (legend.ramp || []).join(',');
      const ticks = (legend.ticks || []).map((t) => `<span>${esc(t)}</span>`).join('');
      box.innerHTML = caption +
        `<div class="legend-ramp" style="background:linear-gradient(90deg,${stops})"></div>` +
        `<div class="legend-ticks">${ticks}</div>` + note;
      return;
    }
    if (legend.type === 'cat') {
      const rows = (legend.items || []).map((e) =>
        `<div class="legend-row"><i class="legend-sw" style="background:${esc(e.colour)}"></i>` +
        `<span>${esc(e.label)}</span></div>`).join('');
      box.innerHTML = caption + `<div class="legend-block">${rows}</div>` +
        (legend.more ? `<p class="legend-more">+ ${fmtInt(legend.more)} more</p>` : '') + note;
      return;
    }
    box.innerHTML = caption + note;
  }
  app.updateColourLegend = updateColourLegend;

  function fillComponentSelects() {
    const comps = graph.components || [];
    for (const id of ['lod-component', 'rr-component']) {
      const sel = $(id);
      if (!sel) continue;
      const previous = sel.value;
      const options = ['<option value="">All components</option>'];
      for (const c of comps.slice(0, 400)) {
        options.push(
          `<option value="${esc(c.id)}">#${esc(c.id)} — ${fmtInt(c.segs.length)} seg, ${esc(fmtBp(c.length))}</option>`,
        );
      }
      sel.innerHTML = options.join('');
      if (previous && [...sel.options].some((o) => o.value === previous)) sel.value = previous;
    }
  }

  /* ---- find nodes ---- */

  /**
   * Segment names matching a query.
   *
   * Bandage matches either the whole name or any part of it, over a
   * comma-separated list of terms; this is the same, over the segments actually
   * drawn. A name the current scope excluded cannot be found here, which is
   * why the no-match note points back at the scope panel.
   */
  function matchNames(query, exact) {
    const terms = String(query || '').split(/[,\s]+/).filter(Boolean);
    if (!terms.length) return [];
    // Deduplicated, or a name typed twice is counted and listed twice while
    // only ever selecting one segment.
    if (exact) return [...new Set(terms)].filter((t) => graph.byName.has(t));
    const lowered = terms.map((t) => t.toLowerCase());
    const out = [];
    for (const seg of graph.segments) {
      const name = seg.name.toLowerCase();
      if (lowered.some((t) => name.includes(t))) out.push(seg.name);
    }
    return out;
  }

  /** Find by name, select the matches and frame them. Returns the names. */
  app.findNodes = (query, { exact = checked('find-exact') } = {}) => {
    const box = $('find-results');
    const text = String(query || '').trim();
    if (!text) {
      if (box) box.innerHTML = '';
      return [];
    }
    const names = matchNames(text, exact);
    if (!names.length) {
      if (box) {
        box.innerHTML = '<p class="result-note">No drawn segment has that name. ' +
          'If it was filtered out, widen the graph scope and redraw.</p>';
      }
      setStatus(`No segment matches “${text}”`);
      return [];
    }
    app.selectSegments(names, { focus: true });
    if (box) {
      box.innerHTML =
        `<p class="result-note">${fmtInt(names.length)} match${names.length === 1 ? '' : 'es'}</p>` +
        `<div class="dl-row name-chips">${names.slice(0, 40).map((n) =>
          `<a href="#" class="pill" data-goto="${esc(n)}">${esc(n)}</a>`).join('')}</div>` +
        (names.length > 40 ? `<p class="legend-more">+ ${fmtInt(names.length - 40)} more, all selected</p>` : '');
      wireGotoLinks(box);
    }
    setStatus(`${fmtInt(names.length)} segment(s) found`);
    return names;
  };

  /* ---- selection ---- */

  /** Every "click this name to go there" link in a rendered block. */
  function wireGotoLinks(box) {
    box.querySelectorAll('[data-goto]').forEach((b) =>
      b.addEventListener('click', (ev) => {
        ev.preventDefault();
        app.selectSegments([b.dataset.goto], { focus: true });
      }));
  }

  app.selectSegments = (names, { add = false, focus = false } = {}) => {
    renderer.selectNames(names, add);
    updateCounts();
    showSelectionDetail();
    if (focus) {
      const idx = (names || [])
        .map((n) => graph.segmentByName(n))
        .filter(Boolean)
        .map((s) => s.idx);
      if (idx.length) renderer.fitToView(idx);
    }
  };

  /** Frame whatever is selected, without changing the selection. */
  app.frameSelection = () => {
    const idx = renderer.selectedNames()
      .map((n) => graph.segmentByName(n))
      .filter(Boolean)
      .map((s) => s.idx);
    if (!idx.length) { app.toast('Nothing is selected', 'warn'); return; }
    renderer.fitToView(idx);
  };

  /**
   * Count, total length and mean depth of the selection.
   *
   * Depth is a per-base quantity, so the only average that means anything is
   * length-weighted — the same one Bandage divides by when it sizes a node.
   */
  function updateSelectionStats() {
    const segs = renderer.selectedNames().map((n) => graph.segmentByName(n)).filter(Boolean);
    for (const id of ['sel-frame', 'sel-copy', 'sel-fasta', 'sel-clear']) {
      const b = $(id);
      if (b) b.disabled = segs.length === 0;
    }
    const box = $('selection-stats');
    if (!box) return;
    if (!segs.length) {
      box.innerHTML = '<p class="result-note">Nothing selected.</p>';
      return;
    }
    let total = 0;
    let depthSum = 0;
    let depthLen = 0;
    for (const s of segs) {
      total += s.length;
      if (s.depth !== null) { depthSum += s.depth * s.length; depthLen += s.length; }
    }
    const mean = depthLen > 0 ? depthSum / depthLen : null;
    const cell = (label, value) =>
      `<div class="stat"><span class="stat-k">${esc(label)}</span><b class="stat-v">${esc(value)}</b></div>`;
    box.innerHTML =
      cell('Segments', fmtInt(segs.length)) +
      cell('Total length', fmtBp(total)) +
      cell('Mean depth', mean === null ? '—' : `${fmtNum(mean, 2)}×`);
  }

  /**
   * How many sequences the fallback below will fetch one at a time before it
   * gives up and asks for a narrower selection.
   */
  const FASTA_FALLBACK_LIMIT = 400;

  /**
   * FASTA for the named segments.
   *
   * The route is POST /api/fasta `{names, wrap}` -> `{fasta}` (a name may carry
   * a `+`/`-` suffix, and `-` gives the reverse complement). A server that
   * answers it with a real complaint — too many bases, no sequences in this
   * graph — has that complaint shown; one that does not have the endpoint at
   * all replies with something that is not our error shape, and the segment
   * endpoint every server has is used instead.
   */
  async function fastaFor(names) {
    try {
      const res = await fetch(`${API_BASE}/fasta`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
        body: JSON.stringify({ names, wrap: 60 }),
      });
      const text = await res.text();
      let data = null;
      try { data = JSON.parse(text); } catch { /* not our JSON */ }
      if (res.ok && typeof data?.fasta === 'string') return data.fasta;
      if (res.ok && text.trim().startsWith('>')) return text;
      if (data?.error) throw new ApiError(data.error, data.detail || '', res.status, '/fasta');
    } catch (err) {
      if (err instanceof ApiError) throw err;
      /* no such endpoint here; fall through */
    }

    if (names.length > FASTA_FALLBACK_LIMIT) {
      throw new Error(
        `This server has no bulk FASTA endpoint, so ${fmtInt(names.length)} sequences would be ` +
        `${fmtInt(names.length)} requests. Select at most ${FASTA_FALLBACK_LIMIT}, or export the graph.`);
    }
    const parts = [];
    for (const name of names) {
      const detail = await api.segment(name);
      const seq = detail.sequence || '';
      const lines = [];
      for (let i = 0; i < seq.length; i += 60) lines.push(seq.slice(i, i + 60));
      parts.push(`>${name} length=${detail.length}` +
        (detail.depth === null || detail.depth === undefined ? '' : ` depth=${detail.depth}`) +
        `\n${lines.join('\n')}\n`);
    }
    return parts.join('');
  }

  /** The selection as FASTA, or null with the reason already reported. */
  app.selectionFasta = async () => {
    const names = renderer.selectedNames();
    if (!names.length) { app.toast('Nothing is selected', 'warn'); return null; }
    return guard('Collecting sequences…', async () => {
      const fasta = await fastaFor(names);
      if (!fasta.trim()) {
        app.toast('The selected segments carry no sequence — this graph has none', 'warn', 6000);
        return null;
      }
      return fasta;
    });
  };

  async function showSelectionDetail() {
    updateSelectionStats();
    const box = $('selection-detail');
    if (!box) return;
    const names = renderer.selectedNames();
    if (!names.length) {
      box.innerHTML = '<p class="result-note">Click a segment to see its details. ' +
        'Shift-click adds to the selection; drag on empty space to box-select.</p>';
      return;
    }
    if (names.length > 1) {
      const segs = names.map((n) => graph.segmentByName(n)).filter(Boolean);
      box.innerHTML =
        `<h4>Selected contigs</h4><ul class="sel-list">` +
        segs.slice(0, 60).map((s) =>
          `<li><a href="#" data-goto="${esc(s.name)}"><code>${esc(s.name)}</code></a>` +
          `<span class="muted">${esc(fmtBp(s.length))}` +
          (s.depth === null ? '' : ` · ${fmtNum(s.depth, 1)}×`) + '</span></li>').join('') +
        (segs.length > 60 ? `<li class="muted">…and ${fmtInt(segs.length - 60)} more</li>` : '') +
        '</ul>';
      wireGotoLinks(box);
      return;
    }

    const name = names[0];
    const seg = graph.segmentByName(name);
    box.innerHTML = `<p class="result-note"><code>${esc(name)}</code> loading…</p>`;
    let detail;
    try {
      detail = await api.segment(name);
    } catch (err) {
      box.innerHTML = `<p class="result-note err">${esc(err.message)}</p>`;
      return;
    }

    const rows = [
      ['Length', fmtBp(detail.length)],
      ['Depth', detail.depth === null ? '—' : `${fmtNum(detail.depth, 2)}×`],
      ['GC', detail.gc === null ? '—' : `${fmtNum(detail.gc * 100, 2)}%`],
      ['Links', `${detail.deg_start} at start · ${detail.deg_end} at end`],
      ['Component', seg ? `#${seg.component}` : '—'],
      ['Circular', detail.circular ? 'yes' : 'no'],
    ];
    const hits = (detail.ref_hits || []).map((h) =>
      `<li><code>${esc(h.ref)}</code> ${fmtInt(h.r_st)}–${fmtInt(h.r_en)} ` +
      `<span class="muted">${h.strand > 0 ? '+' : '−'} ${fmtNum((h.identity || 0) * 100, 2)}%` +
      `${h.is_primary ? '' : ' (secondary)'}</span></li>`).join('');

    box.innerHTML =
      `<div class="results"><h4><code>${esc(name)}</code></h4>` +
      `<table><tbody>${rows.map(([k, v]) =>
        `<tr><th>${esc(k)}</th><td>${esc(v)}</td></tr>`).join('')}</tbody></table>` +
      (hits ? `<h4>Reference hits</h4><ul>${hits}</ul>` : '') +
      (detail.neighbours?.length
        ? `<h4>Neighbours</h4><div class="dl-row">${detail.neighbours.slice(0, 24).map((n) =>
            `<a href="#" class="pill" data-goto="${esc(n)}">${esc(n)}</a>`).join('')}</div>`
        : '') +
      (detail.sequence
        ? `<h4>Sequence</h4><textarea class="seq-box" readonly rows="4" aria-label="Sequence of ${esc(name)}">${esc(detail.sequence)}</textarea>` +
          `<div class="btn-row"><button class="btn btn-sm" id="copy-seq">Copy sequence</button></div>`
        : '') + '</div>';

    wireGotoLinks(box);
    box.querySelector('#copy-seq')?.addEventListener('click', async () => {
      try {
        await navigator.clipboard.writeText(detail.sequence);
        app.toast('Sequence copied', 'ok', 2000);
      } catch {
        box.querySelector('.seq-box')?.select();
        app.toast('Press Ctrl+C to copy', 'info');
      }
    });
  }
  app.showSelectionDetail = showSelectionDetail;

  /* ---- layout ---- */

  app.layout = new LayoutController(graph, {
    onTick: () => {
      graph.updateBounds();
      if (app._fitAfterLayout) renderer.fitToView(null);
      renderer.requestDraw();
    },
    onDone: (info) => {
      graph.updateBounds();
      if (app._fitAfterLayout) { app._fitAfterLayout = false; renderer.fitToView(null); }
      renderer.requestDraw();
      setLayoutRunning(false);
      $('status-layout').textContent = info?.iter ? `layout: ${fmtInt(info.iter)} iterations` : 'layout: idle';
      $('rr-status').textContent = info?.iter ? `Settled after ${fmtInt(info.iter)} iterations.` : '';
    },
    onError: (msg) => { setLayoutRunning(false); app.toast(String(msg), 'error'); },
    onState: (state) => {
      const el = $('status-layout');
      if (!el) return;
      el.textContent = typeof state === 'boolean'
        ? (state ? 'layout: running…' : 'layout: idle')
        : `layout: ${state}`;
    },
  });

  function setLayoutRunning(running) {
    const stopBtn = $('btn-stop-layout');
    if (stopBtn) stopBtn.hidden = !running;
    $('btn-rearrange')?.classList.toggle('running', running);
    const stop = $('rr-stop');
    if (stop) stop.disabled = !running;
    const apply = $('rr-apply');
    if (apply) apply.disabled = running;
    document.body.classList.toggle('layout-running', running);
  }

  function subsetFor(scope) {
    if (scope === 'selection') {
      const idx = renderer.selectedNames()
        .map((n) => graph.segmentByName(n))
        .filter(Boolean)
        .map((s) => s.idx);
      return idx.length ? idx : null;
    }
    if (scope === 'component') {
      const raw = $('rr-component')?.value;
      if (raw === '' || raw === undefined || raw === null) return null;
      // componentIndex is keyed by number; the select gives a string.
      const id = Number.isNaN(Number(raw)) ? raw : Number(raw);
      const members = graph.segmentsInComponent(id) || [];
      const idx = members.map((s) => (typeof s === 'number' ? s : s.idx))
        .filter((n) => Number.isFinite(n));
      return idx.length ? idx : null;
    }
    return null;
  }

  function startLayout({ mode = 'force', scope = 'all' } = {}) {
    if (graph.isEmpty) { app.toast('Load an assembly first', 'warn'); return; }
    const subset = subsetFor(scope);
    if (scope === 'selection' && !subset) {
      app.toast('Nothing is selected — select segments first, or choose "whole graph"', 'warn');
      return;
    }
    if (scope === 'component' && !subset) {
      app.toast('Pick a component first, or choose "whole graph"', 'warn');
      return;
    }
    const params = {
      ...DEFAULT_PARAMS,
      repulsion: num('rr-repulsion', DEFAULT_PARAMS.repulsion),
      linkStrength: num('rr-linkstr', DEFAULT_PARAMS.linkStrength),
      gravity: num('rr-gravity', DEFAULT_PARAMS.gravity),
      maxIter: int('rr-iters', DEFAULT_PARAMS.maxIter),
    };
    setLayoutRunning(true);
    $('rr-status').textContent = `Running ${mode} layout…`;
    const started = app.layout.start({ mode, subset, params });
    if (!started) setLayoutRunning(false);
  }
  app.startLayout = startLayout;

  app.stopLayout = () => { app.layout.stop(); setLayoutRunning(false); };

  /* ---- reference ---- */

  function clearReferenceResults() {
    const box = $('ref-results');
    if (box) box.innerHTML = '<p class="result-note">No reference loaded. Align one to get genome fraction, misassemblies, NGA50, and reference-guided scaffolding.</p>';
  }

  function renderReferenceResults(report) {
    const box = $('ref-results');
    if (!box) return;
    if (!report) { clearReferenceResults(); return; }
    const stat = (label, value, tone = '') =>
      `<div class="kv"><span>${esc(label)}</span>` +
      `<b class="${esc(tone)}">${esc(value)}</b></div>`;

    const gfTone = report.genome_fraction >= 90 ? 'ok' : (report.genome_fraction >= 70 ? '' : 'err');
    const misTone = report.num_misassemblies === 0 ? 'ok' : 'err';

    const coverage = Object.entries(report.per_reference_coverage || {}).map(([ref, pct]) => {
      const blocks = (report.coverage_blocks || {})[ref] || [];
      const span = blocks.length ? Math.max(...blocks.map((b) => b[1])) : 0;
      // The real length, not one inferred from span/percentage -- that
      // inference is short by exactly the uncovered tail, so a replicon
      // covered to 95.81% drew its bar as if it were 100% full.
      const total = (report.reference_sizes || {})[ref]
        || (pct > 0 ? (span / (pct / 100)) : span);
      const bars = blocks.map(([a, b]) => {
        const left = total ? (a / total) * 100 : 0;
        const width = total ? Math.max(0.25, ((b - a) / total) * 100) : 0;
        return `<i style="left:${left.toFixed(3)}%;width:${width.toFixed(3)}%"></i>`;
      }).join('');
      return `<div class="covbar-wrap"><div class="covbar-label"><code>${esc(ref)}</code>` +
        `<span>${fmtNum(pct, 2)}% covered</span></div>` +
        `<div class="covbar">${bars}</div></div>`;
    }).join('');

    // Two headline numbers and the coverage bars. The full statistics table
    // lives in Assembly QC; printing all of it twice made the sidebar longer
    // without telling anyone anything new.
    box.innerHTML =
      '<div class="results">' +
      stat('Genome fraction', `${fmtNum(report.genome_fraction, 2)}%`, gfTone) +
      stat('Misassemblies', fmtInt(report.num_misassemblies), misTone) +
      (coverage ? `<h4>Reference coverage</h4>${coverage}` : '') +
      (report.misassemblies?.length
        ? `<h4>Misassemblies</h4><ul>${report.misassemblies.slice(0, 40).map((m) =>
            `<li><a href="#" data-goto="${esc(m.contig)}"><code>${esc(m.contig)}</code></a> ` +
            `<span class="gap-tag ${esc(m.kind === 'local' ? 'manual' : 'reference')}">${esc(m.kind)}</span> ` +
            `<span class="muted">${esc(m.description)}</span></li>`).join('')}</ul>`
        : '') + '</div>';

    wireGotoLinks(box);
  }
  app.renderReferenceResults = renderReferenceResults;

  app.loadReference = () => guard('Aligning to reference — this can take a moment…', async () => {
    const path = $('ref-path').value.trim();
    if (!path) { app.toast('Choose a reference FASTA first', 'warn'); return; }
    const result = await api.addReference({
      path,
      preset: $('ref-preset').value,
      min_identity: num('ref-minid', 0),
      min_length: int('ref-minlen', 200),
      threads: int('ref-threads', 8),
    });
    renderReferenceResults(result.report);
    await app.refreshStatus();
    await app.reloadGraph({ keepPositions: true });
    // Reference colouring is the reason people load a reference; switch to it.
    $('colour-mode').value = 'reference';
    renderer.setColourMode('reference');
    updateColourLegend();
    await app.loadReport();
    app.toast(
      `Aligned: ${fmtNum(result.report.genome_fraction, 2)}% genome fraction, ` +
      `${result.report.num_misassemblies} misassemblies`, 'ok', 7000,
    );
  });

  /* ---- QC report ---- */

  app.loadReport = async () => {
    const data = await api.report().catch(() => null);
    if (!data) return null;
    app.report = data;
    renderQc(data);
    return data;
  };

  function renderQc(data) {
    const box = $('qc-tables');
    if (!box) return;
    const m = data.metrics || {};
    const r = data.reference;
    const row = (label, value) => `<tr><th>${esc(label)}</th><td>${esc(value)}</td></tr>`;

    let html = '<div class="results"><table><tbody>' +
      row('# contigs', fmtInt(m.num_contigs)) +
      row('# ≥ 1 kb', fmtInt(m.num_contigs_ge_1kb)) +
      row('Total length', fmtBp(m.total_length)) +
      row('Largest contig', fmtBp(m.largest_contig)) +
      row('N50', fmtBp(m.n50)) +
      row('L50', fmtInt(m.l50)) +
      row('N75', fmtBp(m.n75)) +
      row('auN', fmtInt(Math.round(m.auN || 0))) +
      (m.ng50 ? row('NG50', fmtBp(m.ng50)) : '') +
      (m.gc_percent !== null && m.gc_percent !== undefined ? row('GC', `${fmtNum(m.gc_percent, 2)}%`) : '') +
      (m.mean_depth !== null && m.mean_depth !== undefined ? row('Mean depth', `${fmtNum(m.mean_depth, 1)}×`) : '') +
      row('Links', fmtInt(m.num_links)) +
      row('Components', fmtInt(m.num_components)) +
      row('Dead ends', fmtInt(m.dead_ends)) +
      row('Circular contigs', fmtInt(m.num_circular)) +
      '</tbody></table>';

    if (r) {
      html += '<h4>Reference-based</h4><table><tbody>' +
        row('Genome fraction', `${fmtNum(r.genome_fraction, 2)}%`) +
        row('Duplication ratio', fmtNum(r.duplication_ratio, 3)) +
        row('NA50', fmtBp(r.na50)) +
        (r.nga50 ? row('NGA50', fmtBp(r.nga50)) : '') +
        row('Mismatches / 100 kb', fmtNum(r.mismatches_per_100kb, 2)) +
        row('Indels / 100 kb', fmtNum(r.indels_per_100kb, 2)) +
        row('Misassemblies', fmtInt(r.num_misassemblies)) +
        row('  relocations', fmtInt(r.num_relocations)) +
        row('  inversions', fmtInt(r.num_inversions)) +
        row('  translocations', fmtInt(r.num_translocations)) +
        row('Local misassemblies', fmtInt(r.num_local_misassemblies)) +
        row('Unaligned contigs', fmtInt(r.unaligned_contigs)) +
        '</tbody></table>';
    }
    box.innerHTML = html + '</div>';

    drawSeries($('plot-cumulative'), m.cumulative_curve, {
      xLabel: 'contig rank', yLabel: 'cumulative bp', kind: 'line',
    });
    drawSeries($('plot-nx'), m.nx_curve, {
      xLabel: 'x (%)', yLabel: 'Nx (bp)', kind: 'step',
    });
    drawSeries($('plot-hist'), m.length_histogram, {
      xLabel: 'length (bp)', yLabel: 'contigs', kind: 'bars', logX: true,
    });
  }
  app.renderQc = renderQc;

  /** Minimal inline-SVG plotting, so the panel needs no chart library. */
  function drawSeries(host, series, { xLabel = '', yLabel = '', kind = 'line', logX = false } = {}) {
    if (!host) return;
    if (!Array.isArray(series) || series.length < 2) {
      host.innerHTML = '<p class="result-note">no data</p>';
      return;
    }
    const W = 300, H = 150, padL = 44, padB = 26, padT = 8, padR = 8;
    const xs = series.map((p) => p[0]);
    const ys = series.map((p) => p[1]);
    const xMin = Math.min(...xs), xMax = Math.max(...xs);
    const yMax = Math.max(...ys, 1);
    const tx = (x) => {
      if (logX) {
        const lo = Math.log10(Math.max(xMin, 1)), hi = Math.log10(Math.max(xMax, 10));
        const t = hi > lo ? (Math.log10(Math.max(x, 1)) - lo) / (hi - lo) : 0.5;
        return padL + t * (W - padL - padR);
      }
      const t = xMax > xMin ? (x - xMin) / (xMax - xMin) : 0.5;
      return padL + t * (W - padL - padR);
    };
    const ty = (y) => H - padB - (y / yMax) * (H - padB - padT);

    let d = '';
    if (kind === 'bars') {
      const bw = Math.max(1.5, (W - padL - padR) / series.length - 1.5);
      d = series.map((p) =>
        `<rect x="${(tx(p[0]) - bw / 2).toFixed(2)}" y="${ty(p[1]).toFixed(2)}" ` +
        `width="${bw.toFixed(2)}" height="${Math.max(0, H - padB - ty(p[1])).toFixed(2)}" ` +
        `fill="var(--accent, #5b9bd5)" opacity="0.85"/>`).join('');
    } else {
      const cmds = [];
      series.forEach((p, i) => {
        const X = tx(p[0]).toFixed(2), Y = ty(p[1]).toFixed(2);
        if (i === 0) cmds.push(`M${X},${Y}`);
        else if (kind === 'step') cmds.push(`H${X}`, `V${Y}`);
        else cmds.push(`L${X},${Y}`);
      });
      d = `<path d="${cmds.join(' ')}" fill="none" stroke="var(--accent, #5b9bd5)" stroke-width="1.8"/>`;
    }

    host.className = 'plot';
    host.innerHTML =
      (yLabel ? `<div class="plot-title">${esc(yLabel)}</div>` : '') +
      `<svg viewBox="0 0 ${W} ${H}" role="img">` +
      `<line x1="${padL}" y1="${padT}" x2="${padL}" y2="${H - padB}" stroke="currentColor" opacity="0.35"/>` +
      `<line x1="${padL}" y1="${H - padB}" x2="${W - padR}" y2="${H - padB}" stroke="currentColor" opacity="0.35"/>` +
      `<text x="${padL - 4}" y="${padT + 8}" text-anchor="end" font-size="8" fill="currentColor" opacity="0.7">` +
      `${esc(kind === 'bars' ? fmtInt(yMax) : fmtBp(yMax))}</text>` +
      `<text x="${padL - 4}" y="${H - padB}" text-anchor="end" font-size="8" fill="currentColor" opacity="0.7">0</text>` +
      `<text x="${(W + padL) / 2}" y="${H - 4}" text-anchor="middle" font-size="8" fill="currentColor" opacity="0.7">${esc(xLabel)}</text>` +
      d + '</svg>';
  }

  /* ---- graph operations ---- */

  app.runOp = (op, args, label) => guard(label || 'Applying…', async () => {
    const result = await api.op(op, args);
    await app.refreshStatus();
    await app.reloadGraph({ keepPositions: true });
    await app.loadReport();
    app.scaffoldPanel?.clear();
    app.toast(`${op}: ${fmtInt(result.count)} change(s)`, 'ok');
    return result;
  });

  /* ---- search ---- */

  app.runSearch = () => guard('Searching…', async () => {
    const query = $('search-query').value.trim();
    if (!query) { app.toast('Paste a sequence, or pick a FASTA file', 'warn'); return; }
    const result = await api.search(query, num('search-minid', 0.8));
    graph.setSearchHits(result.hits);
    renderer.setColourMode('search');
    $('colour-mode').value = 'search';
    updateColourLegend();
    renderer.requestDraw();

    const box = $('search-results');
    if (result.hits.length) {
      box.innerHTML =
        `<div class="results"><p class="result-note">${fmtInt(result.hits.length)} hit(s) via ${esc(result.backend)}</p>` +
        `<ul>${result.hits.slice(0, 100).map((h) =>
          `<li><a href="#" data-goto="${esc(h.segment)}"><code>${esc(h.segment)}</code></a> ` +
          `<span class="muted">${fmtNum(h.identity * 100, 1)}% · ${fmtBp(h.length)} · ` +
          `${fmtInt(h.s_st)}–${fmtInt(h.s_en)} ${h.strand > 0 ? '+' : '−'}</span></li>`).join('')}</ul></div>`;
      wireGotoLinks(box);
      app.selectSegments(result.hits.map((h) => h.segment));
    } else {
      box.innerHTML = '<p class="result-note">No hits.</p>';
    }
  });

  app.clearSearch = () => {
    graph.clearSearchHits();
    $('search-results').innerHTML = '';
    $('search-query').value = '';
    if ($('colour-mode').value === 'search') {
      $('colour-mode').value = 'uniform';
      renderer.setColourMode('uniform');
    }
    updateColourLegend();
    renderer.requestDraw();
  };

  /* ---- export ---- */

  app.runExport = () => guard('Writing files…', async () => {
    const what = [...$('export-what').querySelectorAll('input[type=checkbox]:checked')]
      .map((c) => c.value);
    if (!what.length) { app.toast('Pick at least one thing to export', 'warn'); return; }
    const outdir = $('export-outdir').value.trim() || 'plastr_out';
    const result = await api.exportFiles(outdir, what);
    $('export-results').innerHTML =
      `<div class="results"><ul>${result.written.map((w) =>
        `<li><code>${esc(w.path)}</code> <span class="muted">${fmtInt(w.bytes)} bytes</span></li>`).join('')}</ul></div>`;
    // Scaffolds and AGP are skipped server-side when there is no plan, so a
    // request for them alone came back empty and still toasted success.
    const wrote = new Set(result.written.map((w) => w.kind));
    const skipped = what.filter((k) => !wrote.has(k));
    if (skipped.length) {
      const needsPlan = skipped.filter((k) => k === 'scaffolds' || k === 'agp');
      app.toast(
        `Wrote ${result.written.length} file(s); ${skipped.join(', ')} not written`
        + (needsPlan.length ? ' — build a scaffold plan first' : ''),
        result.written.length ? 'warn' : 'error',
        9000,
      );
      return;
    }
    app.toast(`Wrote ${result.written.length} file(s) to ${outdir}`, 'ok', 6000);
  });

  function buildDownloadLinks() {
    const box = $('download-links');
    if (!box) return;
    const kinds = [
      ['scaffolds', 'scaffolds.fasta'], ['agp', 'scaffolds.agp'], ['gfa', 'graph.gfa'],
      ['csv', 'segments.csv'], ['report', 'report.html'], ['session', 'session.json'],
    ];
    box.innerHTML = `<div class="dl-row">${kinds.map(([kind, label]) =>
      `<a href="${esc(api.downloadUrl(kind))}" data-kind="${esc(kind)}" download>${esc(label)}</a>`).join('')}</div>`;
    // A plain <a download> never sees a non-2xx: the browser either cancels or
    // saves the JSON error body under a .fasta name, with nothing shown in the
    // app. Fetch it instead so the server's message reaches the user.
    box.querySelectorAll('a[data-kind]').forEach((a) => {
      a.addEventListener('click', async (ev) => {
        ev.preventDefault();
        try {
          const res = await fetch(a.href);
          if (!res.ok) {
            let message = `HTTP ${res.status}`;
            try { message = (await res.json()).error || message; } catch { /* not JSON */ }
            app.showError(message);
            return;
          }
          download(await res.blob(), a.textContent.trim());
        } catch (err) {
          app.showError((err && err.message) || String(err));
        }
      });
    });
  }

  /* ---- image export ---- */

  function download(blobOrText, filename, type) {
    const blob = blobOrText instanceof Blob ? blobOrText : new Blob([blobOrText], { type });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 4000);
  }
  app.download = download;

  app.savePNG = async () => {
    try {
      const blob = await renderer.toPNG(2);
      download(blob, 'plastr-graph.png');
      app.toast('PNG saved', 'ok', 2500);
    } catch (err) { app.showError(err); }
  };

  app.saveSVG = () => {
    try {
      download(renderer.toSVG(), 'plastr-graph.svg', 'image/svg+xml');
      app.toast('SVG saved', 'ok', 2500);
    } catch (err) { app.showError(err); }
  };

  /* ---- file browser modal ---- */

  const browser = { target: null, current: '', chooseDir: false };

  app.openBrowser = async (targetInputId, { chooseDir = false, title = 'Choose a file' } = {}) => {
    browser.target = targetInputId;
    browser.chooseDir = chooseDir;
    $('browse-title').textContent = title;
    $('modal-browse').hidden = false;
    $('browse-choose-dir').hidden = !chooseDir;
    // Seeding with a directory that does not exist yet (the default output
    // folder) would open the picker on a 404.
    const start = $(targetInputId)?.value || '';
    await browseTo(start.includes('/') ? start : '');
  };

  const closeBrowser = () => { $('modal-browse').hidden = true; };

  async function browseTo(path) {
    try {
      const data = await api.browse(path || '');
      browser.current = data.path;
      $('browse-path').value = data.path;
      $('browse-up').disabled = !data.parent;
      $('browse-up').dataset.parent = data.parent || '';
      $('browse-sel').textContent = '';
      const list = $('browse-list');
      list.innerHTML = data.entries.map((e) =>
        `<li class="${e.is_dir ? 'dir' : ''}" data-name="${esc(e.name)}" data-dir="${e.is_dir ? 1 : 0}">` +
        `<span class="bi">${e.is_dir ? '▸' : '·'}</span>` +
        `<span class="bname">${esc(e.name)}</span>` +
        `<span class="bsize">${e.is_dir ? '' : fmtInt(e.size) + ' B'}</span></li>`).join('')
        || '<li class="muted">nothing here that Plastr can read</li>';

      list.querySelectorAll('li[data-name]').forEach((b) => {
        b.addEventListener('click', () => {
          const full = data.path.replace(/\/$/, '') + '/' + b.dataset.name;
          if (b.dataset.dir === '1') { browseTo(full); return; }
          list.querySelectorAll('li').forEach((x) => x.classList.remove('sel'));
          b.classList.add('sel');
          $('browse-sel').textContent = full;
        });
        b.addEventListener('dblclick', () => {
          if (b.dataset.dir === '1') return;
          $('browse-sel').textContent = data.path.replace(/\/$/, '') + '/' + b.dataset.name;
          $('browse-choose').click();
        });
      });
    } catch (err) {
      app.showError(err);
    }
  }

  on('browse-close', 'click', closeBrowser);
  on('browse-up', 'click', (e) => browseTo(e.currentTarget.dataset.parent || ''));
  on('browse-go', 'click', () => browseTo($('browse-path').value));
  on('browse-path', 'keydown', (e) => { if (e.key === 'Enter') browseTo($('browse-path').value); });
  on('browse-choose-dir', 'click', () => {
    if (browser.target) $(browser.target).value = browser.current;
    closeBrowser();
  });
  on('browse-choose', 'click', () => {
    const chosen = $('browse-sel').textContent.trim();
    if (!chosen) { app.toast('Pick a file first', 'warn'); return; }
    if (browser.target) $(browser.target).value = chosen;
    closeBrowser();
  });
  $('modal-browse')?.addEventListener('click', (e) => {
    if (e.target.id === 'modal-browse' || e.target.classList.contains('modal-backdrop')) closeBrowser();
  });

  /* ---- settings ---- */

  /**
   * Drawing tunables, Bandage's settings dialog cut down to the ones that
   * change the picture rather than the algorithm.
   *
   * There is deliberately no antialiasing toggle. Bandage has one because Qt
   * can be told to stop antialiasing a QPainterPath, which on a large graph is
   * a real speed-up; the Canvas 2D API has no equivalent -- `imageSmoothing`
   * governs drawn images, not path edges -- so the control could only ever
   * have been a switch that did nothing.
   */
  const SETTINGS_KEY = 'plastr-settings';
  const SETTINGS_DEFAULTS = Object.freeze({
    nodeWidth: 5,
    depthPower: 0.5,
    depthEffect: 0.5,
    edgeWidth: 1.1,
    textSize: 11,
    outline: true,
    autoLength: true,
    perMegabase: 1000,
    labelName: false,
    labelLength: false,
    labelDepth: false,
  });
  let settings = { ...SETTINGS_DEFAULTS };

  const SETTINGS_CONTROLS = [
    ['nodeWidth', 'set-nodewidth', 'set-nodewidth-out', 1],
    ['depthPower', 'set-depthpower', 'set-depthpower-out', 2],
    ['depthEffect', 'set-deptheffect', 'set-deptheffect-out', 2],
    ['edgeWidth', 'set-edgewidth', 'set-edgewidth-out', 1],
    ['textSize', 'set-textsize', 'set-textsize-out', 0],
  ];

  function saveSettings() {
    try { localStorage.setItem(SETTINGS_KEY, JSON.stringify(settings)); } catch { /* private mode */ }
  }

  /** Push the numbers at the renderer. */
  function applyRenderSettings() {
    renderer.setOption('baseWidth', settings.nodeWidth);
    renderer.setOption('depthPower', settings.depthPower);
    renderer.setOption('depthEffectOnWidth', settings.depthEffect);
    renderer.setOption('edgeWidth', settings.edgeWidth);
    renderer.setOption('labelTextSize', settings.textSize);
    renderer.setOption('outlineNodes', settings.outline);
    renderer.setOption('labelName', settings.labelName);
    renderer.setOption('labelLength', settings.labelLength);
    renderer.setOption('labelDepth', settings.labelDepth);
    // The renderer's master switch: with it off nothing is drawn whatever the
    // three fields say, so it is the union of them.
    renderer.setOption('showLabels', settings.labelName || settings.labelLength || settings.labelDepth);
  }

  /**
   * Re-map every drawn length after the calibration rule changed.
   *
   * `calibrate` needs the segment count and the base count it was given at load
   * time; they are recovered from the model rather than remembered, because a
   * graph operation can change both between two visits to this panel.
   */
  function applyLengthSettings() {
    graph.geom.autoLength = settings.autoLength;
    graph.geom.unitsPerMegabase = settings.perMegabase;
    const perMb = $('set-permb');
    if (perMb) perMb.disabled = settings.autoLength;
    if (graph.isEmpty) return;
    let bases = 0;
    for (const seg of graph.segments) bases += seg.length;
    graph.calibrate(graph.segments.length, bases);
    graph.rescale(graph.geom.lengthScale || 1);
    graph.updateBounds();
    renderer.updateStyle();
    app.layout?.markDirty();
    renderer.requestDraw();
  }

  /** Write the settings into their controls, then apply them. */
  function syncSettingsControls() {
    for (const [key, slider, out, dp] of SETTINGS_CONTROLS) {
      const el = $(slider);
      if (el) el.value = String(settings[key]);
      if ($(out)) $(out).value = Number(settings[key]).toFixed(dp);
    }
    const pairs = [
      ['set-outline', 'outline'],
      ['set-autolength', 'autoLength'],
      ['lab-name', 'labelName'], ['lab-length', 'labelLength'], ['lab-depth', 'labelDepth'],
    ];
    for (const [id, key] of pairs) { const el = $(id); if (el) el.checked = !!settings[key]; }
    const perMb = $('set-permb');
    if (perMb) perMb.value = String(settings.perMegabase);
    applyRenderSettings();
    applyLengthSettings();
  }

  function loadSettings() {
    try {
      const saved = JSON.parse(localStorage.getItem(SETTINGS_KEY) || '{}') || {};
      for (const key of Object.keys(SETTINGS_DEFAULTS)) {
        if (saved[key] !== undefined) settings[key] = saved[key];
      }
    } catch { /* private mode, or a stale shape */ }
    syncSettingsControls();
  }
  app.applySettings = syncSettingsControls;

  /* ---- wiring ---- */

  function wire() {
    // Load panel
    on('load-go', 'click', () => app.loadFile($('load-path').value.trim(), $('load-format').value));
    on('load-path', 'keydown', (e) => { if (e.key === 'Enter') $('load-go').click(); });
    on('load-browse', 'click', () => app.openBrowser('load-path', { title: 'Choose an assembly' }));
    on('load-refresh', 'click', () => guard('Reloading…', () => app.reloadGraph({ keepPositions: true })));
    on('paths-browse', 'click', () => app.openBrowser('paths-path', { title: 'Choose contigs.paths' }));
    on('paths-go', 'click', () => guard('Attaching paths…', async () => {
      const r = await api.loadPaths($('paths-path').value.trim());
      await app.reloadGraph({ keepPositions: true });
      app.toast(`Attached ${fmtInt(r.added)} path(s)`, 'ok');
    }));

    // Graph scope
    const syncScopeArgs = () => {
      const scope = $('scope-mode')?.value || 'entire';
      for (const [name, id] of [['around', 'scope-args-around'], ['depth', 'scope-args-depth'],
        ['component', 'scope-args-component']]) {
        const el = $(id);
        if (el) el.hidden = scope !== name;
      }
    };
    on('scope-mode', 'change', syncScopeArgs);
    syncScopeArgs();
    on('scope-nodes', 'keydown', (e) => { if (e.key === 'Enter') $('lod-apply').click(); });
    on('scope-from-sel', 'click', () => {
      const names = renderer.selectedNames();
      if (!names.length) { app.toast('Nothing is selected', 'warn'); return; }
      $('scope-nodes').value = names.join(', ');
      // Names taken from the graph are exact by construction.
      const exact = document.querySelector('input[name="scope-match"][value="exact"]');
      if (exact) exact.checked = true;
    });
    on('scope-reset', 'click', () => {
      $('scope-mode').value = 'entire';
      syncScopeArgs();
      guard('Redrawing…', () => app.reloadGraph({ keepPositions: true }));
    });
    on('lod-apply', 'click', () => guard('Redrawing…', () => app.reloadGraph({ keepPositions: true })));

    // Find nodes
    on('find-go', 'click', () => app.findNodes($('find-input').value));
    on('find-input', 'keydown', (e) => {
      if (e.key === 'Enter') { e.preventDefault(); app.findNodes(e.target.value); }
    });
    on('find-exact', 'change', () => {
      const q = $('find-input')?.value.trim();
      if (q) app.findNodes(q);
    });

    // Selection
    on('sel-frame', 'click', () => app.frameSelection());
    on('sel-clear', 'click', () => renderer.clearSelection());
    on('sel-copy', 'click', async () => {
      const fasta = await app.selectionFasta();
      if (!fasta) return;
      try {
        await navigator.clipboard.writeText(fasta);
        app.toast(`Copied ${fmtInt(renderer.selectedNames().length)} sequence(s)`, 'ok', 2500);
      } catch {
        // Clipboard access needs a secure context; saving is the way out.
        app.download(fasta, 'plastr-selection.fasta', 'text/x-fasta');
        app.toast('The clipboard is unavailable here, so the FASTA was saved instead', 'info', 6000);
      }
    });
    on('sel-fasta', 'click', async () => {
      const fasta = await app.selectionFasta();
      if (!fasta) return;
      app.download(fasta, 'plastr-selection.fasta', 'text/x-fasta');
      app.toast('FASTA saved', 'ok', 2500);
    });

    // Node labels and settings
    for (const [id, key] of [['lab-name', 'labelName'], ['lab-length', 'labelLength'],
      ['lab-depth', 'labelDepth'], ['set-outline', 'outline']]) {
      on(id, 'change', () => { settings[key] = checked(id); applyRenderSettings(); saveSettings(); });
    }
    for (const [key, slider, out, dp] of SETTINGS_CONTROLS) {
      on(slider, 'input', (e) => {
        settings[key] = Number(e.target.value);
        if ($(out)) $(out).value = settings[key].toFixed(dp);
        applyRenderSettings();
        saveSettings();
      });
    }
    on('set-autolength', 'change', () => {
      settings.autoLength = checked('set-autolength');
      applyLengthSettings();
      saveSettings();
    });
    on('set-permb', 'change', () => {
      settings.perMegabase = Math.max(1, num('set-permb', SETTINGS_DEFAULTS.perMegabase));
      applyLengthSettings();
      saveSettings();
    });
    on('set-reset', 'click', () => {
      settings = { ...SETTINGS_DEFAULTS };
      syncSettingsControls();
      saveSettings();
      app.toast('Drawing settings restored', 'ok', 2500);
    });

    // Colour / style
    on('colour-mode', 'change', (e) => {
      renderer.setColourMode(e.target.value);
      updateColourLegend();
    });
    const styleToggle = (id, key) => on(id, 'change', (e) => {
      renderer.setOption(key, e.target.checked);
    });
    styleToggle('opt-depth-width', 'depthWidth');
    styleToggle('opt-arrows', 'showArrows');
    styleToggle('opt-links', 'showLinks');
    on('width-scale', 'input', (e) => {
      $('width-scale-out').value = Number(e.target.value).toFixed(1);
      renderer.setOption('widthScale', Number(e.target.value));
    });
    on('size-scale', 'input', (e) => {
      const v = Number(e.target.value) || 1;
      $('size-scale-out').value = v.toFixed(1);
      // Node length is a property of the model, not a render option: the
      // polyline geometry and the layout's rest lengths both derive from it.
      graph.rescale(v);
      graph.updateBounds();
      renderer.updateStyle();
      app.layout?.markDirty();
      renderer.requestDraw();
    });

    // Rearrange
    const pop = $('rearrange-pop');
    /** Place the popover under the control that opened it, kept on screen. */
    const anchorPop = () => {
      const anchor = $('btn-rearrange-menu') || $('btn-rearrange');
      if (!pop || !anchor) return;
      const r = anchor.getBoundingClientRect();
      const width = pop.offsetWidth || 300;
      const left = Math.max(8, Math.min(r.right - width, window.innerWidth - width - 8));
      pop.style.left = `${left}px`;
      pop.style.right = 'auto';
      pop.style.top = `${Math.min(r.bottom + 6, window.innerHeight - 80)}px`;
    };
    const togglePop = (show) => {
      if (!pop) return;
      const next = show === undefined ? pop.hidden : show;
      pop.hidden = !next;
      if (next) anchorPop();
    };
    window.addEventListener('resize', () => { if (pop && !pop.hidden) anchorPop(); });
    on('btn-rearrange', 'click', () => startLayout({ mode: radio('rr-mode') || 'force', scope: radio('rr-scope') || 'all' }));
    on('btn-rearrange-menu', 'click', () => togglePop());
    on('rearrange-close', 'click', () => togglePop(false));
    on('rr-apply', 'click', () => startLayout({ mode: radio('rr-mode') || 'force', scope: radio('rr-scope') || 'all' }));
    on('rr-stop', 'click', () => app.stopLayout());
    on('btn-stop-layout', 'click', () => app.stopLayout());
    // Decimals must match each slider's step, or the readout disagrees with the
    // value actually used -- repulsion steps by 0.005 and read '0.1' at 0.05.
    //
    // The sliders are also *seeded* from DEFAULT_PARAMS rather than trusting the
    // value in the markup. `startLayout` reads the physics out of the DOM, so a
    // stale `value=` attribute silently overrides the engine's defaults and
    // tuning the engine appears to do nothing at all.
    for (const [slider, out, dp, key] of [
      ['rr-repulsion', 'rr-rep-out', 3, 'repulsion'],
      ['rr-linkstr', 'rr-link-out', 2, 'linkStrength'],
      ['rr-gravity', 'rr-grav-out', 3, 'gravity'],
    ]) {
      const el = $(slider);
      if (el && DEFAULT_PARAMS[key] !== undefined) {
        el.value = String(DEFAULT_PARAMS[key]);
        if ($(out)) $(out).value = Number(el.value).toFixed(dp);
      }
      on(slider, 'input', (e) => { $(out).value = Number(e.target.value).toFixed(dp); });
    }
    const iters = $('rr-iters');
    if (iters) iters.value = String(DEFAULT_PARAMS.maxIter);
    document.addEventListener('click', (e) => {
      if (!pop || pop.hidden) return;
      if (pop.contains(e.target) || e.target.closest('#btn-rearrange-menu')) return;
      togglePop(false);
    });

    // View
    on('btn-fit', 'click', () => renderer.fitToView());
    on('btn-save-png', 'click', () => app.savePNG());
    on('btn-save-svg', 'click', () => app.saveSVG());
    on('btn-theme', 'click', () => {
      const root = document.documentElement;
      const current = root.dataset.theme
        || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
      const next = current === 'dark' ? 'light' : 'dark';
      root.dataset.theme = next;
      try { localStorage.setItem('plastr-theme', next); } catch { /* private mode */ }
      renderer.refreshTheme();
      updateColourLegend();
    });
    on('btn-help', 'click', () => { $('help-overlay').hidden = false; });
    on('help-close', 'click', () => { $('help-overlay').hidden = true; });
    $('help-overlay')?.addEventListener('click', (e) => {
      if (e.target.id === 'help-overlay' || e.target.classList.contains('modal-backdrop')) {
        $('help-overlay').hidden = true;
      }
    });

    // Reference
    on('ref-browse', 'click', () => app.openBrowser('ref-path', { title: 'Choose a reference FASTA' }));
    on('ref-go', 'click', () => app.loadReference());
    on('ref-drop', 'click', () => guard('Removing reference…', async () => {
      await api.dropReference();
      clearReferenceResults();
      await app.refreshStatus();
      await app.reloadGraph({ keepPositions: true });
      await app.loadReport();
      if ($('colour-mode').value === 'reference') {
        $('colour-mode').value = 'uniform';
        renderer.setColourMode('uniform');
        updateColourLegend();
      }
    }));

    // QC
    on('qc-refresh', 'click', () => guard('Computing…', () => app.loadReport()));

    // Operations
    on('op-delete', 'click', () => {
      const names = renderer.selectedNames();
      if (!names.length) { app.toast('Nothing is selected', 'warn'); return; }
      app.runOp('delete', { names }, `Deleting ${names.length} segment(s)…`);
    });
    on('op-simplify', 'click', () => app.runOp('simplify', {}, 'Merging unbranching paths…'));
    on('op-break', 'click', () =>
      app.runOp('break_misassemblies', { extensive_only: checked('op-break-extensive') },
        'Breaking misassembled contigs…'));
    on('op-reverse', 'click', () => {
      const names = renderer.selectedNames();
      if (names.length !== 1) { app.toast('Select exactly one segment to reverse', 'warn'); return; }
      app.runOp('reverse', { name: names[0] }, 'Reversing…');
    });
    on('op-filter', 'click', () => app.runOp('filter', {
      min_length: int('flt-minlen', 0),
      min_depth: $('flt-mindepth').value === '' ? null : num('flt-mindepth', 0),
      max_depth: $('flt-maxdepth').value === '' ? null : num('flt-maxdepth', 0),
      keep_components: $('flt-keepcomp').value === '' ? null : int('flt-keepcomp', 0),
      min_component_length: int('flt-mincomplen', 0),
    }, 'Filtering…'));
    on('op-split', 'click', () => {
      const name = $('split-name').value.trim() || renderer.selectedNames()[0];
      if (!name) { app.toast('Name a segment to split', 'warn'); return; }
      const positions = $('split-pos').value.split(/[,\s]+/).map((s) => parseInt(s, 10)).filter(Number.isFinite);
      if (!positions.length) { app.toast('Give at least one split position', 'warn'); return; }
      app.runOp('split', { name, positions }, 'Splitting…');
    });
    on('merge-from-sel', 'click', () => {
      $('merge-steps').value = renderer.selectedNames().map((n) => `${n}+`).join(',');
    });
    on('op-merge', 'click', () => {
      const steps = $('merge-steps').value.split(/[,\s]+/).filter(Boolean);
      if (steps.length < 2) { app.toast('Give at least two steps, e.g. a+,b+', 'warn'); return; }
      app.runOp('merge', { steps }, 'Merging…');
    });
    on('op-undo', 'click', () => guard('Undoing…', async () => {
      await api.undo();
      await app.refreshStatus();
      await app.reloadGraph({ keepPositions: true });
      await app.loadReport();
      app.scaffoldPanel?.clear();
      app.toast('Undone', 'ok', 2000);
    }));

    // Search
    on('search-go', 'click', () => app.runSearch());
    on('search-clear', 'click', () => app.clearSearch());
    on('search-browse', 'click', () => app.openBrowser('search-query', { title: 'Choose a query FASTA' }));

    // Export
    on('export-go', 'click', () => app.runExport());
    on('export-browse', 'click', () => app.openBrowser('export-outdir', { chooseDir: true, title: 'Choose an output directory' }));
    on('session-save', 'click', () => guard('Saving session…', async () => {
      await api.saveSession({
        layout: graph.exportLayout(),
        settings: {
          colourMode: $('colour-mode').value,
          widthScale: num('width-scale', 1),
          showLinks: checked('opt-links'),
          labels: {
            name: checked('lab-name'), length: checked('lab-length'), depth: checked('lab-depth'),
          },
        },
      });
      app.toast('Session saved on the server — use the session.json link to keep a copy', 'ok', 6000);
    }));
    on('session-load', 'click', () => guard('Restoring session…', async () => {
      const data = await api.getSession();
      if (data.layout) { graph.importLayout(data.layout); graph.updateBounds(); renderer.requestDraw(); }
      const s = data.settings || {};
      if (s.colourMode) { $('colour-mode').value = s.colourMode; renderer.setColourMode(s.colourMode); }
      if (s.labels) {
        settings.labelName = !!s.labels.name;
        settings.labelLength = !!s.labels.length;
        settings.labelDepth = !!s.labels.depth;
        syncSettingsControls();
        saveSettings();
      }
      updateColourLegend();
      app.toast('Session restored', 'ok');
    }));

    // Canvas interaction
    renderer.attach({
      // The renderer emits ('hover', segmentIndex, pointerEvent).
      hover: (idx, event) => {
        const tip = $('hover-tip');
        if (!tip) return;
        const seg = idx >= 0 ? graph.segments[idx] : null;
        if (!seg) { tip.hidden = true; return; }
        const hit = bestRefHit(seg.refHits);
        tip.innerHTML =
          `<strong>${esc(seg.name)}</strong><br>${esc(fmtBp(seg.length))}` +
          (seg.depth !== null ? ` · ${fmtNum(seg.depth, 1)}×` : '') +
          (seg.gc !== null ? ` · GC ${fmtNum(seg.gc * 100, 1)}%` : '') +
          (hit ? `<br><span class="muted">${esc(hit.ref)} ${fmtInt(hit.r_st)}–${fmtInt(hit.r_en)}</span>` : '');
        tip.hidden = false;
        if (event && typeof event.clientX === 'number') {
          const stage = $('stage').getBoundingClientRect();
          // Flip to the other side of the cursor near the right/bottom edge so
          // the tip is never clipped away by the stage's overflow:hidden.
          const x = event.clientX - stage.left;
          const y = event.clientY - stage.top;
          const w = tip.offsetWidth || 180;
          const h = tip.offsetHeight || 48;
          tip.style.left = `${Math.max(0, Math.min(x + 14, stage.width - w - 8))}px`;
          tip.style.top = `${Math.max(0, Math.min(y + 14, stage.height - h - 8))}px`;
        }
      },
      selection: () => { updateCounts(); showSelectionDetail(); },
      viewchange: () => {
        const z = $('status-zoom');
        if (z) z.textContent = `zoom ${(renderer.view.scale * 100).toFixed(0)}%`;
      },
      dragend: () => { graph.updateBounds(); app.layout?.syncPositions(); },
    });

    /** Open the Find panel if it is collapsed, then put the caret in its box. */
    const focusFind = () => {
      const panel = $('panel-find');
      if (panel?.classList.contains('collapsed')) panel.querySelector('.panel-head')?.click();
      const box = $('find-input');
      box?.focus();
      box?.select();
    };

    // The canvas is a focusable control, so it has to answer the keyboard: a
    // mouse is not the only way to reach a part of the drawing.
    canvas?.addEventListener('keydown', (e) => {
      if (e.ctrlKey || e.metaKey || e.altKey) return;
      const step = e.shiftKey ? 140 : 40;
      switch (e.key) {
        case 'ArrowLeft': renderer.panBy(step, 0); break;
        case 'ArrowRight': renderer.panBy(-step, 0); break;
        case 'ArrowUp': renderer.panBy(0, step); break;
        case 'ArrowDown': renderer.panBy(0, -step); break;
        case '+': case '=':
          renderer.zoomAt(canvas.clientWidth / 2, canvas.clientHeight / 2, 1.25); break;
        case '-': case '_':
          renderer.zoomAt(canvas.clientWidth / 2, canvas.clientHeight / 2, 1 / 1.25); break;
        default: return;
      }
      e.preventDefault();
    });

    // Keyboard shortcuts
    window.addEventListener('keydown', (e) => {
      const typing = ['INPUT', 'TEXTAREA', 'SELECT'].includes(e.target.tagName) || e.target.isContentEditable;
      if (e.key === 'Escape') {
        if ($('help-overlay')) $('help-overlay').hidden = true;
        if ($('modal-browse')) $('modal-browse').hidden = true;
        if ($('rearrange-pop')) $('rearrange-pop').hidden = true;
        if (!typing) renderer.clearSelection();
        return;
      }
      if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'f') {
        // The browser's own find cannot see anything painted on a canvas.
        e.preventDefault(); focusFind(); return;
      }
      if (typing) return;
      if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'z') {
        e.preventDefault(); $('op-undo')?.click(); return;
      }
      if ((e.ctrlKey || e.metaKey) && e.key.toLowerCase() === 'a') {
        e.preventDefault(); renderer.selectAllVisible(); updateCounts(); showSelectionDetail(); return;
      }
      if (e.ctrlKey || e.metaKey || e.altKey) return;
      switch (e.key.toLowerCase()) {
        case 'r': e.preventDefault(); $('btn-rearrange').click(); break;
        case 'f': e.preventDefault(); renderer.fitToView(null); break;
        case 'c': e.preventDefault(); $('sel-copy')?.click(); break;
        case '/': e.preventDefault(); focusFind(); break;
        case '?': e.preventDefault(); if ($('help-overlay')) $('help-overlay').hidden = false; break;
        case 'a':
          e.preventDefault();
          renderer.selectAllVisible();
          updateCounts();
          showSelectionDetail();
          break;
        case 'delete': case 'backspace': e.preventDefault(); $('op-delete')?.click(); break;
        default: break;
      }
    });

    // Collapsible sidebar sections. Which ones are open is remembered, because
    // the sidebar is taller than any screen and re-opening the same three
    // panels on every launch is a chore.
    const PANEL_KEY = 'plastr-panels';
    let panelState = {};
    try { panelState = JSON.parse(localStorage.getItem(PANEL_KEY) || '{}') || {}; } catch { panelState = {}; }
    document.querySelectorAll('.panel').forEach((panel) => {
      const head = panel.querySelector('.panel-head');
      if (!head) return;
      const id = panel.id || '';
      if (Object.prototype.hasOwnProperty.call(panelState, id)) {
        panel.classList.toggle('collapsed', !panelState[id]);
      }
      const button = head.querySelector('.panel-toggle');
      const sync = () => {
        const open = !panel.classList.contains('collapsed');
        if (button) button.setAttribute('aria-expanded', String(open));
        panelState[id] = open;
        try { localStorage.setItem(PANEL_KEY, JSON.stringify(panelState)); } catch { /* private mode */ }
      };
      sync();
      head.addEventListener('click', () => { panel.classList.toggle('collapsed'); sync(); });
    });

    window.addEventListener('resize', () => renderer.resize());
    buildDownloadLinks();
    clearReferenceResults();
    loadSettings();
    updateSelectionStats();
    showSelectionDetail();
  }

  app.wire = wire;
  return app;
}

/* --------------------------------------------------------------- bootstrap */

export async function start() {
  try {
    const saved = localStorage.getItem('plastr-theme');
    if (saved) document.documentElement.dataset.theme = saved;
  } catch { /* ignore */ }

  const app = createApp();
  app.wire();
  app.scaffoldPanel = initScaffoldPanel(app);

  try {
    const status = await app.refreshStatus();
    if (status.loaded) {
      await app.reloadGraph({ relayout: true, keepPositions: false });
      await app.loadReport();
      if (status.reference) {
        const data = await api.report().catch(() => null);
        if (data?.reference) app.renderReferenceResults(data.reference);
      }
      app.setStatus('Ready');
    } else {
      app.setStatus('Open an assembly to begin');
      const hint = $('empty-hint');
      if (hint) hint.hidden = false;
    }
    if (status.align_backend === 'none') {
      app.toast(
        'minimap2/mappy not found: reference alignment, misassembly detection and ' +
        'reference scaffolding are unavailable. Install the conda environment to enable them.',
        'warn', 12000,
      );
    }
  } catch (err) {
    app.showError(err);
  }

  window.plastr = app; // handy for debugging from the console
  return app;
}

export default start;
