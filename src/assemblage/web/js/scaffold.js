/**
 * scaffold.js — the scaffold plan editor.
 *
 * The plan returned by `POST /api/scaffold` is plain data, so the manual edits
 * made here and the automatic method share one representation. Every edit is
 * pushed back with `PUT /api/scaffold`, which returns a fresh `preview` (total
 * length, gap bases, bridged bases, N50) that we show immediately.
 *
 * Supported edits: reorder members within and between scaffolds (HTML5 drag &
 * drop, including dragging unplaced contigs in), flip a member's orientation,
 * change a gap size, remove a member, rename a scaffold, and create a new empty
 * scaffold.
 */

import { api, ApiError } from './api.js';
import { fmtBp, fmtInt, fmtNum } from './graph.js';

/* ------------------------------------------------------------- helpers */

function el(tag, attrs, ...children) {
  const n = document.createElement(tag);
  if (attrs) {
    for (const [k, v] of Object.entries(attrs)) {
      if (v === null || v === undefined || v === false) continue;
      if (k === 'class') n.className = v;
      else if (k === 'text') n.textContent = v;
      else if (k === 'html') n.innerHTML = v;
      else if (k.startsWith('on') && typeof v === 'function') n.addEventListener(k.slice(2), v);
      else n.setAttribute(k, v === true ? '' : String(v));
    }
  }
  for (const c of children) {
    if (c === null || c === undefined || c === false) continue;
    n.appendChild(typeof c === 'string' || typeof c === 'number' ? document.createTextNode(String(c)) : c);
  }
  return n;
}

const $ = (id) => document.getElementById(id);

/** Gap evidence -> tag colour class. graph/adjacent gaps hold real sequence. */
const EVIDENCE_CLASSES = new Set(['graph', 'adjacent', 'reference', 'manual', 'default']);
function evidenceClass(e) {
  const v = String(e || 'default').toLowerCase();
  return EVIDENCE_CLASSES.has(v) ? v : 'default';
}

/* -------------------------------------------------------------- module */

/** Drag payload lives here: dataTransfer.getData() is unreadable on dragover. */
let dragSrc = null;

export function initScaffoldPanel(app) {
  const editorEl = $('scaffold-editor');
  const previewEl = $('sc-preview');
  const notesEl = $('sc-notes');
  const unplacedEl = $('sc-unplaced-list');

  const ctrl = {
    plan: null,
    preview: null,
    /** True while a PUT is in flight, so rapid edits don't race each other. */
    saving: false,
    pending: false,
  };

  /* ------------------------------------------------------ plan helpers */

  function ensurePlan() {
    if (!ctrl.plan) {
      ctrl.plan = {
        method: 'manual', scaffold_count: 0, placed_count: 0,
        unplaced: [], redundant: [], notes: [], scaffolds: [],
      };
    }
    if (!Array.isArray(ctrl.plan.scaffolds)) ctrl.plan.scaffolds = [];
    if (!Array.isArray(ctrl.plan.unplaced)) ctrl.plan.unplaced = [];
    for (const s of ctrl.plan.scaffolds) if (!Array.isArray(s.members)) s.members = [];
    return ctrl.plan;
  }

  /** Refresh derived counters and clear the trailing gap of every scaffold. */
  function normalise() {
    const plan = ensurePlan();
    let placed = 0;
    for (const s of plan.scaffolds) {
      const m = s.members;
      for (let i = 0; i < m.length; i++) {
        m[i].orientation = m[i].orientation === '-' ? '-' : '+';
        const g = Number(m[i].gap_after);
        m[i].gap_after = Number.isFinite(g) ? Math.max(0, Math.round(g)) : 0;
        if (i === m.length - 1) m[i].gap_after = 0;
      }
      placed += m.length;
    }
    plan.scaffold_count = plan.scaffolds.length;
    plan.placed_count = placed;
    return plan;
  }

  /** Approximate scaffold length from the graph, for instant feedback. */
  function scaffoldLength(s) {
    let total = 0;
    let unknown = false;
    for (const m of s.members) {
      const seg = app.graph.segmentByName(m.segment);
      if (seg) total += seg.length; else unknown = true;
      total += Number(m.gap_after) || 0;
    }
    return { total, unknown };
  }

  /* --------------------------------------------------------- rendering */

  function renderPreview() {
    previewEl.textContent = '';
    const p = ctrl.preview;
    if (!p) return;
    const rows = [
      ['Scaffolds', fmtInt(p.scaffold_count)],
      ['Total length', fmtBp(p.total_length)],
      ['Gaps', `${fmtInt(p.num_gaps)} · ${fmtBp(p.gap_bases)} N`],
      ['Bridged from graph', fmtBp(p.bridged_bases)],
      ['N50', fmtBp(p.n50)],
      ['Largest', fmtBp(p.largest)],
    ];
    const table = el('table');
    for (const [k, v] of rows) {
      table.appendChild(el('tr', null, el('th', { text: k }), el('td', { text: v })));
    }
    previewEl.appendChild(el('h4', { text: 'Preview' }));
    previewEl.appendChild(table);
    if (Array.isArray(p.warnings) && p.warnings.length) {
      for (const w of p.warnings) {
        previewEl.appendChild(el('div', { class: 'warnbox', text: String(w) }));
      }
    }
  }

  function renderNotes() {
    notesEl.textContent = '';
    const plan = ctrl.plan;
    if (!plan) return;
    const bits = [];
    if (plan.method) bits.push('method: ' + plan.method);
    if (Number.isFinite(Number(plan.placed_count))) bits.push(fmtInt(plan.placed_count) + ' placed');
    if (Array.isArray(plan.redundant) && plan.redundant.length) {
      bits.push(fmtInt(plan.redundant.length) + ' redundant');
    }
    if (bits.length) notesEl.appendChild(el('p', { class: 'muted small', text: bits.join(' · ') }));
    if (Array.isArray(plan.notes)) {
      for (const n of plan.notes) notesEl.appendChild(el('div', { class: 'warnbox', text: String(n) }));
    }
  }

  function renderUnplaced() {
    unplacedEl.textContent = '';
    const plan = ctrl.plan;
    if (!plan || !plan.unplaced || !plan.unplaced.length) return;
    unplacedEl.appendChild(el('h4', { text: `Unplaced (${plan.unplaced.length})` }));
    const wrap = el('ul', { class: 'sc-members' });
    plan.unplaced.forEach((name, i) => {
      const li = el('li', {
        class: 'sc-member',
        draggable: 'true',
        title: 'Drag into a scaffold to place this contig',
      },
        el('span', { class: 'sc-grip', text: '⠿' }),
        el('span', { class: 'sc-name', text: String(name) }),
        el('span', { class: 'muted small', text: lengthLabel(name) }),
        el('span', null),
        el('span', null));
      li.addEventListener('dragstart', (ev) => {
        dragSrc = { kind: 'unplaced', index: i, name: String(name) };
        li.classList.add('dragging');
        try { ev.dataTransfer.setData('text/plain', String(name)); } catch { /* optional */ }
        ev.dataTransfer.effectAllowed = 'move';
      });
      li.addEventListener('dragend', () => { li.classList.remove('dragging'); dragSrc = null; });
      li.addEventListener('click', () => app.selectSegments([String(name)], false));
      wrap.appendChild(li);
    });
    unplacedEl.appendChild(wrap);
  }

  function lengthLabel(name) {
    const seg = app.graph.segmentByName(name);
    return seg ? fmtBp(seg.length) : '—';
  }

  function render() {
    editorEl.textContent = '';
    const plan = ctrl.plan;
    renderNotes();
    renderPreview();
    renderUnplaced();
    if (!plan) {
      editorEl.appendChild(el('p', { class: 'muted small', text: 'No plan yet — build one above.' }));
      return;
    }
    if (!plan.scaffolds.length) {
      editorEl.appendChild(el('p', { class: 'muted small', text: 'The plan contains no scaffolds.' }));
      return;
    }

    plan.scaffolds.forEach((s, si) => {
      const { total, unknown } = scaffoldLength(s);
      const box = el('div', { class: 'sc-scaffold', 'data-s': si });

      const nameInput = el('input', {
        class: 'sc-sc-name', type: 'text', value: s.name || `scaffold_${si + 1}`,
        title: 'Scaffold name (used in the exported FASTA/AGP)',
      });
      nameInput.addEventListener('change', () => {
        s.name = nameInput.value.trim() || `scaffold_${si + 1}`;
        save();
      });

      const meta = el('span', {
        class: 'sc-sc-meta',
        text: `${s.members.length} · ${fmtBp(total)}${unknown ? '?' : ''}`,
        title: s.reference ? 'reference: ' + s.reference : (s.source || ''),
      });

      const del = el('button', {
        class: 'sc-rm', title: 'Delete this scaffold (members become unplaced)', text: '✕',
        onclick: () => {
          for (const m of s.members) if (!plan.unplaced.includes(m.segment)) plan.unplaced.push(m.segment);
          plan.scaffolds.splice(si, 1);
          save();
        },
      });

      const selectAll = el('button', {
        class: 'btn btn-sm', text: 'select', title: 'Select this scaffold on the canvas',
        onclick: () => app.selectSegments(s.members.map((m) => m.segment), false),
      });

      box.appendChild(el('div', { class: 'sc-sc-head' }, nameInput, meta, selectAll, del));

      const ul = el('ul', { class: 'sc-members', 'data-s': si });
      if (!s.members.length) {
        ul.appendChild(el('li', { class: 'sc-empty', text: 'empty — drag contigs here' }));
      }
      s.members.forEach((m, mi) => ul.appendChild(renderMember(plan, s, si, m, mi)));
      attachListDnD(ul, si);
      box.appendChild(ul);
      editorEl.appendChild(box);
    });
  }

  function renderMember(plan, s, si, m, mi) {
    const isLast = mi === s.members.length - 1;
    const li = el('li', {
      class: 'sc-member', draggable: 'true', 'data-s': si, 'data-m': mi,
      title: memberTitle(m),
    });

    li.appendChild(el('span', { class: 'sc-grip', text: '⠿' }));
    li.appendChild(el('span', { class: 'sc-name', text: String(m.segment) }));

    const orient = el('button', {
      class: 'sc-orient' + (m.orientation === '-' ? ' minus' : ''),
      text: m.orientation === '-' ? '−' : '+',
      title: 'Flip orientation',
    });
    orient.addEventListener('click', (ev) => {
      ev.stopPropagation();
      m.orientation = m.orientation === '-' ? '+' : '-';
      save();
    });
    li.appendChild(orient);

    if (isLast) {
      li.appendChild(el('span', { class: 'muted small', text: lengthLabel(m.segment) }));
    } else {
      const gapIn = el('input', {
        type: 'number', min: '0', step: '10',
        value: String(Number(m.gap_after) || 0),
        title: 'Gap after this contig, in bases',
      });
      gapIn.addEventListener('click', (ev) => ev.stopPropagation());
      gapIn.addEventListener('change', () => {
        const v = Math.max(0, Math.round(Number(gapIn.value) || 0));
        if (v !== Number(m.gap_after)) {
          m.gap_after = v;
          // A hand-edited gap is no longer supported by reference/graph evidence.
          m.gap_evidence = 'manual';
        }
        save();
      });
      const tagText = evidenceLabel(m);
      const tag = el('span', {
        class: 'gap-tag ' + evidenceClass(m.gap_evidence),
        text: tagText,
        title: evidenceTitle(m),
      });
      li.appendChild(el('span', { class: 'sc-gap' }, gapIn, tag));
    }

    li.appendChild(el('button', {
      class: 'sc-rm', text: '✕', title: 'Remove from this scaffold',
      onclick: (ev) => {
        ev.stopPropagation();
        s.members.splice(mi, 1);
        if (!plan.unplaced.includes(m.segment)) plan.unplaced.push(m.segment);
        save();
      },
    }));

    li.addEventListener('click', () => app.selectSegments([String(m.segment)], false));

    li.addEventListener('dragstart', (ev) => {
      dragSrc = { kind: 'member', s: si, m: mi, name: String(m.segment) };
      li.classList.add('dragging');
      try { ev.dataTransfer.setData('text/plain', String(m.segment)); } catch { /* optional */ }
      ev.dataTransfer.effectAllowed = 'move';
    });
    li.addEventListener('dragend', () => { li.classList.remove('dragging'); dragSrc = null; clearDropMarks(); });

    li.addEventListener('dragover', (ev) => {
      if (!dragSrc) return;
      ev.preventDefault();
      ev.dataTransfer.dropEffect = 'move';
      const r = li.getBoundingClientRect();
      const after = (ev.clientY - r.top) > r.height / 2;
      clearDropMarks();
      li.classList.add(after ? 'drop-after' : 'drop-before');
    });
    li.addEventListener('drop', (ev) => {
      if (!dragSrc) return;
      ev.preventDefault();
      ev.stopPropagation();
      const r = li.getBoundingClientRect();
      const after = (ev.clientY - r.top) > r.height / 2;
      clearDropMarks();
      moveInto(si, mi + (after ? 1 : 0));
    });

    return li;
  }

  function evidenceLabel(m) {
    const e = evidenceClass(m.gap_evidence);
    const bridged = Number(m.bridge_sequence_length);
    if ((e === 'graph' || e === 'adjacent') && Number.isFinite(bridged) && bridged > 0) {
      return e + ' ' + fmtBp(bridged);
    }
    return e;
  }

  function evidenceTitle(m) {
    const e = evidenceClass(m.gap_evidence);
    const path = Array.isArray(m.bridge_path) && m.bridge_path.length
      ? '\npath: ' + m.bridge_path.join(' → ') : '';
    switch (e) {
      case 'graph': return 'Gap filled with a real path through the assembly graph' + path;
      case 'adjacent': return 'Contigs are adjacent in the graph — the gap holds real sequence' + path;
      case 'reference': return 'Gap size estimated from the reference coordinates';
      case 'manual': return 'Gap size set by hand';
      default: return 'Default gap size (no evidence)';
    }
  }

  function memberTitle(m) {
    const bits = [String(m.segment)];
    const seg = app.graph.segmentByName(m.segment);
    if (seg) bits.push(fmtBp(seg.length));
    if (m.ref) {
      const st = Number(m.ref_start), en = Number(m.ref_end);
      bits.push(`${m.ref}:${Number.isFinite(st) ? fmtInt(st) : '?'}-${Number.isFinite(en) ? fmtInt(en) : '?'}`);
    }
    if (Number.isFinite(Number(m.identity))) bits.push('id ' + fmtNum(Number(m.identity) * 100, 2) + '%');
    if (m.overlaps_previous) bits.push('overlaps previous');
    return bits.join(' · ');
  }

  function clearDropMarks() {
    for (const n of editorEl.querySelectorAll('.drop-before, .drop-after')) {
      n.classList.remove('drop-before', 'drop-after');
    }
    for (const n of editorEl.querySelectorAll('.drop-target')) n.classList.remove('drop-target');
  }

  /** Allow dropping onto the empty area of a scaffold's member list. */
  function attachListDnD(ul, si) {
    ul.addEventListener('dragover', (ev) => {
      if (!dragSrc) return;
      ev.preventDefault();
      ev.dataTransfer.dropEffect = 'move';
      const box = ul.closest('.sc-scaffold');
      if (box) box.classList.add('drop-target');
    });
    ul.addEventListener('dragleave', () => {
      const box = ul.closest('.sc-scaffold');
      if (box) box.classList.remove('drop-target');
    });
    ul.addEventListener('drop', (ev) => {
      if (!dragSrc) return;
      ev.preventDefault();
      clearDropMarks();
      const plan = ensurePlan();
      const target = plan.scaffolds[si];
      moveInto(si, target ? target.members.length : 0);
    });
  }

  /* ------------------------------------------------------------- edits */

  /** Remove every placement of `name`, discarding any scaffold left empty.
   *  A contig belongs in exactly one place. */
  function dropSegmentFromScaffolds(plan, name) {
    for (const scaffold of plan.scaffolds) {
      for (let i = scaffold.members.length - 1; i >= 0; i--) {
        if (scaffold.members[i].segment === name) scaffold.members.splice(i, 1);
      }
    }
    for (let i = plan.scaffolds.length - 1; i >= 0; i--) {
      if (!plan.scaffolds[i].members.length) plan.scaffolds.splice(i, 1);
    }
  }

  /** Move the dragged item into scaffold `dstS` at position `dstIndex`. */
  function moveInto(dstS, dstIndex) {
    const plan = ensurePlan();
    const dst = plan.scaffolds[dstS];
    if (!dst || !dragSrc) { dragSrc = null; return; }
    let member = null;

    if (dragSrc.kind === 'unplaced') {
      const name = dragSrc.name;
      const at = plan.unplaced.indexOf(name);
      if (at >= 0) plan.unplaced.splice(at, 1);
      // With "include unplaced" on, the server also gives each unplaced contig
      // its own singleton scaffold. Without this the contig ends up in two
      // scaffolds and is written to the FASTA and AGP twice.
      dropSegmentFromScaffolds(plan, name);
      member = { segment: name, orientation: '+', gap_after: 100, gap_evidence: 'manual' };
    } else {
      const src = plan.scaffolds[dragSrc.s];
      if (!src) { dragSrc = null; return; }
      const removed = src.members.splice(dragSrc.m, 1);
      if (!removed.length) { dragSrc = null; return; }
      member = removed[0];
      // Removing an earlier element in the same list shifts the target.
      if (dragSrc.s === dstS && dragSrc.m < dstIndex) dstIndex--;
    }
    // dropSegmentFromScaffolds may have removed scaffolds, so re-find the
    // destination by identity rather than trusting the old index.
    const target = plan.scaffolds.includes(dst) ? dst : plan.scaffolds[dstS];
    if (!target) { dragSrc = null; render(); return; }
    dstIndex = Math.max(0, Math.min(target.members.length, dstIndex));
    target.members.splice(dstIndex, 0, member);
    dragSrc = null;
    save();
  }

  function addScaffold() {
    const plan = ensurePlan();
    let n = plan.scaffolds.length + 1;
    const names = new Set(plan.scaffolds.map((s) => s.name));
    while (names.has('scaffold_' + n)) n++;
    plan.scaffolds.push({ name: 'scaffold_' + n, source: 'manual', reference: null, members: [] });
    save();
  }

  /* -------------------------------------------------------------- I/O */

  async function save() {
    normalise();
    render();
    if (ctrl.saving) { ctrl.pending = true; return; }
    ctrl.saving = true;
    try {
      const res = await api.putScaffold(ctrl.plan);
      if (res && res.plan) ctrl.plan = res.plan;
      ctrl.preview = (res && res.preview) || null;
      normalise();
      render();
      app.status('Scaffold plan updated');
    } catch (e) {
      app.showError(e, 'Could not save the scaffold plan');
    } finally {
      ctrl.saving = false;
      if (ctrl.pending) { ctrl.pending = false; save(); }
    }
  }

  async function build() {
    const body = {
      method: $('sc-method').value,
      min_identity: Number($('sc-minid').value) || 0,
      min_query_coverage: Number($('sc-mincov').value) || 0,
      min_align_length: Number($('sc-minalign').value) || 0,
      min_gap: Number($('sc-mingap').value) || 0,
      fill_gaps_from_graph: $('sc-fillgaps').checked,
      include_unplaced: $('sc-unplaced').checked,
      break_misassemblies_first: $('sc-breakmis').checked,
    };
    app.setBusy(true, 'building scaffolds…');
    try {
      const res = await api.buildScaffold(body);
      ctrl.plan = (res && res.plan) || null;
      ctrl.preview = (res && res.preview) || null;
      if (ctrl.plan) normalise();
      render();
      const n = ctrl.plan ? ctrl.plan.scaffolds.length : 0;
      app.toast(`Built ${n} scaffold${n === 1 ? '' : 's'}`, 'ok');
      app.refreshStatus();
    } catch (e) {
      app.showError(e, 'Scaffolding failed');
    } finally {
      app.setBusy(false);
    }
  }

  async function reload() {
    app.setBusy(true, 'loading plan…');
    try {
      const res = await api.getScaffold();
      // The endpoint returns the plan directly, or null when none exists.
      const plan = res && res.plan ? res.plan : res;
      if (!plan || !Array.isArray(plan.scaffolds)) {
        ctrl.plan = null;
        ctrl.preview = null;
        render();
        app.toast('No scaffold plan on the server yet', 'info');
        return;
      }
      ctrl.plan = plan;
      ctrl.preview = (res && res.preview) || null;
      normalise();
      render();
    } catch (e) {
      if (e instanceof ApiError && e.status === 404) {
        ctrl.plan = null; ctrl.preview = null; render();
        app.toast('No scaffold plan on the server yet', 'info');
      } else {
        app.showError(e, 'Could not fetch the scaffold plan');
      }
    } finally {
      app.setBusy(false);
    }
  }

  /* ------------------------------------------------------------- wiring */

  $('sc-build').addEventListener('click', build);
  $('sc-reload').addEventListener('click', reload);
  $('sc-add').addEventListener('click', addScaffold);
  editorEl.addEventListener('dragend', () => { dragSrc = null; clearDropMarks(); });

  render();

  const panel = {
    build,
    reload,
    render,
    get plan() { return ctrl.plan; },
    /** Called when the graph is replaced: names may no longer exist. */
    onGraphChanged() { render(); },
    clear() { ctrl.plan = null; ctrl.preview = null; render(); },
  };
  app.scaffoldPanel = panel;
  return panel;
}

export default initScaffoldPanel;
