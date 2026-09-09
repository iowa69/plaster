/**
 * The layout's contract.
 *
 * These are the properties that make a drawing readable rather than a hairball:
 * joined contigs end up adjacent, separate components do not sit on top of one
 * another, and a circular contig does not stop the simulation settling.
 */
import test from 'node:test';
import assert from 'node:assert/strict';

import { GraphModel } from '../../src/plastr/web/js/graph.js';
import { LayoutEngine } from '../../src/plastr/web/js/layout.js';

function build(payload) {
  const g = new GraphModel();
  g.setData(payload, { keepPositions: false });
  g.seedPositions();
  const eng = new LayoutEngine();
  const d = g.toLayoutArrays();
  d.px = g.px;
  d.py = g.py;
  eng.init(d);
  return { g, eng };
}

function run(eng, iterations = null) {
  eng.begin('force', null, {});
  const limit = iterations || eng.maxIter;
  for (let i = 0; i < limit; i++) if (eng._stepForce().done) break;
}

/** `comps` chains, each a separate connected component. */
function payload(comps, { circular = [] } = {}) {
  const segments = [];
  const links = [];
  comps.forEach((lengths, c) => {
    lengths.forEach((L, i) => {
      segments.push({
        name: `c${c}_${i}`, length: L, depth: 30, gc: 0.5, component: c,
        deg_start: 1, deg_end: 1, circular: false, ref_hits: [],
      });
      if (i) {
        links.push({
          from: `c${c}_${i - 1}`, to: `c${c}_${i}`,
          from_orient: '+', to_orient: '+', overlap: 0,
        });
      }
    });
    if (circular.includes(c)) {
      links.push({
        from: `c${c}_0`, to: `c${c}_0`, from_orient: '+', to_orient: '+', overlap: 0,
      });
    }
  });
  return { segments, links };
}

/** Median distance between the two ends a link joins. */
function linkGaps(g, eng) {
  const gaps = [];
  for (let l = 0; l < eng.nLinks; l++) {
    if (eng.linkFrom[l] === eng.linkTo[l]) continue;
    const a = eng.particleOf(eng.linkFrom[l], eng.linkFromEnd[l]);
    const b = eng.particleOf(eng.linkTo[l], eng.linkToEnd[l]);
    gaps.push(Math.hypot(eng.px[b] - eng.px[a], eng.py[b] - eng.py[a]));
  }
  gaps.sort((x, y) => x - y);
  return gaps;
}

test('joined contigs come to rest end to end, not on long leaders', () => {
  // The regression this guards: links used to rest at an absolute 24 units
  // against 5-unit contigs, so the drawing was dots joined by long lines.
  for (const chain of [new Array(12).fill(8_000), new Array(120).fill(5_000)]) {
    const { g, eng } = build(payload([chain]));
    run(eng);
    const gaps = linkGaps(g, eng);
    const median = gaps[gaps.length >> 1];
    const meanNode = g.segments.reduce((a, s) => a + s.drawLen, 0) / g.segments.length;
    assert.ok(
      median < meanNode,
      `median link gap ${median.toFixed(1)} exceeds the mean contig length ${meanNode.toFixed(1)}`,
    );
  }
});

test('separate components do not overlap once the layout settles', () => {
  const { g, eng } = build(payload([
    new Array(40).fill(9_000),
    new Array(8).fill(6_000),
    new Array(3).fill(4_000),
  ]));
  run(eng);

  const boxes = new Map();
  for (let i = 0; i < eng.nParticles; i++) {
    const c = g.particleSeg[i] >= 0 ? g.segments[g.particleSeg[i]].component : 0;
    const b = boxes.get(c) || { x0: Infinity, y0: Infinity, x1: -Infinity, y1: -Infinity };
    b.x0 = Math.min(b.x0, eng.px[i]); b.x1 = Math.max(b.x1, eng.px[i]);
    b.y0 = Math.min(b.y0, eng.py[i]); b.y1 = Math.max(b.y1, eng.py[i]);
    boxes.set(c, b);
  }
  const list = [...boxes.values()];
  for (let i = 0; i < list.length; i++) {
    for (let j = i + 1; j < list.length; j++) {
      const a = list[i], b = list[j];
      const overlaps = a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1;
      assert.ok(!overlaps, `component bounding boxes ${i} and ${j} overlap`);
    }
  }
});

test('a circular contig does not destabilise the layout', () => {
  // A segment's self-link used to be given a spring that fought the polyline's
  // own distance constraints, so the pair sat at full speed forever.
  const { eng } = build(payload([new Array(10).fill(7_000)], { circular: [0] }));
  run(eng);
  for (let i = 0; i < eng.nParticles; i++) {
    assert.ok(Number.isFinite(eng.px[i]) && Number.isFinite(eng.py[i]), 'position went non-finite');
    assert.ok(Math.abs(eng.vx[i]) < eng.unit, 'a particle is still moving at layout scale per step');
  }
});

test('positions stay finite for degenerate graphs', () => {
  for (const p of [
    payload([[5_000]]),
    payload([[1_000, 1_000]]),
    { segments: [], links: [] },
  ]) {
    const { eng } = build(p);
    if (!eng.ready) continue;
    run(eng, 120);
    for (let i = 0; i < eng.nParticles; i++) {
      assert.ok(Number.isFinite(eng.px[i]) && Number.isFinite(eng.py[i]));
    }
  }
});

test('a link never rests further apart than the contigs it joins are long', () => {
  const { g, eng } = build(payload([[200_000, 200_000, 200_000, 200_000]]));
  run(eng);
  const gaps = linkGaps(g, eng);
  const longest = Math.max(...g.segments.map((s) => s.drawLen));
  for (const gap of gaps) assert.ok(gap < longest, `gap ${gap.toFixed(1)} vs contig ${longest.toFixed(1)}`);
});

test('a contig that closes on itself is drawn as a ring', () => {
  // A complete circular plasmid assembles into one contig whose ends join. The
  // whole point of looking at it is to see that it closes, so it has to come
  // out round rather than as a straight bar with a loop tacked on.
  const segments = [{
    name: 'plasmid', length: 40_000, depth: 30, gc: 0.5, component: 0,
    deg_start: 1, deg_end: 1, circular: true, ref_hits: [],
  }];
  const links = [{ from: 'plasmid', to: 'plasmid', from_orient: '+', to_orient: '+', overlap: 0 }];
  const { g, eng } = build({ segments, links });

  const seg = g.segments[0];
  assert.ok(seg.closed, 'a self-linked contig must be marked closed');
  assert.ok(seg.k >= 12, `a ring needs vertices to be round; got ${seg.k}`);

  run(eng);

  const pts = [];
  for (let i = 0; i < seg.k; i++) pts.push([eng.px[seg.p0 + i], eng.py[seg.p0 + i]]);
  const cx = pts.reduce((a, p) => a + p[0], 0) / pts.length;
  const cy = pts.reduce((a, p) => a + p[1], 0) / pts.length;
  const rs = pts.map((p) => Math.hypot(p[0] - cx, p[1] - cy));
  const mean = rs.reduce((a, b) => a + b, 0) / rs.length;
  const cv = Math.sqrt(rs.reduce((a, b) => a + (b - mean) ** 2, 0) / rs.length) / mean;
  assert.ok(cv < 0.15, `not round: radial CV ${cv.toFixed(3)}`);

  // And the ring must close: the seam is as tight as any other join.
  const seam = Math.hypot(
    eng.px[seg.p0] - eng.px[seg.p0 + seg.k - 1],
    eng.py[seg.p0] - eng.py[seg.p0 + seg.k - 1],
  );
  const step = seg.drawLen / seg.k;
  assert.ok(seam < step * 2, `ring left open: seam ${seam.toFixed(1)} vs step ${step.toFixed(1)}`);
});

test('contigs joined nose to tail come out as a circle', () => {
  // The other shape circularity arrives in: a closed molecule assembled in
  // several pieces. A ring of contigs is neutrally stable, so without a shape
  // prior it crumples into a blob and hides the fact that it closes.
  const n = 8;
  const segments = Array.from({ length: n }, (_, i) => ({
    name: 'r' + i, length: 9_000, depth: 30, gc: 0.5, component: 0,
    deg_start: 1, deg_end: 1, circular: false, ref_hits: [],
  }));
  const links = Array.from({ length: n }, (_, i) => ({
    from: 'r' + i, to: 'r' + ((i + 1) % n), from_orient: '+', to_orient: '+', overlap: 0,
  }));
  const { g, eng } = build({ segments, links });
  run(eng);

  const pts = [];
  for (let i = 0; i < eng.nParticles; i++) pts.push([eng.px[i], eng.py[i]]);
  const cx = pts.reduce((a, p) => a + p[0], 0) / pts.length;
  const cy = pts.reduce((a, p) => a + p[1], 0) / pts.length;
  const rs = pts.map((p) => Math.hypot(p[0] - cx, p[1] - cy));
  const mean = rs.reduce((a, b) => a + b, 0) / rs.length;
  const cv = Math.sqrt(rs.reduce((a, b) => a + (b - mean) ** 2, 0) / rs.length) / mean;
  assert.ok(cv < 0.2, `ring of contigs did not stay circular: radial CV ${cv.toFixed(3)}`);
  assert.ok(mean > 0, 'ring collapsed to a point');
});

test('a component that is not a closed molecule is left alone', () => {
  // The ring prior must not reshape ordinary graphs: a branching component has
  // a real shape of its own and inventing a circle for it would be a lie.
  const segments = ['a', 'b', 'c', 'd', 'e'].map((nm) => ({
    name: nm, length: 9_000, depth: 30, gc: 0.5, component: 0,
    deg_start: 1, deg_end: 1, circular: false, ref_hits: [],
  }));
  // A star: every contig hangs off 'a'. Nothing here is a cycle.
  const links = ['b', 'c', 'd', 'e'].map((nm) => ({
    from: 'a', to: nm, from_orient: '+', to_orient: '+', overlap: 0,
  }));
  const { g } = build({ segments, links });
  const walk = g._walkComponent(g.components[0]);
  assert.equal(walk.ring, false, 'a tree must not be treated as a ring');
  for (const seg of g.segments) assert.ok(!seg.closed);
});

test('smoothing the joins does not stretch the drawing into a line', () => {
  // Straightening every join looks like the way to make a chain of contigs
  // flow, and it is a trap: the term has no length limit, so chains pull
  // straight and the whole graph extends into rays. Measured once at full
  // strength on every join, the demo graph went from 952 units across to 6522
  // and from square to a 10:1 streak. The relaxation is thresholded so it only
  // touches real kinks; this pins that it stays bounded.
  const n = 60;
  const segments = Array.from({ length: n }, (_, i) => ({
    name: 'c' + i, length: 6_000, depth: 30, gc: 0.5, component: 0,
    deg_start: 1, deg_end: 1, circular: false, ref_hits: [],
  }));
  const links = Array.from({ length: n - 1 }, (_, i) => ({
    from: 'c' + i, to: 'c' + (i + 1), from_orient: '+', to_orient: '+', overlap: 0,
  }));

  const extentOf = (params) => {
    const { eng } = build({ segments, links });
    eng.begin('force', null, params);
    for (let i = 0; i < eng.maxIter; i++) if (eng._stepForce().done) break;
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    for (let i = 0; i < eng.nParticles; i++) {
      x0 = Math.min(x0, eng.px[i]); x1 = Math.max(x1, eng.px[i]);
      y0 = Math.min(y0, eng.py[i]); y1 = Math.max(y1, eng.py[i]);
    }
    const w = x1 - x0, h = y1 - y0;
    return { extent: Math.max(w, h), aspect: Math.max(w, h) / Math.max(1, Math.min(w, h)) };
  };

  const off = extentOf({ jointStrength: 0 });
  const on = extentOf({});
  assert.ok(on.extent < off.extent * 3,
    `joint relaxation stretched the layout ${(on.extent / off.extent).toFixed(1)}x`);
  assert.ok(on.aspect < 6, `layout collapsed towards a line: aspect ${on.aspect.toFixed(1)}`);
});

test('joint relaxation leaves gentle turns alone', () => {
  // It must only spend force where the eye catches a corner. A chain already
  // laid out straight has nothing to relax, so the term should be inert on it
  // rather than quietly rearranging a layout that was already fine.
  const { g, eng } = build(payload([new Array(10).fill(8_000)]));
  // Lay the chain out perfectly straight, then check a step barely moves it.
  let x = 0;
  for (const seg of g.segments) {
    for (let i = 0; i < seg.k; i++) {
      eng.px[seg.p0 + i] = x + (seg.drawLen / Math.max(1, seg.k - 1)) * i;
      eng.py[seg.p0 + i] = 0;
    }
    x += seg.drawLen + 5;
  }
  eng.begin('force', null, { repulsion: 0, gravity: 0, curvature: 0, linkStrength: 0 });
  const before = Array.from(eng.py);
  eng._stepForce();
  let moved = 0;
  for (let i = 0; i < eng.nParticles; i++) moved = Math.max(moved, Math.abs(eng.py[i] - before[i]));
  assert.ok(moved < 1e-3, `a straight chain was disturbed by ${moved.toFixed(4)}`);
});
