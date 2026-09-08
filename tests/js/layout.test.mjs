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
