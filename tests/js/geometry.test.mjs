/**
 * The drawing's geometry contract.
 *
 * Every assertion here stands for a bug that shipped: the scale that was never
 * calibrated, the polylines that could never bend, and the links that were
 * longer than the contigs they joined. They are written as properties of the
 * model rather than as fixed numbers, so tuning the constants is allowed and
 * breaking the contract is not.
 */
import test from 'node:test';
import assert from 'node:assert/strict';

import {
  DEFAULT_GEOM, GraphModel, calibrateScale, drawLengthFor, particleCountFor,
  endParticle, innerParticle, END_START, END_END,
} from '../../src/plastr/web/js/graph.js';

/** A chain of `n` contigs of the given lengths, each linked to the next. */
function chainPayload(lengths, { depths = null, circular = false } = {}) {
  const segments = lengths.map((L, i) => ({
    name: 'ctg' + i,
    length: L,
    depth: depths ? depths[i] : 30,
    gc: 0.5,
    component: 0,
    deg_start: i === 0 ? 0 : 1,
    deg_end: i === lengths.length - 1 ? 0 : 1,
    circular: false,
    ref_hits: [],
  }));
  const links = [];
  for (let i = 0; i + 1 < lengths.length; i++) {
    links.push({ from: 'ctg' + i, to: 'ctg' + (i + 1), from_orient: '+', to_orient: '+', overlap: 0 });
  }
  if (circular && lengths.length) {
    links.push({ from: 'ctg0', to: 'ctg0', from_orient: '+', to_orient: '+', overlap: 0 });
  }
  return { segments, links };
}

const meanDrawLen = (g) =>
  g.segments.reduce((a, s) => a + s.drawLen, 0) / (g.segments.length || 1);

test('calibration puts the mean contig near the target length at any assembly size', () => {
  // A fixed units-per-megabase cannot do this, which is what made every contig
  // in a real assembly collapse onto the minimum length together.
  const target = DEFAULT_GEOM.meanNodeLength;
  for (const [n, each] of [[10, 4_000], [300, 5_000], [500, 10_000], [50, 100_000]]) {
    const g = new GraphModel();
    g.setData(chainPayload(new Array(n).fill(each)), { keepPositions: false });
    const mean = meanDrawLen(g);
    assert.ok(
      mean > target * 0.5 && mean < target * 2,
      `mean drawn length ${mean.toFixed(1)} is not near ${target} for ${n} x ${each} bp`,
    );
  }
});

test('calibration is independent of the units-per-megabase fallback', () => {
  const payload = chainPayload(new Array(200).fill(6_000));
  const a = new GraphModel();
  const b = new GraphModel();
  b.geom = { ...DEFAULT_GEOM, unitsPerMegabase: 999_999 };
  a.setData(payload, { keepPositions: false });
  b.setData(payload, { keepPositions: false });
  assert.equal(meanDrawLen(a).toFixed(6), meanDrawLen(b).toFixed(6));
});

test('drawn length keeps the spread of the underlying contig lengths', () => {
  // The floor must not swallow the distribution: if it does, the picture stops
  // carrying length and every contig reads as the same stub.
  const lengths = [500, 1_000, 5_000, 20_000, 80_000, 300_000];
  const g = new GraphModel();
  g.setData(chainPayload(lengths), { keepPositions: false });
  const drawn = g.segments.map((s) => s.drawLen);
  const ratio = Math.max(...drawn) / Math.min(...drawn);
  assert.ok(ratio > 20, `drawn lengths span only ${ratio.toFixed(1)}x for a 600x length span`);
  // Monotonic in length, so a longer contig is never drawn shorter.
  for (let i = 1; i < drawn.length; i++) assert.ok(drawn[i] >= drawn[i - 1]);
});

test('a long contig earns enough vertices to bend', () => {
  const g = new GraphModel();
  g.setData(chainPayload([2_000, 2_000, 2_000, 200_000]), { keepPositions: false });
  const longest = g.segments[g.segments.length - 1];
  assert.ok(
    longest.k >= 5,
    `a 100x-mean contig got ${longest.k} vertices; it cannot curve with fewer than a handful`,
  );
});

test('particleCountFor never returns a degenerate polyline', () => {
  for (const len of [0, 1, 5, 50, 5_000, 1e6]) {
    const k = particleCountFor(len, DEFAULT_GEOM);
    assert.ok(Number.isInteger(k) && k >= 2, `k=${k} for drawLen=${len}`);
    assert.ok(k <= DEFAULT_GEOM.maxParticles);
  }
});

test('calibrateScale is finite on degenerate input', () => {
  for (const [n, total] of [[0, 0], [1, 0], [0, 1000], [5, -1]]) {
    const s = calibrateScale(n, total, DEFAULT_GEOM);
    assert.ok(Number.isFinite(s) && s > 0, `scale ${s} for n=${n} total=${total}`);
  }
  assert.ok(Number.isFinite(drawLengthFor(NaN, DEFAULT_GEOM, 1e-3)));
});

test('each end of a segment exposes an inward neighbour for the edge tangent', () => {
  // Without a distinct inward vertex an edge has no direction to continue, and
  // links fall back to straight chords.
  const g = new GraphModel();
  g.setData(chainPayload([5_000, 50_000]), { keepPositions: false });
  for (const seg of g.segments) {
    for (const end of [END_START, END_END]) {
      const e = endParticle(seg, end);
      const i = innerParticle(seg, end);
      assert.ok(e >= seg.p0 && e < seg.p0 + seg.k);
      assert.ok(i >= seg.p0 && i < seg.p0 + seg.k);
      if (seg.k >= 2) assert.notEqual(e, i, 'inward vertex must differ from the end vertex');
    }
  }
});

test('links carry their tangent vertices and self-loops are marked', () => {
  const g = new GraphModel();
  g.setData(chainPayload([9_000, 9_000], { circular: true }), { keepPositions: false });
  assert.ok(g.links.length >= 2);
  for (const l of g.links) {
    assert.equal(typeof l.paIn, 'number');
    assert.equal(typeof l.pbIn, 'number');
    assert.equal(typeof l.selfLoop, 'boolean');
  }
  assert.equal(g.links.filter((l) => l.selfLoop).length, 1);
});

test('particles are allocated contiguously and map back to their segment', () => {
  const g = new GraphModel();
  g.setData(chainPayload([1_000, 40_000, 3_000, 120_000]), { keepPositions: false });
  let expected = 0;
  for (const seg of g.segments) {
    assert.equal(seg.p0, expected, 'segments must tile the particle array without gaps');
    expected += seg.k;
    for (let i = 0; i < seg.k; i++) assert.equal(g.particleSeg[seg.p0 + i], seg.idx);
  }
  assert.equal(g.nParticles, expected);
  assert.equal(g.px.length, expected);
});

test('rescale keeps the drawing instead of collapsing it to the origin', () => {
  // The node-length slider calls this without re-seeding; zeroing the positions
  // would drop the whole graph onto the origin until the next layout run.
  const g = new GraphModel();
  g.setData(chainPayload([4_000, 40_000, 9_000]), { keepPositions: false });
  g.seedPositions();
  const before = g.segments.map((s) => [g.px[s.p0], g.py[s.p0]]);
  const spreadBefore = Math.max(...before.map(([x, y]) => Math.hypot(x, y)));
  assert.ok(spreadBefore > 0, 'seeding should place segments away from the origin');

  g.rescale(2.5);

  let moved = 0;
  g.segments.forEach((s, i) => {
    assert.ok(Number.isFinite(g.px[s.p0]) && Number.isFinite(g.py[s.p0]));
    moved = Math.max(moved, Math.hypot(g.px[s.p0] - before[i][0], g.py[s.p0] - before[i][1]));
  });
  const spreadAfter = Math.max(...g.segments.map((s) => Math.hypot(g.px[s.p0], g.py[s.p0])));
  assert.ok(spreadAfter > spreadBefore * 0.25, 'rescale must not collapse the layout');
  // And the vertex counts must follow the new lengths.
  for (const seg of g.segments) assert.equal(seg.k, particleCountFor(seg.drawLen, g.geom));
});

test('an empty graph is handled without throwing', () => {
  const g = new GraphModel();
  g.setData({ segments: [], links: [] }, { keepPositions: false });
  assert.equal(g.segments.length, 0);
  assert.equal(g.nParticles, 0);
  g.seedPositions();
  g.rescale(2);
});
