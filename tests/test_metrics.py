"""Known-answer tests for the reference-free metrics.

The worked example throughout is ``[100, 200, 300, 400]`` (total 1000 bp), for
which every statistic can be checked by hand:

    sorted descending    400  300  200  100
    cumulative           400  700  900  1000
    N50 -> first length whose cumulative reaches 500  -> 300 (L50 = 2)
    N75 -> first length whose cumulative reaches 750  -> 200 (L75 = 3)
    auN  = (400^2 + 300^2 + 200^2 + 100^2) / 1000     -> 300.0
"""

from __future__ import annotations

import pytest

from plastr.core.analysis.metrics import (
    AssemblyMetrics,
    auN,
    compare_metrics,
    compute_metrics,
    nx_stat,
)
from plastr.core.model import AssemblyGraph, Link, Segment

LENGTHS = [100, 200, 300, 400]


@pytest.fixture
def four_contig_graph() -> AssemblyGraph:
    """Four unlinked, sequence-free segments of 100/200/300/400 bp."""
    graph = AssemblyGraph("four")
    for length in LENGTHS:
        graph.add_segment(Segment(name=f"s{length}", sequence=None, length=length))
    return graph


# ---------------------------------------------------------------------------
# nx_stat / auN by hand
# ---------------------------------------------------------------------------


class TestNxStat:
    def test_n50_and_l50(self):
        assert nx_stat(LENGTHS, 0.5) == (300, 2)

    def test_n75_and_l75(self):
        assert nx_stat(LENGTHS, 0.75) == (200, 3)

    def test_n100_is_the_smallest_contig(self):
        assert nx_stat(LENGTHS, 1.0) == (100, 4)

    def test_input_order_does_not_matter(self):
        assert nx_stat([300, 100, 400, 200], 0.5) == nx_stat(LENGTHS, 0.5)

    def test_single_contig(self):
        assert nx_stat([1234], 0.5) == (1234, 1)

    def test_empty_input_returns_zeros(self):
        assert nx_stat([], 0.5) == (0, 0)
        assert nx_stat([], 0.75) == (0, 0)

    def test_ng50_with_an_explicit_genome_size(self):
        # Genome size 1000 == the assembly total, so NG50 == N50.
        assert nx_stat(LENGTHS, 0.5, total=1000) == (300, 2)
        # Genome size 2000: the target is 1000, only reached by the last contig.
        assert nx_stat(LENGTHS, 0.5, total=2000) == (100, 4)
        # Genome size 800: the target is 400, reached immediately.
        assert nx_stat(LENGTHS, 0.5, total=800) == (400, 1)

    def test_ng50_is_zero_when_the_assembly_cannot_reach_the_target(self):
        # 1000 bp assembled of a 4 Mb genome: half of it is never covered.
        assert nx_stat(LENGTHS, 0.5, total=4000) == (0, 0)


class TestAuN:
    def test_known_value(self):
        assert auN(LENGTHS) == pytest.approx(300.0)

    def test_equals_the_contig_length_for_a_single_contig(self):
        assert auN([777]) == pytest.approx(777.0)

    def test_is_at_least_the_mean_and_at_most_the_maximum(self):
        assert max(LENGTHS) >= auN(LENGTHS) >= sum(LENGTHS) / len(LENGTHS)

    def test_with_an_explicit_genome_size(self):
        assert auN(LENGTHS, total=2000) == pytest.approx(150.0)

    def test_empty_and_degenerate_inputs(self):
        assert auN([]) == 0.0
        assert auN(LENGTHS, total=0) == 0.0
        assert auN(LENGTHS, total=-5) == 0.0


# ---------------------------------------------------------------------------
# compute_metrics
# ---------------------------------------------------------------------------


class TestComputeMetrics:
    def test_all_length_statistics(self, four_contig_graph):
        m = compute_metrics(four_contig_graph)
        assert m.num_contigs == 4
        assert m.total_length == 1000
        assert m.largest_contig == 400
        assert m.smallest_contig == 100
        assert m.mean_contig == pytest.approx(250.0)
        assert m.median_contig == 250  # (200 + 300) // 2
        assert (m.n50, m.l50) == (300, 2)
        assert (m.n75, m.l75) == (200, 3)
        assert m.auN == pytest.approx(300.0)
        assert m.num_contigs_ge_1kb == 0
        assert m.total_length_ge_1kb == 0

    def test_size_buckets(self):
        graph = AssemblyGraph()
        for length in (999, 1000, 9999, 10_000, 49_999, 50_000):
            graph.add_segment(Segment(f"s{length}", None, length))
        m = compute_metrics(graph)
        assert m.num_contigs_ge_1kb == 5
        assert m.num_contigs_ge_10kb == 3
        assert m.num_contigs_ge_50kb == 1
        assert m.total_length_ge_1kb == 1000 + 9999 + 10_000 + 49_999 + 50_000

    def test_ng50_needs_a_genome_size(self, four_contig_graph):
        without = compute_metrics(four_contig_graph)
        assert without.ng50 is None and without.lg50 is None and without.genome_size is None

        with_size = compute_metrics(four_contig_graph, genome_size=2000)
        assert with_size.genome_size == 2000
        assert (with_size.ng50, with_size.lg50) == (100, 4)
        # NG75 needs 1500 bp of the 2000 bp genome, which 1000 bp cannot reach.
        assert with_size.ng75 == 0

    def test_ng50_equals_n50_when_the_genome_size_matches(self, four_contig_graph):
        m = compute_metrics(four_contig_graph, genome_size=1000)
        assert m.ng50 == m.n50 == 300
        assert m.lg50 == m.l50 == 2

    def test_min_length_filters_like_quast_min_contig(self, four_contig_graph):
        m = compute_metrics(four_contig_graph, min_length=250)
        assert m.num_contigs == 2
        assert m.total_length == 700
        assert (m.n50, m.l50) == (400, 1)

    def test_empty_graph_returns_zeros_without_raising(self):
        m = compute_metrics(AssemblyGraph())
        assert isinstance(m, AssemblyMetrics)
        assert m.num_contigs == 0
        assert m.total_length == 0
        assert (m.n50, m.l50, m.n75, m.l75) == (0, 0, 0, 0)
        assert m.auN == 0.0
        assert m.largest_contig == 0 and m.smallest_contig == 0
        assert m.mean_contig == 0.0 and m.median_contig == 0
        assert m.gc_percent is None and m.mean_depth is None
        assert m.nx_curve == [] and m.cumulative_curve == [] and m.length_histogram == []
        assert m.to_dict()["num_contigs"] == 0

    def test_everything_filtered_out_returns_zeros(self, four_contig_graph):
        m = compute_metrics(four_contig_graph, min_length=100_000)
        assert m.num_contigs == 0
        assert m.total_length == 0
        assert m.auN == 0.0

    def test_gc_and_n_content_are_length_weighted(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("gc", "GGGGAAAA"))  # 8 bp, 50% GC
        graph.add_segment(Segment("ns", "NNNNNNNN"))  # 8 bp, 0% GC, 8 Ns
        m = compute_metrics(graph)
        assert m.gc_percent == pytest.approx(25.0)  # 100 * (0.5*8 + 0*8) / 16
        assert m.n_per_100kb == pytest.approx(50_000.0)  # 100000 * 8 / 16

    def test_depth_mean_and_cv(self):
        graph = AssemblyGraph()
        for i, depth in enumerate((10.0, 20.0, 30.0)):
            graph.add_segment(Segment(f"s{i}", "ACGT" * 25, depth=depth))
        m = compute_metrics(graph)
        assert m.mean_depth == pytest.approx(20.0)
        # sample stdev of (10,20,30) is 10 -> CV = 0.5
        assert m.depth_cv == pytest.approx(0.5)

    def test_graph_statistics(self, tiny_graph):
        m = compute_metrics(tiny_graph)
        assert m.num_links == 3
        assert m.num_components == 3
        assert m.dead_ends == 4
        assert m.num_circular == 1
        assert m.total_length == 310
        assert m.num_contigs == 5

    def test_demo_graph_metrics(self, demo_graph):
        m = compute_metrics(demo_graph, genome_size=169_000)
        assert m.num_contigs == 10
        assert m.total_length == 137_400
        assert m.num_links == 9
        assert m.num_components == 2
        assert 40.0 < m.gc_percent < 60.0  # random ACGT sequence
        assert m.n_per_100kb == 0.0
        assert m.n50 >= m.n75
        assert m.ng50 is not None and m.ng50 <= m.n50


# ---------------------------------------------------------------------------
# Plot series
# ---------------------------------------------------------------------------


class TestCurves:
    def test_cumulative_curve_is_exact_and_monotonic(self, four_contig_graph):
        curve = compute_metrics(four_contig_graph).cumulative_curve
        assert curve == [(1, 400), (2, 700), (3, 900), (4, 1000)]
        xs = [x for x, _ in curve]
        ys = [y for _, y in curve]
        assert xs == sorted(xs) and len(set(xs)) == len(xs)
        assert all(b > a for a, b in zip(ys, ys[1:]))
        assert ys[-1] == 1000

    def test_nx_curve_shape(self, four_contig_graph):
        curve = compute_metrics(four_contig_graph).nx_curve
        assert len(curve) == 100
        assert [x for x, _ in curve] == list(range(1, 101))
        ys = [y for _, y in curve]
        # Nx never increases as x grows.
        assert all(b <= a for a, b in zip(ys, ys[1:]))
        assert ys[0] == 400  # N1 is the largest contig
        assert ys[-1] == 100  # N100 is the smallest
        assert curve[49] == (50, 300)  # N50 agrees with nx_stat

    def test_nx_curve_midpoint_matches_n50_on_the_demo_graph(self, demo_graph):
        m = compute_metrics(demo_graph)
        assert dict(m.nx_curve)[50] == m.n50

    def test_length_histogram_covers_every_contig(self, demo_graph):
        m = compute_metrics(demo_graph)
        assert sum(count for _edge, count in m.length_histogram) == m.num_contigs
        edges = [edge for edge, _ in m.length_histogram]
        assert edges == sorted(edges)

    def test_histogram_of_identical_lengths_is_a_single_bin(self):
        graph = AssemblyGraph()
        for i in range(4):
            graph.add_segment(Segment(f"s{i}", None, 500))
        assert compute_metrics(graph).length_histogram == [(500, 4)]

    def test_cumulative_curve_is_downsampled_for_huge_assemblies(self):
        graph = AssemblyGraph()
        for i in range(5000):
            graph.add_segment(Segment(f"s{i}", None, 100 + i))
        curve = compute_metrics(graph).cumulative_curve
        assert len(curve) <= 2001
        assert curve[-1][0] == 5000  # the final point is always kept


# ---------------------------------------------------------------------------
# compare_metrics
# ---------------------------------------------------------------------------


class TestCompareMetrics:
    def test_side_by_side_table(self, four_contig_graph):
        a = compute_metrics(four_contig_graph)
        b = compute_metrics(four_contig_graph, min_length=250)
        table = compare_metrics({"all": a, "big": b})
        assert table["columns"] == ["all", "big"]
        rows = {r["label"]: r["values"] for r in table["rows"]}
        assert rows["# contigs"] == ["4", "2"]
        assert rows["N50"] == ["300", "400"]
        assert rows["Total length"] == ["1,000", "700"]
        # GC is unknown for sequence-free segments, so that row is dropped.
        assert "GC (%)" not in rows


class TestWeightedMedianDepth:
    """Median depth is length-weighted, which is what Bandage reports.

    A plain median over segments is dominated by the many tiny nodes in a real
    graph, and the mean is dragged upwards by collapsed repeats. Validated
    against Bandage on two real assemblies (15.4764 and 49.5349).
    """

    @staticmethod
    def _graph(pairs):
        from plastr.core.model import AssemblyGraph, Segment

        g = AssemblyGraph("depth")
        for i, (depth, length) in enumerate(pairs):
            g.add_segment(Segment(f"s{i}", "A" * length, length, depth))
        return g

    def test_many_short_nodes_do_not_dominate(self):
        from plastr.core.analysis.metrics import compute_metrics

        # Nine 100 bp nodes at 100x, one 100 kb node at 10x. Most *bases* are
        # at 10x, so that is the median; a plain median over nodes would say 100.
        pairs = [(100.0, 100)] * 9 + [(10.0, 100_000)]
        m = compute_metrics(self._graph(pairs))
        assert m.median_depth == 10.0
        assert m.mean_depth == pytest.approx(91.0)

    def test_a_single_segment_is_its_own_median(self):
        from plastr.core.analysis.metrics import compute_metrics

        m = compute_metrics(self._graph([(33.5, 5000)]))
        assert m.median_depth == 33.5

    def test_no_depth_information_gives_none(self):
        from plastr.core.analysis.metrics import compute_metrics
        from plastr.core.model import AssemblyGraph, Segment

        g = AssemblyGraph("nodepth")
        g.add_segment(Segment("a", "ACGT" * 100))
        m = compute_metrics(g)
        assert m.median_depth is None
        assert m.mean_depth is None


class TestFilteredTopology:
    """With --min-contig set, topology describes the same segments as the lengths.

    Reporting 63 contigs next to 125 components reads as a bug.
    """

    @staticmethod
    def _chain():
        from plastr.core.model import AssemblyGraph, Link, Segment

        g = AssemblyGraph("chain")
        # big -- tiny -- big : filtering the tiny node splits the component
        g.add_segment(Segment("big1", "A" * 5000))
        g.add_segment(Segment("tiny", "C" * 100))
        g.add_segment(Segment("big2", "G" * 5000))
        g.add_link(Link("big1", "+", "tiny", "+"))
        g.add_link(Link("tiny", "+", "big2", "+"))
        return g

    def test_unfiltered_topology_covers_the_whole_graph(self):
        from plastr.core.analysis.metrics import compute_metrics

        m = compute_metrics(self._chain())
        assert m.num_contigs == 3
        assert m.num_links == 2
        assert m.num_components == 1
        assert m.dead_ends == 2  # the two outer ends

    def test_filtering_out_a_connector_splits_the_component(self):
        from plastr.core.analysis.metrics import compute_metrics

        m = compute_metrics(self._chain(), min_length=1000)
        assert m.num_contigs == 2
        assert m.num_links == 0, "a link to a filtered-out segment is not a connection"
        assert m.num_components == 2
        assert m.dead_ends == 4  # both ends of both surviving segments
