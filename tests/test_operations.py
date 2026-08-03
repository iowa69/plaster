"""Graph editing operations.

The subtle one is :func:`split_segment`: links on the original segment's *start*
end must land on the first piece and links on its *end* end on the last piece,
whichever orientation they were written in.
"""

from __future__ import annotations

import pytest

from plastr.core.analysis.operations import (
    delete_segments,
    filter_graph,
    merge_path,
    restore,
    reverse_segment,
    simplify,
    split_segment,
)
from plastr.core.errors import GraphOperationError
from plastr.core.model import AssemblyGraph, Link, Segment
from plastr.core.sequence import revcomp


def snapshot(graph: AssemblyGraph):
    """Everything an undo must restore."""
    segments = {n: (s.length, s.sequence, s.depth) for n, s in graph.segments.items()}
    links = {
        key: (lk.from_name, lk.from_orient, lk.to_name, lk.to_orient, lk.overlap)
        for key, lk in graph.links.items()
    }
    return segments, links


# ---------------------------------------------------------------------------
# delete / restore
# ---------------------------------------------------------------------------


class TestDeleteAndRestore:
    def test_delete_removes_the_segment_and_its_links(self, tiny_graph):
        record = delete_segments(tiny_graph, ["B"])
        assert record["kind"] == "delete_segments"
        assert record["count"] == 1
        assert [s.name for s in record["segments"]] == ["B"]
        assert len(record["links"]) == 2
        assert "B" not in tiny_graph
        assert tiny_graph.link_count == 1
        assert tiny_graph.successors("A", "+") == []
        assert tiny_graph.predecessors("C", "+") == []

    def test_delete_several_at_once(self, tiny_graph):
        record = delete_segments(tiny_graph, ["A", "D", "R"])
        assert record["count"] == 3
        assert sorted(tiny_graph.segments) == ["B", "C"]
        assert tiny_graph.link_count == 1

    def test_unknown_names_are_skipped_quietly(self, tiny_graph):
        record = delete_segments(tiny_graph, ["ghost", "B"])
        assert record["count"] == 1
        assert len(tiny_graph.segments) == 4

    def test_restore_puts_everything_back(self, tiny_graph):
        before = snapshot(tiny_graph)
        record = delete_segments(tiny_graph, ["B", "R"])
        assert snapshot(tiny_graph) != before
        restore(tiny_graph, record)
        assert snapshot(tiny_graph) == before
        assert tiny_graph.walk_sequence([("A", "+"), ("B", "+"), ("C", "+")])
        assert tiny_graph.is_circular("R") is True

    def test_restore_of_a_split_removes_the_pieces_again(self, tiny_graph):
        before = snapshot(tiny_graph)
        record = split_segment(tiny_graph, "B", [30])
        assert "B_part1" in tiny_graph
        restore(tiny_graph, record)
        assert snapshot(tiny_graph) == before
        assert "B_part1" not in tiny_graph


# ---------------------------------------------------------------------------
# filter
# ---------------------------------------------------------------------------


class TestFilterGraph:
    @pytest.fixture
    def sized(self) -> AssemblyGraph:
        graph = AssemblyGraph()
        for name, length, depth in (
            ("tiny", 100, 2.0),
            ("small", 500, 10.0),
            ("medium", 5_000, 30.0),
            ("large", 50_000, 90.0),
        ):
            graph.add_segment(Segment(name, "A" * length, depth=depth))
        return graph

    def test_min_length(self, sized):
        record = filter_graph(sized, min_length=1_000)
        assert record["count"] == 2
        assert sorted(sized.segments) == ["large", "medium"]

    def test_min_length_keeps_the_boundary_value(self, sized):
        filter_graph(sized, min_length=500)
        assert sorted(sized.segments) == ["large", "medium", "small"]

    def test_depth_window(self, sized):
        filter_graph(sized, min_depth=5.0, max_depth=50.0)
        assert sorted(sized.segments) == ["medium", "small"]

    def test_segments_without_depth_fail_a_min_depth_filter(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("has", "A" * 100, depth=10.0))
        graph.add_segment(Segment("hasnt", "A" * 100))
        filter_graph(graph, min_depth=1.0)
        assert list(graph.segments) == ["has"]

    def test_keep_components_keeps_the_largest(self, tiny_graph):
        # Components are 240 bp (A,B,C), 40 bp (D) and 30 bp (R).
        record = filter_graph(tiny_graph, keep_components=1)
        assert record["count"] == 2
        assert sorted(tiny_graph.segments) == ["A", "B", "C"]
        assert tiny_graph.link_count == 2

    def test_keep_components_two(self, tiny_graph):
        filter_graph(tiny_graph, keep_components=2)
        assert sorted(tiny_graph.segments) == ["A", "B", "C", "D"]

    def test_min_component_length(self, tiny_graph):
        filter_graph(tiny_graph, min_component_length=100)
        assert sorted(tiny_graph.segments) == ["A", "B", "C"]

    def test_min_component_length_drops_everything_when_too_high(self, tiny_graph):
        filter_graph(tiny_graph, min_component_length=10_000)
        assert tiny_graph.segments == {}
        assert tiny_graph.link_count == 0

    def test_a_no_op_filter_changes_nothing(self, tiny_graph):
        before = snapshot(tiny_graph)
        record = filter_graph(tiny_graph)
        assert record["count"] == 0
        assert snapshot(tiny_graph) == before

    def test_filter_is_undoable(self, tiny_graph):
        before = snapshot(tiny_graph)
        record = filter_graph(tiny_graph, min_length=50)
        assert sorted(tiny_graph.segments) == ["A", "B", "C"]
        restore(tiny_graph, record)
        assert snapshot(tiny_graph) == before


# ---------------------------------------------------------------------------
# split_segment
# ---------------------------------------------------------------------------


@pytest.fixture
def split_target() -> AssemblyGraph:
    """X (100 bp) with one external link on each of its four possible spellings.

        A.end   -> X.start     Link(A, +, X, +)
        X.end   -> B.start     Link(X, +, B, +)
        X.start -> C.start     Link(X, -, C, +)
        D.end   -> X.end       Link(D, +, X, -)
    """
    graph = AssemblyGraph()
    graph.add_segment(Segment("X", "".join("ACGT"[i % 4] for i in range(100))))
    for name in ("A", "B", "C", "D"):
        graph.add_segment(Segment(name, "T" * 20))
    graph.add_link(Link("A", "+", "X", "+"))
    graph.add_link(Link("X", "+", "B", "+"))
    graph.add_link(Link("X", "-", "C", "+"))
    graph.add_link(Link("D", "+", "X", "-"))
    return graph


class TestSplitSegment:
    def test_pieces_carry_the_right_sub_sequences(self, split_target):
        original = split_target["X"].sequence
        record = split_segment(split_target, "X", [40])

        assert record["kind"] == "split_segment"
        assert record["added_segments"] == ["X_part1", "X_part2"]
        assert record["count"] == 2
        assert "X" not in split_target
        assert split_target["X_part1"].sequence == original[:40]
        assert split_target["X_part2"].sequence == original[40:]
        assert split_target["X_part1"].length == 40
        assert split_target["X_part2"].length == 60
        assert (
            split_target["X_part1"].sequence + split_target["X_part2"].sequence == original
        )

    def test_consecutive_pieces_are_joined(self, split_target):
        split_segment(split_target, "X", [40])
        assert split_target.has_link("X_part1", "+", "X_part2", "+")

    def test_external_links_reattach_to_the_correct_end_piece(self, split_target):
        split_segment(split_target, "X", [40])

        # A joined X's start, so it must now join the *first* piece's start.
        assert split_target.has_link("A", "+", "X_part1", "+")
        assert not split_target.has_link("A", "+", "X_part2", "+")

        # B was on X's end, so it hangs off the *last* piece.
        assert split_target.has_link("X_part2", "+", "B", "+")
        assert not split_target.has_link("X_part1", "+", "B", "+")

        # (X,-) -> (C,+) leaves X's start: first piece.
        assert split_target.has_link("X_part1", "-", "C", "+")
        assert not split_target.has_link("X_part2", "-", "C", "+")

        # (D,+) -> (X,-) arrives at X's end: last piece.
        assert split_target.has_link("D", "+", "X_part2", "-")
        assert not split_target.has_link("D", "+", "X_part1", "-")

        # Four external links plus the one internal join, and nothing else.
        assert split_target.link_count == 5

    def test_the_ends_of_the_whole_still_behave_like_the_original(self, split_target):
        before_start = sorted(split_target.successors("X", "-"))
        before_end = sorted(split_target.successors("X", "+"))
        split_segment(split_target, "X", [40])
        assert sorted(split_target.successors("X_part1", "-")) == before_start
        assert sorted(split_target.successors("X_part2", "+")) == before_end

    def test_several_cut_points(self, split_target):
        original = split_target["X"].sequence
        record = split_segment(split_target, "X", [25, 50, 75])
        assert record["count"] == 4
        names = ["X_part1", "X_part2", "X_part3", "X_part4"]
        assert record["added_segments"] == names
        assert [split_target[n].length for n in names] == [25, 25, 25, 25]
        assert "".join(split_target[n].sequence for n in names) == original
        for left, right in zip(names, names[1:]):
            assert split_target.has_link(left, "+", right, "+")
        # Outer links still go to the outermost pieces.
        assert split_target.has_link("A", "+", "X_part1", "+")
        assert split_target.has_link("X_part4", "+", "B", "+")

    def test_out_of_range_and_duplicate_positions(self, split_target):
        record = split_segment(split_target, "X", [40, 40, 0, 100, 5_000, -3])
        assert record["count"] == 2  # only 40 is usable

    def test_no_usable_position_raises(self, split_target):
        with pytest.raises(GraphOperationError, match="no usable split position"):
            split_segment(split_target, "X", [0, 100, 500])

    def test_unknown_segment_raises(self, split_target):
        with pytest.raises(GraphOperationError, match="no such segment"):
            split_segment(split_target, "ghost", [10])

    def test_a_segment_without_sequence_cannot_be_split(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("s", None, 1_000))
        with pytest.raises(GraphOperationError, match="cannot be split"):
            split_segment(graph, "s", [500])

    def test_self_loops_are_dropped_rather_than_duplicated(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("R", "ACGT" * 25))
        graph.add_link(Link("R", "+", "R", "+"))
        split_segment(graph, "R", [50])
        assert sorted(graph.segments) == ["R_part1", "R_part2"]
        assert graph.link_count == 1  # just the internal join
        assert graph.has_link("R_part1", "+", "R_part2", "+")

    def test_depth_and_tags_are_inherited(self, split_target):
        split_target["X"].depth = 42.0
        split_target["X"].tags = {"note": "keep me"}
        split_segment(split_target, "X", [40])
        for name in ("X_part1", "X_part2"):
            assert split_target[name].depth == 42.0
            assert split_target[name].tags == {"note": "keep me"}


# ---------------------------------------------------------------------------
# merge_path
# ---------------------------------------------------------------------------


class TestMergePath:
    def test_concatenates_the_sequences_and_rewires_the_neighbours(self, chain_graph):
        seq_a = chain_graph["A"].sequence
        seq_b = chain_graph["B"].sequence
        record = merge_path(chain_graph, [("A", "+"), ("B", "+")])

        assert record["kind"] == "merge_path"
        assert record["count"] == 2
        new_name = record["new_name"]
        assert new_name == "merged_A"
        assert chain_graph[new_name].sequence == seq_a + seq_b
        assert chain_graph[new_name].length == 20
        assert "A" not in chain_graph and "B" not in chain_graph

        # P fed A and C followed B; both must now attach to the merged segment.
        assert chain_graph.has_link("P", "+", new_name, "+")
        assert chain_graph.has_link(new_name, "+", "C", "+")
        assert chain_graph.has_link("C", "+", "Q", "+")
        assert chain_graph.link_count == 3

    def test_overlaps_are_trimmed_when_merging(self, tiny_graph):
        record = merge_path(tiny_graph, [("A", "+"), ("B", "+"), ("C", "+")])
        merged = tiny_graph[record["new_name"]]
        from conftest import TINY

        assert merged.sequence == TINY.walk_abc
        assert merged.length == TINY.walk_abc_length == 220

    def test_reverse_oriented_steps_are_revcomped(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("a", "AAAACCCCGGGGTTTA"))
        graph.add_segment(Segment("b", "ACGTACGTACGTACGA"))
        graph.add_link(Link("a", "+", "b", "-"))
        record = merge_path(graph, [("a", "+"), ("b", "-")])
        assert graph[record["new_name"]].sequence == "AAAACCCCGGGGTTTA" + revcomp(
            "ACGTACGTACGTACGA"
        )

    def test_depth_is_length_weighted(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("a", "A" * 100, depth=10.0))
        graph.add_segment(Segment("b", "C" * 300, depth=30.0))
        graph.add_link(Link("a", "+", "b", "+"))
        record = merge_path(graph, [("a", "+"), ("b", "+")])
        # (10*100 + 30*300) / 400
        assert graph[record["new_name"]].depth == pytest.approx(25.0)

    def test_an_explicit_name_is_used(self, chain_graph):
        record = merge_path(chain_graph, [("A", "+"), ("B", "+")], new_name="contig_1")
        assert record["new_name"] == "contig_1"
        assert "contig_1" in chain_graph

    def test_a_clashing_name_is_made_unique(self, chain_graph):
        chain_graph.add_segment(Segment("contig_1", "A" * 5))
        record = merge_path(chain_graph, [("A", "+"), ("B", "+")], new_name="contig_1")
        assert record["new_name"] == "contig_1_2"

    def test_too_few_steps_raises(self, chain_graph):
        with pytest.raises(GraphOperationError, match="at least two"):
            merge_path(chain_graph, [("A", "+")])

    def test_unknown_segment_raises(self, chain_graph):
        with pytest.raises(GraphOperationError, match="unknown segment"):
            merge_path(chain_graph, [("A", "+"), ("ghost", "+")])

    def test_total_length_is_preserved_for_overlap_free_merges(self, chain_graph):
        before = chain_graph.total_length
        merge_path(chain_graph, [("A", "+"), ("B", "+"), ("C", "+")])
        assert chain_graph.total_length == before


# ---------------------------------------------------------------------------
# reverse_segment
# ---------------------------------------------------------------------------


class TestReverseSegment:
    def test_sequence_is_revcomped_and_links_flip(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("P", "AAAACCCCGGGGTTTA"))
        graph.add_segment(Segment("X", "ACGTACGTACGTACGA"))
        graph.add_segment(Segment("Q", "TTTTGGGGCCCCAAAT"))
        graph.add_link(Link("P", "+", "X", "+"))
        graph.add_link(Link("X", "+", "Q", "+"))
        original = graph["X"].sequence

        record = reverse_segment(graph, "X")
        assert record == {"kind": "reverse_segment", "name": "X", "count": 1}
        assert graph["X"].sequence == revcomp(original)
        assert graph["X"].sequence != original

        # P still joins X's (now flipped) 5' end, expressed as X-.
        assert graph.has_link("P", "+", "X", "-")
        assert not graph.has_link("P", "+", "X", "+")
        assert graph.has_link("X", "-", "Q", "+")
        assert not graph.has_link("X", "+", "Q", "+")
        assert graph.link_count == 2

    def test_reversing_twice_is_a_no_op(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("P", "AAAACCCCGGGGTTTA"))
        graph.add_segment(Segment("X", "ACGTACGTACGTACGA"))
        graph.add_link(Link("P", "+", "X", "+"))
        before = snapshot(graph)
        reverse_segment(graph, "X")
        reverse_segment(graph, "X")
        assert snapshot(graph) == before

    def test_the_walk_through_a_reversed_segment_is_unchanged(self, chain_graph):
        before = chain_graph.walk_sequence([("P", "+"), ("A", "+"), ("B", "+")])
        reverse_segment(chain_graph, "A")
        after = chain_graph.walk_sequence([("P", "+"), ("A", "-"), ("B", "+")])
        assert after == before

    def test_unknown_segment_raises(self, chain_graph):
        with pytest.raises(GraphOperationError, match="no such segment"):
            reverse_segment(chain_graph, "ghost")

    def test_a_segment_without_sequence_raises(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("s", None, 100))
        with pytest.raises(GraphOperationError, match="no sequence"):
            reverse_segment(graph, "s")


# ---------------------------------------------------------------------------
# simplify
# ---------------------------------------------------------------------------


class TestSimplify:
    def test_a_linear_chain_collapses_to_one_segment(self, chain_graph):
        expected = chain_graph.walk_sequence(
            [("P", "+"), ("A", "+"), ("B", "+"), ("C", "+"), ("Q", "+")]
        )
        record = simplify(chain_graph)

        assert record["kind"] == "simplify"
        assert record["count"] == 1
        assert len(chain_graph.segments) == 1
        assert chain_graph.link_count == 0
        (only,) = chain_graph.segments.values()
        assert only.sequence == expected
        assert only.length == 50

    def test_simplify_stops_at_a_branch(self):
        graph = AssemblyGraph()
        for name in ("a", "b", "c", "d"):
            graph.add_segment(Segment(name, "ACGT" * 10))
        graph.add_link(Link("a", "+", "b", "+"))
        graph.add_link(Link("b", "+", "c", "+"))
        graph.add_link(Link("b", "+", "d", "+"))
        record = simplify(graph)
        assert record["count"] == 1
        assert sorted(graph.segments) == ["c", "d", "merged_a"]
        assert graph["merged_a"].length == 80
        assert graph.has_link("merged_a", "+", "c", "+")
        assert graph.has_link("merged_a", "+", "d", "+")

    def test_nothing_to_do(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("lonely", "ACGT" * 10))
        record = simplify(graph)
        assert record["count"] == 0
        assert list(graph.segments) == ["lonely"]

    def test_total_length_is_preserved(self, chain_graph):
        before = chain_graph.total_length
        simplify(chain_graph)
        assert chain_graph.total_length == before

    def test_simplify_respects_overlaps(self, tiny_graph):
        from conftest import TINY

        simplify(tiny_graph)
        merged = [s for s in tiny_graph.segments.values() if s.length == TINY.walk_abc_length]
        assert len(merged) == 1
        assert merged[0].sequence == TINY.walk_abc
        # D and R were not part of any chain.
        assert "D" in tiny_graph and "R" in tiny_graph


class TestMergeSafety:
    """merge_path must not fabricate joins or lose overlaps."""

    @staticmethod
    def _graph(overlap=55):
        import random

        from plastr.core.model import AssemblyGraph, Link, Segment

        rng = random.Random(4)
        seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        g = AssemblyGraph("m")
        a = seq(500)
        b = a[-overlap:] + seq(500 - overlap)
        c = b[-overlap:] + seq(500 - overlap)
        lonely = seq(400)
        for name, s in (("a", a), ("b", b), ("c", c), ("lonely", lonely)):
            g.add_segment(Segment(name, s))
        g.add_link(Link("a", "+", "b", "+", overlap, f"{overlap}M"))
        g.add_link(Link("b", "+", "c", "+", overlap, f"{overlap}M"))
        g.overlap_default = overlap
        return g

    def test_merging_unconnected_segments_is_refused(self):
        from plastr.core.analysis.operations import merge_path
        from plastr.core.errors import GraphOperationError

        g = self._graph()
        # walk_sequence would fall back to overlap_default and silently delete
        # 55 real bases from the second contig.
        with pytest.raises(GraphOperationError, match="not connected"):
            merge_path(g, [("a", "+"), ("lonely", "+")])

    def test_a_boundary_link_keeps_its_overlap_after_a_merge(self):
        from plastr.core.analysis.operations import merge_path
        from plastr.core.model import Link, Segment

        g = self._graph()
        g.add_segment(Segment("d", "A" * 300))
        g.add_link(Link("c", "+", "d", "+", 55, "55M"))
        record = merge_path(g, [("a", "+"), ("b", "+")])
        merged = record["new_name"]
        onward = [lk for lk in g.links.values() if merged in (lk.from_name, lk.to_name)]
        assert onward, "the merge dropped the boundary link entirely"
        assert all(lk.overlap == 55 for lk in onward), (
            "boundary link lost its overlap; the exported GFA would declare a "
            "blunt join and later merges would duplicate bases"
        )

    def test_the_merged_sequence_is_the_walk_not_the_concatenation(self):
        from plastr.core.analysis.operations import merge_path

        g = self._graph()
        expected = g.walk_sequence([("a", "+"), ("b", "+"), ("c", "+")])
        record = merge_path(g, [("a", "+"), ("b", "+"), ("c", "+")])
        assert g.segments[record["new_name"]].sequence == expected
        assert len(expected) == 3 * 500 - 2 * 55


class TestReverseIsUndoable:
    """undo after reverse was a silent no-op presented as success."""

    @staticmethod
    def _project():
        from plastr.core.model import AssemblyGraph, Link, Segment
        from plastr.core.project import Project

        g = AssemblyGraph("r")
        g.add_segment(Segment("a", "AAAACCCGGGTTT"))
        g.add_segment(Segment("b", "GGGGTTTTAAAA"))
        g.add_link(Link("a", "+", "b", "+", 0, "*"))
        p = Project()
        p.graph = g
        return p

    def test_undo_restores_the_sequence_and_the_links(self):
        p = self._project()
        g = p.graph
        before_seq = g.segments["a"].sequence
        before_links = {lk.key() for lk in g.links.values()}

        p.apply_operation("reverse", {"name": "a"})
        assert g.segments["a"].sequence != before_seq

        p.undo()
        assert g.segments["a"].sequence == before_seq
        assert {lk.key() for lk in g.links.values()} == before_links
        assert p.undo_stack == []
