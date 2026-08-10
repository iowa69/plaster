"""Tests for the double-stranded graph model.

The convention (from the module docstring and docs/API.md) is the subtle part of
Plastr, so it is pinned down here in detail: a link and its reverse form are
one physical connection, and the physical ends a link joins follow from the two
orientations.
"""

from __future__ import annotations

import pytest

from plastr.core.errors import GraphOperationError
from plastr.core.model import (
    AssemblyGraph,
    Link,
    Path,
    Segment,
    depth_from_name,
    depth_from_tags,
    length_from_name,
)

from conftest import TINY

# The table in docs/API.md ("Link end convention").
END_TABLE = [
    ("+", "+", ("A", "end"), ("B", "start")),
    ("+", "-", ("A", "end"), ("B", "end")),
    ("-", "+", ("A", "start"), ("B", "start")),
    ("-", "-", ("A", "start"), ("B", "end")),
]


# ---------------------------------------------------------------------------
# Segment
# ---------------------------------------------------------------------------


class TestSegment:
    def test_length_is_derived_from_the_sequence(self):
        assert Segment("s", "ACGTACGT").length == 8

    def test_explicit_length_wins_when_there_is_no_sequence(self):
        seg = Segment("s", None, 1234)
        assert seg.length == 1234
        assert seg.has_sequence is False

    def test_gc_is_none_without_a_sequence(self):
        assert Segment("s", None, 100).gc is None
        assert Segment("s", "GGCC").gc == pytest.approx(1.0)

    def test_seq_oriented(self):
        seg = Segment("s", "ACCTGA")
        assert seg.seq_oriented("+") == "ACCTGA"
        assert seg.seq_oriented("-") == "TCAGGT"

    def test_seq_oriented_without_sequence_raises(self):
        with pytest.raises(GraphOperationError, match="no sequence"):
            Segment("s", None, 10).seq_oriented("+")

    def test_to_dict(self):
        seg = Segment("s", "GGCC", depth=12.5, tags={"xx": 1})
        d = seg.to_dict()
        assert d["name"] == "s" and d["length"] == 4 and d["depth"] == 12.5
        assert d["gc"] == pytest.approx(1.0)
        assert "sequence" not in d
        assert seg.to_dict(include_sequence=True)["sequence"] == "GGCC"


# ---------------------------------------------------------------------------
# Link canonicalisation and the end table
# ---------------------------------------------------------------------------


class TestLinkCanonicalisation:
    def test_a_link_and_its_reverse_share_a_canonical_key(self):
        forward = Link("A", "+", "B", "+")
        backward = Link("B", "-", "A", "-")
        assert forward.reverse() == backward
        assert backward.reverse() == forward
        assert forward.canonical_key() == backward.canonical_key()
        assert forward.key() != backward.key()

    @pytest.mark.parametrize("from_or,to_or", [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")])
    def test_reverse_is_an_involution(self, from_or, to_or):
        link = Link("A", from_or, "B", to_or, 7, "7M")
        assert link.reverse().reverse() == link

    def test_adding_a_link_twice_in_opposite_forms_yields_one_link(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("A", "A" * 10))
        graph.add_segment(Segment("B", "C" * 10))

        assert graph.add_link(Link("A", "+", "B", "+")) is True
        # The same physical connection, written from B's side.
        assert graph.add_link(Link("B", "-", "A", "-")) is False

        assert graph.link_count == 1
        assert len(graph.links) == 1
        # Both spellings are still queryable.
        assert graph.has_link("A", "+", "B", "+")
        assert graph.has_link("B", "-", "A", "-")

    def test_the_opposite_form_of_a_different_link_is_not_deduplicated(self):
        graph = AssemblyGraph()
        for name in "AB":
            graph.add_segment(Segment(name, "A" * 10))
        assert graph.add_link(Link("A", "+", "B", "+")) is True
        # (A,+)->(B,-) is a *different* physical connection (A.end -> B.end).
        assert graph.add_link(Link("A", "+", "B", "-")) is True
        assert graph.link_count == 2

    def test_dangling_links_are_dropped_but_can_be_made_fatal(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("A", "A" * 10))
        assert graph.add_link(Link("A", "+", "ghost", "+")) is False
        assert graph.link_count == 0
        with pytest.raises(GraphOperationError, match="unknown segment"):
            graph.add_link(Link("A", "+", "ghost", "+"), strict=True)


class TestLinkEnds:
    @pytest.mark.parametrize("from_or,to_or,tail,head", END_TABLE)
    def test_ends_match_the_documented_table(self, from_or, to_or, tail, head):
        assert Link("A", from_or, "B", to_or).ends() == (tail, head)

    @pytest.mark.parametrize("from_or,to_or,tail,head", END_TABLE)
    def test_the_reverse_form_joins_the_same_two_ends(self, from_or, to_or, tail, head):
        """A link and its reverse describe the same physical pair of ends."""
        reverse_ends = Link("A", from_or, "B", to_or).reverse().ends()
        assert set(reverse_ends) == {tail, head}


# ---------------------------------------------------------------------------
# Traversal consistency
# ---------------------------------------------------------------------------


def _four_orientation_graph() -> AssemblyGraph:
    """One graph containing all four orientation combinations, on disjoint pairs."""
    graph = AssemblyGraph()
    for name in ("A", "B", "C", "D", "E", "F", "G", "H"):
        graph.add_segment(Segment(name, "ACGT" * 5))
    graph.add_link(Link("A", "+", "B", "+"))
    graph.add_link(Link("C", "+", "D", "-"))
    graph.add_link(Link("E", "-", "F", "+"))
    graph.add_link(Link("G", "-", "H", "-"))
    return graph


class TestSuccessorsPredecessors:
    def test_successor_and_predecessor_are_two_views_of_one_link(self, tiny_graph):
        assert tiny_graph.successors("A", "+") == [("B", "+")]
        assert tiny_graph.predecessors("B", "+") == [("A", "+")]
        # ...and the same statement read on the other strand.
        assert tiny_graph.successors("B", "-") == [("A", "-")]
        assert tiny_graph.predecessors("A", "-") == [("B", "-")]

    def test_consistency_over_every_orientation_combination(self):
        graph = _four_orientation_graph()
        checked = 0
        for name in graph.segments:
            for orient in ("+", "-"):
                for succ in graph.successors(name, orient):
                    assert (name, orient) in graph.predecessors(*succ), (
                        f"{name}{orient} -> {succ} is not mirrored by a predecessor"
                    )
                    assert graph.has_link(name, orient, succ[0], succ[1])
                    checked += 1
                for pred in graph.predecessors(name, orient):
                    assert (name, orient) in graph.successors(*pred)
                    assert graph.has_link(pred[0], pred[1], name, orient)
                    checked += 1
        # 4 links x 2 strands x (successor + predecessor view)
        assert checked == 16

    def test_successors_of_a_reverse_oriented_link(self):
        graph = _four_orientation_graph()
        # (C,+) -> (D,-): joins C.end to D.end.
        assert graph.successors("C", "+") == [("D", "-")]
        assert graph.successors("D", "+") == [("C", "-")]
        assert graph.predecessors("D", "-") == [("C", "+")]

    def test_neighbours_excludes_self(self, tiny_graph):
        assert tiny_graph.neighbours("B") == {"A", "C"}
        assert tiny_graph.neighbours("D") == set()
        assert tiny_graph.neighbours("R") == set()  # self-loop is not a neighbour


class TestDegreeAndDeadEnds:
    def test_degree_is_start_end_ordered(self, tiny_graph):
        assert tiny_graph.degree("A") == (0, 1)  # nothing on the start end
        assert tiny_graph.degree("B") == (1, 1)
        assert tiny_graph.degree("C") == (1, 0)
        assert tiny_graph.degree("D") == (0, 0)
        assert tiny_graph.degree("R") == (1, 1)  # the self loop uses both ends

    def test_is_dead_end(self, tiny_graph):
        assert tiny_graph.is_dead_end("A") is True
        assert tiny_graph.is_dead_end("B") is False
        assert tiny_graph.is_dead_end("C") is True
        assert tiny_graph.is_dead_end("D") is True
        assert tiny_graph.is_dead_end("R") is False

    def test_dead_ends_and_count(self, tiny_graph):
        assert set(tiny_graph.dead_ends()) == TINY.dead_end_names
        # D is isolated and so contributes two dead ends.
        assert tiny_graph.dead_end_count() == TINY.dead_end_count == 4


class TestCircular:
    def test_a_self_looping_segment_is_circular(self, tiny_graph):
        assert tiny_graph.is_circular("R") is True

    def test_a_chain_member_is_not_circular(self, tiny_graph):
        for name in ("A", "B", "C", "D"):
            assert tiny_graph.is_circular(name) is False

    def test_circularity_needs_a_true_self_loop(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("X", "ACGT" * 10))
        graph.add_segment(Segment("Y", "ACGT" * 10))
        # A two-node cycle: neither node is a self-contained circle.
        graph.add_link(Link("X", "+", "Y", "+"))
        graph.add_link(Link("Y", "+", "X", "+"))
        assert graph.is_circular("X") is False
        assert graph.is_circular("Y") is False

    def test_a_self_loop_plus_another_link_is_not_circular(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("X", "ACGT" * 10))
        graph.add_segment(Segment("Y", "ACGT" * 10))
        graph.add_link(Link("X", "+", "X", "+"))
        graph.add_link(Link("X", "+", "Y", "+"))
        assert graph.is_circular("X") is False


# ---------------------------------------------------------------------------
# Components
# ---------------------------------------------------------------------------


class TestConnectedComponents:
    def test_three_components_sorted_largest_first(self, tiny_graph):
        comps = tiny_graph.connected_components()
        assert len(comps) == 3
        assert [sorted(c) for c in comps] == TINY.components
        lengths = [sum(tiny_graph[n].length for n in c) for c in comps]
        assert lengths == TINY.component_lengths == [240, 40, 30]
        assert lengths == sorted(lengths, reverse=True)

    def test_component_map_indexes_from_the_largest(self, tiny_graph):
        mapping = tiny_graph.component_map()
        assert mapping["A"] == mapping["B"] == mapping["C"] == 0
        assert mapping["D"] == 1
        assert mapping["R"] == 2

    def test_every_segment_appears_exactly_once(self, tiny_graph):
        flat = [n for comp in tiny_graph.connected_components() for n in comp]
        assert sorted(flat) == sorted(TINY.names)

    def test_ordering_is_by_total_length_not_by_segment_count(self):
        graph = AssemblyGraph()
        # Component of three tiny segments (30 bp) vs one big one (1000 bp).
        for name in ("t1", "t2", "t3"):
            graph.add_segment(Segment(name, "A" * 10))
        graph.add_link(Link("t1", "+", "t2", "+"))
        graph.add_link(Link("t2", "+", "t3", "+"))
        graph.add_segment(Segment("big", "C" * 1000))
        comps = graph.connected_components()
        assert comps[0] == ["big"]
        assert sorted(comps[1]) == ["t1", "t2", "t3"]


# ---------------------------------------------------------------------------
# Walks
# ---------------------------------------------------------------------------


class TestWalkSequence:
    def test_empty_walk(self, tiny_graph):
        assert tiny_graph.walk_sequence([]) == ""

    def test_single_step_is_the_segment_itself(self, tiny_graph):
        assert tiny_graph.walk_sequence([("A", "+")]) == tiny_graph["A"].sequence

    def test_overlaps_are_trimmed_exactly(self, tiny_graph):
        walk = tiny_graph.walk_sequence([("A", "+"), ("B", "+"), ("C", "+")])
        # 100 + (80 - 10) + (60 - 10)
        assert len(walk) == TINY.walk_abc_length == 220
        assert walk == TINY.walk_abc
        # The shared 10-mers appear once each, not twice.
        assert walk.count("GATTACAGAT") == 1
        assert walk.count("TTCCGGAACC") == 1

    def test_reverse_orientation_revcomps_before_trimming(self, tiny_graph):
        from plastr.core.sequence import revcomp

        # C- -> B- -> A- is the same physical walk read the other way round.
        walk = tiny_graph.walk_sequence([("C", "-"), ("B", "-"), ("A", "-")])
        assert len(walk) == 220
        assert walk == revcomp(TINY.walk_abc)

    def test_unknown_segment_raises(self, tiny_graph):
        with pytest.raises(GraphOperationError, match="unknown segment"):
            tiny_graph.walk_sequence([("A", "+"), ("ghost", "+")])

    def test_overlap_between_is_symmetric_and_falls_back_to_the_default(self, tiny_graph):
        assert tiny_graph.overlap_between("A", "+", "B", "+") == 10
        assert tiny_graph.overlap_between("B", "-", "A", "-") == 10
        # No such link -> the graph default (0 here).
        assert tiny_graph.overlap_default == 0
        assert tiny_graph.overlap_between("A", "+", "D", "+") == 0
        tiny_graph.overlap_default = 55
        assert tiny_graph.overlap_between("A", "+", "D", "+") == 55
        # An existing link still reports its own overlap.
        assert tiny_graph.overlap_between("A", "+", "B", "+") == 10


class TestFindPaths:
    @staticmethod
    def _diamond() -> AssemblyGraph:
        """S -> {M1(300), M2(50)} -> E."""
        graph = AssemblyGraph()
        for name, length in (("S", 100), ("M1", 300), ("M2", 50), ("E", 100)):
            graph.add_segment(Segment(name, "A" * length))
        for a, b in (("S", "M1"), ("M1", "E"), ("S", "M2"), ("M2", "E")):
            graph.add_link(Link(a, "+", b, "+"))
        return graph

    def test_returns_the_intermediate_walk_only(self):
        graph = self._diamond()
        paths = graph.find_paths(("S", "+"), ("E", "+"), max_length=1000)
        assert sorted(paths) == [[("M1", "+")], [("M2", "+")]]

    def test_max_length_prunes_long_routes(self):
        graph = self._diamond()
        # M1 is 300 bp so only the 50 bp route fits.
        paths = graph.find_paths(("S", "+"), ("E", "+"), max_length=100)
        assert paths == [[("M2", "+")]]
        assert graph.find_paths(("S", "+"), ("E", "+"), max_length=10) == []

    def test_max_nodes_prunes_deep_routes(self):
        graph = self._diamond()
        assert graph.find_paths(("S", "+"), ("E", "+"), 1000, max_nodes=0) == []
        assert len(graph.find_paths(("S", "+"), ("E", "+"), 1000, max_nodes=1)) == 2

    def test_max_nodes_bounds_a_longer_chain(self):
        graph = AssemblyGraph()
        for name in ("S", "m1", "m2", "m3", "E"):
            graph.add_segment(Segment(name, "A" * 10))
        for a, b in (("S", "m1"), ("m1", "m2"), ("m2", "m3"), ("m3", "E")):
            graph.add_link(Link(a, "+", b, "+"))
        assert graph.find_paths(("S", "+"), ("E", "+"), 1000, max_nodes=3) == [
            [("m1", "+"), ("m2", "+"), ("m3", "+")]
        ]
        assert graph.find_paths(("S", "+"), ("E", "+"), 1000, max_nodes=2) == []

    def test_direct_adjacency_yields_one_empty_walk(self):
        graph = AssemblyGraph()
        for name in ("S", "E"):
            graph.add_segment(Segment(name, "A" * 10))
        graph.add_link(Link("S", "+", "E", "+"))
        assert graph.find_paths(("S", "+"), ("E", "+"), 1000) == [[]]

    def test_max_paths_caps_the_result(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("S", "A" * 10))
        graph.add_segment(Segment("E", "A" * 10))
        for i in range(6):
            name = f"m{i}"
            graph.add_segment(Segment(name, "A" * 10))
            graph.add_link(Link("S", "+", name, "+"))
            graph.add_link(Link(name, "+", "E", "+"))
        assert len(graph.find_paths(("S", "+"), ("E", "+"), 1000, max_paths=3)) == 3


# ---------------------------------------------------------------------------
# Mutation and export
# ---------------------------------------------------------------------------


class TestMutation:
    def test_remove_segment_drops_its_links_and_returns_an_inverse(self, tiny_graph):
        inverse = tiny_graph.remove_segment("B")
        assert "B" not in tiny_graph
        assert tiny_graph.link_count == 1  # only the R self-loop survives
        assert inverse["op"] == "add_segment"
        assert inverse["segment"].name == "B"
        assert len(inverse["links"]) == 2
        assert tiny_graph.successors("A", "+") == []

    def test_remove_unknown_segment_raises(self, tiny_graph):
        with pytest.raises(GraphOperationError, match="no such segment"):
            tiny_graph.remove_segment("ghost")

    def test_remove_link_accepts_either_spelling(self, tiny_graph):
        tiny_graph.remove_link(Link("B", "-", "A", "-"))
        assert tiny_graph.link_count == 2
        assert not tiny_graph.has_link("A", "+", "B", "+")
        with pytest.raises(GraphOperationError, match="no such link"):
            tiny_graph.remove_link(Link("A", "+", "B", "+"))

    def test_duplicate_segment_names_raise_unless_replacing(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("A", "AAAA"))
        with pytest.raises(GraphOperationError, match="duplicate segment"):
            graph.add_segment(Segment("A", "CCCC"))
        graph.add_segment(Segment("A", "CCCC"), replace=True)
        assert graph["A"].sequence == "CCCC"


class TestExport:
    def test_len_contains_getitem(self, tiny_graph):
        assert len(tiny_graph) == 5
        assert "A" in tiny_graph and "ghost" not in tiny_graph
        assert tiny_graph["A"].length == 100

    def test_total_length_and_link_count(self, tiny_graph):
        assert tiny_graph.total_length == TINY.total_length == 310
        assert tiny_graph.link_count == TINY.link_count == 3

    def test_segment_names_by_length(self, tiny_graph):
        assert tiny_graph.segment_names_by_length() == ["A", "B", "C", "D", "R"]

    def test_links_of(self, tiny_graph):
        assert len(tiny_graph.links_of("B")) == 2
        assert len(tiny_graph.links_of("D")) == 0
        assert len(tiny_graph.links_of("R")) == 1

    def test_summary(self, tiny_graph):
        assert tiny_graph.summary() == {
            "segments": 5,
            "links": 3,
            "paths": 0,
            "total_length": 310,
            "components": 3,
            "dead_ends": 4,
            "has_sequences": True,
        }

    def test_to_dict_carries_component_degree_and_circularity(self, tiny_graph):
        tiny_graph.add_path(Path("p1", [("A", "+"), ("B", "+")]))
        data = tiny_graph.to_dict()
        by_name = {s["name"]: s for s in data["segments"]}
        assert by_name["A"]["component"] == 0
        assert (by_name["A"]["deg_start"], by_name["A"]["deg_end"]) == (0, 1)
        assert by_name["R"]["circular"] is True
        assert by_name["A"]["circular"] is False
        assert len(data["links"]) == 3
        assert data["paths"] == [{"name": "p1", "steps": ["A+", "B+"], "gaps": []}]
        assert "sequence" not in by_name["A"]
        assert tiny_graph.to_dict(include_sequence=True)["segments"][0]["sequence"]


# ---------------------------------------------------------------------------
# Parser helpers
# ---------------------------------------------------------------------------


class TestDepthHelpers:
    @pytest.mark.parametrize("key", ["dp", "DP"])
    def test_direct_depth_tags(self, key):
        assert depth_from_tags({key: 12.5}, 1000) == pytest.approx(12.5)

    @pytest.mark.parametrize("key", ["KC", "kc", "RC", "rc", "FC", "fc"])
    def test_count_tags_are_divided_by_length(self, key):
        # KC, RC and FC are all counts in the GFA1 spec -- k-mer, read and
        # fragment. FC used to be read as a depth, inflating it by exactly the
        # segment length on every graph that carries it.
        assert depth_from_tags({key: 1000}, 100) == pytest.approx(10.0)

    def test_zero_depth_is_kept_not_treated_as_unknown(self):
        assert depth_from_tags({"DP": 0}, 100) == 0.0
        assert depth_from_tags({"KC": 0}, 100) == 0.0

    def test_dp_wins_over_kc(self):
        assert depth_from_tags({"dp": 7.0, "KC": 100_000}, 100) == pytest.approx(7.0)

    def test_no_usable_tag(self):
        assert depth_from_tags({}, 100) is None
        assert depth_from_tags({"KC": 1000}, 0) is None
        assert depth_from_tags({"dp": "not-a-number"}, 100) is None

    def test_depth_and_length_from_spades_names(self):
        name = "NODE_1_length_5000_cov_12.34"
        assert depth_from_name(name) == pytest.approx(12.34)
        assert length_from_name(name) == 5000
        assert depth_from_name("plain_name") is None
        assert length_from_name("plain_name") is None
