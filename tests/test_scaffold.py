"""Scaffolding: plan, graph bridging, and the FASTA/AGP contract.

The central test here is :meth:`TestFastaAgpContract.test_every_agp_row_describes
_its_slice_of_the_fasta`. AGP is only useful if downstream tools can trust it, so
every row is checked against the sequence it claims to describe.
"""

from __future__ import annotations

from collections import defaultdict, deque

import pytest

from plastr.core.analysis.align import align_graph, alignment_backend
from plastr.core.analysis.misassembly import evaluate_against_reference
from plastr.core.analysis.operations import break_misassemblies
from plastr.core.io.agp import COMPONENT_CONTIG, GAP_KNOWN, GAP_UNKNOWN, read_agp
from plastr.core.io.fasta import read_fasta
from plastr.core.model import AssemblyGraph, Link, Segment
from plastr.core.scaffold.builder import build_scaffolds, write_scaffolds
from plastr.core.scaffold.graph_bridge import (
    _bridge_sequence,
    bridge_gap,
    bridge_pieces,
    find_unbranching_paths,
)
from plastr.core.scaffold.plan import (
    GAP_ADJACENT,
    GAP_DEFAULT,
    GAP_GRAPH,
    GAP_REFERENCE,
    Scaffold,
    ScaffoldMember,
    ScaffoldPlan,
)
from plastr.core.scaffold.reference_guided import scaffold_by_reference
from plastr.core.sequence import revcomp

pytestmark = pytest.mark.skipif(
    alignment_backend() == "none",
    reason="neither the mappy module nor the minimap2 binary is available",
)

GAP_TYPES = (GAP_KNOWN, GAP_UNKNOWN)


@pytest.fixture
def demo_plan(demo_graph, demo_alignments) -> ScaffoldPlan:
    return scaffold_by_reference(demo_graph, demo_alignments)


@pytest.fixture
def demo_built(demo_graph, demo_plan):
    return build_scaffolds(demo_graph, demo_plan)


def expected_bridges(graph: AssemblyGraph, plan: ScaffoldPlan):
    """component_id -> the bridging sequences the builder should have emitted.

    Recomputed straight from the graph so the check does not simply trust the
    string the builder happened to write.
    """
    out: dict[str, deque] = defaultdict(deque)
    for scaffold in plan.scaffolds:
        members = [m for m in scaffold.members if m.segment in graph.segments]
        for index, member in enumerate(members[:-1]):
            if not (member.bridge_path or member.gap_after == 0):
                continue
            nxt = members[index + 1]
            seq = _bridge_sequence(
                graph,
                (member.segment, member.orientation),
                list(member.bridge_path),
                (nxt.segment, nxt.orientation),
            )
            if not seq:
                continue
            key = "+".join(f"{n}{o}" for n, o in member.bridge_path) or (
                f"{member.segment}_bridge"
            )
            out[key].append(seq)
    return out


# ---------------------------------------------------------------------------
# The FASTA / AGP contract
# ---------------------------------------------------------------------------


class TestFastaAgpContract:
    def test_the_demo_data_builds_without_warnings(self, demo_built):
        assert demo_built.warnings == []
        assert len(demo_built.records) == 2
        assert {name for name, _ in demo_built.records} == {
            "scaffold_chromosome",
            "scaffold_plasmid",
        }

    def test_every_agp_row_describes_its_slice_of_the_fasta(self, demo_graph, demo_built):
        """The contract that must never break.

        Every component row names a segment that exists and slices out exactly
        the sub-range it declares -- including the rows covering sequence
        recovered from the graph, which are written per segment of the walk
        rather than as one synthetic component.
        """
        sequences = dict(demo_built.records)
        contig_rows = gap_rows = 0

        for row in demo_built.agp_rows:
            assert row.object_name in sequences, row.object_name
            piece = sequences[row.object_name][row.object_beg - 1 : row.object_end]
            # The slice must be inside the sequence, not silently truncated.
            assert len(piece) == row.object_end - row.object_beg + 1

            if row.component_type == COMPONENT_CONTIG:
                assert len(piece) == row.component_end - row.component_beg + 1
                assert row.component_id in demo_graph.segments, (
                    f"AGP names component {row.component_id!r}, which is not a "
                    "segment -- no consumer could resolve it"
                )
                segment = demo_graph[row.component_id]
                assert 1 <= row.component_beg <= row.component_end <= segment.length
                # AGP component coordinates are in the forward contig, whatever
                # the orientation column says; the slice is reverse-complemented
                # afterwards. Slicing the already-oriented sequence instead
                # agrees only when nothing was trimmed.
                forward = segment.sequence[row.component_beg - 1 : row.component_end]
                expected = revcomp(forward) if row.orientation == "-" else forward
                assert piece == expected, (
                    f"{row.object_name} part {row.part_number} does not match "
                    f"{row.component_id}{row.orientation}"
                )
                assert "N" not in piece
                contig_rows += 1
            else:
                assert row.component_type in GAP_TYPES
                assert set(piece) == {"N"}
                assert len(piece) == row.gap_length
                gap_rows += 1

        assert contig_rows == 10  # 9 placed contigs + 1 graph bridge
        assert gap_rows == demo_built.num_gaps == 4

    def test_the_rows_tile_each_scaffold_with_no_gap_or_overlap(self, demo_built):
        by_object = defaultdict(list)
        for row in demo_built.agp_rows:
            by_object[row.object_name].append(row)
        for name, sequence in demo_built.records:
            rows = sorted(by_object[name], key=lambda r: r.part_number)
            assert [r.part_number for r in rows] == list(range(1, len(rows) + 1))
            assert rows[0].object_beg == 1
            assert rows[-1].object_end == len(sequence)
            rebuilt = "".join(sequence[r.object_beg - 1 : r.object_end] for r in rows)
            assert rebuilt == sequence

    def test_the_graph_bridge_carries_the_real_repeat_sequence(self, demo_built, demo_graph):
        """The bridged gap is written as the repeat segment it walks through.

        It used to be one row whose component id was the synthesised string
        'ctg_repeat+' -- a name no FASTA contains, with an orientation of '+'
        whatever direction the walk took.
        """
        sequences = dict(demo_built.records)
        # The bridge sits between the first two members of the chromosome
        # scaffold, i.e. AGP part 2, and covers the repeat in full.
        row = next(
            r
            for r in demo_built.agp_rows
            if r.object_name == "scaffold_chromosome" and r.part_number == 2
        )
        assert row.component_id == "ctg_repeat"
        assert row.orientation == "+"
        assert (row.component_beg, row.component_end) == (1, 2_400)
        piece = sequences[row.object_name][row.object_beg - 1 : row.object_end]
        assert piece == demo_graph["ctg_repeat"].sequence
        assert len(piece) == 2_400

    def test_totals_add_up(self, demo_graph, demo_plan, demo_built):
        member_bases = sum(
            demo_graph[m.segment].length
            for s in demo_plan.scaffolds
            for m in s.members
            if m.segment in demo_graph.segments
        )
        assert demo_built.bridged_bases == 2_400
        assert demo_built.gap_bases == sum(
            r.gap_length for r in demo_built.agp_rows if r.component_type in GAP_TYPES
        )
        assert (
            demo_built.total_length
            == member_bases + demo_built.gap_bases + demo_built.bridged_bases
        )
        assert demo_built.total_length == sum(demo_built.lengths)
        assert demo_built.total_length == sum(len(s) for _n, s in demo_built.records)

    def test_written_files_agree_with_the_in_memory_build(self, demo_graph, demo_built, tmp_path):
        fasta_path = tmp_path / "scaffolds.fasta"
        agp_path = tmp_path / "scaffolds.agp"
        write_scaffolds(demo_built, str(fasta_path), str(agp_path), comments=["test run"])

        from_disk = {name: seq for name, _desc, seq in read_fasta(str(fasta_path))}
        assert from_disk == dict(demo_built.records)

        rows = list(read_agp(str(agp_path)))
        assert rows == demo_built.agp_rows

        # And the contract holds against the files, not just the objects.
        for row in rows:
            piece = from_disk[row.object_name][row.object_beg - 1 : row.object_end]
            if row.component_type == COMPONENT_CONTIG and row.component_id in demo_graph.segments:
                assert piece == demo_graph[row.component_id].seq_oriented(row.orientation)
            elif row.component_type in GAP_TYPES:
                assert piece == "N" * row.gap_length

    def test_reverse_oriented_members_are_revcomped_in_the_fasta(self):
        from plastr.core.sequence import revcomp

        seq_a = "AAAACCCCGGGGTTTA"
        seq_b = "ACGTACGTACGTACGA"
        # Both must be asymmetric or the test would pass without revcomping.
        assert revcomp(seq_a) != seq_a and revcomp(seq_b) != seq_b

        graph = AssemblyGraph()
        graph.add_segment(Segment("a", seq_a))
        graph.add_segment(Segment("b", seq_b))
        plan = ScaffoldPlan(
            scaffolds=[
                Scaffold(
                    name="s1",
                    members=[
                        ScaffoldMember("a", "-", gap_after=10, gap_evidence=GAP_REFERENCE),
                        ScaffoldMember("b", "-"),
                    ],
                )
            ]
        )
        built = build_scaffolds(graph, plan, min_gap=10)
        assert built.warnings == []
        name, sequence = built.records[0]
        assert name == "s1"
        assert sequence == revcomp(seq_a) + "N" * 10 + revcomp(seq_b)
        rows = built.agp_rows
        assert [r.component_type for r in rows] == [COMPONENT_CONTIG, GAP_KNOWN, COMPONENT_CONTIG]
        assert rows[0].orientation == "-"
        assert sequence[rows[0].object_beg - 1 : rows[0].object_end] == graph["a"].seq_oriented("-")
        assert sequence[rows[2].object_beg - 1 : rows[2].object_end] == graph["b"].seq_oriented("-")

    def test_min_gap_is_enforced(self):
        graph = AssemblyGraph()
        for name in ("a", "b"):
            graph.add_segment(Segment(name, "ACGT" * 10))
        plan = ScaffoldPlan(
            scaffolds=[
                Scaffold(
                    name="s1",
                    members=[
                        ScaffoldMember("a", "+", gap_after=5, gap_evidence=GAP_DEFAULT),
                        ScaffoldMember("b", "+"),
                    ],
                )
            ]
        )
        built = build_scaffolds(graph, plan, min_gap=100)
        assert built.gap_bases == 100
        assert built.agp_rows[1].gap_length == 100
        # gap_evidence 'default' is an unknown-size gap.
        assert built.agp_rows[1].component_type == GAP_UNKNOWN
        assert built.records[0][1].count("N") == 100

    def test_a_missing_segment_is_warned_about_not_silently_dropped(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("a", "ACGT" * 10))
        plan = ScaffoldPlan(
            scaffolds=[
                Scaffold(
                    name="s1",
                    members=[ScaffoldMember("a", "+"), ScaffoldMember("ghost", "+")],
                )
            ]
        )
        built = build_scaffolds(graph, plan)
        assert any("ghost" in w for w in built.warnings)
        assert built.records[0][1] == "ACGT" * 10

    def test_a_segment_without_sequence_is_a_hard_error(self):
        from plastr.core.errors import GraphOperationError

        graph = AssemblyGraph()
        graph.add_segment(Segment("a", None, 100))
        plan = ScaffoldPlan(scaffolds=[Scaffold(name="s1", members=[ScaffoldMember("a", "+")])])
        with pytest.raises(GraphOperationError, match="no sequence"):
            build_scaffolds(graph, plan)


# ---------------------------------------------------------------------------
# The plan itself
# ---------------------------------------------------------------------------


class TestScaffoldPlan:
    def test_demo_plan_shape(self, demo_plan):
        assert demo_plan.method == "reference"
        assert len(demo_plan.scaffolds) == 2
        assert demo_plan.placed_count == 9
        assert demo_plan.unplaced == ["ctg_foreign"]
        names = [s.name for s in demo_plan.scaffolds]
        assert names == ["scaffold_chromosome", "scaffold_plasmid"]
        assert demo_plan.scaffolds[0].reference == "chromosome"

    def test_members_are_ordered_along_the_reference(self, demo_plan):
        for scaffold in demo_plan.scaffolds:
            starts = [m.ref_start for m in scaffold.members]
            assert starts == sorted(starts), scaffold.name

    def test_member_index(self, demo_plan):
        index = demo_plan.member_index()
        assert index["ctg_clean_1"] == ("scaffold_chromosome", 0)
        assert index["ctg_plasmid"] == ("scaffold_plasmid", 0)

    def test_gap_evidence_values_are_from_the_documented_set(self, demo_plan):
        allowed = {GAP_REFERENCE, GAP_GRAPH, GAP_ADJACENT, GAP_DEFAULT, "manual"}
        for scaffold in demo_plan.scaffolds:
            for member in scaffold.members:
                assert member.gap_evidence in allowed

    def test_the_repeat_gap_was_filled_from_the_graph(self, demo_plan):
        first = demo_plan.scaffolds[0].members[0]
        assert first.segment == "ctg_clean_1"
        assert first.gap_evidence == GAP_GRAPH
        assert first.bridge_path == [("ctg_repeat", "+")]
        assert first.bridge_sequence_length == 2_400
        assert first.gap_after == 0
        assert any("recovered from the assembly graph" in n for n in demo_plan.notes)

    def test_to_dict_from_dict_round_trips_losslessly(self, demo_plan):
        first = demo_plan.to_dict()
        second = ScaffoldPlan.from_dict(first).to_dict()
        assert second == first

    def test_round_trip_preserves_bridge_paths_and_orientations(self):
        plan = ScaffoldPlan(
            method="manual",
            unplaced=["u1"],
            redundant=["r1"],
            notes=["a note"],
            scaffolds=[
                Scaffold(
                    name="s1",
                    source="manual",
                    reference="chr1",
                    members=[
                        ScaffoldMember(
                            segment="c1",
                            orientation="-",
                            gap_after=250,
                            gap_evidence=GAP_GRAPH,
                            bridge_path=[("m1", "+"), ("m2", "-")],
                            bridge_sequence_length=500,
                            ref="chr1",
                            ref_start=10,
                            ref_end=1010,
                            identity=0.987,
                            overlaps_previous=True,
                        ),
                        ScaffoldMember(segment="c2"),
                    ],
                )
            ],
        )
        restored = ScaffoldPlan.from_dict(plan.to_dict())
        assert restored.to_dict() == plan.to_dict()
        member = restored.scaffolds[0].members[0]
        assert member.bridge_path == [("m1", "+"), ("m2", "-")]
        assert member.orientation == "-"
        assert member.identity == pytest.approx(0.987)
        assert member.overlaps_previous is True
        assert restored.unplaced == ["u1"] and restored.redundant == ["r1"]

    def test_from_dict_accepts_bridge_paths_as_pairs(self):
        data = {
            "method": "manual",
            "scaffolds": [
                {"name": "s", "members": [{"segment": "c", "bridge_path": [["m", "+"]]}]}
            ],
        }
        plan = ScaffoldPlan.from_dict(data)
        assert plan.scaffolds[0].members[0].bridge_path == [("m", "+")]

    def test_an_empty_plan_round_trips(self):
        plan = ScaffoldPlan()
        assert ScaffoldPlan.from_dict(plan.to_dict()).to_dict() == plan.to_dict()

    def test_include_unplaced_adds_singleton_scaffolds(self, demo_graph, demo_alignments):
        plan = scaffold_by_reference(demo_graph, demo_alignments, include_unplaced=True)
        names = [s.name for s in plan.scaffolds]
        assert "scaffold_unplaced_ctg_foreign" in names
        built = build_scaffolds(demo_graph, plan)
        assert built.warnings == []
        # Every contig now appears in the output.
        assert sum(built.lengths) >= demo_graph.total_length


# ---------------------------------------------------------------------------
# graph_bridge
# ---------------------------------------------------------------------------


class TestBridgeGap:
    def test_directly_linked_contigs_close_to_nothing(self, demo_graph):
        # ctg_clean_4 -> ctg_inversion is a real link in the demo graph.
        assert demo_graph.has_link("ctg_clean_4", "+", "ctg_inversion", "+")
        assert bridge_gap(demo_graph, ("ctg_clean_4", "+"), ("ctg_inversion", "+"), 5_000) == (
            [],
            "",
        )

    def test_no_path_returns_none(self, demo_graph):
        assert bridge_gap(demo_graph, ("ctg_clean_1", "+"), ("ctg_foreign", "+"), 1_000) is None

    def test_unknown_segments_return_none(self, demo_graph):
        assert bridge_gap(demo_graph, ("ghost", "+"), ("ctg_repeat", "+"), 100) is None
        assert bridge_gap(demo_graph, ("ctg_repeat", "+"), ("ghost", "+"), 100) is None

    def test_a_graph_without_links_returns_none(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("a", "A" * 100))
        graph.add_segment(Segment("b", "C" * 100))
        assert bridge_gap(graph, ("a", "+"), ("b", "+"), 100) is None

    def test_the_repeat_path_is_recovered(self, demo_graph):
        result = bridge_gap(demo_graph, ("ctg_clean_1", "+"), ("ctg_clean_2", "+"), 2_400)
        assert result is not None
        path, sequence = result
        assert path == [("ctg_repeat", "+")]
        assert sequence == demo_graph["ctg_repeat"].sequence
        assert len(sequence) == 2_400

    def test_a_path_of_the_wrong_length_is_rejected(self, demo_graph):
        # The only route is 2400 bp; asking for a 20 kb gap must not accept it.
        assert bridge_gap(demo_graph, ("ctg_clean_1", "+"), ("ctg_clean_2", "+"), 20_000) is None

    def test_ambiguous_gaps_are_refused(self):
        """Two different sequences of the same length must not be guessed at."""
        graph = AssemblyGraph()
        graph.add_segment(Segment("left", "A" * 500))
        graph.add_segment(Segment("right", "T" * 500))
        graph.add_segment(Segment("m1", "C" * 300))
        graph.add_segment(Segment("m2", "G" * 300))
        for middle in ("m1", "m2"):
            graph.add_link(Link("left", "+", middle, "+"))
            graph.add_link(Link(middle, "+", "right", "+"))
        assert bridge_gap(graph, ("left", "+"), ("right", "+"), 300) is None

    def test_overlaps_are_trimmed_out_of_the_bridge(self):
        graph = AssemblyGraph()
        graph.add_segment(Segment("left", "A" * 100))
        graph.add_segment(Segment("mid", "C" * 100))
        graph.add_segment(Segment("right", "G" * 100))
        graph.add_link(Link("left", "+", "mid", "+", 10, "10M"))
        graph.add_link(Link("mid", "+", "right", "+", 10, "10M"))
        result = bridge_gap(graph, ("left", "+"), ("right", "+"), 90)
        assert result is not None
        path, sequence = result
        assert path == [("mid", "+")]
        # 100 bp of 'mid' minus the 10 bp it shares with 'left'.
        assert sequence == "C" * 90


class TestFindUnbranchingPaths:
    def test_a_linear_chain_is_one_maximal_path(self, chain_graph):
        chains = find_unbranching_paths(chain_graph)
        assert len(chains) == 1
        assert chains[0] == [("P", "+"), ("A", "+"), ("B", "+"), ("C", "+"), ("Q", "+")]

    def test_a_branch_stops_the_chain(self):
        graph = AssemblyGraph()
        for name in ("a", "b", "c", "d"):
            graph.add_segment(Segment(name, "ACGT" * 10))
        graph.add_link(Link("a", "+", "b", "+"))
        graph.add_link(Link("b", "+", "c", "+"))
        graph.add_link(Link("b", "+", "d", "+"))
        chains = find_unbranching_paths(graph)
        assert chains == [[("a", "+"), ("b", "+")]]

    def test_min_length_filters_short_chains(self, chain_graph):
        assert find_unbranching_paths(chain_graph, min_length=10) == []


# ---------------------------------------------------------------------------
# Breaking misassemblies feeds back into a clean evaluation
# ---------------------------------------------------------------------------


class TestBreakMisassembliesEndToEnd:
    def test_breaking_the_demo_contigs_removes_every_misassembly(
        self, demo_graph, demo_alignments, demo_reference_lengths, demo_paths
    ):
        contig_lengths = {n: s.length for n, s in demo_graph.segments.items()}
        before = evaluate_against_reference(
            demo_alignments, demo_reference_lengths, contig_lengths
        )
        assert before.num_misassemblies == 3

        record = break_misassemblies(demo_graph, before.misassemblies)
        assert record["count"] == 3
        assert record["pieces"] == 6
        assert sorted(record["added_segments"]) == [
            "ctg_inversion_part1",
            "ctg_inversion_part2",
            "ctg_relocation_part1",
            "ctg_relocation_part2",
            "ctg_translocation_part1",
            "ctg_translocation_part2",
        ]
        for name in ("ctg_inversion", "ctg_relocation", "ctg_translocation"):
            assert name not in demo_graph

        after_alignments = align_graph(demo_graph, str(demo_paths["reference"]))
        after = evaluate_against_reference(
            after_alignments,
            demo_reference_lengths,
            {n: s.length for n, s in demo_graph.segments.items()},
        )
        assert after.num_misassemblies == 0
        assert after.num_inversions == 0
        assert after.num_relocations == 0
        assert after.num_translocations == 0
        assert after.misassembled_contigs == []
        # Breaking must not lose sequence or reference coverage.
        assert demo_graph.total_length == 137_400
        assert after.genome_fraction == pytest.approx(before.genome_fraction, abs=0.5)

    def test_scaffolding_the_broken_graph_still_honours_the_agp_contract(
        self, demo_graph, demo_alignments, demo_reference_lengths, demo_paths
    ):
        contig_lengths = {n: s.length for n, s in demo_graph.segments.items()}
        report = evaluate_against_reference(
            demo_alignments, demo_reference_lengths, contig_lengths
        )
        break_misassemblies(demo_graph, report.misassemblies)

        alignments = align_graph(demo_graph, str(demo_paths["reference"]))
        plan = scaffold_by_reference(demo_graph, alignments)
        built = build_scaffolds(demo_graph, plan)

        assert built.warnings == []
        sequences = dict(built.records)
        for row in built.agp_rows:
            piece = sequences[row.object_name][row.object_beg - 1 : row.object_end]
            if row.component_type in GAP_TYPES:
                assert piece == "N" * row.gap_length
            elif row.component_id in demo_graph.segments:
                assert piece == demo_graph[row.component_id].seq_oriented(row.orientation)


class TestStaleBridges:
    """A hand-edited plan must never splice in a walk the graph cannot support.

    Reordering or flipping a member leaves the previously-computed ``bridge_path``
    describing a walk that no longer exists. Splicing its sequence in anyway would
    invent a join the reads never supported, so the bridge is dropped and an
    honest gap is written instead.
    """

    @staticmethod
    def _plan():
        from plastr.core.scaffold.plan import Scaffold, ScaffoldMember, ScaffoldPlan

        plan = ScaffoldPlan(method="manual")
        plan.scaffolds.append(
            Scaffold(
                name="s1",
                members=[
                    ScaffoldMember("A", "+", gap_after=0, bridge_path=[("M", "+")]),
                    ScaffoldMember("B", "+"),
                ],
            )
        )
        return plan

    @staticmethod
    def _graph():
        from plastr.core.model import AssemblyGraph, Link, Segment

        g = AssemblyGraph("bridge")
        g.add_segment(Segment("A", "A" * 500))
        g.add_segment(Segment("M", "C" * 300))
        g.add_segment(Segment("B", "G" * 500))
        g.add_link(Link("A", "+", "M", "+"))
        g.add_link(Link("M", "+", "B", "+"))
        return g

    def test_a_valid_bridge_is_still_spliced_in(self):
        from plastr.core.scaffold.builder import build_scaffolds

        graph = self._graph()
        built = build_scaffolds(graph, self._plan())
        assert built.warnings == []
        assert built.bridged_bases == 300
        assert built.records[0][1] == "A" * 500 + "C" * 300 + "G" * 500

    def test_flipping_a_member_drops_the_bridge_and_leaves_a_gap(self):
        from plastr.core.scaffold.builder import build_scaffolds

        graph = self._graph()
        plan = self._plan()
        plan.scaffolds[0].members[0].orientation = "-"  # no link from A- to M+

        built = build_scaffolds(graph, plan)
        assert any("does not connect" in w for w in built.warnings)
        assert built.bridged_bases == 0
        assert built.gap_bases >= 100
        sequence = built.records[0][1]
        assert "C" * 300 not in sequence, "spliced in a walk the graph does not support"
        assert "N" * 100 in sequence

    def test_the_agp_still_describes_the_fasta_after_a_dropped_bridge(self):
        from plastr.core.io.agp import GAP_KNOWN, GAP_UNKNOWN
        from plastr.core.scaffold.builder import build_scaffolds

        graph = self._graph()
        plan = self._plan()
        plan.scaffolds[0].members[0].orientation = "-"
        built = build_scaffolds(graph, plan)

        sequence = dict(built.records)["s1"]
        cursor = 1
        for row in sorted(built.agp_rows, key=lambda r: r.part_number):
            assert row.object_beg == cursor
            piece = sequence[row.object_beg - 1 : row.object_end]
            if row.component_type in (GAP_KNOWN, GAP_UNKNOWN):
                assert set(piece) == {"N"}
                assert len(piece) == row.gap_length
            else:
                assert piece == graph[row.component_id].seq_oriented(row.orientation)
            cursor = row.object_end + 1
        assert cursor - 1 == len(sequence)


class TestScaffoldReconstructsTheGenome:
    """The test that matters: does the exported sequence equal the truth?

    Checking that the AGP agrees with the FASTA is necessary but not sufficient
    -- both are generated from one coordinate counter, so they agree even when
    the sequence is wrong. This suite once passed while every graph-closed gap
    duplicated the link's overlap, making the scaffold k-1 bases too long per
    join on exactly the SPAdes/Unicycler graphs the tool targets.
    """

    @staticmethod
    def _cut_genome(overlap, seed=3, length=3000, cuts=(0, 800, 1500, 2300)):
        """A known genome cut into pieces that share `overlap` bases, as a k-mer
        assembler emits them."""
        import random

        from plastr.core.model import AssemblyGraph, Link, Segment

        rng = random.Random(seed)
        genome = "".join(rng.choice("ACGT") for _ in range(length))
        bounds = [*cuts, length]
        g = AssemblyGraph("cut")
        for i in range(len(bounds) - 1):
            start = bounds[i]
            end = bounds[i + 1] + (overlap if i < len(bounds) - 2 else 0)
            g.add_segment(Segment(f"s{i}", genome[start:end]))
            if i:
                g.add_link(
                    Link(f"s{i - 1}", "+", f"s{i}", "+", overlap,
                         f"{overlap}M" if overlap else "*")
                )
        return genome, g

    def _scaffold(self, graph):
        from plastr.core.project import Project
        from plastr.core.scaffold.builder import build_scaffolds

        p = Project()
        p.graph = graph
        p.build_plan(method="graph")
        built = build_scaffolds(graph, p.plan)
        return built, max((s for _n, s in built.records), key=len)

    @pytest.mark.parametrize("overlap", [0, 21, 55, 127])
    def test_the_scaffold_equals_the_genome_it_came_from(self, overlap):
        genome, graph = self._cut_genome(overlap)
        built, scaffold = self._scaffold(graph)
        assert built.warnings == []
        assert len(scaffold) == len(genome), (
            f"scaffold is {len(scaffold) - len(genome):+d} bp off the truth; "
            f"3 joins x {overlap} bp overlap"
        )
        assert scaffold == genome

    def test_a_bridged_gap_and_a_reversed_member_still_reconstruct_the_genome(self):
        import random

        from plastr.core.io.agp import GAP_KNOWN, GAP_UNKNOWN
        from plastr.core.model import AssemblyGraph, Link, Segment
        from plastr.core.scaffold.builder import build_scaffolds
        from plastr.core.scaffold.plan import Scaffold, ScaffoldMember, ScaffoldPlan
        from plastr.core.sequence import revcomp

        rng = random.Random(9)
        overlap = 40
        genome = "".join(rng.choice("ACGT") for _ in range(2000))
        left, middle, right = (
            genome[0 : 800 + overlap],
            genome[800 : 1400 + overlap],
            genome[1400:2000],
        )
        g = AssemblyGraph("bridge")
        g.add_segment(Segment("L", left))
        g.add_segment(Segment("M", middle))
        g.add_segment(Segment("Rrc", revcomp(right)))  # stored on the other strand
        g.add_link(Link("L", "+", "M", "+", overlap, f"{overlap}M"))
        g.add_link(Link("M", "+", "Rrc", "-", overlap, f"{overlap}M"))

        plan = ScaffoldPlan(method="manual")
        plan.scaffolds.append(
            Scaffold(
                name="s1",
                members=[
                    ScaffoldMember("L", "+", gap_after=0, bridge_path=[("M", "+")]),
                    ScaffoldMember("Rrc", "-"),
                ],
            )
        )
        built = build_scaffolds(g, plan)
        scaffold = dict(built.records)["s1"]
        assert built.warnings == []
        assert scaffold == genome

        # The AGP must describe the trimmed component, in original contig
        # coordinates -- a reversed member is trimmed at its *end*.
        rows = sorted(built.agp_rows, key=lambda r: r.part_number)
        cursor = 1
        for row in rows:
            assert row.object_beg == cursor
            piece = scaffold[row.object_beg - 1 : row.object_end]
            if row.component_type in (GAP_KNOWN, GAP_UNKNOWN):
                assert set(piece) == {"N"}
            elif row.component_id in g.segments:
                source = g.segments[row.component_id].sequence
                span = source[row.component_beg - 1 : row.component_end]
                assert piece == (revcomp(span) if row.orientation == "-" else span)
            cursor = row.object_end + 1
        assert cursor - 1 == len(scaffold)

        reversed_row = next(r for r in rows if r.component_id == "Rrc")
        assert reversed_row.component_beg == 1
        assert reversed_row.component_end == len(right) - overlap


def test_a_contig_placed_twice_is_reported():
    """An editing slip that places one contig in two scaffolds writes it twice.

    A repeat legitimately recurs, so this is a warning rather than an error --
    but it must not pass silently.
    """
    from plastr.core.model import AssemblyGraph, Segment
    from plastr.core.scaffold.builder import build_scaffolds
    from plastr.core.scaffold.plan import Scaffold, ScaffoldMember, ScaffoldPlan

    g = AssemblyGraph("dup")
    g.add_segment(Segment("a", "A" * 500))
    g.add_segment(Segment("b", "C" * 500))
    plan = ScaffoldPlan()
    plan.scaffolds.append(
        Scaffold(
            name="s1",
            members=[ScaffoldMember("a", gap_after=100), ScaffoldMember("b")],
        )
    )
    plan.scaffolds.append(Scaffold(name="s2", members=[ScaffoldMember("a")]))

    built = build_scaffolds(g, plan)
    assert any("placed 2 times" in w and "a" in w for w in built.warnings)
    # and it really is written twice, which is what the warning is about
    assert sum(seq.count("A" * 500) for _n, seq in built.records) == 2


def test_a_single_placement_produces_no_warning():
    from plastr.core.model import AssemblyGraph, Segment
    from plastr.core.scaffold.builder import build_scaffolds
    from plastr.core.scaffold.plan import Scaffold, ScaffoldMember, ScaffoldPlan

    g = AssemblyGraph("ok")
    g.add_segment(Segment("a", "A" * 500))
    g.add_segment(Segment("b", "C" * 500))
    plan = ScaffoldPlan()
    plan.scaffolds.append(
        Scaffold(
            name="s1",
            members=[ScaffoldMember("a", gap_after=100), ScaffoldMember("b")],
        )
    )
    assert build_scaffolds(g, plan).warnings == []


# ---------------------------------------------------------------------------
# Regressions from reviewing real assemblies
# ---------------------------------------------------------------------------


class TestBridgeAgpNamesRealSegments:
    """A graph-bridged join used to be one row with a synthesised component id.

    ``"+".join(f"{n}{o}" ...)`` produced 'seg0146-' for a one-step walk and
    'a+b-' for a two-step one: names present in no FASTA, coordinates that were
    not coordinates in them, and orientation '+' whichever way the walk ran.
    """

    def test_pieces_reassemble_the_bridge_exactly(self, demo_graph, demo_plan):
        for scaffold in demo_plan.scaffolds:
            members = [m for m in scaffold.members if m.segment in demo_graph.segments]
            for index, member in enumerate(members[:-1]):
                if not member.bridge_path:
                    continue
                nxt = members[index + 1]
                bridge = _bridge_sequence(
                    demo_graph,
                    (member.segment, member.orientation),
                    list(member.bridge_path),
                    (nxt.segment, nxt.orientation),
                )
                pieces = bridge_pieces(
                    demo_graph,
                    (member.segment, member.orientation),
                    list(member.bridge_path),
                )
                assert pieces is not None
                rebuilt = ""
                for name, orient, beg, end in pieces:
                    forward = demo_graph[name].sequence[beg - 1 : end]
                    rebuilt += revcomp(forward) if orient == "-" else forward
                assert rebuilt == bridge

    def test_no_component_id_carries_an_orientation_suffix(self, demo_built, demo_graph):
        for row in demo_built.agp_rows:
            if row.component_type != COMPONENT_CONTIG:
                continue
            assert not row.component_id.endswith(("+", "-"))
            assert row.component_id in demo_graph.segments


class TestPlacedCountIsHonest:
    def test_unplaced_singletons_do_not_count_as_placed(self):
        plan = ScaffoldPlan(method="reference")
        plan.scaffolds.append(
            Scaffold(
                name="s1",
                source="reference",
                members=[
                    ScaffoldMember(segment="a", orientation="+"),
                    ScaffoldMember(segment="b", orientation="+"),
                ],
            )
        )
        plan.scaffolds.append(
            Scaffold(
                name="s_unplaced_c",
                source="unplaced",
                members=[ScaffoldMember(segment="c", orientation="+")],
            )
        )
        plan.unplaced = ["c"]

        assert plan.placed_count == 2
        assert plan.scaffold_count == 1
        # placed + unplaced must partition the assembly, never exceed it
        assert plan.placed_count + len(plan.unplaced) == 3
        data = plan.to_dict()
        assert data["placed_count"] == 2
        assert data["scaffold_count"] == 1
        assert data["record_count"] == 2


class TestScaffoldNamesAreUnique:
    def test_references_that_sanitise_alike_get_distinct_names(self, demo_graph):
        from plastr.core.scaffold.reference_guided import _unique

        used: set[str] = set()
        assert _unique("scaffold_plasmid_1", used) == "scaffold_plasmid_1"
        assert _unique("scaffold_plasmid_1", used) == "scaffold_plasmid_1_2"
        assert _unique("scaffold_plasmid_1", used) == "scaffold_plasmid_1_3"
