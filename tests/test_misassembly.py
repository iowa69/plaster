"""Reference-based evaluation: the ground-truth test.

``examples/make_demo_data.py`` plants exactly one inversion, one relocation, one
translocation and one unalignable contig, so the counts below are not
approximations -- they are what the generator built.
"""

from __future__ import annotations

import pytest

from assemblage.core.analysis.align import Alignment, alignment_backend
from assemblage.core.analysis.misassembly import (
    INVERSION,
    LOCAL,
    RELOCATION,
    TRANSLOCATION,
    ReferenceReport,
    classify_misassemblies,
    evaluate_against_reference,
)

pytestmark = pytest.mark.skipif(
    alignment_backend() == "none",
    reason="neither the mappy module nor the minimap2 binary is available",
)


def block(
    q_st: int,
    q_en: int,
    r_st: int,
    r_en: int,
    *,
    query: str = "ctg",
    ref: str = "chr1",
    strand: int = 1,
    q_len: int = 20_000,
    nm: int = 0,
) -> Alignment:
    span = q_en - q_st
    return Alignment(
        query=query,
        q_len=q_len,
        q_st=q_st,
        q_en=q_en,
        strand=strand,
        ref=ref,
        r_len=1_000_000,
        r_st=r_st,
        r_en=r_en,
        matches=span,
        block_len=span,
        mapq=60,
        nm=nm,
    )


# ---------------------------------------------------------------------------
# The ground-truth test
# ---------------------------------------------------------------------------


class TestDemoGroundTruth:
    def test_the_generator_recorded_what_it_planted(self, demo_truth):
        assert demo_truth == {
            "inversions": 1,
            "relocations": 1,
            "translocations": 1,
            "unaligned_contigs": 1,
            "reference_length": 169_000,
        }

    def test_misassembly_counts_are_exactly_the_planted_ones(self, demo_report, demo_truth):
        assert demo_report.num_inversions == demo_truth["inversions"] == 1
        assert demo_report.num_relocations == demo_truth["relocations"] == 1
        assert demo_report.num_translocations == demo_truth["translocations"] == 1
        assert demo_report.num_misassemblies == 3
        assert demo_report.num_local_misassemblies == 0

    def test_the_right_contigs_are_blamed(self, demo_report):
        assert demo_report.misassembled_contigs == [
            "ctg_inversion",
            "ctg_relocation",
            "ctg_translocation",
        ]
        kinds = {m.contig: m.kind for m in demo_report.misassemblies}
        assert kinds == {
            "ctg_inversion": INVERSION,
            "ctg_relocation": RELOCATION,
            "ctg_translocation": TRANSLOCATION,
        }

    def test_each_break_is_near_the_place_it_was_planted(self, demo_report):
        by_contig = {m.contig: m for m in demo_report.misassemblies}
        # The inversion contig is 8 kb forward + 8 kb reverse complemented.
        assert by_contig["ctg_inversion"].contig_pos == pytest.approx(8_000, abs=200)
        # The relocation contig is 6 kb + 7 kb, skipping 14 kb of chromosome.
        assert by_contig["ctg_relocation"].contig_pos == pytest.approx(6_000, abs=200)
        assert by_contig["ctg_relocation"].ref_gap == pytest.approx(14_000, abs=200)
        # The translocation contig is 6 kb of chromosome + 5 kb of plasmid.
        assert by_contig["ctg_translocation"].contig_pos == pytest.approx(6_000, abs=200)
        assert by_contig["ctg_translocation"].left_ref == "chromosome"
        assert by_contig["ctg_translocation"].right_ref == "plasmid"

    def test_all_reported_misassemblies_are_extensive(self, demo_report):
        assert all(m.is_extensive for m in demo_report.misassemblies)

    def test_exactly_one_contig_is_unalignable(self, demo_report, demo_truth):
        assert demo_report.unaligned_contigs == demo_truth["unaligned_contigs"] == 1
        assert demo_report.aligned_contigs == 9
        assert demo_report.fully_unaligned_length == 3_500  # ctg_foreign

    def test_genome_fraction(self, demo_report):
        # The contigs deliberately leave parts of the chromosome uncovered.
        assert 80.0 < demo_report.genome_fraction < 84.0
        assert demo_report.reference_length == 169_000
        assert demo_report.reference_sequences == 2
        assert demo_report.covered_bases == pytest.approx(
            demo_report.genome_fraction * 169_000 / 100, abs=1
        )

    def test_duplication_ratio_is_near_one(self, demo_report):
        assert demo_report.duplication_ratio == pytest.approx(1.0, abs=0.02)

    def test_mismatch_rate_reflects_the_simulated_snp_rate(self, demo_report):
        # The generator mutates at 0.1% == 100 per 100 kb.
        assert 60.0 < demo_report.mismatches_per_100kb < 140.0
        assert demo_report.total_mismatches > 0

    def test_no_indels_were_simulated(self, demo_report):
        assert demo_report.total_indels == 0
        assert demo_report.indels_per_100kb == 0.0

    def test_aligned_block_statistics(self, demo_report):
        assert demo_report.largest_alignment == 32_000  # ctg_clean_2
        assert demo_report.na50 > 0
        assert demo_report.nga50 is not None
        assert demo_report.nga50 <= demo_report.na50

    def test_per_reference_coverage(self, demo_report):
        assert set(demo_report.per_reference_coverage) == {"chromosome", "plasmid"}
        for ref, pct in demo_report.per_reference_coverage.items():
            assert 0.0 < pct <= 100.0, ref
        assert demo_report.per_reference_coverage["plasmid"] == pytest.approx(100.0, abs=1.0)

    def test_coverage_blocks_are_merged_and_sorted(self, demo_report):
        for ref, blocks in demo_report.coverage_blocks.items():
            assert blocks == sorted(blocks), ref
            for (_a, b), (c, _d) in zip(blocks, blocks[1:]):
                assert c > b, f"{ref} blocks overlap or touch"

    def test_report_serialises(self, demo_report):
        d = demo_report.to_dict()
        assert d["num_misassemblies"] == 3
        assert isinstance(d["misassemblies"], list)
        assert d["misassemblies"][0]["kind"] in {INVERSION, RELOCATION, TRANSLOCATION}

    def test_an_explicit_genome_size_changes_only_the_normalised_numbers(
        self, _demo_alignments_session, demo_reference_lengths, demo_graph
    ):
        contig_lengths = {n: s.length for n, s in demo_graph.segments.items()}
        doubled = evaluate_against_reference(
            _demo_alignments_session,
            demo_reference_lengths,
            contig_lengths,
            genome_size=338_000,
        )
        assert doubled.covered_bases > 0
        assert doubled.genome_fraction == pytest.approx(
            100.0 * doubled.covered_bases / 338_000
        )
        assert doubled.num_misassemblies == 3


class TestEvaluateEdgeCases:
    def test_no_alignments_at_all(self):
        report = evaluate_against_reference([], {"chr1": 1000}, {"c1": 500, "c2": 400})
        assert isinstance(report, ReferenceReport)
        assert report.genome_fraction == 0.0
        assert report.covered_bases == 0
        assert report.duplication_ratio == 0.0
        assert report.aligned_contigs == 0
        assert report.unaligned_contigs == 2
        assert report.fully_unaligned_length == 900
        assert report.num_misassemblies == 0
        assert report.mismatches_per_100kb == 0.0

    def test_a_perfectly_covered_reference(self):
        alns = [block(0, 1_000, 0, 1_000, query="c1")]
        report = evaluate_against_reference(alns, {"chr1": 1_000}, {"c1": 1_000})
        assert report.genome_fraction == pytest.approx(100.0)
        assert report.duplication_ratio == pytest.approx(1.0)
        assert report.unaligned_contigs == 0
        assert report.partially_unaligned == 0

    def test_duplication_shows_up_as_a_ratio_above_one(self):
        # Two contigs claiming the same 1 kb of reference.
        alns = [
            block(0, 1_000, 0, 1_000, query="c1"),
            block(0, 1_000, 0, 1_000, query="c2"),
        ]
        report = evaluate_against_reference(alns, {"chr1": 1_000}, {"c1": 1_000, "c2": 1_000})
        assert report.covered_bases == 1_000
        assert report.total_aligned_length == 2_000
        assert report.duplication_ratio == pytest.approx(2.0)

    def test_partially_unaligned_contigs_are_counted(self):
        alns = [block(0, 500, 0, 500, query="c1", q_len=5_000)]
        report = evaluate_against_reference(alns, {"chr1": 5_000}, {"c1": 5_000})
        assert report.aligned_contigs == 1
        assert report.partially_unaligned == 1

    def test_blocks_below_min_block_are_ignored(self):
        alns = [block(0, 150, 0, 150, query="c1")]
        report = evaluate_against_reference(alns, {"chr1": 1_000}, {"c1": 150})
        assert report.aligned_contigs == 0
        assert report.covered_bases == 0


# ---------------------------------------------------------------------------
# classify_misassemblies, one category at a time
# ---------------------------------------------------------------------------


class TestClassifyMisassemblies:
    def test_translocation(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000, ref="chrA"),
                block(5_000, 10_000, 0, 5_000, ref="chrB"),
            ]
        )
        assert len(events) == 1
        event = events[0]
        assert event.kind == TRANSLOCATION
        assert event.is_extensive is True
        assert event.contig_pos == 5_000
        assert (event.left_ref, event.right_ref) == ("chrA", "chrB")
        assert event.ref_gap is None
        assert "different reference" in event.description

    def test_inversion(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000, strand=1),
                block(5_000, 10_000, 5_000, 10_000, strand=-1),
            ]
        )
        assert len(events) == 1
        event = events[0]
        assert event.kind == INVERSION
        assert event.is_extensive is True
        assert event.contig_pos == 5_000
        assert "strand flips" in event.description

    def test_relocation(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000),
                block(5_000, 10_000, 20_000, 25_000),  # 15 kb forward jump
            ]
        )
        assert len(events) == 1
        event = events[0]
        assert event.kind == RELOCATION
        assert event.is_extensive is True
        assert event.ref_gap == 15_000
        assert event.contig_pos == 5_000
        assert "15,000 bp" in event.description
        assert "forward" in event.description

    def test_backward_relocation(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 20_000, 25_000),
                block(5_000, 10_000, 0, 5_000),  # jumps 20 kb backwards
            ]
        )
        assert len(events) == 1
        assert events[0].kind == RELOCATION
        assert events[0].ref_gap == -25_000  # 0 - 25 000
        assert "backward" in events[0].description

    def test_local_misassembly_from_a_backward_shift(self):
        # 600 bp is above LOCAL_THRESHOLD (200) and below EXTENSIVE (1000).
        events = classify_misassemblies(
            [
                block(0, 5_000, 5_000, 10_000),
                block(5_000, 10_000, 9_400, 14_400),
            ]
        )
        assert len(events) == 1
        event = events[0]
        assert event.kind == LOCAL
        assert event.is_extensive is False
        assert event.ref_gap == -600
        assert "local shift of 600 bp" in event.description

    def test_local_misassembly_from_a_500bp_forward_shift(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000),
                block(5_000, 10_000, 5_500, 10_500),  # 500 bp forward shift
            ]
        )
        assert [e.kind for e in events] == [LOCAL]
        assert events[0].ref_gap == 500

    def test_a_five_kb_forward_relocation_is_extensive(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000),
                block(5_000, 10_000, 10_000, 15_000),  # 5 kb forward jump
            ]
        )
        assert [e.kind for e in events] == [RELOCATION]

    def test_collinear_blocks_are_not_a_misassembly(self):
        events = classify_misassemblies(
            [block(0, 5_000, 0, 5_000), block(5_000, 10_000, 5_000, 10_000)]
        )
        assert events == []

    def test_a_single_block_can_never_be_a_misassembly(self):
        assert classify_misassemblies([block(0, 5_000, 0, 5_000)]) == []

    def test_blocks_shorter_than_min_block_are_dropped(self):
        events = classify_misassemblies(
            [block(0, 5_000, 0, 5_000, ref="chrA"), block(5_000, 5_100, 0, 100, ref="chrB")]
        )
        assert events == []

    def test_a_repeat_hitting_many_places_is_not_a_pile_of_misassemblies(self):
        # One 5 kb contig aligning to four different places on the reference:
        # the blocks overlap on the contig, so only the best one is kept.
        hits = [
            block(0, 5_000, offset, offset + 5_000)
            for offset in (0, 50_000, 100_000, 150_000)
        ]
        assert classify_misassemblies(hits) == []

    def test_two_contigs_are_classified_independently(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000, query="a", ref="chrA"),
                block(5_000, 10_000, 0, 5_000, query="a", ref="chrB"),
                block(0, 5_000, 0, 5_000, query="b", strand=1),
                block(5_000, 10_000, 5_000, 10_000, query="b", strand=-1),
            ]
        )
        assert [(e.contig, e.kind) for e in events] == [
            ("a", TRANSLOCATION),
            ("b", INVERSION),
        ]

    def test_three_blocks_give_two_breakpoints(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000, ref="chrA"),
                block(5_000, 10_000, 0, 5_000, ref="chrB"),
                block(10_000, 15_000, 0, 5_000, ref="chrC"),
            ]
        )
        assert [e.kind for e in events] == [TRANSLOCATION, TRANSLOCATION]
        assert [e.contig_pos for e in events] == [5_000, 10_000]

    def test_results_are_sorted_by_contig_then_position(self):
        events = classify_misassemblies(
            [
                block(0, 5_000, 0, 5_000, query="z", ref="chrA"),
                block(5_000, 10_000, 0, 5_000, query="z", ref="chrB"),
                block(0, 5_000, 0, 5_000, query="a", ref="chrA"),
                block(5_000, 10_000, 0, 5_000, query="a", ref="chrB"),
            ]
        )
        assert [e.contig for e in events] == ["a", "z"]
