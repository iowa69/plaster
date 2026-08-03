"""Tests for plastr.core.analysis.align.

The hand-built alignments here are deliberately arithmetic-friendly so every
coordinate can be checked exactly.
"""

from __future__ import annotations

import pytest

from plastr.core.analysis.align import (
    Alignment,
    align_graph,
    alignment_backend,
    best_placements,
    merge_collinear,
    parse_cigar,
    reference_lengths,
    split_all,
    split_at_long_indels,
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
    query: str = "c1",
    ref: str = "chr1",
    strand: int = 1,
    q_len: int = 20_000,
    matches: int | None = None,
    block_len: int | None = None,
    is_primary: bool = True,
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
        matches=span if matches is None else matches,
        block_len=span if block_len is None else block_len,
        mapq=60,
        nm=nm,
        is_primary=is_primary,
    )


# ---------------------------------------------------------------------------
# Alignment properties
# ---------------------------------------------------------------------------


class TestAlignmentProperties:
    def test_spans_identity_and_coverage(self):
        aln = block(0, 5_000, 100, 5_100, matches=4_950, block_len=5_000)
        assert aln.q_span == 5_000
        assert aln.r_span == 5_000
        assert aln.identity == pytest.approx(0.99)
        assert aln.query_coverage == pytest.approx(0.25)  # 5000 of 20000

    def test_degenerate_values_do_not_divide_by_zero(self):
        aln = block(0, 0, 0, 0, q_len=0, matches=0, block_len=0)
        assert aln.identity == 0.0
        assert aln.query_coverage == 0.0

    def test_indel_counts_and_mismatch_estimate(self):
        aln = block(0, 1_000, 0, 1_010, nm=35)
        aln.cigar = [(500, "M"), (10, "D"), (400, "M"), (5, "I"), (100, "M")]
        assert aln.indel_counts() == (1, 1, 15)
        # 35 edits, 15 of which are indel bases -> 20 substitutions.
        assert aln.mismatch_estimate() == 20

    def test_mismatch_estimate_never_goes_negative(self):
        aln = block(0, 1_000, 0, 1_010, nm=5)
        aln.cigar = [(500, "M"), (100, "D"), (500, "M")]
        assert aln.mismatch_estimate() == 0

    def test_parse_cigar(self):
        assert parse_cigar("10M5D3I") == [(10, "M"), (5, "D"), (3, "I")]
        assert parse_cigar("") == []
        assert parse_cigar(None) == []

    def test_to_dict(self):
        aln = block(0, 100, 5, 105, matches=99, block_len=100)
        d = aln.to_dict()
        assert d["query"] == "c1" and d["ref"] == "chr1"
        assert d["q_st"] == 0 and d["r_en"] == 105
        assert d["identity"] == pytest.approx(0.99)


# ---------------------------------------------------------------------------
# split_at_long_indels
# ---------------------------------------------------------------------------


def _spanning_deletion(strand: int) -> Alignment:
    """One alignment that carries straight across a 5 kb reference deletion.

    5000M 5000D 5000M: 10 000 query bases against 15 000 reference bases,
    starting at reference position 1000. ``nm`` is 5050, i.e. the 5000 deleted
    bases plus 50 real substitutions.
    """
    aln = block(
        1_000 if strand < 0 else 0,
        11_000 if strand < 0 else 10_000,
        1_000,
        16_000,
        strand=strand,
        q_len=12_000 if strand < 0 else 10_000,
        matches=9_950,
        block_len=15_000,
        nm=5_050,
    )
    aln.cigar = [(5_000, "M"), (5_000, "D"), (5_000, "M")]
    return aln


class TestSplitAtLongIndels:
    def test_forward_strand_splits_into_two_pieces(self):
        pieces = split_at_long_indels(_spanning_deletion(1), min_indel=1_000)
        assert len(pieces) == 2

        left, right = pieces
        # Query coordinates tile the original range with no gap.
        assert (left.q_st, left.q_en) == (0, 5_000)
        assert (right.q_st, right.q_en) == (5_000, 10_000)
        # Reference coordinates skip exactly the 5 kb deletion.
        assert (left.r_st, left.r_en) == (1_000, 6_000)
        assert (right.r_st, right.r_en) == (11_000, 16_000)
        assert right.r_st - left.r_en == 5_000

    def test_the_deleted_bases_are_no_longer_counted_as_covered(self):
        pieces = split_at_long_indels(_spanning_deletion(1), min_indel=1_000)
        covered = sum(p.r_span for p in pieces)
        assert covered == 10_000  # not the original 15 000
        assert sum(p.block_len for p in pieces) == 10_000

    def test_nm_excludes_the_deletion_length(self):
        pieces = split_at_long_indels(_spanning_deletion(1), min_indel=1_000)
        # 5050 - 5000 deleted bases = 50 real edits, shared between the pieces.
        assert [p.nm for p in pieces] == [25, 25]
        assert sum(p.nm for p in pieces) == 50
        assert all(p.nm < 100 for p in pieces)

    def test_identity_stays_high(self):
        pieces = split_at_long_indels(_spanning_deletion(1), min_indel=1_000)
        assert [p.matches for p in pieces] == [4_975, 4_975]
        for piece in pieces:
            assert piece.identity == pytest.approx(0.995)
            assert piece.identity > 0.99

    def test_reverse_strand_query_coordinates_map_back_correctly(self):
        aln = _spanning_deletion(-1)
        pieces = split_at_long_indels(aln, min_indel=1_000)
        assert len(pieces) == 2

        first, second = pieces  # in reference order
        assert all(p.strand == -1 for p in pieces)
        # On the minus strand the CIGAR runs along the reference, so the first
        # reference piece corresponds to the *last* query bases.
        assert (first.r_st, first.r_en) == (1_000, 6_000)
        assert (first.q_st, first.q_en) == (6_000, 11_000)
        assert (second.r_st, second.r_en) == (11_000, 16_000)
        assert (second.q_st, second.q_en) == (1_000, 6_000)
        # Together they tile exactly the original query interval.
        assert min(p.q_st for p in pieces) == aln.q_st
        assert max(p.q_en for p in pieces) == aln.q_en
        assert sum(p.q_span for p in pieces) == aln.q_en - aln.q_st
        # Every coordinate stays inside the query.
        assert all(0 <= p.q_st < p.q_en <= aln.q_len for p in pieces)

    def test_reverse_strand_keeps_the_same_nm_and_matches_split(self):
        pieces = split_at_long_indels(_spanning_deletion(-1), min_indel=1_000)
        assert [p.nm for p in pieces] == [25, 25]
        assert [p.matches for p in pieces] == [4_975, 4_975]

    def test_a_threshold_above_the_indel_leaves_it_alone(self):
        aln = _spanning_deletion(1)
        assert split_at_long_indels(aln, min_indel=6_000) == [aln]

    def test_short_indels_are_never_split(self):
        aln = block(0, 1_000, 0, 1_100, matches=990, block_len=1_100, nm=110)
        aln.cigar = [(500, "M"), (100, "D"), (500, "M")]
        assert split_at_long_indels(aln, min_indel=1_000) == [aln]

    def test_an_alignment_without_a_cigar_is_untouched(self):
        aln = block(0, 1_000, 0, 1_000)
        assert split_at_long_indels(aln, min_indel=1_000) == [aln]

    def test_a_long_insertion_splits_the_query_not_the_reference(self):
        aln = block(0, 15_000, 0, 10_000, matches=9_950, block_len=15_000, nm=5_050)
        aln.cigar = [(5_000, "M"), (5_000, "I"), (5_000, "M")]
        pieces = split_at_long_indels(aln, min_indel=1_000)
        assert len(pieces) == 2
        assert (pieces[0].q_st, pieces[0].q_en) == (0, 5_000)
        assert (pieces[1].q_st, pieces[1].q_en) == (10_000, 15_000)
        assert (pieces[0].r_st, pieces[0].r_en) == (0, 5_000)
        assert (pieces[1].r_st, pieces[1].r_en) == (5_000, 10_000)

    def test_three_way_split(self):
        aln = block(0, 15_000, 0, 25_000, matches=14_900, block_len=25_000, nm=10_100)
        aln.cigar = [(5_000, "M"), (5_000, "D"), (5_000, "M"), (5_000, "D"), (5_000, "M")]
        pieces = split_at_long_indels(aln, min_indel=1_000)
        assert len(pieces) == 3
        assert [(p.r_st, p.r_en) for p in pieces] == [
            (0, 5_000),
            (10_000, 15_000),
            (20_000, 25_000),
        ]
        assert [(p.q_st, p.q_en) for p in pieces] == [
            (0, 5_000),
            (5_000, 10_000),
            (10_000, 15_000),
        ]
        # 10 000 of the 10 100 edits were the two deletions; only ~100 remain
        # (apportioned by block, hence the rounding slack).
        assert 95 <= sum(p.nm for p in pieces) <= 105

    def test_split_all_passes_untouched_alignments_through(self):
        plain = block(0, 1_000, 0, 1_000)
        out = split_all([plain, _spanning_deletion(1)], min_indel=1_000)
        assert len(out) == 3
        assert plain in out


# ---------------------------------------------------------------------------
# merge_collinear
# ---------------------------------------------------------------------------


class TestMergeCollinear:
    def test_adjacent_blocks_are_joined(self):
        a = block(0, 1_000, 0, 1_000)
        b = block(1_000, 2_000, 1_000, 2_000)
        merged = merge_collinear([a, b])
        assert len(merged) == 1
        joined = merged[0]
        assert (joined.q_st, joined.q_en) == (0, 2_000)
        assert (joined.r_st, joined.r_en) == (0, 2_000)
        assert joined.matches == 2_000
        assert joined.block_len == 2_000

    def test_a_small_gap_is_still_collinear(self):
        merged = merge_collinear([block(0, 1_000, 0, 1_000), block(1_050, 2_000, 1_060, 2_010)])
        assert len(merged) == 1
        assert merged[0].q_en == 2_000

    def test_a_big_reference_jump_is_refused(self):
        a = block(0, 1_000, 0, 1_000)
        b = block(1_000, 2_000, 50_000, 51_000)  # 49 kb jump
        merged = merge_collinear([a, b])
        assert len(merged) == 2
        assert sorted((m.r_st for m in merged)) == [0, 50_000]

    def test_a_big_query_jump_is_refused(self):
        merged = merge_collinear([block(0, 1_000, 0, 1_000), block(50_000, 51_000, 1_000, 2_000)])
        assert len(merged) == 2

    def test_a_large_backwards_step_is_refused(self):
        # r_gap of -5000 is well below the -500 tolerance.
        merged = merge_collinear([block(0, 1_000, 10_000, 11_000), block(1_000, 2_000, 5_000, 6_000)])
        assert len(merged) == 2

    def test_max_gap_is_configurable(self):
        blocks = [block(0, 1_000, 0, 1_000), block(1_000, 2_000, 20_000, 21_000)]
        assert len(merge_collinear(blocks, max_gap=10_000)) == 2
        assert len(merge_collinear(blocks, max_gap=30_000)) == 1

    def test_different_references_never_merge(self):
        blocks = [block(0, 1_000, 0, 1_000, ref="chr1"), block(1_000, 2_000, 1_000, 2_000, ref="chr2")]
        assert len(merge_collinear(blocks)) == 2

    def test_different_strands_never_merge(self):
        blocks = [block(0, 1_000, 0, 1_000, strand=1), block(1_000, 2_000, 1_000, 2_000, strand=-1)]
        assert len(merge_collinear(blocks)) == 2

    def test_different_queries_stay_separate(self):
        blocks = [block(0, 1_000, 0, 1_000, query="a"), block(1_000, 2_000, 1_000, 2_000, query="b")]
        assert len(merge_collinear(blocks)) == 2

    def test_reverse_strand_blocks_merge_in_descending_reference_order(self):
        # On the minus strand the reference runs backwards as the query advances.
        a = block(0, 1_000, 5_000, 6_000, strand=-1)
        b = block(1_000, 2_000, 4_000, 5_000, strand=-1)
        merged = merge_collinear([a, b])
        assert len(merged) == 1
        assert (merged[0].r_st, merged[0].r_en) == (4_000, 6_000)

    def test_inputs_are_not_mutated(self):
        a = block(0, 1_000, 0, 1_000)
        b = block(1_000, 2_000, 1_000, 2_000)
        merge_collinear([a, b])
        assert (a.q_en, a.r_en, a.matches) == (1_000, 1_000, 1_000)
        assert (b.q_st, b.r_st) == (1_000, 1_000)

    def test_empty_input(self):
        assert merge_collinear([]) == []


# ---------------------------------------------------------------------------
# best_placements
# ---------------------------------------------------------------------------


class TestBestPlacements:
    def test_picks_the_longest_primary_hit(self):
        short_primary = block(0, 3_000, 0, 3_000, ref="chrA")
        long_primary = block(5_000, 11_000, 0, 6_000, ref="chrB")
        longer_secondary = block(0, 9_000, 0, 9_000, ref="chrC", is_primary=False)

        best = best_placements([short_primary, longer_secondary, long_primary])
        assert set(best) == {"c1"}
        assert best["c1"].ref == "chrB"
        assert best["c1"].q_span == 6_000

    def test_a_secondary_hit_wins_only_when_no_primary_exists(self):
        best = best_placements([block(0, 9_000, 0, 9_000, ref="chrC", is_primary=False)])
        assert best["c1"].ref == "chrC"

    def test_blocks_are_merged_before_the_longest_is_chosen(self):
        # One contig broken into two collinear blocks on chrA (8 kb together)
        # against a single 6 kb block on chrB.
        hits = [
            block(0, 4_000, 0, 4_000, ref="chrA"),
            block(4_000, 8_000, 4_000, 8_000, ref="chrA"),
            block(8_000, 14_000, 0, 6_000, ref="chrB"),
        ]
        best = best_placements(hits)
        assert best["c1"].ref == "chrA"
        assert best["c1"].q_span == 8_000

    def test_min_query_coverage_filter(self):
        hits = [block(0, 6_000, 0, 6_000, q_len=20_000)]  # coverage 0.30
        assert best_placements(hits, min_query_coverage=0.2)["c1"].q_span == 6_000
        assert best_placements(hits, min_query_coverage=0.5) == {}

    def test_min_identity_filter(self):
        poor = block(0, 6_000, 0, 6_000, matches=3_000, block_len=6_000)  # 50%
        assert poor.identity == pytest.approx(0.5)
        assert best_placements([poor], min_identity=0.9) == {}
        assert set(best_placements([poor], min_identity=0.4)) == {"c1"}

    def test_one_entry_per_query(self):
        hits = [
            block(0, 5_000, 0, 5_000, query="a"),
            block(0, 5_000, 0, 5_000, query="b"),
            block(6_000, 9_000, 20_000, 23_000, query="b"),
        ]
        assert set(best_placements(hits)) == {"a", "b"}

    def test_empty_input(self):
        assert best_placements([]) == {}


# ---------------------------------------------------------------------------
# Real alignments against the demo reference
# ---------------------------------------------------------------------------


class TestAgainstTheDemoReference:
    def test_reference_lengths(self, demo_paths, demo_reference_lengths):
        assert demo_reference_lengths == {"chromosome": 160_000, "plasmid": 9_000}
        assert sum(demo_reference_lengths.values()) == 169_000

    def test_every_contig_but_the_foreign_one_aligns(self, demo_alignments):
        aligned = {a.query for a in demo_alignments}
        assert "ctg_foreign" not in aligned
        assert "ctg_clean_1" in aligned
        assert len(aligned) == 9

    def test_a_clean_contig_lands_where_it_was_cut_from(self, demo_alignments):
        hits = [a for a in demo_alignments if a.query == "ctg_clean_1" and a.is_primary]
        assert len(hits) == 1
        hit = hits[0]
        assert hit.ref == "chromosome"
        assert hit.strand == 1
        assert hit.r_st == 0
        assert hit.r_en == 19_000
        # The generator mutates at 0.1%, so identity should be ~0.999.
        assert hit.identity > 0.99

    def test_the_inversion_contig_produces_both_strands(self, demo_alignments):
        strands = {a.strand for a in demo_alignments if a.query == "ctg_inversion"}
        assert strands == {1, -1}

    def test_the_translocation_contig_hits_both_references(self, demo_alignments):
        refs = {a.ref for a in demo_alignments if a.query == "ctg_translocation"}
        assert refs == {"chromosome", "plasmid"}

    def test_min_length_skips_short_queries(self, demo_graph, demo_paths):
        hits = align_graph(demo_graph, str(demo_paths["reference"]), min_length=20_000)
        # Only ctg_clean_2 (32 kb) and ctg_clean_4 (25 kb) are that long.
        assert {a.query for a in hits} == {"ctg_clean_2", "ctg_clean_4"}

    def test_reference_lengths_helper_matches_the_fasta(self, demo_paths):
        assert reference_lengths(str(demo_paths["reference"])) == {
            "chromosome": 160_000,
            "plasmid": 9_000,
        }
