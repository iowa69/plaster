"""Tests for plastr.core.sequence."""

from __future__ import annotations

import pytest

from plastr.core.sequence import flip, gc_content, n_count, oriented, revcomp


class TestRevcomp:
    def test_plain_acgt(self):
        assert revcomp("ACGT") == "ACGT"
        assert revcomp("AAAA") == "TTTT"
        assert revcomp("ACCTGA") == "TCAGGT"

    def test_empty(self):
        assert revcomp("") == ""

    def test_is_an_involution(self):
        seq = "ACGTTGCAnnRYSWKM-ACGT"
        assert revcomp(revcomp(seq)) == seq

    def test_iupac_ambiguity_codes(self):
        # R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N, then reversed.
        assert revcomp("RYSWKMBDHVN") == "NBDHVKMWSRY"

    def test_lowercase_is_preserved(self):
        assert revcomp("aAcCgGtT") == "AaCcGgTt"
        assert revcomp("acgt") == "acgt"

    def test_uracil_becomes_adenine(self):
        assert revcomp("ACGU") == "ACGT"
        assert revcomp("acgu") == "acgt"

    def test_gap_characters_survive(self):
        assert revcomp("A-C") == "G-T"

    def test_length_is_never_changed(self):
        seq = "ACGTNRYKMacgtn-U"
        assert len(revcomp(seq)) == len(seq)


class TestGcContent:
    def test_empty_sequence_is_zero_not_an_error(self):
        assert gc_content("") == 0.0

    def test_all_n_is_zero_not_an_error(self):
        assert gc_content("NNNNNNNN") == 0.0
        assert gc_content("nnnn") == 0.0

    @pytest.mark.parametrize(
        "seq,expected",
        [
            ("GC", 1.0),
            ("AT", 0.0),
            ("GCAT", 0.5),
            ("gcat", 0.5),
            ("GgCcAaTt", 0.5),
            ("G", 1.0),
            ("A", 0.0),
        ],
    )
    def test_known_fractions(self, seq, expected):
        assert gc_content(seq) == pytest.approx(expected)

    def test_ns_are_excluded_from_the_denominator(self):
        # 2 of the 4 unambiguous bases are G/C; the Ns must not dilute it.
        assert gc_content("GCATNNNN") == pytest.approx(0.5)

    def test_s_code_counts_as_gc(self):
        # 'S' means "strong" (G or C), so it counts in both numerator and
        # denominator; this is the documented behaviour of the _GC set.
        assert gc_content("S") == pytest.approx(1.0)
        assert gc_content("SSAT") == pytest.approx(0.5)


class TestNCount:
    def test_counts_both_cases(self):
        assert n_count("ACGT") == 0
        assert n_count("NNACGTnn") == 4
        assert n_count("") == 0


class TestOrientedAndFlip:
    def test_oriented_forward_returns_the_input(self):
        assert oriented("ACCTGA", "+") == "ACCTGA"

    def test_oriented_reverse_is_revcomp(self):
        assert oriented("ACCTGA", "-") == revcomp("ACCTGA")
        assert oriented("ACCTGA", "-") == "TCAGGT"

    def test_oriented_twice_returns_the_original(self):
        seq = "ACGTTGCA"
        assert oriented(oriented(seq, "-"), "-") == seq

    def test_flip(self):
        assert flip("+") == "-"
        assert flip("-") == "+"
        assert flip(flip("+")) == "+"
