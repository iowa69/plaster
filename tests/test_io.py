"""Tests for the readers and writers: GFA, FASTA, FASTG, AGP."""

from __future__ import annotations

import gzip
import textwrap

import pytest

from assemblage.core.errors import AssemblageFormatError
from assemblage.core.io.agp import (
    COMPONENT_CONTIG,
    GAP_KNOWN,
    AgpRow,
    read_agp,
    write_agp,
)
from assemblage.core.io.fasta import (
    fasta_string,
    read_fasta,
    read_fasta_dict,
    write_fasta,
)
from assemblage.core.io.fastg import read_fastg
from assemblage.core.io.gfa import (
    cigar_overlap,
    detect_format,
    read_gfa,
    write_gfa,
)
from assemblage.core.io.loader import load_graph, read_contigs_fasta
from assemblage.core.sequence import revcomp


def write(path, text: str):
    path.write_text(textwrap.dedent(text).lstrip("\n"))
    return path


def _fingerprint(graph):
    """Everything a GFA round trip must preserve."""
    segments = {
        name: (seg.length, seg.sequence, seg.depth)
        for name, seg in graph.segments.items()
    }
    links = {
        key: (lk.from_name, lk.from_orient, lk.to_name, lk.to_orient, lk.overlap)
        for key, lk in graph.links.items()
    }
    return segments, links


# ---------------------------------------------------------------------------
# GFA
# ---------------------------------------------------------------------------


class TestGfaRoundTrip:
    def test_demo_graph_survives_read_write_read(self, demo_paths, tmp_path):
        first = read_gfa(str(demo_paths["gfa"]))
        out = tmp_path / "rt.gfa"
        write_gfa(first, str(out))
        second = read_gfa(str(out))

        assert _fingerprint(first) == _fingerprint(second)
        assert len(second.segments) == 10
        assert second.link_count == 9
        # Sequences really are there, not just equal-and-empty.
        assert all(s.sequence for s in second.segments.values())
        assert second.total_length == first.total_length == 137_400

    def test_round_trip_of_a_handwritten_file(self, tmp_path):
        src = write(
            tmp_path / "in.gfa",
            """
            H\tVN:Z:1.0
            S\ts1\tACGTACGTAA\tdp:f:12.5
            S\ts2\tTTTTGGGGCC\tdp:f:8
            S\ts3\t*\tLN:i:4321
            L\ts1\t+\ts2\t-\t4M
            L\ts2\t+\ts3\t+\t*
            P\tp1\ts1+,s2-\t*
            """,
        )
        first = read_gfa(str(src))
        out = tmp_path / "out.gfa"
        write_gfa(first, str(out))
        second = read_gfa(str(out))

        assert _fingerprint(first) == _fingerprint(second)
        assert second["s1"].sequence == "ACGTACGTAA"
        assert second["s3"].sequence is None and second["s3"].length == 4321
        assert second.has_link("s1", "+", "s2", "-")
        assert second.overlap_between("s1", "+", "s2", "-") == 4
        assert list(second.paths) == ["p1"]
        assert second.paths["p1"].steps == [("s1", "+"), ("s2", "-")]

    def test_reverse_written_link_reads_back_as_the_same_connection(self, tmp_path):
        src = write(
            tmp_path / "in.gfa",
            """
            H\tVN:Z:1.0
            S\ta\tACGTACGTAC
            S\tb\tTTTTGGGGCC
            L\ta\t+\tb\t+\t0M
            L\tb\t-\ta\t-\t0M
            """,
        )
        graph = read_gfa(str(src))
        assert graph.link_count == 1


class TestGfaDepthParsing:
    @pytest.mark.parametrize(
        "tag,expected",
        [
            ("dp:f:31.5", 31.5),
            ("DP:f:31.5", 31.5),
            ("KC:i:3150", 31.5),  # 3150 / 100 bp
            ("RC:i:3150", 31.5),
        ],
    )
    def test_depth_from_tags(self, tmp_path, tag, expected):
        src = write(
            tmp_path / "d.gfa",
            f"H\tVN:Z:1.0\nS\tseg\t{'ACGT' * 25}\t{tag}\n",
        )
        graph = read_gfa(str(src))
        assert graph["seg"].length == 100
        assert graph["seg"].depth == pytest.approx(expected)

    def test_depth_from_a_spades_style_name(self, tmp_path):
        src = write(
            tmp_path / "d.gfa",
            f"H\tVN:Z:1.0\nS\tNODE_7_length_100_cov_42.75\t{'ACGT' * 25}\n",
        )
        graph = read_gfa(str(src))
        assert graph["NODE_7_length_100_cov_42.75"].depth == pytest.approx(42.75)

    def test_tag_beats_name(self, tmp_path):
        src = write(
            tmp_path / "d.gfa",
            f"H\tVN:Z:1.0\nS\tNODE_7_length_100_cov_42.75\t{'ACGT' * 25}\tdp:f:3.5\n",
        )
        assert read_gfa(str(src))["NODE_7_length_100_cov_42.75"].depth == pytest.approx(3.5)

    def test_no_depth_information_leaves_it_none(self, tmp_path):
        src = write(tmp_path / "d.gfa", f"H\tVN:Z:1.0\nS\tplain\t{'ACGT' * 25}\n")
        assert read_gfa(str(src))["plain"].depth is None


class TestGfaPlaceholderSequences:
    def test_star_sequence_with_ln_tag(self, tmp_path):
        src = write(
            tmp_path / "s.gfa",
            """
            H\tVN:Z:1.0
            S\tbig\t*\tLN:i:250000\tdp:f:5.5
            S\treal\tACGTACGT
            """,
        )
        graph = read_gfa(str(src))
        assert graph["big"].sequence is None
        assert graph["big"].has_sequence is False
        assert graph["big"].length == 250_000
        assert graph["big"].depth == pytest.approx(5.5)
        assert graph.total_length == 250_008

    def test_star_sequence_with_kc_tag_uses_the_ln_length(self, tmp_path):
        src = write(
            tmp_path / "s.gfa",
            "H\tVN:Z:1.0\nS\tbig\t*\tLN:i:1000\tKC:i:25000\n",
        )
        assert read_gfa(str(src))["big"].depth == pytest.approx(25.0)

    def test_star_sequence_falls_back_to_the_length_in_the_name(self, tmp_path):
        src = write(
            tmp_path / "s.gfa",
            "H\tVN:Z:1.0\nS\tNODE_2_length_777_cov_3\t*\n",
        )
        graph = read_gfa(str(src))
        assert graph["NODE_2_length_777_cov_3"].length == 777

    def test_star_sequence_is_written_back_with_ln(self, tmp_path):
        src = write(tmp_path / "s.gfa", "H\tVN:Z:1.0\nS\tbig\t*\tLN:i:4242\n")
        out = tmp_path / "o.gfa"
        write_gfa(read_gfa(str(src)), str(out))
        assert "LN:i:4242" in out.read_text()
        assert read_gfa(str(out))["big"].length == 4242


class TestGfaErrors:
    def test_bad_link_orientation(self, tmp_path):
        src = write(
            tmp_path / "bad.gfa",
            """
            H\tVN:Z:1.0
            S\ta\tACGTACGT
            S\tb\tACGTACGT
            L\ta\tx\tb\t+\t0M
            """,
        )
        with pytest.raises(AssemblageFormatError) as exc:
            read_gfa(str(src))
        assert exc.value.line_no == 4
        assert "orientation" in exc.value.message
        assert str(src) in str(exc.value)

    def test_s_line_too_short(self, tmp_path):
        src = write(
            tmp_path / "bad.gfa",
            """
            H\tVN:Z:1.0
            S\tlonely
            """,
        )
        with pytest.raises(AssemblageFormatError) as exc:
            read_gfa(str(src))
        assert exc.value.line_no == 2
        assert "S line" in exc.value.message

    def test_l_line_too_short(self, tmp_path):
        src = write(
            tmp_path / "bad.gfa",
            """
            H\tVN:Z:1.0
            S\ta\tACGT
            L\ta\t+\tb
            """,
        )
        with pytest.raises(AssemblageFormatError) as exc:
            read_gfa(str(src))
        assert exc.value.line_no == 3
        assert "L line" in exc.value.message

    def test_no_s_lines(self, tmp_path):
        src = write(
            tmp_path / "bad.gfa",
            """
            H\tVN:Z:1.0
            L\ta\t+\tb\t+\t0M
            """,
        )
        with pytest.raises(AssemblageFormatError) as exc:
            read_gfa(str(src))
        assert "no S (segment) lines" in exc.value.message
        # There is no single offending line for this one, but the file is named.
        assert exc.value.path == str(src)
        assert exc.value.line_no is None

    def test_unrecognised_record_only_fails_in_strict_mode(self, tmp_path):
        src = write(
            tmp_path / "odd.gfa",
            """
            H\tVN:Z:1.0
            S\ta\tACGTACGT
            Z\tsomething\telse
            """,
        )
        assert len(read_gfa(str(src)).segments) == 1
        with pytest.raises(AssemblageFormatError) as exc:
            read_gfa(str(src), strict=True)
        assert exc.value.line_no == 3


class TestCigarOverlap:
    @pytest.mark.parametrize(
        "cigar,expected",
        [("*", 0), ("", 0), ("55M", 55), ("10M5D10M", 25), ("10M5I10M", 20), ("3=2X", 5)],
    )
    def test_overlap_length(self, cigar, expected):
        assert cigar_overlap(cigar) == expected

    def test_modal_overlap_becomes_the_graph_default(self, tmp_path):
        src = write(
            tmp_path / "ov.gfa",
            """
            H\tVN:Z:1.0
            S\ta\tACGTACGTAC
            S\tb\tACGTACGTAC
            S\tc\tACGTACGTAC
            S\td\tACGTACGTAC
            L\ta\t+\tb\t+\t55M
            L\tb\t+\tc\t+\t55M
            L\tc\t+\td\t+\t3M
            """,
        )
        graph = read_gfa(str(src))
        assert graph.overlap_default == 55
        assert graph.overlap_between("c", "+", "d", "+") == 3


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------


class TestFastaReading:
    def test_multiline_records(self, tmp_path):
        src = write(
            tmp_path / "m.fa",
            """
            >first this is the description
            ACGTACGTAC
            GGGGTTTTAA
            CC
            >second
            TTTT
            """,
        )
        records = list(read_fasta(str(src)))
        assert records == [
            ("first", "this is the description", "ACGTACGTACGGGGTTTTAACC"),
            ("second", "", "TTTT"),
        ]

    def test_blank_lines_and_semicolon_comments_are_ignored(self, tmp_path):
        src = write(
            tmp_path / "c.fa",
            """
            >a

            ; an ancient FASTA comment
            ACGT

            ACGT
            """,
        )
        assert list(read_fasta(str(src))) == [("a", "", "ACGTACGT")]

    def test_gzipped_input(self, tmp_path):
        src = tmp_path / "g.fa.gz"
        with gzip.open(src, "wt") as handle:
            handle.write(">a desc\nACGTACGT\nTTTT\n>b\nGGGG\n")
        assert list(read_fasta(str(src))) == [
            ("a", "desc", "ACGTACGTTTTT"),
            ("b", "", "GGGG"),
        ]

    def test_empty_header_raises(self, tmp_path):
        src = write(tmp_path / "e.fa", ">\nACGT\n")
        with pytest.raises(AssemblageFormatError) as exc:
            list(read_fasta(str(src)))
        assert "empty header" in exc.value.message
        assert exc.value.line_no == 1

    def test_sequence_before_any_header_raises(self, tmp_path):
        src = write(tmp_path / "e.fa", "ACGT\n>a\nACGT\n")
        with pytest.raises(AssemblageFormatError) as exc:
            list(read_fasta(str(src)))
        assert "before any" in exc.value.message
        assert exc.value.line_no == 1

    def test_a_file_without_records_raises(self, tmp_path):
        src = write(tmp_path / "e.fa", "\n\n")
        with pytest.raises(AssemblageFormatError, match="no FASTA records"):
            list(read_fasta(str(src)))

    def test_read_fasta_dict_rejects_duplicates(self, tmp_path):
        src = write(tmp_path / "d.fa", ">a\nACGT\n>a\nTTTT\n")
        with pytest.raises(AssemblageFormatError, match="duplicate sequence name"):
            read_fasta_dict(str(src))

    def test_read_fasta_dict(self, tmp_path):
        src = write(tmp_path / "d.fa", ">a\nACGT\n>b\nTTTT\n")
        assert read_fasta_dict(str(src)) == {"a": "ACGT", "b": "TTTT"}


class TestFastaWriting:
    def test_default_wrap_is_60(self):
        text = fasta_string([("s", "A" * 145)])
        lines = text.splitlines()
        assert lines[0] == ">s"
        assert [len(x) for x in lines[1:]] == [60, 60, 25]
        assert "".join(lines[1:]) == "A" * 145

    @pytest.mark.parametrize("wrap", [1, 7, 10, 80])
    def test_arbitrary_wrap_widths(self, wrap):
        seq = "ACGT" * 30  # 120 bp
        lines = fasta_string([("s", seq)], wrap=wrap).splitlines()[1:]
        assert all(len(line) <= wrap for line in lines)
        assert all(len(line) == wrap for line in lines[:-1])
        assert "".join(lines) == seq

    def test_wrap_zero_writes_one_line(self):
        seq = "ACGT" * 30
        lines = fasta_string([("s", seq)], wrap=0).splitlines()
        assert lines == [">s", seq]

    def test_round_trip_through_a_file(self, tmp_path):
        records = [("a", "ACGT" * 17), ("b", "TTTTGGGG")]
        out = tmp_path / "o.fa"
        assert write_fasta(str(out), records, wrap=11) == 2
        assert [(n, s) for n, _d, s in read_fasta(str(out))] == records


# ---------------------------------------------------------------------------
# FASTG
# ---------------------------------------------------------------------------

E1 = "EDGE_1_length_10_cov_5.0"
E2 = "EDGE_2_length_8_cov_4.0"
E3 = "EDGE_3_length_6_cov_3.0"
E4 = "EDGE_4_length_5_cov_2.0"
S1, S2, S3, S4 = "AAAACCCCGG", "TTTTGGGG", "ACGTAC", "GGCTA"


@pytest.fixture
def spades_fastg(tmp_path):
    """A hand-built SPAdes-style FASTG.

    ``EDGE_1`` connects forward to the reverse copy of ``EDGE_2`` and to the
    forward copy of ``EDGE_3``; ``EDGE_2``'s record states the same connection
    from its own side, which must fold into one link. ``EDGE_3`` and ``EDGE_4``
    appear only as their apostrophe (reverse-complement) records.
    """
    path = tmp_path / "assembly_graph.fastg"
    path.write_text(
        f">{E1}:{E2}',{E3};\n{S1}\n"
        f">{E2}:{E1}';\n{S2}\n"
        f">{E3}';\n{revcomp(S3)}\n"
        f">{E4}';\n{S4}\n"
    )
    return path


class TestFastg:
    def test_one_segment_per_edge(self, spades_fastg):
        graph = read_fastg(str(spades_fastg))
        assert sorted(graph.segments) == sorted([E1, E2, E3, E4])
        assert graph.source_format == "fastg"

    def test_forward_records_keep_their_sequence(self, spades_fastg):
        graph = read_fastg(str(spades_fastg))
        assert graph[E1].sequence == S1
        assert graph[E2].sequence == S2

    def test_apostrophe_records_are_folded_back_by_revcomp(self, spades_fastg):
        graph = read_fastg(str(spades_fastg))
        # The file only held EDGE_3' and EDGE_4', so the stored + strand is the
        # reverse complement of what was written.
        assert graph[E3].sequence == S3
        assert graph[E4].sequence == revcomp(S4)
        assert graph[E4].length == len(S4)

    def test_links_and_their_orientations(self, spades_fastg):
        graph = read_fastg(str(spades_fastg))
        # EDGE_1 -> EDGE_2' and EDGE_1 -> EDGE_3; the EDGE_2 record repeats the
        # first of those from the other side and must not create a second link.
        assert graph.link_count == 2
        assert graph.has_link(E1, "+", E2, "-")
        assert graph.has_link(E2, "+", E1, "-")  # the same physical link
        assert graph.has_link(E1, "+", E3, "+")
        assert not graph.has_link(E1, "+", E2, "+")

    def test_depth_comes_from_the_name(self, spades_fastg):
        graph = read_fastg(str(spades_fastg))
        assert graph[E1].depth == pytest.approx(5.0)
        assert graph[E4].depth == pytest.approx(2.0)

    def test_plain_fasta_is_rejected(self, tmp_path):
        src = write(tmp_path / "plain.fastg", ">edge_1\nACGT\n>edge_2\nTTTT\n")
        with pytest.raises(AssemblageFormatError, match="looks like plain FASTA"):
            read_fastg(str(src))


# ---------------------------------------------------------------------------
# AGP
# ---------------------------------------------------------------------------


class TestAgp:
    ROWS = [
        AgpRow("scaf1", 1, 1000, 1, COMPONENT_CONTIG, "ctgA", 1, 1000, "+"),
        AgpRow(
            "scaf1",
            1001,
            1100,
            2,
            GAP_KNOWN,
            gap_length=100,
            gap_type="scaffold",
            linkage="yes",
            linkage_evidence="align_genus",
        ),
        AgpRow("scaf1", 1101, 1600, 3, COMPONENT_CONTIG, "ctgB", 1, 500, "-"),
    ]

    def test_write_then_read_preserves_every_value(self, tmp_path):
        out = tmp_path / "s.agp"
        assert write_agp(str(out), self.ROWS) == 3
        back = list(read_agp(str(out)))
        assert back == self.ROWS

    def test_header_and_comments_are_written_and_skipped_on_read(self, tmp_path):
        out = tmp_path / "s.agp"
        write_agp(str(out), self.ROWS, comments=["built by tests", "two lines\nhere"])
        text = out.read_text()
        assert text.startswith("##agp-version\t2.1\n")
        assert "# built by tests\n" in text
        assert "# two lines\n# here\n" in text
        assert list(read_agp(str(out))) == self.ROWS

    def test_row_columns(self):
        assert self.ROWS[0].to_line() == "scaf1\t1\t1000\t1\tW\tctgA\t1\t1000\t+"
        assert self.ROWS[1].to_line() == "scaf1\t1001\t1100\t2\tN\t100\tscaffold\tyes\talign_genus"

    def test_short_row_raises_with_a_line_number(self, tmp_path):
        out = tmp_path / "bad.agp"
        out.write_text("##agp-version\t2.1\nscaf1\t1\t1000\t1\tW\tctgA\n")
        with pytest.raises(AssemblageFormatError) as exc:
            list(read_agp(str(out)))
        assert exc.value.line_no == 2
        assert "9 columns" in exc.value.message

    def test_non_integer_coordinate_raises(self, tmp_path):
        out = tmp_path / "bad.agp"
        out.write_text("scaf1\tone\t1000\t1\tW\tctgA\t1\t1000\t+\n")
        with pytest.raises(AssemblageFormatError, match="must be integers"):
            list(read_agp(str(out)))


# ---------------------------------------------------------------------------
# Format detection and the loader
# ---------------------------------------------------------------------------


class TestDetectFormat:
    def test_gfa(self, demo_paths):
        assert detect_format(str(demo_paths["gfa"])) == "gfa"

    def test_fasta(self, demo_paths):
        assert detect_format(str(demo_paths["reference"])) == "fasta"
        assert detect_format(str(demo_paths["contigs"])) == "fasta"

    def test_fastg(self, spades_fastg):
        assert detect_format(str(spades_fastg)) == "fastg"

    def test_gzipped_fasta(self, tmp_path):
        src = tmp_path / "z.fasta.gz"
        with gzip.open(src, "wt") as handle:
            handle.write(">a\nACGTACGT\n")
        assert detect_format(str(src)) == "fasta"

    def test_gfa2(self, tmp_path):
        src = write(
            tmp_path / "x.gfa",
            """
            H\tVN:Z:2.0
            S\ts1\t10\tACGTACGTAC
            """,
        )
        assert detect_format(str(src)) == "gfa2"

    def test_unknown(self, tmp_path):
        src = write(tmp_path / "notes.txt", "just some prose\nwith no markers\n")
        assert detect_format(str(src)) == "unknown"


class TestLoader:
    def test_load_graph_sniffs_gfa(self, demo_paths):
        graph = load_graph(str(demo_paths["gfa"]))
        assert graph.source_format == "gfa"
        assert len(graph.segments) == 10
        assert graph.link_count == 9

    def test_load_graph_sniffs_fasta_into_an_edge_free_graph(self, demo_paths):
        graph = load_graph(str(demo_paths["contigs"]))
        assert graph.source_format == "fasta"
        assert len(graph.segments) == 10
        assert graph.link_count == 0
        assert graph.total_length == 137_400

    def test_contig_fasta_and_gfa_agree_on_sequence(self, demo_paths):
        from_gfa = load_graph(str(demo_paths["gfa"]))
        from_fa = read_contigs_fasta(str(demo_paths["contigs"]))
        assert {n: s.sequence for n, s in from_gfa.segments.items()} == {
            n: s.sequence for n, s in from_fa.segments.items()
        }

    def test_missing_file(self, tmp_path):
        with pytest.raises(AssemblageFormatError, match="file not found"):
            load_graph(str(tmp_path / "nope.gfa"))

    def test_unknown_format(self, tmp_path):
        src = write(tmp_path / "notes.txt", "just some prose\n")
        with pytest.raises(AssemblageFormatError, match="could not determine the format"):
            load_graph(str(src))

    def test_explicit_format_overrides_sniffing(self, demo_paths):
        graph = load_graph(str(demo_paths["contigs"]), fmt="fasta")
        assert graph.source_format == "fasta"


class TestOverlapInference:
    """FASTG states no overlap, but SPAdes edges really do overlap by k-1 bases.

    Taking them as blunt joins silently corrupts every merge and every scaffold
    built through such a link: the shared bases appear twice. Validated against
    a real SPAdes FASTG where Bandage reports a 127 bp overlap.
    """

    @staticmethod
    def _overlapping_graph(overlap_bp, n=6):
        """A chain whose neighbours share ``overlap_bp`` bases, links marked blunt."""
        import random

        from assemblage.core.model import AssemblyGraph, Link, Segment

        rng = random.Random(11)
        g = AssemblyGraph("inferred")
        previous = None
        for i in range(n):
            if previous is None:
                seq = "".join(rng.choice("ACGT") for _ in range(400))
            else:
                seq = previous[-overlap_bp:] + "".join(
                    rng.choice("ACGT") for _ in range(400 - overlap_bp)
                )
            g.add_segment(Segment(f"s{i}", seq))
            if i:
                g.add_link(Link(f"s{i - 1}", "+", f"s{i}", "+", 0, "*"))
            previous = seq
        return g

    def test_exact_overlap_finds_the_shared_bases(self):
        from assemblage.core.io.overlaps import exact_overlap

        assert exact_overlap("AAAACCCGGG", "CCCGGGTTTT") == 6
        assert exact_overlap("AAAA", "TTTT") == 0
        assert exact_overlap("ACGT", "ACGT") == 4

    def test_a_real_overlap_is_measured_and_applied(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps

        g = self._overlapping_graph(127)
        assert all(link.overlap == 0 for link in g.links.values())
        assert apply_inferred_overlaps(g) == 127
        assert {link.overlap for link in g.links.values()} == {127}
        assert g.overlap_default == 127

    def test_the_merged_sequence_no_longer_duplicates_the_junction(self):
        from assemblage.core.analysis.operations import merge_path
        from assemblage.core.io.overlaps import apply_inferred_overlaps

        g = self._overlapping_graph(127, n=3)
        apply_inferred_overlaps(g)
        merged = merge_path(g, [("s0", "+"), ("s1", "+"), ("s2", "+")])
        segment = g.segments[merged["new_name"]]
        # 3 x 400 bp sharing 127 bp at each of 2 junctions
        assert segment.length == 3 * 400 - 2 * 127

    def test_genuinely_blunt_joins_are_left_alone(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps

        g = self._overlapping_graph(0)
        assert apply_inferred_overlaps(g) == 0
        assert {link.overlap for link in g.links.values()} == {0}

    def test_a_coincidental_short_match_is_not_treated_as_an_overlap(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps

        # Neighbours share only a couple of bases by chance.
        g = self._overlapping_graph(3)
        assert apply_inferred_overlaps(g) == 0

    def test_explicit_cigars_are_never_overridden(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps
        from assemblage.core.model import Link

        g = self._overlapping_graph(127)
        for key, link in list(g.links.items()):
            g.links[key] = Link(
                link.from_name, link.from_orient, link.to_name, link.to_orient, 55, "55M"
            )
        assert apply_inferred_overlaps(g) == 0
        assert {link.overlap for link in g.links.values()} == {55}

    def test_a_graph_without_sequences_cannot_be_measured(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps
        from assemblage.core.model import AssemblyGraph, Link, Segment

        g = AssemblyGraph("noseq")
        g.add_segment(Segment("a", None, length=500))
        g.add_segment(Segment("b", None, length=500))
        g.add_link(Link("a", "+", "b", "+", 0, "*"))
        assert apply_inferred_overlaps(g) == 0


class TestMixedOverlapGraphs:
    """A measured overlap must be checked against each link, not stamped on all.

    A graph where some joins overlap and others are blunt would otherwise have
    its blunt links given a phantom overlap, silently deleting real bases from
    every merge and scaffold through them.
    """

    @staticmethod
    def _mixed(n_overlapping=20, n_blunt=6, overlap=55):
        import random

        from assemblage.core.model import AssemblyGraph, Link, Segment

        rng = random.Random(5)
        seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        g = AssemblyGraph("mixed")
        previous = None
        for i in range(n_overlapping):
            s = seq(400) if previous is None else previous[-overlap:] + seq(400 - overlap)
            g.add_segment(Segment(f"o{i}", s))
            if i:
                g.add_link(Link(f"o{i - 1}", "+", f"o{i}", "+", 0, "*"))
            previous = s
        for i in range(n_blunt):
            g.add_segment(Segment(f"b{i}", seq(400)))
            if i:
                g.add_link(Link(f"b{i - 1}", "+", f"b{i}", "+", 0, "*"))
        return g

    def test_only_the_links_that_really_overlap_get_the_overlap(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps

        g = self._mixed()
        assert apply_inferred_overlaps(g) == 55
        overlapping = [lk for k, lk in g.links.items() if k[0].startswith("o")]
        blunt = [lk for k, lk in g.links.items() if k[0].startswith("b")]
        assert {lk.overlap for lk in overlapping} == {55}
        assert {lk.overlap for lk in blunt} == {0}

    def test_walking_a_blunt_join_keeps_every_base(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps

        g = self._mixed()
        apply_inferred_overlaps(g)
        assert len(g.walk_sequence([("b0", "+"), ("b1", "+")])) == 800

    def test_a_single_coincidental_match_does_not_become_the_default(self):
        from assemblage.core.io.overlaps import apply_inferred_overlaps
        from assemblage.core.model import AssemblyGraph, Link, Segment

        # One measurable link among many that carry no sequence.
        g = AssemblyGraph("sparse")
        shared = "ACGTACGTACGTAC"  # 14 bp
        g.add_segment(Segment("a", "TTTT" * 25 + shared))
        g.add_segment(Segment("b", shared + "GGGG" * 25))
        g.add_link(Link("a", "+", "b", "+", 0, "*"))
        for i in range(20):
            g.add_segment(Segment(f"n{i}", None, length=500))
            if i:
                g.add_link(Link(f"n{i - 1}", "+", f"n{i}", "+", 0, "*"))
        apply_inferred_overlaps(g)
        assert g.overlap_default == 0, (
            "one verifiable link must not set the graph-wide fallback overlap"
        )
