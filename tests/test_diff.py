"""Comparing two assemblies by content."""

from __future__ import annotations

import random

import pytest

from plastr.core.analysis import align as align_mod
from plastr.core.analysis.diff import (
    _merged_length,
    diff_graphs,
)
from plastr.core.errors import PlastrError
from plastr.core.model import AssemblyGraph, Link, Segment

needs_aligner = pytest.mark.skipif(
    align_mod.alignment_backend() == "none",
    reason="needs minimap2 or mappy",
)


def _seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mutate(rng: random.Random, seq: str, rate: float) -> str:
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def _graph(segments, links=()) -> AssemblyGraph:
    g = AssemblyGraph()
    for name, seq in segments:
        g.add_segment(Segment(name=name, sequence=seq, depth=30.0))
    for a, b in links:
        g.add_link(Link(from_name=a, from_orient="+", to_name=b, to_orient="+"))
    return g


# --------------------------------------------------------------- interval maths


def test_merged_length_counts_overlaps_once():
    # The whole point of measuring coverage per base: six alignments piled on
    # the same third of a contig must not add up to the contig being present.
    assert _merged_length([]) == 0
    assert _merged_length([(0, 100)]) == 100
    assert _merged_length([(0, 100), (0, 100), (0, 100)]) == 100
    assert _merged_length([(0, 50), (50, 100)]) == 100
    assert _merged_length([(0, 60), (40, 100)]) == 100
    assert _merged_length([(0, 10), (90, 100)]) == 20
    # Order must not matter.
    assert _merged_length([(90, 100), (0, 10), (5, 12)]) == 22


# ------------------------------------------------------------------ guard rails


def test_rejects_a_graph_with_no_sequences():
    a = AssemblyGraph()
    rng = random.Random(0)
    b = (_graph([("b", _seq(rng, 500))]))
    with pytest.raises(PlastrError, match="no sequences"):
        diff_graphs(a, b)


def test_rejects_nonsense_thresholds():
    rng = random.Random(0)
    a = (_graph([("a", _seq(rng, 500))]))
    with pytest.raises(PlastrError, match="greater than"):
        diff_graphs(a, a, shared_fraction=0.1, present_fraction=0.9)


# -------------------------------------------------------------- the comparison


@needs_aligner
def test_identical_assemblies_are_wholly_shared():
    rng = random.Random(1)
    segs = [("c1", _seq(rng, 20_000)), ("c2", _seq(rng, 12_000))]
    a = _graph(segs)
    b = _graph(segs)

    d = diff_graphs(a, b, "first", "second", preset="asm5")

    assert [c.status for c in d.a.contigs] == ["shared", "shared"]
    assert [c.status for c in d.b.contigs] == ["shared", "shared"]
    assert d.a.missing_bases == 0
    assert d.b.missing_bases == 0
    assert d.a.shared_fraction > 0.99


@needs_aligner
def test_a_contig_absent_from_the_other_assembly_is_missing():
    rng = random.Random(2)
    shared = _seq(rng, 25_000)
    only_in_a = _seq(rng, 9_000)

    a = (_graph([("shared", shared), ("extra", only_in_a)]))
    b = (_graph([("shared", shared)]))

    d = diff_graphs(a, b, "A", "B", preset="asm5")

    by_name = {c.name: c for c in d.a.contigs}
    assert by_name["shared"].status == "shared"
    assert by_name["extra"].status == "missing"
    assert by_name["extra"].covered == 0
    assert by_name["extra"].best_match is None

    assert d.a.missing_bases == 9_000
    # B is wholly contained in A, so nothing is missing the other way.
    assert d.b.missing_bases == 0
    assert d.a.counts() == {"shared": 1, "partial": 0, "missing": 1}


@needs_aligner
def test_a_whole_plasmid_missing_shows_up_as_a_missing_component():
    # The headline case: two isolates of one strain, one carrying a plasmid.
    # No table of N50s shows this; a component that never aligns does.
    rng = random.Random(3)
    chrom = [(f"chr{i}", _seq(rng, 30_000)) for i in range(4)]
    chrom_links = [(f"chr{i}", f"chr{i+1}") for i in range(3)]
    plasmid = [("plasmid", _seq(rng, 40_000))]

    a = (_graph(chrom + plasmid, chrom_links))
    b = (_graph(chrom, chrom_links))

    d = diff_graphs(a, b, "with_plasmid", "without", preset="asm5")

    missing = d.a.missing_components()
    assert len(missing) == 1
    assert missing[0].contigs == ["plasmid"]
    assert missing[0].status == "missing"
    assert missing[0].length == 40_000
    # And the chromosome component is not flagged.
    shared = [c for c in d.a.components if c.status == "shared"]
    assert sum(len(c.contigs) for c in shared) == 4


@needs_aligner
def test_a_contig_half_present_is_partial_not_shared():
    # A contig whose second half the other assembly lacks. Splitting alignments
    # at long indels would count the gap as covered and call this shared, which
    # is why the comparison does not split them.
    rng = random.Random(4)
    left = _seq(rng, 20_000)
    right = _seq(rng, 20_000)

    a = (_graph([("joined", left + right)]))
    b = (_graph([("left_only", left)]))

    d = diff_graphs(a, b, "A", "B", preset="asm5")
    contig = d.a.contigs[0]
    assert contig.status == "partial"
    assert 0.35 < contig.covered_fraction < 0.65
    assert contig.uncovered > 15_000


@needs_aligner
def test_divergent_copies_still_count_as_shared():
    # Same molecule, 1% different. It is the same plasmid and must not be
    # reported as newly gained.
    rng = random.Random(5)
    original = _seq(rng, 30_000)
    a = (_graph([("p", original)]))
    b = (_graph([("p", _mutate(rng, original, 0.01))]))

    d = diff_graphs(a, b, "A", "B", preset="asm10")
    assert d.a.contigs[0].status == "shared"
    assert d.a.contigs[0].best_identity is not None
    assert d.a.contigs[0].best_identity > 0.95


@needs_aligner
def test_direction_matters_and_both_are_reported():
    rng = random.Random(6)
    common = _seq(rng, 25_000)
    a = (_graph([("common", common), ("only_a", _seq(rng, 8_000))]))
    b = (_graph([("common", common), ("only_b", _seq(rng, 6_000))]))

    d = diff_graphs(a, b, "A", "B", preset="asm5")

    assert [c.name for c in d.a.missing_contigs()] == ["only_a"]
    assert [c.name for c in d.b.missing_contigs()] == ["only_b"]
    assert d.a.missing_bases == 8_000
    assert d.b.missing_bases == 6_000


@needs_aligner
def test_to_dict_is_serialisable_and_complete():
    import json

    rng = random.Random(7)
    a = (_graph([("x", _seq(rng, 12_000))]))
    b = (_graph([("y", _seq(rng, 12_000))]))
    d = diff_graphs(a, b, "A", "B", preset="asm5")

    blob = json.loads(json.dumps(d.to_dict()))
    assert blob["a"]["label"] == "A"
    assert blob["a"]["counts"]["missing"] == 1
    assert blob["preset"] == "asm5"
    assert "contigs" in blob["b"] and "components" in blob["b"]


# ------------------------------------------------------------ span arithmetic


def test_invert_spans_is_the_complement():
    from plastr.core.analysis.diff import invert_spans

    assert invert_spans([], 100) == [(0, 100)]
    assert invert_spans([(0, 100)], 100) == []
    assert invert_spans([(10, 20)], 100) == [(0, 10), (20, 100)]
    assert invert_spans([(0, 30), (70, 100)], 100) == [(30, 70)]
    # Overlapping and unsorted input must still give a clean complement.
    assert invert_spans([(70, 100), (0, 30), (10, 25)], 100) == [(30, 70)]
    # Spans past the end are clamped rather than producing negative gaps.
    assert invert_spans([(90, 500)], 100) == [(0, 90)]


# ------------------------------------------------------------- recovering it


@needs_aligner
def test_write_outputs_recovers_the_missing_sequence(tmp_path):
    rng = random.Random(8)
    shared = _seq(rng, 30_000)
    plasmid = _seq(rng, 12_000)

    a = _graph([("chr", shared), ("plasmid", plasmid)])
    b = _graph([("chr", shared)])
    d = diff_graphs(a, b, "A", "B", preset="asm5")

    from plastr.core.analysis.diff import write_outputs

    written = write_outputs(d, a, b, tmp_path)
    kinds = {kind for kind, _, _ in written}
    assert kinds == {"only-contigs", "only-regions", "contigs-tsv", "components-tsv"}

    # The plasmid comes back, sequence intact, so it can go straight into BLAST.
    text = (tmp_path / "A.only-contigs.fasta").read_text()
    assert ">plasmid" in text
    assert "".join(text.split("\n")[1:]).strip() == plasmid
    # And nothing from the chromosome, which B has.
    assert ">chr" not in text

    # B lost nothing, so its file exists and is empty rather than missing.
    assert (tmp_path / "B.only-contigs.fasta").read_text().strip() == ""

    # The regions file names its coordinates.
    regions = (tmp_path / "A.only-regions.fasta").read_text()
    assert "plasmid:1-12000" in regions

    rows = (tmp_path / "A.contigs.tsv").read_text().strip().split("\n")
    assert rows[0].startswith("contig\tlength\tstatus")
    assert len(rows) == 3  # header plus two contigs
    assert any(r.startswith("plasmid\t12000\tmissing") for r in rows)


@needs_aligner
def test_only_regions_carries_the_novel_part_of_a_partial_contig(tmp_path):
    # The point of the regions file: a contig that is mostly shared contributes
    # only its novel insert, not the whole contig.
    rng = random.Random(9)
    left = _seq(rng, 20_000)
    insert = _seq(rng, 5_000)

    a = _graph([("withinsert", left + insert)])
    b = _graph([("plain", left)])
    d = diff_graphs(a, b, "A", "B", preset="asm5")

    from plastr.core.analysis.diff import write_outputs

    write_outputs(d, a, b, tmp_path)
    # Partial, so it is not in only-contigs ...
    assert (tmp_path / "A.only-contigs.fasta").read_text().strip() == ""
    # ... but its uncovered tail is in only-regions, and it is the insert.
    regions = (tmp_path / "A.only-regions.fasta").read_text()
    body = "".join(line for line in regions.split("\n") if not line.startswith(">"))
    assert 4_000 < len(body) < 6_000
    assert body in (left + insert)


@needs_aligner
def test_write_outputs_refuses_to_clobber(tmp_path):
    rng = random.Random(10)
    a = _graph([("x", _seq(rng, 9_000))])
    b = _graph([("y", _seq(rng, 9_000))])
    d = diff_graphs(a, b, "A", "B", preset="asm5")

    from plastr.core.analysis.diff import write_outputs

    write_outputs(d, a, b, tmp_path)
    with pytest.raises(PlastrError, match="already exists"):
        write_outputs(d, a, b, tmp_path)
    write_outputs(d, a, b, tmp_path, overwrite=True)


@needs_aligner
def test_matches_are_recorded_for_the_map(tmp_path):
    rng = random.Random(11)
    shared = _seq(rng, 25_000)
    a = _graph([("a1", shared), ("a2", _seq(rng, 7_000))])
    b = _graph([("b1", shared)])
    d = diff_graphs(a, b, "A", "B", preset="asm5")

    pairs = {(q, r) for q, r, _ in d.a.matches}
    assert ("a1", "b1") in pairs
    # The unmatched contig has no edge at all, which is what leaves it hanging
    # unconnected in the entanglement map.
    assert not any(q == "a2" for q, _, _ in d.a.matches)
    assert all(w > 0 for _, _, w in d.a.matches)
