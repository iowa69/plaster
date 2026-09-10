"""Compare two assemblies by content: what is in one and not the other.

The side-by-side statistics in :mod:`metrics` answer "which assembly is
better". They cannot answer "what did I gain, and what did I lose", because two
assemblies with identical N50 and identical total length can still disagree
about a whole plasmid. That question is the one people actually ask when they
reassemble an isolate, change an assembler, or look for mobile elements, and it
can only be answered by aligning the sequence.

Each contig is aligned to the *other* assembly and scored by how much of it the
other one covers. The interesting output is not a percentage but a list: the
contigs present in one assembly and missing from the other, and -- more useful
still -- the connected components entirely missing, since a plasmid absent from
one of a pair is a whole component that never aligns.

Coverage is measured per base and not per alignment. A contig hit by six
partial alignments that together cover it is present, and one hit six times
over the same third of itself is not; counting alignments or summing their
lengths gets both of those wrong.
"""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass, field
from typing import Sequence

from ..errors import PlastrError
from ..io.fasta import write_fasta
from ..model import AssemblyGraph
from . import align as align_mod

__all__ = [
    "invert_spans",
    "ContigDiff",
    "ComponentDiff",
    "AssemblySide",
    "GraphDiff",
    "diff_graphs",
    "write_outputs",
    "SHARED_FRACTION",
    "PRESENT_FRACTION",
]

#: A contig this much covered by the other assembly counts as shared outright.
SHARED_FRACTION = 0.95

#: Below this it is reported as missing. Between the two it is partial: real
#: sequence in common, but enough of it unaccounted for to be worth a look.
PRESENT_FRACTION = 0.10


@dataclass
class ContigDiff:
    """One contig, scored against the whole of the other assembly."""

    name: str
    length: int
    depth: float | None
    component: int
    circular: bool
    covered: int
    best_match: str | None
    best_identity: float | None
    status: str  # "shared" | "partial" | "missing"
    #: Half-open [start, end) runs of this contig that nothing in the other
    #: assembly aligns to. This is the sequence that is actually new, and it is
    #: what gets written out -- a partial contig's novel insert is the
    #: interesting part of it, not the whole contig.
    uncovered_spans: list[tuple[int, int]] = field(default_factory=list)

    @property
    def covered_fraction(self) -> float:
        return (self.covered / self.length) if self.length else 0.0

    @property
    def uncovered(self) -> int:
        return max(0, self.length - self.covered)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "length": self.length,
            "depth": self.depth,
            "component": self.component,
            "circular": self.circular,
            "covered": self.covered,
            "covered_fraction": self.covered_fraction,
            "best_match": self.best_match,
            "best_identity": self.best_identity,
            "status": self.status,
            "uncovered_spans": [list(sp) for sp in self.uncovered_spans],
        }


@dataclass
class ComponentDiff:
    """A connected component summarised by how much of it the other side has.

    A component every one of whose contigs is missing is the headline result:
    on a bacterial pair that is a plasmid one isolate carries and the other does
    not, and it is invisible in any table of N50s.
    """

    index: int
    contigs: list[str]
    length: int
    covered: int
    circular: bool
    status: str  # "shared" | "partial" | "missing"

    @property
    def covered_fraction(self) -> float:
        return (self.covered / self.length) if self.length else 0.0

    def to_dict(self) -> dict:
        return {
            "index": self.index,
            "contigs": self.contigs,
            "length": self.length,
            "covered": self.covered,
            "covered_fraction": self.covered_fraction,
            "circular": self.circular,
            "status": self.status,
        }


@dataclass
class AssemblySide:
    """One assembly's view of the comparison."""

    label: str
    path: str | None
    total_length: int
    num_contigs: int
    contigs: list[ContigDiff] = field(default_factory=list)
    components: list[ComponentDiff] = field(default_factory=list)
    #: Contigs with no sequence, which cannot be compared at all.
    without_sequence: int = 0
    #: (this contig, other contig, aligned bases), heaviest first. The edges of
    #: the entanglement map: which contig over there accounts for which contig
    #: over here.
    matches: list[tuple[str, str, int]] = field(default_factory=list)

    def _bases(self, status: str) -> int:
        return sum(c.length for c in self.contigs if c.status == status)

    @property
    def shared_bases(self) -> int:
        """Bases of this assembly the other one covers, summed per base."""
        return sum(c.covered for c in self.contigs)

    @property
    def unique_bases(self) -> int:
        return sum(c.uncovered for c in self.contigs)

    @property
    def missing_bases(self) -> int:
        """Bases held in contigs the other assembly does not have at all."""
        return self._bases("missing")

    @property
    def partial_bases(self) -> int:
        return self._bases("partial")

    @property
    def shared_fraction(self) -> float:
        return (self.shared_bases / self.total_length) if self.total_length else 0.0

    def counts(self) -> dict[str, int]:
        out = {"shared": 0, "partial": 0, "missing": 0}
        for c in self.contigs:
            out[c.status] = out.get(c.status, 0) + 1
        return out

    def missing_contigs(self) -> list[ContigDiff]:
        """Missing and partial contigs, longest first -- the things to look at."""
        rows = [c for c in self.contigs if c.status != "shared"]
        rows.sort(key=lambda c: (c.status == "partial", -c.length))
        return rows

    def missing_components(self) -> list[ComponentDiff]:
        rows = [c for c in self.components if c.status != "shared"]
        rows.sort(key=lambda c: (c.status == "partial", -c.length))
        return rows

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "path": self.path,
            "total_length": self.total_length,
            "num_contigs": self.num_contigs,
            "without_sequence": self.without_sequence,
            "shared_bases": self.shared_bases,
            "unique_bases": self.unique_bases,
            "missing_bases": self.missing_bases,
            "partial_bases": self.partial_bases,
            "shared_fraction": self.shared_fraction,
            "counts": self.counts(),
            "contigs": [c.to_dict() for c in self.contigs],
            "components": [c.to_dict() for c in self.components],
            "matches": [list(m) for m in self.matches],
        }


@dataclass
class GraphDiff:
    """Both directions of the comparison. Neither alone is the answer."""

    a: AssemblySide
    b: AssemblySide
    preset: str
    min_identity: float
    shared_fraction: float
    present_fraction: float

    def to_dict(self) -> dict:
        return {
            "a": self.a.to_dict(),
            "b": self.b.to_dict(),
            "preset": self.preset,
            "min_identity": self.min_identity,
            "shared_fraction": self.shared_fraction,
            "present_fraction": self.present_fraction,
        }


def _merge_spans(spans: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
    """Overlapping intervals collapsed into disjoint ones, in order."""
    if not spans:
        return []
    ordered = sorted(spans)
    out = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start > out[-1][1]:
            out.append([start, end])
        else:
            out[-1][1] = max(out[-1][1], end)
    return [(a, b) for a, b in out]


def _merged_length(spans: Sequence[tuple[int, int]]) -> int:
    """Total length of a set of intervals, counting overlaps once."""
    return sum(b - a for a, b in _merge_spans(spans))


def invert_spans(spans: Sequence[tuple[int, int]], length: int) -> list[tuple[int, int]]:
    """The gaps between a set of intervals, over [0, length)."""
    out: list[tuple[int, int]] = []
    cursor = 0
    for start, end in _merge_spans(spans):
        start = max(0, min(start, length))
        end = max(0, min(end, length))
        if start > cursor:
            out.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < length:
        out.append((cursor, length))
    return out


def _write_contigs(graph: AssemblyGraph, path: str) -> int:
    records = [
        (name, seg.sequence)
        for name, seg in graph.segments.items()
        if seg.sequence
    ]
    write_fasta(path, records)
    return len(records)


def _classify(fraction: float, shared_at: float, present_at: float) -> str:
    if fraction >= shared_at:
        return "shared"
    if fraction >= present_at:
        return "partial"
    return "missing"


def _score_side(
    graph: AssemblyGraph,
    label: str,
    path: str | None,
    other_fasta: str,
    preset: str,
    min_identity: float,
    threads: int,
    shared_at: float,
    present_at: float,
) -> AssemblySide:
    """Align every contig of ``graph`` to ``other_fasta`` and score coverage."""
    queries = [(n, s.sequence) for n, s in graph.segments.items() if s.sequence]
    hits = []
    if queries:
        hits = align_mod.align_sequences(
            queries,
            other_fasta,
            preset=preset,
            threads=threads,
            # Splitting at long indels would count the span across a deletion as
            # covered, which is the whole point of not doing it here: a contig
            # whose middle is absent from the other assembly must read as
            # partial, not as shared.
            split_long_indels=0,
        )
    if min_identity > 0:
        hits = [h for h in hits if h.identity >= min_identity]

    spans: dict[str, list[tuple[int, int]]] = {}
    best: dict[str, tuple[float, str]] = {}
    weights: dict[tuple[str, str], int] = {}
    for h in hits:
        spans.setdefault(h.query, []).append((h.q_st, h.q_en))
        key = (h.query, h.ref)
        weights[key] = weights.get(key, 0) + max(0, h.q_span)
        score = h.q_span * max(h.identity, 0.0)
        prev = best.get(h.query)
        if prev is None or score > prev[0]:
            best[h.query] = (score, h.ref)

    identity_of: dict[str, float] = {}
    for h in hits:
        cur = identity_of.get(h.query)
        if cur is None or h.identity > cur:
            identity_of[h.query] = h.identity

    components = graph.connected_components()
    comp_index = graph.component_map()

    contigs: list[ContigDiff] = []
    without_sequence = 0
    for name, seg in graph.segments.items():
        if not seg.sequence:
            without_sequence += 1
        merged = _merge_spans(spans.get(name, []))
        covered = min(seg.length, sum(b - a for a, b in merged))
        fraction = (covered / seg.length) if seg.length else 0.0
        contigs.append(
            ContigDiff(
                name=name,
                length=seg.length,
                depth=seg.depth,
                component=comp_index.get(name, 0),
                circular=graph.is_circular(name),
                covered=covered,
                best_match=best.get(name, (0.0, None))[1],
                best_identity=identity_of.get(name),
                status=_classify(fraction, shared_at, present_at),
                uncovered_spans=invert_spans(merged, seg.length),
            )
        )
    contigs.sort(key=lambda c: -c.length)

    by_name = {c.name: c for c in contigs}
    comp_rows: list[ComponentDiff] = []
    for i, members in enumerate(components):
        length = sum(by_name[m].length for m in members if m in by_name)
        covered = sum(by_name[m].covered for m in members if m in by_name)
        fraction = (covered / length) if length else 0.0
        comp_rows.append(
            ComponentDiff(
                index=i,
                contigs=sorted(members),
                length=length,
                covered=covered,
                circular=any(by_name[m].circular for m in members if m in by_name),
                status=_classify(fraction, shared_at, present_at),
            )
        )
    comp_rows.sort(key=lambda c: -c.length)

    return AssemblySide(
        label=label,
        path=path,
        total_length=sum(c.length for c in contigs),
        num_contigs=len(contigs),
        contigs=contigs,
        components=comp_rows,
        without_sequence=without_sequence,
        matches=sorted(
            ((q, r, w) for (q, r), w in weights.items()),
            key=lambda m: -m[2],
        ),
    )


def diff_graphs(
    graph_a: AssemblyGraph,
    graph_b: AssemblyGraph,
    label_a: str = "A",
    label_b: str = "B",
    path_a: str | None = None,
    path_b: str | None = None,
    preset: str = "asm10",
    min_identity: float = 0.0,
    threads: int = 4,
    shared_fraction: float = SHARED_FRACTION,
    present_fraction: float = PRESENT_FRACTION,
) -> GraphDiff:
    """Compare two assemblies by content, in both directions.

    Both directions are computed because neither answers the question alone: a
    contig missing from B tells you nothing about what B has that A does not,
    and a collapsed repeat can be fully covered one way and only partly the
    other.
    """
    if shared_fraction <= present_fraction:
        raise PlastrError(
            "shared_fraction must be greater than present_fraction "
            f"(got {shared_fraction} and {present_fraction})"
        )
    for graph, label in ((graph_a, label_a), (graph_b, label_b)):
        if not any(s.sequence for s in graph.segments.values()):
            raise PlastrError(
                f"{label} has no sequences, so there is nothing to compare. "
                "Load a GFA with S-line sequences, or a FASTA."
            )

    with tempfile.TemporaryDirectory(prefix="plastr-diff-") as tmp:
        fa_a = os.path.join(tmp, "a.fasta")
        fa_b = os.path.join(tmp, "b.fasta")
        _write_contigs(graph_a, fa_a)
        _write_contigs(graph_b, fa_b)

        side_a = _score_side(
            graph_a, label_a, path_a, fa_b, preset, min_identity, threads,
            shared_fraction, present_fraction,
        )
        side_b = _score_side(
            graph_b, label_b, path_b, fa_a, preset, min_identity, threads,
            shared_fraction, present_fraction,
        )

    return GraphDiff(
        a=side_a,
        b=side_b,
        preset=preset,
        min_identity=min_identity,
        shared_fraction=shared_fraction,
        present_fraction=present_fraction,
    )


def write_outputs(
    diff: "GraphDiff",
    graph_a: AssemblyGraph,
    graph_b: AssemblyGraph,
    outdir: str | os.PathLike[str],
    overwrite: bool = False,
    min_span: int = 200,
) -> list[tuple[str, str, int]]:
    """Write the comparison out as files you can feed to the next tool.

    Four things per assembly, because "what is missing" has two useful readings
    and both are wanted:

      * ``<label>.only-contigs.fasta`` -- whole contigs the other assembly does
        not have. This is the plasmid, the phage, the insertion element.
      * ``<label>.only-regions.fasta`` -- just the stretches nothing over there
        aligns to, cut out of their contigs and named with their coordinates.
        A contig that is 90% shared contributes only its novel 10% here, which
        is the part worth blasting.
      * ``<label>.contigs.tsv`` -- every contig with its verdict, coverage, best
        match and depth, for sorting in a spreadsheet or reading with pandas.
      * ``<label>.components.tsv`` -- the same by connected component, which is
        where a whole missing molecule shows up.

    Returns ``(kind, path, count)`` for each file written.
    """
    outdir = str(outdir)
    os.makedirs(outdir, exist_ok=True)
    written: list[tuple[str, str, int]] = []

    def _open(name: str):
        path = os.path.join(outdir, name)
        if os.path.exists(path) and not overwrite:
            raise PlastrError(
                f"{path} already exists; pass overwrite to replace it"
            )
        return path

    def _safe(label: str) -> str:
        return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in label)

    for side, graph in ((diff.a, graph_a), (diff.b, graph_b)):
        stem = _safe(side.label)

        whole = [
            (c.name, graph.segments[c.name].sequence)
            for c in side.contigs
            if c.status == "missing" and graph.segments.get(c.name) is not None
            and graph.segments[c.name].sequence
        ]
        path = _open(f"{stem}.only-contigs.fasta")
        write_fasta(path, whole)
        written.append(("only-contigs", path, len(whole)))

        regions: list[tuple[str, str]] = []
        for c in side.contigs:
            seg = graph.segments.get(c.name)
            if seg is None or not seg.sequence:
                continue
            for start, end in c.uncovered_spans:
                if end - start < min_span:
                    continue
                # Coordinates in the name, 1-based inclusive, so the record can
                # be traced back to the contig it was cut from.
                regions.append(
                    (f"{c.name}:{start + 1}-{end} len={end - start}", seg.sequence[start:end])
                )
        path = _open(f"{stem}.only-regions.fasta")
        write_fasta(path, regions)
        written.append(("only-regions", path, len(regions)))

        path = _open(f"{stem}.contigs.tsv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(
                "contig\tlength\tstatus\tcovered_bp\tcovered_fraction\tuncovered_bp"
                "\tbest_match\tbest_identity\tdepth\tcomponent\tcircular\n"
            )
            for c in side.contigs:
                fh.write(
                    f"{c.name}\t{c.length}\t{c.status}\t{c.covered}\t"
                    f"{c.covered_fraction:.4f}\t{c.uncovered}\t"
                    f"{c.best_match or ''}\t"
                    f"{'' if c.best_identity is None else f'{c.best_identity:.4f}'}\t"
                    f"{'' if c.depth is None else f'{c.depth:.2f}'}\t"
                    f"{c.component}\t{int(c.circular)}\n"
                )
        written.append(("contigs-tsv", path, len(side.contigs)))

        path = _open(f"{stem}.components.tsv")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("component\tcontigs\tlength\tstatus\tcovered_fraction\tcircular\tmembers\n")
            for c in side.components:
                fh.write(
                    f"{c.index}\t{len(c.contigs)}\t{c.length}\t{c.status}\t"
                    f"{c.covered_fraction:.4f}\t{int(c.circular)}\t{','.join(c.contigs)}\n"
                )
        written.append(("components-tsv", path, len(side.components)))

    return written
