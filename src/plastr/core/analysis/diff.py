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
    "ContigDiff",
    "ComponentDiff",
    "AssemblySide",
    "GraphDiff",
    "diff_graphs",
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


def _merged_length(spans: Sequence[tuple[int, int]]) -> int:
    """Total length of a set of intervals, counting overlaps once."""
    if not spans:
        return 0
    ordered = sorted(spans)
    total = 0
    cur_start, cur_end = ordered[0]
    for start, end in ordered[1:]:
        if start > cur_end:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    total += cur_end - cur_start
    return total


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
    for h in hits:
        spans.setdefault(h.query, []).append((h.q_st, h.q_en))
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
        covered = min(seg.length, _merged_length(spans.get(name, [])))
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
