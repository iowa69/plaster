"""Reference-guided ordering and orientation of contigs.

The method is the standard one (align, order by reference coordinate, orient by
strand, size gaps from reference distance) with two refinements that matter in
practice:

* Contigs whose reference placements overlap substantially are not silently
  stacked. The longer alignment keeps the position; the other is reported as
  redundant, because a contig that lands on top of another is usually a
  collapsed repeat or a duplicated haplotype and quietly interleaving them
  produces a scaffold nobody can trust.
* Where the assembly graph connects two neighbouring contigs by a path of about
  the right length, the gap is filled with that real sequence instead of Ns.
  See :mod:`assemblage.core.scaffold.graph_bridge`.
"""

from __future__ import annotations

from typing import Sequence

from ..analysis.align import Alignment, best_placements
from ..model import AssemblyGraph
from .graph_bridge import bridge_gap
from .plan import (
    GAP_ADJACENT,
    GAP_DEFAULT,
    GAP_GRAPH,
    GAP_REFERENCE,
    Scaffold,
    ScaffoldMember,
    ScaffoldPlan,
)

DEFAULT_MIN_GAP = 100


def scaffold_by_reference(
    graph: AssemblyGraph,
    alignments: Sequence[Alignment],
    min_identity: float = 0.80,
    min_query_coverage: float = 0.30,
    min_align_length: int = 500,
    min_gap: int = DEFAULT_MIN_GAP,
    max_overlap_fraction: float = 0.5,
    fill_gaps_from_graph: bool = True,
    max_bridge_length: int = 50_000,
    include_unplaced: bool = False,
    scaffold_prefix: str = "scaffold",
) -> ScaffoldPlan:
    """Build a scaffold plan from contig-to-reference alignments."""
    plan = ScaffoldPlan(method="reference")

    usable = [a for a in alignments if a.q_span >= min_align_length]
    placements = best_placements(
        usable, min_query_coverage=min_query_coverage, min_identity=min_identity
    )

    all_contigs = set(graph.segments)
    placed_names = set(placements)
    unplaced = sorted(all_contigs - placed_names)

    # Group by reference sequence and order along it.
    by_ref: dict[str, list[Alignment]] = {}
    for aln in placements.values():
        by_ref.setdefault(aln.ref, []).append(aln)

    for ref in sorted(by_ref, key=lambda r: -sum(a.q_span for a in by_ref[r])):
        hits = sorted(by_ref[ref], key=lambda a: (a.r_st, -a.q_span))
        kept: list[Alignment] = []
        for aln in hits:
            if kept:
                previous = kept[-1]
                overlap = previous.r_en - aln.r_st
                shorter = min(previous.r_span, aln.r_span)
                if shorter > 0 and overlap > max_overlap_fraction * shorter:
                    # Same reference span claimed twice: keep the better one.
                    loser = aln if aln.q_span <= previous.q_span else previous
                    winner = previous if loser is aln else aln
                    plan.redundant.append(loser.query)
                    if loser is previous:
                        kept[-1] = winner
                    continue
            kept.append(aln)

        if not kept:
            continue

        scaffold = Scaffold(
            name=f"{scaffold_prefix}_{_safe(ref)}",
            source="reference",
            reference=ref,
        )
        for i, aln in enumerate(kept):
            member = ScaffoldMember(
                segment=aln.query,
                orientation="+" if aln.strand > 0 else "-",
                ref=ref,
                ref_start=aln.r_st,
                ref_end=aln.r_en,
                identity=round(aln.identity, 5),
            )
            if i + 1 < len(kept):
                nxt = kept[i + 1]
                raw_gap = nxt.r_st - aln.r_en
                if raw_gap <= 0:
                    member.overlaps_previous = True
                    member.gap_after = min_gap
                    member.gap_evidence = GAP_DEFAULT
                else:
                    member.gap_after = max(raw_gap, min_gap)
                    member.gap_evidence = GAP_REFERENCE
            scaffold.members.append(member)

        plan.scaffolds.append(scaffold)

    # Optionally keep unplaced contigs as single-contig scaffolds so that the
    # exported FASTA is a complete assembly rather than a subset.
    if include_unplaced:
        for name in unplaced:
            plan.scaffolds.append(
                Scaffold(
                    name=f"{scaffold_prefix}_unplaced_{_safe(name)}",
                    source="unplaced",
                    members=[ScaffoldMember(segment=name, orientation="+")],
                )
            )
    plan.unplaced = unplaced

    if fill_gaps_from_graph and graph.link_count:
        filled = apply_graph_bridges(graph, plan, max_bridge_length=max_bridge_length)
        if filled:
            plan.notes.append(
                f"{filled} gap(s) replaced with real sequence recovered from the assembly graph"
            )

    plan.notes.append(
        f"{plan.placed_count} contig(s) placed on {len(plan.scaffolds)} scaffold(s); "
        f"{len(unplaced)} unplaced, {len(plan.redundant)} redundant"
    )
    return plan


def apply_graph_bridges(
    graph: AssemblyGraph,
    plan: ScaffoldPlan,
    max_bridge_length: int = 50_000,
    length_tolerance: float = 0.5,
) -> int:
    """Replace Ns with graph-derived sequence wherever the graph supports it."""
    filled = 0
    for scaffold in plan.scaffolds:
        for i in range(len(scaffold.members) - 1):
            member = scaffold.members[i]
            nxt = scaffold.members[i + 1]
            if member.gap_after <= 0:
                continue
            result = bridge_gap(
                graph,
                (member.segment, member.orientation),
                (nxt.segment, nxt.orientation),
                expected_gap=member.gap_after,
                tolerance=length_tolerance,
                max_length=min(max_bridge_length, max(2000, member.gap_after * 3)),
            )
            if result is None:
                continue
            path, sequence = result
            member.bridge_path = path
            member.bridge_sequence_length = len(sequence)
            member.gap_after = 0
            member.gap_evidence = GAP_ADJACENT if not path else GAP_GRAPH
            filled += 1
    return filled


def _safe(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)
