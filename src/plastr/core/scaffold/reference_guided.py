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
  See :mod:`plastr.core.scaffold.graph_bridge`.
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
    # Scaffold names double as FASTA ids and AGP object names, so two reference
    # sequences that sanitise to the same string must not collide: that wrote
    # duplicate FASTA headers and folded two objects into one non-tiling AGP.
    used_names: set[str] = set()

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
            name=_unique(f"{scaffold_prefix}_{_safe(ref)}", used_names),
            source="reference",
            reference=ref,
        )
        overlap_flags: set[int] = set()
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
                    # The overlap belongs to the contig on the right of the
                    # pair -- it is the one that starts before the previous one
                    # ended. Flagging the left member put the mark one row too
                    # early, including on the first member of a scaffold, which
                    # has no previous member at all.
                    pending_overlap = True
                    member.gap_after = min_gap
                    member.gap_evidence = GAP_DEFAULT
                    # Writing both contigs whole and then a run of Ns emits the
                    # shared bases twice and invents a gap where the reference
                    # says the two placements touch -- a false tandem repeat at
                    # every such junction. Trim instead, but only when the ends
                    # really do match; a reference-estimated overlap is not
                    # evidence enough to delete sequence.
                    shared = _verified_overlap(
                        graph,
                        (member.segment, member.orientation),
                        (nxt.query, "+" if nxt.strand > 0 else "-"),
                        -raw_gap,
                    )
                    if shared:
                        member.trim_next = shared
                        member.gap_after = 0
                        member.gap_evidence = GAP_ADJACENT
                else:
                    pending_overlap = False
                    member.gap_after = max(raw_gap, min_gap)
                    member.gap_evidence = GAP_REFERENCE
            else:
                pending_overlap = False
            scaffold.members.append(member)
            if pending_overlap:
                overlap_flags.add(len(scaffold.members))  # index of the next member

        for index in overlap_flags:
            if index < len(scaffold.members):
                scaffold.members[index].overlaps_previous = True
        overlap_flags.clear()

        plan.scaffolds.append(scaffold)

    # Optionally keep unplaced contigs as single-contig scaffolds so that the
    # exported FASTA is a complete assembly rather than a subset.
    if include_unplaced:
        for name in unplaced:
            plan.scaffolds.append(
                Scaffold(
                    name=_unique(
                        f"{scaffold_prefix}_unplaced_{_safe(name)}", used_names
                    ),
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
        # A segment spliced into a bridge is already written inside a scaffold,
        # so keeping its own pass-through record too writes that sequence
        # twice: 16-26 kb per genome on real data, inflating assembly length
        # and producing duplicate hits downstream.
        consumed = {
            name
            for scaffold in plan.scaffolds
            for member in scaffold.members
            for name, _orient in member.bridge_path
        }
        if consumed:
            before = len(plan.scaffolds)
            plan.scaffolds = [
                s
                for s in plan.scaffolds
                if not (
                    s.source == "unplaced"
                    and len(s.members) == 1
                    and s.members[0].segment in consumed
                )
            ]
            plan.unplaced = [n for n in plan.unplaced if n not in consumed]
            unplaced = plan.unplaced
            dropped = before - len(plan.scaffolds)
            if dropped:
                plan.notes.append(
                    f"{dropped} contig(s) already written inside a graph bridge are "
                    "not repeated as separate records"
                )

    carried = len(plan.scaffolds) - plan.scaffold_count
    note = (
        f"{plan.placed_count} contig(s) placed on {plan.scaffold_count} scaffold(s); "
        f"{len(unplaced)} unplaced, {len(plan.redundant)} redundant"
    )
    if carried:
        note += f" ({carried} unplaced contig(s) carried through as single records)"
    plan.notes.append(note)
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


def _verified_overlap(
    graph: AssemblyGraph,
    left: tuple[str, str],
    right: tuple[str, str],
    estimate: int,
    slack: int = 30,
    minimum: int = 20,
) -> int:
    """Longest real suffix/prefix match between two oriented contig ends.

    The reference says these two placements overlap by roughly ``estimate``
    bases. That is an estimate from coordinates, not a sequence alignment, so
    before deleting anything the bases have to agree exactly. Returns 0 when
    they do not, which leaves an honest gap.
    """
    if estimate < minimum:
        return 0
    a = graph.segments.get(left[0])
    b = graph.segments.get(right[0])
    if a is None or b is None or a.sequence is None or b.sequence is None:
        return 0
    tail = a.seq_oriented(left[1])
    head = b.seq_oriented(right[1])
    hi = min(estimate + slack, len(tail), len(head))
    lo = max(minimum, estimate - slack)
    # Prefer the longest match, so a repeat boundary is cut in the same place a
    # merge would cut it.
    for size in range(hi, lo - 1, -1):
        if tail[-size:] == head[:size]:
            return size
    return 0


def _unique(name: str, used: set[str]) -> str:
    """``name``, suffixed if something already claimed it."""
    candidate = name
    suffix = 2
    while candidate in used:
        candidate = f"{name}_{suffix}"
        suffix += 1
    used.add(candidate)
    return candidate
