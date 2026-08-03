"""Turning a scaffold plan into sequence and the AGP that describes it.

FASTA and AGP are produced in one pass from the same coordinate counter, so the
AGP always describes the FASTA exactly. This is worth insisting on: a mismatched
AGP is worse than no AGP, because downstream tools trust it.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Iterable

from ..errors import GraphOperationError
from ..io.agp import COMPONENT_CONTIG, GAP_KNOWN, GAP_UNKNOWN, AgpRow, write_agp
from ..io.fasta import write_fasta
from ..model import AssemblyGraph
from .plan import GAP_GRAPH, GAP_REFERENCE, ScaffoldPlan


@dataclass
class BuiltScaffolds:
    records: list[tuple[str, str]] = field(default_factory=list)
    agp_rows: list[AgpRow] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    total_length: int = 0
    gap_bases: int = 0
    bridged_bases: int = 0
    num_gaps: int = 0

    @property
    def lengths(self) -> list[int]:
        return [len(seq) for _n, seq in self.records]


def build_scaffolds(
    graph: AssemblyGraph,
    plan: ScaffoldPlan,
    min_gap: int = 100,
) -> BuiltScaffolds:
    """Render the plan into sequences plus matching AGP rows."""
    out = BuiltScaffolds()

    for scaffold in plan.scaffolds:
        pieces: list[str] = []
        rows: list[AgpRow] = []
        position = 1  # AGP is 1-based inclusive
        part = 1

        members = [m for m in scaffold.members if m.segment in graph.segments]
        if len(members) != len(scaffold.members):
            missing = {m.segment for m in scaffold.members} - set(graph.segments)
            out.warnings.append(
                f"scaffold {scaffold.name}: dropped missing segment(s) "
                + ", ".join(sorted(missing))
            )
        if not members:
            continue

        for index, member in enumerate(members):
            segment = graph.segments[member.segment]
            if not segment.has_sequence:
                raise GraphOperationError(
                    f"segment {member.segment!r} has no sequence, so it cannot be "
                    "written to a scaffold. Load a graph that includes sequences."
                )
            seq = segment.seq_oriented(member.orientation)
            pieces.append(seq)
            end = position + len(seq) - 1
            rows.append(
                AgpRow(
                    object_name=scaffold.name,
                    object_beg=position,
                    object_end=end,
                    part_number=part,
                    component_type=COMPONENT_CONTIG,
                    component_id=member.segment,
                    component_beg=1,
                    component_end=len(seq),
                    orientation=member.orientation,
                )
            )
            position = end + 1
            part += 1

            if index == len(members) - 1:
                break

            # Either real bridging sequence from the graph, or a run of Ns.
            if member.bridge_path or member.gap_after == 0:
                bridge = _bridge_sequence_for(graph, members, index, member)
                if bridge:
                    pieces.append(bridge)
                    end = position + len(bridge) - 1
                    rows.append(
                        AgpRow(
                            object_name=scaffold.name,
                            object_beg=position,
                            object_end=end,
                            part_number=part,
                            component_type=COMPONENT_CONTIG,
                            component_id="+".join(
                                f"{n}{o}" for n, o in member.bridge_path
                            )
                            or f"{member.segment}_bridge",
                            component_beg=1,
                            component_end=len(bridge),
                            orientation="+",
                        )
                    )
                    out.bridged_bases += len(bridge)
                    position = end + 1
                    part += 1
                continue

            gap = max(int(member.gap_after), min_gap)
            pieces.append("N" * gap)
            end = position + gap - 1
            known = member.gap_evidence in (GAP_REFERENCE, GAP_GRAPH)
            rows.append(
                AgpRow(
                    object_name=scaffold.name,
                    object_beg=position,
                    object_end=end,
                    part_number=part,
                    component_type=GAP_KNOWN if known else GAP_UNKNOWN,
                    gap_length=gap,
                    gap_type="scaffold",
                    linkage="yes",
                    linkage_evidence=(
                        "align_genus" if member.gap_evidence == GAP_REFERENCE else "unspecified"
                    ),
                )
            )
            out.gap_bases += gap
            out.num_gaps += 1
            position = end + 1
            part += 1

        sequence = "".join(pieces)
        out.records.append((scaffold.name, sequence))
        out.agp_rows.extend(rows)
        out.total_length += len(sequence)

    _verify(out)
    return out


def _bridge_sequence_for(
    graph: AssemblyGraph,
    members: list,
    index: int,
    member,
) -> str:
    """Recompute the bridging sequence so it always matches the current graph."""
    from .graph_bridge import _bridge_sequence

    nxt = members[index + 1]
    try:
        seq = _bridge_sequence(
            graph,
            (member.segment, member.orientation),
            list(member.bridge_path),
            (nxt.segment, nxt.orientation),
        )
    except Exception:
        return ""
    return seq or ""


def _verify(built: BuiltScaffolds) -> None:
    """Check every AGP row against the sequence it claims to describe."""
    by_object: dict[str, list[AgpRow]] = {}
    for row in built.agp_rows:
        by_object.setdefault(row.object_name, []).append(row)
    lengths = {name: len(seq) for name, seq in built.records}
    for name, rows in by_object.items():
        expected = lengths.get(name)
        if expected is None:
            built.warnings.append(f"AGP describes {name} but no sequence was written")
            continue
        last_end = max(r.object_end for r in rows)
        if last_end != expected:
            built.warnings.append(
                f"AGP for {name} ends at {last_end} but the sequence is {expected} bp"
            )
        ordered = sorted(rows, key=lambda r: r.part_number)
        cursor = 1
        for row in ordered:
            if row.object_beg != cursor:
                built.warnings.append(
                    f"AGP for {name} part {row.part_number} starts at "
                    f"{row.object_beg}, expected {cursor}"
                )
                break
            cursor = row.object_end + 1


def write_scaffolds(
    built: BuiltScaffolds,
    fasta_path: str | os.PathLike[str],
    agp_path: str | os.PathLike[str] | None = None,
    wrap: int = 60,
    comments: Iterable[str] = (),
) -> None:
    write_fasta(fasta_path, built.records, wrap=wrap)
    if agp_path:
        write_agp(agp_path, built.agp_rows, comments=comments)
