"""The scaffold plan: an editable description of how contigs become scaffolds.

Both the automatic reference-guided method and the manual drag-to-reorder UI
produce a :class:`ScaffoldPlan`. Everything downstream -- FASTA, AGP, the QC
report -- is generated from the plan, so a hand-edited scaffold and a
computed one cannot behave differently.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

# How a gap between two members was decided
GAP_REFERENCE = "reference"  # size estimated from reference coordinates
GAP_GRAPH = "graph"  # the assembly graph supplied the real sequence
GAP_ADJACENT = "adjacent"  # the two contigs are directly linked in the graph
GAP_MANUAL = "manual"  # the user typed a size
GAP_DEFAULT = "default"  # nothing better was known


@dataclass
class ScaffoldMember:
    """One contig placed in a scaffold."""

    segment: str
    orientation: str = "+"
    #: Bases of N inserted *after* this member. Ignored on the last member.
    gap_after: int = 0
    gap_evidence: str = GAP_DEFAULT
    #: Real sequence recovered from the graph to bridge the gap after this member.
    bridge_path: list[tuple[str, str]] = field(default_factory=list)
    bridge_sequence_length: int = 0
    # Where the reference put this contig, when a reference was used.
    ref: str | None = None
    ref_start: int | None = None
    ref_end: int | None = None
    identity: float | None = None
    #: True when the reference placement overlapped its neighbour rather than
    #: leaving a gap; useful to surface in the UI because it usually means a
    #: collapsed repeat.
    overlaps_previous: bool = False

    def to_dict(self) -> dict:
        d = asdict(self)
        d["bridge_path"] = [f"{n}{o}" for n, o in self.bridge_path]
        return d


@dataclass
class Scaffold:
    name: str
    members: list[ScaffoldMember] = field(default_factory=list)
    source: str = "reference"
    reference: str | None = None

    @property
    def contig_count(self) -> int:
        return len(self.members)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "source": self.source,
            "reference": self.reference,
            "members": [m.to_dict() for m in self.members],
        }


@dataclass
class ScaffoldPlan:
    scaffolds: list[Scaffold] = field(default_factory=list)
    unplaced: list[str] = field(default_factory=list)
    #: Contigs excluded because another contig claimed the same reference span.
    redundant: list[str] = field(default_factory=list)
    method: str = "reference"
    notes: list[str] = field(default_factory=list)

    @property
    def placed_count(self) -> int:
        """Contigs actually anchored to something, not carried through.

        ``include_unplaced`` appends one single-member scaffold per unplaced
        contig so the exported FASTA is a complete assembly. Counting those as
        placed made the summary contradict itself -- "217 placed ... 138
        unplaced" out of a 217-contig assembly.
        """
        return sum(
            len(s.members) for s in self.scaffolds if s.source != "unplaced"
        )

    @property
    def scaffold_count(self) -> int:
        """Scaffolds that join something, excluding pass-through singletons."""
        return sum(1 for s in self.scaffolds if s.source != "unplaced")

    def member_index(self) -> dict[str, tuple[str, int]]:
        """segment -> (scaffold name, position)."""
        out: dict[str, tuple[str, int]] = {}
        for scaffold in self.scaffolds:
            for i, member in enumerate(scaffold.members):
                out[member.segment] = (scaffold.name, i)
        return out

    def to_dict(self) -> dict:
        return {
            "method": self.method,
            "scaffolds": [s.to_dict() for s in self.scaffolds],
            "unplaced": list(self.unplaced),
            "redundant": list(self.redundant),
            "notes": list(self.notes),
            "placed_count": self.placed_count,
            "scaffold_count": self.scaffold_count,
            #: Every record the export will write, including the single-contig
            #: pass-throughs. `scaffold_count` counts only real joins.
            "record_count": len(self.scaffolds),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ScaffoldPlan":
        plan = cls(
            method=data.get("method", "manual"),
            unplaced=list(data.get("unplaced", [])),
            redundant=list(data.get("redundant", [])),
            notes=list(data.get("notes", [])),
        )
        for s in data.get("scaffolds", []):
            scaffold = Scaffold(
                name=s["name"],
                source=s.get("source", "manual"),
                reference=s.get("reference"),
            )
            for m in s.get("members", []):
                bridge = []
                for token in m.get("bridge_path", []) or []:
                    if isinstance(token, str) and len(token) > 1 and token[-1] in "+-":
                        bridge.append((token[:-1], token[-1]))
                    elif isinstance(token, (list, tuple)) and len(token) == 2:
                        bridge.append((token[0], token[1]))
                scaffold.members.append(
                    ScaffoldMember(
                        segment=m["segment"],
                        orientation=m.get("orientation", "+"),
                        gap_after=int(m.get("gap_after", 0) or 0),
                        gap_evidence=m.get("gap_evidence", GAP_MANUAL),
                        bridge_path=bridge,
                        bridge_sequence_length=int(m.get("bridge_sequence_length", 0) or 0),
                        ref=m.get("ref"),
                        ref_start=m.get("ref_start"),
                        ref_end=m.get("ref_end"),
                        identity=m.get("identity"),
                        overlaps_previous=bool(m.get("overlaps_previous", False)),
                    )
                )
            plan.scaffolds.append(scaffold)
        return plan
