"""AGP 2.1 writing and reading.

AGP describes how components (contigs) and gaps are laid out along an object
(scaffold). Plastr always writes the AGP from the same scaffold plan that
produced the FASTA, so the two cannot drift apart.

Columns:
    object  object_beg  object_end  part_number  component_type
    then either  component_id  component_beg  component_end  orientation
    or           gap_length    gap_type        linkage        linkage_evidence
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import IO, Iterable, Iterator

from ..errors import PlastrFormatError

# Component types: W = WGS contig, N = gap of known size, U = gap of unknown size
COMPONENT_CONTIG = "W"
GAP_KNOWN = "N"
GAP_UNKNOWN = "U"


@dataclass(slots=True)
class AgpRow:
    object_name: str
    object_beg: int
    object_end: int
    part_number: int
    component_type: str
    # contig rows
    component_id: str | None = None
    component_beg: int | None = None
    component_end: int | None = None
    orientation: str | None = None
    # gap rows
    gap_length: int | None = None
    gap_type: str = "scaffold"
    linkage: str = "yes"
    linkage_evidence: str = "align_genus"

    def to_line(self) -> str:
        head = [
            self.object_name,
            str(self.object_beg),
            str(self.object_end),
            str(self.part_number),
            self.component_type,
        ]
        if self.component_type in (GAP_KNOWN, GAP_UNKNOWN):
            tail = [
                str(self.gap_length),
                self.gap_type,
                self.linkage,
                self.linkage_evidence,
            ]
        else:
            tail = [
                str(self.component_id),
                str(self.component_beg),
                str(self.component_end),
                str(self.orientation),
            ]
        return "\t".join(head + tail)


def write_agp(
    path: "str | os.PathLike[str] | IO[str]",
    rows: Iterable[AgpRow],
    comments: Iterable[str] = (),
) -> int:
    handle: IO[str]
    should_close = False
    if isinstance(path, (str, os.PathLike)):
        handle = open(path, "w")
        should_close = True
    else:
        handle = path
    count = 0
    try:
        handle.write("##agp-version\t2.1\n")
        for comment in comments:
            for line in str(comment).splitlines():
                handle.write(f"# {line}\n")
        for row in rows:
            handle.write(row.to_line() + "\n")
            count += 1
    finally:
        if should_close:
            handle.close()
    return count


def read_agp(path: str | os.PathLike[str]) -> Iterator[AgpRow]:
    with open(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                raise PlastrFormatError(
                    "AGP row needs 9 columns", str(path), line_no, raw
                )
            ctype = f[4]
            try:
                common = dict(
                    object_name=f[0],
                    object_beg=int(f[1]),
                    object_end=int(f[2]),
                    part_number=int(f[3]),
                    component_type=ctype,
                )
            except ValueError as exc:
                raise PlastrFormatError(
                    "AGP coordinate columns must be integers", str(path), line_no, raw
                ) from exc
            if ctype in (GAP_KNOWN, GAP_UNKNOWN):
                yield AgpRow(
                    **common,
                    gap_length=int(f[5]),
                    gap_type=f[6],
                    linkage=f[7],
                    linkage_evidence=f[8],
                )
            else:
                yield AgpRow(
                    **common,
                    component_id=f[5],
                    component_beg=int(f[6]),
                    component_end=int(f[7]),
                    orientation=f[8],
                )
