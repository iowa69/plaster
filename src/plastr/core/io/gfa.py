"""GFA 1.0 (and pragmatic GFA 2.0) reading and writing.

The reader is deliberately forgiving: assemblers emit slightly different dialects
and a viewer that refuses to open a file is useless. Structural problems raise
``PlastrFormatError`` with a line number; cosmetic ones are ignored.
"""

from __future__ import annotations

import os
import re
from typing import Iterator

from ..errors import PlastrFormatError
from ..model import (
    AssemblyGraph,
    Link,
    Path,
    Segment,
    depth_from_name,
    depth_from_tags,
    length_from_name,
)
from .fasta import open_text

_CIGAR_M = re.compile(r"(\d+)([MIDNSHP=X])")
_TAG = re.compile(r"^([A-Za-z][A-Za-z0-9]):([AifZJHB]):(.*)$")


def _parse_tag(field: str) -> tuple[str, object] | None:
    m = _TAG.match(field)
    if not m:
        return None
    key, typ, value = m.groups()
    if typ == "i":
        try:
            return key, int(value)
        except ValueError:
            return key, value
    if typ == "f":
        try:
            return key, float(value)
        except ValueError:
            return key, value
    return key, value


def _parse_tags(fields: list[str]) -> dict[str, object]:
    tags: dict[str, object] = {}
    for f in fields:
        parsed = _parse_tag(f)
        if parsed:
            tags[parsed[0]] = parsed[1]
    return tags


def cigar_overlap(cigar: str) -> int:
    """Length of the overlap a link CIGAR describes.

    Assemblers write overlaps as e.g. ``55M``. We sum the operations that consume
    both sequences, which is the length shared between the two segments.
    """
    if not cigar or cigar == "*":
        return 0
    total = 0
    for count, op in _CIGAR_M.findall(cigar):
        if op in "M=XD":
            total += int(count)
    return total


def detect_format(path: str | os.PathLike[str]) -> str:
    """Sniff the file type: 'gfa', 'gfa2', 'fastg', 'fasta', or 'unknown'."""
    name = str(path).lower()
    for suffix in (".gz",):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    with open_text(path) as handle:
        for _ in range(200):
            line = handle.readline()
            if not line:
                break
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            head = line.split("\t", 1)[0]
            if head in {"H", "S", "L", "P", "C", "W"} and "\t" in line:
                if head == "H" and "VN:Z:2" in line:
                    return "gfa2"
                if head == "S":
                    parts = line.split("\t")
                    # GFA2 S lines are: S <sid> <slen> <sequence>
                    if len(parts) >= 4 and parts[2].isdigit():
                        return "gfa2"
                return "gfa"
            if head in {"E", "G", "O", "U"} and "\t" in line:
                return "gfa2"
            if line.startswith(">"):
                return "fastg" if name.endswith(".fastg") or ":" in line else "fasta"
    if name.endswith(".fastg"):
        return "fastg"
    if name.endswith((".fa", ".fasta", ".fna", ".ffn", ".contigs")):
        return "fasta"
    if name.endswith((".gfa", ".gfa1")):
        return "gfa"
    return "unknown"


def read_gfa(path: str | os.PathLike[str], strict: bool = False) -> AssemblyGraph:
    """Parse a GFA 1.0 file into an :class:`AssemblyGraph`."""
    graph = AssemblyGraph(name=os.path.basename(str(path)))
    graph.source_path = str(path)
    graph.source_format = "gfa"
    pending_links: list[Link] = []
    pending_paths: list[Path] = []
    overlaps_seen: list[int] = []

    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            rec = fields[0]

            if rec == "H":
                continue

            if rec == "S":
                if len(fields) < 3:
                    raise PlastrFormatError(
                        "S line needs at least a name and a sequence field",
                        str(path),
                        line_no,
                        raw,
                    )
                name = fields[1]
                seq_field = fields[2]
                tags = _parse_tags(fields[3:])
                sequence = None if seq_field == "*" else seq_field.upper()
                length = len(sequence) if sequence is not None else 0
                if sequence is None:
                    for key in ("LN", "ln"):
                        if key in tags:
                            try:
                                length = int(tags[key])  # type: ignore[arg-type]
                            except (TypeError, ValueError):
                                pass
                            break
                    else:
                        guessed = length_from_name(name)
                        if guessed:
                            length = guessed
                depth = depth_from_tags(tags, length) or depth_from_name(name)
                graph.add_segment(
                    Segment(name=name, sequence=sequence, length=length, depth=depth, tags=tags),
                    replace=not strict,
                )
                continue

            if rec == "L":
                if len(fields) < 5:
                    raise PlastrFormatError(
                        "L line needs from, from-orient, to, to-orient",
                        str(path),
                        line_no,
                        raw,
                    )
                f_name, f_or, t_name, t_or = fields[1], fields[2], fields[3], fields[4]
                if f_or not in "+-" or t_or not in "+-":
                    raise PlastrFormatError(
                        f"link orientation must be + or -, got {f_or!r}/{t_or!r}",
                        str(path),
                        line_no,
                        raw,
                    )
                cigar = fields[5] if len(fields) > 5 else "*"
                ov = cigar_overlap(cigar)
                overlaps_seen.append(ov)
                pending_links.append(Link(f_name, f_or, t_name, t_or, ov, cigar))
                continue

            if rec == "P":
                if len(fields) < 3:
                    continue
                pname = fields[1]
                steps = []
                for token in fields[2].split(","):
                    token = token.strip()
                    if len(token) < 2 or token[-1] not in "+-":
                        continue
                    steps.append((token[:-1], token[-1]))
                pending_paths.append(Path(name=pname, steps=steps))
                continue

            if rec == "W":
                # rGFA/GFA1.1 walk line: W sample hap seqid start end walk
                if len(fields) < 7:
                    continue
                walk = fields[6]
                steps = []
                for m in re.finditer(r"([<>])([^<>]+)", walk):
                    direction, seg = m.groups()
                    steps.append((seg, "+" if direction == ">" else "-"))
                if steps:
                    pending_paths.append(Path(name=f"{fields[1]}_{fields[3]}", steps=steps))
                continue

            if rec in {"C", "A"}:
                continue  # containment / read alignment: not needed for the view

            if strict and rec not in {"H", "S", "L", "P", "C", "W", "A"}:
                raise PlastrFormatError(
                    f"unrecognised record type {rec!r}", str(path), line_no, raw
                )

    if not graph.segments:
        raise PlastrFormatError("no S (segment) lines found", str(path))

    for link in pending_links:
        graph.add_link(link, strict=strict)

    for p in pending_paths:
        if all(seg in graph.segments for seg, _ in p.steps):
            graph.add_path(p)

    if overlaps_seen:
        # The modal overlap is a good default for links that omitted a CIGAR.
        counts: dict[int, int] = {}
        for ov in overlaps_seen:
            counts[ov] = counts.get(ov, 0) + 1
        graph.overlap_default = max(counts, key=lambda k: counts[k])

    # Some writers emit '*' for every link even though the segments do overlap.
    # This is a no-op when the CIGARs were present or the joins really are blunt.
    from .overlaps import apply_inferred_overlaps

    apply_inferred_overlaps(graph)
    return graph


def read_gfa2(path: str | os.PathLike[str]) -> AssemblyGraph:
    """Parse GFA 2.0. Only the records that matter for a graph view are read."""
    graph = AssemblyGraph(name=os.path.basename(str(path)))
    graph.source_path = str(path)
    graph.source_format = "gfa2"
    pending: list[Link] = []

    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            rec = fields[0]
            if rec == "S":
                # S <sid> <slen> <sequence> <tag>*
                if len(fields) < 4:
                    raise PlastrFormatError(
                        "GFA2 S line needs sid, slen, sequence", str(path), line_no, raw
                    )
                name = fields[1]
                try:
                    slen = int(fields[2])
                except ValueError as exc:
                    raise PlastrFormatError(
                        "GFA2 S line length is not an integer", str(path), line_no, raw
                    ) from exc
                seq = None if fields[3] == "*" else fields[3].upper()
                tags = _parse_tags(fields[4:])
                depth = depth_from_tags(tags, slen)
                graph.add_segment(
                    Segment(name, seq, slen if seq is None else len(seq), depth, tags),
                    replace=True,
                )
            elif rec == "E":
                # E <eid> <sid1><or> <sid2><or> beg1 end1 beg2 end2 <alignment>
                if len(fields) < 8:
                    continue
                a, b = fields[2], fields[3]
                if not a or not b or a[-1] not in "+-" or b[-1] not in "+-":
                    continue
                try:
                    beg1, end1 = int(fields[4].rstrip("$")), int(fields[5].rstrip("$"))
                except ValueError:
                    beg1 = end1 = 0
                pending.append(Link(a[:-1], a[-1], b[:-1], b[-1], max(0, end1 - beg1)))
            elif rec in {"O", "U"}:
                if len(fields) < 3:
                    continue
                steps = []
                for token in fields[2].split():
                    if token and token[-1] in "+-":
                        steps.append((token[:-1], token[-1]))
                if steps:
                    graph.add_path(Path(fields[1], steps))

    if not graph.segments:
        raise PlastrFormatError("no S (segment) lines found", str(path))
    for link in pending:
        graph.add_link(link)
    return graph


def write_gfa(
    graph: AssemblyGraph,
    path: str | os.PathLike[str],
    include_paths: bool = True,
) -> None:
    """Write the graph back out as GFA 1.0."""
    with open(path, "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for name in graph.segments:
            seg = graph.segments[name]
            seq = seg.sequence if seg.sequence else "*"
            parts = [f"S\t{name}\t{seq}"]
            if seg.sequence is None and seg.length:
                parts.append(f"LN:i:{seg.length}")
            if seg.depth is not None:
                parts.append(f"dp:f:{seg.depth:.6g}")
            for key, value in seg.tags.items():
                if key in {"LN", "ln", "dp", "DP"}:
                    continue
                if isinstance(value, int):
                    parts.append(f"{key}:i:{value}")
                elif isinstance(value, float):
                    parts.append(f"{key}:f:{value:.6g}")
                else:
                    parts.append(f"{key}:Z:{value}")
            out.write("\t".join(parts) + "\n")
        for link in graph.links.values():
            cigar = link.cigar if link.cigar and link.cigar != "*" else (
                f"{link.overlap}M" if link.overlap else "*"
            )
            out.write(
                f"L\t{link.from_name}\t{link.from_orient}\t"
                f"{link.to_name}\t{link.to_orient}\t{cigar}\n"
            )
        if include_paths:
            for p in graph.paths.values():
                if not p.steps:
                    continue
                steps = ",".join(f"{n}{o}" for n, o in p.steps)
                overlaps = "*" if not p.overlaps else ",".join(f"{o}M" for o in p.overlaps)
                out.write(f"P\t{p.name}\t{steps}\t{overlaps}\n")


def gfa_records(graph: AssemblyGraph) -> Iterator[str]:
    """Yield GFA lines without touching the filesystem (used by the download API)."""
    yield "H\tVN:Z:1.0\n"
    for name, seg in graph.segments.items():
        seq = seg.sequence if seg.sequence else "*"
        extra = "" if seg.depth is None else f"\tdp:f:{seg.depth:.6g}"
        yield f"S\t{name}\t{seq}{extra}\n"
    for link in graph.links.values():
        cigar = f"{link.overlap}M" if link.overlap else "*"
        yield (
            f"L\t{link.from_name}\t{link.from_orient}\t"
            f"{link.to_name}\t{link.to_orient}\t{cigar}\n"
        )
