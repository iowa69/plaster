"""Editing operations on the graph.

Every operation returns a short record describing what changed, which the server
turns into an undo entry and the UI shows in the activity log.
"""

from __future__ import annotations

from typing import Sequence

from ..errors import GraphOperationError
from ..model import AssemblyGraph, Link, Segment
from ..sequence import flip
from .misassembly import Misassembly


def delete_segments(graph: AssemblyGraph, names: Sequence[str]) -> dict:
    """Remove segments and every link touching them."""
    inverse = graph.remove_segments(names)
    return {
        "kind": "delete_segments",
        "segments": inverse["segments"],
        "links": inverse["links"],
        "count": len(inverse["segments"]),
    }


def restore(graph: AssemblyGraph, record: dict) -> None:
    """Undo an operation, given the record it returned."""
    # Reversing is its own inverse, and it edits a segment in place rather than
    # removing anything -- so it needs handling of its own, or undo silently
    # does nothing and the edit is presented as reversible when it is not.
    if record.get("kind") == "reverse_segment":
        name = record.get("name")
        if name in graph.segments:
            reverse_segment(graph, name)
        return

    for segment in record.get("segments", []):
        graph.add_segment(segment, replace=True)
    for link in record.get("links", []):
        graph.add_link(link)
    for name in record.get("added_segments", []):
        if name in graph.segments:
            graph.remove_segment(name)


def filter_graph(
    graph: AssemblyGraph,
    min_length: int = 0,
    min_depth: float | None = None,
    max_depth: float | None = None,
    keep_components: int | None = None,
    min_component_length: int = 0,
) -> dict:
    """Drop segments that fail the given thresholds.

    ``keep_components`` keeps only the N largest connected components, which is
    the usual way to strip the cloud of tiny fragments off a bacterial graph.
    """
    doomed: set[str] = set()
    for name, segment in graph.segments.items():
        if segment.length < min_length:
            doomed.add(name)
        if min_depth is not None and (segment.depth is None or segment.depth < min_depth):
            doomed.add(name)
        if max_depth is not None and segment.depth is not None and segment.depth > max_depth:
            doomed.add(name)

    components = graph.connected_components()
    if keep_components is not None and keep_components > 0:
        for comp in components[keep_components:]:
            doomed.update(comp)
    if min_component_length > 0:
        for comp in components:
            if sum(graph.segments[n].length for n in comp) < min_component_length:
                doomed.update(comp)

    return delete_segments(graph, sorted(doomed))


def split_segment(graph: AssemblyGraph, name: str, positions: Sequence[int]) -> dict:
    """Cut a segment at the given 0-based offsets.

    Links on the original segment's start end move to the first piece and links
    on its end end move to the last piece; consecutive pieces are joined.
    """
    segment = graph.segments.get(name)
    if segment is None:
        raise GraphOperationError(f"no such segment {name!r}")
    if not segment.has_sequence:
        raise GraphOperationError(
            f"segment {name!r} has no sequence, so it cannot be split"
        )
    cuts = sorted({p for p in positions if 0 < p < segment.length})
    if not cuts:
        raise GraphOperationError(
            f"no usable split position inside {name!r} (length {segment.length})"
        )

    bounds = [0, *cuts, segment.length]
    assert segment.sequence is not None
    pieces: list[Segment] = []
    for i in range(len(bounds) - 1):
        start, end = bounds[i], bounds[i + 1]
        piece_name = _unique_name(graph, f"{name}_part{i + 1}")
        pieces.append(
            Segment(
                name=piece_name,
                sequence=segment.sequence[start:end],
                length=end - start,
                depth=segment.depth,
                tags=dict(segment.tags),
            )
        )

    old_links = graph.links_of(name)
    inverse = graph.remove_segment(name)

    for piece in pieces:
        graph.add_segment(piece, replace=True)
    for left, right in zip(pieces, pieces[1:]):
        graph.add_link(Link(left.name, "+", right.name, "+", 0, "*"))

    # Reattach the original external links to the correct end piece.
    first, last = pieces[0].name, pieces[-1].name
    for link in old_links:
        if link.from_name == name and link.to_name == name:
            continue  # a self-loop on a split segment is not meaningful
        if link.from_name == name:
            anchor = last if link.from_orient == "+" else first
            new_from_orient = link.from_orient
            graph.add_link(
                Link(anchor, new_from_orient, link.to_name, link.to_orient, link.overlap, link.cigar)
            )
        else:
            anchor = first if link.to_orient == "+" else last
            graph.add_link(
                Link(link.from_name, link.from_orient, anchor, link.to_orient, link.overlap, link.cigar)
            )

    return {
        "kind": "split_segment",
        "segments": [inverse["segment"]],
        "links": inverse["links"],
        "added_segments": [p.name for p in pieces],
        "count": len(pieces),
    }


def break_misassemblies(
    graph: AssemblyGraph,
    events: Sequence[Misassembly],
    extensive_only: bool = True,
    min_piece: int = 200,
) -> dict:
    """Split contigs at detected misassembly breakpoints.

    This is the step that makes a reference-based evaluation actionable: instead
    of only reporting that a contig is chimeric, cut it there so the pieces can
    be scaffolded into their correct places.
    """
    by_contig: dict[str, list[int]] = {}
    for event in events:
        if extensive_only and not event.is_extensive:
            continue
        by_contig.setdefault(event.contig, []).append(event.contig_pos)

    all_segments: list[Segment] = []
    all_links: list[Link] = []
    added: list[str] = []
    broken = 0
    for contig, positions in by_contig.items():
        segment = graph.segments.get(contig)
        if segment is None or not segment.has_sequence:
            continue
        usable = [
            p
            for p in sorted(set(positions))
            if min_piece <= p <= segment.length - min_piece
        ]
        if not usable:
            continue
        record = split_segment(graph, contig, usable)
        all_segments.extend(record["segments"])
        all_links.extend(record["links"])
        added.extend(record["added_segments"])
        broken += 1

    return {
        "kind": "break_misassemblies",
        "segments": all_segments,
        "links": all_links,
        "added_segments": added,
        "count": broken,
        "pieces": len(added),
    }


def merge_path(graph: AssemblyGraph, steps: Sequence[tuple[str, str]], new_name: str | None = None) -> dict:
    """Collapse an unbranching walk into a single segment."""
    if len(steps) < 2:
        raise GraphOperationError("a merge needs at least two segments")
    names = [n for n, _ in steps]
    missing = [n for n in names if n not in graph.segments]
    if missing:
        raise GraphOperationError(f"unknown segment(s): {', '.join(missing)}")

    # Every consecutive pair must actually be linked. walk_sequence falls back
    # to the graph's default overlap when it finds no link, so merging two
    # unconnected contigs would silently delete k-1 real bases from the second.
    for (a_name, a_or), (b_name, b_or) in zip(steps, steps[1:]):
        if not graph.has_link(a_name, a_or, b_name, b_or):
            raise GraphOperationError(
                f"{a_name}{a_or} and {b_name}{b_or} are not connected in the graph, "
                "so they cannot be merged. Use scaffolding to join contigs that the "
                "assembly does not link."
            )

    sequence = graph.walk_sequence(steps)
    lengths = [graph.segments[n].length for n in names]
    depths = [graph.segments[n].depth for n in names]
    weighted = [
        (d * l) for d, l in zip(depths, lengths) if d is not None
    ]
    total_len = sum(l for d, l in zip(depths, lengths) if d is not None)
    depth = (sum(weighted) / total_len) if total_len else None

    merged_name = _unique_name(graph, new_name or f"merged_{names[0]}")

    # Remember what the outside world was attached to before we delete anything.
    head, head_orient = steps[0]
    tail, tail_orient = steps[-1]
    # Keep each boundary link's overlap: rebuilding them with the Link defaults
    # would declare blunt joins where the assembler declared an overlap, and the
    # next merge or scaffold through that link would duplicate k-1 bases.
    inside = set(names)
    incoming = [
        (n, o, graph.overlap_between(n, o, head, head_orient))
        for n, o in graph.predecessors(head, head_orient)
        if n not in inside
    ]
    outgoing = [
        (n, o, graph.overlap_between(tail, tail_orient, n, o))
        for n, o in graph.successors(tail, tail_orient)
        if n not in inside
    ]

    inverse = graph.remove_segments(names)
    removed_segments: list[Segment] = inverse["segments"]
    removed_links: list[Link] = inverse["links"]

    graph.add_segment(
        Segment(merged_name, sequence, len(sequence), depth), replace=True
    )
    for name, orient, overlap in incoming:
        graph.add_link(
            Link(name, orient, merged_name, "+", overlap, f"{overlap}M" if overlap else "*")
        )
    for name, orient, overlap in outgoing:
        graph.add_link(
            Link(merged_name, "+", name, orient, overlap, f"{overlap}M" if overlap else "*")
        )

    return {
        "kind": "merge_path",
        "segments": removed_segments,
        "links": removed_links,
        "added_segments": [merged_name],
        "count": len(names),
        "new_name": merged_name,
    }


def simplify(graph: AssemblyGraph, min_chain: int = 2) -> dict:
    """Merge every unbranching chain in the graph.

    The graph-only analogue of scaffolding: wherever the assembler left a chain
    of segments with no choice at either end, join them.
    """
    from ..scaffold.graph_bridge import find_unbranching_paths

    chains = find_unbranching_paths(graph, min_length=min_chain)
    all_segments: list[Segment] = []
    all_links: list[Link] = []
    added: list[str] = []
    merged = 0
    for chain in chains:
        if any(n not in graph.segments for n, _ in chain):
            continue
        try:
            record = merge_path(graph, chain)
        except GraphOperationError:
            continue
        all_segments.extend(record["segments"])
        all_links.extend(record["links"])
        added.extend(record["added_segments"])
        merged += 1
    return {
        "kind": "simplify",
        "segments": all_segments,
        "links": all_links,
        "added_segments": added,
        "count": merged,
    }


def reverse_segment(graph: AssemblyGraph, name: str) -> dict:
    """Flip a segment's stored strand, rewriting its links to match."""
    from ..sequence import revcomp

    segment = graph.segments.get(name)
    if segment is None:
        raise GraphOperationError(f"no such segment {name!r}")
    if not segment.has_sequence:
        raise GraphOperationError(f"segment {name!r} has no sequence to reverse")
    assert segment.sequence is not None
    segment.sequence = revcomp(segment.sequence)

    rewritten: list[Link] = []
    for link in list(graph.links.values()):
        if name not in (link.from_name, link.to_name):
            continue
        graph.links.pop(link.canonical_key(), None)
        new = Link(
            link.from_name,
            flip(link.from_orient) if link.from_name == name else link.from_orient,
            link.to_name,
            flip(link.to_orient) if link.to_name == name else link.to_orient,
            link.overlap,
            link.cigar,
        )
        rewritten.append(new)
    for link in rewritten:
        graph.links[link.canonical_key()] = link
    graph._reindex()
    return {"kind": "reverse_segment", "name": name, "count": 1}


def _unique_name(graph: AssemblyGraph, base: str) -> str:
    if base not in graph.segments:
        return base
    i = 2
    while f"{base}_{i}" in graph.segments:
        i += 1
    return f"{base}_{i}"
