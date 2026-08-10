"""Filling scaffold gaps with real sequence taken from the assembly graph.

Reference-guided scaffolding estimates a gap size from reference coordinates and
writes that many Ns. But the assembler often *did* assemble the intervening
sequence -- it just could not decide how to attach it, so it sits in the graph
as one or more unplaced segments. If the graph contains a path between the two
flanking contigs whose length is close to the estimated gap, that path is almost
certainly the missing sequence, and using it turns a gap of Ns into real bases.

The check is deliberately conservative: an ambiguous gap (several plausible
paths of similar length) is left as Ns rather than guessed at.
"""

from __future__ import annotations

from ..model import AssemblyGraph


def bridge_gap(
    graph: AssemblyGraph,
    left: tuple[str, str],
    right: tuple[str, str],
    expected_gap: int,
    tolerance: float = 0.5,
    max_length: int = 50_000,
    max_nodes: int = 20,
) -> tuple[list[tuple[str, str]], str] | None:
    """Find sequence in the graph that bridges ``left`` -> ``right``.

    Returns ``(intermediate_steps, sequence)`` or ``None`` when the graph offers
    no confident answer. An empty step list with an empty sequence means the two
    contigs are directly linked, so the gap closes to nothing.
    """
    if left[0] not in graph.segments or right[0] not in graph.segments:
        return None

    # Direct adjacency: the two contigs are already joined in the graph.
    if graph.has_link(left[0], left[1], right[0], right[1]):
        return [], ""

    if not graph.link_count:
        return None

    lo = max(0, int(expected_gap * (1 - tolerance)))
    hi = int(expected_gap * (1 + tolerance)) + 100
    search_limit = min(max_length, max(hi * 2, 2000))

    candidates = graph.find_paths(
        left, right, max_length=search_limit, max_nodes=max_nodes, max_paths=16
    )
    if not candidates:
        return None

    scored: list[tuple[int, list[tuple[str, str]], str]] = []
    for path in candidates:
        try:
            sequence = _bridge_sequence(graph, left, path, right)
        except Exception:
            continue
        if sequence is None:
            continue
        if not (lo <= len(sequence) <= hi):
            continue
        scored.append((abs(len(sequence) - expected_gap), path, sequence))

    if not scored:
        return None

    scored.sort(key=lambda item: item[0])
    best_delta, best_path, best_sequence = scored[0]

    # Refuse to choose between two similarly good but different answers.
    if len(scored) > 1:
        runner_up = scored[1]
        if runner_up[2] != best_sequence and runner_up[0] <= best_delta * 1.25:
            return None

    return best_path, best_sequence


def bridge_pieces(
    graph: AssemblyGraph,
    left: tuple[str, str],
    middle: list[tuple[str, str]],
) -> list[tuple[str, str, int, int]] | None:
    """Break a bridge into the per-segment slices it is actually made of.

    Returns ``(name, orientation, component_beg, component_end)`` per step, in
    1-based inclusive coordinates *within that segment*, so the AGP can name a
    real segment and a real sub-range instead of a synthesised component id.

    The bridge sequence is the concatenation of the middle segments, each with
    its leading link overlap trimmed, so this decomposition is exact by
    construction -- and checked against ``_bridge_sequence`` in the tests.
    """
    out: list[tuple[str, str, int, int]] = []
    previous = left
    for name, orient in middle:
        segment = graph.segments.get(name)
        if segment is None or not segment.has_sequence:
            return None
        overlap = graph.overlap_between(previous[0], previous[1], name, orient)
        trim = min(overlap, segment.length)
        if trim < segment.length:
            # Trimming the front of the oriented sequence is the back of the
            # contig when the step is reverse-complemented.
            if orient == "-":
                beg, end = 1, segment.length - trim
            else:
                beg, end = trim + 1, segment.length
            out.append((name, orient, beg, end))
        previous = (name, orient)
    return out


def _bridge_sequence(
    graph: AssemblyGraph,
    left: tuple[str, str],
    middle: list[tuple[str, str]],
    right: tuple[str, str],
) -> str | None:
    """Sequence strictly between ``left`` and ``right`` along the walk.

    ``walk_sequence`` trims each link's overlap from the following segment, so
    the bridge is the full walk with the leading left contig and the trailing
    (already overlap-trimmed) right contig removed.
    """
    walk = [left] + middle + [right]
    for name, _ in walk:
        segment = graph.segments.get(name)
        if segment is None or not segment.has_sequence:
            return None

    full = graph.walk_sequence(walk)
    left_len = graph.segments[left[0]].length
    penultimate = middle[-1] if middle else left
    right_overlap = graph.overlap_between(
        penultimate[0], penultimate[1], right[0], right[1]
    )
    right_len = max(0, graph.segments[right[0]].length - right_overlap)

    start = left_len
    end = len(full) - right_len
    if end < start:
        return None
    return full[start:end]


def find_unbranching_paths(
    graph: AssemblyGraph, min_length: int = 2
) -> list[list[tuple[str, str]]]:
    """Maximal chains where each junction has exactly one way in and one way out.

    These are the stretches an assembler could have merged but did not, and
    collapsing them is the graph-only equivalent of scaffolding.
    """
    visited: set[str] = set()
    chains: list[list[tuple[str, str]]] = []

    def single_successor(node: tuple[str, str]) -> tuple[str, str] | None:
        succ = graph.successors(*node)
        if len(succ) != 1:
            return None
        nxt = succ[0]
        # The next node must also have exactly one way back, or this is a join.
        if len(graph.predecessors(*nxt)) != 1:
            return None
        return nxt

    for name in graph.segments:
        if name in visited:
            continue
        for orient in ("+", "-"):
            start = (name, orient)
            # Only begin a chain at a node that is not itself a simple extension.
            preds = graph.predecessors(*start)
            if len(preds) == 1 and len(graph.successors(*preds[0])) == 1:
                continue
            chain = [start]
            seen_local = {name}
            node = start
            while True:
                nxt = single_successor(node)
                if nxt is None or nxt[0] in seen_local:
                    break
                chain.append(nxt)
                seen_local.add(nxt[0])
                node = nxt
            if len(chain) >= min_length:
                chains.append(chain)
                visited.update(n for n, _ in chain)
                break
    return chains
