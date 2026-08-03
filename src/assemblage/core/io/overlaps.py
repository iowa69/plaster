"""Working out link overlaps when the file does not state them.

GFA records an overlap per link as a CIGAR. FASTG has no such field, yet SPAdes
FASTG edges really do overlap by k-1 bases -- the same overlap the GFA would have
written as ``127M``. Taking those links at face value as blunt joins is not a
cosmetic problem: every merge or scaffold built through such a link duplicates
the shared bases, silently corrupting exported sequence.

So when a graph arrives with no overlap information we measure it, by finding the
longest exact suffix/prefix match across a sample of links and requiring the
answer to be consistent. Assemblers use one k for the whole graph, so a real
overlap shows up as an overwhelming consensus; a graph of genuinely blunt joins
produces no consensus and is left alone.
"""

from __future__ import annotations

from collections import Counter

from ..model import AssemblyGraph, Link

#: No assembler uses a k-mer overlap longer than this.
MAX_OVERLAP = 512

#: Short "overlaps" are coincidence -- two contigs sharing a few bases at a
#: junction by chance. Real k-mer overlaps are far longer.
MIN_CREDIBLE_OVERLAP = 10

#: Fraction of sampled links that must agree before the value is trusted.
CONSENSUS = 0.6


def exact_overlap(left: str, right: str, cap: int = MAX_OVERLAP) -> int:
    """Longest k where the last k bases of ``left`` equal the first k of ``right``."""
    limit = min(cap, len(left), len(right))
    for k in range(limit, 0, -1):
        if left[-k:] == right[:k]:
            return k
    return 0


def infer_overlap(
    graph: AssemblyGraph,
    cap: int = MAX_OVERLAP,
    sample: int = 120,
    min_overlap: int = MIN_CREDIBLE_OVERLAP,
) -> int:
    """The overlap the links of ``graph`` appear to share, or 0 if unclear."""
    links = list(graph.links.values())
    if not links:
        return 0

    step = max(1, len(links) // sample)
    counts: Counter[int] = Counter()
    for link in links[::step][:sample]:
        left = graph.segments.get(link.from_name)
        right = graph.segments.get(link.to_name)
        if left is None or right is None or not left.has_sequence or not right.has_sequence:
            continue
        counts[
            exact_overlap(
                left.seq_oriented(link.from_orient),
                right.seq_oriented(link.to_orient),
                cap,
            )
        ] += 1

    if not counts:
        return 0
    value, hits = counts.most_common(1)[0]
    if value < min_overlap:
        return 0
    if hits < CONSENSUS * sum(counts.values()):
        return 0
    return value


def apply_inferred_overlaps(graph: AssemblyGraph) -> int:
    """Measure and record link overlaps when the file did not supply them.

    Returns the overlap applied, or 0 if none was. Only runs when every link
    claims a blunt join, so an explicit CIGAR is never overridden.
    """
    if not graph.links:
        return 0
    if any(link.overlap for link in graph.links.values()):
        return 0
    if not any(s.has_sequence for s in graph.segments.values()):
        return 0

    overlap = infer_overlap(graph)
    if not overlap:
        return 0

    # The sampled consensus is only a candidate. Apply it link by link, and only
    # where the two sequences really do share those bases -- a graph where some
    # joins overlap and others are blunt would otherwise have its blunt links
    # given a phantom overlap, silently deleting real bases from every merge and
    # scaffold through them.
    applied = 0
    for key, link in list(graph.links.items()):
        left = graph.segments.get(link.from_name)
        right = graph.segments.get(link.to_name)
        if left is None or right is None or not left.has_sequence or not right.has_sequence:
            continue
        a_seq = left.seq_oriented(link.from_orient)
        b_seq = right.seq_oriented(link.to_orient)
        if len(a_seq) < overlap or len(b_seq) < overlap:
            continue
        if a_seq[-overlap:] != b_seq[:overlap]:
            continue
        graph.links[key] = Link(
            link.from_name,
            link.from_orient,
            link.to_name,
            link.to_orient,
            overlap,
            f"{overlap}M",
        )
        applied += 1

    if not applied:
        return 0
    # overlap_default is the fallback used when a walk step has no stored link,
    # so only adopt it when the measured overlap describes most of the graph --
    # not when a single coincidental match was the only thing verifiable.
    if applied >= max(2, 0.5 * len(graph.links)):
        graph.overlap_default = overlap
    return overlap
