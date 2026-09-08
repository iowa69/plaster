"""The assembly graph model.

Plastr follows Bandage's double-stranded convention. A *segment* is a single
piece of sequence; it is drawn as one polyline with two ends:

    start ●━━━━━━━━━━━━━━━━━━━━● end
          5' of the + strand    3' of the + strand

A link is stored as ``(from_name, from_orient) -> (to_name, to_orient)``. Which
physical ends that joins follows from the orientations:

    (A,+) -> (B,+)   A.end   -> B.start
    (A,+) -> (B,-)   A.end   -> B.end
    (A,-) -> (B,+)   A.start -> B.start
    (A,-) -> (B,-)   A.start -> B.end

Every link has an equivalent reverse form -- ``(A,+)->(B,+)`` is the same physical
connection as ``(B,-)->(A,-)`` -- so links are canonicalised on insert and stored
once. Traversal exposes both directions.
"""

from __future__ import annotations

import re
from collections import defaultdict, deque
from dataclasses import dataclass, field
from typing import Iterable, Sequence

from .errors import GraphOperationError
from .sequence import flip, gc_content, oriented

# ---------------------------------------------------------------------------
# Segments and links
# ---------------------------------------------------------------------------

# SPAdes/Velvet style names carry coverage: NODE_1_length_5000_cov_12.34
_SPADES_COV = re.compile(r"(?:_cov_|_cov)([0-9]+(?:\.[0-9]+)?)")
_SPADES_LEN = re.compile(r"(?:_length_|_length)([0-9]+)")


@dataclass(slots=True)
class Segment:
    """One node of the assembly graph."""

    name: str
    sequence: str | None = None
    length: int = 0
    depth: float | None = None
    tags: dict[str, object] = field(default_factory=dict)
    # Populated by analysis; kept on the segment so the renderer can read it cheaply.
    ref_hits: list[dict] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.sequence is not None and self.length == 0:
            self.length = len(self.sequence)

    @property
    def has_sequence(self) -> bool:
        return bool(self.sequence)

    @property
    def gc(self) -> float | None:
        return gc_content(self.sequence) if self.sequence else None

    def seq_oriented(self, orientation: str) -> str:
        if self.sequence is None:
            raise GraphOperationError(
                f"segment {self.name!r} has no sequence (graph loaded without sequences)"
            )
        return oriented(self.sequence, orientation)

    def to_dict(self, include_sequence: bool = False) -> dict:
        d: dict[str, object] = {
            "name": self.name,
            "length": self.length,
            "depth": self.depth,
            "gc": self.gc,
        }
        if self.tags:
            d["tags"] = dict(self.tags)
        if self.ref_hits:
            d["ref_hits"] = self.ref_hits
        if include_sequence and self.sequence:
            d["sequence"] = self.sequence
        return d


@dataclass(frozen=True, slots=True)
class Link:
    """A directed connection between two oriented segments."""

    from_name: str
    from_orient: str
    to_name: str
    to_orient: str
    overlap: int = 0
    cigar: str = "*"

    def reverse(self) -> "Link":
        """The same physical connection expressed from the other side."""
        return Link(
            self.to_name,
            flip(self.to_orient),
            self.from_name,
            flip(self.from_orient),
            self.overlap,
            self.cigar,
        )

    def key(self) -> tuple[str, str, str, str]:
        return (self.from_name, self.from_orient, self.to_name, self.to_orient)

    def canonical_key(self) -> tuple[str, str, str, str]:
        """A key that is identical for a link and its reverse form."""
        a = self.key()
        b = self.reverse().key()
        return min(a, b)

    def ends(self) -> tuple[tuple[str, str], tuple[str, str]]:
        """Which physical ends this link joins, as ((seg, 'start'|'end'), ...)."""
        tail = (self.from_name, "end" if self.from_orient == "+" else "start")
        head = (self.to_name, "start" if self.to_orient == "+" else "end")
        return tail, head

    def to_dict(self) -> dict:
        return {
            "from": self.from_name,
            "from_orient": self.from_orient,
            "to": self.to_name,
            "to_orient": self.to_orient,
            "overlap": self.overlap,
        }


@dataclass(slots=True)
class Path:
    """A named ordered walk through the graph (GFA `P` line, SPAdes contigs.paths)."""

    name: str
    steps: list[tuple[str, str]] = field(default_factory=list)  # (segment, orient)
    overlaps: list[int] = field(default_factory=list)
    #: Indices ``i`` such that the join between ``steps[i]`` and ``steps[i+1]``
    #: is a scaffold gap rather than a graph edge. SPAdes writes these as ``;``
    #: in a P line or contigs.paths; the two sides are not adjacent in the graph
    #: and no sequence should be spliced across the join.
    gaps: list[int] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "steps": [f"{n}{o}" for n, o in self.steps],
            "gaps": list(self.gaps),
        }


# ---------------------------------------------------------------------------
# Graph
# ---------------------------------------------------------------------------


class AssemblyGraph:
    """Segments plus links, with the queries the rest of Plastr needs.

    Links are stored canonically (one entry per physical connection) in
    ``self.links``. Adjacency is maintained in both directions so traversal never
    has to scan.
    """

    def __init__(self, name: str = "assembly") -> None:
        self.name = name
        self.segments: dict[str, Segment] = {}
        self.links: dict[tuple[str, str, str, str], Link] = {}
        self.paths: dict[str, Path] = {}
        self.overlap_default: int = 0
        self.source_path: str | None = None
        self.source_format: str | None = None
        # (segment, orient) -> list of (segment, orient) reachable by walking forward
        self._out: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
        # segment -> the canonical keys of every link touching it, so that
        # removing a segment costs its degree rather than a scan of all links.
        self._seg_links: dict[str, set[tuple[str, str, str, str]]] = defaultdict(set)

    # -- construction -------------------------------------------------------

    def add_segment(self, segment: Segment, replace: bool = False) -> Segment:
        if segment.name in self.segments and not replace:
            raise GraphOperationError(f"duplicate segment name {segment.name!r}")
        self.segments[segment.name] = segment
        return segment

    def add_link(self, link: Link, strict: bool = False) -> bool:
        """Add a link. Returns False if it was a duplicate or dangled.

        Dangling links (referring to an unknown segment) are dropped rather than
        raising, because real-world GFA from assemblers occasionally contains
        them; ``strict=True`` turns that into an error instead.
        """
        if link.from_name not in self.segments or link.to_name not in self.segments:
            if strict:
                raise GraphOperationError(
                    f"link references unknown segment: {link.from_name} -> {link.to_name}"
                )
            return False
        key = link.canonical_key()
        if key in self.links:
            return False
        self.links[key] = link
        self._index_link(link)
        return True

    def _index_link(self, link: Link) -> None:
        key = link.canonical_key()
        self._seg_links[link.from_name].add(key)
        self._seg_links[link.to_name].add(key)
        self._out[(link.from_name, link.from_orient)].append((link.to_name, link.to_orient))
        rev = link.reverse()
        # A hairpin (a link joining one end of a segment to itself) is its own
        # reverse; indexing it twice would double its degree.
        if rev.key() != link.key():
            self._out[(rev.from_name, rev.from_orient)].append((rev.to_name, rev.to_orient))

    def _unindex_link(self, link: Link) -> None:
        key = link.canonical_key()
        self._seg_links.get(link.from_name, set()).discard(key)
        self._seg_links.get(link.to_name, set()).discard(key)
        _discard(self._out.get((link.from_name, link.from_orient)), (link.to_name, link.to_orient))
        rev = link.reverse()
        if rev.key() != link.key():
            _discard(
                self._out.get((rev.from_name, rev.from_orient)),
                (rev.to_name, rev.to_orient),
            )

    def _reindex(self) -> None:
        self._out = defaultdict(list)
        self._seg_links = defaultdict(set)
        for link in self.links.values():
            self._index_link(link)

    def add_path(self, path: Path) -> None:
        self.paths[path.name] = path

    # -- basic queries ------------------------------------------------------

    def __len__(self) -> int:
        return len(self.segments)

    def __contains__(self, name: object) -> bool:
        return name in self.segments

    def __getitem__(self, name: str) -> Segment:
        return self.segments[name]

    @property
    def total_length(self) -> int:
        return sum(s.length for s in self.segments.values())

    @property
    def link_count(self) -> int:
        return len(self.links)

    def successors(self, name: str, orient: str = "+") -> list[tuple[str, str]]:
        """Oriented segments reachable by leaving the 3' end of ``(name, orient)``."""
        return list(self._out.get((name, orient), ()))

    def predecessors(self, name: str, orient: str = "+") -> list[tuple[str, str]]:
        """Oriented segments that lead into the 5' end of ``(name, orient)``."""
        return [(n, flip(o)) for n, o in self._out.get((name, flip(orient)), ())]

    def neighbours(self, name: str) -> set[str]:
        out = {n for n, _ in self._out.get((name, "+"), ())}
        out |= {n for n, _ in self._out.get((name, "-"), ())}
        out.discard(name)
        return out

    def degree(self, name: str) -> tuple[int, int]:
        """(links on the start end, links on the end end)."""
        return (
            len(self._out.get((name, "-"), ())),
            len(self._out.get((name, "+"), ())),
        )

    def is_dead_end(self, name: str) -> bool:
        left, right = self.degree(name)
        return left == 0 or right == 0

    def dead_ends(self) -> list[str]:
        """Segments with at least one unconnected end. Isolated segments count twice."""
        return [n for n in self.segments if self.is_dead_end(n)]

    def dead_end_count(self) -> int:
        total = 0
        for name in self.segments:
            left, right = self.degree(name)
            total += (left == 0) + (right == 0)
        return total

    def links_of(self, name: str) -> list[Link]:
        return [
            self.links[key]
            for key in self._seg_links.get(name, ())
            if key in self.links
        ]

    # -- components ---------------------------------------------------------

    def connected_components(self) -> list[list[str]]:
        """Weakly connected components, largest first."""
        seen: set[str] = set()
        components: list[list[str]] = []
        for start in self.segments:
            if start in seen:
                continue
            stack = [start]
            seen.add(start)
            comp = []
            while stack:
                node = stack.pop()
                comp.append(node)
                for nb in self.neighbours(node):
                    if nb not in seen:
                        seen.add(nb)
                        stack.append(nb)
            components.append(comp)
        components.sort(key=lambda c: sum(self.segments[n].length for n in c), reverse=True)
        return components

    def component_map(self) -> dict[str, int]:
        """segment name -> component index (0 is the largest component)."""
        mapping: dict[str, int] = {}
        for idx, comp in enumerate(self.connected_components()):
            for name in comp:
                mapping[name] = idx
        return mapping

    def is_circular(self, name: str) -> bool:
        """True if the segment forms a self-contained loop (a circular contig)."""
        succ = self._out.get((name, "+"), ())
        if len(succ) != 1 or succ[0] != (name, "+"):
            return False
        return len(self._out.get((name, "-"), ())) == 1

    # -- scope (which segments a view should show) --------------------------

    def match_names(
        self, queries: Iterable[str], exact: bool = True
    ) -> tuple[list[str], list[str]]:
        """Resolve typed-in segment names, following Bandage.

        Exact matching tolerates a trailing ``+``/``-`` because Bandage names
        each strand separately and people paste what they see. Partial matching
        takes every segment whose name *contains* the query, which is how one
        pulls a whole ``NODE_12_length_...`` family out with ``NODE_12_``.

        Returns ``(matched names in graph order, queries that matched nothing)``.
        """
        matched: dict[str, None] = {}
        missing: list[str] = []
        for raw in queries:
            query = (raw or "").strip()
            if not query:
                continue
            if exact:
                if query in self.segments:
                    hits = [query]
                elif query[-1] in "+-" and query[:-1] in self.segments:
                    hits = [query[:-1]]
                else:
                    hits = []
            else:
                hits = [name for name in self.segments if query in name]
            if hits:
                matched.update(dict.fromkeys(hits))
            else:
                missing.append(query)
        return list(matched), missing

    def within_distance(
        self,
        seeds: Iterable[str],
        distance: int,
        allowed: set[str] | None = None,
    ) -> set[str]:
        """The seeds plus every segment at most ``distance`` steps from one.

        ``allowed`` walls the walk in: a segment outside it is neither returned
        nor stepped through, so a length or component filter cannot be tunnelled
        under by a two-step expansion.
        """
        frontier = {
            n
            for n in seeds
            if n in self.segments and (allowed is None or n in allowed)
        }
        seen = set(frontier)
        for _ in range(max(0, distance)):
            nxt: set[str] = set()
            for name in frontier:
                for nb in self.neighbours(name):
                    if nb in seen or (allowed is not None and nb not in allowed):
                        continue
                    nxt.add(nb)
            if not nxt:
                break
            seen |= nxt
            frontier = nxt
        return seen

    def in_depth_range(
        self, minimum: float | None = None, maximum: float | None = None
    ) -> list[str]:
        """Segments whose depth lies in ``[minimum, maximum]``; either end may be open.

        A segment of unknown depth is in no range at all -- it cannot be shown to
        be, and letting it through would quietly widen the filter it was meant
        to obey.
        """
        out = []
        for name, seg in self.segments.items():
            depth = seg.depth
            if depth is None:
                continue
            if minimum is not None and depth < minimum:
                continue
            if maximum is not None and depth > maximum:
                continue
            out.append(name)
        return out

    # -- traversal ----------------------------------------------------------

    def walk_sequence(self, steps: Sequence[tuple[str, str]]) -> str:
        """Concatenate a walk, trimming each link's overlap from the following piece."""
        if not steps:
            return ""
        pieces: list[str] = []
        for i, (name, orient) in enumerate(steps):
            seg = self.segments.get(name)
            if seg is None:
                raise GraphOperationError(f"walk references unknown segment {name!r}")
            seq = seg.seq_oriented(orient)
            if i > 0:
                prev_name, prev_orient = steps[i - 1]
                ov = self.overlap_between(prev_name, prev_orient, name, orient)
                if ov:
                    seq = seq[ov:]
            pieces.append(seq)
        return "".join(pieces)

    def overlap_between(
        self, from_name: str, from_orient: str, to_name: str, to_orient: str
    ) -> int:
        link = Link(from_name, from_orient, to_name, to_orient)
        stored = self.links.get(link.canonical_key())
        if stored is not None:
            return stored.overlap
        return self.overlap_default

    def has_link(
        self, from_name: str, from_orient: str, to_name: str, to_orient: str
    ) -> bool:
        return Link(from_name, from_orient, to_name, to_orient).canonical_key() in self.links

    def find_paths(
        self,
        start: tuple[str, str],
        end: tuple[str, str],
        max_length: int,
        max_nodes: int = 30,
        max_paths: int = 24,
    ) -> list[list[tuple[str, str]]]:
        """Bounded search for walks from ``start`` to ``end``.

        Used by graph-aware gap filling: we want the real sequence between two
        scaffolded contigs if the graph provides one of roughly the right length.
        The intermediate walk is returned *without* the start and end segments.
        """
        results: list[list[tuple[str, str]]] = []
        # (current node, path so far, length of path so far)
        queue: deque[tuple[tuple[str, str], list[tuple[str, str]], int]] = deque()
        queue.append((start, [], 0))
        while queue and len(results) < max_paths:
            node, walk, length = queue.popleft()
            for nxt in self.successors(*node):
                if nxt == end:
                    results.append(list(walk))
                    if len(results) >= max_paths:
                        break
                    continue
                if len(walk) >= max_nodes:
                    continue
                seg = self.segments.get(nxt[0])
                if seg is None:
                    continue
                new_len = length + seg.length
                if new_len > max_length:
                    continue
                # Allow a segment to repeat only once, to survive small repeats
                # without exploding on tandem arrays.
                if sum(1 for n, _ in walk if n == nxt[0]) >= 2:
                    continue
                queue.append((nxt, walk + [nxt], new_len))
        return results

    # -- mutation (each returns an inverse for undo) -------------------------

    def remove_segment(self, name: str) -> dict:
        if name not in self.segments:
            raise GraphOperationError(f"no such segment {name!r}")
        seg = self.segments.pop(name)
        removed = self.links_of(name)
        for lk in removed:
            self.links.pop(lk.canonical_key(), None)
            self._unindex_link(lk)
        self._seg_links.pop(name, None)
        self._out.pop((name, "+"), None)
        self._out.pop((name, "-"), None)
        return {
            "op": "add_segment",
            "segment": seg,
            "links": removed,
        }

    def remove_segments(self, names: Iterable[str]) -> dict:
        """Bulk removal. Cheaper than repeated :meth:`remove_segment`."""
        doomed = [n for n in dict.fromkeys(names) if n in self.segments]
        segments = []
        links: list[Link] = []
        seen: set[tuple[str, str, str, str]] = set()
        for name in doomed:
            segments.append(self.segments.pop(name))
            for link in self.links_of(name):
                key = link.canonical_key()
                if key in seen:
                    continue
                seen.add(key)
                links.append(link)
                self.links.pop(key, None)
                self._unindex_link(link)
            self._seg_links.pop(name, None)
            self._out.pop((name, "+"), None)
            self._out.pop((name, "-"), None)
        return {"op": "add_segment", "segments": segments, "links": links}

    def remove_link(self, link: Link) -> dict:
        key = link.canonical_key()
        if key not in self.links:
            raise GraphOperationError("no such link")
        stored = self.links.pop(key)
        self._unindex_link(stored)
        return {"op": "add_link", "links": [stored]}

    def merge_simple_path(self, names: Sequence[str]) -> None:  # pragma: no cover - see simplify
        raise NotImplementedError("use plastr.core.analysis.simplify.merge_unbranching")

    # -- export helpers -----------------------------------------------------

    def segment_names_by_length(self) -> list[str]:
        return sorted(self.segments, key=lambda n: self.segments[n].length, reverse=True)

    def to_dict(self, include_sequence: bool = False) -> dict:
        comp = self.component_map()
        segs = []
        for name, seg in self.segments.items():
            d = seg.to_dict(include_sequence)
            d["component"] = comp.get(name, 0)
            left, right = self.degree(name)
            d["deg_start"] = left
            d["deg_end"] = right
            d["circular"] = self.is_circular(name)
            segs.append(d)
        return {
            "name": self.name,
            "segments": segs,
            "links": [lk.to_dict() for lk in self.links.values()],
            "paths": [p.to_dict() for p in self.paths.values()],
            "source_format": self.source_format,
            "source_path": self.source_path,
        }

    def summary(self) -> dict:
        comps = self.connected_components()
        return {
            "segments": len(self.segments),
            "links": len(self.links),
            "paths": len(self.paths),
            "total_length": self.total_length,
            "components": len(comps),
            "dead_ends": self.dead_end_count(),
            "has_sequences": any(s.has_sequence for s in self.segments.values()),
        }


# ---------------------------------------------------------------------------
# Graph scope
# ---------------------------------------------------------------------------

#: The ways a view can be narrowed, spelled as Bandage's graph scope offers them.
SCOPES = ("entire", "around", "depth", "component")

#: Values for :attr:`Scope.rule`.
NO_TRUNCATION = "none"
BREADTH_FIRST = "breadth-first from the scope seeds"


@dataclass(slots=True)
class Scope:
    """The answer to "which segments should this view contain?"."""

    names: list[str]
    #: How many were in scope before ``max_nodes`` was applied.
    total: int
    seeds: list[str]
    #: Queries that matched no segment. Bandage warns about these by name.
    missing: list[str]
    rule: str = NO_TRUNCATION

    @property
    def dropped(self) -> int:
        return self.total - len(self.names)


def scope_segments(
    graph: AssemblyGraph,
    scope: str = "entire",
    *,
    names: Iterable[str] = (),
    distance: int = 0,
    exact: bool = True,
    min_depth: float | None = None,
    max_depth: float | None = None,
    component: int | None = None,
    min_length: int = 0,
    max_nodes: int | None = None,
    components: dict[str, int] | None = None,
) -> Scope:
    """Pick the segments a drawing should contain, as Bandage's graph scope does.

    ``min_length`` and ``component`` narrow the pool the scope works *within*
    rather than being applied to its result, so "everything two steps from node
    5" never counts a step through a segment the caller already excluded.
    ``distance`` is ignored outside the ``around`` scope, which is what Bandage
    does with its node-distance setting.
    """
    if scope not in SCOPES:
        raise GraphOperationError(
            f"unknown graph scope {scope!r}; choose one of {', '.join(SCOPES)}"
        )
    if component is not None and components is None:
        components = graph.component_map()
    pool = {
        name
        for name, seg in graph.segments.items()
        if seg.length >= min_length
        and (component is None or (components or {}).get(name, 0) == component)
    }

    missing: list[str] = []
    if scope == "around":
        matched, missing = graph.match_names(names, exact)
        seeds = [n for n in matched if n in pool]
        keep = graph.within_distance(seeds, distance, pool)
    elif scope == "depth":
        keep = {n for n in graph.in_depth_range(min_depth, max_depth) if n in pool}
        seeds = []
    else:
        keep = pool
        seeds = []

    total = len(keep)
    if max_nodes is None or total <= max_nodes:
        return Scope(sorted(keep), total, seeds, missing)
    return Scope(
        sorted(_grow(graph, seeds, keep, max_nodes)), total, seeds, missing, BREADTH_FIRST
    )


def _grow(
    graph: AssemblyGraph, seeds: Iterable[str], allowed: set[str], limit: int
) -> list[str]:
    """Breadth-first from the seeds until ``limit`` segments are collected.

    Keeping the *longest* segments instead is easy and useless: the longest
    segments of a large assembly are rarely neighbours, so the drawing becomes a
    field of unconnected sticks. Growing outward keeps whatever is returned
    joined up. With no seeds (the scopes that name none) the largest segment
    anchors the first blob, and each further blob starts from the largest
    segment not yet reached, so the big structures still arrive first.
    """
    def biggest_first(name: str) -> tuple[int, str]:
        return (-graph.segments[name].length, name)

    # Every seed is queued at once. A name the caller typed outranks a neighbour
    # of some other seed, so a tight cap must not spend the whole budget on the
    # first seed's blob and drop a segment that was asked for by name.
    queue: deque[str] = deque(
        sorted((n for n in seeds if n in allowed), key=biggest_first)
    )
    starts = deque(sorted(allowed, key=biggest_first))
    seen: set[str] = set(queue)
    kept: list[str] = []
    while len(kept) < limit:
        if not queue:
            while starts and starts[0] in seen:
                starts.popleft()
            if not starts:
                break
            start = starts.popleft()
            seen.add(start)
            queue.append(start)
        name = queue.popleft()
        kept.append(name)
        for nb in sorted(graph.neighbours(name), key=biggest_first):
            if nb in allowed and nb not in seen:
                seen.add(nb)
                queue.append(nb)
    return kept


def _discard(items: list | None, value) -> None:
    """Remove one occurrence of ``value`` from ``items`` if present."""
    if not items:
        return
    try:
        items.remove(value)
    except ValueError:
        pass


# ---------------------------------------------------------------------------
# Helpers used by parsers
# ---------------------------------------------------------------------------


def depth_from_tags(tags: dict[str, object], length: int) -> float | None:
    """Work out per-base depth from the usual GFA tag spellings.

    ``dp``/``DP`` are depth directly. ``KC``/``RC``/``FC`` are k-mer, read and
    fragment *counts* -- all three become depth once divided by segment length,
    which is what the GFA1 spec defines and what Bandage does.
    """
    for key in ("dp", "DP"):
        if key in tags:
            try:
                return float(tags[key])  # type: ignore[arg-type]
            except (TypeError, ValueError):
                pass
    for key in ("KC", "kc", "RC", "rc", "FC", "fc"):
        if key in tags and length:
            try:
                return float(tags[key]) / length  # type: ignore[arg-type]
            except (TypeError, ValueError, ZeroDivisionError):
                pass
    return None


def depth_from_name(name: str) -> float | None:
    """SPAdes-style names embed coverage: NODE_1_length_5000_cov_12.34."""
    m = _SPADES_COV.search(name)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None


def length_from_name(name: str) -> int | None:
    m = _SPADES_LEN.search(name)
    if m:
        try:
            return int(m.group(1))
        except ValueError:
            return None
    return None
