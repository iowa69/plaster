"""Reference-free assembly statistics.

These are the QUAST "basic statistics" numbers, computed from the graph's
segments. NG50/LG50 need an expected genome size; without one they are omitted
rather than guessed.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Sequence

from ..model import AssemblyGraph
from ..sequence import gc_content, n_count


def nx_stat(lengths: Sequence[int], fraction: float, total: int | None = None) -> tuple[int, int]:
    """Return ``(Nx, Lx)`` for the given fraction (0.5 for N50).

    ``Nx`` is the length of the contig at which the cumulative sum first reaches
    ``fraction`` of ``total``; ``Lx`` is how many contigs that took. Passing an
    explicit ``total`` (an expected genome size) gives the NGx/LGx variant.
    """
    if not lengths:
        return 0, 0
    ordered = sorted(lengths, reverse=True)
    target = (total if total is not None else sum(ordered)) * fraction
    running = 0
    for i, length in enumerate(ordered, 1):
        running += length
        if running >= target:
            return length, i
    return 0, 0


def auN(lengths: Sequence[int], total: int | None = None) -> float:
    """Area under the Nx curve -- a single number that is less arbitrary than N50."""
    if not lengths:
        return 0.0
    denom = total if total is not None else sum(lengths)
    if denom <= 0:
        return 0.0
    return sum(length * length for length in lengths) / denom


@dataclass
class AssemblyMetrics:
    """Reference-free assembly statistics."""

    num_contigs: int = 0
    num_contigs_ge_1kb: int = 0
    num_contigs_ge_10kb: int = 0
    num_contigs_ge_50kb: int = 0
    total_length: int = 0
    total_length_ge_1kb: int = 0
    largest_contig: int = 0
    smallest_contig: int = 0
    mean_contig: float = 0.0
    median_contig: int = 0
    gc_percent: float | None = None
    n_per_100kb: float = 0.0
    n50: int = 0
    l50: int = 0
    n75: int = 0
    l75: int = 0
    auN: float = 0.0
    ng50: int | None = None
    lg50: int | None = None
    ng75: int | None = None
    genome_size: int | None = None
    # graph-specific
    num_links: int = 0
    num_components: int = 0
    dead_ends: int = 0
    num_circular: int = 0
    mean_depth: float | None = None
    depth_cv: float | None = None
    # plotting series
    nx_curve: list[tuple[float, int]] = field(default_factory=list)
    cumulative_curve: list[tuple[int, int]] = field(default_factory=list)
    length_histogram: list[tuple[int, int]] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _median(values: Sequence[int]) -> int:
    if not values:
        return 0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) // 2


def _nx_curve(lengths: Sequence[int], total: int | None = None, steps: int = 100):
    """Nx for x = 1..100, for the step plot in the report."""
    if not lengths:
        return []
    ordered = sorted(lengths, reverse=True)
    denom = total if total is not None else sum(ordered)
    if denom <= 0:
        return []
    curve: list[tuple[float, int]] = []
    running = 0
    idx = 0
    for step in range(1, steps + 1):
        target = denom * step / steps
        while running < target and idx < len(ordered):
            running += ordered[idx]
            idx += 1
        curve.append((step, ordered[idx - 1] if idx else 0))
    return curve


def _cumulative(lengths: Sequence[int]) -> list[tuple[int, int]]:
    ordered = sorted(lengths, reverse=True)
    out: list[tuple[int, int]] = []
    running = 0
    for i, length in enumerate(ordered, 1):
        running += length
        out.append((i, running))
    # Down-sample so the JSON stays small for huge assemblies.
    if len(out) > 2000:
        stride = len(out) // 2000 + 1
        sampled = out[::stride]
        if sampled[-1] != out[-1]:
            sampled.append(out[-1])
        return sampled
    return out


def _histogram(lengths: Sequence[int], bins: int = 30) -> list[tuple[int, int]]:
    """Log-spaced length histogram: [(bin_lower_bound, count), ...]."""
    if not lengths:
        return []
    lo = max(1, min(lengths))
    hi = max(lengths)
    if hi <= lo:
        return [(lo, len(lengths))]
    log_lo, log_hi = math.log10(lo), math.log10(hi)
    edges = [10 ** (log_lo + (log_hi - log_lo) * i / bins) for i in range(bins + 1)]
    counts = [0] * bins
    for length in lengths:
        pos = min(
            bins - 1,
            int((math.log10(max(length, 1)) - log_lo) / (log_hi - log_lo) * bins),
        )
        counts[max(0, pos)] += 1
    return [(int(edges[i]), counts[i]) for i in range(bins)]


def compute_metrics(
    graph: AssemblyGraph,
    genome_size: int | None = None,
    min_length: int = 0,
) -> AssemblyMetrics:
    """Compute reference-free statistics over the segments of ``graph``.

    ``min_length`` mirrors QUAST's ``--min-contig`` and defaults to counting
    everything, because in a graph view short segments are real structure rather
    than noise.
    """
    segments = [s for s in graph.segments.values() if s.length >= min_length]
    lengths = [s.length for s in segments]
    m = AssemblyMetrics()
    m.genome_size = genome_size
    m.num_contigs = len(lengths)
    if not lengths:
        return m

    m.total_length = sum(lengths)
    m.num_contigs_ge_1kb = sum(1 for x in lengths if x >= 1000)
    m.num_contigs_ge_10kb = sum(1 for x in lengths if x >= 10_000)
    m.num_contigs_ge_50kb = sum(1 for x in lengths if x >= 50_000)
    m.total_length_ge_1kb = sum(x for x in lengths if x >= 1000)
    m.largest_contig = max(lengths)
    m.smallest_contig = min(lengths)
    m.mean_contig = m.total_length / len(lengths)
    m.median_contig = _median(lengths)

    m.n50, m.l50 = nx_stat(lengths, 0.5)
    m.n75, m.l75 = nx_stat(lengths, 0.75)
    m.auN = auN(lengths)
    if genome_size and genome_size > 0:
        m.ng50, m.lg50 = nx_stat(lengths, 0.5, genome_size)
        m.ng75, _ = nx_stat(lengths, 0.75, genome_size)

    with_seq = [s for s in segments if s.has_sequence]
    if with_seq:
        joined_gc = 0.0
        weighted = 0
        ns = 0
        for s in with_seq:
            assert s.sequence is not None
            joined_gc += gc_content(s.sequence) * s.length
            weighted += s.length
            ns += n_count(s.sequence)
        if weighted:
            m.gc_percent = 100.0 * joined_gc / weighted
            m.n_per_100kb = 100_000.0 * ns / weighted

    depths = [s.depth for s in segments if s.depth is not None]
    if depths:
        mean = sum(depths) / len(depths)
        m.mean_depth = mean
        if mean > 0 and len(depths) > 1:
            var = sum((d - mean) ** 2 for d in depths) / (len(depths) - 1)
            m.depth_cv = math.sqrt(var) / mean

    m.num_links = graph.link_count
    m.num_components = len(graph.connected_components())
    m.dead_ends = graph.dead_end_count()
    m.num_circular = sum(1 for n in graph.segments if graph.is_circular(n))

    m.nx_curve = _nx_curve(lengths)
    m.cumulative_curve = _cumulative(lengths)
    m.length_histogram = _histogram(lengths)
    return m


def compare_metrics(named: dict[str, AssemblyMetrics]) -> dict:
    """Shape several metric sets into the side-by-side table the report uses."""
    rows = [
        ("# contigs", "num_contigs", "{:,}"),
        ("# contigs >= 1 kb", "num_contigs_ge_1kb", "{:,}"),
        ("# contigs >= 10 kb", "num_contigs_ge_10kb", "{:,}"),
        ("Total length", "total_length", "{:,}"),
        ("Largest contig", "largest_contig", "{:,}"),
        ("GC (%)", "gc_percent", "{:.2f}"),
        ("N50", "n50", "{:,}"),
        ("L50", "l50", "{:,}"),
        ("N75", "n75", "{:,}"),
        ("auN", "auN", "{:,.0f}"),
        ("NG50", "ng50", "{:,}"),
        ("LG50", "lg50", "{:,}"),
        ("# N's per 100 kb", "n_per_100kb", "{:.2f}"),
        ("# links", "num_links", "{:,}"),
        ("# components", "num_components", "{:,}"),
        ("# dead ends", "dead_ends", "{:,}"),
    ]
    names = list(named)
    table = []
    for label, key, fmt in rows:
        values = []
        any_value = False
        for name in names:
            raw = getattr(named[name], key, None)
            if raw is None:
                values.append("-")
            else:
                any_value = True
                try:
                    values.append(fmt.format(raw))
                except (ValueError, TypeError):
                    values.append(str(raw))
        if any_value:
            table.append({"label": label, "values": values})
    return {"columns": names, "rows": table}
