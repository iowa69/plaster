"""Standalone HTML QC report.

Everything the report needs is baked into one file: the CSS lives in a
``<style>`` block, the charts are inline SVG generated here in Python, and the
only JavaScript is a few dozen lines of vanilla code for tab switching and
table sorting. There are no CDN links and no chart library, so the report can
be emailed, archived next to the assembly, or opened on a cluster login node
with no network.

The public API is deliberately small::

    html = build_report(metrics, reference_report=..., plan=..., built=...)
    path = write_report("qc.html", metrics=metrics)

Every argument except ``metrics`` is optional and ``None`` is handled
everywhere -- a missing value prints as an en dash rather than raising.
"""

from __future__ import annotations

import html
import math
import os
from datetime import datetime, timezone
from typing import Any, Mapping, Sequence

from .analysis.metrics import AssemblyMetrics, compare_metrics

__all__ = ["build_report", "write_report", "render_charts"]

DASH = "–"

TEMPLATE_NAME = "report.html.j2"

# Colours used both by the CSS and by the generated SVG, kept here so the two
# never drift apart.
INK = "#1b1f23"
MUTED = "#6a737d"
GRID = "#e3e6ea"
AXIS = "#aab1b8"
ACCENT = "#2563a8"
ACCENT_FILL = "#cbdcf0"
COVERED = "#3d8f5a"

KIND_COLOURS = {
    "relocation": "#d97706",
    "inversion": "#b3357a",
    "translocation": "#c02626",
    "local": "#6a737d",
}


# --------------------------------------------------------------------------
# number formatting
# --------------------------------------------------------------------------


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return True
    return False


def fmt_int(value: Any) -> str:
    """``12345`` -> ``12,345``; ``None`` -> en dash."""
    if _is_missing(value):
        return DASH
    try:
        return f"{int(round(float(value))):,}"
    except (TypeError, ValueError):
        return str(value)


def fmt_float(value: Any, dp: int = 2) -> str:
    if _is_missing(value):
        return DASH
    try:
        return f"{float(value):,.{dp}f}"
    except (TypeError, ValueError):
        return str(value)


def fmt_pct(value: Any, dp: int = 2) -> str:
    if _is_missing(value):
        return DASH
    try:
        return f"{float(value):.{dp}f}%"
    except (TypeError, ValueError):
        return str(value)


def fmt_bp(value: Any) -> str:
    """Human-scaled base counts: ``1,240,000`` -> ``1.24 Mb``."""
    if _is_missing(value):
        return DASH
    try:
        n = float(value)
    except (TypeError, ValueError):
        return str(value)
    sign = "-" if n < 0 else ""
    n = abs(n)
    if n >= 1_000_000_000:
        return f"{sign}{n / 1_000_000_000:.2f} Gb"
    if n >= 1_000_000:
        return f"{sign}{n / 1_000_000:.2f} Mb"
    if n >= 10_000:
        return f"{sign}{n / 1_000:.1f} kb"
    if n >= 1_000:
        return f"{sign}{n / 1_000:.2f} kb"
    return f"{sign}{int(round(n)):,} bp"


def fmt_span(value: Any) -> str:
    """Exact count plus a human-scaled hint, e.g. ``169,000 (169.0 kb)``."""
    if _is_missing(value):
        return DASH
    try:
        n = int(round(float(value)))
    except (TypeError, ValueError):
        return str(value)
    if abs(n) < 1000:
        return f"{n:,}"
    return f"{n:,} ({fmt_bp(n)})"


def _axis_label(value: float) -> str:
    """Compact tick label: 1.5M, 24k, 900."""
    n = abs(value)
    if n >= 1_000_000_000:
        text = f"{value / 1_000_000_000:g}G"
    elif n >= 1_000_000:
        text = f"{value / 1_000_000:g}M"
    elif n >= 1_000:
        text = f"{value / 1_000:g}k"
    else:
        text = f"{value:g}"
    return text


def _e(text: Any) -> str:
    return html.escape(str(text), quote=True)


def _round(value: float, places: int = 2) -> float:
    """Coordinate rounding that can never emit ``nan`` or ``inf`` into the SVG."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return 0.0
    if math.isnan(v) or math.isinf(v):
        return 0.0
    return round(v, places)


# --------------------------------------------------------------------------
# tiny SVG plotting kit
# --------------------------------------------------------------------------


class _Plot:
    """A minimal cartesian plotting surface that emits SVG source."""

    def __init__(
        self,
        width: int = 720,
        height: int = 300,
        margin: tuple[int, int, int, int] = (18, 18, 46, 74),
        title: str = "",
    ) -> None:
        self.width = width
        self.height = height
        self.top, self.right, self.bottom, self.left = margin
        self.title = title
        self.parts: list[str] = []
        self.x0, self.x1 = 0.0, 1.0
        self.y0, self.y1 = 0.0, 1.0
        self.log_x = False

    # -- geometry ---------------------------------------------------------
    @property
    def plot_w(self) -> float:
        return max(1.0, self.width - self.left - self.right)

    @property
    def plot_h(self) -> float:
        return max(1.0, self.height - self.top - self.bottom)

    def set_domain(
        self,
        x0: float,
        x1: float,
        y0: float,
        y1: float,
        log_x: bool = False,
    ) -> None:
        self.log_x = log_x
        if log_x:
            x0 = max(x0, 1e-9)
            x1 = max(x1, x0 * 10)
            self.x0, self.x1 = math.log10(x0), math.log10(x1)
        else:
            if x1 <= x0:
                x1 = x0 + 1
            self.x0, self.x1 = float(x0), float(x1)
        if y1 <= y0:
            y1 = y0 + 1
        self.y0, self.y1 = float(y0), float(y1)

    def sx(self, value: float) -> float:
        if self.log_x:
            value = math.log10(max(float(value), 1e-9))
        frac = (float(value) - self.x0) / (self.x1 - self.x0)
        return _round(self.left + frac * self.plot_w)

    def sy(self, value: float) -> float:
        frac = (float(value) - self.y0) / (self.y1 - self.y0)
        return _round(self.top + (1.0 - frac) * self.plot_h)

    # -- decoration -------------------------------------------------------
    def frame(
        self,
        x_ticks: Sequence[float],
        y_ticks: Sequence[float],
        x_label: str = "",
        y_label: str = "",
        x_fmt=None,
        y_fmt=None,
    ) -> None:
        x_fmt = x_fmt or _axis_label
        y_fmt = y_fmt or _axis_label
        bottom = self.sy(self.y0)
        left = self.sx(10**self.x0 if self.log_x else self.x0)

        for tick in y_ticks:
            y = self.sy(tick)
            self.parts.append(
                f'<line class="grid" x1="{self.left}" y1="{y}" '
                f'x2="{_round(self.left + self.plot_w)}" y2="{y}"/>'
            )
            self.parts.append(
                f'<text class="tick" x="{self.left - 8}" y="{_round(y + 3.5)}" '
                f'text-anchor="end">{_e(y_fmt(tick))}</text>'
            )
        for tick in x_ticks:
            x = self.sx(tick)
            self.parts.append(
                f'<line class="grid" x1="{x}" y1="{self.top}" '
                f'x2="{x}" y2="{bottom}"/>'
            )
            self.parts.append(
                f'<text class="tick" x="{x}" y="{_round(bottom + 17)}" '
                f'text-anchor="middle">{_e(x_fmt(tick))}</text>'
            )

        self.parts.append(
            f'<line class="axis" x1="{self.left}" y1="{bottom}" '
            f'x2="{_round(self.left + self.plot_w)}" y2="{bottom}"/>'
        )
        self.parts.append(
            f'<line class="axis" x1="{left}" y1="{self.top}" '
            f'x2="{left}" y2="{bottom}"/>'
        )
        if x_label:
            self.parts.append(
                f'<text class="axlabel" x="{_round(self.left + self.plot_w / 2)}" '
                f'y="{self.height - 6}" text-anchor="middle">{_e(x_label)}</text>'
            )
        if y_label:
            cy = _round(self.top + self.plot_h / 2)
            self.parts.append(
                f'<text class="axlabel" x="14" y="{cy}" text-anchor="middle" '
                f'transform="rotate(-90 14 {cy})">{_e(y_label)}</text>'
            )

    def add(self, markup: str) -> None:
        self.parts.append(markup)

    def render(self, css_class: str = "chart") -> str:
        body = "".join(self.parts)
        # ``max-width`` pins the chart to its natural size so the text inside it
        # keeps a consistent size instead of scaling with the column width.
        return (
            f'<svg class="{css_class}" viewBox="0 0 {self.width} {self.height}" '
            f'width="100%" style="max-width:{self.width}px" '
            f'preserveAspectRatio="xMidYMid meet" '
            f'role="img" aria-label="{_e(self.title)}" '
            f'xmlns="http://www.w3.org/2000/svg">{body}</svg>'
        )


def _nice_ticks(
    lo: float,
    hi: float,
    count: int = 5,
    integer: bool = False,
) -> list[float]:
    """A handful of round numbers spanning ``[lo, hi]``.

    With ``integer=True`` the step is never smaller than 1, so a count axis
    cannot end up labelled 0, 0, 0, 1 after rounding.
    """
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return [lo] if math.isfinite(lo) else [0.0]
    raw = (hi - lo) / max(1, count)
    magnitude = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1.0
    step = magnitude
    for mult in (1, 2, 2.5, 5, 10):
        step = magnitude * mult
        if raw <= step:
            break
    if integer:
        step = max(1.0, round(step))
    start = math.ceil(lo / step) * step
    ticks: list[float] = []
    value = start
    while value <= hi + step * 1e-9 and len(ticks) < 40:
        ticks.append(round(value, 10))
        value += step
    if not ticks:
        ticks = [lo, hi]
    return ticks


def _log_ticks(lo: float, hi: float) -> list[float]:
    """Decade ticks, subdivided at 2x/5x when the range spans few decades."""
    lo = max(lo, 1.0)
    hi = max(hi, lo * 10)
    decades = math.log10(hi) - math.log10(lo)
    multipliers = (1,) if decades > 4 else ((1, 5) if decades > 2 else (1, 2, 5))
    ticks: list[float] = []
    exp = math.floor(math.log10(lo))
    while 10**exp <= hi * 1.0000001 and len(ticks) < 24:
        for mult in multipliers:
            value = float(10**exp) * mult
            if lo * 0.999 <= value <= hi * 1.0000001:
                ticks.append(value)
        exp += 1
    return ticks or [lo, hi]


def _empty_chart(message: str, height: int = 160, width: int = 500) -> str:
    return (
        f'<svg class="chart empty" viewBox="0 0 {width} {height}" width="100%" '
        f'style="max-width:{width}px" role="img" '
        f'aria-label="{_e(message)}" xmlns="http://www.w3.org/2000/svg">'
        f'<rect x="0.5" y="0.5" width="{width - 1}" height="{height - 1}" fill="none" '
        f'stroke="{GRID}" stroke-dasharray="4 4"/>'
        f'<text x="{width // 2}" y="{height // 2}" text-anchor="middle" class="tick">'
        f"{_e(message)}</text></svg>"
    )


# --------------------------------------------------------------------------
# the individual charts
# --------------------------------------------------------------------------


def svg_cumulative(curve: Sequence[Sequence[int]] | None) -> str:
    """Cumulative assembly length against contig rank."""
    points = [(float(a), float(b)) for a, b in (curve or [])]
    if len(points) < 1:
        return _empty_chart("no contigs to plot")
    if len(points) == 1:
        points = [(0.0, 0.0), points[0]]

    plot = _Plot(width=500, height=300, title="Cumulative length")
    max_x = max(p[0] for p in points)
    max_y = max(p[1] for p in points)
    plot.set_domain(0, max_x, 0, max_y * 1.05)
    plot.frame(
        _nice_ticks(0, max_x, 5, integer=True),
        _nice_ticks(0, max_y * 1.05, 5),
        x_label="Contig rank (longest first)",
        y_label="Cumulative length (bp)",
        x_fmt=lambda v: f"{int(v):,}" if max_x < 1000 else _axis_label(v),
    )

    line = " ".join(f"{plot.sx(x)},{plot.sy(y)}" for x, y in points)
    area = (
        f"{plot.sx(points[0][0])},{plot.sy(0)} "
        + line
        + f" {plot.sx(points[-1][0])},{plot.sy(0)}"
    )
    plot.add(f'<polygon class="area" points="{area}"/>')
    plot.add(f'<polyline class="series" points="{line}"/>')

    # Mark the total.
    tx, ty = plot.sx(points[-1][0]), plot.sy(points[-1][1])
    plot.add(f'<circle class="marker" cx="{tx}" cy="{ty}" r="3"/>')
    anchor = "end" if tx > plot.left + plot.plot_w * 0.6 else "start"
    dx = -8 if anchor == "end" else 8
    plot.add(
        f'<text class="note" x="{_round(tx + dx)}" y="{_round(ty - 8)}" '
        f'text-anchor="{anchor}">total {_e(fmt_bp(points[-1][1]))}</text>'
    )
    return plot.render()


def svg_nx(curve: Sequence[Sequence[float]] | None, n50: int | None = None) -> str:
    """The Nx step curve, x running 1..100."""
    points = [(float(a), float(b)) for a, b in (curve or [])]
    if not points:
        return _empty_chart("no contigs to plot")

    plot = _Plot(width=500, height=300, title="Nx curve")
    max_y = max(p[1] for p in points) or 1.0
    plot.set_domain(0, 100, 0, max_y * 1.05)
    plot.frame(
        [0, 20, 40, 60, 80, 100],
        _nice_ticks(0, max_y * 1.05, 5),
        x_label="x (% of total assembly length)",
        y_label="Nx contig length (bp)",
        x_fmt=lambda v: f"{int(v)}",
    )

    # Step plot: hold each value across its bin.
    step: list[str] = []
    previous_x = 0.0
    for x, y in points:
        step.append(f"{plot.sx(previous_x)},{plot.sy(y)}")
        step.append(f"{plot.sx(x)},{plot.sy(y)}")
        previous_x = x
    line = " ".join(step)
    plot.add(
        f'<polygon class="area" points="{plot.sx(0)},{plot.sy(0)} {line} '
        f'{plot.sx(previous_x)},{plot.sy(0)}"/>'
    )
    plot.add(f'<polyline class="series" points="{line}"/>')

    if n50:
        y = plot.sy(float(n50))
        x = plot.sx(50)
        plot.add(
            f'<line class="guide" x1="{plot.left}" y1="{y}" '
            f'x2="{x}" y2="{y}"/>'
        )
        plot.add(
            f'<line class="guide" x1="{x}" y1="{y}" x2="{x}" '
            f'y2="{plot.sy(plot.y0)}"/>'
        )
        plot.add(
            f'<text class="note" x="{_round(x + 8)}" y="{_round(y - 7)}">'
            f"N50 {_e(fmt_bp(n50))}</text>"
        )
    return plot.render()


def svg_histogram(
    histogram: Sequence[Sequence[int]] | None, largest: int | None = None
) -> str:
    """Contig length histogram on a log x axis.

    ``largest`` is the longest contig. Without it the axis maximum has to be
    reverse-engineered from the bin lower bounds, which the metrics emit as
    truncated integers -- so the recovered ratio came out slightly small and
    the last bin claimed an upper bound below the contig it was counting.
    """
    bins = [(float(a), int(b)) for a, b in (histogram or [])]
    bins = [(max(lo, 1.0), count) for lo, count in bins]
    if not bins:
        return _empty_chart("no contigs to plot", 200, 1000)

    lo = bins[0][0]
    if len(bins) > 1:
        ratio = bins[1][0] / bins[0][0] if bins[0][0] > 0 else 2.0
        ratio = ratio if ratio > 1.0 else 2.0
        hi = bins[-1][0] * ratio
    else:
        ratio = 10.0
        hi = lo * 10
    if largest:
        hi = max(hi, float(largest))

    plot = _Plot(
        width=1000, height=320, margin=(18, 22, 48, 78),
        title="Contig length histogram",
    )
    max_count = max(c for _lo, c in bins) or 1
    plot.set_domain(lo, hi, 0, max_count * 1.12, log_x=True)
    plot.frame(
        _log_ticks(lo, hi),
        _nice_ticks(0, max_count * 1.12, 5, integer=True),
        x_label="Contig length (bp, log scale)",
        y_label="Number of contigs",
        x_fmt=fmt_bp,
        y_fmt=lambda v: f"{int(v):,}",
    )

    edges = [b[0] for b in bins] + [hi]
    baseline = plot.sy(0)
    for i, (edge, count) in enumerate(bins):
        if count <= 0:
            continue
        x_left = plot.sx(edge)
        x_right = plot.sx(edges[i + 1])
        width = max(1.0, x_right - x_left - 1.5)
        y = plot.sy(count)
        height = max(1.0, baseline - y)
        plot.add(
            f'<rect class="bar" x="{_round(x_left + 0.75)}" y="{y}" '
            f'width="{_round(width)}" height="{_round(height)}">'
            f"<title>{_e(fmt_int(count))} contigs between "
            f"{_e(fmt_bp(edge))} and {_e(fmt_bp(edges[i + 1]))}</title></rect>"
        )
    return plot.render()


def svg_ideogram(reference_report: Any) -> str:
    """Reference coverage bars with misassembly breakpoints marked.

    The reference sequence lengths are not stored on the report, so they are
    recovered from ``per_reference_coverage`` (a percentage) and the covered
    span; where that is impossible the furthest block end is used instead.
    """
    if reference_report is None:
        return ""
    blocks: Mapping[str, Sequence[Sequence[int]]] = (
        getattr(reference_report, "coverage_blocks", None) or {}
    )
    coverage: Mapping[str, float] = (
        getattr(reference_report, "per_reference_coverage", None) or {}
    )
    names = sorted(set(blocks) | set(coverage))
    if not names:
        return _empty_chart("no alignments to the reference", 120, 1000)

    lengths: dict[str, float] = {}
    for name in names:
        spans = [(int(a), int(b)) for a, b in (blocks.get(name) or [])]
        covered = sum(b - a for a, b in spans)
        pct = float(coverage.get(name) or 0.0)
        furthest = max((b for _a, b in spans), default=0)
        if pct > 0 and covered > 0:
            estimate = 100.0 * covered / pct
        else:
            estimate = float(furthest)
        lengths[name] = max(estimate, float(furthest), 1.0)

    longest = max(lengths.values())

    # Misassembly breakpoints, grouped per reference sequence.
    marks: dict[str, list[tuple[int, str, str]]] = {}
    for event in getattr(reference_report, "misassemblies", None) or []:
        kind = str(getattr(event, "kind", "") or "")
        contig = str(getattr(event, "contig", "") or "")
        left_ref = getattr(event, "left_ref", None)
        right_ref = getattr(event, "right_ref", None)
        left_end = getattr(event, "left_end", None)
        right_start = getattr(event, "right_start", None)
        if left_ref in lengths and left_end is not None:
            marks.setdefault(str(left_ref), []).append((int(left_end), kind, contig))
        if right_ref in lengths and right_ref != left_ref and right_start is not None:
            marks.setdefault(str(right_ref), []).append(
                (int(right_start), kind, contig)
            )

    row_h = 46
    label_w = 150
    right_pad = 84
    width = 1000
    top = 14
    height = top + row_h * len(names) + 44
    track_w = width - label_w - right_pad
    bar_h = 16

    parts: list[str] = []
    for i, name in enumerate(names):
        y = top + i * row_h
        length = lengths[name]
        span = max(1.0, track_w * (length / longest)) if longest else track_w

        parts.append(
            f'<text class="idlabel" x="{label_w - 10}" y="{_round(y + bar_h - 3)}" '
            f'text-anchor="end">{_e(name)}</text>'
        )
        parts.append(
            f'<rect class="idtrack" x="{label_w}" y="{y}" width="{_round(span)}" '
            f'height="{bar_h}" rx="2"/>'
        )
        for a, b in blocks.get(name) or []:
            x = label_w + span * (float(a) / length if length else 0.0)
            w = span * (max(0.0, float(b) - float(a)) / length if length else 0.0)
            parts.append(
                f'<rect class="idcov" x="{_round(x)}" y="{y}" '
                f'width="{_round(max(w, 0.6))}" height="{bar_h}">'
                f"<title>{_e(name)}:{_e(fmt_int(a))}"
                f"{DASH}{_e(fmt_int(b))}</title></rect>"
            )
        for position, kind, contig in sorted(marks.get(name, [])):
            x = _round(label_w + span * (float(position) / length if length else 0.0))
            colour = KIND_COLOURS.get(kind, MUTED)
            parts.append(
                f'<line x1="{x}" y1="{_round(y - 5)}" x2="{x}" '
                f'y2="{_round(y + bar_h + 5)}" stroke="{colour}" '
                f'stroke-width="1.6"/>'
            )
            parts.append(
                f'<polygon points="{x},{_round(y - 5)} {_round(x - 4)},'
                f'{_round(y - 12)} {_round(x + 4)},{_round(y - 12)}" '
                f'fill="{colour}"><title>{_e(kind)} in {_e(contig)} at '
                f"{_e(name)}:{_e(fmt_int(position))}</title></polygon>"
            )

        pct = coverage.get(name)
        parts.append(
            f'<text class="idpct" x="{_round(label_w + span + 8)}" '
            f'y="{_round(y + bar_h - 3)}">{_e(fmt_pct(pct))}</text>'
        )
        parts.append(
            f'<text class="idsize" x="{label_w}" y="{_round(y + bar_h + 15)}">'
            f"{_e(fmt_bp(length))}</text>"
        )

    legend_y = top + row_h * len(names) + 16
    x = label_w
    for kind, colour in KIND_COLOURS.items():
        parts.append(
            f'<rect x="{_round(x)}" y="{_round(legend_y - 8)}" width="9" '
            f'height="9" fill="{colour}"/>'
        )
        parts.append(
            f'<text class="tick" x="{_round(x + 14)}" y="{_round(legend_y)}">'
            f"{_e(kind)}</text>"
        )
        x += 30 + 7.2 * len(kind)
    parts.append(
        f'<rect x="{_round(x)}" y="{_round(legend_y - 8)}" width="9" height="9" '
        f'fill="{COVERED}"/>'
    )
    parts.append(
        f'<text class="tick" x="{_round(x + 14)}" y="{_round(legend_y)}">'
        f"covered</text>"
    )

    body = "".join(parts)
    return (
        f'<svg class="chart ideogram" viewBox="0 0 {width} {_round(height)}" '
        f'width="100%" style="max-width:{width}px" '
        f'preserveAspectRatio="xMidYMid meet" '
        f'role="img" aria-label="Reference coverage ideogram" '
        f'xmlns="http://www.w3.org/2000/svg">{body}</svg>'
    )


def svg_scaffold_n50(before: int | None, after: int | None) -> str:
    """A two-bar before/after comparison of contig N50 vs scaffold N50."""
    values = [("Contig N50", before or 0), ("Scaffold N50", after or 0)]
    top = max(v for _l, v in values)
    if top <= 0:
        return _empty_chart("no scaffold lengths to compare", 160, 500)

    # Same natural width as the other paired charts, so every chart in the
    # report shrinks by the same factor on paper and the type stays consistent.
    plot = _Plot(width=500, height=270, margin=(18, 18, 46, 76), title="N50 change")
    plot.set_domain(0, 1, 0, top * 1.18)
    plot.frame([], _nice_ticks(0, top * 1.18, 4), y_label="N50 (bp)")

    baseline = plot.sy(0)
    slot = plot.plot_w / len(values)
    for i, (label, value) in enumerate(values):
        centre = plot.left + slot * (i + 0.5)
        bar_w = min(88.0, slot * 0.5)
        y = plot.sy(value)
        plot.add(
            f'<rect class="bar" x="{_round(centre - bar_w / 2)}" y="{y}" '
            f'width="{_round(bar_w)}" height="{_round(max(1.0, baseline - y))}"/>'
        )
        plot.add(
            f'<text class="note" x="{_round(centre)}" y="{_round(y - 7)}" '
            f'text-anchor="middle">{_e(fmt_bp(value))}</text>'
        )
        plot.add(
            f'<text class="tick" x="{_round(centre)}" y="{_round(baseline + 17)}" '
            f'text-anchor="middle">{_e(label)}</text>'
        )
    return plot.render()


def render_charts(
    metrics: AssemblyMetrics,
    reference_report: Any = None,
    scaffold_n50: int | None = None,
) -> dict[str, str]:
    """All the SVG the report needs, keyed by name."""
    return {
        "cumulative": svg_cumulative(getattr(metrics, "cumulative_curve", None)),
        "nx": svg_nx(getattr(metrics, "nx_curve", None), getattr(metrics, "n50", None)),
        "histogram": svg_histogram(
            getattr(metrics, "length_histogram", None),
            getattr(metrics, "largest_contig", None),
        ),
        "ideogram": svg_ideogram(reference_report) if reference_report else "",
        "scaffold_n50": (
            svg_scaffold_n50(getattr(metrics, "n50", None), scaffold_n50)
            if scaffold_n50 is not None
            else ""
        ),
    }


# --------------------------------------------------------------------------
# table shaping
# --------------------------------------------------------------------------


def _row(label: str, value: str, hint: str = "") -> dict[str, str]:
    return {"label": label, "value": value, "hint": hint}


def _basic_rows(m: AssemblyMetrics) -> list[dict[str, str]]:
    g = lambda name: getattr(m, name, None)  # noqa: E731
    rows = [
        _row("# contigs", fmt_int(g("num_contigs"))),
        _row("# contigs >= 1 kb", fmt_int(g("num_contigs_ge_1kb"))),
        _row("# contigs >= 10 kb", fmt_int(g("num_contigs_ge_10kb"))),
        _row("# contigs >= 50 kb", fmt_int(g("num_contigs_ge_50kb"))),
        _row("Total length", fmt_span(g("total_length"))),
        _row("Total length (>= 1 kb)", fmt_span(g("total_length_ge_1kb"))),
        _row("Largest contig", fmt_span(g("largest_contig"))),
        _row("Smallest contig", fmt_span(g("smallest_contig"))),
        _row("Mean contig length", fmt_float(g("mean_contig"), 1)),
        _row("Median contig length", fmt_int(g("median_contig"))),
        _row("GC content", fmt_pct(g("gc_percent"))),
        _row("# N's per 100 kb", fmt_float(g("n_per_100kb"))),
        _row("N50", fmt_span(g("n50")), "half the assembly sits in contigs this long or longer"),
        _row("L50", fmt_int(g("l50")), "contigs needed to reach half the assembly"),
        _row("N75", fmt_span(g("n75"))),
        _row("L75", fmt_int(g("l75"))),
        _row("auN", fmt_float(g("auN"), 0), "area under the Nx curve"),
    ]
    if g("genome_size"):
        rows += [
            _row("Expected genome size", fmt_span(g("genome_size"))),
            _row("NG50", fmt_span(g("ng50"))),
            _row("LG50", fmt_int(g("lg50"))),
            _row("NG75", fmt_span(g("ng75"))),
        ]
    rows += [
        _row("# links", fmt_int(g("num_links"))),
        _row("# connected components", fmt_int(g("num_components"))),
        _row("# dead ends", fmt_int(g("dead_ends"))),
        _row("# circular contigs", fmt_int(g("num_circular"))),
        _row("Mean depth", fmt_float(g("mean_depth"))),
        _row("Depth CV", fmt_float(g("depth_cv"))),
    ]
    return rows


def _reference_rows(r: Any) -> list[dict[str, str]]:
    g = lambda name: getattr(r, name, None)  # noqa: E731
    return [
        _row("Reference length", fmt_span(g("reference_length"))),
        _row("# reference sequences", fmt_int(g("reference_sequences"))),
        _row("Genome fraction", fmt_pct(g("genome_fraction")), "reference bases covered by any alignment"),
        _row("Covered reference bases", fmt_span(g("covered_bases"))),
        _row("Duplication ratio", fmt_float(g("duplication_ratio"), 3), "aligned contig bases / covered reference bases"),
        _row("# aligned contigs", fmt_int(g("aligned_contigs"))),
        _row("# unaligned contigs", fmt_int(g("unaligned_contigs"))),
        _row("# partially unaligned", fmt_int(g("partially_unaligned"))),
        _row("Fully unaligned length", fmt_span(g("fully_unaligned_length"))),
        _row("Total aligned length", fmt_span(g("total_aligned_length"))),
        _row("Largest alignment", fmt_span(g("largest_alignment"))),
        _row("NA50", fmt_span(g("na50")), "N50 over aligned blocks rather than contigs"),
        _row("NGA50", fmt_span(g("nga50"))),
        _row("LA50", fmt_int(g("la50"))),
        _row("# mismatches per 100 kb", fmt_float(g("mismatches_per_100kb"))),
        _row("# indels per 100 kb", fmt_float(g("indels_per_100kb"))),
        _row("Total mismatches", fmt_int(g("total_mismatches"))),
        _row("Total indels", fmt_int(g("total_indels"))),
        _row("# misassemblies (extensive)", fmt_int(g("num_misassemblies"))),
        _row("# relocations", fmt_int(g("num_relocations"))),
        _row("# inversions", fmt_int(g("num_inversions"))),
        _row("# translocations", fmt_int(g("num_translocations"))),
        _row("# local misassemblies", fmt_int(g("num_local_misassemblies"))),
        _row("# misassembled contigs", fmt_int(len(g("misassembled_contigs") or []))),
        _row("Misassembled contigs length", fmt_span(g("misassembled_contigs_length"))),
    ]


def _misassembly_rows(r: Any) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for event in getattr(r, "misassemblies", None) or []:
        kind = str(getattr(event, "kind", "") or "")
        left_ref = getattr(event, "left_ref", None)
        right_ref = getattr(event, "right_ref", None)
        left_end = getattr(event, "left_end", None)
        right_start = getattr(event, "right_start", None)
        if left_ref and right_ref and left_ref != right_ref:
            location = (
                f"{left_ref}:{fmt_int(left_end)} {DASH}> {right_ref}:{fmt_int(right_start)}"
            )
        elif left_ref:
            location = f"{left_ref}:{fmt_int(left_end)} {DASH}> {fmt_int(right_start)}"
        else:
            location = DASH
        out.append(
            {
                "contig": str(getattr(event, "contig", "") or DASH),
                "kind": kind or DASH,
                "colour": KIND_COLOURS.get(kind, MUTED),
                "extensive": bool(getattr(event, "is_extensive", False)),
                "contig_pos": getattr(event, "contig_pos", None),
                "contig_pos_text": fmt_int(getattr(event, "contig_pos", None)),
                "location": location,
                "ref_gap": fmt_int(getattr(event, "ref_gap", None)),
                "description": str(getattr(event, "description", "") or ""),
            }
        )
    return out


def _scaffold_context(plan: Any, built: Any, metrics: AssemblyMetrics) -> dict[str, Any]:
    """Everything the scaffolding section shows."""
    lengths: list[int] = []
    if built is not None:
        try:
            lengths = list(built.lengths)
        except Exception:
            lengths = []

    scaffold_n50 = None
    if lengths:
        from .analysis.metrics import nx_stat

        scaffold_n50 = nx_stat(lengths, 0.5)[0]

    length_of = {}
    if built is not None:
        for name, seq in getattr(built, "records", None) or []:
            length_of[name] = len(seq)

    scaffolds: list[dict[str, Any]] = []
    for scaffold in getattr(plan, "scaffolds", None) or []:
        members: list[dict[str, Any]] = []
        member_list = list(getattr(scaffold, "members", None) or [])
        for i, member in enumerate(member_list):
            last = i == len(member_list) - 1
            bridge = list(getattr(member, "bridge_path", None) or [])
            members.append(
                {
                    "segment": str(getattr(member, "segment", "") or DASH),
                    "orientation": str(getattr(member, "orientation", "+") or "+"),
                    "ref": getattr(member, "ref", None) or DASH,
                    "ref_span": (
                        f"{fmt_int(getattr(member, 'ref_start', None))}"
                        f"{DASH}{fmt_int(getattr(member, 'ref_end', None))}"
                        if getattr(member, "ref_start", None) is not None
                        else DASH
                    ),
                    "identity": (
                        fmt_pct(100.0 * getattr(member, "identity"))
                        if getattr(member, "identity", None) is not None
                        else DASH
                    ),
                    "gap_after": DASH if last else fmt_int(getattr(member, "gap_after", 0)),
                    "gap_evidence": DASH if last else str(getattr(member, "gap_evidence", "") or DASH),
                    "bridge": (
                        " ".join(f"{n}{o}" for n, o in bridge) if bridge else DASH
                    ),
                    "bridge_length": (
                        fmt_int(getattr(member, "bridge_sequence_length", 0))
                        if bridge
                        else DASH
                    ),
                    "overlaps_previous": bool(getattr(member, "overlaps_previous", False)),
                }
            )
        name = str(getattr(scaffold, "name", "") or DASH)
        scaffolds.append(
            {
                "name": name,
                "source": str(getattr(scaffold, "source", "") or DASH),
                "reference": getattr(scaffold, "reference", None) or DASH,
                "members": members,
                "member_count": len(members),
                "length": fmt_span(length_of.get(name)),
            }
        )

    all_records = len(getattr(plan, "scaffolds", None) or [])
    joined = getattr(plan, "scaffold_count", all_records)
    summary = [
        _row("Method", str(getattr(plan, "method", None) or DASH)),
        _row("# scaffolds", fmt_int(joined), "contigs joined to at least one other"),
        _row("# contigs placed", fmt_int(getattr(plan, "placed_count", None))),
        _row("# contigs unplaced", fmt_int(len(getattr(plan, "unplaced", None) or [])),
             "kept as single-contig records so the FASTA is a complete assembly"),
        _row("# contigs redundant", fmt_int(len(getattr(plan, "redundant", None) or [])),
             "another contig claimed the same reference span; not written"),
        _row("# records written", fmt_int(all_records),
             "scaffolds plus the unplaced contigs carried through"),
    ]
    if built is not None:
        summary += [
            _row("Total scaffold length", fmt_span(getattr(built, "total_length", None))),
            _row("# gaps", fmt_int(getattr(built, "num_gaps", None))),
            _row("Gap bases (N)", fmt_span(getattr(built, "gap_bases", None))),
            _row("Bridged bases (from graph)", fmt_span(getattr(built, "bridged_bases", None)),
                 "real sequence recovered from the graph instead of Ns"),
            _row("Largest scaffold", fmt_span(max(lengths) if lengths else None)),
        ]
    summary += [
        _row("Contig N50 (before)", fmt_span(getattr(metrics, "n50", None))),
        _row("Scaffold N50 (after)", fmt_span(scaffold_n50)),
    ]

    return {
        "summary": summary,
        "scaffolds": scaffolds,
        "unplaced": list(getattr(plan, "unplaced", None) or []),
        "redundant": list(getattr(plan, "redundant", None) or []),
        "notes": list(getattr(plan, "notes", None) or []),
        "warnings": list(getattr(built, "warnings", None) or []) if built is not None else [],
        "scaffold_n50": scaffold_n50,
    }


def _stat_cards(
    metrics: AssemblyMetrics,
    reference_report: Any,
    scaffold: dict[str, Any] | None,
) -> list[dict[str, str]]:
    cards = [
        {"label": "Contigs", "value": fmt_int(getattr(metrics, "num_contigs", None)), "note": f"{fmt_int(getattr(metrics, 'num_contigs_ge_1kb', None))} >= 1 kb", "tone": ""},
        {"label": "Total length", "value": fmt_bp(getattr(metrics, "total_length", None)), "note": f"largest {fmt_bp(getattr(metrics, 'largest_contig', None))}", "tone": ""},
        {"label": "N50", "value": fmt_bp(getattr(metrics, "n50", None)), "note": f"L50 {fmt_int(getattr(metrics, 'l50', None))}", "tone": ""},
        {"label": "GC", "value": fmt_pct(getattr(metrics, "gc_percent", None)), "note": f"{fmt_float(getattr(metrics, 'n_per_100kb', None))} N / 100 kb", "tone": ""},
        {"label": "Graph links", "value": fmt_int(getattr(metrics, "num_links", None)), "note": f"{fmt_int(getattr(metrics, 'num_components', None))} components, {fmt_int(getattr(metrics, 'dead_ends', None))} dead ends", "tone": ""},
    ]
    if reference_report is not None:
        frac = getattr(reference_report, "genome_fraction", None)
        mis = getattr(reference_report, "num_misassemblies", None) or 0
        cards.append(
            {
                "label": "Genome fraction",
                "value": fmt_pct(frac),
                "note": f"NGA50 {fmt_bp(getattr(reference_report, 'nga50', None))}",
                "tone": "good" if (frac or 0) >= 90 else ("warn" if (frac or 0) >= 70 else "bad"),
            }
        )
        cards.append(
            {
                "label": "Misassemblies",
                "value": fmt_int(mis),
                "note": f"{fmt_int(getattr(reference_report, 'num_local_misassemblies', None))} local",
                "tone": "good" if mis == 0 else ("warn" if mis <= 5 else "bad"),
            }
        )
        cards.append(
            {
                "label": "Mismatches / 100 kb",
                "value": fmt_float(getattr(reference_report, "mismatches_per_100kb", None)),
                "note": f"{fmt_float(getattr(reference_report, 'indels_per_100kb', None))} indels / 100 kb",
                "tone": "",
            }
        )
    if scaffold:
        cards.append(
            {
                "label": "Scaffolds",
                "value": fmt_int(len(scaffold["scaffolds"])),
                "note": f"N50 {fmt_bp(scaffold['scaffold_n50'])}",
                "tone": "",
            }
        )
    return cards


# --------------------------------------------------------------------------
# rendering
# --------------------------------------------------------------------------


def _template_dir() -> str:
    """Locate ``plastr/templates`` whether installed or run from source."""
    try:
        from importlib import resources

        path = resources.files("plastr").joinpath("templates")
        if os.path.isdir(str(path)):
            return str(path)
    except Exception:
        pass
    return os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "templates")


def _environment():
    import jinja2

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader(_template_dir()),
        autoescape=jinja2.select_autoescape(["html", "xml", "j2", "html.j2"], default=True),
        undefined=jinja2.StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    env.filters["fmt_int"] = fmt_int
    env.filters["fmt_float"] = fmt_float
    env.filters["fmt_pct"] = fmt_pct
    env.filters["fmt_bp"] = fmt_bp
    env.filters["fmt_span"] = fmt_span
    return env


def build_report(
    metrics: AssemblyMetrics,
    reference_report: Any = None,
    plan: Any = None,
    built: Any = None,
    comparison: Mapping[str, AssemblyMetrics] | None = None,
    title: str = "Plastr report",
    source_path: str | os.PathLike[str] | None = None,
    reference_path: str | os.PathLike[str] | None = None,
    subject: str | None = None,
) -> str:
    """Render a complete, standalone HTML QC report.

    Only ``metrics`` is required. Each optional argument unlocks one more
    section; passing nothing but metrics still produces a valid document.
    """
    if metrics is None:
        raise ValueError("build_report() needs an AssemblyMetrics instance")

    scaffold = None
    if plan is not None or built is not None:
        scaffold = _scaffold_context(plan, built, metrics)

    charts = render_charts(
        metrics,
        reference_report,
        scaffold["scaffold_n50"] if scaffold else None,
    )

    comparison_table = None
    if comparison:
        try:
            comparison_table = compare_metrics(dict(comparison))
        except Exception:
            comparison_table = None

    now = datetime.now().astimezone()
    context: dict[str, Any] = {
        "title": title,
        "generated": now.strftime("%Y-%m-%d %H:%M:%S %Z").strip(),
        "generated_iso": now.isoformat(timespec="seconds"),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
        "source_path": str(source_path) if source_path else None,
        "reference_path": str(reference_path) if reference_path else None,
        "metrics": metrics,
        "basic_rows": _basic_rows(metrics),
        "cards": _stat_cards(metrics, reference_report, scaffold),
        "reference": reference_report,
        "reference_rows": _reference_rows(reference_report) if reference_report is not None else [],
        "misassemblies": _misassembly_rows(reference_report) if reference_report is not None else [],
        "per_reference": (
            sorted((getattr(reference_report, "per_reference_coverage", None) or {}).items())
            if reference_report is not None
            else []
        ),
        "misassembled_contigs": (
            list(getattr(reference_report, "misassembled_contigs", None) or [])
            if reference_report is not None
            else []
        ),
        "scaffold": scaffold,
        "comparison": comparison_table,
        #: When several assemblies are compared, every tab except Comparison
        #: describes just this one, and says so rather than leaving the reader
        #: to assume the numbers cover all of them.
        "subject": subject,
        "charts": charts,
        "kind_colours": KIND_COLOURS,
        "dash": DASH,
        "version": _version(),
    }
    return _environment().get_template(TEMPLATE_NAME).render(**context)


def write_report(path: str | os.PathLike[str], **kwargs: Any) -> str:
    """Render a report and write it to ``path``; returns the path."""
    html_text = build_report(**kwargs)
    target = os.fspath(path)
    parent = os.path.dirname(os.path.abspath(target))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(target, "w", encoding="utf-8") as fh:
        fh.write(html_text)
    return target


def _version() -> str:
    try:
        from importlib.metadata import version

        return version("plastr")
    except Exception:
        return ""
