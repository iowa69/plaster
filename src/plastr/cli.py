"""Command line interface.

``plastr view`` opens the interactive studio; the other subcommands run the
same operations headlessly so the tool fits into a pipeline as well as a
browser.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import webbrowser

from . import __version__
from .core.analysis import align as align_mod
from .core.analysis import search as search_mod
from .core.errors import PlastrError
from .core.project import Project

BANNER = r"""
   ___  __         __
  / _ \/ /__ ____ / /____ ____
 / ___/ / _ `(_-</ __/ -_) __/     ,----------------------------.
/_/  /_/\_,_/___/\__/\__/_/       |  ::::  [~~~~~~~~~~]  ::::  |
                                   `----------------------------'
 patching up assembly graphs -- visualise, evaluate, scaffold
"""


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _colour_enabled() -> bool:
    """Colour only when a person is actually looking at a terminal.

    Honours NO_COLOR (https://no-color.org) and FORCE_COLOR, and stays off when
    stdout is a pipe so redirected output stays clean for grep and for logs.
    """
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    if os.environ.get("TERM") == "dumb":
        return False
    return sys.stdout.isatty()


_COLOUR = None


def _c(text: str, code: str) -> str:
    global _COLOUR
    if _COLOUR is None:
        _COLOUR = _colour_enabled()
    return f"\033[{code}m{text}\033[0m" if _COLOUR else text


def _dim(text: str) -> str:
    return _c(text, "2")


def _bold(text: str) -> str:
    return _c(text, "1")


def _good(text: str) -> str:
    return _c(text, "32")


def _warn(text: str) -> str:
    return _c(text, "33")


def _bad(text: str) -> str:
    return _c(text, "31")


def _plural(n: int, singular: str, plural: str | None = None) -> str:
    return singular if n == 1 else (plural or singular + "s")


def _human(n: float | int | None) -> str:
    if n is None:
        return "-"
    n = float(n)
    for unit, size in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if abs(n) >= size:
            return f"{n / size:.2f} {unit}"
    return f"{n:,.0f} bp"


def _print_metrics(project: Project) -> None:
    m = project.metrics()
    rows = [
        ("contigs", f"{m.num_contigs:,}"),
        ("contigs >= 1 kb", f"{m.num_contigs_ge_1kb:,}"),
        ("total length", _human(m.total_length)),
        ("total length (no ovl)", _human(m.total_length_no_overlaps)),
        ("largest contig", _human(m.largest_contig)),
        ("lower quartile", _human(m.q1_length)),
        ("upper quartile", _human(m.q3_length)),
        ("N50", _human(m.n50)),
        ("L50", f"{m.l50:,}"),
        ("N75", _human(m.n75)),
        ("auN", f"{m.auN:,.0f}"),
    ]
    if m.ng50 is not None:
        rows.append(("NG50", _human(m.ng50)))
    if m.gc_percent is not None:
        rows.append(("GC", f"{m.gc_percent:.2f} %"))
    if m.median_depth is not None:
        rows.append(("median depth", f"{m.median_depth:.1f} x"))
    if m.mean_depth is not None:
        rows.append(("mean depth", f"{m.mean_depth:.1f} x"))
    if m.estimated_sequence_length is not None:
        rows.append(("estimated sequence", _human(m.estimated_sequence_length)))
    if m.overlap_min is None:
        overlaps = "-"
    elif m.overlap_min == m.overlap_max:
        overlaps = f"{m.overlap_min:,} bp"
    else:
        overlaps = f"{m.overlap_min:,} to {m.overlap_max:,} bp"
    rows += [
        ("links", f"{m.num_links:,}"),
        ("edge overlap range", overlaps),
        ("components", f"{m.num_components:,}"),
        ("largest component", _human(m.largest_component_length) + (
            f" ({m.largest_component_percent:.1f} %)"
            if m.largest_component_percent is not None
            else ""
        )),
        ("orphaned length", _human(m.orphaned_length)),
        ("dead ends", f"{m.dead_ends:,}"),
        ("percentage dead ends", f"{m.dead_end_percent:.2f} %"),
        ("circular contigs", f"{m.num_circular:,}"),
    ]
    print("\n  " + _bold("Assembly statistics"))
    print(_dim("  " + "-" * 40))
    for label, value in rows:
        print(f"  {_dim(label.ljust(22))}{value:>18}")


def _print_reference(project: Project) -> None:
    r = project.reference_report
    if r is None:
        return
    print("\n  " + _bold("Reference-based statistics"))
    print(_dim("  " + "-" * 40))
    def gf_tone(text: str) -> str:
        if r.genome_fraction >= 90:
            return _good(text)
        return _warn(text) if r.genome_fraction >= 70 else _bad(text)

    def mis_tone(text: str) -> str:
        return _good(text) if r.num_misassemblies == 0 else _bad(text)

    # (label, value, tone). Colour is applied after padding, because an ANSI
    # escape is characters too and would push every column out of line.
    rows = [
        ("genome fraction", f"{r.genome_fraction:.2f} %", gf_tone),
        ("duplication ratio", f"{r.duplication_ratio:.3f}", None),
        ("largest alignment", _human(r.largest_alignment), None),
        ("NA50", _human(r.na50), None),
        ("NGA50", _human(r.nga50) if r.nga50 else "-", None),
        ("mismatches / 100 kb", f"{r.mismatches_per_100kb:.2f}", None),
        ("indels / 100 kb", f"{r.indels_per_100kb:.2f}", None),
        ("misassemblies", f"{r.num_misassemblies:,}", mis_tone),
        ("  relocations", f"{r.num_relocations:,}", None),
        ("  inversions", f"{r.num_inversions:,}", None),
        ("  translocations", f"{r.num_translocations:,}", None),
        ("local misassemblies", f"{r.num_local_misassemblies:,}", None),
        ("unaligned contigs", f"{r.unaligned_contigs:,}", None),
    ]
    for label, value, tone in rows:
        cell = value.rjust(18)
        print(f"  {_dim(label.ljust(22))}{tone(cell) if tone else cell}")
    if r.misassemblies:
        print("\n  " + _bold("Misassembly detail"))
        print(_dim("  " + "-" * 40))
        for event in r.misassemblies[:25]:
            flag = _bad("!") if event.is_extensive else " "
            print(f"  {flag} {event.contig:<24} {event.kind:<14} {_dim(event.description)}")
        if len(r.misassemblies) > 25:
            print(_dim(f"    ... and {len(r.misassemblies) - 25} more"))


def _load(args) -> Project:
    project = Project()
    project.load(args.assembly, getattr(args, "format", None))
    if getattr(args, "paths", None):
        added = project.attach_paths(args.paths)
        print(f"  attached {added} path(s) from {args.paths}")
    if getattr(args, "genome_size", None):
        project.genome_size = args.genome_size
    project.min_contig = int(getattr(args, "min_contig", 0) or 0)
    project.circular_references = not getattr(args, "linear_references", False)
    project.primary_only = not getattr(args, "count_ambiguous", False)
    summary = project.graph.summary() if project.graph else {}
    print(
        f"  loaded {args.assembly}: {summary.get('segments', 0):,} segments, "
        f"{summary.get('links', 0):,} links, {_human(summary.get('total_length'))}"
    )
    if not summary.get("has_sequences"):
        print("  note: this file has no sequences; export and scaffolding are unavailable")
    return project


def _maybe_reference(project: Project, args) -> None:
    if not getattr(args, "reference", None):
        return
    backend = align_mod.alignment_backend()
    if backend == "none":
        print("  ! minimap2/mappy not found -- skipping reference analysis", file=sys.stderr)
        return
    print(f"  aligning to {args.reference} (preset {args.preset}, backend {backend}) ...")
    project.set_reference(
        args.reference, preset=args.preset, threads=args.threads
    )
    print(f"  {len(project.alignments):,} alignment blocks")


# ---------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------


def cmd_view(args) -> int:
    import uvicorn

    from .server.app import create_app

    project = Project()
    if args.assembly:
        project.load(args.assembly, args.format)
        if args.paths:
            project.attach_paths(args.paths)
        if args.genome_size:
            project.genome_size = args.genome_size
        project.min_contig = int(getattr(args, "min_contig", 0) or 0)
        project.circular_references = not getattr(args, "linear_references", False)
        project.primary_only = not getattr(args, "count_ambiguous", False)
        summary = project.graph.summary() if project.graph else {}
        print(
            f"  loaded {args.assembly}: {summary.get('segments', 0):,} segments, "
            f"{summary.get('links', 0):,} links"
        )
        if args.reference:
            _maybe_reference(project, args)
        if getattr(args, "diff", None):
            print(f"    comparing against {args.diff} ...")
            result = project.set_comparison(
                args.diff, preset=args.preset, threads=args.threads
            )
            counts = result.a.counts()
            print(f"    {counts['shared']:,} shared, {counts['partial']:,} partial, "
                  f"{counts['missing']:,} not in {os.path.basename(args.diff)}"
                  "   (colour by 'comparison' to see them)")

    if not 1 <= args.port <= 65535:
        print(f"  error: --port must be between 1 and 65535, got {args.port}", file=sys.stderr)
        return 2

    app = create_app(project, threads=args.threads)
    url = f"http://{args.host}:{args.port}"

    config = uvicorn.Config(app, host=args.host, port=args.port, log_level=args.log_level)
    server = uvicorn.Server(config)
    # Claim the port before announcing anything. Printing "Plastr is running
    # at ..." and opening a browser first meant that when the port was already
    # taken, the banner still appeared and the browser opened on somebody
    # else's session -- a different assembly, presented as yours.
    try:
        # uvicorn logs the OSError and raises SystemExit rather than letting the
        # error out, so catch both and say something useful instead.
        sockets = [config.bind_socket()]
    except (OSError, SystemExit) as exc:
        detail = exc if isinstance(exc, OSError) else "address already in use"
        print(
            f"\nerror: cannot listen on {args.host}:{args.port}: {detail}."
            "\n       Another Plastr may already be running there; try --port.\n",
            file=sys.stderr,
        )
        return 1

    print(BANNER)
    print(f"  Plastr is running at {url}")
    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print()
        print(f"  !  Bound to {args.host}, so this is reachable from other machines.")
        print("  !  There is no authentication: anyone who can reach this port can")
        print("  !  browse your filesystem, read files, and write exports as you.")
        print("  !  Use the default 127.0.0.1 unless you trust the whole network.")
    print("  press Ctrl+C to stop\n")

    if not args.no_browser:
        timer = threading.Timer(1.0, lambda: webbrowser.open(url))
        timer.daemon = True
        timer.start()

    server.run(sockets=sockets)
    return 0


def cmd_qc(args) -> int:
    project = _load(args)
    _maybe_reference(project, args)
    _print_metrics(project)
    _print_reference(project)

    if args.html:
        from .core.report import build_report

        html = build_report(
            project.metrics(),
            reference_report=project.reference_report,
            source_path=project.source_path,
            reference_path=project.reference_path,
            graph=project.graph,
            min_contig=project.min_contig,
        )
        # The report always contains an en dash and declares UTF-8, so it must
        # be written as UTF-8 whatever the host locale says.
        with open(args.html, "w", encoding="utf-8") as fh:
            fh.write(html)
        print(f"\n  report written to {args.html}")
    if args.json:
        payload = {
            "metrics": project.metrics().to_dict(),
            "reference": (
                project.reference_report.to_dict() if project.reference_report else None
            ),
        }
        with open(args.json, "w") as fh:
            json.dump(payload, fh, indent=2)
        print(f"  json written to {args.json}")
    print()
    return 0


def cmd_scaffold(args) -> int:
    project = _load(args)
    if args.method == "reference":
        if not args.reference:
            print("  error: --reference is required for reference-guided scaffolding", file=sys.stderr)
            return 2
        _maybe_reference(project, args)
        if not project.alignments:
            print("  error: alignment produced no results", file=sys.stderr)
            return 1
    elif args.reference:
        # Graph scaffolding does not use the reference to order contigs, but
        # aligning anyway lets --break-misassemblies work and fills in the QC
        # section of the report. Silently ignoring -r would be worse.
        _maybe_reference(project, args)
    elif args.break_misassemblies:
        print(
            "  error: --break-misassemblies needs a reference; pass -r, or drop the flag",
            file=sys.stderr,
        )
        return 2

    if args.break_misassemblies and project.reference_report:
        before = project.reference_report.num_misassemblies
        record = project.apply_operation("break_misassemblies", {"extensive_only": True})
        print(
            f"  broke {record['count']} misassembled contig(s) into "
            f"{record['pieces']} pieces (was {before} misassemblies)"
        )
        project.refresh_reference(threads=args.threads)
        if project.reference_report:
            print(f"  after breaking: {project.reference_report.num_misassemblies} misassemblies")

    project.build_plan(
        method=args.method,
        min_identity=args.min_identity,
        min_query_coverage=args.min_coverage,
        min_align_length=args.min_align_length,
        min_gap=args.min_gap,
        fill_gaps_from_graph=not args.no_graph_gapfill,
        include_unplaced=not args.drop_unplaced,
        threads=args.threads,
    )
    plan = project.plan
    assert plan is not None
    for note in plan.notes:
        print(f"  {note}")

    preview = project.preview()
    carried = preview["scaffold_count"] - plan.scaffold_count
    print("\n  " + _bold("Scaffolds"))
    print(_dim("  " + "-" * 40))
    rows = [
        ("scaffolds", f"{plan.scaffold_count:,}"),
        ("  contigs joined", f"{plan.placed_count:,}"),
        ("unplaced, kept as is", f"{carried:,}"),
        ("records written", f"{preview['scaffold_count']:,}"),
        ("total length", _human(preview["total_length"])),
        ("largest", _human(preview["largest"])),
        ("N50", _human(preview["n50"])),
        ("gap bases (N)", _human(preview["gap_bases"])),
        ("bases from graph", _human(preview["bridged_bases"])),
    ]
    for label, value in rows:
        print(f"  {_dim(label.ljust(22))}{value:>18}")
    for warning in preview["warnings"]:
        print(f"  {_warn('!')} {warning}", file=sys.stderr)

    written = project.export(args.outdir, args.export, overwrite=args.overwrite)
    print(f"\n  written to {os.path.abspath(args.outdir)}/")
    for item in written:
        print(f"    {item['kind']:<12} {os.path.basename(item['path']):<20} {item['bytes']:>12,} bytes")
    print()
    return 0


def cmd_export(args) -> int:
    project = _load(args)
    _maybe_reference(project, args)
    if "scaffolds" in args.export or "agp" in args.export:
        if project.reference_path:
            project.build_plan(method="reference", threads=args.threads)
        else:
            project.build_plan(method="graph")
    written = project.export(args.outdir, args.export, overwrite=args.overwrite)
    print(f"  written to {os.path.abspath(args.outdir)}/")
    for item in written:
        print(f"    {item['kind']:<12} {os.path.basename(item['path']):<20} {item['bytes']:>12,} bytes")
    return 0


def cmd_search(args) -> int:
    project = _load(args)
    backend, hits = search_mod.search_graph(
        project.require_graph(),
        args.query,
        min_identity=args.min_identity,
        threads=args.threads,
    )
    shown = min(len(hits), args.max_hits)
    extra = "" if shown == len(hits) else f" (showing the top {shown})"
    print(f"  backend: {backend}; {len(hits)} hit(s){extra}\n")
    print(f"  {'segment':<28}{'ident':>8}{'len':>9}{'start':>10}{'end':>10}{'strand':>8}")
    for hit in hits[: args.max_hits]:
        print(
            f"  {hit.segment:<28}{hit.identity * 100:>7.2f}%{hit.length:>9,}"
            f"{hit.s_st:>10,}{hit.s_en:>10,}{'+' if hit.strand > 0 else '-':>8}"
        )
    return 0


def cmd_doctor(_args) -> int:
    """Report what is installed and what each missing piece would cost you."""
    from pathlib import Path

    print(BANNER)
    ok = True

    def line(label: str, good: bool, detail: str) -> None:
        mark = "ok  " if good else "MISS"
        print(f"  [{mark}] {label:<22} {detail}")

    print("  Environment\n  " + "-" * 60)
    line("python", True, sys.version.split()[0])
    line("plastr", True, __version__)

    backend = align_mod.alignment_backend()
    if backend == "none":
        ok = False
        line(
            "minimap2 / mappy", False,
            "no reference alignment: no genome fraction, misassemblies, or "
            "reference scaffolding",
        )
    else:
        line("minimap2 / mappy", True, f"using {backend}")

    if search_mod.blast_available():
        line("blast", True, "sequence search enabled")
    else:
        line("blast", True, "not found -- search will fall back to minimap2")

    try:
        import uvicorn  # noqa: F401
        import fastapi  # noqa: F401

        line("web service", True, "fastapi + uvicorn available")
    except ImportError:
        ok = False
        line("web service", False, "'plastr view' will not start")

    web_dir = Path(__file__).resolve().parent / "web"
    has_ui = (web_dir / "index.html").exists()
    if has_ui:
        line("web interface", True, str(web_dir))
    else:
        ok = False
        line("web interface", False, f"missing files under {web_dir}")

    try:
        from .core.report import build_report  # noqa: F401

        line("html reports", True, "available")
    except ImportError as exc:
        ok = False
        line("html reports", False, str(exc))

    print()
    if ok:
        print("  Everything needed is installed.\n")
        print("  Try:  python examples/make_demo_data.py")
        print("        plastr view examples/demo/assembly.gfa -r examples/demo/reference.fasta\n")
    else:
        print("  Some pieces are missing. The conda environment installs all of them:\n")
        print("        conda env create -f environment.yml && conda activate plastr\n")
    return 0 if ok else 1


def _print_diff_side(side, other_label: str) -> None:
    counts = side.counts()
    print(f"\n  {side.label}")
    print("  " + "-" * 56)
    print(f"    contigs                        {side.num_contigs:>12,}")
    print(f"    total length                   {_human(side.total_length):>12}")
    print(f"    also in {other_label[:20]:<22}{_human(side.shared_bases):>12}"
          f"  ({100.0 * side.shared_fraction:.2f} %)")
    print(f"    not in {other_label[:20]:<23}{_human(side.unique_bases):>12}")
    print(f"    contigs shared / partial / missing"
          f"{counts['shared']:>7,} /{counts['partial']:>6,} /{counts['missing']:>6,}")
    if side.without_sequence:
        print(f"    contigs without sequence       {side.without_sequence:>12,}"
              "   (not compared)")

    gone = side.missing_components()
    if gone:
        print(f"\n    Components not fully in {other_label}")
        print("    " + "-" * 54)
        for comp in gone[:12]:
            shape = "circular" if comp.circular else "linear"
            print(f"    ! {len(comp.contigs):>3} contig(s)  {_human(comp.length):>10}"
                  f"  {shape:<9} {comp.status:<8}"
                  f" {100.0 * comp.covered_fraction:5.1f} % covered")
        if len(gone) > 12:
            print(f"      ... and {len(gone) - 12} more")

    rows = side.missing_contigs()
    if rows:
        print(f"\n    Longest contigs not fully in {other_label}")
        print("    " + "-" * 54)
        for c in rows[:15]:
            print(f"    {c.status:<8} {c.name[:26]:<26} {_human(c.length):>10}"
                  f" {100.0 * c.covered_fraction:5.1f} % covered")
        if len(rows) > 15:
            print(f"      ... and {len(rows) - 15} more")


def cmd_diff(args) -> int:
    """Compare two assemblies by content: what is in one and not the other."""
    from .core.analysis.diff import diff_graphs

    projects = []
    for path in (args.a, args.b):
        project = Project()
        project.load(path)
        summary = project.graph.summary() if project.graph else {}
        print(f"  loaded {path}: {summary.get('segments', 0):,} segments, "
              f"{_human(summary.get('total_length'))}")
        projects.append(project)

    label_a = args.label_a or os.path.splitext(os.path.basename(args.a))[0]
    label_b = args.label_b or os.path.splitext(os.path.basename(args.b))[0]
    if label_a == label_b:
        label_a, label_b = f"{label_a} (1)", f"{label_b} (2)"

    print(f"\n  comparing sequence, {args.preset} ...")
    result = diff_graphs(
        projects[0].graph, projects[1].graph,
        label_a=label_a, label_b=label_b,
        path_a=args.a, path_b=args.b,
        preset=args.preset,
        min_identity=args.min_identity,
        threads=args.threads,
        shared_fraction=args.shared_fraction,
        present_fraction=args.present_fraction,
    )

    _print_diff_side(result.a, label_b)
    _print_diff_side(result.b, label_a)

    if args.json:
        import json
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(result.to_dict(), fh, indent=2)
        print(f"\n  json written to {args.json}")

    if args.html:
        from .core.report import build_diff_report
        html = build_diff_report(result)
        with open(args.html, "w", encoding="utf-8") as fh:
            fh.write(html)
        print(f"  report written to {args.html}")

    return 0


def cmd_compare(args) -> int:
    """Evaluate several assemblies side by side, QUAST-style."""
    from .core.analysis.metrics import compare_metrics

    projects: dict[str, Project] = {}
    for path in args.assemblies:
        label = os.path.splitext(os.path.basename(path))[0]
        suffix = 2
        while label in projects:
            label = f"{os.path.splitext(os.path.basename(path))[0]}_{suffix}"
            suffix += 1
        project = Project()
        project.load(path)
        if args.genome_size:
            project.genome_size = args.genome_size
        project.min_contig = int(getattr(args, "min_contig", 0) or 0)
        summary = project.graph.summary() if project.graph else {}
        print(
            f"  loaded {path}: {summary.get('segments', 0):,} segments, "
            f"{_human(summary.get('total_length'))}"
        )
        if args.reference and align_mod.alignment_backend() != "none":
            print(f"    aligning to {args.reference} ...")
            project.set_reference(args.reference, preset=args.preset, threads=args.threads)
        projects[label] = project

    named = {label: p.metrics() for label, p in projects.items()}
    table = compare_metrics(named)

    col_width = max(14, max((len(c) for c in table["columns"]), default=14) + 2)
    print("\n  Comparison")
    print("  " + "-" * (30 + col_width * len(table["columns"])))
    header = f"  {'metric':<28}" + "".join(f"{c:>{col_width}}" for c in table["columns"])
    print(header)
    for row in table["rows"]:
        print(
            f"  {row['label']:<28}"
            + "".join(f"{v:>{col_width}}" for v in row["values"])
        )

    if args.reference and any(p.reference_report for p in projects.values()):
        print("\n  Reference-based")
        print("  " + "-" * (30 + col_width * len(table["columns"])))
        ref_rows = [
            ("genome fraction (%)", lambda r: f"{r.genome_fraction:.2f}"),
            ("duplication ratio", lambda r: f"{r.duplication_ratio:.3f}"),
            ("NGA50", lambda r: _human(r.nga50) if r.nga50 else "-"),
            ("misassemblies", lambda r: f"{r.num_misassemblies:,}"),
            ("mismatches / 100 kb", lambda r: f"{r.mismatches_per_100kb:.1f}"),
            ("unaligned contigs", lambda r: f"{r.unaligned_contigs:,}"),
        ]
        for label, fmt in ref_rows:
            cells = []
            for name in table["columns"]:
                report = projects[name].reference_report
                cells.append(fmt(report) if report else "-")
            print(f"  {label:<28}" + "".join(f"{v:>{col_width}}" for v in cells))

    if args.html:
        from .core.report import build_report

        first_label = next(iter(projects))
        first = projects[first_label]
        html = build_report(
            first.metrics(),
            reference_report=first.reference_report,
            comparison=named,
            title="Plastr comparison",
            # Only the Comparison tab covers every assembly; the rest describe
            # this one, so name it rather than listing them all in the header.
            source_path=first.source_path,
            reference_path=args.reference,
            subject=first_label,
            graph=first.graph,
            min_contig=first.min_contig,
        )
        with open(args.html, "w", encoding="utf-8") as fh:
            fh.write(html)
        print(f"\n  report written to {args.html}")
    print()
    return 0


def cmd_info(args) -> int:
    project = _load(args)
    _print_metrics(project)
    graph = project.require_graph()
    components = graph.connected_components()
    print("\n  " + _bold("Largest components"))
    print(_dim("  " + "-" * 40))
    for i, comp in enumerate(components[:10]):
        total = sum(graph.segments[n].length for n in comp)
        label = f"{len(comp):,} {_plural(len(comp), 'segment')}"
        print(f"  {_dim(f'{i:>3}')}  {label:<18}{_human(total):>14}")
    if len(components) > 10:
        print(_dim(f"       ... and {len(components) - 10:,} more"))
    print()
    return 0


# ---------------------------------------------------------------------------
# argument parsing
# ---------------------------------------------------------------------------


def _add_common(parser: argparse.ArgumentParser, assembly_required: bool = True) -> None:
    parser.add_argument(
        "assembly",
        nargs=None if assembly_required else "?",
        help="assembly graph (GFA/GFA2/FASTG) or contig FASTA",
    )
    parser.add_argument("--format", choices=["gfa", "gfa2", "fastg", "fasta"], help="override format detection")
    parser.add_argument("--paths", help="SPAdes contigs.paths to attach")
    parser.add_argument("--genome-size", type=int, help="expected genome size, for NG50/NGA50")
    parser.add_argument(
        "--min-contig",
        type=int,
        default=0,
        metavar="BP",
        help="ignore segments shorter than this in the statistics "
        "(default 0, count everything; use 500 to match QUAST)",
    )
    parser.add_argument(
        "--linear-references",
        action="store_true",
        help="treat reference sequences as linear; by default they are assumed "
        "circular, so a contig spanning a replicon's origin is not a misassembly",
    )
    parser.add_argument(
        "--count-ambiguous",
        action="store_true",
        help="include secondary alignments in the statistics, so a repeat counts "
        "at every place it maps (default: primary alignments only, as QUAST does)",
    )
    parser.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4)


def _add_overwrite(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace files that already exist in the output directory "
        "(by default the export stops rather than destroying them)",
    )


def _add_reference(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-r", "--reference", help="reference FASTA")
    parser.add_argument(
        "--preset",
        choices=list(align_mod.PRESETS),
        default=align_mod.DEFAULT_PRESET,
        help="minimap2 preset: asm5 same strain, asm10 same species, asm20 related species",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plastr",
        description="Interactive assembly graph studio: visualise, evaluate, and scaffold.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  plastr view assembly.gfa\n"
            "  plastr view assembly.gfa -r reference.fasta\n"
            "  plastr qc assembly.gfa -r reference.fasta --html report.html\n"
            "  plastr scaffold assembly.gfa -r reference.fasta -o out/ --break-misassemblies\n"
            "  plastr search assembly.gfa --query gene.fasta\n"
            "  plastr compare spades.gfa mine.gfa -r reference.fasta --html compare.html\n"
        ),
    )
    parser.add_argument("--version", action="version", version=f"plastr {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    p_view = sub.add_parser("view", help="open the interactive studio in a browser")
    _add_common(p_view, assembly_required=False)
    _add_reference(p_view)
    p_view.add_argument(
        "--host",
        default="127.0.0.1",
        help="address to bind. The default keeps the studio on this machine. "
        "Binding elsewhere exposes an unauthenticated API that can read and "
        "write files as you -- only do it on a network you trust",
    )
    p_view.add_argument("--diff", metavar="OTHER.gfa",
                        help="compare against another assembly and colour by what it lacks")
    p_view.add_argument("--port", type=int, default=8781, help="port to listen on")
    p_view.add_argument("--no-browser", action="store_true", help="do not open a browser")
    p_view.add_argument(
        "--log-level",
        default="warning",
        choices=["critical", "error", "warning", "info", "debug", "trace"],
    )
    p_view.set_defaults(func=cmd_view)

    p_qc = sub.add_parser("qc", help="assembly statistics, with reference evaluation")
    _add_common(p_qc)
    _add_reference(p_qc)
    p_qc.add_argument("--html", help="write an HTML report here")
    p_qc.add_argument("--json", help="write the numbers as JSON here")
    p_qc.set_defaults(func=cmd_qc)

    p_sc = sub.add_parser("scaffold", help="build scaffolds and export FASTA + AGP")
    _add_common(p_sc)
    _add_reference(p_sc)
    p_sc.add_argument("-o", "--outdir", default="plastr_out")
    p_sc.add_argument("--method", choices=["reference", "graph"], default="reference")
    p_sc.add_argument("--min-identity", type=float, default=0.80)
    p_sc.add_argument("--min-coverage", type=float, default=0.30, help="minimum fraction of a contig that must align")
    p_sc.add_argument("--min-align-length", type=int, default=500)
    p_sc.add_argument("--min-gap", type=int, default=100)
    p_sc.add_argument("--no-graph-gapfill", action="store_true", help="do not fill gaps with graph sequence")
    p_sc.add_argument("--drop-unplaced", action="store_true", help="omit contigs that could not be placed")
    p_sc.add_argument(
        "--break-misassemblies",
        action="store_true",
        help="split contigs at detected misassemblies before scaffolding",
    )
    p_sc.add_argument(
        "--export",
        nargs="+",
        default=["scaffolds", "agp", "gfa", "csv", "report"],
        choices=["scaffolds", "agp", "gfa", "csv", "report", "session"],
    )
    _add_overwrite(p_sc)
    p_sc.set_defaults(func=cmd_scaffold)

    p_ex = sub.add_parser("export", help="write artefacts without the GUI")
    _add_common(p_ex)
    _add_reference(p_ex)
    p_ex.add_argument("-o", "--outdir", default="plastr_out")
    p_ex.add_argument(
        "--export",
        nargs="+",
        default=["gfa", "csv", "report"],
        choices=["scaffolds", "agp", "gfa", "csv", "report", "session"],
    )
    _add_overwrite(p_ex)
    p_ex.set_defaults(func=cmd_export)

    p_se = sub.add_parser("search", help="find a sequence in the assembly")
    _add_common(p_se)
    p_se.add_argument("-q", "--query", required=True, help="sequence, FASTA text, or a FASTA path")
    p_se.add_argument("--min-identity", type=float, default=0.8)
    p_se.add_argument("--max-hits", type=int, default=50)
    p_se.set_defaults(func=cmd_search)

    p_cmp = sub.add_parser("compare", help="evaluate several assemblies side by side")
    p_cmp.add_argument("assemblies", nargs="+", help="two or more assemblies to compare")
    _add_reference(p_cmp)
    p_cmp.add_argument("--genome-size", type=int, help="expected genome size, for NG50/NGA50")
    p_cmp.add_argument("--min-contig", type=int, default=0, metavar="BP",
                       help="ignore segments shorter than this (use 500 to match QUAST)")
    p_cmp.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4)
    p_cmp.add_argument("--html", help="write a comparison report here")
    p_cmp.set_defaults(func=cmd_compare)

    p_diff = sub.add_parser(
        "diff", help="compare two assemblies by content: what is in one and not the other")
    p_diff.add_argument("a", help="first assembly")
    p_diff.add_argument("b", help="second assembly")
    p_diff.add_argument("--label-a", help="name for the first assembly in the report")
    p_diff.add_argument("--label-b", help="name for the second assembly in the report")
    p_diff.add_argument("--preset", default="asm10", help="minimap2 preset (default asm10)")
    p_diff.add_argument("--min-identity", type=float, default=0.0, metavar="F",
                        help="ignore alignments below this identity, 0-1")
    p_diff.add_argument("--shared-fraction", type=float, default=0.95, metavar="F",
                        help="covered at least this much counts as shared (default 0.95)")
    p_diff.add_argument("--present-fraction", type=float, default=0.10, metavar="F",
                        help="below this counts as missing rather than partial (default 0.10)")
    p_diff.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4)
    p_diff.add_argument("--html", help="write a comparison report here")
    p_diff.add_argument("--json", help="write the full comparison as JSON here")
    p_diff.set_defaults(func=cmd_diff)

    p_doc = sub.add_parser("doctor", help="check that everything is installed")
    p_doc.set_defaults(func=cmd_doctor)

    p_in = sub.add_parser("info", help="quick summary of a graph")
    _add_common(p_in)
    p_in.set_defaults(func=cmd_info)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except PlastrError as exc:
        print(f"\nerror: {exc}\n", file=sys.stderr)
        return 1
    except (OSError, RuntimeError, ValueError) as exc:
        # A missing reference, an unwritable output path or a failed aligner
        # should read as an error, not as a stack trace.
        print(f"\nerror: {exc}\n", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
