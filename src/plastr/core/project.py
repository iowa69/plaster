"""The Project: everything a working session holds.

Both the web server and the CLI drive a Project, so a workflow done by clicking
and a workflow done with subcommands produce identical output.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any

from .analysis import align as align_mod
from .analysis import misassembly as misassembly_mod
from .analysis import operations as ops_mod
from .analysis.align import Alignment
from .analysis.metrics import AssemblyMetrics, compute_metrics, nx_stat
from .analysis.misassembly import ReferenceReport
from .errors import PlastrError, GraphOperationError
from .io import gfa as gfa_mod
from .io.loader import load_graph, load_spades_paths
from .model import AssemblyGraph
from .scaffold.builder import BuiltScaffolds, build_scaffolds, write_scaffolds
from .scaffold.plan import Scaffold, ScaffoldMember, ScaffoldPlan
from .scaffold.reference_guided import scaffold_by_reference

MAX_UNDO = 40


@dataclass
class Project:
    graph: AssemblyGraph | None = None
    source_path: str | None = None
    reference_path: str | None = None
    reference_lengths: dict[str, int] = field(default_factory=dict)
    alignments: list[Alignment] = field(default_factory=list)
    reference_report: ReferenceReport | None = None
    plan: ScaffoldPlan | None = None
    genome_size: int | None = None
    align_preset: str = align_mod.DEFAULT_PRESET
    #: Segments shorter than this are excluded from statistics. Plastr
    #: counts everything by default, because in a graph a short segment is real
    #: structure; QUAST's equivalent default is 500.
    min_contig: int = 0
    #: Exclude secondary alignments from the statistics (QUAST's default
    #: handling of a repeat that maps to several places).
    primary_only: bool = True
    #: Treat reference sequences as circular, so a contig spanning the origin of
    #: a replicon is not reported as a relocation.
    circular_references: bool = True
    #: Shortest run of Ns written for a gap. AGP forbids a zero-length gap.
    min_gap: int = 100
    undo_stack: list[dict] = field(default_factory=list)
    layout: dict[str, Any] = field(default_factory=dict)
    settings: dict[str, Any] = field(default_factory=dict)

    # -- loading ------------------------------------------------------------

    def load(self, path: str, fmt: str | None = None) -> AssemblyGraph:
        graph = load_graph(path, fmt)
        self.graph = graph
        self.source_path = str(path)
        self.reference_path = None
        self.reference_lengths = {}
        self.alignments = []
        self.reference_report = None
        self.plan = None
        self.undo_stack = []
        self.layout = {}
        return graph

    def require_graph(self) -> AssemblyGraph:
        if self.graph is None:
            raise PlastrError("no assembly has been loaded yet")
        return self.graph

    def attach_paths(self, path: str) -> int:
        return load_spades_paths(self.require_graph(), path)

    # -- reference ----------------------------------------------------------

    def set_reference(
        self,
        path: str,
        preset: str = align_mod.DEFAULT_PRESET,
        min_identity: float = 0.0,
        min_length: int = 200,
        threads: int = 4,
    ) -> ReferenceReport:
        graph = self.require_graph()
        if not os.path.exists(path):
            raise PlastrError(f"reference not found: {path}")
        if os.path.isdir(path):
            raise PlastrError(f"{path} is a directory, not a reference FASTA")
        if not os.access(path, os.R_OK):
            raise PlastrError(f"cannot read {path}: permission denied")
        if not any(s.has_sequence for s in graph.segments.values()):
            raise PlastrError(
                "this graph has no sequences, so it cannot be aligned to a reference"
            )
        alignments = align_mod.align_graph(
            graph,
            path,
            preset=preset,
            min_length=max(min_length, self.min_contig),
            min_identity=min_identity,
            threads=threads,
        )
        self.reference_path = str(path)
        self.align_preset = preset
        self.reference_lengths = align_mod.reference_lengths(path)
        self.alignments = alignments
        align_mod.annotate_graph(graph, alignments)
        self.reference_report = misassembly_mod.evaluate_against_reference(
            alignments,
            self.reference_lengths,
            {
                n: s.length
                for n, s in graph.segments.items()
                if s.length >= self.min_contig
            },
            genome_size=self.genome_size,
            primary_only=self.primary_only,
            circular_references=self.circular_references,
        )
        return self.reference_report

    def clear_reference(self) -> None:
        self.reference_path = None
        self.reference_lengths = {}
        self.alignments = []
        self.reference_report = None
        if self.graph:
            for segment in self.graph.segments.values():
                segment.ref_hits = []

    def refresh_reference(self, threads: int = 4) -> ReferenceReport | None:
        """Re-align after the graph has been edited."""
        if not self.reference_path:
            return None
        return self.set_reference(
            self.reference_path, preset=self.align_preset, threads=threads
        )

    # -- metrics ------------------------------------------------------------

    def metrics(self) -> AssemblyMetrics:
        return compute_metrics(
            self.require_graph(),
            genome_size=self.effective_genome_size,
            min_length=self.min_contig,
        )

    @property
    def effective_genome_size(self) -> int | None:
        if self.genome_size:
            return self.genome_size
        if self.reference_lengths:
            return sum(self.reference_lengths.values())
        return None

    # -- graph operations ---------------------------------------------------

    def apply_operation(self, op: str, args: dict | None = None) -> dict:
        graph = self.require_graph()
        args = args or {}
        if op == "delete":
            record = ops_mod.delete_segments(graph, args.get("names", []))
        elif op == "filter":
            record = ops_mod.filter_graph(
                graph,
                min_length=int(args.get("min_length", 0) or 0),
                min_depth=args.get("min_depth"),
                max_depth=args.get("max_depth"),
                keep_components=args.get("keep_components"),
                min_component_length=int(args.get("min_component_length", 0) or 0),
            )
        elif op == "simplify":
            record = ops_mod.simplify(graph)
        elif op == "break_misassemblies":
            if not self.reference_report:
                raise GraphOperationError(
                    "load a reference first -- breaking contigs needs detected misassemblies"
                )
            record = ops_mod.break_misassemblies(
                graph,
                self.reference_report.misassemblies,
                extensive_only=bool(args.get("extensive_only", True)),
            )
        elif op == "split":
            record = ops_mod.split_segment(
                graph, args["name"], [int(p) for p in args.get("positions", [])]
            )
        elif op == "merge":
            steps = [
                (token[:-1], token[-1])
                for token in args.get("steps", [])
                if isinstance(token, str) and len(token) > 1 and token[-1] in "+-"
            ]
            record = ops_mod.merge_path(graph, steps, args.get("new_name"))
        elif op == "reverse":
            record = ops_mod.reverse_segment(graph, args["name"])
        else:
            raise GraphOperationError(f"unknown operation {op!r}")

        self.undo_stack.append(record)
        del self.undo_stack[:-MAX_UNDO]
        # Edits invalidate any plan built from the old graph.
        self.plan = None
        return record

    def undo(self) -> dict:
        if not self.undo_stack:
            raise GraphOperationError("nothing to undo")
        record = self.undo_stack.pop()
        ops_mod.restore(self.require_graph(), record)
        self.plan = None
        return record

    # -- scaffolding --------------------------------------------------------

    def build_plan(
        self,
        method: str = "reference",
        min_identity: float = 0.80,
        min_query_coverage: float = 0.30,
        min_align_length: int = 500,
        min_gap: int = 100,
        fill_gaps_from_graph: bool = True,
        include_unplaced: bool = True,
        break_misassemblies_first: bool = False,
        threads: int = 4,
    ) -> ScaffoldPlan:
        graph = self.require_graph()
        self.min_gap = max(1, int(min_gap))
        if method == "graph":
            self.plan = self._plan_from_graph()
            return self.plan

        if not self.alignments:
            raise PlastrError(
                "reference-guided scaffolding needs a reference -- load one first"
            )
        if break_misassemblies_first and self.reference_report:
            self.apply_operation("break_misassemblies", {"extensive_only": True})
            self.refresh_reference(threads=threads)

        self.plan = scaffold_by_reference(
            graph,
            self.alignments,
            min_identity=min_identity,
            min_query_coverage=min_query_coverage,
            min_align_length=min_align_length,
            min_gap=min_gap,
            fill_gaps_from_graph=fill_gaps_from_graph,
            include_unplaced=include_unplaced,
        )
        return self.plan

    def _plan_from_graph(self) -> ScaffoldPlan:
        """Scaffold using only the graph: every unbranching chain becomes one."""
        from .scaffold.graph_bridge import find_unbranching_paths

        graph = self.require_graph()
        plan = ScaffoldPlan(method="graph")
        used: set[str] = set()
        for index, chain in enumerate(find_unbranching_paths(graph, min_length=2), 1):
            scaffold = Scaffold(name=f"scaffold_graph_{index}", source="graph")
            for step_index, (name, orient) in enumerate(chain):
                member = ScaffoldMember(segment=name, orientation=orient)
                if step_index + 1 < len(chain):
                    member.gap_after = 0
                    member.gap_evidence = "adjacent"
                scaffold.members.append(member)
                used.add(name)
            plan.scaffolds.append(scaffold)

        for name in graph.segment_names_by_length():
            if name in used:
                continue
            plan.scaffolds.append(
                Scaffold(
                    name=f"scaffold_single_{name}",
                    source="graph",
                    members=[ScaffoldMember(segment=name, orientation="+")],
                )
            )
        plan.notes.append(
            f"{len(plan.scaffolds)} scaffold(s) from unbranching graph paths"
        )
        return plan

    def set_plan(self, plan: ScaffoldPlan) -> ScaffoldPlan:
        self.plan = plan
        return plan

    def build(self, min_gap: int | None = None) -> BuiltScaffolds:
        if self.plan is None:
            raise PlastrError("no scaffold plan has been built yet")
        return build_scaffolds(
            self.require_graph(),
            self.plan,
            min_gap=max(1, self.min_gap if min_gap is None else min_gap),
        )

    def preview(self) -> dict:
        built = self.build()
        lengths = built.lengths
        n50, _l50 = nx_stat(lengths, 0.5)
        return {
            "scaffold_count": len(built.records),
            "total_length": built.total_length,
            "gap_bases": built.gap_bases,
            "bridged_bases": built.bridged_bases,
            "num_gaps": built.num_gaps,
            "n50": n50,
            "largest": max(lengths) if lengths else 0,
            "warnings": built.warnings,
        }

    # -- export -------------------------------------------------------------

    def export(
        self,
        outdir: str,
        what: list[str] | None = None,
        overwrite: bool = True,
    ) -> list[dict]:
        """Write the requested artefacts into ``outdir``.

        ``overwrite=False`` refuses rather than replacing an existing file. The
        server passes False by default: an export is the one operation that can
        destroy data the user did not create, and the directory is chosen by
        whatever is driving the API.
        """
        graph = self.require_graph()
        what = what or ["scaffolds", "agp", "gfa", "csv", "report"]
        os.makedirs(outdir, exist_ok=True)

        if not overwrite:
            names = {
                "scaffolds": "scaffolds.fasta", "agp": "scaffolds.agp",
                "gfa": "graph.gfa", "csv": "segments.csv",
                "report": "report.html", "session": "session.json",
            }
            clashes = sorted(
                names[k] for k in what
                if k in names and os.path.exists(os.path.join(outdir, names[k]))
            )
            if clashes:
                raise PlastrError(
                    f"{outdir} already contains {', '.join(clashes)}. "
                    "Choose another directory, or pass --overwrite to replace them."
                )
        written: list[dict] = []

        def record(kind: str, path: str) -> None:
            written.append(
                {"kind": kind, "path": os.path.abspath(path), "bytes": os.path.getsize(path)}
            )

        if ("scaffolds" in what or "agp" in what) and self.plan is not None:
            built = self.build()
            fasta_path = os.path.join(outdir, "scaffolds.fasta")
            agp_path = os.path.join(outdir, "scaffolds.agp") if "agp" in what else None
            write_scaffolds(
                built,
                fasta_path,
                agp_path,
                comments=[
                    f"generated by Plastr from {self.source_path}",
                    f"reference: {self.reference_path or 'none'}",
                    *built.warnings,
                ],
            )
            if "scaffolds" in what:
                record("scaffolds", fasta_path)
            if agp_path:
                record("agp", agp_path)

        if "gfa" in what:
            path = os.path.join(outdir, "graph.gfa")
            gfa_mod.write_gfa(graph, path)
            record("gfa", path)

        if "csv" in what:
            path = os.path.join(outdir, "segments.csv")
            self.write_csv(path)
            record("csv", path)

        if "report" in what:
            from .report import build_report

            path = os.path.join(outdir, "report.html")
            html = build_report(
                self.metrics(),
                reference_report=self.reference_report,
                plan=self.plan,
                built=self.build() if self.plan else None,
                source_path=self.source_path,
                reference_path=self.reference_path,
            )
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(html)
            record("report", path)

        if "session" in what:
            path = os.path.join(outdir, "session.json")
            with open(path, "w") as fh:
                json.dump(self.session_dict(), fh, indent=2)
            record("session", path)

        return written

    def write_csv(self, path) -> None:
        """Write the per-segment table to a path or an open text handle."""
        import contextlib
        import csv
        import os as _os

        graph = self.require_graph()
        components = graph.component_map()
        placement = self.plan.member_index() if self.plan else {}
        opened = isinstance(path, (str, _os.PathLike))
        handle = open(path, "w", newline="") if opened else path
        with contextlib.nullcontext(handle) if not opened else handle as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "name", "length", "depth", "gc", "component",
                    "deg_start", "deg_end", "circular",
                    "ref", "ref_start", "ref_end", "ref_identity", "ref_strand",
                    "scaffold", "scaffold_position",
                ]
            )
            for name, segment in graph.segments.items():
                left, right = graph.degree(name)
                hit = segment.ref_hits[0] if segment.ref_hits else {}
                scaffold, position = placement.get(name, ("", ""))
                writer.writerow(
                    [
                        name,
                        segment.length,
                        f"{segment.depth:.4f}" if segment.depth is not None else "",
                        f"{segment.gc:.5f}" if segment.gc is not None else "",
                        components.get(name, 0),
                        left,
                        right,
                        int(graph.is_circular(name)),
                        hit.get("ref", ""),
                        hit.get("r_st", ""),
                        hit.get("r_en", ""),
                        f"{hit['identity']:.5f}" if "identity" in hit else "",
                        hit.get("strand", ""),
                        scaffold,
                        position,
                    ]
                )

    # -- session ------------------------------------------------------------

    def session_dict(self) -> dict:
        return {
            "version": 1,
            "source_path": self.source_path,
            "reference_path": self.reference_path,
            "align_preset": self.align_preset,
            "genome_size": self.genome_size,
            "layout": self.layout,
            "settings": self.settings,
            "plan": self.plan.to_dict() if self.plan else None,
        }

    def load_session(self, data: dict, reload_files: bool = True) -> None:
        if reload_files and data.get("source_path"):
            self.load(data["source_path"])
        self.genome_size = data.get("genome_size")
        self.align_preset = data.get("align_preset", align_mod.DEFAULT_PRESET)
        self.layout = data.get("layout") or {}
        self.settings = data.get("settings") or {}
        if reload_files and data.get("reference_path"):
            try:
                self.set_reference(data["reference_path"], preset=self.align_preset)
            except PlastrError:
                pass
        if data.get("plan"):
            self.plan = ScaffoldPlan.from_dict(data["plan"])

    def status(self) -> dict:
        graph = self.graph
        return {
            "loaded": graph is not None,
            "graph": graph.summary() if graph else None,
            "source_path": self.source_path,
            "source_format": graph.source_format if graph else None,
            "reference": (
                {
                    "path": self.reference_path,
                    "sequences": len(self.reference_lengths),
                    "length": sum(self.reference_lengths.values()),
                }
                if self.reference_path
                else None
            ),
            "has_alignments": bool(self.alignments),
            "has_plan": self.plan is not None,
            "align_backend": align_mod.alignment_backend(),
            "undo_depth": len(self.undo_stack),
            "settings": {
                "min_contig": self.min_contig,
                "circular_references": self.circular_references,
                "primary_only": self.primary_only,
                "genome_size": self.genome_size,
            },
        }
