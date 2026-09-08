"""The Plastr HTTP service.

One project lives in memory; the browser drives it. See ``docs/API.md`` for the
contract. There is no biology in this module -- it validates input, calls the
core library, and shapes JSON.
"""

from __future__ import annotations

import io
import json
import os
import stat
import threading
from pathlib import Path
from typing import Any

from fastapi import Body, FastAPI, HTTPException, Query, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from starlette.exceptions import HTTPException as StarletteHTTPException

from ..core.analysis import align as align_mod
from ..core.analysis import search as search_mod
from ..core.errors import PlastrError, PlastrFormatError, MissingDependencyError
from ..core.io import gfa as gfa_mod
from ..core.project import Project
from ..core.scaffold.plan import ScaffoldPlan

WEB_DIR = Path(__file__).resolve().parent.parent / "web"

# The project is mutated from request handlers; a lock keeps a slow alignment
# from racing a graph edit.
_lock = threading.RLock()


def _num(body: dict, key: str, default, cast, low=None, high=None, label=None):
    """Read a numeric field, honouring an explicit zero.

    ``body.get(key, default) or default`` looks harmless and silently turns
    every explicit 0 into the default -- a request for ``min_gap: 0`` came back
    byte-identical to ``min_gap: 100``.
    """
    value = body.get(key)
    if value is None or value == "":
        return default
    try:
        out = cast(value)
    except (TypeError, ValueError):
        raise HTTPException(
            status_code=400,
            detail={"error": f"{label or key} must be a number, got {value!r}", "detail": ""},
        ) from None
    if (low is not None and out < low) or (high is not None and out > high):
        span = f"between {low} and {high}" if low is not None and high is not None else (
            f"at least {low}" if high is None else f"at most {high}"
        )
        raise HTTPException(
            status_code=400,
            detail={"error": f"{label or key} must be {span}, got {out}", "detail": ""},
        )
    return out


def create_app(project: Project | None = None, threads: int = 4) -> FastAPI:
    app = FastAPI(title="Plastr", version="1.0.0", docs_url="/api/docs")
    app.state.project = project or Project()
    app.state.threads = threads

    # ---------------------------------------------------------------- errors

    @app.exception_handler(PlastrFormatError)
    async def _format_error(_request: Request, exc: PlastrFormatError):
        return JSONResponse(
            status_code=400, content={"error": str(exc), "detail": "format error"}
        )

    @app.exception_handler(MissingDependencyError)
    async def _missing_dep(_request: Request, exc: MissingDependencyError):
        return JSONResponse(
            status_code=503, content={"error": str(exc), "detail": "missing dependency"}
        )

    @app.exception_handler(PlastrError)
    async def _general(_request: Request, exc: PlastrError):
        return JSONResponse(status_code=400, content={"error": str(exc), "detail": ""})

    @app.exception_handler(HTTPException)
    async def _http_error(_request: Request, exc: HTTPException):
        detail = exc.detail
        if isinstance(detail, dict) and "error" in detail:
            content = {"error": detail["error"], "detail": detail.get("detail", "")}
        else:
            content = {"error": str(detail), "detail": ""}
        return JSONResponse(status_code=exc.status_code, content=content)

    @app.exception_handler(RequestValidationError)
    async def _validation(_request: Request, exc: RequestValidationError):
        # docs/API.md promises {"error", "detail"} on every 4xx. FastAPI's own
        # 422 body is a list under "detail", which the UI rendered as
        # "HTTP 422 Unprocessable Entity" or "[object Object]".
        first = (exc.errors() or [{}])[0]
        where = ".".join(str(p) for p in first.get("loc", ()) if p != "body") or "request"
        return JSONResponse(
            status_code=400,
            content={
                "error": f"{where}: {first.get('msg', 'invalid value')}",
                "detail": "validation error",
            },
        )

    @app.exception_handler(StarletteHTTPException)
    async def _starlette_error(_request: Request, exc: StarletteHTTPException):
        detail = exc.detail
        if isinstance(detail, dict) and "error" in detail:
            content = {"error": detail["error"], "detail": detail.get("detail", "")}
        else:
            content = {"error": str(detail), "detail": ""}
        return JSONResponse(status_code=exc.status_code, content=content)

    @app.exception_handler(Exception)
    async def _unexpected(_request: Request, exc: Exception):
        # An unhandled OSError would otherwise surface as a plain-text 500 that
        # the UI shows as "HTTP 500 Internal Server Error".
        return JSONResponse(
            status_code=500,
            content={"error": f"{type(exc).__name__}: {exc}", "detail": "unexpected error"},
        )

    def proj() -> Project:
        return app.state.project

    def fail(message: str, status: int = 400) -> None:
        # docs/API.md promises {"error", "detail"} and the UI reads `error`.
        # FastAPI's default HTTPException body is {"detail": ...}, which would
        # demote the real message behind a bare "HTTP 400".
        raise HTTPException(status_code=status, detail={"error": message, "detail": ""})

    # ---------------------------------------------------------------- status

    @app.get("/api/status")
    def status() -> dict:
        with _lock:
            out = proj().status()
        out["blast_available"] = search_mod.blast_available()
        return out

    @app.post("/api/settings")
    def settings(body: dict = Body(default={})) -> dict:
        """Change the QC options that govern metrics and reference evaluation."""
        body = body or {}
        with _lock:
            p = proj()
            if "min_contig" in body:
                # A negative length is meaningless rather than wrong, so clamp
                # it; only a non-number is rejected.
                p.min_contig = max(0, _num(body, "min_contig", 0, int))
            if "circular_references" in body:
                p.circular_references = bool(body["circular_references"])
            if "primary_only" in body:
                p.primary_only = bool(body["primary_only"])
            if "genome_size" in body:
                value = body["genome_size"]
                p.genome_size = _num(body, "genome_size", None, int, 1) if value else None
            return _qc_settings(p)

    @app.post("/api/load")
    def load(body: dict = Body(...)) -> dict:
        path = (body or {}).get("path")
        if not path:
            fail("no path given")
        fmt = (body or {}).get("format") or None
        with _lock:
            proj().load(str(path), fmt)
            out = proj().status()
        out["blast_available"] = search_mod.blast_available()
        return out

    @app.post("/api/load-paths")
    def load_paths(body: dict = Body(...)) -> dict:
        path = (body or {}).get("path")
        if not path:
            fail("no path given")
        with _lock:
            return {"added": proj().attach_paths(str(path))}

    # ----------------------------------------------------------------- graph

    @app.get("/api/graph")
    def graph(
        min_length: int = Query(0, ge=0),
        max_nodes: int = Query(15000, ge=1, le=500000),
        component: int | None = Query(None),
    ) -> dict:
        with _lock:
            p = proj()
            g = p.require_graph()
            components = g.component_map()

            names = [
                name
                for name, segment in g.segments.items()
                if segment.length >= min_length
                and (component is None or components.get(name, 0) == component)
            ]
            total = len(names)
            truncated = False
            if len(names) > max_nodes:
                # Keep the longest segments -- they are the ones worth looking at.
                names = sorted(
                    names, key=lambda n: g.segments[n].length, reverse=True
                )[:max_nodes]
                truncated = True
            keep = set(names)

            segments = []
            for name in keep:
                segment = g.segments[name]
                left, right = g.degree(name)
                segments.append(
                    {
                        "name": name,
                        "length": segment.length,
                        "depth": segment.depth,
                        "gc": segment.gc,
                        "component": components.get(name, 0),
                        "deg_start": left,
                        "deg_end": right,
                        "circular": g.is_circular(name),
                        "ref_hits": segment.ref_hits,
                    }
                )
            links = [
                link.to_dict()
                for link in g.links.values()
                if link.from_name in keep and link.to_name in keep
            ]
            paths = [
                p_.to_dict()
                for p_ in g.paths.values()
                if all(step in keep for step, _ in p_.steps)
            ]
            references = sorted(p.reference_lengths)

        return {
            "segments": segments,
            "links": links,
            "paths": paths,
            "truncated": truncated,
            "shown": len(segments),
            "total": total,
            "references": references,
        }

    @app.get("/api/segment/{name}")
    def segment(name: str) -> dict:
        with _lock:
            g = proj().require_graph()
            seg = g.segments.get(name)
            if seg is None:
                fail(f"no segment named {name!r}", 404)
            left, right = g.degree(name)
            return {
                "name": name,
                "length": seg.length,
                "depth": seg.depth,
                "gc": seg.gc,
                "circular": g.is_circular(name),
                "deg_start": left,
                "deg_end": right,
                "tags": {k: v for k, v in seg.tags.items()},
                "ref_hits": seg.ref_hits,
                "neighbours": sorted(g.neighbours(name)),
                "sequence": seg.sequence or "",
            }

    # ------------------------------------------------------------- reference

    @app.post("/api/reference")
    def set_reference(body: dict = Body(...)) -> dict:
        body = body or {}
        path = body.get("path")
        if not path:
            fail("no reference path given")
        if not os.path.exists(str(path)):
            fail(f"reference not found: {path}")
        with _lock:
            p = proj()
            # QC options are sticky: set once here, they also govern later
            # re-alignments and the metrics endpoint.
            if "min_contig" in body:
                p.min_contig = max(0, int(body.get("min_contig") or 0))
            if "circular_references" in body:
                p.circular_references = bool(body["circular_references"])
            if "primary_only" in body:
                p.primary_only = bool(body["primary_only"])
            preset = body.get("preset") or align_mod.DEFAULT_PRESET
            if preset not in align_mod.PRESETS:
                fail(
                    f"unknown preset {preset!r}; choose one of "
                    f"{', '.join(align_mod.PRESETS)}"
                )
            report = p.set_reference(
                str(path),
                preset=preset,
                min_identity=_num(body, "min_identity", 0.0, float, 0.0, 1.0),
                min_length=_num(body, "min_length", 200, int, 0),
                threads=_num(body, "threads", app.state.threads, int, 1, 512),
            )
            return {
                "status": "aligned",
                "report": report.to_dict(),
                "metrics": p.metrics().to_dict(),
            }

    @app.delete("/api/reference")
    def clear_reference() -> dict:
        with _lock:
            proj().clear_reference()
            return proj().status()

    @app.get("/api/report")
    def report() -> dict:
        with _lock:
            p = proj()
            p.require_graph()
            return {
                "metrics": p.metrics().to_dict(),
                "reference": p.reference_report.to_dict() if p.reference_report else None,
            }

    # ------------------------------------------------------------ operations

    @app.post("/api/op")
    def operation(body: dict = Body(...)) -> dict:
        body = body or {}
        op = body.get("op")
        if not op:
            fail("no operation given")
        args = body.get("args") or {}
        if not isinstance(args, dict):
            fail("'args' must be an object")
        with _lock:
            p = proj()
            try:
                record = p.apply_operation(str(op), args)
            except (KeyError, TypeError) as exc:
                fail(f"operation {op!r} is missing or misuses an argument: {exc}")
            except ValueError as exc:
                fail(f"operation {op!r} got a bad argument: {exc}")
            return {
                "applied": op,
                "count": record.get("count", 0),
                "graph": p.require_graph().summary(),
                "undo_depth": len(p.undo_stack),
            }

    @app.post("/api/undo")
    def undo() -> dict:
        with _lock:
            p = proj()
            record = p.undo()
            return {
                "applied": f"undo:{record.get('kind', '')}",
                "count": record.get("count", 0),
                "graph": p.require_graph().summary(),
                "undo_depth": len(p.undo_stack),
            }

    # ------------------------------------------------------------- scaffolds

    @app.post("/api/scaffold")
    def make_plan(body: dict = Body(default={})) -> dict:
        body = body or {}
        with _lock:
            p = proj()
            method = body.get("method") or "reference"
            if method not in ("reference", "graph"):
                fail(f"unknown method {method!r}; choose 'reference' or 'graph'")
            plan = p.build_plan(
                method=method,
                min_identity=_num(body, "min_identity", 0.80, float, 0.0, 1.0),
                min_query_coverage=_num(body, "min_query_coverage", 0.30, float, 0.0, 1.0),
                min_align_length=_num(body, "min_align_length", 500, int, 0),
                min_gap=_num(body, "min_gap", 100, int, 0),
                fill_gaps_from_graph=bool(body.get("fill_gaps_from_graph", True)),
                include_unplaced=bool(body.get("include_unplaced", True)),
                break_misassemblies_first=bool(body.get("break_misassemblies_first", False)),
                threads=app.state.threads,
            )
            return {"plan": plan.to_dict(), "preview": p.preview()}

    @app.get("/api/scaffold")
    def get_plan() -> dict:
        with _lock:
            p = proj()
            if p.plan is None:
                return {"plan": None, "preview": None}
            return {"plan": p.plan.to_dict(), "preview": p.preview()}

    @app.put("/api/scaffold")
    def put_plan(body: dict = Body(...)) -> dict:
        body = body or {}
        data = body.get("plan", body)
        if not isinstance(data, dict) or "scaffolds" not in data:
            fail("body must be a scaffold plan with a 'scaffolds' list")
        if not isinstance(data["scaffolds"], list):
            fail("'scaffolds' must be a list")
        with _lock:
            p = proj()
            try:
                plan = ScaffoldPlan.from_dict(data)
            except (KeyError, TypeError, ValueError) as exc:
                fail(f"malformed scaffold plan: {exc}")
            p.set_plan(plan)
            return {"plan": p.plan.to_dict(), "preview": p.preview()}

    # ---------------------------------------------------------------- search

    @app.post("/api/search")
    def search(body: dict = Body(...)) -> dict:
        body = body or {}
        query = body.get("query")
        if not query:
            fail("no query given")
        with _lock:
            g = proj().require_graph()
            backend, hits = search_mod.search_graph(
                g,
                str(query),
                min_identity=_num(body, "min_identity", 0.8, float, 0.0, 1.0),
                min_length=_num(body, "min_length", 0, int, 0),
                threads=app.state.threads,
            )
        return {"backend": backend, "hits": [h.to_dict() for h in hits]}

    # ---------------------------------------------------------------- export

    @app.post("/api/export")
    def export(body: dict = Body(default={})) -> dict:
        body = body or {}
        outdir = str(body.get("outdir") or "plastr_out")
        what = body.get("what") or ["scaffolds", "agp", "gfa", "csv", "report"]
        if not isinstance(what, list) or not all(isinstance(k, str) for k in what):
            fail("'what' must be a list of artefact names")
        with _lock:
            try:
                written = proj().export(
                    outdir, list(what), overwrite=bool(body.get("overwrite", False))
                )
            except PermissionError:
                fail(f"permission denied writing to {outdir}", 403)
            except OSError as exc:
                fail(f"cannot write to {outdir}: {exc.strerror or exc}")
        return {"written": written}

    @app.get("/api/download/{kind}")
    def download(kind: str):
        with _lock:
            p = proj()
            p.require_graph()
            if kind == "gfa":
                data = "".join(gfa_mod.gfa_records(p.require_graph()))
                return _text_download(data, "graph.gfa")
            if kind == "session":
                return _text_download(
                    json.dumps(p.session_dict(), indent=2), "session.json"
                )
            if kind in ("scaffolds", "agp"):
                if p.plan is None:
                    fail("build a scaffold plan first")
                built = p.build()
                if kind == "scaffolds":
                    from ..core.io.fasta import fasta_string

                    return _text_download(
                        fasta_string(built.records), "scaffolds.fasta"
                    )
                from ..core.io.agp import write_agp

                buf = io.StringIO()
                write_agp(buf, built.agp_rows, comments=built.warnings)
                return _text_download(buf.getvalue(), "scaffolds.agp")
            if kind == "csv":
                # Built in memory rather than via a temp file, which previously
                # leaked one copy of the table into /tmp per download.
                buf = io.StringIO()
                p.write_csv(buf)
                return _text_download(buf.getvalue(), "segments.csv", "text/csv")
            if kind == "report":
                from ..core.report import build_report

                html = build_report(
                    p.metrics(),
                    reference_report=p.reference_report,
                    plan=p.plan,
                    built=p.build() if p.plan else None,
                    source_path=p.source_path,
                    reference_path=p.reference_path,
                    graph=p.graph,
                )
                return _text_download(html, "report.html", "text/html")
        fail(f"unknown download kind {kind!r}", 404)

    # --------------------------------------------------------------- session

    @app.get("/api/session")
    def get_session() -> dict:
        with _lock:
            return proj().session_dict()

    @app.post("/api/session")
    def post_session(body: dict = Body(...)) -> dict:
        body = body or {}
        with _lock:
            p = proj()
            if "layout" in body:
                p.layout = body["layout"] or {}
            if "settings" in body:
                p.settings = body["settings"] or {}
            if body.get("plan"):
                p.set_plan(ScaffoldPlan.from_dict(body["plan"]))
            if body.get("restore"):
                p.load_session(body["restore"], reload_files=True)
            return {"saved": True}

    # ---------------------------------------------------------------- browse

    @app.get("/api/browse")
    def browse(path: str = Query(default="")) -> dict:
        target = Path(path).expanduser() if path else Path.cwd()
        try:
            target = target.resolve()
            if target.is_file():
                target = target.parent
            if not target.is_dir():
                fail(f"not a directory: {target}", 404)
        except PermissionError:
            fail(f"permission denied: {target}", 403)
        except OSError as exc:
            fail(f"cannot read {path}: {exc.strerror or exc}")

        interesting = {
            ".gfa", ".gfa1", ".gfa2", ".fastg", ".fa", ".fasta", ".fna",
            ".ffn", ".gz", ".paths", ".json", ".agp",
        }
        entries: list[dict[str, Any]] = []
        try:
            for child in sorted(
                target.iterdir(), key=lambda c: (not c.is_dir(), c.name.lower())
            ):
                if child.name.startswith("."):
                    continue
                is_dir = child.is_dir()
                if not is_dir and child.suffix.lower() not in interesting:
                    continue
                try:
                    st = child.stat()
                except OSError:
                    continue
                # Never offer a FIFO or device: opening one blocks forever, and
                # the file picker is the one place a user can pick blind.
                if not is_dir and not stat.S_ISREG(st.st_mode):
                    continue
                size = st.st_size if not is_dir else 0
                entries.append({"name": child.name, "is_dir": is_dir, "size": size})
        except PermissionError:
            fail(f"permission denied: {target}", 403)

        return {
            "path": str(target),
            "parent": str(target.parent) if target.parent != target else None,
            "entries": entries,
        }

    # ------------------------------------------------------------ static UI

    if WEB_DIR.is_dir():
        app.mount("/", _NoCacheStatic(directory=str(WEB_DIR), html=True), name="web")

    return app


class _NoCacheStatic(StaticFiles):
    """Serve the UI with revalidation forced.

    The interface is a handful of local files, so caching saves nothing, and a
    stale cached ``index.html`` or ES module after an upgrade produces a broken
    page that a normal reload does not fix.
    """

    def file_response(self, *args, **kwargs):  # type: ignore[override]
        response = super().file_response(*args, **kwargs)
        response.headers["Cache-Control"] = "no-cache, must-revalidate"
        return response


def _qc_settings(project: Project) -> dict:
    return {
        "min_contig": project.min_contig,
        "circular_references": project.circular_references,
        "primary_only": project.primary_only,
        "genome_size": project.genome_size,
    }


def _text_download(text: str, filename: str, media_type: str = "text/plain"):
    data = text.encode()
    return StreamingResponse(
        io.BytesIO(data),
        media_type=media_type,
        headers={
            "Content-Disposition": f'attachment; filename="{filename}"',
            "Content-Length": str(len(data)),
        },
    )
