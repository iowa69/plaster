"""Integration tests for the HTTP API.

These check the contract in docs/API.md: the shapes the browser depends on, and
that errors come back as readable messages rather than stack traces.
"""

from __future__ import annotations

import json

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from plastr.core.project import Project  # noqa: E402

from conftest import requires_aligner  # noqa: E402
from plastr.server.app import create_app  # noqa: E402


@pytest.fixture
def client(demo_paths):
    app = create_app(Project(), threads=2)
    with TestClient(app) as c:
        c.post("/api/load", json={"path": str(demo_paths["gfa"])})
        yield c


@pytest.fixture
def empty_client():
    app = create_app(Project(), threads=2)
    with TestClient(app) as c:
        yield c


# --------------------------------------------------------------------- status


def test_status_before_loading(empty_client):
    body = empty_client.get("/api/status").json()
    assert body["loaded"] is False
    assert body["graph"] is None
    assert "align_backend" in body
    assert "blast_available" in body


def test_load_reports_graph_summary(client):
    body = client.get("/api/status").json()
    assert body["loaded"] is True
    assert body["graph"]["segments"] == 10
    assert body["graph"]["links"] == 9
    assert body["graph"]["has_sequences"] is True
    assert body["source_format"] == "gfa"


def test_load_missing_file_is_a_readable_400(empty_client):
    r = empty_client.post("/api/load", json={"path": "/definitely/not/here.gfa"})
    assert r.status_code == 400
    assert "not found" in r.json()["error"]


def test_operations_before_load_do_not_crash(empty_client):
    r = empty_client.get("/api/graph")
    assert r.status_code == 400
    assert "no assembly" in r.json()["error"].lower()


# ---------------------------------------------------------------------- graph


def test_graph_payload_shape(client):
    body = client.get("/api/graph").json()
    assert body["shown"] == body["total"] == 10
    assert body["truncated"] is False
    segment = body["segments"][0]
    for key in ("name", "length", "depth", "gc", "component", "deg_start", "deg_end", "circular"):
        assert key in segment
    link = body["links"][0]
    assert set(link) == {"from", "from_orient", "to", "to_orient", "overlap"}
    assert link["from_orient"] in "+-"


def test_graph_min_length_filter(client):
    body = client.get("/api/graph", params={"min_length": 10_000}).json()
    assert body["shown"] < 10
    assert all(s["length"] >= 10_000 for s in body["segments"])


def test_graph_truncation_keeps_longest(client):
    body = client.get("/api/graph", params={"max_nodes": 3}).json()
    assert body["truncated"] is True
    assert body["shown"] == 3
    assert body["total"] == 10
    lengths = [s["length"] for s in body["segments"]]
    assert min(lengths) >= 11_000  # the three longest in the demo set


def test_links_never_reference_a_hidden_segment(client):
    body = client.get("/api/graph", params={"max_nodes": 4}).json()
    shown = {s["name"] for s in body["segments"]}
    for link in body["links"]:
        assert link["from"] in shown and link["to"] in shown


def test_segment_detail_includes_sequence(client):
    body = client.get("/api/segment/ctg_clean_1").json()
    assert body["length"] == 19_000
    assert len(body["sequence"]) == 19_000
    assert body["neighbours"] == ["ctg_repeat"]


def test_unknown_segment_is_404(client):
    assert client.get("/api/segment/nope").status_code == 404


# ------------------------------------------------------------------ reference


@requires_aligner
def test_reference_alignment_and_report(client, demo_paths):
    r = client.post(
        "/api/reference",
        json={"path": str(demo_paths["reference"]), "preset": "asm5"},
    )
    assert r.status_code == 200
    report = r.json()["report"]
    assert 78 < report["genome_fraction"] < 81
    assert report["num_relocations"] == 1
    assert report["num_inversions"] == 1
    assert report["num_translocations"] == 1

    graph = client.get("/api/graph").json()
    assert graph["references"] == ["chromosome", "plasmid"]
    assert any(s["ref_hits"] for s in graph["segments"])

    client.delete("/api/reference")
    graph = client.get("/api/graph").json()
    assert all(not s["ref_hits"] for s in graph["segments"])


def test_missing_reference_file_is_400(client):
    r = client.post("/api/reference", json={"path": "/nope.fasta"})
    assert r.status_code == 400


# ----------------------------------------------------------------- operations


def test_delete_and_undo_restores_exactly(client):
    before = client.get("/api/status").json()["graph"]
    r = client.post("/api/op", json={"op": "delete", "args": {"names": ["ctg_repeat"]}})
    assert r.status_code == 200
    assert r.json()["undo_depth"] == 1
    assert r.json()["graph"]["segments"] == before["segments"] - 1

    r = client.post("/api/undo")
    assert r.json()["graph"] == before
    assert r.json()["undo_depth"] == 0


def test_undo_with_empty_stack_is_400(client):
    assert client.post("/api/undo").status_code == 400


def test_unknown_operation_is_400(client):
    r = client.post("/api/op", json={"op": "nonsense"})
    assert r.status_code == 400
    assert "unknown operation" in r.json()["error"]


def test_break_misassemblies_needs_a_reference(client):
    r = client.post("/api/op", json={"op": "break_misassemblies", "args": {}})
    assert r.status_code == 400
    assert "reference" in r.json()["error"].lower()


# ------------------------------------------------------------------ scaffolds


def test_scaffold_requires_a_reference(client):
    r = client.post("/api/scaffold", json={"method": "reference"})
    assert r.status_code == 400


def test_graph_method_scaffolding_needs_no_reference(client):
    r = client.post("/api/scaffold", json={"method": "graph"})
    assert r.status_code == 200
    plan = r.json()["plan"]
    assert plan["method"] == "graph"
    assert plan["scaffold_count"] >= 1
    assert r.json()["preview"]["total_length"] > 0


@requires_aligner
def test_reference_scaffold_round_trip_and_manual_edit(client, demo_paths):
    client.post("/api/reference", json={"path": str(demo_paths["reference"]), "preset": "asm5"})
    r = client.post("/api/scaffold", json={"method": "reference", "include_unplaced": True})
    assert r.status_code == 200
    plan = r.json()["plan"]
    preview = r.json()["preview"]
    assert plan["scaffold_count"] >= 2
    assert preview["warnings"] == []

    # GET returns what POST built.
    assert client.get("/api/scaffold").json()["plan"] == plan

    # A manual edit: drop the last member of the first scaffold.
    target = next(s for s in plan["scaffolds"] if len(s["members"]) > 1)
    removed = target["members"].pop()
    r = client.put("/api/scaffold", json={"plan": plan})
    assert r.status_code == 200
    updated = r.json()["plan"]
    edited = next(s for s in updated["scaffolds"] if s["name"] == target["name"])
    assert removed["segment"] not in [m["segment"] for m in edited["members"]]
    assert r.json()["preview"]["total_length"] < preview["total_length"]


def test_put_rejects_a_body_that_is_not_a_plan(client):
    assert client.put("/api/scaffold", json={"plan": {"nope": 1}}).status_code == 400


# --------------------------------------------------------------------- export


def test_downloads_have_content_and_filenames(client):
    for kind, needle in (("gfa", "S\t"), ("csv", "name,length")):
        r = client.get(f"/api/download/{kind}")
        assert r.status_code == 200, kind
        assert "attachment" in r.headers["content-disposition"]
        assert needle in r.text


def test_scaffold_download_requires_a_plan(client):
    assert client.get("/api/download/scaffolds").status_code == 400


def test_unknown_download_kind_is_404(client):
    assert client.get("/api/download/banana").status_code == 404


def test_export_writes_files(client, tmp_path):
    client.post("/api/scaffold", json={"method": "graph"})
    r = client.post(
        "/api/export",
        json={"outdir": str(tmp_path / "out"), "what": ["scaffolds", "agp", "gfa", "csv"]},
    )
    assert r.status_code == 200
    written = {item["kind"]: item for item in r.json()["written"]}
    assert set(written) == {"scaffolds", "agp", "gfa", "csv"}
    for item in written.values():
        assert item["bytes"] > 0


# -------------------------------------------------------------------- session


def test_session_round_trip(client):
    r = client.post(
        "/api/session",
        json={"layout": {"zoom": 2.5}, "settings": {"colour": "depth"}},
    )
    assert r.json()["saved"] is True
    body = client.get("/api/session").json()
    assert body["layout"] == {"zoom": 2.5}
    assert body["settings"]["colour"] == "depth"
    assert json.dumps(body)  # must be serialisable for the download endpoint


# --------------------------------------------------------------------- browse


def test_browse_lists_assembly_files(client, demo_paths):
    body = client.get("/api/browse", params={"path": str(demo_paths["dir"])}).json()
    names = {e["name"] for e in body["entries"]}
    assert "assembly.gfa" in names
    assert "reference.fasta" in names
    assert body["parent"] is not None


def test_browse_rejects_a_missing_directory(client):
    assert client.get("/api/browse", params={"path": "/no/such/dir"}).status_code == 404


# --------------------------------------------------------------------- search


def test_search_finds_a_segment_by_its_own_sequence(client, demo_paths):
    sequence = client.get("/api/segment/ctg_clean_3").json()["sequence"][:1500]
    r = client.post("/api/search", json={"query": sequence, "min_identity": 0.95})
    if r.status_code == 503:
        pytest.skip("no blast or minimap2 available")
    assert r.status_code == 200
    hits = r.json()["hits"]
    assert any(h["segment"] == "ctg_clean_3" for h in hits)


def test_search_rejects_junk(client):
    r = client.post("/api/search", json={"query": "this is not a sequence!!"})
    assert r.status_code == 400


# -------------------------------------------------------------------- settings


def test_status_reports_the_qc_settings(client):
    settings = client.get("/api/status").json()["settings"]
    assert settings == {
        "min_contig": 0,
        "circular_references": True,
        "primary_only": True,
        "genome_size": None,
    }


def test_settings_can_be_changed_and_are_returned(client):
    r = client.post("/api/settings", json={"min_contig": 500, "primary_only": False})
    assert r.status_code == 200
    assert r.json()["min_contig"] == 500
    assert r.json()["primary_only"] is False
    # unlisted keys are left alone
    assert r.json()["circular_references"] is True
    assert client.get("/api/status").json()["settings"]["min_contig"] == 500


def test_min_contig_is_honoured_by_the_report(client):
    before = client.get("/api/report").json()["metrics"]["num_contigs"]
    client.post("/api/settings", json={"min_contig": 5_000})
    after = client.get("/api/report").json()["metrics"]["num_contigs"]
    assert after < before, "min_contig did not filter the statistics"


def test_settings_accepts_an_empty_body(client):
    r = client.post("/api/settings", json={})
    assert r.status_code == 200
    assert r.json()["min_contig"] == 0


def test_a_negative_min_contig_is_clamped(client):
    assert client.post("/api/settings", json={"min_contig": -10}).json()["min_contig"] == 0


def test_genome_size_can_be_cleared(client):
    client.post("/api/settings", json={"genome_size": 5_000_000})
    assert client.get("/api/status").json()["settings"]["genome_size"] == 5_000_000
    client.post("/api/settings", json={"genome_size": None})
    assert client.get("/api/status").json()["settings"]["genome_size"] is None


@requires_aligner
def test_settings_may_be_passed_inline_when_aligning(client, demo_paths):
    r = client.post(
        "/api/reference",
        json={
            "path": str(demo_paths["reference"]),
            "preset": "asm5",
            "circular_references": False,
            "min_contig": 500,
        },
    )
    assert r.status_code == 200
    settings = client.get("/api/status").json()["settings"]
    assert settings["circular_references"] is False
    assert settings["min_contig"] == 500


# --------------------------------------------------------- error contract


def test_every_error_carries_the_documented_error_key(client):
    """docs/API.md promises {"error", "detail"} and the UI reads `error`.

    FastAPI's default HTTPException body is {"detail": ...}, which made the UI
    show a bare "HTTP 400" and demote the real message.
    """
    cases = [
        client.get("/api/segment/nope"),
        client.get("/api/download/bogus"),
        client.post("/api/load", json={}),
        client.post("/api/op", json={}),
        client.post("/api/reference", json={"path": "/nope.fasta"}),
        client.put("/api/scaffold", json={"plan": {"nope": 1}}),
        client.post("/api/undo"),
    ]
    for r in cases:
        body = r.json()
        assert "error" in body, f"{r.url} returned {body}"
        assert body["error"], f"{r.url} returned an empty error"
        assert "detail" in body


def test_foreseeable_io_failures_are_readable_errors_not_500s(client):
    for r in (
        client.post("/api/load-paths", json={"path": "/nope.paths"}),
        client.post("/api/load", json={"path": "/tmp"}),
        client.post("/api/reference", json={"path": "/tmp"}),
    ):
        assert r.status_code == 400, r.text
        assert "error" in r.json()


# ----------------------------------------------------------- export safety


def test_export_refuses_to_overwrite_unless_asked(client, tmp_path):
    client.post("/api/scaffold", json={"method": "graph"})
    out = str(tmp_path / "out")

    first = client.post("/api/export", json={"outdir": out, "what": ["csv"]})
    assert first.status_code == 200

    again = client.post("/api/export", json={"outdir": out, "what": ["csv"]})
    assert again.status_code == 400
    assert "already contains" in again.json()["error"]

    forced = client.post(
        "/api/export", json={"outdir": out, "what": ["csv"], "overwrite": True}
    )
    assert forced.status_code == 200


def test_csv_download_leaves_no_temporary_file(client, tmp_path, monkeypatch):
    import tempfile

    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    before = set(tmp_path.iterdir())
    r = client.get("/api/download/csv")
    assert r.status_code == 200
    assert "name,length" in r.text
    assert set(tmp_path.iterdir()) == before, "download left a temp file behind"
