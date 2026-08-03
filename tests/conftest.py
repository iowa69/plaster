"""Fixtures shared by the Plastr test suite.

Two families of fixture live here.

``tiny_graph`` is a hand-built :class:`AssemblyGraph` whose every interesting
property is written down in :data:`TINY` below, so tests assert exact numbers
instead of "is not None".

The ``demo_*`` family is the ground-truth dataset produced by
``examples/make_demo_data.py``. The generator is seeded, so every number derived
from it is reproducible; it is run once per session into a tmp directory.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from plastr.core.analysis.align import _copy_alignment, alignment_backend
from plastr.core.model import AssemblyGraph, Link, Segment

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEMO_SCRIPT = PROJECT_ROOT / "examples" / "make_demo_data.py"

#: ``"mappy"``, ``"minimap2"`` or ``"none"``. Modules that need real alignments
#: skip themselves when this is ``"none"``.
ALIGNMENT_BACKEND = alignment_backend()

HAVE_ALIGNER = ALIGNMENT_BACKEND != "none"

requires_aligner = pytest.mark.skipif(
    not HAVE_ALIGNER,
    reason="neither the mappy module nor the minimap2 binary is available",
)


# ---------------------------------------------------------------------------
# The hand-built tiny graph
# ---------------------------------------------------------------------------

# Overlap sequences are spelled out so walk_sequence's trimming can be asserted
# against an exact string rather than only a length.
OV_AB = "GATTACAGAT"  # 10 bp shared by A's 3' end and B's 5' end
OV_BC = "TTCCGGAACC"  # 10 bp shared by B's 3' end and C's 5' end

SEQ_A = "A" * 90 + OV_AB  # 100 bp
SEQ_B = OV_AB + "C" * 60 + OV_BC  # 80 bp
SEQ_C = OV_BC + "G" * 50  # 60 bp
SEQ_D = "T" * 40  # 40 bp, isolated
SEQ_R = "ACGT" * 7 + "AC"  # 30 bp, self-looping (circular)


class TINY:
    """Everything true about :func:`tiny_graph`, as literals.

    Layout::

        A(100) --10bp--> B(80) --10bp--> C(60)     component 0, 240 bp
        D(40)                                      component 1, 40 bp, isolated
        R(30) ---self-loop--->                     component 2, 30 bp, circular
    """

    names = ("A", "B", "C", "D", "R")
    lengths = {"A": 100, "B": 80, "C": 60, "D": 40, "R": 30}
    depths = {"A": 10.0, "B": 20.0, "C": 30.0, "D": 5.0, "R": 40.0}
    total_length = 310
    link_count = 3  # A->B, B->C, R self-loop
    components = [["A", "B", "C"], ["D"], ["R"]]
    component_lengths = [240, 40, 30]
    dead_end_names = {"A", "C", "D"}
    dead_end_count = 4  # A one end, C one end, D both ends
    #: The exact result of walking A+ -> B+ -> C+ with the 10 bp overlaps trimmed.
    walk_abc = "A" * 90 + OV_AB + "C" * 60 + OV_BC + "G" * 50
    walk_abc_length = 220  # 100 + (80-10) + (60-10)


@pytest.fixture
def tiny_graph() -> AssemblyGraph:
    """A small graph with fully known properties. See :class:`TINY`."""
    graph = AssemblyGraph(name="tiny")
    for name, seq in (
        ("A", SEQ_A),
        ("B", SEQ_B),
        ("C", SEQ_C),
        ("D", SEQ_D),
        ("R", SEQ_R),
    ):
        graph.add_segment(Segment(name=name, sequence=seq, depth=TINY.depths[name]))
    assert graph.add_link(Link("A", "+", "B", "+", 10, "10M"))
    assert graph.add_link(Link("B", "+", "C", "+", 10, "10M"))
    assert graph.add_link(Link("R", "+", "R", "+", 0, "*"))
    return graph


@pytest.fixture
def chain_graph() -> AssemblyGraph:
    """A plain unbranching chain P+ -> A+ -> B+ -> C+ -> Q+ with no overlaps."""
    graph = AssemblyGraph(name="chain")
    pieces = {
        "P": "AAAAAAAAAA",
        "A": "CCCCCCCCCC",
        "B": "GGGGGGGGGG",
        "C": "TTTTTTTTTT",
        "Q": "ACACACACAC",
    }
    for name, seq in pieces.items():
        graph.add_segment(Segment(name=name, sequence=seq))
    for left, right in (("P", "A"), ("A", "B"), ("B", "C"), ("C", "Q")):
        assert graph.add_link(Link(left, "+", right, "+", 0, "*"))
    return graph


# ---------------------------------------------------------------------------
# The generated demo dataset
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def demo_dir(tmp_path_factory) -> Path:
    """Run ``examples/make_demo_data.py`` once into a session tmp directory."""
    out = tmp_path_factory.mktemp("plastr_demo")
    proc = subprocess.run(
        [sys.executable, str(DEMO_SCRIPT), str(out)],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
    )
    if proc.returncode != 0:
        pytest.fail(
            "make_demo_data.py failed:\n"
            f"stdout: {proc.stdout}\nstderr: {proc.stderr}"
        )
    return out


@pytest.fixture(scope="session")
def demo_paths(demo_dir: Path) -> dict[str, Path]:
    paths = {
        "dir": demo_dir,
        "gfa": demo_dir / "assembly.gfa",
        "reference": demo_dir / "reference.fasta",
        "contigs": demo_dir / "contigs.fasta",
        "truth": demo_dir / "truth.json",
    }
    for key, path in paths.items():
        assert path.exists(), f"demo data is missing {key}: {path}"
    return paths


@pytest.fixture(scope="session")
def demo_truth(demo_paths) -> dict:
    with open(demo_paths["truth"]) as handle:
        return json.load(handle)["expected"]


@pytest.fixture
def demo_graph(demo_paths) -> AssemblyGraph:
    """A freshly loaded demo graph -- function scoped because tests mutate it."""
    from plastr.core.io.loader import load_graph

    return load_graph(str(demo_paths["gfa"]))


@pytest.fixture(scope="session")
def demo_reference_lengths(demo_paths) -> dict[str, int]:
    from plastr.core.analysis.align import reference_lengths

    return reference_lengths(str(demo_paths["reference"]))


@pytest.fixture(scope="session")
def _demo_alignments_session(demo_paths):
    """Align the demo contigs once; every test gets copies of these."""
    if not HAVE_ALIGNER:
        pytest.skip("no alignment backend")
    from plastr.core.analysis.align import align_graph
    from plastr.core.io.loader import load_graph

    graph = load_graph(str(demo_paths["gfa"]))
    return align_graph(graph, str(demo_paths["reference"]))


@pytest.fixture
def demo_alignments(_demo_alignments_session):
    """Private copies of the demo alignments, safe to mutate."""
    return [_copy_alignment(a) for a in _demo_alignments_session]


@pytest.fixture(scope="session")
def demo_report(_demo_alignments_session, demo_reference_lengths, demo_paths):
    from plastr.core.analysis.misassembly import evaluate_against_reference
    from plastr.core.io.loader import load_graph

    graph = load_graph(str(demo_paths["gfa"]))
    contig_lengths = {n: s.length for n, s in graph.segments.items()}
    return evaluate_against_reference(
        _demo_alignments_session, demo_reference_lengths, contig_lengths
    )
