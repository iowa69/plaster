"""Assemblage -- interactive assembly graph studio.

Bandage-style visualisation, QUAST-style evaluation, and reference-guided
scaffolding that exports usable FASTA and AGP.

Typical library use::

    from assemblage import Project

    project = Project()
    project.load("assembly.gfa")
    project.set_reference("reference.fasta", preset="asm5")
    print(project.reference_report.genome_fraction)

    project.build_plan(method="reference")
    project.export("out/", ["scaffolds", "agp", "report"])
"""

from __future__ import annotations

__version__ = "1.0.0"

__all__ = [
    "Project",
    "AssemblyGraph",
    "Segment",
    "Link",
    "load_graph",
    "compute_metrics",
    "__version__",
]


def __getattr__(name: str):
    # Imported lazily so `import assemblage` stays cheap and does not require
    # the optional alignment dependencies to be present.
    if name == "Project":
        from .core.project import Project

        return Project
    if name in ("AssemblyGraph", "Segment", "Link"):
        from .core import model

        return getattr(model, name)
    if name == "load_graph":
        from .core.io.loader import load_graph

        return load_graph
    if name == "compute_metrics":
        from .core.analysis.metrics import compute_metrics

        return compute_metrics
    raise AttributeError(f"module 'assemblage' has no attribute {name!r}")
