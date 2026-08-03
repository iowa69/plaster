"""One entry point for loading any supported assembly file."""

from __future__ import annotations

import os

from ..errors import AssemblageFormatError
from ..model import AssemblyGraph, Path, Segment, depth_from_name
from . import fastg as fastg_mod
from . import gfa as gfa_mod
from .fasta import read_fasta

SUPPORTED = ("gfa", "gfa2", "fastg", "fasta")


def read_contigs_fasta(path: str | os.PathLike[str]) -> AssemblyGraph:
    """Load a plain contig FASTA as an edge-free graph.

    Still useful: QC, reference alignment, and scaffolding all work without links,
    they just cannot fill gaps from graph paths.
    """
    graph = AssemblyGraph(name=os.path.basename(str(path)))
    graph.source_path = str(path)
    graph.source_format = "fasta"
    for name, desc, seq in read_fasta(path):
        depth = depth_from_name(name)
        tags: dict[str, object] = {}
        if desc:
            tags["desc"] = desc
        graph.add_segment(
            Segment(name, seq.upper(), len(seq), depth, tags), replace=True
        )
    if not graph.segments:
        raise AssemblageFormatError("FASTA contained no sequences", str(path))
    return graph


def load_graph(path: str | os.PathLike[str], fmt: str | None = None) -> AssemblyGraph:
    """Load ``path``, sniffing the format unless ``fmt`` is given."""
    if not os.path.exists(path):
        raise AssemblageFormatError(f"file not found: {path}")
    detected = fmt or gfa_mod.detect_format(path)
    if detected == "gfa":
        return gfa_mod.read_gfa(path)
    if detected == "gfa2":
        return gfa_mod.read_gfa2(path)
    if detected == "fastg":
        return fastg_mod.read_fastg(path)
    if detected == "fasta":
        return read_contigs_fasta(path)
    raise AssemblageFormatError(
        f"could not determine the format of {path}. "
        f"Supported: {', '.join(SUPPORTED)}. Pass --format to override."
    )


def load_spades_paths(graph: AssemblyGraph, path: str | os.PathLike[str]) -> int:
    """Attach a SPAdes ``contigs.paths`` file to a graph. Returns paths added.

    The file alternates a name line with a path line; a name ending in ``'``
    denotes the reverse-complement copy, which we skip since it is redundant.
    """
    added = 0
    with open(path) as handle:
        lines = [ln.strip() for ln in handle if ln.strip()]
    i = 0
    while i < len(lines) - 1:
        name = lines[i]
        # A path may be split across lines ending with ';'
        steps_text = lines[i + 1]
        i += 2
        while steps_text.endswith(";") and i < len(lines):
            steps_text = steps_text[:-1] + lines[i]
            i += 1
        if name.endswith("'"):
            continue
        steps: list[tuple[str, str]] = []
        for token in steps_text.replace(";", "").split(","):
            token = token.strip()
            if len(token) < 2 or token[-1] not in "+-":
                continue
            seg, orient = token[:-1], token[-1]
            if seg in graph.segments:
                steps.append((seg, orient))
        if steps:
            graph.add_path(Path(name, steps))
            added += 1
    return added
