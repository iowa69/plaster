"""SPAdes FASTG reading.

FASTG headers from SPAdes look like::

    >EDGE_1_length_500_cov_10.5:EDGE_2_length_300_cov_9.1',EDGE_3_length_20_cov_4;

The name before the colon is this edge; the comma-separated list after it is the
set of edges it connects to. A trailing apostrophe means the reverse-complement
copy. SPAdes emits both strands as separate records, so we fold them back into
one segment per edge and derive orientation from the apostrophes.
"""

from __future__ import annotations

import os
import re

from ..errors import AssemblageFormatError
from ..model import AssemblyGraph, Link, Segment, depth_from_name
from .fasta import read_fasta

_EDGE = re.compile(r"^(.*?)(')?$")


def _split_name(token: str) -> tuple[str, str]:
    """Return (base name, orientation) for a FASTG edge token."""
    token = token.strip().rstrip(";")
    if token.endswith("'"):
        return token[:-1], "-"
    return token, "+"


def read_fastg(path: str | os.PathLike[str]) -> AssemblyGraph:
    graph = AssemblyGraph(name=os.path.basename(str(path)))
    graph.source_path = str(path)
    graph.source_format = "fastg"

    # First pass: collect sequences for the forward copy of every edge.
    records: list[tuple[str, str, str]] = []  # (base, orient, sequence)
    edges: list[tuple[tuple[str, str], tuple[str, str]]] = []
    saw_colon = False

    for header, description, sequence in read_fasta(path):
        full = header if not description else f"{header} {description}"
        full = full.strip().rstrip(";")
        if ":" in full:
            saw_colon = True
            left, right = full.split(":", 1)
            neighbours = [t for t in right.split(",") if t.strip()]
        else:
            left, neighbours = full, []
        base, orient = _split_name(left)
        if not base:
            raise AssemblageFormatError("FASTG record with an empty edge name", str(path))
        records.append((base, orient, sequence.upper()))
        for token in neighbours:
            nb_base, nb_orient = _split_name(token)
            if nb_base:
                edges.append(((base, orient), (nb_base, nb_orient)))

    if not records:
        raise AssemblageFormatError("no FASTG records found", str(path))
    if not saw_colon:
        raise AssemblageFormatError(
            "no ':' found in any header -- this looks like plain FASTA, not FASTG", str(path)
        )

    for base, orient, sequence in records:
        if base in graph.segments:
            continue
        # Only the '+' record's sequence is stored; '-' records are its revcomp.
        if orient == "+":
            graph.add_segment(
                Segment(base, sequence, len(sequence), depth_from_name(base)), replace=True
            )
    # Some files only carry the reverse copy of an edge; fall back to it.
    for base, orient, sequence in records:
        if base not in graph.segments:
            from ..sequence import revcomp

            graph.add_segment(
                Segment(base, revcomp(sequence), len(sequence), depth_from_name(base)),
                replace=True,
            )

    for (a_base, a_or), (b_base, b_or) in edges:
        graph.add_link(Link(a_base, a_or, b_base, b_or, 0, "*"))

    # FASTG has no field for the overlap, but SPAdes edges really do overlap by
    # k-1 bases. Measure it, or every merge through a link duplicates them.
    from .overlaps import apply_inferred_overlaps

    apply_inferred_overlaps(graph)
    return graph
