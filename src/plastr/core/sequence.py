"""Sequence helpers.

Kept tiny and dependency-free so the graph model can use it without pulling in
numpy or biopython.
"""

from __future__ import annotations

_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn-",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn-",
)

_GC = frozenset("GCgcSs")
_ACGT = frozenset("ACGTacgt")


def revcomp(seq: str) -> str:
    """Reverse complement, IUPAC aware."""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """Fraction of G+C among unambiguous bases. Returns 0.0 for an empty/N-only sequence."""
    gc = 0
    total = 0
    for ch in seq:
        if ch in _ACGT:
            total += 1
            if ch in _GC:
                gc += 1
        elif ch in _GC:
            total += 1
            gc += 1
    return gc / total if total else 0.0


def n_count(seq: str) -> int:
    return seq.count("N") + seq.count("n")


def oriented(seq: str, orientation: str) -> str:
    """Return ``seq`` on the requested strand. ``orientation`` is '+' or '-'."""
    return seq if orientation == "+" else revcomp(seq)


def flip(orientation: str) -> str:
    return "-" if orientation == "+" else "+"
