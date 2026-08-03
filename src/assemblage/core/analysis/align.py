"""Aligning assembly segments to a reference.

Prefers the ``mappy`` module (minimap2's Python binding, no temp files); falls
back to the ``minimap2`` binary and PAF parsing when mappy is unavailable. Both
paths produce the same :class:`Alignment` records so nothing downstream cares
which was used.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Iterable, Sequence

from ..errors import MissingDependencyError
from ..io.fasta import read_fasta, write_fasta
from ..model import AssemblyGraph

# minimap2 presets by expected divergence between assembly and reference
PRESETS = {
    "asm5": "≤1% divergence (same strain)",
    "asm10": "≤5% divergence (same species)",
    "asm20": "≤10% divergence (related species)",
}
DEFAULT_PRESET = "asm10"

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(slots=True)
class Alignment:
    """One alignment block of a query segment against a reference sequence."""

    query: str
    q_len: int
    q_st: int
    q_en: int
    strand: int  # +1 forward, -1 reverse
    ref: str
    r_len: int
    r_st: int
    r_en: int
    matches: int = 0
    block_len: int = 0
    mapq: int = 0
    nm: int = 0
    is_primary: bool = True
    cigar: list[tuple[int, str]] = field(default_factory=list)

    @property
    def q_span(self) -> int:
        return self.q_en - self.q_st

    @property
    def r_span(self) -> int:
        return self.r_en - self.r_st

    @property
    def identity(self) -> float:
        return self.matches / self.block_len if self.block_len else 0.0

    @property
    def query_coverage(self) -> float:
        return self.q_span / self.q_len if self.q_len else 0.0

    def indel_counts(self) -> tuple[int, int, int]:
        """(insertion events, deletion events, total indel bases) from the CIGAR."""
        ins = dele = bases = 0
        for count, op in self.cigar:
            if op == "I":
                ins += 1
                bases += count
            elif op == "D":
                dele += 1
                bases += count
        return ins, dele, bases

    def mismatch_estimate(self) -> int:
        """Substitutions = edit distance minus indel bases."""
        _, _, indel_bases = self.indel_counts()
        return max(0, self.nm - indel_bases)

    def to_dict(self) -> dict:
        return {
            "query": self.query,
            "q_len": self.q_len,
            "q_st": self.q_st,
            "q_en": self.q_en,
            "strand": self.strand,
            "ref": self.ref,
            "r_len": self.r_len,
            "r_st": self.r_st,
            "r_en": self.r_en,
            "identity": round(self.identity, 5),
            "mapq": self.mapq,
            "is_primary": self.is_primary,
            "block_len": self.block_len,
        }


def parse_cigar(text: str) -> list[tuple[int, str]]:
    return [(int(n), op) for n, op in _CIGAR_RE.findall(text or "")]


def have_mappy() -> bool:
    try:
        import mappy  # noqa: F401
    except ImportError:
        return False
    return True


def have_minimap2() -> bool:
    return shutil.which("minimap2") is not None


def alignment_backend() -> str:
    if have_mappy():
        return "mappy"
    if have_minimap2():
        return "minimap2"
    return "none"


# ---------------------------------------------------------------------------
# mappy path
# ---------------------------------------------------------------------------


def _align_with_mappy(
    queries: Sequence[tuple[str, str]],
    reference_path: str,
    preset: str,
    min_length: int,
    threads: int,
) -> list[Alignment]:
    import mappy

    aligner = mappy.Aligner(str(reference_path), preset=preset, n_threads=threads)
    if not aligner:
        raise RuntimeError(f"minimap2 failed to index {reference_path}")

    out: list[Alignment] = []
    for name, seq in queries:
        if not seq or len(seq) < min_length:
            continue
        for hit in aligner.map(seq, cs=False, MD=False):
            out.append(
                Alignment(
                    query=name,
                    q_len=len(seq),
                    q_st=hit.q_st,
                    q_en=hit.q_en,
                    strand=hit.strand,
                    ref=hit.ctg,
                    r_len=hit.ctg_len,
                    r_st=hit.r_st,
                    r_en=hit.r_en,
                    matches=hit.mlen,
                    block_len=hit.blen,
                    mapq=hit.mapq,
                    nm=hit.NM,
                    is_primary=bool(hit.is_primary),
                    cigar=[(c, "MIDNSHP=X"[op]) for c, op in hit.cigar],
                )
            )
    return out


# ---------------------------------------------------------------------------
# minimap2 CLI path
# ---------------------------------------------------------------------------


def _parse_paf(lines: Iterable[str]) -> list[Alignment]:
    out: list[Alignment] = []
    for raw in lines:
        line = raw.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        tags = {}
        for extra in f[12:]:
            parts = extra.split(":", 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]
        out.append(
            Alignment(
                query=f[0],
                q_len=int(f[1]),
                q_st=int(f[2]),
                q_en=int(f[3]),
                strand=1 if f[4] == "+" else -1,
                ref=f[5],
                r_len=int(f[6]),
                r_st=int(f[7]),
                r_en=int(f[8]),
                matches=int(f[9]),
                block_len=int(f[10]),
                mapq=int(f[11]),
                nm=int(tags.get("NM", 0) or 0),
                is_primary=tags.get("tp", "P") == "P",
                cigar=parse_cigar(tags.get("cg", "")),
            )
        )
    return out


def _align_with_cli(
    queries: Sequence[tuple[str, str]],
    reference_path: str,
    preset: str,
    min_length: int,
    threads: int,
) -> list[Alignment]:
    with tempfile.TemporaryDirectory(prefix="assemblage_aln_") as tmp:
        query_path = os.path.join(tmp, "query.fa")
        write_fasta(
            query_path, [(n, s) for n, s in queries if s and len(s) >= min_length]
        )
        cmd = [
            "minimap2",
            "-c",
            "-x",
            preset,
            "-t",
            str(threads),
            "--secondary=yes",
            "-N",
            "5",
            str(reference_path),
            query_path,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {proc.stderr.strip()[:500]}")
        return _parse_paf(proc.stdout.splitlines())


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------


def align_sequences(
    queries: Sequence[tuple[str, str]],
    reference_path: str | os.PathLike[str],
    preset: str = DEFAULT_PRESET,
    min_length: int = 200,
    threads: int = 4,
    split_long_indels: int = 1000,
) -> list[Alignment]:
    """Align ``(name, sequence)`` pairs to a reference FASTA.

    Alignments spanning an indel of ``split_long_indels`` bp or more are broken
    into separate blocks; pass 0 to keep minimap2's raw blocks.
    """
    backend = alignment_backend()
    if backend == "none":
        raise MissingDependencyError(
            "minimap2",
            "Install it with 'conda install -c bioconda minimap2 mappy', "
            "or use the assemblage conda environment.",
        )
    if preset not in PRESETS:
        preset = DEFAULT_PRESET
    if backend == "mappy":
        try:
            hits = _align_with_mappy(
                queries, str(reference_path), preset, min_length, threads
            )
        except Exception:
            if not have_minimap2():
                raise
            hits = _align_with_cli(
                queries, str(reference_path), preset, min_length, threads
            )
    else:
        hits = _align_with_cli(queries, str(reference_path), preset, min_length, threads)

    if split_long_indels > 0:
        hits = split_all(hits, split_long_indels)
    return hits


def align_graph(
    graph: AssemblyGraph,
    reference_path: str | os.PathLike[str],
    preset: str = DEFAULT_PRESET,
    min_length: int = 200,
    min_identity: float = 0.0,
    threads: int = 4,
    split_long_indels: int = 1000,
) -> list[Alignment]:
    """Align every segment that has a sequence to the reference."""
    queries = [
        (name, seg.sequence)
        for name, seg in graph.segments.items()
        if seg.sequence
    ]
    hits = align_sequences(
        queries, reference_path, preset, min_length, threads, split_long_indels
    )
    if min_identity > 0:
        hits = [h for h in hits if h.identity >= min_identity]
    return hits


def reference_lengths(reference_path: str | os.PathLike[str]) -> dict[str, int]:
    return {name: len(seq) for name, _d, seq in read_fasta(reference_path)}


def annotate_graph(graph: AssemblyGraph, alignments: Sequence[Alignment]) -> None:
    """Attach a compact per-segment summary of alignments for the renderer."""
    by_query: dict[str, list[Alignment]] = {}
    for aln in alignments:
        by_query.setdefault(aln.query, []).append(aln)
    for name, seg in graph.segments.items():
        hits = by_query.get(name, [])
        if not hits:
            seg.ref_hits = []
            continue
        hits.sort(key=lambda a: a.q_span, reverse=True)
        seg.ref_hits = [h.to_dict() for h in hits[:8]]


def best_placements(
    alignments: Sequence[Alignment],
    min_query_coverage: float = 0.0,
    min_identity: float = 0.0,
) -> dict[str, Alignment]:
    """One representative alignment per query: the longest primary hit.

    Merging collinear blocks first means a contig broken into several alignment
    blocks by small indels is still placed by its true full extent.
    """
    grouped: dict[str, list[Alignment]] = {}
    for aln in alignments:
        if aln.identity < min_identity:
            continue
        grouped.setdefault(aln.query, []).append(aln)

    best: dict[str, Alignment] = {}
    for query, hits in grouped.items():
        merged = merge_collinear(hits)
        if not merged:
            continue
        candidate = max(merged, key=lambda a: (a.is_primary, a.q_span, a.matches))
        if candidate.query_coverage < min_query_coverage:
            continue
        best[query] = candidate
    return best


def merge_collinear(
    hits: Sequence[Alignment],
    max_gap: int = 10_000,
    max_overlap: int | None = None,
) -> list[Alignment]:
    """Join alignment blocks of one query that continue on the same reference.

    Blocks qualify when they share a reference and strand and advance in both
    query and reference coordinates by no more than ``max_gap``, and step
    backwards by no more than ``max_overlap`` (which defaults to ``max_gap``).

    The window matters and callers should choose it deliberately. Placing a
    contig wants a *wide* window, so a contig split into blocks by a deletion
    still yields one placement spanning its true extent. Classifying
    misassemblies wants a *narrow* one, because every gap the merge absorbs is a
    discrepancy that then cannot be reported.
    """
    if not hits:
        return []
    if max_overlap is None:
        # Blocks routinely overlap a little where an alignment was split, but a
        # large step *backwards* on the reference is a rearrangement, not a
        # continuation -- so the backward tolerance stays tight even when the
        # forward window is wide.
        max_overlap = min(500, max_gap)
    groups: dict[tuple[str, str, int], list[Alignment]] = {}
    for h in hits:
        groups.setdefault((h.query, h.ref, h.strand), []).append(h)

    merged: list[Alignment] = []
    for (query, ref, strand), blocks in groups.items():
        blocks.sort(key=lambda a: a.q_st)
        current = _copy_alignment(blocks[0])
        for nxt in blocks[1:]:
            q_gap = nxt.q_st - current.q_en
            if strand > 0:
                r_gap = nxt.r_st - current.r_en
            else:
                r_gap = current.r_st - nxt.r_en
            if -max_overlap <= q_gap <= max_gap and -max_overlap <= r_gap <= max_gap:
                current.q_en = max(current.q_en, nxt.q_en)
                current.r_st = min(current.r_st, nxt.r_st)
                current.r_en = max(current.r_en, nxt.r_en)
                current.matches += nxt.matches
                current.block_len += nxt.block_len
                current.nm += nxt.nm
                current.mapq = max(current.mapq, nxt.mapq)
                current.is_primary = current.is_primary or nxt.is_primary
                current.cigar = current.cigar + nxt.cigar
            else:
                merged.append(current)
                current = _copy_alignment(nxt)
        merged.append(current)
    return merged


def split_at_long_indels(aln: Alignment, min_indel: int = 1000) -> list[Alignment]:
    """Break an alignment wherever its CIGAR contains a long insertion/deletion.

    minimap2 will happily carry a single alignment straight across a 14 kb
    deletion, recording it as one ``D`` operation. Left alone, that hides a real
    structural difference (the reference gap never shows up as a break) *and*
    inflates genome fraction, because the untouched region between ``r_st`` and
    ``r_en`` is counted as covered. Splitting restores both.
    """
    ops = aln.cigar
    if not ops or not any(c >= min_indel and op in "ID" for c, op in ops):
        return [aln]

    pieces: list[Alignment] = []
    # The split-out indels are removed from the alignment entirely, so they must
    # not be part of the denominator when apportioning matches and edit distance.
    total_block = max(
        1,
        sum(
            c
            for c, op in ops
            if op in "M=XID" and not (op in "ID" and c >= min_indel)
        ),
    )
    # NM counts the split-out indels as edit operations. Since those indels are
    # no longer part of any piece, drop them before apportioning, or a single
    # 14 kb deletion turns into thousands of phantom mismatches.
    removed = sum(c for c, op in ops if op in "ID" and c >= min_indel)
    distributable_nm = max(0, aln.nm - removed)

    q_off = 0  # query bases consumed, in CIGAR orientation
    r_pos = aln.r_st
    cur_q_start = 0
    cur_r_start = aln.r_st
    cur_block = 0
    cur_ops: list[tuple[int, str]] = []

    def emit(q_end_cigar: int, r_end: int, block: int, sub_ops: list[tuple[int, str]]) -> None:
        if q_end_cigar <= cur_q_start or r_end <= cur_r_start or block <= 0:
            return
        if aln.strand > 0:
            q_st = aln.q_st + cur_q_start
            q_en = aln.q_st + q_end_cigar
        else:
            q_st = aln.q_en - q_end_cigar
            q_en = aln.q_en - cur_q_start
        share = block / total_block
        pieces.append(
            Alignment(
                query=aln.query,
                q_len=aln.q_len,
                q_st=max(0, q_st),
                q_en=min(aln.q_len, q_en),
                strand=aln.strand,
                ref=aln.ref,
                r_len=aln.r_len,
                r_st=cur_r_start,
                r_en=r_end,
                matches=int(round(aln.matches * share)),
                block_len=block,
                mapq=aln.mapq,
                nm=int(round(distributable_nm * share)),
                is_primary=aln.is_primary,
                cigar=list(sub_ops),
            )
        )

    for count, op in ops:
        if op in "ID" and count >= min_indel:
            emit(q_off, r_pos, cur_block, cur_ops)
            if op == "I":
                q_off += count
            else:
                r_pos += count
            cur_q_start = q_off
            cur_r_start = r_pos
            cur_block = 0
            cur_ops = []
            continue
        if op in "M=X":
            q_off += count
            r_pos += count
            cur_block += count
            cur_ops.append((count, op))
        elif op == "I":
            q_off += count
            cur_block += count
            cur_ops.append((count, op))
        elif op in "SH":
            q_off += count
        elif op in "DN":
            r_pos += count
            cur_block += count
            cur_ops.append((count, op))
    emit(q_off, r_pos, cur_block, cur_ops)

    return pieces or [aln]


def split_all(alignments: Sequence[Alignment], min_indel: int = 1000) -> list[Alignment]:
    out: list[Alignment] = []
    for aln in alignments:
        out.extend(split_at_long_indels(aln, min_indel))
    return out


def _copy_alignment(a: Alignment) -> Alignment:
    return Alignment(
        query=a.query,
        q_len=a.q_len,
        q_st=a.q_st,
        q_en=a.q_en,
        strand=a.strand,
        ref=a.ref,
        r_len=a.r_len,
        r_st=a.r_st,
        r_en=a.r_en,
        matches=a.matches,
        block_len=a.block_len,
        mapq=a.mapq,
        nm=a.nm,
        is_primary=a.is_primary,
        cigar=list(a.cigar),
    )
