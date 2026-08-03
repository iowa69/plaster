"""Searching for a query sequence among the graph's segments.

Uses BLAST when it is installed (better for short and divergent queries, which
is what people usually paste in -- a gene, a primer pair, a resistance cassette)
and falls back to minimap2 otherwise.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass

from ..errors import AssemblageFormatError, MissingDependencyError
from ..io.fasta import read_fasta, write_fasta
from ..model import AssemblyGraph


@dataclass(slots=True)
class SearchHit:
    segment: str
    query: str
    identity: float
    length: int
    q_st: int
    q_en: int
    s_st: int
    s_en: int
    strand: int
    bitscore: float = 0.0
    evalue: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def blast_available() -> bool:
    return shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None


def search_backend() -> str:
    if blast_available():
        return "blast"
    from .align import alignment_backend

    backend = alignment_backend()
    return "minimap2" if backend != "none" else "none"


def parse_query(query: str) -> list[tuple[str, str]]:
    """Accept a raw sequence, pasted FASTA text, or a path to a FASTA file."""
    text = (query or "").strip()
    if not text:
        raise AssemblageFormatError("empty search query")

    if os.path.exists(text) and len(text) < 4096 and "\n" not in text:
        return [(name, seq) for name, _d, seq in read_fasta(text)]

    if text.startswith(">"):
        records: list[tuple[str, str]] = []
        name: str | None = None
        chunks: list[str] = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0] or f"query{len(records) + 1}"
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            records.append((name, "".join(chunks)))
        if not records:
            raise AssemblageFormatError("no sequence found in the pasted FASTA")
        return records

    cleaned = "".join(text.split()).upper()
    invalid = set(cleaned) - set("ACGTURYSWKMBDHVN")
    if invalid:
        raise AssemblageFormatError(
            "query is neither a readable file, FASTA text, nor a nucleotide sequence "
            f"(unexpected characters: {''.join(sorted(invalid))[:20]})"
        )
    return [("query", cleaned)]


def search_graph(
    graph: AssemblyGraph,
    query: str,
    min_identity: float = 0.8,
    min_length: int = 0,
    max_hits: int = 500,
    threads: int = 4,
) -> tuple[str, list[SearchHit]]:
    """Find where ``query`` occurs in the assembly. Returns (backend, hits)."""
    queries = parse_query(query)
    subjects = [
        (name, seg.sequence)
        for name, seg in graph.segments.items()
        if seg.sequence
    ]
    if not subjects:
        raise AssemblageFormatError(
            "this graph was loaded without sequences, so it cannot be searched"
        )

    backend = search_backend()
    if backend == "blast":
        return "blast", _blast(queries, subjects, min_identity, min_length, max_hits, threads)
    if backend == "minimap2":
        return "minimap2", _minimap(queries, subjects, min_identity, min_length, max_hits, threads)
    raise MissingDependencyError(
        "blastn or minimap2",
        "Install with 'conda install -c bioconda blast minimap2'.",
    )


def _blast(
    queries: list[tuple[str, str]],
    subjects: list[tuple[str, str]],
    min_identity: float,
    min_length: int,
    max_hits: int,
    threads: int,
) -> list[SearchHit]:
    with tempfile.TemporaryDirectory(prefix="assemblage_blast_") as tmp:
        db = os.path.join(tmp, "segments.fa")
        qf = os.path.join(tmp, "query.fa")
        write_fasta(db, subjects)
        write_fasta(qf, queries)
        make = subprocess.run(
            ["makeblastdb", "-in", db, "-dbtype", "nucl", "-out", os.path.join(tmp, "db")],
            capture_output=True,
            text=True,
        )
        if make.returncode != 0:
            raise RuntimeError(f"makeblastdb failed: {make.stderr.strip()[:400]}")
        proc = subprocess.run(
            [
                "blastn",
                "-query", qf,
                "-db", os.path.join(tmp, "db"),
                "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send bitscore evalue",
                "-perc_identity", str(max(0.0, min_identity * 100)),
                "-max_target_seqs", str(max_hits),
                "-num_threads", str(threads),
                "-dust", "no",
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(f"blastn failed: {proc.stderr.strip()[:400]}")

        hits: list[SearchHit] = []
        for line in proc.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 10:
                continue
            length = int(f[3])
            if length < min_length:
                continue
            s_st, s_en = int(f[6]), int(f[7])
            strand = 1 if s_en >= s_st else -1
            hits.append(
                SearchHit(
                    segment=f[1],
                    query=f[0],
                    identity=float(f[2]) / 100.0,
                    length=length,
                    q_st=int(f[4]) - 1,
                    q_en=int(f[5]),
                    s_st=min(s_st, s_en) - 1,
                    s_en=max(s_st, s_en),
                    strand=strand,
                    bitscore=float(f[8]),
                    evalue=float(f[9]),
                )
            )
        hits.sort(key=lambda h: h.bitscore, reverse=True)
        return hits[:max_hits]


def _minimap(
    queries: list[tuple[str, str]],
    subjects: list[tuple[str, str]],
    min_identity: float,
    min_length: int,
    max_hits: int,
    threads: int,
) -> list[SearchHit]:
    from .align import align_sequences

    with tempfile.TemporaryDirectory(prefix="assemblage_search_") as tmp:
        subject_path = os.path.join(tmp, "segments.fa")
        write_fasta(subject_path, subjects)
        # Query and subject are swapped relative to reference alignment: here the
        # segments are the "reference" we are searching within.
        alignments = align_sequences(
            queries, subject_path, preset="asm20", min_length=0, threads=threads,
            split_long_indels=0,
        )
    hits = [
        SearchHit(
            segment=a.ref,
            query=a.query,
            identity=a.identity,
            length=a.block_len,
            q_st=a.q_st,
            q_en=a.q_en,
            s_st=a.r_st,
            s_en=a.r_en,
            strand=a.strand,
            bitscore=float(a.matches),
        )
        for a in alignments
        if a.identity >= min_identity and a.block_len >= min_length
    ]
    hits.sort(key=lambda h: h.bitscore, reverse=True)
    return hits[:max_hits]
