"""FASTA reading and writing, with transparent gzip support."""

from __future__ import annotations

import gzip
import io
import os
from pathlib import Path
from typing import IO, Iterator

from ..errors import PlastrFormatError


def open_text(path: str | os.PathLike[str]) -> IO[str]:
    """Open a possibly-gzipped text file."""
    p = Path(path)
    with open(p, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(p, "rt")
    return open(p, "r")


def read_fasta(path: str | os.PathLike[str]) -> Iterator[tuple[str, str, str]]:
    """Yield ``(name, description, sequence)`` for each record.

    ``name`` is the header up to the first whitespace; ``description`` is the rest.
    """
    name: str | None = None
    description = ""
    chunks: list[str] = []
    seen_any = False
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                seen_any = True
                if name is not None:
                    yield name, description, "".join(chunks)
                header = line[1:].strip()
                if not header:
                    raise PlastrFormatError(
                        "FASTA record has an empty header", str(path), line_no, raw
                    )
                parts = header.split(None, 1)
                name = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                chunks = []
            elif line.startswith(";"):
                continue  # ancient FASTA comment
            else:
                if name is None:
                    raise PlastrFormatError(
                        "sequence data before any '>' header", str(path), line_no, raw
                    )
                chunks.append(line.strip())
    if name is not None:
        yield name, description, "".join(chunks)
    elif not seen_any:
        raise PlastrFormatError("file contains no FASTA records", str(path))


def read_fasta_dict(path: str | os.PathLike[str]) -> dict[str, str]:
    """Load a FASTA into ``{name: sequence}``. Duplicate names raise."""
    out: dict[str, str] = {}
    for name, _desc, seq in read_fasta(path):
        if name in out:
            raise PlastrFormatError(f"duplicate sequence name {name!r}", str(path))
        out[name] = seq
    return out


def write_fasta(
    path: str | os.PathLike[str] | IO[str],
    records: "Iterator[tuple[str, str]] | list[tuple[str, str]]",
    wrap: int = 60,
) -> int:
    """Write ``(name, sequence)`` records. Returns the number written.

    ``wrap=0`` writes each sequence on a single line.
    """
    handle: IO[str]
    should_close = False
    if isinstance(path, (str, os.PathLike)):
        handle = open(path, "w")
        should_close = True
    else:
        handle = path
    count = 0
    try:
        for name, seq in records:
            handle.write(f">{name}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    handle.write(seq[i : i + wrap])
                    handle.write("\n")
                if not seq:
                    handle.write("\n")
            else:
                handle.write(seq)
                handle.write("\n")
            count += 1
    finally:
        if should_close:
            handle.close()
    return count


def fasta_string(records: "list[tuple[str, str]]", wrap: int = 60) -> str:
    buf = io.StringIO()
    write_fasta(buf, records, wrap)
    return buf.getvalue()
