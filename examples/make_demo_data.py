#!/usr/bin/env python3
"""Build a small demo dataset with known ground truth.

Produces a reference genome and an assembly graph derived from it, into which
specific errors have been introduced on purpose:

  * one **inversion**    -- a contig whose second half is reverse complemented
  * one **translocation** -- a chimeric contig joining chromosome and plasmid
  * one **relocation**   -- a contig that skips a chunk of the chromosome
  * SNPs at roughly 0.1% so mismatch rates are non-zero
  * a repeat present three times, which becomes a shared graph node

Run:  python examples/make_demo_data.py [outdir]
"""

from __future__ import annotations

import json
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from plastr.core.io.fasta import write_fasta  # noqa: E402
from plastr.core.sequence import revcomp  # noqa: E402

SEED = 20260803
CHROM_LEN = 160_000
PLASMID_LEN = 9_000
REPEAT_LEN = 2_400


def random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng: random.Random, seq: str, rate: float = 0.001) -> str:
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def main(outdir: str = "examples/demo") -> None:
    rng = random.Random(SEED)
    os.makedirs(outdir, exist_ok=True)

    # ---- reference ----
    repeat = random_seq(rng, REPEAT_LEN)
    chrom = random_seq(rng, CHROM_LEN)
    # Plant the repeat three times in the chromosome so the graph has a real
    # branching node rather than a trivial path.
    for pos in (20_000, 70_000, 120_000):
        chrom = chrom[:pos] + repeat + chrom[pos + REPEAT_LEN :]
    plasmid = random_seq(rng, PLASMID_LEN)

    reference = os.path.join(outdir, "reference.fasta")
    write_fasta(reference, [("chromosome", chrom), ("plasmid", plasmid)])

    # ---- contigs, with deliberate errors ----
    truth: dict[str, object] = {"expected": {}}
    contigs: list[tuple[str, str]] = []

    # Clean contigs covering most of the chromosome.
    clean_spans = [(0, 19_000), (23_000, 55_000), (58_000, 69_500), (75_000, 100_000)]
    for i, (a, b) in enumerate(clean_spans, 1):
        contigs.append((f"ctg_clean_{i}", mutate(rng, chrom[a:b])))

    # Inversion: 100k-108k followed by the reverse complement of 108k-116k.
    inv = chrom[100_000:108_000] + revcomp(chrom[108_000:116_000])
    contigs.append(("ctg_inversion", mutate(rng, inv)))

    # Relocation: 125k-131k joined directly to 145k-152k (skips 14 kb).
    reloc = chrom[125_000:131_000] + chrom[145_000:152_000]
    contigs.append(("ctg_relocation", mutate(rng, reloc)))

    # Translocation: chromosome 152k-158k fused to the first 5 kb of the plasmid.
    trans = chrom[152_000:158_000] + plasmid[:5_000]
    contigs.append(("ctg_translocation", mutate(rng, trans)))

    # Remainder of the plasmid, correct.
    contigs.append(("ctg_plasmid", mutate(rng, plasmid[5_000:])))

    # The collapsed repeat, as its own node (this is what makes it a graph).
    contigs.append(("ctg_repeat", mutate(rng, repeat)))

    # A short contig with no reference match at all (contamination).
    contigs.append(("ctg_foreign", random_seq(rng, 3_500)))

    truth["expected"] = {
        "inversions": 1,
        "relocations": 1,
        "translocations": 1,
        "unaligned_contigs": 1,
        "reference_length": len(chrom) + len(plasmid),
    }

    # ---- write GFA with links that reflect the repeat structure ----
    lines = ["H\tVN:Z:1.0"]
    depth_of = {"ctg_repeat": 33.0}
    for name, seq in contigs:
        depth = depth_of.get(name, round(rng.uniform(28.0, 38.0), 2))
        lines.append(f"S\t{name}\t{seq}\tdp:f:{depth}\tLN:i:{len(seq)}")

    links = [
        ("ctg_clean_1", "+", "ctg_repeat", "+"),
        ("ctg_repeat", "+", "ctg_clean_2", "+"),
        ("ctg_clean_3", "+", "ctg_repeat", "+"),
        ("ctg_repeat", "+", "ctg_clean_4", "+"),
        ("ctg_clean_4", "+", "ctg_inversion", "+"),
        ("ctg_inversion", "+", "ctg_relocation", "+"),
        ("ctg_relocation", "+", "ctg_translocation", "+"),
        ("ctg_translocation", "+", "ctg_plasmid", "+"),
        ("ctg_plasmid", "+", "ctg_translocation", "+"),
    ]
    for a, ao, b, bo in links:
        lines.append(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M")

    gfa = os.path.join(outdir, "assembly.gfa")
    with open(gfa, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    write_fasta(os.path.join(outdir, "contigs.fasta"), contigs)
    with open(os.path.join(outdir, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=2)

    print(f"reference : {reference}  ({len(chrom) + len(plasmid):,} bp, 2 sequences)")
    print(f"graph     : {gfa}  ({len(contigs)} segments, {len(links)} links)")
    print(f"contigs   : {os.path.join(outdir, 'contigs.fasta')}")
    print(f"truth     : {truth['expected']}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "examples/demo")
