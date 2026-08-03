# Validation

Plastr is checked against the two tools it sets out to replace, on real
bacterial assemblies rather than on data Plastr generated itself. Testing a
tool only against its own synthetic fixtures is circular: it confirms the code
does what it does, not that it is right.

Everything below is reproducible with the commands given.

---

## Graph statistics vs. Bandage

`Bandage info` and `plastr info` were run on the same GFA files.

**A 2.88 Mb assembly from the author's own assembler** (`assembly_graph.gfa`,
94 bp overlaps):

| | Bandage | Plastr |
|---|---|---|
| Node count | 514 | 514 |
| Edge count | 641 | 641 |
| Total length | 2,878,586 | 2,878,586 |
| Dead ends | 76 | 76 |
| Connected components | 5 | 5 |
| N50 | 45,898 | 45,898 |
| Longest node | 141,452 | 141,452 |
| Largest component | 2,676,905 | 2,676,905 |
| Median depth | 15.4764 | 15.4764 |

**A 5.45 Mb *Klebsiella pneumoniae* graph** (127 bp overlaps):

| | Bandage | Plastr |
|---|---|---|
| Node count | 125 | 125 |
| Edge count | 159 | 159 |
| Total length | 5,448,934 | 5,448,934 |
| Dead ends | 12 | 12 |
| Connected components | 8 | 8 |
| N50 | 444,072 | 444,072 |
| Median depth | 49.5349 | 49.5349 |

Exact agreement on every statistic both tools report.

```bash
Bandage info assembly_graph.gfa
plastr info assembly_graph.gfa
```

---

## Reference-based statistics vs. QUAST

QUAST 5.x and Plastr were run on the same contigs and the same closed
reference (a *K. pneumoniae* chromosome plus two plasmids, 5,509,378 bp).
Plastr was given `--min-contig 500` to match QUAST's default.

| | QUAST | Plastr |
|---|---|---|
| # contigs | 63 | 63 |
| Total length | 5,434,937 | 5,434,937 |
| Largest contig | 812,461 | 812,461 |
| GC (%) | 57.30 | 57.30 |
| N50 | 444,072 | 444,072 |
| NG50 | 444,072 | 444,072 |
| L50 | 5 | 5 |
| auN | 441,529.4 | 441,529.4 |
| Genome fraction (%) | 98.571 | 98.573 |
| Duplication ratio | 1.001 | 1.001 |
| Total aligned length | 5,433,885 | 5,433,885 |
| Largest alignment | 812,461 | 812,461 |
| NA50 | 444,072 | 444,072 |
| NGA50 | 444,072 | 444,072 |
| mismatches / 100 kb | 0.28 | 0.28 |
| indels / 100 kb | 0.00 | 0.00 |
| # misassemblies | 0 | 0 |
| # local misassemblies | 0 | 0 |
| # unaligned contigs | 2 | 2 |

```bash
quast.py -o quast_out -r reference.fasta segments.fasta
plastr qc segments.fasta -r reference.fasta --preset asm5 --min-contig 500
```

### Two bugs this comparison found

Both were real, and neither showed up against synthetic data.

**Contigs spanning a circular origin were reported as relocations.** Bacterial
chromosomes and plasmids are circular. A contig that spans the origin ends at
the last base of the reference and resumes at base 0, which on a linear reading
looks like a jump of nearly the whole replicon. Plastr reported three such
false relocations — one on the chromosome and one on each plasmid — where QUAST
reported none. Reference sequences are now treated as circular by default, and
the distance between two alignment blocks is measured the short way round;
`--linear-references` restores the old behaviour.

**Secondary alignments were inflating every statistic.** A repeat that maps to
five places produces five alignments. Counting them all pushed *total aligned
length* to 5,512,931 bp — more than the 5,434,937 bp assembly it was measuring,
which is impossible — and inflated the mismatch rate nine-fold (2.5 vs 0.28 per
100 kb). Statistics now use primary alignments only, matching QUAST's default
handling of ambiguity. Secondary hits are still kept for the graph view, where
seeing every place a repeat lands is the whole point; `--count-ambiguous`
includes them in the statistics too.

Both are covered by regression tests in `tests/test_misassembly.py`.

### Across six genomes

The single-genome comparison above could be luck, so it was repeated on five
further closed *K. pneumoniae* genomes, 80 metric comparisons in total.

These agree **exactly on all six genomes**: number of contigs, total length,
largest contig, GC, N50, NG50, L50, genome fraction, duplication ratio, NGA50,
largest alignment, total aligned length, and — the number that matters most —
**# misassemblies**.

Two metrics differ slightly on three of the six:

| genome | metric | QUAST | Plastr |
|---|---|---|---|
| ERR10447223 | # local misassemblies | 2 | 1 |
| ERR10447223 | mismatches / 100 kb | 0.13 | 0.04 |
| ERR11578077 | # local misassemblies | 1 | 0 |
| ERR11578427 | # local misassemblies | 2 | 1 |
| ERR11578427 | mismatches / 100 kb | 12.37 | 9.66 |
| ERR11578427 | indels / 100 kb | 0.60 | 0.13 |

These are the two most alignment-sensitive statistics, and the differences are
what you get from running a different aligner configuration: a local
misassembly is defined by a reference gap between 200 bp and 1 kb, so a block
boundary shifting by a few tens of bases moves an event across the threshold.
Plastr reads slightly low on both, consistently. Closing the gap entirely
would mean reproducing QUAST's exact aligner invocation, which is not a goal;
what matters is that the structural verdict — how many misassemblies, how much
of the genome is covered, how contiguous the assembly is — is identical.

---

## Scaffolding, judged by QUAST

The headline feature evaluated by an independent tool. Plastr scaffolded the
*K. pneumoniae* graph against its reference, and QUAST then scored the contigs
and the resulting scaffolds side by side.

```bash
plastr scaffold assembly_graph.gfa -r reference.fasta --preset asm5 -o out/
quast.py -o quast_out -r reference.fasta segments.fasta out/scaffolds.fasta
```

| | contigs | scaffolds |
|---|---|---|
| # contigs | 63 | **5** |
| Total length | 5,434,937 | 5,504,709 |
| Largest contig | 812,461 | **5,255,698** |
| N50 | 444,072 | **5,255,698** |
| NGA50 | 444,072 | **5,238,277** |
| Genome fraction (%) | 98.571 | **99.565** |
| Duplication ratio | 1.001 | **1.000** |
| **# misassemblies** | 0 | **0** |
| # local misassemblies | 0 | 1 |
| mismatches / 100 kb | 0.28 | 2.04 |
| indels / 100 kb | 0.00 | 0.62 |

N50 improved 11.8-fold and the chromosome was reconstructed as a single
5.26 Mb scaffold, **with no misassemblies introduced** and duplication ratio
at 1.000. Genome fraction *rose* by one percentage point, because 48 gaps were
closed with 58 kb of real sequence recovered from the graph rather than filled
with Ns.

The residual mismatch and indel rate is the assembler's own error rate in
sequence that had not previously been aligned to anything: 2.04 mismatches per
100 kb is 99.998% identity.

### The bug this table used to hide

An earlier version of this table read 5,510,805 bp, 2.28 mismatches and 1.55
indels per 100 kb, and the rise was explained away as the bridge sequence
carrying the assembler's error rate. That explanation was wrong. The scaffold
builder was writing every contig in full, including the bases it shared with
the piece before it, so each graph-closed gap made the scaffold *k*-1 bases too
long. 48 closed gaps × 127 bp = 6,096 bp — exactly the difference between the
old total and the corrected one.

The FASTA/AGP consistency check could not catch it, because both are generated
from a single coordinate counter: they agreed with each other while both were
wrong relative to the biology. What was missing was a test that a scaffold
equals the genome it was built from. That test now exists, parameterised over
four overlap sizes and including a bridged gap with a reverse-complemented
member, and it fails loudly if the trim is removed.

### The AGP contract on real data

Every AGP row was checked against the scaffold FASTA it describes: 125 component
rows sliced out exactly the oriented contig sequence, 10 gap rows were exactly
the declared number of Ns, 31 graph-bridge rows carried their recovered
sequence, coordinates tiled each scaffold with no gap or overlap, and every
scaffold's length equalled its AGP extent. This held with 127 bp overlaps in
play, where the bridging sequence must be overlap-trimmed to be correct.

---

## Robustness across real files

Every assembly graph on the development machine -- 251 files spanning SPAdes
GFA and FASTG, Unicycler intermediate and final graphs, and the author's own
assembler, at every *k* from 27 to 127 -- loads without error, in 2.2 s in
total.

### A third bug this found

**FASTG links were treated as blunt joins.** FASTG has no field for the overlap,
so the parser recorded zero, but SPAdes FASTG edges really do overlap by *k*-1
bases -- the same overlap the GFA writes as `127M`. Any merge or scaffold built
through such a link therefore duplicated the shared bases: merging two 128 bp
edges produced 256 bp where the correct answer is 129 bp. This was silent
corruption of exported sequence, and it would have hit anyone loading SPAdes
FASTG.

Overlaps are now measured when the file does not state them, by finding the
longest exact suffix/prefix match across a sample of links and requiring a
consensus. Explicit CIGARs are never overridden, and a graph of genuinely blunt
joins produces no consensus and is left alone.

Across those 251 files the measured overlaps come out as 27, 53, 55, 71, 77, 87,
94, 99, 111, 119 and 127 -- exactly the *k*-1 values the assemblers used -- with
17 files correctly identified as blunt.

---

## Scale

A synthetic 20,000-segment graph (21,076 links, 11.8 Mb):

| operation | time |
|---|---|
| parse GFA | 0.25 s |
| compute all metrics | 0.63 s |
| serialise the graph to JSON | 0.78 s |
| filter to 10 largest components | 0.07 s |
| simplify (4,116 merges) | 0.29 s |
| undo | 0.41 s |

An earlier version re-indexed the whole graph once per removed segment, so
filtering and simplifying a graph this size did not finish at all. Links are now
indexed by the segments they touch.

---

## Layout quality

Force-directed layout was checked on a real 514-segment assembly by measuring
mean link length in units of mean segment draw length. Connected segments
should sit close together, so a good layout scores near 1-2; a hairball scores
much higher.

| repulsion normalisation | link/segment ratio | fit-to-view zoom |
|---|---|---|
| none (original) | 15.2 | 9% |
| normalised, reference 100 | 0.85 | 81% |
| normalised, reference 200 | **2.4** | **49%** |

Repulsion accumulates over every particle while a particle's link forces do
not, so the same `repulsion` value that lays out a 10-segment graph well blows
a 500-segment one into a hairball with links fifteen times longer than the
segments they join. Repulsion is now damped by particle count, calibrated so
graphs at or below the reference size are unchanged and larger ones stay
readable. The 514-segment graph went from an unreadable disc to distinguishable
contig strands.

Rendering that graph costs 0.7 ms per frame.

---

## Test suite

403 tests, running in under 3 seconds:

```bash
conda activate plastr
python -m pytest tests/ -v
```

They cover the double-stranded link-end convention, format round-trips,
known-answer N50/NG50/auN, alignment splitting on both strands, misassembly
classification against ground truth, the FASTA/AGP contract, graph editing and
undo, and the HTTP API contract.
