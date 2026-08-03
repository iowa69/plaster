# Validation

Assemblage is checked against the two tools it sets out to replace, on real
bacterial assemblies rather than on data Assemblage generated itself. Testing a
tool only against its own synthetic fixtures is circular: it confirms the code
does what it does, not that it is right.

Everything below is reproducible with the commands given.

---

## Graph statistics vs. Bandage

`Bandage info` and `assemblage info` were run on the same GFA files.

**A 2.88 Mb assembly from the author's own assembler** (`assembly_graph.gfa`,
94 bp overlaps):

| | Bandage | Assemblage |
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

| | Bandage | Assemblage |
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
assemblage info assembly_graph.gfa
```

---

## Reference-based statistics vs. QUAST

QUAST 5.x and Assemblage were run on the same contigs and the same closed
reference (a *K. pneumoniae* chromosome plus two plasmids, 5,509,378 bp).
Assemblage was given `--min-contig 500` to match QUAST's default.

| | QUAST | Assemblage |
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
assemblage qc segments.fasta -r reference.fasta --preset asm5 --min-contig 500
```

### Two bugs this comparison found

Both were real, and neither showed up against synthetic data.

**Contigs spanning a circular origin were reported as relocations.** Bacterial
chromosomes and plasmids are circular. A contig that spans the origin ends at
the last base of the reference and resumes at base 0, which on a linear reading
looks like a jump of nearly the whole replicon. Assemblage reported three such
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

---

## Scaffolding, judged by QUAST

The headline feature evaluated by an independent tool. Assemblage scaffolded the
*K. pneumoniae* graph against its reference, and QUAST then scored the contigs
and the resulting scaffolds side by side.

```bash
assemblage scaffold assembly_graph.gfa -r reference.fasta --preset asm5 -o out/
quast.py -o quast_out -r reference.fasta segments.fasta out/scaffolds.fasta
```

| | contigs | scaffolds |
|---|---|---|
| # contigs | 63 | **5** |
| Largest contig | 812,461 | **5,259,000** |
| N50 | 444,072 | **5,259,000** |
| NGA50 | 444,072 | **5,241,626** |
| Genome fraction (%) | 98.571 | **99.570** |
| Duplication ratio | 1.001 | 1.001 |
| **# misassemblies** | 0 | **0** |
| # local misassemblies | 0 | 0 |
| mismatches / 100 kb | 0.28 | 2.28 |
| indels / 100 kb | 0.00 | 1.55 |

N50 improved 11.8-fold and the chromosome was reconstructed as a single
5.26 Mb scaffold, **with no misassemblies introduced** and no increase in
duplication. Genome fraction *rose* by one percentage point, because 48 gaps
were closed with 58 kb of real sequence recovered from the graph rather than
filled with Ns.

The small rise in mismatch and indel rate is expected and is not degradation:
the newly incorporated bridge sequence had not previously been aligned to
anything, and it carries the assembler's own error rate. 2.28 mismatches per
100 kb is 99.998% identity.

### The AGP contract on real data

Every AGP row was checked against the scaffold FASTA it describes: 125 component
rows sliced out exactly the oriented contig sequence, 10 gap rows were exactly
the declared number of Ns, 31 graph-bridge rows carried their recovered
sequence, coordinates tiled each scaffold with no gap or overlap, and every
scaffold's length equalled its AGP extent. This held with 127 bp overlaps in
play, where the bridging sequence must be overlap-trimmed to be correct.

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

## Test suite

389 tests, running in under 3 seconds:

```bash
conda activate assemblage
python -m pytest tests/ -v
```

They cover the double-stranded link-end convention, format round-trips,
known-answer N50/NG50/auN, alignment splitting on both strands, misassembly
classification against ground truth, the FASTA/AGP contract, graph editing and
undo, and the HTTP API contract.
