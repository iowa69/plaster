# Assemblage — Interactive Assembly Graph Studio

**Date:** 2026-08-03
**Status:** Approved for implementation (autonomous session; assumptions recorded below)

## Purpose

A single tool that does what Bandage does (interactive assembly-graph visualisation) and what
QUAST does (reference-based assembly quality assessment), plus the things people always wish
Bandage had: a real *rearrange* control, the ability to drop a reference into the view, and
reference-guided scaffolding whose output you can actually export and use downstream.

The immediate consumer is a newly written short-read de novo assembler. Its output is assumed to
be GFA (the de facto standard); FASTG, plain FASTA and SPAdes path files are also supported so the
tool is useful against SPAdes/SKESA/Unicycler output for comparison.

## Assumptions (recorded because this was built autonomously)

1. Input graphs are bacterial-to-small-eukaryote scale (up to ~10^5 segments). Layout runs in the
   browser; larger graphs are handled by component filtering rather than by a GPU renderer.
2. The reference is a FASTA file of one or more sequences (chromosomes/plasmids/scaffolds).
3. `minimap2` (via the `mappy` Python binding) is an acceptable dependency for alignment. BLAST is
   optional and used only for the sequence-search panel when present.
4. "Easy to use" means: `conda env create -f environment.yml`, then one command
   (`assemblage view graph.gfa`) that opens a browser. No build step, no Qt, no X11 quirks.

## Architecture

Three layers, each independently testable.

**Core library (`assemblage.core`)** — pure Python, no I/O beyond files, no web knowledge.
Parsers/writers, the graph model, analysis, and scaffolding all live here and are usable from a
script or notebook. Every feature in the GUI is reachable from this layer, which is what makes the
CLI subcommands possible.

**Service layer (`assemblage.server`)** — a FastAPI app that holds one `Project` in memory and
exposes it over REST/JSON. It owns session state, undo/redo, and file export. It contains no
biology.

**Presentation (`assemblage.web`)** — static HTML/CSS/JS, no build step. A canvas renderer, a
force-directed layout worker, and the side panels. Talks to the service layer over `fetch`.

```
 GFA/FASTG/FASTA ─▶ io.* ─▶ AssemblyGraph ─┬─▶ analysis.metrics ──▶ report.html
                                            ├─▶ analysis.align (mappy) ──▶ misassemblies
                                            ├─▶ scaffold.reference_guided ──▶ FASTA + AGP
                                            └─▶ layout (browser) ──▶ canvas / SVG / PNG
```

## Components

### `core/model.py` — `AssemblyGraph`

Segments and links, using Bandage's double-stranded convention: a segment is drawn as one polyline
with a *start* end (5' of the `+` strand) and an *end* end. A link `(A,+) → (B,-)` therefore joins
`A.end` to `B.end`. Reverse complement is a view, not a copy.

Responsibilities: adjacency, connected components, degree, dead ends, path traversal, and the
mutation operations the GUI needs (delete, merge simple path, reverse, duplicate-into-scaffold).
Anything that changes the graph returns an inverse operation so undo is free.

### `core/io/` — formats

`gfa.py` (GFA 1.0 read/write incl. `P` paths and `dp/DP/KC/RC/FC` depth tags), `fastg.py` (SPAdes),
`fasta.py`, `agp.py` (AGP 2.1 writer), `paths.py` (SPAdes `contigs.paths`). Each parser is a
generator over records so a 5 GB GFA doesn't have to fit in memory twice.

### `core/analysis/`

- `metrics.py` — reference-free: count, total length, N50/L50/N75, NG50/LG50 (needs genome size),
  largest contig, GC, depth distribution, dead ends, connected components, Nx and cumulative-length
  series for plots.
- `align.py` — thin wrapper over `mappy`. Aligns segments to a reference, returns normalised
  alignment blocks. Falls back to the `minimap2` binary if `mappy` is missing.
- `misassembly.py` — QUAST-style classification of adjacent alignment blocks within one contig:
  relocation (>1 kb ref gap, same ref/strand), local misassembly (200 bp–1 kb), inversion (strand
  flip), translocation (different reference sequence). Plus genome fraction, duplication ratio, and
  mismatches/indels per 100 kb derived from CIGAR and NM.
- `search.py` — BLAST (if available) or minimap2 query search over segments, for highlighting.

### `core/scaffold/`

- `reference_guided.py` — the headline feature. Align contigs → pick best placement per contig →
  order and orient by reference coordinate → compute gaps from reference distance. Contigs that
  can't be placed confidently become unplaced singletons rather than being silently dropped.
- `graph_gap_fill.py` — where two adjacent contigs in a scaffold are connected in the assembly
  graph by a path of roughly the right length, splice the real path sequence in instead of a run of
  Ns. This is the part that makes the output better than reference scaffolding alone, because the
  gap sequence comes from the reads, not from the reference.
- `builder.py` — turns an ordered list of (segment, orientation, gap) into scaffold FASTA plus a
  matching AGP 2.1 file, so the two can never disagree.

### `server/` — FastAPI

Endpoints for project load, graph fetch (with level-of-detail filtering), mutation, alignment,
scaffolding, metrics, and export. One in-memory `Project`; an undo stack of inverse operations;
autosave of session JSON next to the input file.

### `web/` — renderer

Canvas 2D with viewport culling. Layout is force-directed, run in a Web Worker so the UI stays
responsive, with **Rearrange** exposed as a first-class control: whole graph, selection only, or
one component, in force-directed / linear / circular / component-grid modes. Colouring by depth,
GC, length, component, reference chromosome, BLAST hit, or a user CSV.

## Data flow for the two headline workflows

**Insert reference:** user picks a FASTA in the sidebar → `POST /api/reference` → `align.py` maps
every segment → alignment blocks stored on the project → nodes gain `ref_hits` → the renderer
recolours by chromosome and the QC panel fills in genome fraction, misassemblies, NGA50.

**Scaffold and export:** user clicks *Build scaffolds* → `scaffold/reference_guided.py` produces an
ordered scaffold plan → the plan is rendered in an editable list (drag to reorder, click to flip,
edit gap size) → *Export* writes `scaffolds.fasta` + `scaffolds.agp` + a QC report for the
scaffolded assembly. The plan is data, so the manual edits and the automatic method share one code
path.

## Error handling

Parsers raise a typed `AssemblageFormatError` carrying file, line number, and the offending text;
the server turns these into 400s that the UI shows verbatim, because a truncated GFA is the single
most common failure and a stack trace helps nobody. Missing optional dependencies (`mappy`, BLAST)
degrade features with a clear message rather than crashing at import. Alignment of a huge reference
is run in a thread with progress reported over a polling endpoint.

## Testing

`pytest` over the core library: round-trip GFA read/write, known-answer N50/NG50, a hand-built
graph whose components and dead ends are known, a synthetic reference/contig pair with a
deliberately introduced inversion and translocation so misassembly classification is checked
against ground truth, and a scaffolding case where FASTA and AGP must agree base-for-base. The web
layer gets a smoke test that boots the server and checks the API contract.

## Out of scope

Long-read/HiC scaffolding, gene annotation, variant calling, multi-user sessions, and a
publication-quality figure editor. Colour/label control is good enough to make a figure, but this
is not Illustrator.
