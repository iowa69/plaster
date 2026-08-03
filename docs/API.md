# Plastr HTTP API contract

All endpoints are under `/api`. Request and response bodies are JSON unless
stated. Errors return `{"error": "...", "detail": "..."}` with a 4xx/5xx status;
the UI shows `error` verbatim.

The server holds exactly one project in memory.

---

## `GET /api/status`

```json
{
  "loaded": true,
  "graph": {"segments": 10, "links": 9, "paths": 0, "total_length": 137400,
            "components": 2, "dead_ends": 5, "has_sequences": true},
  "source_path": "examples/demo/assembly.gfa",
  "source_format": "gfa",
  "reference": {"path": "...", "sequences": 2, "length": 169000},
  "has_alignments": true,
  "has_plan": false,
  "align_backend": "mappy",
  "blast_available": true,
  "undo_depth": 0,
  "settings": {"min_contig": 0, "circular_references": true,
               "primary_only": true, "genome_size": null}
}
```
`reference` is `null` when none is loaded.

---

## `POST /api/settings`
Change the options that govern metrics and reference evaluation. Any subset may
be given; unlisted keys are left alone. Returns the settings now in effect.

```json
{"min_contig": 500, "circular_references": true, "primary_only": true,
 "genome_size": null}
```

| key | default | meaning |
|---|---|---|
| `min_contig` | `0` | ignore segments shorter than this in the statistics (500 matches QUAST) |
| `circular_references` | `true` | treat reference sequences as circular, so a contig spanning a replicon's origin is not a relocation |
| `primary_only` | `true` | exclude secondary alignments from the statistics, as QUAST does |
| `genome_size` | `null` | expected genome size for NG50/NGA50; falls back to the reference length |

These are also reported under `settings` in `GET /api/status`, and may be passed
in the body of `POST /api/reference` to set them at the same time as aligning.

---

## `POST /api/load`
Body: `{"path": "/abs/or/relative/path", "format": null}`
`format` is one of `gfa`, `gfa2`, `fastg`, `fasta`, or null to sniff.
Returns the same shape as `/api/status`.

## `POST /api/load-paths`
Body: `{"path": "contigs.paths"}` — attach SPAdes paths. Returns `{"added": 12}`.

---

## `GET /api/graph`
Query params: `min_length` (int, default 0), `max_nodes` (int, default 15000),
`component` (int or omitted → all).

```json
{
  "segments": [
    {"name": "ctg_1", "length": 19000, "depth": 31.2, "gc": 0.503,
     "component": 0, "deg_start": 0, "deg_end": 1, "circular": false,
     "ref_hits": [{"ref": "chromosome", "r_st": 0, "r_en": 19000,
                   "strand": 1, "identity": 0.9989, "q_st": 0, "q_en": 19000,
                   "mapq": 60, "is_primary": true, "q_len": 19000,
                   "r_len": 160000, "block_len": 19000}]}
  ],
  "links": [{"from": "a", "from_orient": "+", "to": "b", "to_orient": "+", "overlap": 0}],
  "paths": [{"name": "contig_1", "steps": ["1+", "2+", "4+"]}],
  "truncated": false,
  "shown": 10,
  "total": 10,
  "references": ["chromosome", "plasmid"]
}
```

**Link end convention (important for the renderer).** Each segment is drawn as a
polyline with a `start` end and an `end` end. A link joins:

| from_orient | to_orient | joins                |
|-------------|-----------|----------------------|
| `+`         | `+`       | A.end   → B.start    |
| `+`         | `-`       | A.end   → B.end      |
| `-`         | `+`       | A.start → B.start    |
| `-`         | `-`       | A.start → B.end      |

---

## `GET /api/segment/{name}`
Returns full detail including `sequence` (string) and `neighbours`.

---

## `POST /api/reference`
Body:
```json
{"path": "reference.fasta", "preset": "asm10", "min_identity": 0.0,
 "min_length": 200, "threads": 8}
```
`preset` ∈ `asm5` | `asm10` | `asm20`. Runs alignment (may take seconds), then
returns `{"status": "...", "report": <reference report>, "metrics": <metrics>}`.

## `DELETE /api/reference`
Drops the reference and its alignments.

---

## `GET /api/report`
```json
{
  "metrics": { "num_contigs": 10, "total_length": 137400, "n50": 19000,
               "l50": 3, "ng50": 13000, "gc_percent": 50.1, "auN": 20123.4,
               "largest_contig": 32000, "num_links": 9, "num_components": 2,
               "dead_ends": 5, "num_circular": 0, "mean_depth": 32.1,
               "nx_curve": [[1, 32000], ...], "cumulative_curve": [[1, 32000], ...],
               "length_histogram": [[1000, 3], ...] },
  "reference": {
      "genome_fraction": 82.06, "duplication_ratio": 1.0, "na50": 19000,
      "nga50": 11500, "mismatches_per_100kb": 92.3, "indels_per_100kb": 0.0,
      "num_misassemblies": 3, "num_relocations": 1, "num_inversions": 1,
      "num_translocations": 1, "num_local_misassemblies": 0,
      "unaligned_contigs": 1, "misassembled_contigs": ["ctg_x"],
      "per_reference_coverage": {"chromosome": 84.1},
      "coverage_blocks": {"chromosome": [[0, 19000], ...]},
      "misassemblies": [{"contig": "ctg_x", "kind": "inversion",
                         "is_extensive": true, "contig_pos": 8000,
                         "description": "..."}]
  }
}
```
`reference` is `null` when no reference is loaded.

---

## `POST /api/op`
Body: `{"op": "<name>", "args": {...}}`. Returns
`{"applied": "<name>", "count": 3, "graph": <status.graph>, "undo_depth": 1}`.

| op | args |
|----|------|
| `delete` | `{"names": ["a","b"]}` |
| `filter` | `{"min_length": 500, "min_depth": null, "max_depth": null, "keep_components": 5, "min_component_length": 0}` |
| `simplify` | `{}` — merge unbranching chains |
| `break_misassemblies` | `{"extensive_only": true}` — needs a reference |
| `split` | `{"name": "ctg_1", "positions": [5000]}` |
| `merge` | `{"steps": ["a+","b+"]}` |
| `reverse` | `{"name": "ctg_1"}` |

## `POST /api/undo`
Returns the same shape as `/api/op`.

---

## `POST /api/scaffold`
Body:
```json
{"method": "reference", "min_identity": 0.8, "min_query_coverage": 0.3,
 "min_align_length": 500, "min_gap": 100, "fill_gaps_from_graph": true,
 "include_unplaced": true, "break_misassemblies_first": false}
```
`method` ∈ `reference` | `graph`. Returns `{"plan": <plan>, "preview": <built stats>}`.

Plan shape:
```json
{"method": "reference", "scaffold_count": 3, "placed_count": 10,
 "unplaced": ["ctg_foreign"], "redundant": [], "notes": ["..."],
 "scaffolds": [
   {"name": "scaffold_chromosome", "source": "reference", "reference": "chromosome",
    "members": [{"segment": "ctg_1", "orientation": "+", "gap_after": 3000,
                 "gap_evidence": "reference", "bridge_path": ["r+"],
                 "bridge_sequence_length": 2400, "ref": "chromosome",
                 "ref_start": 0, "ref_end": 19000, "identity": 0.9989,
                 "overlaps_previous": false}]}
 ]}
```
`gap_evidence` ∈ `reference` | `graph` | `adjacent` | `manual` | `default`.

## `GET /api/scaffold` — current plan or `null`.
## `PUT /api/scaffold` — body is a plan (same shape); replaces it after manual edits.
   Returns `{"plan": ..., "preview": ...}`.

`preview`:
```json
{"scaffold_count": 3, "total_length": 182904, "gap_bases": 43104,
 "bridged_bases": 2400, "num_gaps": 4, "n50": 175404, "largest": 175404,
 "warnings": []}
```

---

## `POST /api/search`
Body: `{"query": "ACGT... or a path to a FASTA", "min_identity": 0.8}`
Returns `{"backend": "blast"|"minimap2", "hits": [{"segment": "ctg_1",
"identity": 0.99, "q_st": 0, "q_en": 500, "s_st": 100, "s_en": 600,
"strand": 1, "bitscore": 900}]}`

---

## `POST /api/export`
Body: `{"outdir": "plastr_out", "what": ["scaffolds","agp","gfa","csv","report","session"]}`
Returns `{"written": [{"kind": "scaffolds", "path": "...", "bytes": 1234}]}`.

## `GET /api/download/{kind}`
Streams a single artefact directly (`scaffolds`, `agp`, `gfa`, `csv`, `report`,
`session`). `Content-Disposition` is set so the browser saves it.

---

## `GET /api/session` / `POST /api/session`
Save/restore `{"layout": {...}, "settings": {...}, "plan": {...}}`. The server
stores `layout` and `settings` opaquely — they are the renderer's business.

---

## `GET /api/browse?path=/some/dir`
Directory listing for the built-in file picker:
`{"path": "/abs", "parent": "/", "entries": [{"name": "x.gfa", "is_dir": false, "size": 123}]}`
