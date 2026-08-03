"""Reference-based assembly evaluation, in the style of QUAST.

Given alignments of contigs to a reference we derive:

* **Misassemblies** -- adjacent alignment blocks of one contig that cannot be
  explained by a simple continuation of the reference. Classified as
  relocation, inversion, or translocation, following QUAST's definitions.
* **Genome fraction** -- how much of the reference is covered by any alignment.
* **Duplication ratio** -- aligned contig bases divided by covered reference
  bases; >1 means the assembly represents parts of the genome more than once.
* **NA50 / NGA50** -- N50 computed over *aligned blocks* rather than contigs,
  so a contig that is misassembled no longer gets credit for its full length.
* **Mismatch and indel rates** per 100 kb.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Sequence

from .align import Alignment, merge_collinear
from .metrics import nx_stat

# QUAST's defaults
EXTENSIVE_THRESHOLD = 1000  # ref gap above this is an extensive misassembly
LOCAL_THRESHOLD = 200  # between this and EXTENSIVE it is a local misassembly

RELOCATION = "relocation"
INVERSION = "inversion"
TRANSLOCATION = "translocation"
LOCAL = "local"


@dataclass(slots=True)
class Misassembly:
    contig: str
    kind: str
    is_extensive: bool
    # where in the contig the break sits
    contig_pos: int
    left_ref: str
    left_end: int
    right_ref: str
    right_start: int
    ref_gap: int | None
    description: str

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ReferenceReport:
    reference_length: int = 0
    reference_sequences: int = 0
    aligned_contigs: int = 0
    unaligned_contigs: int = 0
    partially_unaligned: int = 0
    fully_unaligned_length: int = 0
    genome_fraction: float = 0.0
    covered_bases: int = 0
    duplication_ratio: float = 0.0
    largest_alignment: int = 0
    total_aligned_length: int = 0
    na50: int = 0
    nga50: int | None = None
    la50: int = 0
    mismatches_per_100kb: float = 0.0
    indels_per_100kb: float = 0.0
    total_mismatches: int = 0
    total_indels: int = 0
    misassemblies: list[Misassembly] = field(default_factory=list)
    num_misassemblies: int = 0
    num_relocations: int = 0
    num_inversions: int = 0
    num_translocations: int = 0
    num_local_misassemblies: int = 0
    misassembled_contigs: list[str] = field(default_factory=list)
    misassembled_contigs_length: int = 0
    per_reference_coverage: dict[str, float] = field(default_factory=dict)
    coverage_blocks: dict[str, list[tuple[int, int]]] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["misassemblies"] = [m.to_dict() if hasattr(m, "to_dict") else m for m in self.misassemblies]
        return d


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= out[-1][1]:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(a, b) for a, b in out]


def circular_distance(gap: int, ref_length: int, circular: bool) -> int:
    """Distance between two reference points, going the short way round.

    Bacterial chromosomes and plasmids are circular, and a contig that spans the
    origin looks, on a linear reference, like a jump of nearly the whole
    replicon. Measuring the distance around the circle turns that back into the
    zero-length continuation it really is.
    """
    distance = abs(gap)
    if circular and ref_length > 0:
        distance = min(distance, abs(ref_length - distance))
    return distance


def classify_misassemblies(
    alignments: Sequence[Alignment],
    min_block: int = 200,
    circular_references: bool = True,
) -> list[Misassembly]:
    """Find misassembly breakpoints within each contig.

    Blocks are first merged where they are collinear, then the remaining
    adjacent pairs (ordered along the contig) are examined. A pair that changes
    reference sequence is a translocation; one that changes strand is an
    inversion; one that jumps more than 1 kb on the same reference is a
    relocation. Smaller jumps are recorded as local misassemblies, which QUAST
    treats as assembly-level indels rather than structural errors.
    """
    by_contig: dict[str, list[Alignment]] = {}
    for aln in alignments:
        if aln.q_span < min_block:
            continue
        by_contig.setdefault(aln.query, []).append(aln)

    results: list[Misassembly] = []
    for contig, hits in by_contig.items():
        # Merge only blocks that are genuinely contiguous -- the ones a small
        # indel split apart. A wider window would absorb the very discrepancies
        # we are here to find: at the placement default of 10 kb, every forward
        # relocation shorter than that becomes invisible.
        blocks = merge_collinear(hits, max_gap=LOCAL_THRESHOLD, max_overlap=LOCAL_THRESHOLD)
        blocks = [b for b in blocks if b.q_span >= min_block]
        if len(blocks) < 2:
            continue
        # Keep only blocks that do not substantially overlap each other on the
        # contig; a repeat aligning to many places should not look like dozens
        # of misassemblies.
        blocks.sort(key=lambda a: (a.q_st, -a.q_span))
        chosen: list[Alignment] = []
        for block in blocks:
            if chosen and block.q_st < chosen[-1].q_en - min_block // 2:
                if block.q_span > chosen[-1].q_span:
                    chosen[-1] = block
                continue
            chosen.append(block)
        if len(chosen) < 2:
            continue

        for left, right in zip(chosen, chosen[1:]):
            break_pos = (left.q_en + right.q_st) // 2
            if left.ref != right.ref:
                results.append(
                    Misassembly(
                        contig=contig,
                        kind=TRANSLOCATION,
                        is_extensive=True,
                        contig_pos=break_pos,
                        left_ref=left.ref,
                        left_end=left.r_en,
                        right_ref=right.ref,
                        right_start=right.r_st,
                        ref_gap=None,
                        description=(
                            f"contig continues on a different reference sequence "
                            f"({left.ref} -> {right.ref})"
                        ),
                    )
                )
                continue
            if left.strand != right.strand:
                results.append(
                    Misassembly(
                        contig=contig,
                        kind=INVERSION,
                        is_extensive=True,
                        contig_pos=break_pos,
                        left_ref=left.ref,
                        left_end=left.r_en,
                        right_ref=right.ref,
                        right_start=right.r_st,
                        ref_gap=None,
                        description="strand flips relative to the reference",
                    )
                )
                continue
            if left.strand > 0:
                gap = right.r_st - left.r_en
            else:
                gap = left.r_st - right.r_en
            distance = circular_distance(
                gap, max(left.r_len, right.r_len), circular_references
            )
            wrapped = circular_references and distance < abs(gap)
            if wrapped and distance <= EXTENSIVE_THRESHOLD:
                # A contig spanning the origin of a circular replicon is correct,
                # not misassembled.
                continue
            if distance > EXTENSIVE_THRESHOLD:
                results.append(
                    Misassembly(
                        contig=contig,
                        kind=RELOCATION,
                        is_extensive=True,
                        contig_pos=break_pos,
                        left_ref=left.ref,
                        left_end=left.r_en,
                        right_ref=right.ref,
                        right_start=right.r_st,
                        ref_gap=gap,
                        description=(
                            f"jumps {distance:,} bp along {left.ref} "
                            f"({'forward' if gap > 0 else 'backward'})"
                        ),
                    )
                )
            elif distance > LOCAL_THRESHOLD:
                results.append(
                    Misassembly(
                        contig=contig,
                        kind=LOCAL,
                        is_extensive=False,
                        contig_pos=break_pos,
                        left_ref=left.ref,
                        left_end=left.r_en,
                        right_ref=right.ref,
                        right_start=right.r_st,
                        ref_gap=gap,
                        description=f"local shift of {distance:,} bp",
                    )
                )
    results.sort(key=lambda m: (m.contig, m.contig_pos))
    return results


def evaluate_against_reference(
    alignments: Sequence[Alignment],
    reference_lengths: dict[str, int],
    contig_lengths: dict[str, int],
    genome_size: int | None = None,
    min_block: int = 200,
    primary_only: bool = True,
    circular_references: bool = True,
) -> ReferenceReport:
    """Build the full reference-based report.

    ``primary_only`` excludes secondary alignments. A repeat that maps to five
    places produces five alignments, and counting them all makes the assembly
    look longer than it is -- total aligned length can even exceed the assembly
    itself. Secondary hits are still worth keeping for the graph view, where
    seeing every place a repeat lands is the point, so they are filtered here
    rather than at alignment time.
    """
    report = ReferenceReport()
    report.reference_length = sum(reference_lengths.values())
    report.reference_sequences = len(reference_lengths)
    genome = genome_size or report.reference_length

    primary = [
        a
        for a in alignments
        if a.q_span >= min_block and (a.is_primary or not primary_only)
    ]

    # --- coverage of the reference ---
    per_ref: dict[str, list[tuple[int, int]]] = {}
    for aln in primary:
        per_ref.setdefault(aln.ref, []).append((aln.r_st, aln.r_en))
    covered = 0
    for ref, intervals in per_ref.items():
        merged = _merge_intervals(intervals)
        report.coverage_blocks[ref] = merged
        ref_covered = sum(b - a for a, b in merged)
        covered += ref_covered
        ref_len = reference_lengths.get(ref, 0)
        report.per_reference_coverage[ref] = (
            100.0 * ref_covered / ref_len if ref_len else 0.0
        )
    report.covered_bases = covered
    report.genome_fraction = 100.0 * covered / genome if genome else 0.0

    # --- aligned / unaligned contigs ---
    aligned_span: dict[str, int] = {}
    for aln in primary:
        aligned_span[aln.query] = aligned_span.get(aln.query, 0) + aln.q_span
    report.aligned_contigs = len(aligned_span)
    report.unaligned_contigs = len(contig_lengths) - len(aligned_span)
    report.fully_unaligned_length = sum(
        length for name, length in contig_lengths.items() if name not in aligned_span
    )
    report.partially_unaligned = sum(
        1
        for name, span in aligned_span.items()
        if contig_lengths.get(name, 0) and span < 0.9 * contig_lengths[name]
    )
    report.total_aligned_length = sum(aligned_span.values())
    report.duplication_ratio = (
        report.total_aligned_length / covered if covered else 0.0
    )

    # --- aligned block statistics (NA50 / NGA50) ---
    block_lengths: list[int] = []
    for contig, hits in _group(primary).items():
        # Aligned blocks are split at extensive misassemblies but not at local
        # ones, so a chimeric contig loses credit for its full length while a
        # contig with a small indel keeps it.
        for block in merge_collinear(
            hits, max_gap=EXTENSIVE_THRESHOLD, max_overlap=EXTENSIVE_THRESHOLD
        ):
            if block.q_span >= min_block:
                block_lengths.append(block.q_span)
    if block_lengths:
        report.largest_alignment = max(block_lengths)
        report.na50, report.la50 = nx_stat(block_lengths, 0.5)
        if genome:
            report.nga50, _ = nx_stat(block_lengths, 0.5, genome)

    # --- mismatch / indel rates ---
    total_mismatch = 0
    total_indel = 0
    aligned_bases = 0
    for aln in primary:
        total_mismatch += aln.mismatch_estimate()
        ins, dele, _bases = aln.indel_counts()
        total_indel += ins + dele
        aligned_bases += aln.block_len
    report.total_mismatches = total_mismatch
    report.total_indels = total_indel
    if aligned_bases:
        report.mismatches_per_100kb = 100_000.0 * total_mismatch / aligned_bases
        report.indels_per_100kb = 100_000.0 * total_indel / aligned_bases

    # --- misassemblies ---
    events = classify_misassemblies(
        primary, min_block=min_block, circular_references=circular_references
    )
    report.misassemblies = events
    report.num_relocations = sum(1 for e in events if e.kind == RELOCATION)
    report.num_inversions = sum(1 for e in events if e.kind == INVERSION)
    report.num_translocations = sum(1 for e in events if e.kind == TRANSLOCATION)
    report.num_local_misassemblies = sum(1 for e in events if e.kind == LOCAL)
    report.num_misassemblies = (
        report.num_relocations + report.num_inversions + report.num_translocations
    )
    bad = sorted({e.contig for e in events if e.is_extensive})
    report.misassembled_contigs = bad
    report.misassembled_contigs_length = sum(contig_lengths.get(c, 0) for c in bad)
    return report


def _group(alignments: Sequence[Alignment]) -> dict[str, list[Alignment]]:
    out: dict[str, list[Alignment]] = {}
    for aln in alignments:
        out.setdefault(aln.query, []).append(aln)
    return out
