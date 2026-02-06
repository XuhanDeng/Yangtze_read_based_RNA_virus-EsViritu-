# try_change_esv

This repository contains a focused set of changes to EsViritu for better
progress visibility and faster performance on large datasets.

## What changed
- Added explicit step-completion logs in the main pipeline.
- Optimized clustering tie-break computation to avoid repeated full-table scans.
- Batched per-contig BAM processing to reduce I/O overhead.
- Added multiprocessing for `bam_coverage_windows` and progress logs.
- Added documentation of all changes.

## Files included
- `src/EsViritu/EsViritu.py`
- `src/EsViritu/esv_funcs.py`
- `docs/code-changes.md`

## Details
See `docs/code-changes.md` for a full, detailed, side-by-side explanation of
original vs changed code.

## EsViritu temp outputs (per sample)
Temporary files are written under `{OUTPUT_DIR}/{SAMPLE}_temp` unless `--temp` is set.

### Read preprocessing
- `{sample}.filtered.fastq`: reads after optional `fastp` quality trimming/filtering. If host/spike-in
  filtering is enabled, this becomes the post-filtered read set used for mapping.
- `{sample}.filtered_R1.fastq` / `{sample}.filtered_R2.fastq`: paired read outputs from host/spike-in
  filtering (if enabled).
- `{sample}.fastp.json` / `{sample}.fastp.html`: `fastp` QC reports for the preprocessing step.

### Mapping and clustering
- `{sample}.initial.filt.bam`: mapped reads passing alignment thresholds (unsorted).
- `{sample}.initial.filt.sorted.bam`: sorted + indexed BAM used as the input to round-1 clustering.
- `{sample}.initial_read_comp_compare.tsv`: read-sharing comparison table for round-1 clustering.
- `{sample}.initial_read_comp_clust.tsv`: read-sharing cluster assignments for round-1.
- `{sample}_clustered_refs.fasta`: clustered reference FASTA after round-1 clustering.

- `{sample}.second.filt.bam`: mapped reads passing alignment thresholds (unsorted) against round-1 refs.
- `{sample}.second.filt.sorted.bam`: sorted + indexed BAM used as the input to round-2 clustering.
- `{sample}.second_read_comp_compare.tsv`: read-sharing comparison table for round-2 clustering.
- `{sample}.second_read_comp_clust.tsv`: read-sharing cluster assignments for round-2.
- `{sample}_second_clustered_refs.fasta`: clustered reference FASTA after round-2 clustering.

- `{sample}.third.filt.bam`: mapped reads passing alignment thresholds (unsorted) against round-2 refs.
- `{sample}.third.filt.sorted.bam`: final sorted + indexed BAM used for all downstream stats.
- `{sample}.third.filt.unmapped.bam`: reads that did not align to the round-2 refs (final map only).
- `{sample}.third.filt.unmapped.sorted.bam`: sorted + indexed unmapped reads from the final map.
- `{sample}.third.filt.unmapped.fastq`: FASTQ of unmapped reads from the final map.
- `{sample}.third.filt.failed.bam`: reads that mapped but failed alignment thresholds (final map only).
- `{sample}.third.filt.failed.sorted.bam`: sorted + indexed failed-threshold reads from the final map.
- `{sample}.third.filt.failed.fastq`: FASTQ of failed-threshold reads from the final map.

### How the three filtered BAMs differ
- `initial.filt.sorted.bam`: first mapping to the full virus DB; filtered by alignment thresholds; sorted.
- `second.filt.bam`: second mapping to round-1 clustered references; filtered by the same thresholds; unsorted.
- `third.filt.sorted.bam`: third mapping to round-2 clustered references; filtered by the same thresholds; sorted; used for final stats and reports.

Generation steps (simplified):
1) `fastp` → filtered FASTQ (optional).
2) `minimap2_f` → `initial.filt.bam` → sort/index → `initial.filt.sorted.bam`.
3) Cluster refs (round 1) → `minimap2_f` → `second.filt.bam` → sort/index → `second.filt.sorted.bam`.
4) Cluster refs (round 2) → `minimap2_f` → `third.filt.bam` → sort/index → `third.filt.sorted.bam`.

### Downstream intermediate tables
- `{sample}_final_coverm.tsv`: coverm-like abundance table computed from final BAM.
- `{sample}.read_ani.per_contig.tsv`: per-contig read ANI table.
