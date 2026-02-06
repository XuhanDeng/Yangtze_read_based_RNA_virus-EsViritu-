# Code Changes vs Original (Detailed)

This document explains exactly what was changed in this working copy, where it was changed, and why.
It is meant as a human-readable changelog so you can compare behavior with the original code.

## High-level summary
- Added explicit “step finished” logs so the pipeline progress is visible in real time.
- Optimized the clustering tie-break logic to avoid repeated full-table scans.
- Reduced per-contig I/O overhead in coverage/contig stats by batching work per process.
- Added progress logs inside `bam_coverage_windows`.
- Added concise docstrings to clarify purpose of helper functions.

---

## Detailed changes by file

### `src/EsViritu/EsViritu.py`

**1) High-level pipeline overview comment**
- **Where:** top of file, right after the shebang.
- **Why:** makes the full pipeline easier to understand at a glance.
- **Behavior impact:** none.

**2) Docstring for `esviritu()`**
- **Where:** inside `def esviritu():`
- **Why:** clarify the function’s role as the CLI entrypoint.
- **Behavior impact:** none.

**3) Step-completion log messages**
- **Where:** after each major stage in `esviritu()`:
  - `trim_filter`
  - `fastp_stats`
  - `various_readstats`
  - `minimap2_f` (initial, second, third)
  - `assembly_read_sharing_table` (initial, second)
  - `cluster_assemblies_by_read_sharing` (initial, second)
  - `clust_record_getter` (initial, second)
  - `bam_to_consensus_fasta`
  - `bam_to_coverm_table`
  - `bam_coverage_windows`
  - `read_ani_from_bam`
  - `assembly_table_maker`
  - `tax_profile`
  - report generation (reactable)
- **Why:** to show exactly where the pipeline is when it appears “stuck”.
- **Behavior impact:** logging only.

**4) Multi-process coverage windows**
- **Where:** call to `bam_coverage_windows`.
- **Change:**
  - Before: `bam_coverage_windows(third_map_bam)`
  - After:  `bam_coverage_windows(third_map_bam, int(args.CPU))`
- **Why:** allow `bam_coverage_windows` to use the requested CPU count.
- **Behavior impact:** same output, faster execution when `--cpu > 1`.

---

### `src/EsViritu/esv_funcs.py`

**1) `cluster_assemblies_by_read_sharing` performance fix**
- **Where:** inside the function’s sorting logic.
- **Problem in original:** the sort key called `df.filter(...)` for each assembly, scanning the entire
  dataframe multiple times inside the loop. With large tables, this is effectively non‑terminating.
- **Change:** precompute a dict `max_avg_identity[assembly]` in two passes:
  - Pass 1: update from `assembly_a` values.
  - Pass 2: update from `assembly_b` values.
- **Why two passes:** preserves the same tie‑break semantics as the original
  (`assembly_a` values are seen before `assembly_b` values).
- **Behavior impact:** identical output ordering, drastically faster.

**Added log lines**
- **Where:** start and end of `cluster_assemblies_by_read_sharing`.
- **Lines:**
  - `iam runging the new script`
  - `new script finished successfully`
- **Behavior impact:** logging only.

**2) `bam_to_coverm_table` batching optimization**
- **Where:** new helper functions and rewritten parallel loop.
- **New helpers:**
  - `_calculate_contig_stats_from_bam(...)`  
    same logic as the original `calculate_contig_stats`, but reuses a single BAM handle.
  - `_calculate_contig_stats_batch(...)`  
    opens BAM once and processes a list of contigs.
- **Why:** the original version opened the BAM once per contig, which is huge overhead at scale.
- **Behavior impact:** same output, much faster for many contigs.

**3) `bam_coverage_windows` multiprocessing + progress**
- **Where:** new `_coverage_windows_batch(...)` and `bam_coverage_windows(...)` signature change.
- **Change:**
  - Added batching with `ProcessPoolExecutor`.
  - Added `max_workers` parameter (default 1 for backwards compatibility).
  - Added progress logs every ~100 contigs.
- **Behavior impact:** same output, faster when `max_workers > 1`.

**4) Docstrings**
- **Where:** `fastp_stats`, `print_esviritu_banner`, `ghost_banner`.
- **Why:** improve readability only.
- **Behavior impact:** none.

---

### `src/EsViritu/utils/aniclust.py`
- Added docstrings for:
  - `parse_seqs`
  - `log_time`
  - `parse_arguments`
- **Behavior impact:** none.

---

### `src/EsViritu/utils/anicalc.py`
- Added docstrings for:
  - `parse_blast`
  - `yield_alignment_blocks`
  - `prune_alns`
  - `compute_ani`
  - `compute_cov`
  - `parse_arguments`
- **Behavior impact:** none.

---

## Output compatibility (important)
- All computational changes preserve input/output and data semantics.
- `cluster_assemblies_by_read_sharing` sorting produces the same results as before
  because the same tie-breaker values are used in the same order.
- `bam_coverage_windows` output remains identical; batch order is preserved.

## Usage reminders
- Use `-t` or `--cpu` to enable multi-process coverage windows:
  - Example: `EsViritu -t 4 ...`
- To run the modified repo code instead of a system install:
  - `PYTHONPATH=/path/to/EsViritu-main/src python -m EsViritu.EsViritu ...`

---

## Original vs changed code blocks

**Note:** snippets are abbreviated for readability. `...` means unchanged code omitted.

### `src/EsViritu/EsViritu.py` — pipeline overview + `esviritu()` docstring

**Original**
```python
#!/usr/bin/env python

import argparse
import time
import sys, os
...
def esviritu():

    pathname = os.path.dirname(__file__)
    esviritu_script_path = os.path.abspath(pathname)
    print(esviritu_script_path)
```

**Changed**
```python
#!/usr/bin/env python
# High-level pipeline overview (single-sample run):
# 1) Parse CLI args, validate DB paths, tools, and output/temp dirs; persist params YAML + log.
# 2) Optional read preprocessing via fastp (quality trim, filter, dedup), plus optional host filter.
# 3) Compute read stats and write a readstats YAML for downstream reporting.
# 4) Map reads to the virus reference DB with minimap2; exit early if no alignments.
# 5) Two rounds of read-sharing clustering to collapse similar references:
#    - build read-sharing tables, cluster assemblies, build clustered reference FASTA,
#      then re-map reads to the clustered references.
# 6) From the final BAM, generate consensus FASTA, coverage windows, read ANI per contig,
#    and coverm-like abundance metrics.
# 7) Write main output tables: per-contig detections, per-assembly summaries, and tax profile.
# 8) If R deps are present, render interactive HTML report; otherwise warn and skip.
# 9) Optionally clean temp files, then log total elapsed time.

import argparse
import time
import sys, os
...
def esviritu():
    """Run the EsViritu single-sample pipeline end-to-end from CLI args."""

    pathname = os.path.dirname(__file__)
    esviritu_script_path = os.path.abspath(pathname)
    print(esviritu_script_path)
```

### `src/EsViritu/EsViritu.py` — step-finished logs (pattern)

**Original**
```python
trim_filter_fn = timed_function(logger=logger)(esvf.trim_filter)
trim_filt_reads = trim_filter_fn(...)

logger.info(f"Main input reads: {trim_filt_reads}")
```

**Changed**
```python
trim_filter_fn = timed_function(logger=logger)(esvf.trim_filter)
trim_filt_reads = trim_filter_fn(...)
logger.info("STEP finished: trim_filter")

logger.info(f"Main input reads: {trim_filt_reads}")
```

> The same pattern is added after each major step listed earlier (fastp_stats, minimap2_f, etc).

### `src/EsViritu/EsViritu.py` — pass CPU to `bam_coverage_windows`

**Original**
```python
windows_cov_df = bam_coverage_windows_fn(
    third_map_bam
)
```

**Changed**
```python
windows_cov_df = bam_coverage_windows_fn(
    third_map_bam,
    int(args.CPU)
)
```

### `src/EsViritu/esv_funcs.py` — clustering sort performance fix

**Original**
```python
while remaining:
    # Sort remaining assemblies by read count (descending), breaking ties by highest avg_read_identity
    def get_avg_identity(a):
        # Try both possible columns for avg_read_identity
        val = None
        if 'avg_read_identity_a' in df.columns:
            # Get max value for this assembly as either a or b
            vals = df.filter((pl.col('assembly_a') == a) | (pl.col('assembly_b') == a))
            vals_a = vals.filter(pl.col('assembly_a') == a)
            vals_b = vals.filter(pl.col('assembly_b') == a)
            id_a = vals_a['avg_read_identity_a'].to_list() if 'avg_read_identity_a' in vals_a.columns else []
            id_b = vals_b['avg_read_identity_b'].to_list() if 'avg_read_identity_b' in vals_b.columns else []
            all_ids = [x for x in id_a + id_b if x is not None]
            if all_ids:
                val = max(all_ids)
        if val is None:
            return float('-inf')
        return val
    remaining.sort(key=lambda a: (total_reads.get(a, 0), get_avg_identity(a)), reverse=True)
```

**Changed**
```python
# Precompute per-assembly max avg_read_identity once (same semantics as prior get_avg_identity).
has_avg = 'avg_read_identity_a' in df.columns
max_avg_identity = {}
def _update_max_identity(asm, val):
    if val is None:
        return
    cur = max_avg_identity.get(asm)
    if cur is None or val > cur:
        max_avg_identity[asm] = val

# First pass: assembly_a identities
for row in df.iter_rows(named=True):
    ...
    if has_avg:
        _update_max_identity(row['assembly_a'], row.get('avg_read_identity_a'))
    ...

# Second pass: assembly_b identities (preserves original order)
if has_avg:
    for row in df.iter_rows(named=True):
        _update_max_identity(row['assembly_b'], row.get('avg_read_identity_b'))

while remaining:
    remaining.sort(
        key=lambda a: (total_reads.get(a, 0), max_avg_identity.get(a, float('-inf'))),
        reverse=True
    )
```

### `src/EsViritu/esv_funcs.py` — clustering log lines

**Original**
```python
def cluster_assemblies_by_read_sharing(df, threshold=0.33) -> pl.DataFrame:
    """
    Cluster assemblies that share at least `threshold` ...
    """
    ...
    return pl.DataFrame(records)
```

**Changed**
```python
def cluster_assemblies_by_read_sharing(df, threshold=0.33) -> pl.DataFrame:
    """
    Cluster assemblies that share at least `threshold` ...
    """
    logger.info("iam runging the new script")
    ...
    logger.info("new script finished successfully")
    return pl.DataFrame(records)
```

### `src/EsViritu/esv_funcs.py` — `bam_to_coverm_table` batching

**Original**
```python
def calculate_contig_stats(bam_path: str, contig: str, include_secondary: bool = False) -> Dict:
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        ...
        return {...}

def bam_to_coverm_table(...):
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        contigs_with_reads = set()
        for read in bamfile.fetch(until_eof=True):
            if not read.is_unmapped:
                contigs_with_reads.add(bamfile.get_reference_name(read.reference_id))

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for contig in contigs_with_reads:
                futures.append(executor.submit(calculate_contig_stats, bam_path, contig, include_secondary))
            records = []
            for future in futures:
                result = future.result()
                if result is not None:
                    records.append(result)
```

**Changed**
```python
def _calculate_contig_stats_from_bam(bamfile: pysam.AlignmentFile, contig: str,
                                     include_secondary: bool = False) -> Dict:
    ...
    return {...}

def calculate_contig_stats(bam_path: str, contig: str, include_secondary: bool = False) -> Dict:
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        return _calculate_contig_stats_from_bam(bamfile, contig, include_secondary)

def _calculate_contig_stats_batch(args) -> list:
    bam_path, contigs, include_secondary = args
    records = []
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        for contig in contigs:
            result = _calculate_contig_stats_from_bam(bamfile, contig, include_secondary)
            if result is not None:
                records.append(result)
    return records

def bam_to_coverm_table(...):
    with pysam.AlignmentFile(bam_path, "rb") as bamfile:
        contigs_with_reads = set()
        for read in bamfile.fetch(until_eof=True):
            if not read.is_unmapped:
                contigs_with_reads.add(bamfile.get_reference_name(read.reference_id))
        contigs_with_reads = list(contigs_with_reads)
        ...
        batches = [...]
        args_iter = [(bam_path, batch, include_secondary) for batch in batches]
        records = []
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            for batch_records in executor.map(_calculate_contig_stats_batch, args_iter):
                records.extend(batch_records)
```

### `src/EsViritu/esv_funcs.py` — `bam_coverage_windows` multiprocessing + progress

**Original**
```python
def bam_coverage_windows(bam_path: str) -> pl.DataFrame:
    bam = pysam.AlignmentFile(bam_path, "rb")
    contigs = list(bam.references)
    records = []
    for idx, contig in enumerate(contigs, start=1):
        if idx == 1 or idx % 100 == 0 or idx == len(contigs):
            logger.info(f"bam_coverage_windows progress: {idx}/{len(contigs)} contigs")
        ...
        records.append({...})
    bam.close()
    return pl.DataFrame(records)
```

**Changed**
```python
def _coverage_windows_batch(args) -> list:
    bam_path, contigs = args
    records = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for contig in contigs:
            ...
            records.append({...})
    return records

def bam_coverage_windows(bam_path: str, max_workers: int = 1) -> pl.DataFrame:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        contigs = list(bam.references)
    if max_workers == 1 or len(contigs) == 1:
        records = _coverage_windows_batch((bam_path, contigs))
        return pl.DataFrame(records)
    ...
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for batch, batch_records in zip(batches, executor.map(_coverage_windows_batch, args_iter)):
            records.extend(batch_records)
            processed += len(batch)
            if processed >= next_log:
                logger.info(f"bam_coverage_windows progress: {processed}/{total_contigs} contigs")
                next_log = processed + 100
    return pl.DataFrame(records)
```

### `src/EsViritu/esv_funcs.py` — docstrings (example)

**Original**
```python
def fastp_stats(...):
    pipeline_stats_fastp_html = ...
```

**Changed**
```python
def fastp_stats(...):
    """Collect total read counts from fastp JSON, generating stats if needed."""
    pipeline_stats_fastp_html = ...
```

### `src/EsViritu/utils/aniclust.py` — docstrings (example)

**Original**
```python
def parse_seqs(path):
    handle = gzip.open(path) if ...
```

**Changed**
```python
def parse_seqs(path):
    """Yield (id, sequence) tuples from a FASTA/FASTA.GZ file."""
    handle = gzip.open(path) if ...
```

### `src/EsViritu/utils/anicalc.py` — docstrings (example)

**Original**
```python
def compute_ani(alns):
    return round(sum(...), 2)
```

**Changed**
```python
def compute_ani(alns):
    """Compute average nucleotide identity across alignment blocks."""
    return round(sum(...), 2)
```

### `docs/code-changes.md` — new file

**Original**
```text
(file did not exist)
```

**Changed**
```text
docs/code-changes.md added to document all modifications.
```
