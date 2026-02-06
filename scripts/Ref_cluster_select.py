#!/usr/bin/env python3
"""
Select representative sequences from each cluster.

Rule 1: If a cluster contains sequences from the EsViritu DB, select one of those.
Rule 2: Otherwise, select a sequence from the earliest-priority input sample.
Rule 3: If all members are from the EsViritu DB, output all members (no collapse).

Inputs:
  - Cluster TSV: cluster_id<TAB>member1,member2,...
  - EsViritu DB sequence table (TSV/CSV) or FASTA with DB IDs

Output:
  - Representative sequence IDs (one per line by default)

This script is intentionally verbose and beginner-friendly. Comments explain:
  - why each helper exists,
  - what each variable means,
  - how each decision is made.
"""

import argparse
import gzip
import os
import re
import sys
from typing import Iterable, List, Optional, Tuple, Set


def eprint(message: str) -> None:
    # Log to stderr so status messages do not mix with main output.
    # This keeps the output file clean (only IDs) while still showing progress.
    print(message, file=sys.stderr)


def open_text(path: str):
    # Open plain text or gzip-compressed text transparently.
    # Many bioinformatics files are gzipped to save space; we handle both.
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def detect_fasta(path: str) -> bool:
    # Look at the first non-empty line; FASTA headers start with ">".
    # If we see ">", treat the file as FASTA, otherwise as a table.
    with open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            return stripped.startswith(">")
    return False


def load_fasta_ids(path: str) -> Set[str]:
    # Collect the sequence IDs from FASTA headers.
    # FASTA format: lines starting with ">" are headers; first token is the ID.
    ids = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                if seq_id:
                    ids.add(seq_id)
    return ids


def detect_delimiter(sample_line: str) -> str:
    # Heuristically detect TSV vs CSV; default to TSV if ambiguous.
    # We only need a simple guess, not a full CSV parser.
    if "\t" in sample_line:
        return "\t"
    if "," in sample_line:
        return ","
    return "\t"


def load_table_ids(path: str, id_col: Optional[str]) -> Set[str]:
    # Read a tabular file and return the ID set from a selected column.
    # This supports both:
    #   - files with a header row (column names),
    #   - files without headers (assume first column is ID).
    header_candidates = [
        "accession",
        "acc",
        "seq_id",
        "sequence_id",
        "id",
        "name",
    ]
    ids: Set[str] = set()

    with open_text(path) as handle:
        # Find the first non-empty line to detect delimiter and header.
        first_line = None
        for line in handle:
            if line.strip():
                first_line = line.rstrip("\n")
                break
        if first_line is None:
            return ids

        delimiter = detect_delimiter(first_line)
        cols = [c.strip() for c in first_line.split(delimiter)]
        lower_cols = [c.lower() for c in cols]
        is_header = any(c in header_candidates for c in lower_cols)

        # Decide which column contains the ID.
        # User can pass --db-id-col:
        #   - a column name (e.g., "Accession")
        #   - or a numeric index (0-based)
        if id_col:
            # If the user gave a numeric string, treat it as a 0-based index.
            if id_col.isdigit():
                id_idx = int(id_col)
            else:
                # If they gave a name, require a header and match case-insensitively.
                if not is_header:
                    raise ValueError(
                        "id column name provided but no header detected in table"
                    )
                if id_col.lower() not in lower_cols:
                    raise ValueError(f"id column '{id_col}' not found in table header")
                id_idx = lower_cols.index(id_col.lower())
        else:
            # Auto-detect: pick the first known header name, else fall back to column 0.
            if is_header:
                id_idx = None
                for candidate in header_candidates:
                    if candidate in lower_cols:
                        id_idx = lower_cols.index(candidate)
                        break
                if id_idx is None:
                    if len(cols) == 1:
                        id_idx = 0
                    else:
                        eprint(
                            "Warning: could not detect ID column, defaulting to first column"
                        )
                        id_idx = 0
            else:
                id_idx = 0

        def add_from_line(raw_line: str) -> None:
            # Extract one ID from a data line and add to the set.
            # We strip whitespace so IDs like "NC_000001.1 " still match.
            parts = raw_line.rstrip("\n").split(delimiter)
            if id_idx < len(parts):
                value = parts[id_idx].strip()
                if value:
                    ids.add(value)

        # If there is a header, skip it; otherwise, the first line is data.
        if is_header:
            for line in handle:
                if line.strip():
                    add_from_line(line)
        else:
            add_from_line(first_line)
            for line in handle:
                if line.strip():
                    add_from_line(line)

    return ids


def load_esviritu_ids(path: str, id_col: Optional[str]) -> Set[str]:
    # Auto-detect the DB format and load IDs accordingly.
    # This allows passing either a FASTA or a metadata table.
    if detect_fasta(path):
        eprint(f"Detected FASTA format for EsViritu DB: {path}")
        return load_fasta_ids(path)
    eprint(f"Detected table format for EsViritu DB: {path}")
    return load_table_ids(path, id_col)


def iter_clusters(path: str) -> Iterable[Tuple[str, List[str]]]:
    # Parse cluster TSV: first column is cluster ID, second is comma-separated members.
    # Example line:
    #   Cluster123<TAB>seqA,seqB,seqC
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                eprint(f"Warning: skipping malformed cluster line: {line}")
                continue
            cluster_id = parts[0].strip()
            members_raw = parts[1].strip()
            members = [m.strip() for m in members_raw.split(",") if m.strip()]
            yield cluster_id, members


def extract_sample_name(seq_id: str) -> Optional[str]:
    # Sample IDs are assumed to be the prefix before the last "_<number>".
    # Example: "SampleA_12" -> "SampleA"
    # If the ID does not match this pattern, return None.
    match = re.match(r"(.+)_\d+$", seq_id)
    if match:
        return match.group(1)
    return None


def infer_sample_order(cluster_path: str) -> List[str]:
    # Infer sample priority by first appearance of each sample in the cluster file.
    # This provides a stable, reproducible ordering without extra inputs.
    order: List[str] = []
    seen = set()
    for _, members in iter_clusters(cluster_path):
        for member in members:
            sample = extract_sample_name(member)
            if sample and sample not in seen:
                seen.add(sample)
                order.append(sample)
    return order


def load_sample_order_from_file(path: str) -> List[str]:
    # Read sample names from a file; allow commas or one-per-line.
    # This lets users control the priority explicitly.
    order: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "," in stripped:
                order.extend([s.strip() for s in stripped.split(",") if s.strip()])
            else:
                order.append(stripped)
    return order


def choose_representative(
    members: List[str],
    db_ids: Set[str],
    sample_order: Optional[List[str]],
) -> Tuple[Optional[str], str]:
    # Rule 1: prefer any member present in the EsViritu DB.
    # We return the first DB member we see.
    for member in members:
        if member in db_ids:
            return member, "db"

    # No members means nothing to select.
    if not members:
        return None, "empty"

    # Rule 2: if sample order is known, pick the earliest-priority sample.
    # We build a ranking dictionary for O(1) lookups.
    if sample_order:
        sample_rank = {sample: idx for idx, sample in enumerate(sample_order)}
        best_member = None
        best_idx = None
        for member in members:
            sample = extract_sample_name(member)
            if not sample:
                continue
            if sample not in sample_rank:
                continue
            idx = sample_rank[sample]
            if best_idx is None or idx < best_idx:
                best_idx = idx
                best_member = member
        if best_member:
            return best_member, "sample"

    # If no explicit sample order, pick the first member that looks like a sample contig.
    # This uses the "_<number>" suffix heuristic.
    for member in members:
        if extract_sample_name(member):
            return member, "sample"

    # Final fallback: just take the first member as a representative.
    # This handles cases where IDs are unusual and do not encode sample names.
    return members[0], "fallback"


def main() -> None:
    # Command-line interface definition.
    # Users specify where the cluster file and DB file are located.
    parser = argparse.ArgumentParser(
        description="Select representative sequences from clusters with EsViritu DB preference"
    )
    parser.add_argument(
        "--cluster",
        required=True,
        help="Cluster TSV (cluster_id<TAB>member1,member2,...)",
    )
    parser.add_argument(
        "--esviritu-db",
        required=True,
        help="EsViritu DB sequence table (TSV/CSV) or FASTA with IDs",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file with representative sequence IDs",
    )
    parser.add_argument(
        "--db-id-col",
        default=None,
        help="ID column name or 0-based index for DB table (ignored for FASTA)",
    )
    parser.add_argument(
        "--sample-order",
        default=None,
        help="Comma-separated sample names in priority order",
    )
    parser.add_argument(
        "--sample-order-file",
        default=None,
        help="File with sample names in priority order (one per line or comma-separated)",
    )
    parser.add_argument(
        "--with-cluster-id",
        action="store_true",
        help="Include cluster ID as first column in output",
    )

    args = parser.parse_args()

    # Basic input validation.
    # We exit early if required files are missing.
    if not os.path.exists(args.cluster):
        eprint(f"Cluster file not found: {args.cluster}")
        sys.exit(1)
    if not os.path.exists(args.esviritu_db):
        eprint(f"EsViritu DB file not found: {args.esviritu_db}")
        sys.exit(1)

    # Prevent ambiguous sample-order inputs.
    # We want exactly one source of ordering.
    if args.sample_order and args.sample_order_file:
        eprint("Provide only one of --sample-order or --sample-order-file")
        sys.exit(1)

    # Load or infer sample order.
    # Priority is: CLI string > file > inferred from cluster file.
    sample_order: Optional[List[str]] = None
    if args.sample_order:
        sample_order = [s.strip() for s in args.sample_order.split(",") if s.strip()]
    elif args.sample_order_file:
        sample_order = load_sample_order_from_file(args.sample_order_file)

    if not sample_order:
        # Infer order by scanning clusters and recording first appearance.
        sample_order = infer_sample_order(args.cluster)
        if sample_order:
            preview = ", ".join(sample_order[:10])
            suffix = "..." if len(sample_order) > 10 else ""
            eprint(f"Inferred sample order ({len(sample_order)}): {preview}{suffix}")
        else:
            sample_order = None
            eprint("No sample order inferred; will use member order fallback")

    # Load EsViritu DB IDs for fast membership checks.
    # We store them in a set for O(1) membership testing.
    db_ids = load_esviritu_ids(args.esviritu_db, args.db_id_col)
    eprint(f"Loaded {len(db_ids)} EsViritu DB IDs")

    # Ensure output directory exists if a directory is specified.
    # This avoids "No such file or directory" when writing output.
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Counters for a simple summary.
    total = 0
    all_db = 0
    from_db = 0
    from_sample = 0
    fallback = 0

    with open(args.output, "w") as out_handle:
        for cluster_id, members in iter_clusters(args.cluster):
            # Each line in the cluster file becomes one cluster in this loop.
            total += 1
            if not members:
                # Guard against empty clusters; we skip them.
                eprint(f"Warning: cluster '{cluster_id}' has no members")
                continue

            # New rule: if *all* members are DB sequences, keep them all.
            # In that case, we output every member as a representative.
            if all(member in db_ids for member in members):
                all_db += 1
                for member in members:
                    if args.with_cluster_id:
                        out_handle.write(f"{cluster_id}\t{member}\n")
                    else:
                        out_handle.write(f"{member}\n")
                continue

            # Otherwise pick a single representative using the previous rules.
            rep_id, reason = choose_representative(members, db_ids, sample_order)
            if rep_id is None:
                eprint(f"Warning: cluster '{cluster_id}' has no representative")
                continue
            if reason == "db":
                from_db += 1
            elif reason == "sample":
                from_sample += 1
            else:
                fallback += 1

            if args.with_cluster_id:
                out_handle.write(f"{cluster_id}\t{rep_id}\n")
            else:
                out_handle.write(f"{rep_id}\n")

    eprint(
        f"Processed {total} clusters: {all_db} all-DB clusters expanded, "
        f"{from_db} from DB, {from_sample} from samples, {fallback} fallback"
    )


if __name__ == "__main__":
    main()
