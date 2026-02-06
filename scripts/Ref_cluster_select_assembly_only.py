#!/usr/bin/env python3
"""
Select representative sequences from each cluster (assembly-only).

Rules:
  1) If sample order is provided, choose the first member from the earliest sample.
  2) Otherwise, choose the first member that looks like a sample contig (<sample>_<number>).
  3) Fallback to the first member.
"""

import argparse
import gzip
import os
import re
import sys
from typing import Iterable, List, Optional, Tuple


def eprint(message: str) -> None:
    print(message, file=sys.stderr)


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def iter_clusters(path: str) -> Iterable[Tuple[str, List[str]]]:
    # Cluster TSV: cluster_id<TAB>member1,member2,...
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
    match = re.match(r"(.+)_\d+$", seq_id)
    if match:
        return match.group(1)
    return None


def infer_sample_order(cluster_path: str) -> List[str]:
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
    sample_order: Optional[List[str]],
) -> Tuple[Optional[str], str]:
    if not members:
        return None, "empty"

    if sample_order:
        sample_rank = {sample: idx for idx, sample in enumerate(sample_order)}
        best_member = None
        best_idx = None
        for member in members:
            sample = extract_sample_name(member)
            if not sample or sample not in sample_rank:
                continue
            idx = sample_rank[sample]
            if best_idx is None or idx < best_idx:
                best_idx = idx
                best_member = member
        if best_member:
            return best_member, "sample"

    for member in members:
        if extract_sample_name(member):
            return member, "sample"

    return members[0], "fallback"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select representative sequences from clusters (assembly-only)"
    )
    parser.add_argument(
        "--cluster",
        required=True,
        help="Cluster TSV (cluster_id<TAB>member1,member2,...)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file with representative sequence IDs",
    )
    parser.add_argument(
        "--db-id-col",
        default=None,
        help="Ignored for assembly-only selection (kept for CLI compatibility)",
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

    if not os.path.exists(args.cluster):
        eprint(f"Cluster file not found: {args.cluster}")
        sys.exit(1)

    if args.sample_order and args.sample_order_file:
        eprint("Provide only one of --sample-order or --sample-order-file")
        sys.exit(1)

    sample_order: Optional[List[str]] = None
    if args.sample_order:
        sample_order = [s.strip() for s in args.sample_order.split(",") if s.strip()]
    elif args.sample_order_file:
        sample_order = load_sample_order_from_file(args.sample_order_file)

    if not sample_order:
        sample_order = infer_sample_order(args.cluster)
        if sample_order:
            preview = ", ".join(sample_order[:10])
            suffix = "..." if len(sample_order) > 10 else ""
            eprint(f"Inferred sample order ({len(sample_order)}): {preview}{suffix}")
        else:
            sample_order = None
            eprint("No sample order inferred; will use member order fallback")

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    total = 0
    from_sample = 0
    fallback = 0

    with open(args.output, "w") as out_handle:
        for cluster_id, members in iter_clusters(args.cluster):
            total += 1
            rep_id, reason = choose_representative(members, sample_order)
            if rep_id is None:
                eprint(f"Warning: cluster '{cluster_id}' has no representative")
                continue
            if reason == "sample":
                from_sample += 1
            else:
                fallback += 1

            if args.with_cluster_id:
                out_handle.write(f"{cluster_id}\t{rep_id}\n")
            else:
                out_handle.write(f"{rep_id}\n")

    eprint(
        f"Processed {total} clusters: {from_sample} from samples, {fallback} fallback"
    )


if __name__ == "__main__":
    main()
