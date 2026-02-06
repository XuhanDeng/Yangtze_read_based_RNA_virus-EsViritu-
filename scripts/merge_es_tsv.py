#!/usr/bin/env python3
# Merge original database metadata with new sequences (id + length).
#
# Inputs:
#   1) Original metadata TSV (with header, includes Accession and Length)
#   2) TSV with at least: <id>\t<length>
#   3) (Optional) ID list (one per line) to control which IDs are added
#   4) (Optional) FASTA for length fallback (used when an ID is missing in file 2)
# Output:
#   - A merged TSV that keeps all original rows and appends new rows for IDs
#     not already present. New rows are populated with non-NA placeholders.

from __future__ import annotations

import argparse
import csv
from typing import Dict, List, Tuple


def read_metadata(path: str) -> Tuple[List[str], List[List[str]], int]:
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows: List[List[str]] = []
        for row in reader:
            if not row:
                continue
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))
            rows.append(row)
    accession_idx = header.index("Accession") if "Accession" in header else 0
    return header, rows, accession_idx


def read_id_length(path: str) -> Dict[str, str]:
    pairs: Dict[str, str] = {}
    with open(path, "r") as handle:
        for i, line in enumerate(handle):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if i == 0:
                first = parts[0].strip().lower()
                second = parts[1].strip().lower() if len(parts) > 1 else ""
                if first in {"accession", "acc", "id", "seq_id", "sequence_id"} or "length" in second:
                    continue
            acc = parts[0].strip()
            length = parts[1].strip() if len(parts) > 1 else ""
            if acc:
                pairs[acc] = length
    return pairs


def read_id_list(path: str) -> List[str]:
    ids: List[str] = []
    with open(path, "r") as handle:
        for line in handle:
            acc = line.strip()
            if acc:
                ids.append(acc)
    return ids


def read_fasta_lengths(path: str) -> Dict[str, str]:
    lengths: Dict[str, str] = {}
    with open(path, "r") as handle:
        current_id = None
        current_len = 0
        for line in handle:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = str(current_len)
                current_id = line[1:].strip().split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_id is not None:
            lengths[current_id] = str(current_len)
    return lengths


def build_new_row(
    header: List[str],
    acc: str,
    length: str,
    assembly_existing: set,
) -> List[str]:
    row = [""] * len(header)
    name_to_idx: Dict[str, int] = {name: i for i, name in enumerate(header)}

    def set_if_present(col: str, value: str) -> None:
        idx = name_to_idx.get(col)
        if idx is not None:
            row[idx] = value

    set_if_present("Accession", acc)

    placeholders = {
        "description": f"Novel sequence {acc}",
        "Name": f"Novel_{acc}",
        "Segment": "unknown",
        "kingdom": "unknown",
        "phylum": "unknown",
        "tclass": "unknown",
        "order": "unknown",
        "family": "unknown",
        "genus": "unknown",
        "species": "unknown",
        "subspecies": "unknown",
        "TaxID": "0",
    }
    for col, value in placeholders.items():
        set_if_present(col, value)

    if length:
        set_if_present("Length", length)
        set_if_present("Asm_length", length)

    if "Assembly" in name_to_idx:
        base = f"acc:{acc}"
        assembly = base
        suffix = 1
        while assembly in assembly_existing:
            assembly = f"{base}_{suffix}"
            suffix += 1
        assembly_existing.add(assembly)
        set_if_present("Assembly", assembly)

    return row


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge metadata TSV with new sequence IDs (id + length)."
    )
    parser.add_argument("--metadata", required=True, help="Original metadata TSV")
    parser.add_argument("--id-length", required=True, help="TSV with <id>\\t<length>")
    parser.add_argument(
        "--ids",
        help="Optional list of IDs to add (one per line). If omitted, all IDs from --id-length are used.",
    )
    parser.add_argument(
        "--fasta",
        help="Optional FASTA to fill missing lengths when not present in --id-length.",
    )
    parser.add_argument("--output", required=True, help="Merged metadata TSV output")
    args = parser.parse_args()

    header, rows, accession_idx = read_metadata(args.metadata)
    existing_ids = {row[accession_idx] for row in rows if row[accession_idx]}

    length_idx = header.index("Length") if "Length" in header else None
    asm_length_idx = header.index("Asm_length") if "Asm_length" in header else None

    assembly_idx = header.index("Assembly") if "Assembly" in header else None
    assembly_existing = set()
    if assembly_idx is not None:
        assembly_existing = {
            row[assembly_idx] for row in rows if len(row) > assembly_idx and row[assembly_idx]
        }

    id_lengths = read_id_length(args.id_length)
    fasta_lengths: Dict[str, str] = {}
    if args.fasta:
        fasta_lengths = read_fasta_lengths(args.fasta)
    # Fill/overwrite Length and Asm_length for existing rows when provided.
    if length_idx is not None or asm_length_idx is not None:
        for row in rows:
            acc = row[accession_idx] if accession_idx < len(row) else ""
            if not acc:
                continue
            length = id_lengths.get(acc) or fasta_lengths.get(acc, "")
            if not length:
                continue
            if length_idx is not None:
                row[length_idx] = length
            if asm_length_idx is not None:
                row[asm_length_idx] = length
    if args.ids:
        ids_to_add = read_id_list(args.ids)
    else:
        ids_to_add = list(id_lengths.keys())
    added = 0
    for acc in ids_to_add:
        if acc in existing_ids:
            continue
        length = id_lengths.get(acc) or fasta_lengths.get(acc, "")
        new_row = build_new_row(header, acc, length, assembly_existing)
        rows.append(new_row)
        existing_ids.add(acc)
        added += 1

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)

    print(f"Added {added} new rows")


if __name__ == "__main__":
    main()
