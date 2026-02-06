#!/usr/bin/env python3
# Build metadata TSV for assembly-only sequences (no original EsViritu rows).
#
# Inputs:
#   1) TSV with at least: <id>\t<length>
#   2) (Optional) ID list (one per line) to control which IDs are added
#   3) (Optional) FASTA for length fallback
# Output:
#   - A TSV with EsViritu-style columns containing only the selected IDs.

from __future__ import annotations

import argparse
import csv
from typing import Dict, List


DEFAULT_HEADER = [
    "Accession",
    "description",
    "Name",
    "Segment",
    "kingdom",
    "phylum",
    "tclass",
    "order",
    "family",
    "genus",
    "species",
    "subspecies",
    "Length",
    "TaxID",
    "Assembly",
    "Asm_length",
]


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
        description="Create metadata TSV for assembly-only sequence IDs."
    )
    parser.add_argument("--id-length", required=True, help="TSV with <id>\\t<length>")
    parser.add_argument(
        "--ids",
        help="Optional list of IDs to include (one per line). If omitted, all IDs from --id-length are used.",
    )
    parser.add_argument(
        "--fasta",
        help="Optional FASTA to fill missing lengths when not present in --id-length.",
    )
    parser.add_argument("--output", required=True, help="Metadata TSV output")
    args = parser.parse_args()

    header = DEFAULT_HEADER
    id_lengths = read_id_length(args.id_length)
    fasta_lengths: Dict[str, str] = {}
    if args.fasta:
        fasta_lengths = read_fasta_lengths(args.fasta)

    if args.ids:
        ids_to_add = read_id_list(args.ids)
    else:
        ids_to_add = list(id_lengths.keys())

    assembly_existing = set()
    rows: List[List[str]] = []
    for acc in ids_to_add:
        length = id_lengths.get(acc) or fasta_lengths.get(acc, "")
        rows.append(build_new_row(header, acc, length, assembly_existing))

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows")


if __name__ == "__main__":
    main()
