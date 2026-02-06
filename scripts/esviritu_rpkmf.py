#!/usr/bin/env python3

import argparse
import pandas as pd


def _pick_column(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Compute RPKMF for EsViritu assembly_summary merged table."
    )
    parser.add_argument("--input", required=True, help="Merged assembly_summary TSV")
    parser.add_argument("--output", required=True, help="Output TSV with RPKMF")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    # Drop duplicated header rows from concatenation
    if len(df.columns) > 0:
        df = df[df[df.columns[0]] != df.columns[0]]

    filtered_col = _pick_column(df.columns, ["filtered_reads_in_sample", "filtered_reads"])
    length_col = _pick_column(
        df.columns,
        ["Asm_length", "Asm_Length", "Assembly_length", "assembly_length", "Length", "Contig_length"],
    )
    read_col = _pick_column(df.columns, ["read_count", "Read_count", "reads", "Read_count"])

    if filtered_col is None or length_col is None or read_col is None:
        missing = [n for n, v in [("filtered", filtered_col), ("length", length_col), ("read_count", read_col)] if v is None]
        raise SystemExit(f"Missing required columns: {', '.join(missing)}")

    # Use filtered_reads_in_sample from the first data row for all rows
    filtered_val = pd.to_numeric(df[filtered_col].dropna().head(1), errors="coerce")
    if filtered_val.empty or filtered_val.iloc[0] is None:
        filtered_reads = 0.0
    else:
        filtered_reads = float(filtered_val.iloc[0])

    df[filtered_col] = filtered_reads

    length = pd.to_numeric(df[length_col], errors="coerce").fillna(0.0)
    reads = pd.to_numeric(df[read_col], errors="coerce").fillna(0.0)

    if filtered_reads == 0:
        rpkmf = reads * 0.0
    else:
        rpkmf = (reads / (length / 1000.0)) / (filtered_reads / 1e6)
        rpkmf = rpkmf.where(length > 0, 0.0)

    if "RPKMF" in df.columns:
        df["RPKMF"] = rpkmf.fillna(0.0)
    else:
        df["RPKMF"] = rpkmf.fillna(0.0)

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
