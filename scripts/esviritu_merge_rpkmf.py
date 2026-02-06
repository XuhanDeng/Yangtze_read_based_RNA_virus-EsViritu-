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
        description="Merge per-sample EsViritu RPKMF tables into a wide matrix."
    )
    parser.add_argument("--inputs", nargs="+", required=True, help="Per-sample RPKMF TSVs")
    parser.add_argument("--output", required=True, help="Output merged TSV")
    args = parser.parse_args()

    merged = None
    meta_df = None
    meta_parts = []
    meta_cols_union = set()
    id_col_name = None

    for path in args.inputs:
        df = pd.read_csv(path, sep="\t", dtype=str)
        if len(df.columns) > 0:
            df = df[df[df.columns[0]] != df.columns[0]]

        id_col = _pick_column(df.columns, ["Assembly", "Accession"])
        if id_col is None:
            raise SystemExit(f"No Assembly/Accession column in {path}")

        if id_col_name is None:
            id_col_name = id_col
        elif id_col_name != id_col:
            # Normalize to the first-seen id column name
            df = df.rename(columns={id_col: id_col_name})

        id_cols = [id_col_name]

        sample_col = _pick_column(df.columns, ["sample_ID", "sample", "Sample"])
        if sample_col is None:
            raise SystemExit(f"No sample_ID column in {path}")

        sample_vals = df[sample_col].dropna().unique()
        sample = sample_vals[0] if len(sample_vals) > 0 else "sample"

        rpkmf_col = _pick_column(df.columns, ["RPKMF", "rpkmf"])
        if rpkmf_col is None:
            raise SystemExit(f"No RPKMF column in {path}")

        # Collect metadata rows across all files so we can fill blanks by assembly
        # Exclude per-sample metrics that should not be carried across samples
        per_sample_cols = {
            sample_col,
            rpkmf_col,
            "filtered_reads_in_sample",
            "read_count",
            "covered_bases",
            "avg_read_identity",
        }
        meta_cols = [c for c in df.columns if c not in per_sample_cols]
        meta_cols_union.update(meta_cols)
        meta_parts.append(df[meta_cols])

        per = df[id_cols + [rpkmf_col]].copy()
        per = per.rename(columns={rpkmf_col: f"{sample}_RPKMF"})

        if merged is None:
            merged = per
        else:
            merged = pd.merge(merged, per, on=id_cols, how="outer")

    merged = merged.fillna(0)

    if meta_parts:
        meta_all = pd.concat(meta_parts, ignore_index=True)
        # Ensure all metadata columns exist before aggregation
        for c in meta_cols_union:
            if c not in meta_all.columns:
                meta_all[c] = pd.NA

        def _first_nonempty(series):
            for v in series:
                if pd.notna(v) and str(v).strip() != "":
                    return v
            return pd.NA

        agg_cols = {c: _first_nonempty for c in meta_cols_union if c not in id_cols}
        meta_df = (
            meta_all.groupby(id_cols, as_index=False)
            .agg(agg_cols)
        )

        merged = pd.merge(meta_df, merged, on=id_cols, how="outer")

    merged.to_csv(args.output, sep="\t", index=False)

    def _default_out(suffix):
        if args.output.endswith(".tsv"):
            return args.output[:-4] + suffix + ".tsv"
        return args.output + suffix + ".tsv"

    def _aggregate_by(col_name, out_path):
        if col_name not in merged.columns:
            return
        df = merged.copy()
        df[col_name] = df[col_name].astype(str)
        val = df[col_name].str.strip()
        is_unknown = val.eq("") | val.str.lower().isin(["unknow", "unknown"])
        rpkmf_cols = [c for c in df.columns if c.endswith("_RPKMF")]
        for c in rpkmf_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
        # Known values: cluster by species/subspecies
        known = df[~is_unknown]
        agg = known.groupby(col_name, as_index=False)[rpkmf_cols].sum()

        # Unknown values: keep each row, add Assembly to distinguish
        unknown = df[is_unknown]
        if not unknown.empty and "Assembly" in unknown.columns:
            unknown = unknown.copy()
            unknown[col_name] = "unknown::" + unknown["Assembly"].astype(str)
            unknown = unknown[[col_name] + rpkmf_cols]
            agg = pd.concat([agg, unknown], ignore_index=True)

        agg.to_csv(out_path, sep="\t", index=False)

    # Always write species/subspecies aggregated outputs
    _aggregate_by("species", _default_out(".species"))
    _aggregate_by("subspecies", _default_out(".subspecies"))


if __name__ == "__main__":
    main()
