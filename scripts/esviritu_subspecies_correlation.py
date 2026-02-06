#!/usr/bin/env python3

import argparse
import csv
import multiprocessing as mp
import os
from typing import Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t


UNKNOWN_LABELS = {"unknow", "unknown", ""}


def _is_unknown(value: str) -> bool:
    v = str(value).strip().lower()
    return v in UNKNOWN_LABELS or v.startswith("unknown::") or v.startswith("unknow::")


def _parse_thresholds(text: str) -> List[float]:
    vals = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        vals.append(float(part))
    if not vals:
        raise SystemExit("No valid thresholds provided.")
    return sorted(set(vals))


def _label_subspecies(row, subspecies_col: str) -> str:
    subs = str(row[subspecies_col]) if subspecies_col in row else ""
    if _is_unknown(subs):
        asm = row.get("Assembly", "")
        if str(asm).strip() != "":
            return f"unknown::{asm}"
        return "unknown"
    return subs


def _aggregate_subspecies(df: pd.DataFrame, rpkmf_cols: List[str]) -> pd.DataFrame:
    if "subspecies" not in df.columns:
        raise SystemExit("No subspecies column found in input.")
    subs = df.copy()
    subs["subspecies"] = subs.apply(lambda r: _label_subspecies(r, "subspecies"), axis=1)
    agg = subs.groupby("subspecies", as_index=False)[rpkmf_cols].sum()
    return agg


_CTX = {}


def _init_worker(
    subspecies: List[str],
    ranks: np.ndarray,
    row_mean: np.ndarray,
    row_std: np.ndarray,
    thresholds: List[float],
    p_threshold: float,
    df_t: int,
):
    _CTX["subspecies"] = subspecies
    _CTX["ranks"] = ranks
    _CTX["row_mean"] = row_mean
    _CTX["row_std"] = row_std
    _CTX["thresholds"] = thresholds
    _CTX["min_thr"] = thresholds[0]
    _CTX["p_threshold"] = p_threshold
    _CTX["df_t"] = df_t


def _worker(payload: Tuple[List[int], str, str]):
    seed_indices, out_all_part, out_filt_part = payload
    subspecies = _CTX["subspecies"]
    ranks = _CTX["ranks"]
    row_mean = _CTX["row_mean"]
    row_std = _CTX["row_std"]
    thresholds = _CTX["thresholds"]
    min_thr = _CTX["min_thr"]
    p_threshold = _CTX["p_threshold"]
    df_t = _CTX["df_t"]

    n_rows, n_samples = ranks.shape
    with open(out_all_part, "w", newline="") as f_all, open(out_filt_part, "w", newline="") as f_filt:
        w_all = csv.writer(f_all, delimiter="\t")
        w_filt = csv.writer(f_filt, delimiter="\t")

        for seed_idx in seed_indices:
            seed = ranks[seed_idx]
            seed_mean = row_mean[seed_idx]
            seed_std = row_std[seed_idx]

            if seed_std == 0:
                r = np.zeros(n_rows, dtype=float)
                p = np.ones(n_rows, dtype=float)
            else:
                cov = ((ranks - row_mean[:, None]) * (seed - seed_mean)).sum(axis=1) / (
                    n_samples - 1
                )
                denom = row_std * seed_std
                r = np.zeros(n_rows, dtype=float)
                valid = denom != 0
                r[valid] = cov[valid] / denom[valid]

                r_clip = np.clip(r, -0.999999999, 0.999999999)
                t_stat = r_clip * np.sqrt(df_t / (1 - r_clip**2))
                p = 2 * t.sf(np.abs(t_stat), df_t)
                p[~valid] = 1.0

            seed_name = subspecies[seed_idx]
            for target_idx, target_name in enumerate(subspecies):
                if seed_idx == target_idx:
                    continue
                r_val = float(r[target_idx])
                p_val = float(p[target_idx])
                w_all.writerow([seed_name, target_name, r_val, p_val])

                if r_val >= min_thr and p_val <= p_threshold:
                    r_group = None
                    for thr in thresholds:
                        if r_val >= thr:
                            r_group = thr
                    w_filt.writerow([seed_name, target_name, r_val, p_val, r_group])


def _write_correlations(
    subspecies: List[str],
    data: np.ndarray,
    seed_mask: np.ndarray,
    rpkmf_cols: List[str],
    out_all: str,
    out_filtered: str,
    thresholds: Iterable[float],
    p_threshold: float,
    threads: int,
) -> None:
    n_rows, n_samples = data.shape
    if n_rows == 0 or n_samples < 3:
        with open(out_all, "w", newline="") as f_all, open(out_filtered, "w", newline="") as f_filt:
            csv.writer(f_all, delimiter="\t").writerow(["seed", "target", "r", "p"])
            csv.writer(f_filt, delimiter="\t").writerow(
                ["seed", "target", "r", "p", "r_group"]
            )
        return

    log_data = np.log1p(data)
    ranks = np.apply_along_axis(rankdata, 1, log_data)
    row_mean = ranks.mean(axis=1)
    row_std = ranks.std(axis=1, ddof=1)

    df_t = max(n_samples - 2, 1)
    thresholds = sorted(thresholds)
    seed_indices = np.where(seed_mask)[0].tolist()

    if not seed_indices:
        with open(out_all, "w", newline="") as f_all, open(out_filtered, "w", newline="") as f_filt:
            csv.writer(f_all, delimiter="\t").writerow(["seed", "target", "r", "p"])
            csv.writer(f_filt, delimiter="\t").writerow(
                ["seed", "target", "r", "p", "r_group"]
            )
        return

    if threads <= 1 or len(seed_indices) == 1:
        _init_worker(subspecies, ranks, row_mean, row_std, thresholds, p_threshold, df_t)
        _worker((seed_indices, out_all + ".part0", out_filtered + ".part0"))
        parts = [(out_all + ".part0", out_filtered + ".part0")]
    else:
        workers = min(threads, len(seed_indices))
        chunks = [seed_indices[i::workers] for i in range(workers)]
        parts = [
            (out_all + f".part{i}", out_filtered + f".part{i}") for i in range(workers)
        ]

        ctx = mp.get_context("fork") if hasattr(mp, "get_context") else mp
        with ctx.Pool(
            processes=workers,
            initializer=_init_worker,
            initargs=(subspecies, ranks, row_mean, row_std, thresholds, p_threshold, df_t),
        ) as pool:
            pool.map(_worker, [(chunk, p_all, p_filt) for chunk, (p_all, p_filt) in zip(chunks, parts)])

    with open(out_all, "w", newline="") as f_all, open(out_filtered, "w", newline="") as f_filt:
        w_all = csv.writer(f_all, delimiter="\t")
        w_filt = csv.writer(f_filt, delimiter="\t")
        w_all.writerow(["seed", "target", "r", "p"])
        w_filt.writerow(["seed", "target", "r", "p", "r_group"])

        for p_all, p_filt in parts:
            with open(p_all, "r", newline="") as f_part:
                f_all.writelines(f_part.readlines())
            with open(p_filt, "r", newline="") as f_part:
                f_filt.writelines(f_part.readlines())

    for p_all, p_filt in parts:
        try:
            os.remove(p_all)
            os.remove(p_filt)
        except OSError:
            pass


def main():
    parser = argparse.ArgumentParser(
        description="Subspecies correlations using log1p(RPKMF) + Spearman."
    )
    parser.add_argument("--input", required=True, help="Merged RPKMF TSV")
    parser.add_argument("--output-subspecies", required=True, help="Filtered subspecies table")
    parser.add_argument("--output-all", required=True, help="All correlations table")
    parser.add_argument(
        "--output-filtered",
        required=True,
        help="Filtered correlations table (r thresholds, includes p-values)",
    )
    parser.add_argument("--min-samples", type=int, default=10, help="Min samples present")
    parser.add_argument(
        "--thresholds",
        default="0.6,0.7,0.8,0.9",
        help="Comma-separated r thresholds",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.05,
        help="p-value threshold for filtered table (all table always includes p-values)",
    )
    parser.add_argument("--threads", type=int, default=1, help="Parallel workers")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    rpkmf_cols = [c for c in df.columns if c.endswith("_RPKMF")]
    if not rpkmf_cols:
        raise SystemExit("No RPKMF columns found.")

    for c in rpkmf_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    subs = _aggregate_subspecies(df, rpkmf_cols)
    presence = (subs[rpkmf_cols] > 0).sum(axis=1)
    subs = subs[presence >= args.min_samples].reset_index(drop=True)
    subs.to_csv(args.output_subspecies, sep="\t", index=False)

    subspecies = subs["subspecies"].astype(str).tolist()
    seed_mask = np.array([not _is_unknown(s) for s in subspecies], dtype=bool)
    data = subs[rpkmf_cols].to_numpy(dtype=float)

    thresholds = _parse_thresholds(args.thresholds)
    _write_correlations(
        subspecies,
        data,
        seed_mask,
        rpkmf_cols,
        args.output_all,
        args.output_filtered,
        thresholds,
        args.p_threshold,
        args.threads,
    )


if __name__ == "__main__":
    main()
