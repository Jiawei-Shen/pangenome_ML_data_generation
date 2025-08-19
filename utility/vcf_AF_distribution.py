#!/usr/bin/env python3
"""
Fast AF distribution from VCF/BCF using pysam.

Features:
- Streams AF from INFO (default: AF)
- Handles multi-allelic sites
- Optional PASS-only filter
- Optional region fetch (requires .tbi/.csi)
- Variant-type filter: snp / indel / all
- Streaming histogram (low-mem) or accumulate mode
- Optional progress bar
- Optional threads (if pysam supports VariantFile.set_threads)

Outputs:
- Histogram PNG + CSV of bins
- (Optional) raw AF CSV
- (Accumulate mode) summary CSV (count/mean/std/quantiles)

Author: (revised for speed and correctness)
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pysam
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute AF distribution from VCF/BCF using pysam"
    )
    p.add_argument("vcf", help="Path to .vcf.gz/.vcf.bgz or .bcf")
    p.add_argument("--field", default="AF", help="INFO field to use (default: AF)")
    p.add_argument("--pass-only", action="store_true", help="Keep only FILTER=PASS")
    p.add_argument("--region", default=None,
                   help="Region like 'chr1:1-100000' (requires index)")
    p.add_argument("--variant-type", choices=["all", "snp", "indel"], default="all",
                   help="Variant type to include (default: all)")
    p.add_argument("--max-af", type=float, default=None,
                   help="Keep AF <= value (e.g., 0.05 to focus on rare variants)")
    p.add_argument("--bins", type=int, default=100, help="Histogram bins (default 100)")
    p.add_argument("--logx", action="store_true",
                   help="Plot on log10(AF); zeros are omitted from log plot")
    p.add_argument("--stream-hist", action="store_true",
                   help="Build histogram on the fly (lower memory, fastest)")
    p.add_argument("--out-prefix", default="af_distribution",
                   help="Output file prefix (default: af_distribution)")
    p.add_argument("--show-progress", action="store_true",
                   help="Display a tqdm progress bar while reading")
    p.add_argument("--threads", type=int, default=0,
                   help="htslib bgzf threads if supported (0=auto/disabled)")
    p.add_argument("--save-values", action="store_true",
                   help="Write raw AF values to CSV (slow/large for big files)")
    return p.parse_args()


def is_pass(rec) -> bool:
    """
    PASS definition for pysam:
      - Many VCFs treat PASS as having an empty filter set.
      - Some explicitly include 'PASS'.
    We accept either empty set OR containing 'PASS'.
    """
    keys = rec.filter.keys()
    # keys can be a set-like; empty => pass in many files
    return (len(keys) == 0) or ("PASS" in keys)


def is_snp_record(rec) -> bool:
    """True if REF and all ALT alleles are length 1."""
    if rec.alts is None:
        return False
    if len(rec.ref) != 1:
        return False
    return all(len(alt) == 1 for alt in rec.alts)


def is_indel_record(rec) -> bool:
    """True if any ALT has length different from REF (includes insertions/deletions)."""
    if rec.alts is None:
        return False
    rlen = len(rec.ref)
    return any(len(alt) != rlen for alt in rec.alts)


def af_iterator(vf: pysam.VariantFile, field: str, pass_only: bool,
                variant_type: str, region: str | None, show_progress: bool):
    """
    Yield AF values (float) for each ALT allele of records passing filters.
    """
    # choose iterator (region or whole file)
    it = vf.fetch(region) if region else vf.fetch()

    if show_progress:
        try:
            from tqdm import tqdm
            it = tqdm(it, unit="variants", mininterval=5, desc="Processing VCF")
        except Exception:
            # tqdm absent; proceed without progress
            pass

    for rec in it:
        # FILTER=PASS constraint
        if pass_only and not is_pass(rec):
            continue

        # Variant type filter
        if variant_type == "snp" and not is_snp_record(rec):
            continue
        if variant_type == "indel" and not is_indel_record(rec):
            continue

        # Extract AF(s)
        info_val = rec.info.get(field, None)
        if info_val is None:
            continue

        if isinstance(info_val, (tuple, list)):
            vals = info_val
        else:
            vals = (info_val,)

        for v in vals:
            if v is None:
                continue
            try:
                fv = float(v)
            except Exception:
                continue
            if math.isfinite(fv):
                yield fv


def save_histogram_and_plot(edges: np.ndarray, counts: np.ndarray, out_prefix: Path,
                            logx: bool, title_note: str = ""):
    """
    Save histogram to CSV and PNG plot.
    """
    # CSV
    csvp = out_prefix.with_suffix(".histogram.csv")
    mids = (edges[:-1] + edges[1:]) / 2
    with open(csvp, "w") as f:
        if logx:
            f.write("bin_left_log10,bin_right_log10,bin_mid_log10,count\n")
        else:
            f.write("bin_left,bin_right,bin_mid,count\n")
        for l, r, m, c in zip(edges[:-1], edges[1:], mids, counts):
            f.write(f"{l},{r},{m},{int(c)}\n")

    # NPY (optional binary artifacts for reuse)
    np.save(out_prefix.with_suffix(".hist_edges.npy"), edges)
    np.save(out_prefix.with_suffix(".hist_counts.npy"), counts)

    # Plot
    plt.figure()
    width = np.diff(edges)
    plt.bar(edges[:-1], counts, width=width, align="edge")
    plt.xlabel("log10(AF)" if logx else "Allele Frequency (AF)")
    plt.ylabel("Count (ALT alleles)")
    ttl = "AF Distribution"
    if title_note:
        ttl += f" — {title_note}"
    plt.title(ttl)
    plt.tight_layout()
    fig_path = out_prefix.with_suffix(".png")
    plt.savefig(fig_path, dpi=180)
    plt.close()


def main():
    args = parse_args()

    vcf_path = args.vcf
    out_prefix = Path(args.out_prefix)

    # Open once (pysam uses htslib BGZF)
    vf = pysam.VariantFile(vcf_path)

    # Try to enable multi-threaded bgzf if supported
    if args.threads and args.threads > 0:
        try:
            vf.set_threads(args.threads)
            print(f"[info] Using htslib threads: {args.threads}")
        except Exception:
            print("[warn] VariantFile.set_threads() not available; proceeding single-threaded.")

    # ── STREAMING HISTOGRAM MODE ───────────────────────────────────────────────
    if args.stream_hist:
        if args.logx:
            # log-space bins in log10(AF): ignore zeros
            lo = -8.0  # lower exponent, adjust if your AFs go smaller
            hi_log = math.log10(args.max_af) if (args.max_af and args.max_af > 0) else 0.0
            edges = np.linspace(lo, hi_log, args.bins + 1, dtype=np.float64)
            counts = np.zeros(args.bins, dtype=np.int64)
            n_kept = 0
            for fv in af_iterator(
                vf, args.field, args.pass_only, args.variant_type, args.region, args.show_progress
            ):
                if fv <= 0.0:
                    continue  # cannot log
                if args.max_af is not None and fv > args.max_af:
                    continue
                x = math.log10(fv)
                if edges[0] <= x <= edges[-1]:
                    # use searchsorted in log scale
                    idx = np.searchsorted(edges, x, side="right") - 1
                    if 0 <= idx < counts.size:
                        counts[idx] += 1
                        n_kept += 1

            note = f"log10, n={n_kept}, type={args.variant_type}"
            save_histogram_and_plot(edges, counts, out_prefix, logx=True, title_note=note)
            print(f"[ok] Streaming histogram (log10) saved for {vcf_path} as {out_prefix}.png/.csv")
            return

        else:
            # linear-space bins in [0, hi]
            hi = args.max_af if (args.max_af is not None and args.max_af > 0) else 1.0
            edges = np.linspace(0.0, hi, args.bins + 1, dtype=np.float64)
            counts = np.zeros(args.bins, dtype=np.int64)
            n_kept = 0
            for fv in af_iterator(
                vf, args.field, args.pass_only, args.variant_type, args.region, args.show_progress
            ):
                if args.max_af is not None and fv > args.max_af:
                    continue
                if 0.0 <= fv <= hi:
                    # fast uniform-binning:
                    idx = int((fv / hi) * args.bins)
                    if idx == args.bins:
                        idx -= 1
                    counts[idx] += 1
                    n_kept += 1

            note = f"n={n_kept}, type={args.variant_type}"
            save_histogram_and_plot(edges, counts, out_prefix, logx=False, title_note=note)
            print(f"[ok] Streaming histogram saved for {vcf_path} as {out_prefix}.png/.csv")
            return

    # ── ACCUMULATE MODE (also writes stats) ────────────────────────────────────
    af_vals = []
    push = af_vals.append
    for fv in af_iterator(
        vf, args.field, args.pass_only, args.variant_type, args.region, args.show_progress
    ):
        if args.max_af is not None and fv > args.max_af:
            continue
        push(fv)

    if len(af_vals) == 0:
        print("[warn] No AF values collected. Check field/filters/region/variant-type.")
        return

    af = np.fromiter(af_vals, dtype=np.float64)

    # Optional: write raw AFs (can be huge)
    if args.save_values:
        raw_csv = out_prefix.with_suffix(".af_values.csv")
        np.savetxt(raw_csv, af, delimiter=",", fmt="%.8g")
        print(f"[ok] Saved AF values -> {raw_csv}")

    # Summary stats
    q = np.quantile(af, [0, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0])
    summary_path = out_prefix.with_suffix(".af_summary.csv")
    with open(summary_path, "w") as f:
        f.write(f"count,{af.size}\n")
        f.write(f"mean,{float(af.mean())}\n")
        f.write(f"std,{float(af.std(ddof=0))}\n")
        f.write(f"min,{float(q[0])}\n")
        f.write(f"q25,{float(q[1])}\n")
        f.write(f"median,{float(q[2])}\n")
        f.write(f"q75,{float(q[3])}\n")
        f.write(f"q95,{float(q[4])}\n")
        f.write(f"q99,{float(q[5])}\n")
        f.write(f"max,{float(q[6])}\n")
    print(f"[ok] Saved summary stats -> {summary_path}")

    # Plot (accumulated)
    plt.figure()
    if args.logx:
        af_pos = af[af > 0]
        if af_pos.size:
            plt.hist(np.log10(af_pos), bins=args.bins)
            plt.xlabel("log10(AF)")
        else:
            plt.hist(af, bins=args.bins)
            plt.xlabel("Allele Frequency (AF)")
    else:
        plt.hist(af, bins=args.bins)
        plt.xlabel("Allele Frequency (AF)")

    plt.ylabel("Count (ALT alleles)")
    plt.title(f"AF Distribution (n={af.size}, type={args.variant_type})")
    plt.tight_layout()
    fig_path = out_prefix.with_suffix(".png")
    plt.savefig(fig_path, dpi=180)
    plt.close()
    print(f"[ok] Saved plot -> {fig_path}")


if __name__ == "__main__":
    main()
