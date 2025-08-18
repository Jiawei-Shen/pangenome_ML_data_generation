#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract AF from VCF with bcftools and plot distribution."
    )
    p.add_argument("vcf", help="Path to .vcf or .vcf.bgz")
    p.add_argument("--field", default="AF", help="INFO field to use (default: AF)")
    p.add_argument("--pass-only", action="store_true",
                   help="Use only FILTER=PASS records")
    p.add_argument("--max-af", type=float, default=None,
                   help="Keep AF <= this value (e.g., 0.05 to focus on rare variants)")
    p.add_argument("--bins", type=int, default=100, help="Histogram bins (default: 100)")
    p.add_argument("--logx", action="store_true",
                   help="Plot AF on log10 scale (for very small AFs)")
    p.add_argument("--out-prefix", default="af_distribution",
                   help="Output prefix for files (default: af_distribution)")
    return p.parse_args()


def bcftools_command(vcf_path, field, pass_only):
    # Build bcftools query command to stream AFs (one line per record).
    # If PASS-only, we filter with bcftools view -f PASS first.
    cmd = []
    if pass_only:
        cmd = ["bcftools", "view", "-f", "PASS", vcf_path]
        query = ["bcftools", "query", "-f", f"%INFO/{field}\n", "-"]
    else:
        query = ["bcftools", "query", "-f", f"%INFO/{field}\n", vcf_path]
    return cmd, query


def stream_afs(vcf_path, field, pass_only):
    # Returns a numpy array of AF values, flattened across multiallelic sites
    view_cmd, query_cmd = bcftools_command(vcf_path, field, pass_only)

    if view_cmd:
        view = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
        q = subprocess.Popen(query_cmd, stdin=view.stdout, stdout=subprocess.PIPE, text=True)
        view.stdout.close()
        stream = q.stdout
    else:
        q = subprocess.Popen(query_cmd, stdout=subprocess.PIPE, text=True)
        stream = q.stdout

    af_values = []
    for line in stream:
        s = line.strip()
        if not s or s == ".":
            continue
        # Multi-allelic AFs are comma-separated: 0.01,0.002
        for tok in s.split(","):
            tok = tok.strip()
            if not tok or tok == ".":
                continue
            try:
                val = float(tok)
                if math.isfinite(val):
                    af_values.append(val)
            except ValueError:
                # skip malformed
                continue

    stream.close()
    q.wait()
    return np.array(af_values, dtype=float)


def main():
    args = parse_args()
    vcf_path = args.vcf
    out_prefix = Path(args.out_prefix)

    print(f"[info] Reading AFs from: {vcf_path}")
    if args.pass_only:
        print("[info] Keeping only FILTER=PASS records")
    print(f"[info] INFO field: {args.field}")

    af = stream_afs(vcf_path, args.field, args.pass_only)
    if af.size == 0:
        print("[warn] No AF values found. Check INFO field name and filters.", file=sys.stderr)
        sys.exit(1)

    print(f"[info] Raw AF count (alleles): {af.size}")

    # Optional cap by AF
    if args.max_af is not None:
        keep = af <= args.max_af
        print(f"[info] Applying AF <= {args.max_af}: kept {keep.sum()} / {af.size}")
        af = af[keep]

    # Save raw AFs (one per ALT allele)
    csv_vals = out_prefix.with_suffix(".af_values.csv")
    np.savetxt(csv_vals, af, delimiter=",", fmt="%.8g")
    print(f"[ok] Saved AF values -> {csv_vals}")

    # Summary stats
    summary_path = out_prefix.with_suffix(".af_summary.csv")
    q = np.quantile(af, [0, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0])
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["count", af.size])
        w.writerow(["mean", float(af.mean())])
        w.writerow(["std", float(af.std(ddof=0))])
        w.writerow(["min", float(q[0])])
        w.writerow(["q25", float(q[1])])
        w.writerow(["median", float(q[2])])
        w.writerow(["q75", float(q[3])])
        w.writerow(["q95", float(q[4])])
        w.writerow(["q99", float(q[5])])
        w.writerow(["max", float(q[6])])
    print(f"[ok] Saved summary stats -> {summary_path}")

    # Plot
    fig_path = out_prefix.with_suffix(".png")

    plt.figure()
    if args.logx:
        # Avoid issues with zeros: keep strictly positive AFs for log scale
        af_pos = af[af > 0]
        if af_pos.size == 0:
            print("[warn] No positive AFs to plot on log scale; falling back to linear.")
            data = af
            xlabel = "Allele Frequency (AF)"
        else:
            data = np.log10(af_pos)
            xlabel = "log10(AF)"
        plt.hist(data, bins=args.bins)
    else:
        plt.hist(af, bins=args.bins)
        xlabel = "Allele Frequency (AF)"

    plt.xlabel(xlabel)
    plt.ylabel("Count (ALT alleles)")
    plt.title("AF Distribution")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=180)
    print(f"[ok] Saved histogram -> {fig_path}")


if __name__ == "__main__":
    main()
