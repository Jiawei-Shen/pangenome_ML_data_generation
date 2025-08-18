#!/usr/bin/env python3
import argparse
import math
from pathlib import Path

import numpy as np
from tqdm import tqdm
import pysam
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Fast AF distribution from VCF/BCF using pysam (INFO/AF by default)."
    )
    p.add_argument("vcf", help="Path to .vcf.gz/.vcf.bgz or .bcf")
    p.add_argument("--field", default="AF", help="INFO field to use (default: AF)")
    p.add_argument("--pass-only", action="store_true", help="Keep only FILTER=PASS")
    p.add_argument("--region", default=None,
                   help="Region like 'chr1:1-100000' (requires index)")
    p.add_argument("--max-af", type=float, default=None,
                   help="Keep AF <= value (e.g., 0.05 to focus on rare variants)")
    p.add_argument("--bins", type=int, default=100, help="Histogram bins (default 100)")
    p.add_argument("--logx", action="store_true",
                   help="Plot log10(AF). Zeros are omitted from log plot.")
    p.add_argument("--stream-hist", action="store_true",
                   help="Build histogram on the fly (lower memory, usually faster).")
    p.add_argument("--out-prefix", default="af_distribution",
                   help="Output file prefix (default: af_distribution)")
    return p.parse_args()


def af_iter(vf, field, pass_only, show_progress=False):
    it = vf.fetch()
    if show_progress:
        # no total known, so tqdm shows rate + count
        it = tqdm(it, unit="variants", mininterval=5, desc="Processing VCF")

    for rec in it:
        if pass_only and rec.filter.keys() not in (set(), {"PASS"}) and "PASS" not in rec.filter.keys():
            continue
        info_val = rec.info.get(field, None)
        if info_val is None:
            continue
        if isinstance(info_val, (tuple, list)):
            for v in info_val:
                if v is not None:
                    yield float(v)
        else:
            yield float(info_val)

def main(args):
    # Open once (pysam uses htslib BGZF; ensure the file is indexed if using --region)
    vf = pysam.VariantFile(args.vcf)

    out_prefix = Path(args.out_prefix)

    # Mode A: Streaming histogram (fast, memory-light)
    if args.stream_hist:
        # Choose edges. For linear AF, [0,1] makes sense. If focusing on rare, trim via --max-af.
        hi = args.max_af if (args.max_af is not None and args.max_af > 0) else 1.0
        if not args.logx:
            edges = np.linspace(0.0, hi, args.bins + 1, dtype=np.float64)
            counts = np.zeros(args.bins, dtype=np.int64)
            n_kept = 0
            for fv in af_iter(vf, args.field, args.pass_only, show_progress=True):
                if args.max_af is not None and fv > args.max_af:
                    continue
                # place into bin; clip at edges[-1]
                if 0.0 <= fv <= hi:
                    # np.searchsorted is ~ok; manual int calc faster for uniform bins
                    bin_idx = int((fv / hi) * args.bins)
                    if bin_idx == args.bins:  # fv==hi edge
                        bin_idx -= 1
                    counts[bin_idx] += 1
                    n_kept += 1
            # Save histogram counts & edges
            np.save(out_prefix.with_suffix(".hist_counts.npy"), counts)
            np.save(out_prefix.with_suffix(".hist_edges.npy"), edges)

            # Plot
            plt.figure()
            # For linear scale, we can pass left edges and weights to bar
            width = np.diff(edges)
            plt.bar(edges[:-1], counts, width=width, align="edge")
            plt.xlabel("Allele Frequency (AF)")
            plt.ylabel("Count (ALT alleles)")
            plt.title(f"AF Distribution (n={n_kept}, streaming)")
            plt.tight_layout()
            plt.savefig(out_prefix.with_suffix(".png"), dpi=180)

            # Save a CSV summary
            mids = (edges[:-1] + edges[1:]) / 2
            csvp = out_prefix.with_suffix(".histogram.csv")
            with open(csvp, "w") as f:
                f.write("bin_left,bin_right,bin_mid,count\n")
                for l, r, m, c in zip(edges[:-1], edges[1:], mids, counts):
                    f.write(f"{l},{r},{m},{int(c)}\n")

            print(f"[ok] Streaming hist saved: {csvp}, {out_prefix.with_suffix('.png')}")
            return

        else:
            # log10 histogram streaming: define edges in log-space, ignore zeros
            lo = -8  # 10^-8 lower bound for bins (adjust as needed)
            hi_log = math.log10(args.max_af) if args.max_af else 0.0  # up to 1 => 0
            edges = np.linspace(lo, hi_log, args.bins + 1, dtype=np.float64)
            counts = np.zeros(args.bins, dtype=np.int64)
            n_kept = 0
            for fv in af_iter(vf, args.field, args.pass_only):
                if fv <= 0.0:
                    continue  # cannot log
                if args.max_af is not None and fv > args.max_af:
                    continue
                x = math.log10(fv)
                if edges[0] <= x <= edges[-1]:
                    bin_idx = np.searchsorted(edges, x, side="right") - 1
                    if 0 <= bin_idx < counts.size:
                        counts[bin_idx] += 1
                        n_kept += 1
            # Save
            np.save(out_prefix.with_suffix(".hist_counts.npy"), counts)
            np.save(out_prefix.with_suffix(".hist_edges.npy"), edges)

            # Plot
            plt.figure()
            width = np.diff(edges)
            plt.bar(edges[:-1], counts, width=width, align="edge")
            plt.xlabel("log10(AF)")
            plt.ylabel("Count (ALT alleles)")
            plt.title(f"AF Distribution (log10, n={n_kept}, streaming)")
            plt.tight_layout()
            plt.savefig(out_prefix.with_suffix(".png"), dpi=180)

            csvp = out_prefix.with_suffix(".histogram.csv")
            with open(csvp, "w") as f:
                f.write("bin_left_log10,bin_right_log10,bin_mid_log10,count\n")
                mids = (edges[:-1] + edges[1:]) / 2
                for l, r, m, c in zip(edges[:-1], edges[1:], mids, counts):
                    f.write(f"{l},{r},{m},{int(c)}\n")

            print(f"[ok] Streaming log-hist saved: {csvp}, {out_prefix.with_suffix('.png')}")
            return

    # Mode B: Accumulate values (simpler; still fast for tens of millions of alleles)
    af_vals = []
    push = af_vals.append
    kept = 0
    for fv in af_iter(vf, args.field, args.pass_only):
        if args.max_af is not None and fv > args.max_af:
            continue
        push(fv)
        kept += 1

    if not af_vals:
        print("[warn] No AF values collected. Check field/filters/region.")
        return

    af = np.fromiter(af_vals, dtype=np.float64)
    # Save raw values (optional, large)
    np.savetxt(out_prefix.with_suffix(".af_values.csv"), af, delimiter=",", fmt="%.8g")

    # Summary
    q = np.quantile(af, [0, 0.25, 0.5, 0.75, 0.95, 0.99, 1.0])
    with open(out_prefix.with_suffix(".af_summary.csv"), "w") as f:
        f.write(f"count,{af.size}\n")
        f.write(f"mean,{float(af.mean())}\n")
        f.write(f"std,{float(af.std(ddof=0))}\n")
        f.write(f"min,{float(q[0])}\nq25,{float(q[1])}\nmedian,{float(q[2])}\n")
        f.write(f"q75,{float(q[3])}\nq95,{float(q[4])}\nq99,{float(q[5])}\nmax,{float(q[6])}\n")

    # Plot
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
    plt.title(f"AF Distribution (n={af.size})")
    plt.tight_layout()
    plt.savefig(out_prefix.with_suffix(".png"), dpi=180)
    print(f"[ok] Saved: {out_prefix.with_suffix('.png')}, "
          f"{out_prefix.with_suffix('.af_summary.csv')}")


if __name__ == "__main__":
    args = parse_args()
    main(args)
