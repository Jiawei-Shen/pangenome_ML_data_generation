#!/usr/bin/env python3
import argparse
import numpy as np
import glob, os, csv

# Pre-defined AF bins
AF_BINS = [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0]
AF_LABELS = [
    "0-1e-6",
    "1e-6-1e-5",
    "1e-5-1e-4",
    "1e-4-1e-3",
    "1e-3-1e-2",
    "1e-2-0.1",
    "0.1-0.5",
    "0.5-1.0"
]

def classify_bins(edges, counts):
    """
    Map fine-grained histogram bins into pre-defined super-bins.
    """
    out_counts = {lab: 0 for lab in AF_LABELS}
    mids = (edges[:-1] + edges[1:]) / 2
    for mid, c in zip(mids, counts):
        if c == 0:
            continue
        for i in range(len(AF_BINS)-1):
            if AF_BINS[i] < mid <= AF_BINS[i+1]:
                out_counts[AF_LABELS[i]] += int(c)
                break
    return out_counts

def main():
    ap = argparse.ArgumentParser(description="Re-bin AF histograms into 8 AF classes")
    ap.add_argument("--indir", required=True, help="Folder with *.hist_edges.npy and *.hist_counts.npy")
    ap.add_argument("--out", default="af_classes_summary.csv", help="Output CSV")
    args = ap.parse_args()

    total_counts = {lab: 0 for lab in AF_LABELS}

    for edge_file in glob.glob(os.path.join(args.indir, "*.hist_edges.npy")):
        base = edge_file.replace(".hist_edges.npy", "")
        count_file = base + ".hist_counts.npy"
        if not os.path.exists(count_file):
            continue

        edges = np.load(edge_file)
        counts = np.load(count_file)

        chr_counts = classify_bins(edges, counts)
        for lab in AF_LABELS:
            total_counts[lab] += chr_counts[lab]

    grand_total = sum(total_counts.values())

    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["AF_Class","Count","Proportion"])
        for lab in AF_LABELS:
            cnt = total_counts[lab]
            prop = cnt / grand_total if grand_total > 0 else 0
            writer.writerow([lab, cnt, f"{prop:.6f}"])

    print(f"[ok] Wrote summary -> {args.out}")

if __name__ == "__main__":
    main()
