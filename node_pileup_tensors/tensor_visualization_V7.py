#!/usr/bin/env python3
"""
Visualize head N shard_XXXXX_data.npy files produced by the pileup/tensor script.

For each selected .npy shard:
  - Load the array of shape (num_variants, 6, H, W)
  - For the first K variants in that shard, save a PNG with 6 heatmaps
    (one per channel: base, BQ, mismatch flags, MAPQ, CIGAR, AF-bin)

Usage example:

  python visualize_pileup_npy.py \
      --data-dir /scratch/jshen/data/Pansoma/COLO829T_ONT/6ch_training_data_P100_SNV_flat \
      --max-files 2 \
      --samples-per-file 3 \
      --output-dir ./npy_previews

"""

import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")  # so it works on headless servers
import matplotlib.pyplot as plt


# Channel semantics (matching your generator script)
CHANNEL_NAMES = [
    "ch1: base encoding",
    "ch2: base quality",
    "ch3: mismatch flags",
    "ch4: mapping quality",
    "ch5: CIGAR encoding",
    "ch6: AF bin"
]


def log(msg: str) -> None:
    print(msg, flush=True)


def find_npy_files(data_dir):
    """Return sorted list of *_data.npy files in the directory."""
    if not os.path.isdir(data_dir):
        raise SystemExit(f"Error: data_dir {data_dir!r} is not a directory")

    files = [f for f in os.listdir(data_dir) if f.endswith(".npy")]
    # Optional: prioritize shard_XXXXX_data.npy
    files.sort()
    return [os.path.join(data_dir, f) for f in files]


def visualize_tensor(sample, out_path_prefix):
    """
    Given a single tensor of shape (6, H, W), save a PNG with 6 heatmaps.

    out_path_prefix: base path (without extension) to save figure, e.g. ".../shard00000_idx000"
    """
    if sample.ndim != 3 or sample.shape[0] != 6:
        raise ValueError(f"Expected tensor shape (6, H, W), got {sample.shape}")

    _, H, W = sample.shape

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    axes = axes.ravel()

    for ch in range(6):
        ax = axes[ch]
        im = ax.imshow(sample[ch], aspect="auto", interpolation="nearest")
        ax.set_title(CHANNEL_NAMES[ch])
        ax.set_xlabel("Window position (W)")
        ax.set_ylabel("Row (H)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(f"Tensor preview (shape={sample.shape})", fontsize=14)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    out_png = out_path_prefix + ".png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    log(f"  Saved visualization -> {out_png}")


def main():
    ap = argparse.ArgumentParser(
        description="Visualize head N NPY shard files produced by pileup tensor generator."
    )
    ap.add_argument(
        "--data-dir",
        required=True,
        help="Directory containing shard_XXXXX_data.npy (or similar) files."
    )
    ap.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write PNG previews."
    )
    ap.add_argument(
        "--max-files",
        type=int,
        default=3,
        help="Visualize only the first N .npy files in data-dir."
    )
    ap.add_argument(
        "--samples-per-file",
        type=int,
        default=3,
        help="Number of tensors (variants) to visualize per file."
    )
    ap.add_argument(
        "--verbose",
        action="store_true",
        help="Print extra info about shapes."
    )

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    all_npy = find_npy_files(args.data_dir)
    if not all_npy:
        raise SystemExit(f"No .npy files found in {args.data_dir!r}")

    selected = all_npy[: max(1, args.max_files)]

    log(f"Found {len(all_npy)} .npy files, visualizing first {len(selected)} of them:")
    for path in selected:
        log(f"  {path}")

    for file_idx, npy_path in enumerate(selected):
        log(f"\n=== File {file_idx} / {len(selected)-1}: {npy_path} ===")

        arr = np.load(npy_path)
        if arr.ndim != 4 or arr.shape[1] != 6:
            log(f"Warning: expected array shape (N, 6, H, W), got {arr.shape}; skipping file.")
            continue

        num_variants, channels, H, W = arr.shape
        if args.verbose:
            log(f"  Shape: {arr.shape}  (variants={num_variants}, channels={channels}, H={H}, W={W})")

        num_to_show = min(args.samples_per_file, num_variants)
        log(f"  Visualizing first {num_to_show} variants from this shard...")

        base_name = os.path.splitext(os.path.basename(npy_path))[0]

        for i in range(num_to_show):
            sample = arr[i]  # (6, H, W)
            out_prefix = os.path.join(
                args.output_dir,
                f"{base_name}_var{i:04d}"
            )
            visualize_tensor(sample, out_prefix)

    log("\nDone.")


if __name__ == "__main__":
    main()
