#!/usr/bin/env python3

import argparse
import os
import shutil
import random
from pathlib import Path

def split_and_copy(src_dir, dest_dir, ratio, label):
    files = sorted(os.listdir(src_dir))
    random.shuffle(files)
    total = len(files)
    n_train = int(total * ratio[0])
    n_val = int(total * ratio[1])
    n_test = total - n_train - n_val

    splits = {
        "train": files[:n_train],
        "val": files[n_train:n_train + n_val],
        "test": files[n_train + n_val:]
    }

    for split, split_files in splits.items():
        split_path = os.path.join(dest_dir, split, label)
        os.makedirs(split_path, exist_ok=True)
        for fname in split_files:
            shutil.copy(os.path.join(src_dir, fname), os.path.join(split_path, fname))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, help="Path to source folder with true/ and false/")
    parser.add_argument("--output", required=True, help="Path to output folder")
    parser.add_argument("--organize", nargs=3, type=float, metavar=('TRAIN', 'VAL', 'TEST'),
                        help="Split ratio, e.g., --organize 0.6 0.2 0.2")

    args = parser.parse_args()
    if args.organize:
        assert abs(sum(args.organize) - 1.0) < 1e-6, "Split ratio must sum to 1.0"

        for label in ['true', 'false']:
            src_subdir = os.path.join(args.source, label)
            assert os.path.isdir(src_subdir), f"Missing directory: {src_subdir}"
            split_and_copy(src_subdir, args.output, args.organize, label)

        print(f"Data organized into {args.output}/train|val|test with true/false splits.")

if __name__ == "__main__":
    main()
