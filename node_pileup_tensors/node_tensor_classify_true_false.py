#!/usr/bin/env python3
"""tensor_vcf_classifier_parallel.py

Speed‑optimized version that classifies *.npy tensor files as TRUE/FALSE variants.
Key accelerators
----------------
1.  **ProcessPoolExecutor** – one worker per node directory (default: min(8, CPU cores)).
2.  **Hard‑link instead of copy** when source & destination are on the same filesystem.
3.  **Pickled VCF cache** – the first run stores a *.pickle containing the parsed variant set;
    subsequent runs reload it instantly.
4.  **Progress message every 1 000 tensors processed** (aggregated across workers).

All previous functionality (summary.json output, folder layout, etc.) is preserved.
"""

from __future__ import annotations
import argparse
import concurrent.futures as cf
import json
import os
import pickle
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple, Dict

# -----------------------------------------------------------------------------
# Helper: bcftools query (with pickle cache)
# -----------------------------------------------------------------------------

def load_vcf_variants(vcf_path: Path, chrom: str) -> set[Tuple[int, str, str]]:
    """Return a set of (pos, REF, ALT) tuples for *chrom*.
    Uses a .pickle cache alongside the VCF for speed on re‑runs."""
    cache_path = vcf_path.with_suffix(vcf_path.suffix + f".{chrom}.pickle")
    if cache_path.is_file():
        print(f"[VCF] Loading cached variants from {cache_path} …")
        with cache_path.open("rb") as fh:
            return pickle.load(fh)

    print(f"[VCF] Extracting variants for {chrom} via bcftools …")
    variants: set[Tuple[int, str, str]] = set()
    cmd = ["bcftools", "view", "-r", chrom, str(vcf_path)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"bcftools failed ({e.returncode}): {e.stderr}")

    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        pos = int(f[1])
        ref = f[3].upper()
        for alt in f[4].upper().split(','):
            variants.add((pos, ref, alt))
    print(f"[VCF] Parsed {len(variants):,} variants; caching …")
    with cache_path.open("wb") as fh:
        pickle.dump(variants, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return variants

# -----------------------------------------------------------------------------
# Helper: load node‑>GRCh38 start positions
# -----------------------------------------------------------------------------

def load_node_starts(json_path: Path) -> Dict[str, int]:
    print(f"[JSON] Loading node positions from {json_path} …")
    with json_path.open() as fh:
        data = json.load(fh)
    result: Dict[str, int] = {}
    for entry in data.get("nodes", []):
        nid = str(entry.get("node_id"))
        start = entry.get("grch38_position_start")
        if nid and isinstance(start, int):
            result[nid] = start
    print(f"[JSON] Loaded {len(result):,} node positions.")
    return result

# -----------------------------------------------------------------------------
# IO helper (hard‑link or copy)
# -----------------------------------------------------------------------------

def link_or_copy(src: Path, dst: Path) -> None:
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)

# -----------------------------------------------------------------------------
# Per‑node worker
# -----------------------------------------------------------------------------

def process_node(args: Tuple[str, Path, Dict[str, int], set, Path, Path]) -> Tuple[List[dict], int, int, int]:
    node_id, node_dir, node_starts, vcf_variants, true_dir, false_dir = args
    summary_path = node_dir / "variant_summary.json"
    try:
        with summary_path.open() as fh:
            summary_json = json.load(fh)
    except Exception:
        return [], 0, 0, 0

    start_pos = node_starts.get(node_id)
    if start_pos is None:
        return [], 0, 0, 0

    recs: List[dict] = []
    t_ct = f_ct = processed = 0

    for variant in summary_json.get("variants_passing_af_filter", []):
        tensor_file = variant.get("tensor_file")
        variant_key = variant.get("variant_key")
        if not tensor_file or not variant_key or not tensor_file.endswith(".npy"):
            continue
        tensor_path = node_dir / tensor_file
        if not tensor_path.is_file():
            continue
        try:
            off, _, ref, alt = variant_key.split("_", 3)
            off = int(off)
            ref = ref.upper()
            alt = alt.upper()
        except Exception:
            continue

        gpos = start_pos + off + 1
        is_true = (gpos, ref, alt) in vcf_variants
        dest_dir = true_dir if is_true else false_dir
        dest_path = dest_dir / f"{node_id}_{tensor_file}"
        link_or_copy(tensor_path, dest_path)

        recs.append({
            "node_id": node_id,
            "tensor_file": tensor_file,
            "variant_key": variant_key,
            "genomic_position": gpos,
            "ref": ref,
            "alt": alt,
            "classification": "true" if is_true else "false",
        })
        processed += 1
        if is_true:
            t_ct += 1
        else:
            f_ct += 1

    return recs, t_ct, f_ct, processed

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Classify tensor .npy files as true/false using a VCF ground truth (parallel).")
    ap.add_argument("tensor_folder_path", type=Path)
    ap.add_argument("vcf_file", type=Path)
    ap.add_argument("node_pos_json", type=Path)
    ap.add_argument("--output_folder", type=Path, default="./classification_results")
    ap.add_argument("--chr", default="chr1")
    ap.add_argument("--workers", type=int, default=min(os.cpu_count() or 1, 8),
                    help="Number of parallel workers (default: min(8, CPU cores))")
    args = ap.parse_args()

    args.output_folder.mkdir(parents=True, exist_ok=True)
    true_dir = args.output_folder / "true"
    false_dir = args.output_folder / "false"
    true_dir.mkdir(exist_ok=True)
    false_dir.mkdir(exist_ok=True)

    vcf_variants = load_vcf_variants(args.vcf_file, args.chr)
    node_starts = load_node_starts(args.node_pos_json)

    node_dirs = [p for p in args.tensor_folder_path.iterdir() if p.is_dir()]

    summary_records: List[dict] = []
    true_ct = false_ct = total_processed = 0

    with cf.ProcessPoolExecutor(max_workers=args.workers) as pool:
        iter_args = [
            (p.name, p, node_starts, vcf_variants, true_dir, false_dir)
            for p in node_dirs
        ]
        for recs, t_ct, f_ct, processed in pool.map(process_node, iter_args):
            summary_records.extend(recs)
            true_ct += t_ct
            false_ct += f_ct
            total_processed += processed
            if total_processed and total_processed % 1000 == 0:
                print(f"[Progress] {total_processed:,} tensors processed …")

    # write summary JSON
    summary_path = args.output_folder / "classification_summary.json"
    with summary_path.open("w") as fh:
        json.dump({"chromosome": args.chr, "results": summary_records}, fh, indent=2)
    print(f"Written summary → {summary_path}")

    print("\n--- FINAL STATS ---")
    print(f"Tensors processed : {total_processed:,}")
    print(f"  TRUE            : {true_ct:,}")
    print(f"  FALSE           : {false_ct:,}")

if __name__ == "__main__":
    if sys.version_info < (3, 8):
        sys.exit("Python 3.8+ required for ProcessPoolExecutor improvements.")
    main()