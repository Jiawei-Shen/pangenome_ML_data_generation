#!/usr/bin/env python3
import argparse, json, os, shutil, sys, time, random
from multiprocessing import Pool, cpu_count
from typing import Dict, Iterable, List, Set, Tuple
import pysam

# ----------------------- worker globals (set via initializer) -----------------
G_VCF = None              # pysam.VariantFile handle
G_CHR = None              # chromosome string
G_NODE_POS: Dict[str,int] = {}
G_TRUE_DIR = None
G_FALSE_DIR = None
G_USE_SYMLINKS = False

def _init_worker(vcf_path: str, chrom: str, node_pos: Dict[str,int],
                 true_dir: str, false_dir: str, use_symlinks: bool):
    """Called once per worker; avoids pickling huge objects per task."""
    global G_VCF, G_CHR, G_NODE_POS, G_TRUE_DIR, G_FALSE_DIR, G_USE_SYMLINKS
    G_VCF = pysam.VariantFile(vcf_path)
    G_CHR = chrom
    G_NODE_POS = node_pos
    G_TRUE_DIR = true_dir
    G_FALSE_DIR = false_dir
    G_USE_SYMLINKS = use_symlinks
# ------------------------------------------------------------------------------

def format_time(seconds: float) -> str:
    seconds = int(seconds)
    h, r = divmod(seconds, 3600); m, s = divmod(r, 60)
    return f"{h:02}:{m:02}:{s:02}"

def list_node_dirs(base: str) -> Iterable[str]:
    for e in os.scandir(base):
        if e.is_dir():
            yield e.path

def load_needed_node_positions(json_path: str, needed_ids: Set[str]) -> Dict[str,int]:
    """Only keep node_ids that exist in the tensor folder."""
    # Try streaming if ijson is available
    try:
        import ijson  # type: ignore
        pos: Dict[str,int] = {}
        with open(json_path, "rb") as f:
            for node in ijson.items(f, "nodes.item"):
                nid = str(node.get("node_id"))
                if nid in needed_ids:
                    p = node.get("grch38_position_start")
                    if isinstance(p, int):
                        pos[nid] = p
        return pos
    except Exception:
        # Fallback: load once, then filter (uses more RAM)
        with open(json_path) as f:
            data = json.load(f)
        pos = {}
        for node in data.get("nodes", []):
            nid = str(node.get("node_id"))
            if nid in needed_ids:
                p = node.get("grch38_position_start")
                if isinstance(p, int):
                    pos[nid] = p
        return pos

def _vcf_subset_for_positions(positions: Set[int]) -> Set[Tuple[int,str,str]]:
    """Fetch a small window from the VCF and build a tiny membership set."""
    if not positions:
        return set()
    mn, mx = min(positions), max(positions)
    out: Set[Tuple[int,str,str]] = set()
    # pysam fetch uses 0-based start, half-open end; rec.pos is 1-based
    for rec in G_VCF.fetch(G_CHR, max(0, mn-1), mx+1):
        if rec.pos in positions and rec.alts:
            r = rec.ref.upper()
            for a in rec.alts:
                out.add((rec.pos, r, a.upper()))
    return out

def _classify_node(node_dir: str):
    """Worker task: classify one node directory; returns (records, true, false)."""
    node_id = os.path.basename(node_dir)
    start_pos = G_NODE_POS.get(node_id)
    if start_pos is None:
        return [], 0, 0

    summary_path = os.path.join(node_dir, "variant_summary.json")
    try:
        with open(summary_path) as f:
            summary = json.load(f)
    except Exception:
        return [], 0, 0

    variants = summary.get("variants_passing_af_filter", [])
    if not variants:
        return [], 0, 0

    # First pass: collect needed genomic positions & parse items we’ll output
    items = []  # (tensor_path, tensor_file, pos, ref, alt)
    pos_needed: Set[int] = set()
    for v in variants:
        tf = v.get("tensor_file"); vk = v.get("variant_key")
        if not tf or not vk or not tf.endswith(".npy"):
            continue
        tpath = os.path.join(node_dir, tf)
        if not os.path.isfile(tpath):
            continue
        try:
            off, ref, alt = _parse_key(vk)
        except Exception:
            continue
        pos = start_pos + off
        items.append((tpath, tf, pos, ref, alt))
        pos_needed.add(pos)

    if not items:
        return [], 0, 0

    # Build a small, node-local VCF membership set
    local_vcf = _vcf_subset_for_positions(pos_needed)

    recs = []
    t = f = 0
    for tpath, tf, pos, ref, alt in items:
        is_match = (pos, ref, alt) in local_vcf
        dest = G_TRUE_DIR if is_match else G_FALSE_DIR
        if is_match: t += 1
        else:        f += 1

        dst = os.path.join(dest, f"{node_id}_{tf}")
        try:
            if G_USE_SYMLINKS:
                if not os.path.lexists(dst):
                    os.symlink(os.path.abspath(tpath), dst)
            else:
                if not os.path.exists(dst):
                    shutil.copyfile(tpath, dst)
        except OSError:
            pass

        recs.append({
            "node_id": node_id, "tensor_file": tf, "genomic_position": pos,
            "ref": ref, "alt": alt, "classification": "true" if is_match else "false"
        })

    return recs, t, f

def _parse_key(variant_key: str) -> Tuple[int,str,str]:
    # variant_key looks like "OFFSET_x_REF_ALT_..." — adapt if needed
    parts = variant_key.split("_")
    offset = int(parts[0])
    ref = parts[2].upper()
    alt = parts[3].upper()
    return offset, ref, alt

def _copy_or_link(src: str, dst: str, use_symlinks: bool) -> None:
    """Create dst as a symlink to src, or copy if symlinks not requested/possible."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if use_symlinks:
        if os.path.lexists(dst):
            return
        try:
            os.symlink(os.path.abspath(src), dst)
            return
        except OSError:
            # fall back to copy if symlink fails on this FS
            pass
    if not os.path.exists(dst):
        shutil.copy2(src, dst)

def organize_classified_data(base_dir: str, ratios, use_symlinks: bool, seed=None):
    """
    Populate base_dir/{train,val,test}/{true,false} by SYMLINKING (or copying) from
    base_dir/{true,false}. Do NOT move or delete the source files.
    """
    if seed is not None:
        random.seed(seed)
    assert len(ratios) == 3 and abs(sum(ratios) - 1.0) < 1e-6

    for label in ("true", "false"):
        src = os.path.join(base_dir, label)
        if not os.path.isdir(src):
            continue

        files = sorted(os.listdir(src))
        random.shuffle(files)

        n = len(files)
        n_tr = int(n * ratios[0])
        n_va = int(n * ratios[1])
        n_te = n - n_tr - n_va

        splits = {
            "train": files[:n_tr],
            "val":   files[n_tr:n_tr + n_va],
            "test":  files[n_tr + n_va:]
        }

        # Create links/copies into split dirs; keep originals intact.
        for sp, lst in splits.items():
            outd = os.path.join(base_dir, sp, label)
            os.makedirs(outd, exist_ok=True)
            for f in lst:
                src_path = os.path.join(src, f)
                dst_path = os.path.join(outd, f)
                _copy_or_link(src_path, dst_path, use_symlinks)

    # IMPORTANT: no deletion of source dirs, no moves.

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tensor_folder_path")
    ap.add_argument("vcf_file")
    ap.add_argument("node_pos_json")
    ap.add_argument("--output_folder", default="./classification_results")
    ap.add_argument("--chr", default="chr1")
    ap.add_argument("--use-symlinks", action="store_true",
                    help="Use symlinks instead of copying for both classification and organize steps.")
    ap.add_argument("-j","--workers", type=int, default=max(1, min(8, cpu_count())))
    ap.add_argument("--chunksize", type=int, default=32)
    ap.add_argument("--maxtasksperchild", type=int, default=500)
    ap.add_argument("--organize", nargs=3, type=float, metavar=("TRAIN","VAL","TEST"),
                    help="Ratios for train/val/test (e.g., --organize 0.8 0.1 0.1).")
    ap.add_argument("--seed", type=int)
    args = ap.parse_args()

    os.makedirs(args.output_folder, exist_ok=True)
    true_dir  = os.path.join(args.output_folder, "true")
    false_dir = os.path.join(args.output_folder, "false")
    os.makedirs(true_dir, exist_ok=True); os.makedirs(false_dir, exist_ok=True)

    # Node dirs and the set of IDs we actually need
    node_dirs = list(list_node_dirs(args.tensor_folder_path))
    need_ids = {os.path.basename(p) for p in node_dirs}
    print(f"Found {len(node_dirs):,} node dirs.")

    # Load ONLY the positions we need (streaming if possible)
    node_pos = load_needed_node_positions(args.node_pos_json, need_ids)
    print(f"Loaded positions for {len(node_pos):,} nodes from JSON.")

    # Stream summary to NDJSON to avoid big in-RAM arrays
    summary_path = os.path.join(args.output_folder, "classification_summary.ndjson")
    summary_f = open(summary_path, "w")

    total_nodes = len(node_dirs); total_true = total_false = 0
    start = time.monotonic()

    with Pool(processes=args.workers,
              initializer=_init_worker,
              initargs=(args.vcf_file, args.chr, node_pos,
                        true_dir, false_dir, args.use_symlinks),
              maxtasksperchild=args.maxtasksperchild) as pool:

        for i, (records, t_cnt, f_cnt) in enumerate(
            pool.imap_unordered(_classify_node, node_dirs, chunksize=args.chunksize)
        ):
            if records:
                for rec in records:
                    summary_f.write(json.dumps(rec) + "\n")

            total_true  += t_cnt
            total_false += f_cnt

            done = i + 1
            elapsed = time.monotonic() - start
            rate = done / elapsed if elapsed > 0 else 0.0
            eta = format_time((total_nodes - done)/rate) if rate > 0 else "..."
            print(f"\rProgress {done}/{total_nodes} "
                  f"({done/total_nodes:.1%})  true={total_true} false={total_false}  "
                  f"{rate:.2f} nodes/s  ETA {eta}    ", end="")

    summary_f.close()
    print("\nDone.")
    print(f"Summary written to: {summary_path}")
    print(f"True: {total_true}  False: {total_false}")

    if args.organize:
        print("Organizing into train/val/test (non-destructive)…")
        organize_classified_data(args.output_folder, args.organize, use_symlinks=args.use_symlinks, seed=args.seed)
        print("Organization complete (sources kept in place).")

if __name__ == "__main__":
    main()
