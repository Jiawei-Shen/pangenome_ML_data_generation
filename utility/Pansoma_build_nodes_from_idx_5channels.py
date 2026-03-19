#!/usr/bin/env python3
"""
Build a JSON list of nodes from:
  - GFA file
  - one or more binary .idx files

Behavior:
  * Read node IDs from all --idx files
  * Take the UNION of node IDs, preserving first-seen order
  * Scan the GFA and collect sequence for every wanted node ID
  * Output ONLY:
        {"node_id": "<id>", "sequence": "<seq>"}
    and, if requested:
        {"node_id": "<id>", "sequence": "<seq>", "chrom": "<chr>"}

Notes:
  * No merged JSON input is needed.
  * No extra fields such as "genomead_af" are ever included.
  * If a node ID from idx is not found in the GFA, its sequence is set to "".

Usage:
  python build_nodes_from_idx_gfa_only.py \
      --gfa graph.gfa[.gz] \
      --idx nodes1.idx nodes2.idx \
      --out new_nodes.json \
      --add-chrom chr1
"""

import argparse
import gzip
import json
import os
import struct
import sys
from typing import Dict, List, Set


# ─────────────────────────────────────────────────────────────────────────────
# I/O helpers

def open_maybe_gzip(path: str, mode: str = "rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


# ─────────────────────────────────────────────────────────────────────────────
# IDX loader
#
# Latest .idx entry (30B): <I Q I I H I I>
#   nid, offset, block_size, n_records, flags, R, C
# We only need nid, but we parse strictly and preserve file order.
# ─────────────────────────────────────────────────────────────────────────────

IDX_ENTRY_PACK_LATEST = struct.Struct("<I Q I I H I I")
IDX_ENTRY_SIZE_LATEST = IDX_ENTRY_PACK_LATEST.size  # 30


def load_index_ids(idx_path: str) -> List[str]:
    """Load node IDs from a binary .idx file, preserving file order."""
    ids: List[str] = []
    try:
        with open(idx_path, "rb") as f:
            hdr = f.read(4)
            if len(hdr) != 4:
                print(f"Error: Could not read blocks_num from {idx_path}.", file=sys.stderr)
                return []

            (blocks_num,) = struct.unpack("<I", hdr)

            # Strict size check
            try:
                f.seek(0, os.SEEK_END)
                size = f.tell()
                expected = 4 + blocks_num * IDX_ENTRY_SIZE_LATEST
                if size != expected:
                    print(
                        f"Error: IDX size mismatch for latest format: file={size}, expected={expected} "
                        f"(count={blocks_num}, entry={IDX_ENTRY_SIZE_LATEST}) in {idx_path}.",
                        file=sys.stderr,
                    )
                    return []
            except Exception:
                pass

            f.seek(4, os.SEEK_SET)
            for i in range(blocks_num):
                rec = f.read(IDX_ENTRY_SIZE_LATEST)
                if len(rec) != IDX_ENTRY_SIZE_LATEST:
                    print(f"Error: Truncated IDX at entry {i+1} in {idx_path}.", file=sys.stderr)
                    return []
                block_id, _off, _blk_sz, _nrec, _flags, _R, _C = IDX_ENTRY_PACK_LATEST.unpack(rec)
                ids.append(str(block_id))

    except FileNotFoundError:
        print(f"Error: IDX file not found: {idx_path}", file=sys.stderr)
        return []
    except Exception as e:
        print(f"Unexpected error while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return []

    return ids


def union_node_ids_from_multiple_idx(idx_paths: List[str]) -> List[str]:
    """
    Union node IDs across multiple idx files, preserving first-seen order:
      - iterate idx_paths in given order
      - within each idx, preserve file order
      - keep only first occurrence of each node ID
    """
    seen: Set[str] = set()
    out: List[str] = []

    for idx_path in idx_paths:
        ids = load_index_ids(idx_path)
        for nid in ids:
            if nid not in seen:
                seen.add(nid)
                out.append(nid)

    return out


# ─────────────────────────────────────────────────────────────────────────────
# GFA scan

def scan_gfa_for_sequences(gfa_path: str, wanted_ids: Set[str]) -> Dict[str, str]:
    """
    Stream the GFA once and collect sequences for wanted node IDs.
    GFA S-line format: S <node_id> <sequence> ...
    """
    seqs: Dict[str, str] = {}
    if not wanted_ids:
        return seqs

    with open_maybe_gzip(gfa_path, "rt") as f:
        for line in f:
            if not line or line[0] != "S":
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue

            nid, seq = parts[1], parts[2]
            if nid in wanted_ids and nid not in seqs:
                seqs[nid] = "" if seq == "*" else seq
                if len(seqs) == len(wanted_ids):
                    break

    return seqs


# ─────────────────────────────────────────────────────────────────────────────
# Main

def main():
    ap = argparse.ArgumentParser(
        description="Build JSON node records from GFA using node IDs loaded from one or more .idx files."
    )
    ap.add_argument("--gfa", required=True, help="Input GFA file (.gfa or .gfa.gz)")
    ap.add_argument(
        "--idx",
        required=True,
        nargs="+",
        help="One or more binary .idx files; node IDs are unioned in first-seen order",
    )
    ap.add_argument("--out", required=True, help="Output JSON file")
    ap.add_argument(
        "--add-chrom",
        default=None,
        help='If set (e.g. "chr1"), add {"chrom": "..."} to every output record',
    )
    ap.add_argument(
        "--id-key",
        default="node_id",
        help='Output node ID key name (default: "node_id")',
    )
    ap.add_argument(
        "--seq-key",
        default="sequence",
        help='Output sequence key name (default: "sequence")',
    )
    ap.add_argument(
        "--chrom-key",
        default="chrom",
        help='Output chrom key name when --add-chrom is used (default: "chrom")',
    )

    args = ap.parse_args()

    # 1) Load node IDs from idx files
    node_ids = union_node_ids_from_multiple_idx(args.idx)
    if not node_ids:
        print("ERROR: No node IDs loaded from --idx inputs.", file=sys.stderr)
        sys.exit(2)

    # 2) Scan GFA once for all wanted node IDs
    wanted_ids = set(node_ids)
    seqs = scan_gfa_for_sequences(args.gfa, wanted_ids)

    # 3) Build output records in union order
    out_nodes = []
    for nid in node_ids:
        rec = {
            args.id_key: nid,
            args.seq_key: seqs.get(nid, ""),
        }
        if args.add_chrom is not None:
            rec[args.chrom_key] = args.add_chrom
        out_nodes.append(rec)

    # 4) Write JSON
    with open(args.out, "w", encoding="utf-8") as fo:
        json.dump(out_nodes, fo, ensure_ascii=False)

    found_seq = sum(1 for nid in node_ids if nid in seqs)
    missing_seq = len(node_ids) - found_seq
    print(
        f"[summary] total:{len(node_ids)} found_seq:{found_seq} missing_seq:{missing_seq}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()