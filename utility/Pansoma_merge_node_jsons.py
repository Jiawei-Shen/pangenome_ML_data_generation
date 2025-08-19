#!/usr/bin/env python3
"""
Build a new JSON from:
  - merged node JSON (list of node dicts OR dict with "nodes"),
  - GFA file,
  - node IDs loaded from a binary .idx file (provided loader).

Rules:
  * If node_id exists in merged JSON, copy its FULL record as-is.
  * Otherwise, read its sequence from GFA S-lines and output ONLY:
        {"node_id": "<id>", "sequence": "<seq>"}
  * Output is a JSON LIST of node records, in the order read from .idx.

Now with progress logging that reports ONLY GFA segments (nodes).
"""

import argparse
import gzip
import json
import os
import struct
import sys
import time
from typing import Any, Dict, List, Set

# ── logging helpers ───────────────────────────────────────────────────────────

def ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

def log(msg: str) -> None:
    print(f"[{ts()}] {msg}", file=sys.stderr, flush=True)

def human_int(n: int) -> str:
    return f"{n:,}"

def rate(count: int, seconds: float, unit: str) -> str:
    if seconds <= 0:
        return f"∞ {unit}/s"
    return f"{count/seconds:,.1f} {unit}/s"

# ── I/O helpers ───────────────────────────────────────────────────────────────

def open_maybe_gzip(path: str, mode: str = "rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def load_merged_map(path: str, id_key: str = "node_id") -> Dict[str, Dict[str, Any]]:
    t0 = time.time()
    log(f"Loading merged JSON: {path}")
    with open_maybe_gzip(path, "rt") as f:
        data = json.load(f)

    nodes = data.get("nodes") if isinstance(data, dict) else data
    if not isinstance(nodes, list):
        raise ValueError("Merged JSON must be a list of nodes or a dict with 'nodes'.")

    out: Dict[str, Dict[str, Any]] = {}
    for n in nodes:
        if isinstance(n, dict) and id_key in n:
            out[str(n[id_key])] = n
    dt = time.time() - t0
    log(f"Merged JSON loaded: {human_int(len(out))} nodes in {dt:.2f}s "
        f"({rate(len(out), dt, 'nodes')})")
    return out

# ── GFA scan (segments-only progress) ─────────────────────────────────────────

def scan_gfa_for_sequences(
    gfa_path: str,
    wanted_ids: Set[str],
    log_every_segments: int = 1_000_000,
    log_every_seconds: float = 5.0,
) -> Dict[str, str]:
    """
    Stream the GFA once and collect sequences for the wanted node IDs.
    Counts **only S-segments** for progress (no edges/other lines).

    GFA S-line format:  S <node_id> <sequence> ...
    """
    seqs: Dict[str, str] = {}
    if not wanted_ids:
        log("No missing IDs to fetch from GFA.")
        return seqs

    t0 = time.time()
    last_log_t = t0
    seg_scanned = 0
    found = 0
    need = len(wanted_ids)
    log(f"Scanning GFA segments for {human_int(need)} missing nodes: {gfa_path}")

    with open_maybe_gzip(gfa_path, "rt") as f:
        for line in f:
            if not line or line[0] != "S":
                # ignore edges/paths/headers entirely
                if log_every_seconds and (time.time() - last_log_t) >= log_every_seconds:
                    dt = time.time() - t0
                    log(f"GFA segments: scanned={human_int(seg_scanned)} "
                        f"found={human_int(found)}/{human_int(need)} "
                        f"| {rate(seg_scanned, dt, 'segments')}, {rate(found, dt, 'found')}")
                    last_log_t = time.time()
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            seg_scanned += 1

            nid, seq = parts[1], parts[2]
            if nid in wanted_ids and nid not in seqs:
                seqs[nid] = "" if seq == "*" else seq
                found += 1

            # segment-count based logging
            if log_every_segments and (seg_scanned % log_every_segments == 0):
                dt = time.time() - t0
                log(f"GFA segments: scanned={human_int(seg_scanned)} "
                    f"found={human_int(found)}/{human_int(need)} "
                    f"| {rate(seg_scanned, dt, 'segments')}, {rate(found, dt, 'found')}")

            # time-based logging
            if log_every_seconds and (time.time() - last_log_t) >= log_every_seconds:
                dt = time.time() - t0
                log(f"GFA segments: scanned={human_int(seg_scanned)} "
                    f"found={human_int(found)}/{human_int(need)} "
                    f"| {rate(seg_scanned, dt, 'segments')}, {rate(found, dt, 'found')}")
                last_log_t = time.time()

            # stop once we've found all requested nodes
            if found == need:
                break

    dt_total = time.time() - t0
    log(f"GFA segment scan done: scanned={human_int(seg_scanned)}, "
        f"found={human_int(found)}/{human_int(need)} in {dt_total:.2f}s "
        f"({rate(seg_scanned, dt_total, 'segments')}, {rate(found, dt_total, 'found')})")
    return seqs

# ── .idx loader (unchanged, with timing) ──────────────────────────────────────

def load_index(idx_path):
    """Loads node IDs from a binary .idx file. (Unchanged)"""
    t0 = time.time()
    log(f"Loading IDX: {idx_path}")
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}.", file=sys.stderr)
                return {}
            blocks_num, = struct.unpack("<I", blocks_num_bytes)
            for _ in range(blocks_num):
                header_data = f.read(22)
                if len(header_data) < 22:
                    break
                block_id, _, _, _, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    f.read(metadata_len)  # Skip metadata
                node_index[block_id] = True
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}
    dt = time.time() - t0
    log(f"IDX loaded: {human_int(len(node_index))} IDs in {dt:.2f}s "
        f"({rate(len(node_index), dt, 'IDs')})")
    return node_index

# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Copy full records if present; otherwise add node_id+sequence from GFA; IDs come from .idx. Progress logs show only GFA segments."
    )
    ap.add_argument("--merged", required=True, help="Merged node JSON (list or dict with 'nodes'). Supports .gz")
    ap.add_argument("--gfa", required=True, help="GFA file (.gfa or .gfa.gz)")
    ap.add_argument("--idx", required=True, help="Binary .idx file containing node IDs")
    ap.add_argument("--out", required=True, help="Output JSON (list of nodes)")
    ap.add_argument("--id-key", default="node_id", help="Key for node id in JSON (default: node_id)")
    ap.add_argument("--seq-key", default="sequence", help="Key for sequence in JSON (default: sequence)")
    # New: segments-based progress intervals
    ap.add_argument("--gfa-log-segments", type=int, default=1_000_000,
                    help="Log every N GFA segments scanned (0=disable)")
    ap.add_argument("--gfa-log-seconds", type=float, default=5.0,
                    help="Log every N seconds during GFA scan (0=disable)")
    args = ap.parse_args()

    t_all = time.time()
    log("=== START ===")

    idx_map = load_index(args.idx)
    node_ids: List[str] = [str(k) for k in idx_map.keys()]
    log(f"Node ID list prepared from IDX: {human_int(len(node_ids))} IDs")

    merged_map = load_merged_map(args.merged, id_key=args.id_key)

    t0 = time.time()
    missing: Set[str] = {nid for nid in node_ids if nid not in merged_map}
    dt = time.time() - t0
    log(f"Missing IDs to fetch from GFA: {human_int(len(missing))} (computed in {dt:.2f}s)")

    gfa_seqs = scan_gfa_for_sequences(
        args.gfa,
        missing,
        log_every_segments=args.gfa_log_segments,
        log_every_seconds=args.gfa_log_seconds,
    )

    t0 = time.time()
    out_nodes: List[Dict[str, Any]] = []
    for nid in node_ids:
        if nid in merged_map:
            out_nodes.append(merged_map[nid])
        else:
            out_nodes.append({args.id_key: nid, args.seq_key: gfa_seqs.get(nid, "")})
    dt_build = time.time() - t0
    log(f"Output list built: {human_int(len(out_nodes))} records in {dt_build:.2f}s "
        f"({rate(len(out_nodes), dt_build, 'records')})")

    t0 = time.time()
    with open(args.out, "w") as fo:
        json.dump(out_nodes, fo, ensure_ascii=False)
    dt_write = time.time() - t0
    size_str = ""
    try:
        sz = os.path.getsize(args.out)
        size_str = f", size={sz/1e6:.2f} MB"
    except Exception:
        pass
    log(f"Wrote output JSON: {args.out} in {dt_write:.2f}s{size_str}")

    copied_full = len(node_ids) - len(missing)
    added_from_gfa = sum(1 for nid in node_ids if nid in gfa_seqs)
    missing_seq = sum(1 for nid in node_ids if (nid not in merged_map and nid not in gfa_seqs))
    t_total = time.time() - t_all

    log(f"[SUMMARY] total={human_int(len(node_ids))} | "
        f"copied_full={human_int(copied_full)} | "
        f"added_from_gfa={human_int(added_from_gfa)} | "
        f"missing_seq={human_int(missing_seq)} | "
        f"elapsed={t_total:.2f}s | "
        f"overall rate={rate(len(node_ids), t_total, 'IDs processed')}")
    log("=== DONE ===")

if __name__ == "__main__":
    main()
