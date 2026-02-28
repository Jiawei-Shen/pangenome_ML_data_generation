#!/usr/bin/env python3
import argparse
import pickle
import sys
from typing import Any, Dict, Tuple, Optional


def load_pickle(pkl_path: str) -> Dict[int, Dict[str, int]]:
    with open(pkl_path, "rb") as f:
        obj = pickle.load(f)

    if not isinstance(obj, dict):
        raise TypeError(f"Pickle root object is {type(obj)}, expected dict")

    # Your writer stores int keys. But in case someone converted keys to str, normalize.
    out: Dict[int, Dict[str, int]] = {}
    for k, v in obj.items():
        try:
            ik = int(k)
        except Exception:
            continue
        if isinstance(v, dict):
            out[ik] = v
    return out


def compute_stats(rec: Dict[str, int]) -> Tuple[int, int, int, float, float]:
    perfect = int(rec.get("perfect", 0))
    not_perfect = int(rec.get("not_perfect", 0))
    total = perfect + not_perfect
    not_perfect_ratio = (not_perfect / total) if total > 0 else 0.0
    perfect_ratio = (perfect / total) if total > 0 else 0.0
    return perfect, not_perfect, total, perfect_ratio, not_perfect_ratio


def print_node(node_id: int, rec: Optional[Dict[str, int]]) -> None:
    if rec is None:
        print(f"[NOT FOUND] node_id={node_id}")
        return

    perfect, not_perfect, total, perfect_ratio, not_perfect_ratio = compute_stats(rec)

    max_read_length = int(rec.get("max_read_length", 0))
    max_cigar_length = int(rec.get("max_cigar_length", 0))

    print(f"node_id: {node_id}")
    print(f"  perfect:        {perfect}")
    print(f"  not_perfect:    {not_perfect}")
    print(f"  total:          {total}")
    print(f"  perfect_ratio:  {perfect_ratio:.6f}")
    print(f"  notperf_ratio:  {not_perfect_ratio:.6f}")
    print(f"  max_read_length (sum to_length on node): {max_read_length}")
    print(f"  max_cigar_length (chars):                {max_cigar_length}")


def top_k(
    data: Dict[int, Dict[str, int]],
    by: str,
    k: int,
    min_total: int,
) -> None:
    rows = []
    for nid, rec in data.items():
        perfect, not_perfect, total, perfect_ratio, not_perfect_ratio = compute_stats(rec)
        if total < min_total:
            continue

        if by == "not_perfect":
            key = not_perfect
        elif by == "perfect":
            key = perfect
        elif by == "total":
            key = total
        elif by == "notperf_ratio":
            key = not_perfect_ratio
        elif by == "perfect_ratio":
            key = perfect_ratio
        else:
            raise ValueError(f"Unknown sort key: {by}")

        rows.append((key, nid, perfect, not_perfect, total, perfect_ratio, not_perfect_ratio, rec))

    rows.sort(key=lambda x: x[0], reverse=True)
    rows = rows[:k]

    header = (
        "rank\tkey\tnode_id\tperfect\tnot_perfect\ttotal\tperfect_ratio\tnotperf_ratio\t"
        "max_read_length\tmax_cigar_length"
    )
    print(header)
    for i, (key, nid, p, np, total, pr, npr, rec) in enumerate(rows, 1):
        mrl = int(rec.get("max_read_length", 0))
        mcl = int(rec.get("max_cigar_length", 0))
        if isinstance(key, float):
            key_str = f"{key:.6f}"
        else:
            key_str = str(key)
        print(
            f"{i}\t{key_str}\t{nid}\t{p}\t{np}\t{total}\t{pr:.6f}\t{npr:.6f}\t{mrl}\t{mcl}"
        )


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Query reads_by_node.pkl produced by your GAM parser (node_id -> counts/maxes)."
    )
    ap.add_argument("pkl", help="Path to reads_by_node.pkl")
    ap.add_argument("--node", type=int, help="Node ID to query")
    ap.add_argument(
        "--top",
        type=int,
        default=0,
        help="If >0, list top K nodes instead of a single node (or in addition, if --node is also set).",
    )
    ap.add_argument(
        "--by",
        choices=["not_perfect", "perfect", "total", "notperf_ratio", "perfect_ratio"],
        default="not_perfect",
        help="Metric to sort by for --top (default: not_perfect).",
    )
    ap.add_argument(
        "--min_total",
        type=int,
        default=1,
        help="Minimum total segments (perfect+not_perfect) to include in --top (default: 1).",
    )
    args = ap.parse_args()

    data = load_pickle(args.pkl)
    print(f"[INFO] Loaded {len(data)} nodes from {args.pkl}")

    if args.node is None and args.top <= 0:
        print("[ERROR] Provide --node <ID> and/or --top <K>.")
        sys.exit(2)

    if args.node is not None:
        print_node(args.node, data.get(args.node))

    if args.top and args.top > 0:
        print()
        top_k(data, by=args.by, k=args.top, min_total=args.min_total)


if __name__ == "__main__":
    main()