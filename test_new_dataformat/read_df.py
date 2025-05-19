#!/usr/bin/env python3
import struct
import argparse
import sys

# ─────────────────────────────────────────────────────────────────────────────
RECORD_STRUCT = struct.Struct("<h150s150s20shc")  # updated structure
RECORD_SIZE = RECORD_STRUCT.size  # 325 bytes

# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path):
    node_index = {}
    with open(idx_path, "rb") as f:
        blocks_num, = struct.unpack("<I", f.read(4))
        for _ in range(blocks_num):
            block_id, block_start, block_size, n_records, metadata_len = struct.unpack("<I Q I I H", f.read(4+8+4+4+2))
            if metadata_len > 0:
                f.read(metadata_len)  # skip metadata
            node_index[block_id] = {
                "start": block_start,
                "size": block_size,
                "n_records": n_records
            }
    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def read_segments(dat_path, node_index, node_id):
    if node_id not in node_index:
        raise ValueError(f"Node ID {node_id} not found in index")

    info = node_index[node_id]
    block_start = info["start"]
    n_records = info["n_records"]

    records = []

    with open(dat_path, "rb") as f:
        block_content_start = block_start + 4 + 4 + 2  # skip block header
        f.seek(block_content_start)

        for _ in range(n_records):
            raw = f.read(RECORD_SIZE)
            if not raw:
                break
            offset, seq, bq, cigar_bytes, rq, strand = RECORD_STRUCT.unpack(raw)

            record = {
                "offset": offset,
                "sequence": seq.rstrip(b'\x00').decode(errors='ignore'),
                "base_quality": bq.rstrip(b'\x00').decode(errors='ignore'),
                "cigar": cigar_bytes.rstrip(b'\x00').decode(errors='ignore'),
                "mapping_quality": rq,
                "strand": strand.decode()
            }
            records.append(record)

    return records

# ─────────────────────────────────────────────────────────────────────────────
def print_node_ids(node_index, n):
    node_ids = list(node_index.keys())
    print(f"Total nodes: {len(node_ids)}")
    print("Listing first {} node IDs:".format(n if n != -1 else len(node_ids)))
    for nid in node_ids[:n if n != -1 else None]:
        print(nid)

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Read segments or print node IDs from binary .dat/.idx files (new format with cigar).")
    parser.add_argument("dat_path", help="Path to .dat file")
    parser.add_argument("idx_path", help="Path to .idx file")
    parser.add_argument("--node_id", type=int, help="Node ID to fetch")
    parser.add_argument("--topn", type=int, default=20, help="Show top N records (default 20, use -1 to show all)")
    parser.add_argument("--print-nodes", type=int, help="Print N node IDs from the index and exit")
    args = parser.parse_args()

    # Load index regardless of operation
    node_index = load_index(args.idx_path)

    # If printing node IDs only
    if args.print_nodes is not None:
        print_node_ids(node_index, args.print_nodes)
        return

    # node_id is required if not using --print-nodes
    if args.node_id is None:
        print("Error: --node_id is required unless --print-nodes is specified.", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    # Proceed to read segments
    records = read_segments(args.dat_path, node_index, args.node_id)
    print(f"Node {args.node_id} has {len(records)} records.")

    if args.topn == -1:
        for r in records:
            print(r)
    else:
        for r in records[:args.topn]:
            print(r)

if __name__ == "__main__":
    main()