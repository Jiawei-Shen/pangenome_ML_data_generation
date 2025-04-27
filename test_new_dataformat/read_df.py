#!/usr/bin/env python3
import struct
import argparse

# ─────────────────────────────────────────────────────────────────────────────
RECORD_STRUCT = struct.Struct("<h150s150shc")  # 结构
RECORD_SIZE = RECORD_STRUCT.size  # 305 bytes

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
            offset, seq, bq, rq, strand = RECORD_STRUCT.unpack(raw)
            record = {
                "offset": offset,
                "seq": seq.rstrip(b'\x00').decode(errors='ignore'),
                "base_quality": bq.rstrip(b'\x00').decode(errors='ignore'),
                "read_quality": rq,
                "strand": strand.decode()
            }
            records.append(record)

    return records

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Read segments for a given node_id from binary .dat/.idx files.")
    parser.add_argument("dat_path", help="Path to .dat file")
    parser.add_argument("idx_path", help="Path to .idx file")
    parser.add_argument("node_id", type=int, help="Node ID to fetch")
    parser.add_argument("--topn", type=int, default=10, help="Show top N records (default 5, use -1 to show all)")
    args = parser.parse_args()

    node_index = load_index(args.idx_path)
    records = read_segments(args.dat_path, node_index, args.node_id)

    print(f"Node {args.node_id} has {len(records)} records.")

    # 使用 --topn 控制输出
    if args.topn == -1:
        for r in records:
            print(r)
    else:
        for r in records[:args.topn]:
            print(r)

if __name__ == "__main__":
    main()
