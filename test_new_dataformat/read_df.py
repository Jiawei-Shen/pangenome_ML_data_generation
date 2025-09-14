#!/usr/bin/env python3
import struct
import argparse
import sys
import os

# ─────────────────────────────────────────────────────────────────────────────
# New-format .dat block header
# Try padded 16B first: <I I H 2x I>  (node_id, n_records, flags, pad, node_length)
BLOCK_HDR_PACK_PADDED = struct.Struct("<I I H 2x I")  # 16 bytes
BLOCK_HDR_PACK_14B    = struct.Struct("<I I H I")     # 14 bytes (no pad)

def make_record_struct(node_length: int) -> struct.Struct:
    """
    Per-read record in the new .dat format (variable-sized by node_length):
      <h {L}s {L}s {L}s h c>
        - i16 offset
        - seq[L] bytes (null-padded ASCII)
        - bq[L] bytes (raw qualities, 0..~60; null-padded)
        - cigar[L] bytes (ASCII, null-padded)
        - i16 rq (MAPQ)
        - char strand ('+' / '-')
    """
    return struct.Struct(f"<h{node_length}s{node_length}s{node_length}shc")

# ─────────────────────────────────────────────────────────────────────────────
def load_index(idx_path):
    """
    Supports:
      • New fixed-size entry: 26B  <I Q I I H I>
      • Older fixed-size entry: 22B <I Q I I H>
      • Legacy variable-size with metadata_len per entry: <I Q I I H> + metadata_len bytes
    Returns: { node_id: {"start": offset, "size": block_size, "n_records": n_records, "flags": flags, "node_length": node_len_or_None} }
    """
    node_index = {}

    with open(idx_path, "rb") as f:
        # First 4 bytes: number of entries
        header = f.read(4)
        if len(header) < 4:
            raise ValueError(f"Index file too small: {idx_path}")
        blocks_num, = struct.unpack("<I", header)

        # Detect layout
        file_size = os.fstat(f.fileno()).st_size
        remaining = file_size - 4

        def read_fixed_26():
            rec = f.read(26)
            if len(rec) != 26:
                raise EOFError("Unexpected EOF in 26B index entry")
            node_id, start, size, nrec, flags, node_len = struct.unpack("<I Q I I H I", rec)
            return node_id, start, size, nrec, flags, node_len

        def read_fixed_22():
            rec = f.read(22)
            if len(rec) != 22:
                raise EOFError("Unexpected EOF in 22B index entry")
            node_id, start, size, nrec, flags = struct.unpack("<I Q I I H", rec)
            return node_id, start, size, nrec, flags, None

        def read_legacy_with_meta():
            hdr = f.read(4+8+4+4+2)
            if len(hdr) != 22:
                raise EOFError("Unexpected EOF in legacy index header")
            node_id, start, size, nrec, meta_len = struct.unpack("<I Q I I H", hdr)
            if meta_len:
                skipped = f.read(meta_len)
                if len(skipped) != meta_len:
                    raise EOFError("Unexpected EOF while skipping metadata")
            return node_id, start, size, nrec, 0, None

        # Choose strategy
        strategy = None
        if blocks_num > 0 and remaining // blocks_num == 26 and remaining % blocks_num == 0:
            strategy = "fixed26"
        elif blocks_num > 0 and remaining // blocks_num == 22 and remaining % blocks_num == 0:
            strategy = "fixed22"
        else:
            strategy = "legacy"

        # Read entries
        for _ in range(blocks_num):
            if strategy == "fixed26":
                node_id, start, size, nrec, flags, node_len = read_fixed_26()
            elif strategy == "fixed22":
                node_id, start, size, nrec, flags, node_len = read_fixed_22()
            else:
                node_id, start, size, nrec, flags, node_len = read_legacy_with_meta()

            node_index[node_id] = {
                "start": start,
                "size": size,
                "n_records": nrec,
                "flags": flags,
                "node_length": node_len,   # may be None for older formats
            }

    return node_index

# ─────────────────────────────────────────────────────────────────────────────
def _read_block_header(dat_file, block_start):
    """
    Read .dat block header at block_start.
    Tries 16B padded format first; falls back to 14B if needed.
    Returns: (node_id, n_records, flags, node_length, header_size_bytes)
    """
    # Try padded 16B
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_PADDED.size)
    if len(hdr) == BLOCK_HDR_PACK_PADDED.size:
        nid, nrec, flags, node_len = BLOCK_HDR_PACK_PADDED.unpack(hdr)
        if 1 <= node_len <= 1_000_000 and 0 <= nrec < 10_000_000:
            return nid, nrec, flags, node_len, BLOCK_HDR_PACK_PADDED.size

    # Fallback: 14B
    dat_file.seek(block_start, os.SEEK_SET)
    hdr = dat_file.read(BLOCK_HDR_PACK_14B.size)
    if len(hdr) != BLOCK_HDR_PACK_14B.size:
        raise RuntimeError(f"Cannot read block header at offset {block_start}")
    nid, nrec, flags, node_len = BLOCK_HDR_PACK_14B.unpack(hdr)
    if not (1 <= node_len <= 1_000_000 and 0 <= nrec < 10_000_000):
        raise RuntimeError(f"Suspicious block header at {block_start}: node_len={node_len}, nrec={nrec}")
    return nid, nrec, flags, node_len, BLOCK_HDR_PACK_14B.size

# ─────────────────────────────────────────────────────────────────────────────
def read_segments(dat_path, node_index, node_id):
    if node_id not in node_index:
        raise ValueError(f"Node ID {node_id} not found in index")

    info = node_index[node_id]
    block_start = info["start"]

    records = []

    with open(dat_path, "rb") as f:
        # Read header from .dat to get node_length and true n_records
        nid, n_records, flags, node_length, hdr_size = _read_block_header(f, block_start)
        if nid != node_id:
            # tolerate mismatch but continue; users sometimes re-pack idx/dat differently
            pass

        rec_struct = make_record_struct(int(node_length))
        rec_size = rec_struct.size

        # Jump to records start
        block_content_start = block_start + hdr_size
        f.seek(block_content_start, os.SEEK_SET)

        for _ in range(n_records):
            raw = f.read(rec_size)
            if len(raw) < rec_size:
                break
            offset, seq, bq, cigar_bytes, rq, strand = rec_struct.unpack(raw)

            record = {
                "offset": int(offset),
                "sequence": seq.rstrip(b'\x00').decode('ascii', errors='ignore'),
                # base qualities are raw bytes; keep as list of ints for correctness
                "base_quality": list(bq.rstrip(b'\x00')),
                "cigar": cigar_bytes.rstrip(b'\x00').decode('ascii', errors='ignore'),
                "mapping_quality": int(rq),
                "strand": strand.decode('ascii', errors='ignore') if isinstance(strand, (bytes, bytearray)) else chr(strand),
            }
            records.append(record)

    return records

# ─────────────────────────────────────────────────────────────────────────────
def print_node_ids(node_index, n):
    node_ids = list(node_index.keys())
    node_ids.sort()
    print(f"Total nodes: {len(node_ids)}")
    print("Listing first {} node IDs:".format(n if n != -1 else len(node_ids)))
    for nid in node_ids[:n if n != -1 else None]:
        print(nid)

# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Read segments or print node IDs from binary .dat/.idx files (supports new 26B .idx and legacy formats)."
    )
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
