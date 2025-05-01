#!/usr/bin/env python3
import struct
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

RECORD_STRUCT = struct.Struct("<h150s150s20shc")
RECORD_SIZE = RECORD_STRUCT.size

def reverse_complement(seq):
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]

def load_node_sequences(gfa_path):
    node_sequences = {}
    with open(gfa_path, "r") as f:
        for line in f:
            if line.startswith("S"):
                fields = line.strip().split('\t')
                if len(fields) >= 3:
                    try:
                        node_id = int(fields[1])
                        node_sequences[node_id] = fields[2]
                    except ValueError:
                        continue
    return node_sequences

def load_segments_from_dat_idx(dat_path, idx_path):
    segments_by_node = defaultdict(list)
    with open(idx_path, "rb") as idx_f, open(dat_path, "rb") as dat_f:
        num_nodes = struct.unpack("<I", idx_f.read(4))[0]
        for _ in range(num_nodes):
            node_id, offset, block_size, n_records, _ = struct.unpack("<I Q I I H", idx_f.read(22))
            dat_f.seek(offset + 10)  # skip block header
            for _ in range(n_records):
                data = dat_f.read(RECORD_SIZE)
                if not data or data == b'\x00' * RECORD_SIZE:
                    continue
                offset_val, seq, bq, cigar, rq, strand = RECORD_STRUCT.unpack(data)
                segments_by_node[node_id].append({
                    'offset': offset_val,
                    'seq': seq.rstrip(b'\x00').decode('ascii', errors='ignore'),
                    'base_quality': bq.rstrip(b'\x00').decode('ascii', errors='ignore'),
                    'cigar': cigar.rstrip(b'\x00').decode('ascii', errors='ignore'),
                    'read_quality': rq,
                    'strand': strand.decode()
                })
    return segments_by_node

def normalize_node_segments(args):
    node_id, segments, node_seq = args
    norm_segments = []
    node_len = len(node_seq)
    for seg in segments:
        read_len = len(seg['seq'])
        if seg['strand'] == '-':
            new_offset = node_len - seg['offset'] - read_len
            new_seq = reverse_complement(seg['seq'])
            new_bq = seg['base_quality'][::-1]
            norm_segments.append({
                'offset': new_offset,
                'seq': new_seq,
                'base_quality': new_bq,
                'read_quality': seg['read_quality'],
                'strand': '+'
            })
        else:
            norm_segments.append(seg)
    return node_id, norm_segments

def normalize_all_segments(segments_by_node, node_sequences, workers=8):
    tasks = []
    for node_id, segs in segments_by_node.items():
        node_seq = node_sequences.get(node_id)
        if node_seq:
            tasks.append((node_id, segs, node_seq))

    norm_dict = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for node_id, norm_segs in pool.map(normalize_node_segments, tasks):
            norm_dict[node_id] = norm_segs
    return norm_dict

def main():
    parser = argparse.ArgumentParser(description="Normalize read segments to + strand from .dat/.idx and GFA")
    parser.add_argument("dat", help="Path to .dat file")
    parser.add_argument("idx", help="Path to .idx file")
    parser.add_argument("gfa", help="Path to .gfa file")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel processes")
    args = parser.parse_args()

    print("🔹 Loading segments from .dat/.idx...")
    segments_by_node = load_segments_from_dat_idx(args.dat, args.idx)

    print("🔹 Loading node sequences from GFA...")
    node_sequences = load_node_sequences(args.gfa)

    print("🔹 Normalizing segments to + strand in parallel...")
    normalized = normalize_all_segments(segments_by_node, node_sequences, workers=args.workers)

    for node_id, segs in normalized.items():
        print(f"\nNode {node_id} — {len(segs)} reads")
        for seg in segs[:3]:  # preview only 3 per node
            print(f"  Offset: {seg['offset']}  Strand: {seg['strand']}  Seq: {seg['seq'][:30]}...")

if __name__ == "__main__":
    main()
