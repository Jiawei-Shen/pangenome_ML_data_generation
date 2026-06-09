#!/usr/bin/env python3

import argparse
import logging
import struct
import time
from collections import defaultdict

GLOBAL_MAGIC = b"MYFMT\x01"
GLOBAL_VER_PACK = struct.Struct("<BBI16s")
GLOBAL_HEADER_SIZE = len(GLOBAL_MAGIC) + GLOBAL_VER_PACK.size

BLOCK_HDR_PACK = struct.Struct("<I I H I I")
BLOCK_HDR_SIZE = BLOCK_HDR_PACK.size

IDX_HEADER_PACK = struct.Struct("<I")
IDX_REC_PACK = struct.Struct("<I Q I I H I I")


def make_record_struct(max_read_len, max_cigar_len):
    return struct.Struct(f"<h{max_read_len}s{max_read_len}s{max_cigar_len}shc")


def setup_logger(log_file=None):
    logger = logging.getLogger("graph_candidate_counter")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger


def read_idx(idx_path):
    records = []

    with open(idx_path, "rb") as f:
        n_blocks = IDX_HEADER_PACK.unpack(f.read(IDX_HEADER_PACK.size))[0]

        for _ in range(n_blocks):
            raw = f.read(IDX_REC_PACK.size)
            if len(raw) != IDX_REC_PACK.size:
                raise RuntimeError("IDX file ended unexpectedly")

            nid, offset, block_size, n_records, flags, max_read_len, max_cigar_len = IDX_REC_PACK.unpack(raw)

            records.append({
                "nid": nid,
                "offset": offset,
                "block_size": block_size,
                "n_records": n_records,
                "flags": flags,
                "max_read_len": max_read_len,
                "max_cigar_len": max_cigar_len,
            })

    return records


def parse_cigar(cigar):
    ops = []
    num = ""

    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if num:
                ops.append((int(num), ch))
                num = ""

    return ops


def trim_null_bytes(x):
    return x.split(b"\x00", 1)[0]


def count_node_candidates(records, min_alt_reads, min_af, min_depth):
    """
    Approximate graph candidate logic:
    - For each node offset, count perfect-match bases as REF-like.
    - Count X/I/D CIGAR events as ALT-like alleles.
    - Candidate example = one allele/event passing alt_count and AF threshold.
    """

    depth_by_pos = defaultdict(int)
    alt_by_pos = defaultdict(lambda: defaultdict(int))

    for rec in records:
        node_offset, seq_b, bq_b, cigar_b, mapq, strand_b = rec

        seq = trim_null_bytes(seq_b).decode(errors="ignore")
        cigar = trim_null_bytes(cigar_b).decode(errors="ignore")

        if not cigar:
            continue

        ref_pos = int(node_offset)
        read_pos = 0

        for length, op in parse_cigar(cigar):
            if op == "M":
                for i in range(length):
                    depth_by_pos[ref_pos + i] += 1
                ref_pos += length
                read_pos += length

            elif op == "X":
                alt_seq = seq[read_pos: read_pos + length]
                for i, base in enumerate(alt_seq):
                    pos = ref_pos + i
                    depth_by_pos[pos] += 1
                    alt_by_pos[pos][f"X:{base}"] += 1
                ref_pos += length
                read_pos += length

            elif op == "I":
                anchor = max(ref_pos - 1, 0)
                ins_seq = seq[read_pos: read_pos + length]
                depth_by_pos[anchor] += 1
                alt_by_pos[anchor][f"I:{ins_seq}"] += 1
                read_pos += length

            elif op == "D":
                anchor = ref_pos
                for i in range(length):
                    depth_by_pos[ref_pos + i] += 1
                alt_by_pos[anchor][f"D:{length}"] += 1
                ref_pos += length

    candidate_sites = 0
    candidate_examples = 0
    snv_examples = 0
    indel_examples = 0

    for pos, allele_counts in alt_by_pos.items():
        depth = depth_by_pos.get(pos, 0)
        if depth < min_depth:
            continue

        site_has_candidate = False

        for allele, alt_count in allele_counts.items():
            vaf = alt_count / depth if depth > 0 else 0.0

            if alt_count >= min_alt_reads and vaf >= min_af:
                candidate_examples += 1
                site_has_candidate = True

                if allele.startswith("X:"):
                    snv_examples += 1
                else:
                    indel_examples += 1

        if site_has_candidate:
            candidate_sites += 1

    return candidate_sites, candidate_examples, snv_examples, indel_examples


def run(args, logger):
    t0 = time.time()

    logger.info("Starting graph candidate counting")
    logger.info(f"DAT: {args.dat}")
    logger.info(f"IDX: {args.idx}")
    logger.info(
        f"Parameters: min_depth={args.min_depth}, "
        f"min_alt_reads={args.min_alt_reads}, min_af={args.min_af}"
    )

    idx_records = read_idx(args.idx)
    logger.info(f"Loaded {len(idx_records):,} node blocks from IDX")

    total_nodes = 0
    processed_records = 0
    total_candidate_sites = 0
    total_candidate_examples = 0
    total_snv_examples = 0
    total_indel_examples = 0

    with open(args.dat, "rb") as dat:
        magic = dat.read(len(GLOBAL_MAGIC))
        if magic != GLOBAL_MAGIC:
            raise RuntimeError(f"Unexpected DAT magic: {magic}")

        for i, info in enumerate(idx_records, 1):
            nid = info["nid"]
            offset = info["offset"]
            n_records = info["n_records"]
            max_read_len = info["max_read_len"]
            max_cigar_len = info["max_cigar_len"]

            dat.seek(offset)

            block_hdr = dat.read(BLOCK_HDR_SIZE)
            if len(block_hdr) != BLOCK_HDR_SIZE:
                raise RuntimeError(f"Cannot read block header for node {nid}")

            hdr_nid, hdr_nrec, flags, hdr_R, hdr_C = BLOCK_HDR_PACK.unpack(block_hdr)

            if hdr_nid != nid:
                raise RuntimeError(f"Node ID mismatch: idx={nid}, dat={hdr_nid}")

            rec_struct = make_record_struct(max_read_len, max_cigar_len)
            node_records = []

            for _ in range(n_records):
                raw = dat.read(rec_struct.size)
                if len(raw) != rec_struct.size:
                    raise RuntimeError(f"Cannot read record for node {nid}")
                node_records.append(rec_struct.unpack(raw))

            sites, examples, snv, indel = count_node_candidates(
                node_records,
                min_alt_reads=args.min_alt_reads,
                min_af=args.min_af,
                min_depth=args.min_depth,
            )

            total_nodes += 1
            processed_records += n_records
            total_candidate_sites += sites
            total_candidate_examples += examples
            total_snv_examples += snv
            total_indel_examples += indel

            if total_nodes % args.log_every_nodes == 0:
                elapsed = time.time() - t0
                logger.info(
                    f"Progress: nodes={total_nodes:,}/{len(idx_records):,}, "
                    f"records={processed_records:,}, "
                    f"candidate_sites={total_candidate_sites:,}, "
                    f"candidate_examples={total_candidate_examples:,}, "
                    f"elapsed={elapsed/60:.2f} min"
                )

    elapsed = time.time() - t0

    logger.info("Finished graph candidate counting")
    logger.info(f"nodes_processed:      {total_nodes:,}")
    logger.info(f"records_processed:    {processed_records:,}")
    logger.info(f"candidate_sites:      {total_candidate_sites:,}")
    logger.info(f"candidate_examples:   {total_candidate_examples:,}")
    logger.info(f"SNV_examples:         {total_snv_examples:,}")
    logger.info(f"INDEL_examples:       {total_indel_examples:,}")
    logger.info(f"total_time_min:       {elapsed/60:.2f}")

    print()
    print("Graph tumor-only candidate summary")
    print("----------------------------------")
    print(f"nodes_processed:      {total_nodes}")
    print(f"records_processed:    {processed_records}")
    print(f"candidate_sites:      {total_candidate_sites}")
    print(f"candidate_examples:   {total_candidate_examples}")
    print(f"SNV_examples:         {total_snv_examples}")
    print(f"INDEL_examples:       {total_indel_examples}")


def main():
    parser = argparse.ArgumentParser(
        description="Count graph tumor-only somatic candidate examples from .dat/.idx node pileup files."
    )

    parser.add_argument("--dat", required=True, help="Input .dat file")
    parser.add_argument("--idx", required=True, help="Input .idx file")

    parser.add_argument("--min-depth", type=int, default=1)
    parser.add_argument("--min-alt-reads", type=int, default=3)
    parser.add_argument("--min-af", type=float, default=0.10)

    parser.add_argument("--log-file", default=None)
    parser.add_argument("--log-every-nodes", type=int, default=100000)

    args = parser.parse_args()
    logger = setup_logger(args.log_file)
    run(args, logger)


if __name__ == "__main__":
    main()