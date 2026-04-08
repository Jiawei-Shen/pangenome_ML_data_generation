#!/usr/bin/env python3
"""
Find where an assembly query position maps in a BAM alignment, with extra context.

Outputs for each matching alignment:
1. the local query sequence context around the target position
2. the local CIGAR fragment around the target position, plus the previous 20 bp of query context
3. if the target falls inside an insertion/query-only segment, report the insertion anchor position on the reference
4. the sequence at the target position (single base by default)

Designed for assembly-to-reference BAMs like dipcall hap1/hap2 BAMs.

Example:
    python find_hap_query_pos_to_ref_context.py \
        --bam HG008N_curatedv6_250714_polished6.2.hap1.bam \
        --query-name chr1_hap1 \
        --query-pos 1195496
"""

import argparse
import pysam
from typing import List, Optional, Tuple


# BAM CIGAR op codes
CIGAR_M = 0
CIGAR_I = 1
CIGAR_D = 2
CIGAR_N = 3
CIGAR_S = 4
CIGAR_H = 5
CIGAR_P = 6
CIGAR_EQ = 7
CIGAR_X = 8

CIGAR_CODE_TO_CHAR = {
    CIGAR_M: "M",
    CIGAR_I: "I",
    CIGAR_D: "D",
    CIGAR_N: "N",
    CIGAR_S: "S",
    CIGAR_H: "H",
    CIGAR_P: "P",
    CIGAR_EQ: "=",
    CIGAR_X: "X",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--query-name", required=True, help="Query/contig name, e.g. chr1_hap1")
    parser.add_argument("--query-pos", required=True, type=int, help="1-based query position")
    parser.add_argument(
        "--context-bp",
        type=int,
        default=20,
        help="Number of query bases to show before the target position (default: 20)",
    )
    return parser.parse_args()


def leading_clipped_bases(cigartuples: Optional[List[Tuple[int, int]]]) -> int:
    """
    Return total leading clipped bases on query from the left side.
    Includes leading H and/or S if present in order.
    """
    if not cigartuples:
        return 0

    idx = 0
    total = 0
    while idx < len(cigartuples) and cigartuples[idx][0] in (CIGAR_H, CIGAR_S):
        total += cigartuples[idx][1]
        idx += 1
    return total


def aligned_query_span_1based(rec: pysam.AlignedSegment) -> Optional[Tuple[int, int]]:
    """
    Return the 1-based inclusive query span that is actually represented by this alignment record.
    Leading hard/soft clipping shifts the query start.
    Only aligned/query-consuming ops count toward the aligned span: M, I, =, X.
    """
    cigartuples = rec.cigartuples
    if cigartuples is None:
        return None

    qstart = leading_clipped_bases(cigartuples) + 1
    q_aligned_len = sum(length for op, length in cigartuples if op in (CIGAR_M, CIGAR_I, CIGAR_EQ, CIGAR_X))
    qend = qstart + q_aligned_len - 1
    return qstart, qend


def build_op_intervals(
    rec: pysam.AlignedSegment,
) -> List[dict]:
    """
    Walk through the CIGAR and build per-op intervals on query/reference.

    Returns a list of dicts, each containing:
      - op, op_char, length
      - query_start, query_end (1-based, inclusive) if query-consuming else None
      - ref_start, ref_end (1-based, inclusive) if ref-consuming else None
      - insertion_anchor_ref (1-based reference anchor position for I, else None)
      - insertion_anchor_between (tuple(left_anchor, right_anchor_or_None)) for I, else None
    """
    cigartuples = rec.cigartuples
    if cigartuples is None:
        return []

    ref_name = rec.reference_name
    ref_pos = rec.reference_start + 1  # 1-based
    q_pos = leading_clipped_bases(cigartuples) + 1  # first aligned query position (1-based)

    intervals = []

    for op, length in cigartuples:
        op_char = CIGAR_CODE_TO_CHAR.get(op, "?")

        query_consuming = op in (CIGAR_M, CIGAR_I, CIGAR_S, CIGAR_EQ, CIGAR_X)
        ref_consuming = op in (CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_EQ, CIGAR_X)

        q_start = q_pos if query_consuming and op != CIGAR_S else None
        q_end = (q_pos + length - 1) if query_consuming and op != CIGAR_S else None

        r_start = ref_pos if ref_consuming else None
        r_end = (ref_pos + length - 1) if ref_consuming else None

        insertion_anchor_ref = None
        insertion_anchor_between = None
        if op == CIGAR_I:
            left_anchor = ref_pos - 1
            right_anchor = ref_pos
            insertion_anchor_ref = left_anchor
            insertion_anchor_between = (left_anchor, right_anchor)

        intervals.append(
            {
                "ref_name": ref_name,
                "op": op,
                "op_char": op_char,
                "length": length,
                "query_start": q_start,
                "query_end": q_end,
                "ref_start": r_start,
                "ref_end": r_end,
                "insertion_anchor_ref": insertion_anchor_ref,
                "insertion_anchor_between": insertion_anchor_between,
            }
        )

        if query_consuming:
            q_pos += length
        if ref_consuming:
            ref_pos += length

    return intervals


def find_target_interval(intervals: List[dict], query_pos_1based: int) -> Optional[dict]:
    """
    Return the CIGAR interval containing the query position.
    Only intervals with query coverage are considered.
    """
    for item in intervals:
        q_start = item["query_start"]
        q_end = item["query_end"]
        if q_start is not None and q_end is not None and q_start <= query_pos_1based <= q_end:
            return item
    return None


def query_pos_to_ref_pos(rec: pysam.AlignedSegment, query_pos_1based: int) -> Optional[Tuple[str, int]]:
    """
    Map a single 1-based query position to a single 1-based reference position.
    Returns None if the query base is query-only (e.g. insertion/softclip).
    """
    for qpos0, rpos0 in rec.get_aligned_pairs(matches_only=False):
        if qpos0 is not None and (qpos0 + 1) == query_pos_1based:
            if rpos0 is None:
                return None
            return rec.reference_name, rpos0 + 1
    return None


def get_query_context(
    rec: pysam.AlignedSegment,
    query_pos_1based: int,
    context_bp: int,
) -> Tuple[str, str, str]:
    """
    Return:
      - previous up-to-context_bp query bases
      - target base
      - next up-to-context_bp query bases
    Query positions are interpreted on the full query sequence represented by the BAM record:
    aligned query bases plus leading soft clips if present.
    Hard-clipped bases are not present in query_sequence and cannot be recovered.
    """
    seq = rec.query_sequence
    cigartuples = rec.cigartuples
    if seq is None or cigartuples is None:
        return ("", "", "")

    left_soft = 0
    idx = 0
    while idx < len(cigartuples) and cigartuples[idx][0] == CIGAR_H:
        idx += 1
    if idx < len(cigartuples) and cigartuples[idx][0] == CIGAR_S:
        left_soft = cigartuples[idx][1]

    # rec.query_sequence includes soft-clipped sequence but not hard-clipped sequence.
    # So query coordinate inside rec.query_sequence starts at 1 + leading_soft.
    seq_pos_1based = query_pos_1based - leading_clipped_bases(cigartuples) + left_soft

    if seq_pos_1based < 1 or seq_pos_1based > len(seq):
        return ("", "", "")

    target_idx0 = seq_pos_1based - 1
    prev_seq = seq[max(0, target_idx0 - context_bp): target_idx0]
    target_base = seq[target_idx0]
    next_seq = seq[target_idx0 + 1: target_idx0 + 1 + context_bp]
    return (prev_seq, target_base, next_seq)


def format_local_cigar_fragment(target_item: dict) -> str:
    """
    For the target position, return the specific covering CIGAR fragment only.
    Example: 70I or 105M
    """
    return f"{target_item['length']}{target_item['op_char']}"


def flag_labels(flag: int) -> str:
    labels = []
    if flag & 0x1:
        labels.append("PAIRED")
    if flag & 0x2:
        labels.append("PROPER_PAIR")
    if flag & 0x4:
        labels.append("UNMAP")
    if flag & 0x8:
        labels.append("MUNMAP")
    if flag & 0x10:
        labels.append("REVERSE")
    if flag & 0x20:
        labels.append("MREVERSE")
    if flag & 0x40:
        labels.append("READ1")
    if flag & 0x80:
        labels.append("READ2")
    if flag & 0x100:
        labels.append("SECONDARY")
    if flag & 0x200:
        labels.append("QCFAIL")
    if flag & 0x400:
        labels.append("DUP")
    if flag & 0x800:
        labels.append("SUPPLEMENTARY")
    return ",".join(labels) if labels else "PRIMARY"


def print_alignment_report(
    rec: pysam.AlignedSegment,
    query_pos_1based: int,
    context_bp: int,
) -> None:
    intervals = build_op_intervals(rec)
    target_item = find_target_interval(intervals, query_pos_1based)
    if target_item is None:
        return

    mapped = query_pos_to_ref_pos(rec, query_pos_1based)
    prev_seq, target_base, next_seq = get_query_context(rec, query_pos_1based, context_bp)
    qspan = aligned_query_span_1based(rec)

    print("MATCHING_ALIGNMENT")
    print(f"query_name\t{rec.query_name}")
    print(f"flag\t{rec.flag}")
    print(f"flag_label\t{flag_labels(rec.flag)}")
    print(f"reference\t{rec.reference_name}")
    print(f"reference_start_1based\t{rec.reference_start + 1}")
    print(f"mapping_quality\t{rec.mapping_quality}")
    if qspan is not None:
        print(f"query_span_1based\t{qspan[0]}-{qspan[1]}")

    # 1. local context and corresponding CIGAR fragment
    print(f"target_query_pos_1based\t{query_pos_1based}")
    print(f"prev_{context_bp}bp_query_seq\t{prev_seq if prev_seq else '.'}")
    print(f"target_base_or_bp\t{target_base if target_base else '.'}")
    print(f"next_{context_bp}bp_query_seq\t{next_seq if next_seq else '.'}")
    print(f"covering_cigar_fragment\t{format_local_cigar_fragment(target_item)}")

    # 2. if insertion, output insertion anchor position
    if target_item["op"] == CIGAR_I:
        left_anchor, right_anchor = target_item["insertion_anchor_between"]
        if left_anchor >= 1:
            print(f"insertion_anchor_reference\t{rec.reference_name}:{left_anchor}")
        else:
            print("insertion_anchor_reference\tREFERENCE_START_BEFORE_1")
        if right_anchor is not None:
            print(f"insertion_anchor_between\t{rec.reference_name}:{left_anchor}|{rec.reference_name}:{right_anchor}")
        else:
            print("insertion_anchor_between\tUNKNOWN")
        print("mapped_reference_pos\tUNALIGNED_OR_INSERTION")
    else:
        if mapped is None:
            print("mapped_reference_pos\tUNALIGNED_OR_INSERTION")
        else:
            print(f"mapped_reference_pos\t{mapped[0]}:{mapped[1]}")

    # extra useful fields
    print(f"target_cigar_op\t{target_item['op_char']}")
    print(f"target_cigar_query_span\t{target_item['query_start']}-{target_item['query_end']}")
    if target_item["ref_start"] is not None and target_item["ref_end"] is not None:
        print(f"target_cigar_ref_span\t{target_item['ref_start']}-{target_item['ref_end']}")
    else:
        print("target_cigar_ref_span\t.")
    print()


def main() -> None:
    args = parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")
    found = False

    for rec in bam.fetch(until_eof=True):
        if rec.query_name != args.query_name:
            continue

        qspan = aligned_query_span_1based(rec)
        if qspan is None:
            continue

        qstart, qend = qspan
        if qstart <= args.query_pos <= qend:
            intervals = build_op_intervals(rec)
            target_item = find_target_interval(intervals, args.query_pos)
            if target_item is None:
                continue

            found = True
            print_alignment_report(rec, args.query_pos, args.context_bp)

    bam.close()

    if not found:
        print(f"No alignment covering {args.query_name}:{args.query_pos} was found.")


if __name__ == "__main__":
    main()