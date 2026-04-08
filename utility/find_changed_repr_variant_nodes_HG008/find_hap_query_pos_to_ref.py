#!/usr/bin/env python3
import pysam

BAM = "HG008N_curatedv6_250714_polished6.2.hap1.bam"
TARGET_QNAME = "chr1_hap1"
TARGET_QPOS_1BASED = 1195496  # assembly coordinate, 1-based

def get_query_span_1based(rec):
    """
    Return query span covered by this alignment on the query sequence, 1-based inclusive.
    Uses hard clip + aligned query-consuming ops.
    """
    cigartuples = rec.cigartuples
    if cigartuples is None:
        return None

    # BAM CIGAR ops:
    # 0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X
    left_hard = 0
    if cigartuples and cigartuples[0][0] == 5:
        left_hard = cigartuples[0][1]

    query_consuming = {0, 1, 4, 7, 8}
    qlen = sum(length for op, length in cigartuples if op in query_consuming)

    qstart = left_hard + 1
    qend = left_hard + qlen
    return qstart, qend

def map_query_pos_to_ref_pos(rec, target_qpos_1based):
    """
    Map a 1-based query position to a 1-based reference position using aligned pairs.
    Returns None if the query base is not aligned to a reference base.
    """
    # get_aligned_pairs returns 0-based query/ref positions
    for qpos0, rpos0 in rec.get_aligned_pairs(matches_only=False):
        if qpos0 is not None and (qpos0 + 1) == target_qpos_1based:
            if rpos0 is None:
                return None
            return rec.reference_name, rpos0 + 1
    return None

def main():
    bam = pysam.AlignmentFile(BAM, "rb")
    found = False

    for rec in bam.fetch(until_eof=True):
        if rec.query_name != TARGET_QNAME:
            continue

        span = get_query_span_1based(rec)
        if span is None:
            continue

        qstart, qend = span
        if qstart <= TARGET_QPOS_1BASED <= qend:
            found = True
            mapped = map_query_pos_to_ref_pos(rec, TARGET_QPOS_1BASED)

            print("MATCHING_ALIGNMENT")
            print(f"query_name\t{rec.query_name}")
            print(f"flag\t{rec.flag}")
            print(f"reference\t{rec.reference_name}")
            print(f"reference_start_1based\t{rec.reference_start + 1}")
            print(f"mapping_quality\t{rec.mapping_quality}")
            print(f"query_span_1based\t{qstart}-{qend}")
            print(f"cigar\t{rec.cigarstring}")

            if mapped is None:
                print(f"target_query_pos_1based\t{TARGET_QPOS_1BASED}")
                print("mapped_reference_pos\tUNALIGNED_OR_INSERTION")
            else:
                rname, rpos = mapped
                print(f"target_query_pos_1based\t{TARGET_QPOS_1BASED}")
                print(f"mapped_reference_pos\t{rname}:{rpos}")
            print()

    bam.close()

    if not found:
        print(f"No alignment for {TARGET_QNAME}:{TARGET_QPOS_1BASED} was found.")

if __name__ == "__main__":
    main()