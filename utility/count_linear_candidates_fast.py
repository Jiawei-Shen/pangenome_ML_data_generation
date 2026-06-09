#!/usr/bin/env python3

import argparse
import logging
import re
import time
from collections import defaultdict

import pysam

MD_TOKEN_RE = re.compile(r"(\d+|\^[A-Z]+|[A-Z])")


def setup_logger(log_file=None):
    logger = logging.getLogger("linear_candidate_counter")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    return logger


def parse_region(region):
    if ":" not in region:
        raise ValueError(
            f"Invalid region: {region}. Use autosomes, all, or chr1:1-248956422"
        )
    chrom, rest = region.split(":")
    start, end = rest.replace(",", "").split("-")
    return chrom, int(start) - 1, int(end)


def get_fetch_regions(bam, region_arg):
    if region_arg == "autosomes":
        regions = []
        for i in range(1, 23):
            c1 = f"chr{i}"
            c2 = str(i)

            if c1 in bam.references:
                chrom = c1
            elif c2 in bam.references:
                chrom = c2
            else:
                continue

            regions.append((chrom, 0, bam.get_reference_length(chrom)))

        if not regions:
            raise RuntimeError("No autosomes chr1-chr22 or 1-22 found in BAM header")

        return regions

    if region_arg == "all":
        return [(c, 0, bam.get_reference_length(c)) for c in bam.references]

    return [parse_region(region_arg)]


def passes_read_filters(read, min_mapq):
    if read.is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if read.is_duplicate:
        return False
    if read.mapping_quality < min_mapq:
        return False
    return True


def parse_md_mismatches_and_deletions(read, ref_name, min_baseq):
    events = []

    if not read.has_tag("MD"):
        return events

    md = read.get_tag("MD")
    tokens = MD_TOKEN_RE.findall(md)

    ref_pos = read.reference_start
    query_pos = 0

    cigartuples = read.cigartuples or []
    if not cigartuples:
        return events

    cigar_index = 0
    cigar_op = cigartuples[0][0]
    cigar_remaining = cigartuples[0][1]

    def advance_to_match_block():
        nonlocal cigar_index, cigar_op, cigar_remaining, query_pos, ref_pos

        while cigar_index < len(cigartuples):
            if cigar_op in (0, 7, 8) and cigar_remaining > 0:
                return True

            if cigar_op in (1, 4):
                query_pos += cigar_remaining
            elif cigar_op in (2, 3):
                ref_pos += cigar_remaining

            cigar_index += 1

            if cigar_index >= len(cigartuples):
                return False

            cigar_op = cigartuples[cigar_index][0]
            cigar_remaining = cigartuples[cigar_index][1]

        return False

    for tok in tokens:
        if tok.isdigit():
            n = int(tok)

            while n > 0:
                if not advance_to_match_block():
                    return events

                step = min(n, cigar_remaining)
                ref_pos += step
                query_pos += step
                cigar_remaining -= step
                n -= step

        elif tok.startswith("^"):
            deleted = tok[1:]
            events.append((ref_name, ref_pos, deleted, "-", "DEL"))
            ref_pos += len(deleted)

        else:
            if not advance_to_match_block():
                return events

            ref_base = tok

            if query_pos >= read.query_length:
                ref_pos += 1
                cigar_remaining -= 1
                continue

            if read.query_qualities is not None:
                if read.query_qualities[query_pos] < min_baseq:
                    ref_pos += 1
                    query_pos += 1
                    cigar_remaining -= 1
                    continue

            alt_base = read.query_sequence[query_pos].upper()

            if alt_base in {"A", "C", "G", "T"} and alt_base != ref_base:
                events.append((ref_name, ref_pos, ref_base, alt_base, "SNV"))

            ref_pos += 1
            query_pos += 1
            cigar_remaining -= 1

    return events


def parse_insertions_from_cigar(read, ref_name, min_baseq):
    events = []

    if read.cigartuples is None:
        return events

    ref_pos = read.reference_start
    query_pos = 0

    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            ref_pos += length
            query_pos += length

        elif op == 1:
            ins_seq = read.query_sequence[query_pos: query_pos + length].upper()
            ok = True

            if read.query_qualities is not None:
                qs = read.query_qualities[query_pos: query_pos + length]
                if len(qs) == 0 or min(qs) < min_baseq:
                    ok = False

            if ok and ins_seq:
                anchor_pos = max(ref_pos - 1, read.reference_start)
                events.append((ref_name, anchor_pos, "-", ins_seq, "INS"))

            query_pos += length

        elif op in (2, 3):
            ref_pos += length

        elif op == 4:
            query_pos += length

        elif op in (5, 6):
            pass

    return events


def collect_alt_events(args, logger):
    bam = pysam.AlignmentFile(args.bam, "rb")
    fetch_regions = get_fetch_regions(bam, args.region)

    logger.info(f"Fetch regions: {len(fetch_regions)}")
    for chrom, start, end in fetch_regions:
        logger.info(f"  {chrom}:{start + 1}-{end}")

    alt_counts = defaultdict(int)
    candidate_positions = set()

    total_reads = 0
    used_reads = 0
    skipped_nm0 = 0
    reads_without_md = 0

    t0 = time.time()

    for chrom, start, end in fetch_regions:
        logger.info(f"Scanning region: {chrom}:{start + 1}-{end}")

        for read in bam.fetch(chrom, start, end):
            total_reads += 1

            if not passes_read_filters(read, args.min_mapq):
                continue

            if args.skip_nm0 and read.has_tag("NM") and read.get_tag("NM") == 0:
                skipped_nm0 += 1
                continue

            used_reads += 1
            read_chrom = bam.get_reference_name(read.reference_id)

            events = []

            if read.has_tag("MD"):
                events.extend(
                    parse_md_mismatches_and_deletions(
                        read, read_chrom, args.min_baseq
                    )
                )
            else:
                reads_without_md += 1

            events.extend(parse_insertions_from_cigar(read, read_chrom, args.min_baseq))

            for event_chrom, pos0, ref, alt, typ in events:
                alt_counts[(event_chrom, pos0, ref, alt, typ)] += 1
                candidate_positions.add((event_chrom, pos0))

            if total_reads % args.log_every_reads == 0:
                elapsed = time.time() - t0
                rate = total_reads / elapsed if elapsed > 0 else 0

                logger.info(
                    f"Event scan progress: reads={total_reads:,}, "
                    f"used_reads={used_reads:,}, skipped_NM0={skipped_nm0:,}, "
                    f"reads_without_MD={reads_without_md:,}, "
                    f"raw_alt_alleles={len(alt_counts):,}, "
                    f"candidate_positions={len(candidate_positions):,}, "
                    f"rate={rate:,.1f} reads/sec, elapsed={elapsed/60:.2f} min"
                )

    bam.close()

    logger.info(
        f"Event scan finished: reads={total_reads:,}, "
        f"used_reads={used_reads:,}, skipped_NM0={skipped_nm0:,}, "
        f"reads_without_MD={reads_without_md:,}, "
        f"raw_alt_alleles={len(alt_counts):,}, "
        f"candidate_positions={len(candidate_positions):,}"
    )

    if reads_without_md > 0:
        logger.warning(
            "Some reads do not have MD tags. SNV/DEL detection may be incomplete. "
            "Run: samtools calmd -b input.bam ref.fa > input.MD.bam"
        )

    return alt_counts, candidate_positions


def count_depths_at_positions(args, candidate_positions, logger):
    bam = pysam.AlignmentFile(args.bam, "rb")

    positions_by_chrom = defaultdict(list)
    for chrom, pos0 in candidate_positions:
        positions_by_chrom[chrom].append(pos0)

    depth_counts = {}
    done = 0
    total = len(candidate_positions)
    t0 = time.time()

    for chrom in sorted(positions_by_chrom.keys()):
        positions = sorted(set(positions_by_chrom[chrom]))
        logger.info(f"Depth calculation for {chrom}: {len(positions):,} positions")

        for pos0 in positions:
            depth = 0

            for col in bam.pileup(
                chrom,
                pos0,
                pos0 + 1,
                truncate=True,
                stepper="samtools",
                min_mapping_quality=args.min_mapq,
                ignore_overlaps=False,
                ignore_orphans=False,
                max_depth=args.max_depth,
            ):
                if col.reference_pos != pos0:
                    continue

                for pr in col.pileups:
                    read = pr.alignment

                    if not passes_read_filters(read, args.min_mapq):
                        continue
                    if pr.is_refskip:
                        continue

                    if pr.is_del:
                        depth += 1
                        continue

                    qpos = pr.query_position
                    if qpos is None:
                        continue

                    if read.query_qualities is not None:
                        if read.query_qualities[qpos] < args.min_baseq:
                            continue

                    depth += 1

            depth_counts[(chrom, pos0)] = depth
            done += 1

            if done % args.log_every_positions == 0:
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total - done) / rate if rate > 0 else 0

                logger.info(
                    f"Depth progress: positions={done:,}/{total:,}, "
                    f"rate={rate:,.1f} pos/sec, "
                    f"ETA={eta/60:.2f} min, elapsed={elapsed/60:.2f} min"
                )

    bam.close()
    logger.info(f"Depth calculation finished: {len(depth_counts):,} positions")
    return depth_counts


def summarize_candidates(args, alt_counts, depth_counts, logger):
    candidate_sites = set()
    candidate_examples = 0
    snv_examples = 0
    indel_examples = 0
    raw_snv_alleles = 0
    raw_indel_alleles = 0

    out = open(args.output_tsv, "w") if args.output_tsv else None
    if out:
        out.write("chrom\tpos\tref\talt\ttype\tdepth\talt_count\tvaf\n")

    for (chrom, pos0, ref, alt, typ), alt_count in alt_counts.items():
        if typ == "SNV":
            raw_snv_alleles += 1
        else:
            raw_indel_alleles += 1

        depth = depth_counts.get((chrom, pos0), 0)
        if depth < args.min_depth:
            continue

        vaf = alt_count / depth if depth > 0 else 0.0

        if alt_count < args.min_alt_reads:
            continue
        if vaf < args.min_af:
            continue

        candidate_examples += 1
        candidate_sites.add((chrom, pos0))

        if typ == "SNV":
            snv_examples += 1
        else:
            indel_examples += 1

        if out:
            out.write(
                f"{chrom}\t{pos0 + 1}\t{ref}\t{alt}\t{typ}\t"
                f"{depth}\t{alt_count}\t{vaf:.6f}\n"
            )

    if out:
        out.close()

    summary = {
        "candidate_sites": len(candidate_sites),
        "candidate_examples": candidate_examples,
        "SNV_examples": snv_examples,
        "INDEL_examples": indel_examples,
        "raw_snv_alleles": raw_snv_alleles,
        "raw_indel_alleles": raw_indel_alleles,
    }

    logger.info("Candidate summary")
    for k, v in summary.items():
        logger.info(f"{k}: {v:,}")

    return summary


def run(args):
    logger = setup_logger(args.log_file)
    t0 = time.time()

    logger.info("Starting fast linear tumor-only candidate counting")
    logger.info(f"Reference FASTA: {args.ref}")
    logger.info(f"Tumor BAM/CRAM:  {args.bam}")
    logger.info(f"Region:          {args.region}")
    logger.info(
        f"Parameters: min_mapq={args.min_mapq}, min_baseq={args.min_baseq}, "
        f"min_depth={args.min_depth}, min_alt_reads={args.min_alt_reads}, "
        f"min_af={args.min_af}, skip_nm0={args.skip_nm0}"
    )

    fasta = pysam.FastaFile(args.ref)
    fasta.close()

    alt_counts, candidate_positions = collect_alt_events(args, logger)
    depth_counts = count_depths_at_positions(args, candidate_positions, logger)
    summary = summarize_candidates(args, alt_counts, depth_counts, logger)

    elapsed = time.time() - t0
    logger.info(f"total_time_min: {elapsed/60:.2f}")

    print()
    print("Fast linear tumor-only candidate summary")
    print("----------------------------------------")
    print(f"candidate_sites:    {summary['candidate_sites']}")
    print(f"candidate_examples: {summary['candidate_examples']}")
    print(f"SNV_examples:       {summary['SNV_examples']}")
    print(f"INDEL_examples:     {summary['INDEL_examples']}")
    print(f"raw_snv_alleles:    {summary['raw_snv_alleles']}")
    print(f"raw_indel_alleles:  {summary['raw_indel_alleles']}")
    print(f"total_time_min:     {elapsed/60:.2f}")


def main():
    parser = argparse.ArgumentParser(
        description="Fast DeepSomatic-style tumor-only linear candidate counter."
    )

    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--bam", required=True, help="Tumor BAM/CRAM")

    parser.add_argument(
        "--region",
        default="autosomes",
        help='Default: autosomes = chr1-chr22. Use "all" for all references, or chr1:1-248956422.',
    )

    parser.add_argument("--min-mapq", type=int, default=5)
    parser.add_argument("--min-baseq", type=int, default=10)
    parser.add_argument("--min-depth", type=int, default=1)
    parser.add_argument("--min-alt-reads", type=int, default=3)
    parser.add_argument("--min-af", type=float, default=0.10)
    parser.add_argument("--max-depth", type=int, default=1000000)

    parser.add_argument("--skip-nm0", action="store_true", default=True)
    parser.add_argument("--no-skip-nm0", dest="skip_nm0", action="store_false")

    parser.add_argument("--output-tsv", default=None)
    parser.add_argument("--log-file", default=None)

    parser.add_argument("--log-every-reads", type=int, default=1000000)
    parser.add_argument("--log-every-positions", type=int, default=100000)

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()