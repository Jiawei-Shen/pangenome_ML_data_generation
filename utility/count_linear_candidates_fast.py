#!/usr/bin/env python3

import argparse
import logging
import re
import time
import gc
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
        raise ValueError(f"Invalid region: {region}")
    chrom, rest = region.split(":")
    start, end = rest.replace(",", "").split("-")
    return chrom, int(start) - 1, int(end)


def get_fetch_regions(bam, region_arg):
    if region_arg == "autosomes":
        regions = []
        for i in range(1, 23):
            c1, c2 = f"chr{i}", str(i)
            if c1 in bam.references:
                chrom = c1
            elif c2 in bam.references:
                chrom = c2
            else:
                continue
            regions.append((chrom, 0, bam.get_reference_length(chrom)))
        return regions

    if region_arg == "all":
        return [(c, 0, bam.get_reference_length(c)) for c in bam.references]

    return [parse_region(region_arg)]


def passes_read_filters(read, min_mapq):
    return (
        not read.is_unmapped
        and not read.is_secondary
        and not read.is_supplementary
        and not read.is_duplicate
        and read.mapping_quality >= min_mapq
    )


def parse_md_events(read, chrom, min_baseq):
    events = []
    if not read.has_tag("MD") or not read.cigartuples:
        return events

    tokens = MD_TOKEN_RE.findall(read.get_tag("MD"))
    ref_pos = read.reference_start
    query_pos = 0

    cigartuples = read.cigartuples
    ci = 0
    op, rem = cigartuples[0]

    def advance_to_match():
        nonlocal ci, op, rem, ref_pos, query_pos
        while ci < len(cigartuples):
            if op in (0, 7, 8) and rem > 0:
                return True
            if op in (1, 4):
                query_pos += rem
            elif op in (2, 3):
                ref_pos += rem
            ci += 1
            if ci >= len(cigartuples):
                return False
            op, rem = cigartuples[ci]
        return False

    for tok in tokens:
        if tok.isdigit():
            n = int(tok)
            while n > 0:
                if not advance_to_match():
                    return events
                step = min(n, rem)
                ref_pos += step
                query_pos += step
                rem -= step
                n -= step

        elif tok.startswith("^"):
            deleted = tok[1:]
            events.append((chrom, ref_pos, deleted, "-", "DEL"))
            ref_pos += len(deleted)

        else:
            if not advance_to_match():
                return events

            if query_pos >= read.query_length:
                ref_pos += 1
                rem -= 1
                continue

            if read.query_qualities is not None and read.query_qualities[query_pos] < min_baseq:
                ref_pos += 1
                query_pos += 1
                rem -= 1
                continue

            ref_base = tok
            alt_base = read.query_sequence[query_pos].upper()

            if alt_base in {"A", "C", "G", "T"} and alt_base != ref_base:
                events.append((chrom, ref_pos, ref_base, alt_base, "SNV"))

            ref_pos += 1
            query_pos += 1
            rem -= 1

    return events


def parse_insertions(read, chrom, min_baseq):
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
                anchor = max(ref_pos - 1, read.reference_start)
                events.append((chrom, anchor, "-", ins_seq, "INS"))

            query_pos += length

        elif op in (2, 3):
            ref_pos += length
        elif op == 4:
            query_pos += length

    return events


def collect_alt_events_one_region(args, chrom, start, end, logger):
    bam = pysam.AlignmentFile(args.bam, "rb")

    alt_counts = defaultdict(int)
    candidate_positions = set()

    total_reads = 0
    used_reads = 0
    skipped_nm0 = 0
    no_md = 0
    t0 = time.time()

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
            events.extend(parse_md_events(read, read_chrom, args.min_baseq))
        else:
            no_md += 1

        events.extend(parse_insertions(read, read_chrom, args.min_baseq))

        for e_chrom, pos0, ref, alt, typ in events:
            alt_counts[(e_chrom, pos0, ref, alt, typ)] += 1
            candidate_positions.add((e_chrom, pos0))

        if total_reads % args.log_every_reads == 0:
            elapsed = time.time() - t0
            logger.info(
                f"{chrom} event scan: reads={total_reads:,}, used={used_reads:,}, "
                f"skipped_NM0={skipped_nm0:,}, no_MD={no_md:,}, "
                f"raw_alt_alleles={len(alt_counts):,}, positions={len(candidate_positions):,}, "
                f"elapsed={elapsed/60:.2f} min"
            )

    bam.close()

    logger.info(
        f"{chrom} event scan done: reads={total_reads:,}, used={used_reads:,}, "
        f"skipped_NM0={skipped_nm0:,}, no_MD={no_md:,}, "
        f"raw_alt_alleles={len(alt_counts):,}, positions={len(candidate_positions):,}"
    )

    return alt_counts, candidate_positions


def count_depths_one_region(args, candidate_positions, logger, chrom_label):
    bam = pysam.AlignmentFile(args.bam, "rb")
    depth_counts = {}

    positions_by_chrom = defaultdict(list)
    for chrom, pos0 in candidate_positions:
        positions_by_chrom[chrom].append(pos0)

    total = sum(len(set(v)) for v in positions_by_chrom.values())
    done = 0
    t0 = time.time()

    for chrom, positions in positions_by_chrom.items():
        for pos0 in sorted(set(positions)):
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

                    if read.query_qualities is not None and read.query_qualities[qpos] < args.min_baseq:
                        continue

                    depth += 1

            depth_counts[(chrom, pos0)] = depth
            done += 1

            if done % args.log_every_positions == 0:
                elapsed = time.time() - t0
                rate = done / elapsed if elapsed > 0 else 0
                eta = (total - done) / rate if rate > 0 else 0
                logger.info(
                    f"{chrom_label} depth: {done:,}/{total:,}, "
                    f"rate={rate:,.1f} pos/sec, ETA={eta/60:.2f} min"
                )

    bam.close()
    return depth_counts


def summarize_one_region(args, alt_counts, depth_counts):
    candidate_sites = set()
    candidate_examples = 0
    snv_examples = 0
    indel_examples = 0
    raw_snv_alleles = 0
    raw_indel_alleles = 0

    for (chrom, pos0, ref, alt, typ), alt_count in alt_counts.items():
        if typ == "SNV":
            raw_snv_alleles += 1
        else:
            raw_indel_alleles += 1

        depth = depth_counts.get((chrom, pos0), 0)
        if depth < args.min_depth:
            continue

        vaf = alt_count / depth if depth > 0 else 0.0
        if alt_count < args.min_alt_reads or vaf < args.min_af:
            continue

        candidate_examples += 1
        candidate_sites.add((chrom, pos0))

        if typ == "SNV":
            snv_examples += 1
        else:
            indel_examples += 1

    return {
        "candidate_sites": len(candidate_sites),
        "candidate_examples": candidate_examples,
        "SNV_examples": snv_examples,
        "INDEL_examples": indel_examples,
        "raw_snv_alleles": raw_snv_alleles,
        "raw_indel_alleles": raw_indel_alleles,
    }


def add_summary(total, part):
    for k, v in part.items():
        total[k] += v


def run(args):
    logger = setup_logger(args.log_file)
    t0 = time.time()

    logger.info("Starting memory-optimized linear candidate counting")
    logger.info(f"Reference FASTA: {args.ref}")
    logger.info(f"Tumor BAM/CRAM:  {args.bam}")
    logger.info(f"Region:          {args.region}")

    pysam.FastaFile(args.ref).close()

    bam_tmp = pysam.AlignmentFile(args.bam, "rb")
    regions = get_fetch_regions(bam_tmp, args.region)
    bam_tmp.close()

    total_summary = defaultdict(int)

    logger.info(f"Total regions: {len(regions)}")

    for idx, (chrom, start, end) in enumerate(regions, 1):
        chr_t0 = time.time()
        logger.info(f"Processing region {idx}/{len(regions)}: {chrom}:{start + 1}-{end}")

        alt_counts, candidate_positions = collect_alt_events_one_region(
            args, chrom, start, end, logger
        )

        depth_counts = count_depths_one_region(
            args, candidate_positions, logger, chrom
        )

        summary = summarize_one_region(args, alt_counts, depth_counts)
        add_summary(total_summary, summary)

        logger.info(
            f"Region summary {chrom}: "
            f"candidate_sites={summary['candidate_sites']:,}, "
            f"candidate_examples={summary['candidate_examples']:,}, "
            f"SNV={summary['SNV_examples']:,}, "
            f"INDEL={summary['INDEL_examples']:,}, "
            f"time={(time.time() - chr_t0)/60:.2f} min"
        )

        del alt_counts
        del candidate_positions
        del depth_counts
        gc.collect()

    elapsed = time.time() - t0

    logger.info("Final summary")
    for k in [
        "candidate_sites",
        "candidate_examples",
        "SNV_examples",
        "INDEL_examples",
        "raw_snv_alleles",
        "raw_indel_alleles",
    ]:
        logger.info(f"{k}: {total_summary[k]:,}")
    logger.info(f"total_time_min: {elapsed/60:.2f}")

    print()
    print("Memory-optimized linear tumor-only candidate summary")
    print("----------------------------------------------------")
    print(f"candidate_sites:    {total_summary['candidate_sites']}")
    print(f"candidate_examples: {total_summary['candidate_examples']}")
    print(f"SNV_examples:       {total_summary['SNV_examples']}")
    print(f"INDEL_examples:     {total_summary['INDEL_examples']}")
    print(f"raw_snv_alleles:    {total_summary['raw_snv_alleles']}")
    print(f"raw_indel_alleles:  {total_summary['raw_indel_alleles']}")
    print(f"total_time_min:     {elapsed/60:.2f}")


def main():
    parser = argparse.ArgumentParser(
        description="Memory-optimized DeepSomatic-style tumor-only linear candidate counter."
    )

    parser.add_argument("--ref", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument(
        "--region",
        default="autosomes",
        help='Default: autosomes = chr1-chr22. Use "all" or chr1:1-248956422.',
    )

    parser.add_argument("--min-mapq", type=int, default=5)
    parser.add_argument("--min-baseq", type=int, default=10)
    parser.add_argument("--min-depth", type=int, default=1)
    parser.add_argument("--min-alt-reads", type=int, default=3)
    parser.add_argument("--min-af", type=float, default=0.10)
    parser.add_argument("--max-depth", type=int, default=1000000)

    parser.add_argument("--skip-nm0", action="store_true", default=True)
    parser.add_argument("--no-skip-nm0", dest="skip_nm0", action="store_false")

    parser.add_argument("--log-file", default=None)
    parser.add_argument("--log-every-reads", type=int, default=1000000)
    parser.add_argument("--log-every-positions", type=int, default=100000)

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()