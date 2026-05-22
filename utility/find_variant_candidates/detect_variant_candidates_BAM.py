#!/usr/bin/env python3

import argparse
import time
from collections import Counter
import pysam


def log(msg):
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {msg}", flush=True)


def write_header(summary_txt, args):
    with open(summary_txt, "w") as s:
        s.write("Variant candidate summary by chromosome\n")
        s.write("=" * 80 + "\n")
        s.write(f"Input BAM: {args.bam}\n")
        s.write(f"Input FASTA: {args.fasta}\n")
        s.write(f"Minimum AF: {args.min_af}\n")
        s.write(f"Minimum alt count: {args.min_count}\n")
        s.write(f"Minimum base quality: {args.min_bq}\n")
        s.write(f"Minimum mapping quality: {args.min_mq}\n")
        s.write(f"Log interval: {args.log_interval}\n")
        s.write("=" * 80 + "\n\n")
        s.write(
            "chrom\tchrom_length\tpileup_positions_scanned\tSNV_candidates\t"
            "INDEL_candidates\tTOTAL_candidates\telapsed_minutes\n"
        )


def append_chrom_summary(
    summary_txt,
    chrom,
    chrom_len,
    pileup_sites,
    snv_n,
    indel_n,
    elapsed_min,
):
    total_n = snv_n + indel_n

    with open(summary_txt, "a") as s:
        s.write(
            f"{chrom}\t{chrom_len}\t{pileup_sites}\t{snv_n}\t"
            f"{indel_n}\t{total_n}\t{elapsed_min:.2f}\n"
        )
        s.flush()


def append_final_summary(
    summary_txt,
    total_pileup_sites,
    total_snv,
    total_indel,
    total_elapsed_min,
):
    total_n = total_snv + total_indel

    with open(summary_txt, "a") as s:
        s.write("\n")
        s.write("=" * 80 + "\n")
        s.write("Final total summary\n")
        s.write("=" * 80 + "\n")
        s.write(f"Total pileup positions scanned: {total_pileup_sites}\n")
        s.write(f"Total SNV candidates: {total_snv}\n")
        s.write(f"Total INDEL candidates: {total_indel}\n")
        s.write(f"Total candidates: {total_n}\n")
        s.write(f"Total elapsed minutes: {total_elapsed_min:.2f}\n")


def detect_candidates(
    bam_path,
    fasta_path,
    summary_txt,
    min_af,
    min_count,
    min_bq,
    min_mq,
    log_interval,
):
    log("Starting variant candidate counting")
    log(f"Input BAM: {bam_path}")
    log(f"Input FASTA: {fasta_path}")
    log(f"Output summary TXT: {summary_txt}")
    log(f"Threshold: AF > {min_af}, count > {min_count}")

    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)

    total_snv = 0
    total_indel = 0
    total_pileup_sites = 0

    start_time = time.time()

    for chrom in bam.references:
        chrom_len = bam.get_reference_length(chrom)
        chrom_start_time = time.time()

        chrom_snv = 0
        chrom_indel = 0
        chrom_pileup_sites = 0

        log("=" * 80)
        log(f"Processing chromosome: {chrom}")
        log(f"Chromosome length: {chrom_len:,} bp")

        last_logged_pos = 0

        for col in bam.pileup(
            chrom,
            0,
            chrom_len,
            truncate=True,
            stepper="samtools",
            min_base_quality=min_bq,
            min_mapping_quality=min_mq,
        ):
            pos0 = col.reference_pos
            pos1 = pos0 + 1

            chrom_pileup_sites += 1
            total_pileup_sites += 1

            if pos1 - last_logged_pos >= log_interval:
                pct = (pos1 / chrom_len) * 100 if chrom_len > 0 else 0
                elapsed = time.time() - chrom_start_time

                log(
                    f"{chrom}: position {pos1:,}/{chrom_len:,} "
                    f"({pct:.2f}%), scanned={chrom_pileup_sites:,}, "
                    f"SNV={chrom_snv:,}, INDEL={chrom_indel:,}, "
                    f"elapsed={elapsed / 60:.2f} min"
                )

                last_logged_pos = pos1

            ref_base = fasta.fetch(chrom, pos0, pos0 + 1).upper()
            if ref_base not in {"A", "C", "G", "T"}:
                continue

            depth = 0
            snv_counts = Counter()
            indel_counts = Counter()

            for pr in col.pileups:
                aln = pr.alignment

                if (
                    aln.is_unmapped
                    or aln.is_duplicate
                    or aln.is_secondary
                    or aln.is_supplementary
                ):
                    continue

                if aln.mapping_quality < min_mq:
                    continue

                if not pr.is_del and not pr.is_refskip and pr.query_position is not None:
                    qpos = pr.query_position

                    if aln.query_qualities and aln.query_qualities[qpos] < min_bq:
                        continue

                    base = aln.query_sequence[qpos].upper()

                    if base in {"A", "C", "G", "T"}:
                        depth += 1
                        if base != ref_base:
                            snv_counts[base] += 1

                if pr.indel != 0 and pr.query_position is not None:
                    qpos = pr.query_position

                    if pr.indel > 0:
                        ins_seq = aln.query_sequence[
                            qpos + 1 : qpos + 1 + pr.indel
                        ].upper()

                        ref = ref_base
                        alt = ref_base + ins_seq
                        indel_counts[("INS", ref, alt)] += 1

                    elif pr.indel < 0:
                        del_len = abs(pr.indel)
                        deleted_seq = fasta.fetch(
                            chrom,
                            pos0 + 1,
                            pos0 + 1 + del_len,
                        ).upper()

                        ref = ref_base + deleted_seq
                        alt = ref_base
                        indel_counts[("DEL", ref, alt)] += 1

            for alt, count in snv_counts.items():
                af = count / depth if depth > 0 else 0

                if af > min_af and count > min_count:
                    chrom_snv += 1

            for key, count in indel_counts.items():
                af = count / depth if depth > 0 else 0

                if af > min_af and count > min_count:
                    chrom_indel += 1

        chrom_elapsed_min = (time.time() - chrom_start_time) / 60

        total_snv += chrom_snv
        total_indel += chrom_indel

        append_chrom_summary(
            summary_txt=summary_txt,
            chrom=chrom,
            chrom_len=chrom_len,
            pileup_sites=chrom_pileup_sites,
            snv_n=chrom_snv,
            indel_n=chrom_indel,
            elapsed_min=chrom_elapsed_min,
        )

        log(
            f"Finished {chrom}: scanned={chrom_pileup_sites:,}, "
            f"SNV={chrom_snv:,}, INDEL={chrom_indel:,}, "
            f"elapsed={chrom_elapsed_min:.2f} min"
        )
        log(f"Updated summary TXT: {summary_txt}")

    bam.close()
    fasta.close()

    total_elapsed_min = (time.time() - start_time) / 60

    append_final_summary(
        summary_txt=summary_txt,
        total_pileup_sites=total_pileup_sites,
        total_snv=total_snv,
        total_indel=total_indel,
        total_elapsed_min=total_elapsed_min,
    )

    log("=" * 80)
    log("Finished all chromosomes")
    log(f"Total pileup positions scanned: {total_pileup_sites:,}")
    log(f"Total SNV candidates: {total_snv:,}")
    log(f"Total INDEL candidates: {total_indel:,}")
    log(f"Total candidates: {total_snv + total_indel:,}")
    log(f"Final summary written to: {summary_txt}")


def main():
    parser = argparse.ArgumentParser(
        description="Count SNV and INDEL candidates from BAM by chromosome."
    )

    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA file")

    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output summary TXT file",
    )

    parser.add_argument(
        "--min-af",
        type=float,
        default=0.1,
        help="Minimum allele frequency. Default: 0.1",
    )

    parser.add_argument(
        "--min-count",
        type=int,
        default=2,
        help="Minimum alternative-supporting read count. Default: 2",
    )

    parser.add_argument(
        "--min-bq",
        type=int,
        default=13,
        help="Minimum base quality. Default: 13",
    )

    parser.add_argument(
        "--min-mq",
        type=int,
        default=20,
        help="Minimum mapping quality. Default: 20",
    )

    parser.add_argument(
        "--log-interval",
        type=int,
        default=1_000_000,
        help="Print progress every N reference bases per chromosome. Default: 1,000,000",
    )

    args = parser.parse_args()

    write_header(args.out, args)

    detect_candidates(
        bam_path=args.bam,
        fasta_path=args.fasta,
        summary_txt=args.out,
        min_af=args.min_af,
        min_count=args.min_count,
        min_bq=args.min_bq,
        min_mq=args.min_mq,
        log_interval=args.log_interval,
    )


if __name__ == "__main__":
    main()