#!/usr/bin/env python3

import argparse
import time
from collections import Counter
import pysam


def log(msg):
    now = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {msg}", flush=True)


def detect_candidates(
    bam_path,
    fasta_path,
    out_tsv,
    summary_txt,
    min_af,
    min_count,
    min_bq,
    min_mq,
    log_interval,
):
    log("Starting variant candidate detection")
    log(f"Input BAM: {bam_path}")
    log(f"Input FASTA: {fasta_path}")
    log(f"Output TSV: {out_tsv}")
    log(f"Output summary: {summary_txt}")
    log(f"Threshold: AF > {min_af}, count > {min_count}")
    log(f"Filters: min base quality = {min_bq}, min mapping quality = {min_mq}")

    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)

    snv_n = 0
    indel_n = 0
    total_pileup_positions = 0

    start_time = time.time()

    with open(out_tsv, "w") as out:
        out.write("chrom\tpos\tvariant_type\tref\talt\tdepth\talt_count\tAF\n")

        for chrom in bam.references:
            chrom_len = bam.get_reference_length(chrom)
            chrom_start_time = time.time()
            chrom_pileup_positions = 0

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

                total_pileup_positions += 1
                chrom_pileup_positions += 1

                if pos1 - last_logged_pos >= log_interval:
                    pct = (pos1 / chrom_len) * 100 if chrom_len > 0 else 0
                    elapsed = time.time() - chrom_start_time

                    log(
                        f"{chrom}: position {pos1:,}/{chrom_len:,} "
                        f"({pct:.2f}%), pileup sites scanned: {chrom_pileup_positions:,}, "
                        f"SNV: {snv_n:,}, INDEL: {indel_n:,}, "
                        f"elapsed: {elapsed / 60:.2f} min"
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

                    # SNV counting
                    if not pr.is_del and not pr.is_refskip and pr.query_position is not None:
                        qpos = pr.query_position

                        if aln.query_qualities and aln.query_qualities[qpos] < min_bq:
                            continue

                        base = aln.query_sequence[qpos].upper()

                        if base in {"A", "C", "G", "T"}:
                            depth += 1
                            if base != ref_base:
                                snv_counts[base] += 1

                    # INDEL counting
                    if pr.indel != 0 and pr.query_position is not None:
                        qpos = pr.query_position

                        if pr.indel > 0:
                            ins_seq = aln.query_sequence[
                                qpos + 1 : qpos + 1 + pr.indel
                            ].upper()

                            ref = ref_base
                            alt = ref_base + ins_seq
                            key = ("INS", ref, alt)
                            indel_counts[key] += 1

                        elif pr.indel < 0:
                            del_len = abs(pr.indel)
                            deleted_seq = fasta.fetch(
                                chrom, pos0 + 1, pos0 + 1 + del_len
                            ).upper()

                            ref = ref_base + deleted_seq
                            alt = ref_base
                            key = ("DEL", ref, alt)
                            indel_counts[key] += 1

                # Output SNVs
                for alt, count in snv_counts.items():
                    af = count / depth if depth > 0 else 0

                    if af > min_af and count > min_count:
                        out.write(
                            f"{chrom}\t{pos1}\tSNV\t{ref_base}\t{alt}\t"
                            f"{depth}\t{count}\t{af:.6f}\n"
                        )
                        snv_n += 1

                # Output INDELs
                for (indel_type, ref, alt), count in indel_counts.items():
                    af = count / depth if depth > 0 else 0

                    if af > min_af and count > min_count:
                        out.write(
                            f"{chrom}\t{pos1}\t{indel_type}\t{ref}\t{alt}\t"
                            f"{depth}\t{count}\t{af:.6f}\n"
                        )
                        indel_n += 1

            chrom_elapsed = time.time() - chrom_start_time

            log(
                f"Finished {chrom}: scanned {chrom_pileup_positions:,} pileup sites, "
                f"elapsed {chrom_elapsed / 60:.2f} min"
            )
            log(f"Candidates so far: SNV={snv_n:,}, INDEL={indel_n:,}")

    bam.close()
    fasta.close()

    total_n = snv_n + indel_n
    total_elapsed = time.time() - start_time

    with open(summary_txt, "w") as s:
        s.write(f"Input BAM: {bam_path}\n")
        s.write(f"Input FASTA: {fasta_path}\n")
        s.write(f"Output TSV: {out_tsv}\n")
        s.write(f"Minimum AF: {min_af}\n")
        s.write(f"Minimum alt count: {min_count}\n")
        s.write(f"Minimum base quality: {min_bq}\n")
        s.write(f"Minimum mapping quality: {min_mq}\n")
        s.write(f"Log interval: {log_interval}\n")
        s.write("\n")
        s.write(f"Pileup positions scanned: {total_pileup_positions}\n")
        s.write(f"SNV candidates: {snv_n}\n")
        s.write(f"INDEL candidates: {indel_n}\n")
        s.write(f"Total candidates: {total_n}\n")
        s.write(f"Total elapsed minutes: {total_elapsed / 60:.2f}\n")

    log("=" * 80)
    log("Finished variant candidate detection")
    log(f"Pileup positions scanned: {total_pileup_positions:,}")
    log(f"SNV candidates: {snv_n:,}")
    log(f"INDEL candidates: {indel_n:,}")
    log(f"Total candidates: {total_n:,}")
    log(f"Total elapsed: {total_elapsed / 60:.2f} min")
    log(f"Summary written to: {summary_txt}")


def main():
    parser = argparse.ArgumentParser(
        description="Detect SNV and INDEL candidates from BAM using AF and count thresholds."
    )

    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-f", "--fasta", required=True, help="Reference FASTA file")
    parser.add_argument("-o", "--out", required=True, help="Output candidate TSV file")

    parser.add_argument(
        "--summary",
        default=None,
        help="Output summary TXT file. Default: <out>.summary.txt",
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

    summary_txt = args.summary
    if summary_txt is None:
        summary_txt = args.out + ".summary.txt"

    detect_candidates(
        bam_path=args.bam,
        fasta_path=args.fasta,
        out_tsv=args.out,
        summary_txt=summary_txt,
        min_af=args.min_af,
        min_count=args.min_count,
        min_bq=args.min_bq,
        min_mq=args.min_mq,
        log_interval=args.log_interval,
    )


if __name__ == "__main__":
    main()