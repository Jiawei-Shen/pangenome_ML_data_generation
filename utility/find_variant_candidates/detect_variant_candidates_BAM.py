#!/usr/bin/env python3

import argparse
from collections import Counter
import pysam


def detect_candidates(
    bam_path,
    fasta_path,
    out_tsv,
    summary_txt,
    min_af,
    min_count,
    min_bq,
    min_mq,
):
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path)

    snv_n = 0
    indel_n = 0

    with open(out_tsv, "w") as out:
        out.write("chrom\tpos\tvariant_type\tref\talt\tdepth\talt_count\tAF\n")

        for chrom in bam.references:
            chrom_len = bam.get_reference_length(chrom)

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

                for alt, count in snv_counts.items():
                    af = count / depth if depth > 0 else 0

                    if af > min_af and count > min_count:
                        out.write(
                            f"{chrom}\t{pos1}\tSNV\t{ref_base}\t{alt}\t"
                            f"{depth}\t{count}\t{af:.6f}\n"
                        )
                        snv_n += 1

                for (indel_type, ref, alt), count in indel_counts.items():
                    af = count / depth if depth > 0 else 0

                    if af > min_af and count > min_count:
                        out.write(
                            f"{chrom}\t{pos1}\t{indel_type}\t{ref}\t{alt}\t"
                            f"{depth}\t{count}\t{af:.6f}\n"
                        )
                        indel_n += 1

    bam.close()
    fasta.close()

    total_n = snv_n + indel_n

    with open(summary_txt, "w") as s:
        s.write(f"Input BAM: {bam_path}\n")
        s.write(f"Input FASTA: {fasta_path}\n")
        s.write(f"Output TSV: {out_tsv}\n")
        s.write(f"Minimum AF: {min_af}\n")
        s.write(f"Minimum alt count: {min_count}\n")
        s.write(f"Minimum base quality: {min_bq}\n")
        s.write(f"Minimum mapping quality: {min_mq}\n")
        s.write("\n")
        s.write(f"SNV candidates: {snv_n}\n")
        s.write(f"INDEL candidates: {indel_n}\n")
        s.write(f"Total candidates: {total_n}\n")

    print(f"SNV candidates:   {snv_n}")
    print(f"INDEL candidates: {indel_n}")
    print(f"Total candidates: {total_n}")
    print(f"Summary written to: {summary_txt}")


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
    )


if __name__ == "__main__":
    main()