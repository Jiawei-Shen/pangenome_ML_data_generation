import pysam
import json
import time


def read_vcf_pysam(file_path):
    vcf_file = pysam.VariantFile(file_path)
    variants = []

    for record in vcf_file:
        variant_data = {
            "chromosome": record.chrom,
            "position": record.pos,
            "reference": record.ref,
            "alternative": record.alts,
            "id": record.id,
            "quality": record.qual,
            "samples": {sample: record.samples[sample].alleles for sample in record.samples}
        }
        variants.append(variant_data)

    return variants


if __name__ == "__main__":
    start_time = time.time()  # Start time

    vcf_path = "./Genome_Bottle_VCF/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    vcf_data = read_vcf_pysam(vcf_path)

    # Save to JSON file
    with open("HG005_variants.json", "w") as f:
        json.dump(vcf_data, f, indent=4)

    # Print first 5 entries
    print(vcf_data[:5])

    end_time = time.time()  # End time

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.6f} seconds")
