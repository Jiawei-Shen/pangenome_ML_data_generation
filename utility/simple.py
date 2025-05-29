import json
import argparse
import subprocess

def load_vcf_variants(vcf_file, target_chr):
    """Load variants from a specific chromosome in the VCF using bcftools."""
    variants = set()
    try:
        cmd = ["bcftools", "view", "-r", target_chr, "-H", vcf_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.strip().split('\t')
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            for allele in alt.split(','):
                variants.add((chrom, int(pos), ref, allele))
    except Exception as e:
        print(f"Error reading VCF with bcftools: {e}")
    return variants

def load_json_variants(json_file):
    """Extract (chrom, pos, ref, alt) from raw VCF strings in vcf_query_results."""
    variants = set()
    try:
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)

        for node in data.get("nodes", []):
            for line in node.get("vcf_query_results", []):
                if isinstance(line, str):
                    fields = line.strip().split('\t')
                    if len(fields) >= 5:
                        chrom = fields[0]
                        pos = int(fields[1])
                        ref = fields[3]
                        alt = fields[4]
                        for allele in alt.split(','):
                            variants.add((chrom, pos, ref, allele))
    except Exception as e:
        print(f"Error reading JSON file: {e}")
    return variants

def compare_and_output_missing(vcf_variants, json_variants, output_file):
    missing = sorted(vcf_variants - json_variants)
    print(f"Total variants in VCF: {len(vcf_variants)}")
    print(f"Total variants in JSON: {len(json_variants)}")
    print(f"Variants present in VCF but missing in JSON: {len(missing)}")

    with open(output_file, 'w', encoding='utf-8') as out_f:
        for chrom, pos, ref, alt in missing:
            out_f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")

    print(f"Missing variants written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare VCF and JSON vcf_query_results for a given chromosome")
    parser.add_argument("vcf_file", help="Path to the VCF(.gz) file")
    parser.add_argument("json_file", help="Path to the JSON file with vcf_query_results")
    parser.add_argument("--chr", required=True, help="Chromosome to compare (e.g., chr1)")
    parser.add_argument("--out", default="missing_variants.tsv", help="Output file for missing variants")
    args = parser.parse_args()

    vcf_vars = load_vcf_variants(args.vcf_file, args.chr)
    json_vars = load_json_variants(args.json_file)
    compare_and_output_missing(vcf_vars, json_vars, args.out)
