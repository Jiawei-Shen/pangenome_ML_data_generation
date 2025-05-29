import json
import argparse
import subprocess

def load_vcf_variants(vcf_file):
    """Load variants from a VCF using bcftools."""
    variants = set()
    try:
        cmd = ["bcftools", "view", "-H", vcf_file]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            fields = line.strip().split('\t')
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            for allele in alt.split(','):  # handle multi-allelics
                variants.add((chrom, int(pos), ref, allele))
    except Exception as e:
        print(f"Error reading VCF with bcftools: {e}")
    return variants

def load_json_variants(json_file):
    """Extract variants from JSON file under vcf_query_results."""
    variants = set()
    try:
        with open(json_file, 'r', encoding='utf-8') as f:
            data = json.load(f)
        nodes = data.get("nodes", [])
        for node in nodes:
            for result in node.get("vcf_query_results", []):
                chrom = result.get("chrom") or result.get("CHROM")
                pos = result.get("pos") or result.get("POS")
                ref = result.get("ref") or result.get("REF")
                alt = result.get("alt") or result.get("ALT")
                if chrom and pos and ref and alt:
                    for allele in alt.split(','):
                        variants.add((str(chrom), int(pos), str(ref), str(allele)))
    except Exception as e:
        print(f"Error reading JSON file: {e}")
    return variants

def compare_variants(vcf_variants, json_variants):
    missing = vcf_variants - json_variants
    print(f"\nTotal variants in VCF: {len(vcf_variants)}")
    print(f"Total variants in JSON: {len(json_variants)}")
    print(f"Variants present in VCF but missing in JSON: {len(missing)}")

    if missing:
        print("\nExamples (up to 10):")
        for i, var in enumerate(sorted(missing)[:10], 1):
            print(f"{i}. {var}")
    return missing

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare variants between VCF and JSON vcf_query_results")
    parser.add_argument("vcf_file", help="Path to the VCF or VCF.GZ file")
    parser.add_argument("json_file", help="Path to the JSON file with vcf_query_results")
    args = parser.parse_args()

    vcf_vars = load_vcf_variants(args.vcf_file)
    json_vars = load_json_variants(args.json_file)
    compare_variants(vcf_vars, json_vars)
