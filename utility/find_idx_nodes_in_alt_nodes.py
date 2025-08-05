#!/usr/bin/env python3

import argparse
import struct
import sys
import re
from collections import defaultdict
import pysam  # We need pysam to read VCF files

# This regex helper is useful for parsing VCF AT tags
_SPLIT = re.compile(r"[><]+")


def split_nodes(trav_str):
    """Parses a vg traversal string like '>1>2>3' into a list of ints."""
    if not trav_str:
        return []
    return [int(tok) for tok in _SPLIT.split(trav_str) if tok]


def load_index(idx_path):
    """Loads node IDs from a binary .idx file. (Unchanged)"""
    node_index = {}
    try:
        with open(idx_path, "rb") as f:
            blocks_num_bytes = f.read(4)
            if not blocks_num_bytes or len(blocks_num_bytes) < 4:
                print(f"Error: Could not read blocks_num from {idx_path}.", file=sys.stderr)
                return {}
            blocks_num, = struct.unpack("<I", blocks_num_bytes)
            for i in range(blocks_num):
                header_data = f.read(22)
                if len(header_data) < 22:
                    break
                block_id, _, _, _, metadata_len = struct.unpack("<I Q I I H", header_data)
                if metadata_len > 0:
                    f.read(metadata_len)
                node_index[block_id] = True
    except FileNotFoundError:
        print(f"Error: IDX file not found at {idx_path}", file=sys.stderr)
        return {}
    except Exception as e:
        print(f"An unexpected error occurred while reading IDX file {idx_path}: {e}", file=sys.stderr)
        return {}
    return node_index


def build_vcf_node_map(vcf_path):
    """
    Scans a VCF file and builds a map from alternate node IDs to their variant info.

    Returns:
        A dictionary mapping {node_id: [(chrom, pos, ref, alt), ...]}
    """
    print("Building map from VCF alternate nodes to variants...", file=sys.stderr)
    node_map = defaultdict(list)
    try:
        vcf = pysam.VariantFile(vcf_path)
        if "AT" not in vcf.header.info:
            print("Error: VCF header lacks required INFO/AT tag.", file=sys.stderr)
            return None
    except FileNotFoundError:
        print(f"Error: VCF file not found at {vcf_path}", file=sys.stderr)
        return None
    except ValueError as e:
        print(f"Error: Could not parse VCF file {vcf_path}. It may be corrupted or not a valid VCF. Details: {e}",
              file=sys.stderr)
        return None

    for rec in vcf:
        if "AT" not in rec.info:
            continue

        ref_nodes = set(split_nodes(rec.info["AT"][0]))
        for i, alt_allele in enumerate(rec.alts):
            # The alt path in AT corresponds to the i-th alt allele
            alt_path_str = rec.info["AT"][i + 1]
            alt_nodes = set(split_nodes(alt_path_str))

            # Find nodes unique to this alternate path
            unique_alt_nodes = alt_nodes - ref_nodes

            for node in unique_alt_nodes:
                variant_info = (rec.contig, rec.pos, rec.ref, alt_allele)
                node_map[node].append(variant_info)

    print(f"Mapped {len(node_map)} unique nodes from VCF.", file=sys.stderr)
    return node_map


def verify_matches(tsv_path, idx_nodes, vcf_node_map):
    """
    Verifies that nodes from a TSV file exist in the index and correspond
    to the correct variant and locus in the VCF file.
    """
    print("Verifying TSV records against index and VCF map...", file=sys.stderr)
    results = {
        "total_nodes_in_tsv": 0,
        "nodes_in_idx": 0,
        "nodes_not_in_vcf_map": 0,
        "locus_variant_matches": 0,
        "locus_variant_mismatches": 0,
    }
    line_num = 0

    try:
        with open(tsv_path, 'r') as f:
            next(f, None)  # Skip header line
            for line in f:
                line_num += 1
                if not line.strip():
                    continue

                parts = line.strip().split()
                if len(parts) < 8:
                    print(f"Warning: Skipping malformed line {line_num}: expected at least 8 columns.", file=sys.stderr)
                    continue

                # Extract data from TSV based on the format from your other script
                # CHROM, POS, ID, TYPE, REF_BASE, REF_NODE, ALT_STR, ALT_NODE(S)
                tsv_chrom, tsv_pos, _, _, tsv_ref, _, tsv_alt, alt_nodes_str = parts[:8]
                tsv_pos = int(tsv_pos)

                try:
                    node_ids = [int(n) for n in alt_nodes_str.split(',')]
                except ValueError:
                    print(f"Warning: Non-integer node ID on line {line_num}: '{alt_nodes_str}'", file=sys.stderr)
                    continue

                for node_id in node_ids:
                    results["total_nodes_in_tsv"] += 1

                    # 1. Check if the node is in the index
                    if node_id not in idx_nodes:
                        continue
                    results["nodes_in_idx"] += 1

                    # 2. Check if the node was found in the VCF AT tags
                    if node_id not in vcf_node_map:
                        results["nodes_not_in_vcf_map"] += 1
                        print(f"Info: Node {node_id} (from TSV line {line_num}) is valid but not found in VCF AT tags.",
                              file=sys.stderr)
                        continue

                    # 3. Verify locus and variant correspondence
                    is_match_found = False
                    for vcf_chrom, vcf_pos, vcf_ref, vcf_alt in vcf_node_map[node_id]:
                        if (vcf_chrom == tsv_chrom and vcf_pos == tsv_pos and
                                vcf_ref == tsv_ref and vcf_alt == tsv_alt):
                            is_match_found = True
                            break  # Found a perfect match

                    if is_match_found:
                        results["locus_variant_matches"] += 1
                    else:
                        results["locus_variant_mismatches"] += 1
                        print(f"Mismatch: Node {node_id} (TSV line {line_num}) did not correspond to VCF entry.",
                              file=sys.stderr)
                        print(f"  - TSV says: {tsv_chrom}:{tsv_pos} {tsv_ref}>{tsv_alt}", file=sys.stderr)
                        print(f"  - VCF has:  {vcf_node_map[node_id]}", file=sys.stderr)

    except FileNotFoundError:
        print(f"Error: TSV file not found at {tsv_path}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred while reading TSV file {tsv_path} on line {line_num}: {e}",
              file=sys.stderr)
        return None

    return results


def main():
    """
    Main function to orchestrate the node verification process.
    """
    parser = argparse.ArgumentParser(
        description="Verify that alternate nodes from a TSV file exist in an index and correspond to variants in a VCF."
    )
    parser.add_argument("tsv_file", help="Path to the input TSV file.")
    parser.add_argument("idx_file", help="Path to the .idx index file.")
    parser.add_argument("vcf_file", help="Path to the VCF file with INFO/AT tags.")
    args = parser.parse_args()

    # 1. Load node IDs from the index file
    idx_data = load_index(args.idx_file)
    if not idx_data:
        sys.exit(1)
    idx_nodes = set(idx_data.keys())
    print(f"Loaded {len(idx_nodes)} unique nodes from the index.", file=sys.stderr)

    # 2. Build the node-to-variant map from the VCF
    vcf_node_map = build_vcf_node_map(args.vcf_file)
    if vcf_node_map is None:
        sys.exit(1)

    # 3. Verify TSV records and get results
    results = verify_matches(args.tsv_file, idx_nodes, vcf_node_map)

    # 4. Output the final summary
    if results:
        print("\n--- Verification Summary ---")
        print(f"Total alternate nodes processed from TSV: {results['total_nodes_in_tsv']}")
        print(f"Nodes found in the .idx file:           {results['nodes_in_idx']}")
        print(
            f"  - Corresponding VCF entries found:      {results['locus_variant_matches'] + results['locus_variant_mismatches']}")
        print(f"    - Perfect locus & variant matches:    {results['locus_variant_matches']}")
        print(f"    - Mismatched locus or variant:      {results['locus_variant_mismatches']}")
        print(f"  - No corresponding VCF entry found:     {results['nodes_not_in_vcf_map']}")


if __name__ == "__main__":
    main()