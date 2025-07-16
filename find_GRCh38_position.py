#!/usr/bin/env python3
import json
import argparse
import sys
import re
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', stream=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# Helper Functions
# ─────────────────────────────────────────────────────────────────────────────

def load_af_from_vcf(vcf_path, contig):
    """
    Directly loads allele frequencies from a VCF file into a dictionary using pysam.
    """
    try:
        import pysam
    except ImportError:
        logging.critical("pysam is not installed. Please install it to use VCF annotation: pip install pysam")
        sys.exit(1)

    af_map = {}
    logging.info(f"Loading allele frequency data directly from VCF: {vcf_path} for contig '{contig}'")
    try:
        vcf_file = pysam.VariantFile(vcf_path)
        if contig not in vcf_file.header.contigs:
            logging.warning(f"Contig '{contig}' not found in VCF header. No AF annotation will be added.")
            return {}

        for rec in vcf_file.fetch(contig):
            if 'AF' in rec.info and rec.info['AF'] is not None:
                try:
                    position = rec.pos
                    af = float(rec.info['AF'][0])
                    af_map[position] = round(af, 5)
                except (ValueError, TypeError, IndexError):
                    continue
        logging.info(f"Loaded {len(af_map):,} AF records from VCF.")
    except FileNotFoundError:
        logging.error(f"VCF file not found: {vcf_path}. Cannot perform AF annotation.")
    except Exception as e:
        logging.error(f"Could not load AF data from VCF due to an error: {e}.")

    return af_map


# ─────────────────────────────────────────────────────────────────────────────
# GFA Path Extraction Functionality (Optimized)
# ─────────────────────────────────────────────────────────────────────────────

def parse_w_line_path(w_path_str):
    """Parses a W-line path string into oriented segments."""
    oriented_segments = []
    matches = re.findall(r"([<>])([\w.-]+)", w_path_str)
    for orientation_char, seg_id in matches:
        strand = '+' if orientation_char == '>' else '-'
        oriented_segments.append({'id': seg_id, 'strand': strand})
    return oriented_segments


def parse_p_line_path(p_segments_str):
    """Parses a P-line segment string into oriented segments."""
    oriented_segments = []
    segment_entries = p_segments_str.split(',')
    for seg_orient in segment_entries:
        if len(seg_orient) < 2: continue
        seg_id, strand = seg_orient[:-1], seg_orient[-1]
        if strand in ['+', '-']:
            oriented_segments.append({'id': seg_id, 'strand': strand})
    return oriented_segments


def extract_path_info_from_gfa(gfa_file_path, user_grep_pattern, af_data_map=None):
    """
    Optimized single-pass GFA parser. It reads the GFA file only once to find
    the path and load all sequence segments simultaneously.
    """
    if af_data_map is None:
        af_data_map = {}

    segments = {}
    found_path_line = None
    path_pattern_re = re.compile(user_grep_pattern)

    logging.info(f"Starting single-pass GFA parse. File: '{gfa_file_path}', Pattern: '{user_grep_pattern}'")
    try:
        with open(gfa_file_path, 'r') as f_gfa:
            for line in f_gfa:
                # Load S-line (sequence)
                if line.startswith('S\t'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        segments[parts[1]] = {'seq': parts[2], 'len': len(parts[2])}

                # Find the first matching P or W line
                elif (line.startswith('P\t') or line.startswith('W\t')) and not found_path_line:
                    if path_pattern_re.match(line):
                        found_path_line = line.strip()
                        logging.info(f"Found matching path line: {found_path_line[:100]}...")

    except FileNotFoundError:
        return json.dumps({"error": f"GFA file not found: {gfa_file_path}", "status": "error_file_not_found"})
    except Exception as e:
        return json.dumps({"error": f"Error reading GFA file: {e}", "status": "error_gfa_read"})

    if not found_path_line:
        return json.dumps(
            {"error": f"Path not found using pattern: '{user_grep_pattern}'", "status": "error_path_not_found"})

    # --- Process the found path line ---
    line_parts = found_path_line.split('\t')
    actual_record_type = line_parts[0]
    path_source_type, path_identifier_gfa, path_start_offset, final_path_segments_oriented = None, "", 1, []

    if actual_record_type == "W":
        path_source_type = "W"
        path_identifier_gfa = f"W:{line_parts[1]}/{line_parts[2]}/{line_parts[3]}"
        try:
            path_start_offset = int(line_parts[4]) + 1
        except (ValueError, IndexError):
            path_start_offset = 1
        final_path_segments_oriented = parse_w_line_path(line_parts[6])
    elif actual_record_type == "P":
        path_source_type = "P"
        path_identifier_gfa = f"P:{line_parts[1]}"
        final_path_segments_oriented = parse_p_line_path(line_parts[2])

    # --- Assemble final node list with annotations ---
    output_nodes = []
    current_cumulative_pos = path_start_offset
    for seg_info in final_path_segments_oriented:
        seg_id = seg_info['id']
        if seg_id not in segments:
            return json.dumps({"error": f"Segment '{seg_id}' not in S-records.", "status": "error_segment_not_found"})

        node_data = segments[seg_id]
        node_start_pos = current_cumulative_pos
        node_end_pos = current_cumulative_pos + node_data['len']

        node_af_records = []
        if af_data_map:
            for pos in range(node_start_pos, node_end_pos):
                if pos in af_data_map:
                    node_af_records.append({"pos": pos, "af": af_data_map[pos]})

        output_nodes.append({
            "node_id": seg_id,
            "grch38_position_start": node_start_pos,
            "strand_in_path": seg_info['strand'],
            "length": node_data['len'],
            "gnomad_af": node_af_records,
            "sequence": node_data['seq']
        })
        current_cumulative_pos = node_end_pos

    logging.info(f"Successfully processed {len(output_nodes)} nodes for path '{path_identifier_gfa}'.")
    return json.dumps({
        "path_name_input_pattern": user_grep_pattern,
        "path_identifier_gfa": path_identifier_gfa,
        "path_source_type": path_source_type, "nodes": output_nodes, "status": "success"
    }, indent=4)


# ─────────────────────────────────────────────────────────────────────────────
# Main Execution Logic
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extracts node information for a path from a GFA file, with optional direct annotation from a VCF.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("gfa_file", help="Path to the GFA file.")
    parser.add_argument("path_regex_pattern",
                        help="Python regular expression to find the W or P line.\nExample: '^P\\tmyPathName.*'")

    parser.add_argument("--vcf", help="Optional: Path to VCF file (.vcf.bgz) for AF annotation.")
    parser.add_argument("--vcf-contig",
                        help="The contig/chromosome to use from the VCF file (e.g., 'chr1'). Required if --vcf is used.")

    parser.add_argument("--output_file", "-o",
                        help="Optional path to save JSON output. Prints to stdout if not provided.")
    parser.add_argument("--debug", action="store_true", help="Enable debug level logging.")

    args = parser.parse_args()

    if args.vcf and not args.vcf_contig:
        parser.error("--vcf-contig is required when --vcf is provided.")

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    af_data_map = {}
    if args.vcf:
        af_data_map = load_af_from_vcf(args.vcf, args.vcf_contig)

    json_result_str = extract_path_info_from_gfa(args.gfa_file, args.path_regex_pattern, af_data_map)

    if args.output_file:
        try:
            with open(args.output_file, 'w') as f:
                f.write(json_result_str)
            logging.info(f"Output saved to {args.output_file}")
        except IOError as e:
            logging.error(f"Error writing to '{args.output_file}': {e}")
            sys.exit(1)
    else:
        print(json_result_str)
