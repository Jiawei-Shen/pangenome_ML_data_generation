import re
import argparse
import json # Import the json module

def parse_gfa_to_grch38_positions(gfa_filepath, target_path_prefix="GRCh38", target_chromosome_arg=None):
    """
    Parses a GFA file to map node IDs from specified GRCh38 paths
    to their genomic coordinates, strand, and length.
    Assumes paths representing GRCh38 chromosomes start at position 1.
    Can filter for a specific chromosome.

    Args:
        gfa_filepath (str): Path to the GFA file.
        target_path_prefix (str): Prefix to identify relevant GRCh38 paths.
        target_chromosome_arg (str, optional): Specific chromosome to process (e.g., "1", "chr1", "X").
                                              If None, all relevant GRCh38 paths are processed.

    Returns:
        dict: A dictionary where keys are node IDs and values are dicts
              containing {'strand', 'length', 'grch38_position'}.
              Also returns path_node_details for more granular path-specific info.
    """
    node_lengths = {}
    path_segments_temp = {} # path_name -> list of {'id': node_id, 'orientation': orientation}

    # --- Pass 1: Read S lines for lengths and P/W lines for path segments ---
    try:
        with open(gfa_filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                if not parts: continue
                record_type = parts[0]

                if record_type == 'S':
                    if len(parts) < 3:
                        print(f"Warning: Malformed S line at line {line_num}: {line.strip()}")
                        continue
                    node_id = parts[1]
                    sequence = parts[2]
                    length = len(sequence)
                    for tag in parts[3:]:
                        if tag.startswith("LN:i:"):
                            try: length = int(tag.split(":")[2])
                            except (ValueError, IndexError):
                                print(f"Warning: Could not parse LN tag for segment {node_id} at line {line_num}: {tag}")
                            break
                    node_lengths[node_id] = length

                elif record_type == 'P' or record_type == 'W':
                    if len(parts) < 3:
                        print(f"Warning: Malformed {record_type} line at line {line_num}: {line.strip()}")
                        continue
                    path_name = parts[1]

                    # Make "chr" check in path_name case-insensitive for initial filtering
                    if path_name.startswith(target_path_prefix) and "chr" in path_name.lower():
                        segment_ids_orientations_str = parts[2]
                        if record_type == 'P' or (record_type == 'W' and ',' in segment_ids_orientations_str):
                            segment_definitions = segment_ids_orientations_str.split(',')
                            parsed_segments = []
                            for seg_orient in segment_definitions:
                                if not seg_orient: continue
                                node = seg_orient[:-1]
                                orientation = seg_orient[-1]
                                if orientation not in ['+', '-']:
                                    node = seg_orient # Assume it's just node ID if no +/-
                                    orientation = '?' # Unknown orientation
                                parsed_segments.append({'id': node, 'orientation': orientation})
                            path_segments_temp[path_name] = parsed_segments
                        elif record_type == 'W':
                             print(f"Info: Single-step W line {path_name} at line {line_num} not processed by default path logic.")
    except FileNotFoundError:
        print(f"Error: GFA file not found at {gfa_filepath}")
        return {}, {}
    except Exception as e:
        print(f"An error occurred during GFA parsing: {e}")
        return {}, {}

    # --- Normalize target_chromosome_arg for comparison ---
    normalized_target_chromosome = None
    if target_chromosome_arg:
        temp_target = str(target_chromosome_arg).lower()
        if not temp_target.startswith("chr"):
            normalized_target_chromosome = f"chr{temp_target}"
        else:
            normalized_target_chromosome = temp_target
        # This print is now in main
        # print(f"Filtering results for chromosome: {normalized_target_chromosome}")


    # --- Pass 2: Calculate coordinates and build output data ---
    # This dictionary will be structured for the JSON output
    output_node_data = {}
    # This dictionary can still be useful for detailed console output or other needs
    path_node_details_output = {}


    for path_name, segments in path_segments_temp.items():
        chromosome_name_from_path = None
        match_chr = re.search(r"chr([a-zA-Z0-9_.-]+)", path_name, re.IGNORECASE)
        if match_chr:
            chromosome_name_from_path = f"chr{match_chr.group(1).lower()}"
        else:
            match_direct_chr = re.search(r"(chr[a-zA-Z0-9_.-]+)", path_name, re.IGNORECASE)
            if match_direct_chr:
                chromosome_name_from_path = match_direct_chr.group(1).lower()
            else:
                match_num = re.search(r"([XYM]|(?:[0-9]+))([_.-][a-zA-Z0-9]+)?$", path_name)
                if match_num:
                    chromosome_name_from_path = f"chr{match_num.group(1).lower()}"

        if not chromosome_name_from_path:
            # print(f"Warning: Could not reliably extract chromosome name from path '{path_name}'. Skipping this path.")
            continue

        if normalized_target_chromosome and chromosome_name_from_path != normalized_target_chromosome:
            continue

        current_chromosome_position = 1
        path_specific_node_list = []

        for seg_info in segments:
            node_id = seg_info['id']
            orientation = seg_info['orientation']

            if node_id not in node_lengths:
                print(f"Warning: Node ID '{node_id}' (path '{path_name}') not in S lines. Skipping.")
                continue
            node_len = node_lengths[node_id]
            if node_len == 0:
                print(f"Info: Node ID '{node_id}' (path '{path_name}') has length 0. Skipping.")
                continue

            start_pos = current_chromosome_position
            end_pos = current_chromosome_position + node_len - 1
            grch38_pos_str = f"{chromosome_name_from_path}:{start_pos}-{end_pos}"

            # Populate the main output dictionary for JSON
            output_node_data[node_id] = {
                "strand": orientation,
                "length": node_len,
                "grch38_position": grch38_pos_str
            }

            # Populate detailed list for this path (optional, for other uses)
            path_specific_node_list.append({
                'id': node_id, 'len': node_len, 'strand': orientation,
                'chr': chromosome_name_from_path, 'start': start_pos, 'end': end_pos
            })
            current_chromosome_position = end_pos + 1

        if path_specific_node_list:
            path_node_details_output[path_name] = path_specific_node_list # Retain for potential detailed view

    return output_node_data, path_node_details_output

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Parse a GFA file to map node IDs from GRCh38 paths to genomic coordinates, "
                    "strand, and length. Saves output to JSON. "
                    "Assumes GRCh38 paths start at position 1 of the chromosome. Can filter by chromosome."
    )
    parser.add_argument("gfa_file", help="Path to the GFA input file.")
    parser.add_argument(
        "--prefix", default="GRCh38",
        help="Prefix for GFA path names to identify as reference paths (e.g., 'GRCh38'). Default is 'GRCh38'."
    )
    parser.add_argument(
        "--chromosome", "--chr", metavar="CHR", type=str, default=None,
        help="Optional: Specific chromosome to process (e.g., '1', 'chr1', 'X'). Processes all if not set."
    )
    parser.add_argument(
        "--output_json", "--out", metavar="FILE.json", type=str, default=None,
        help="Optional: Path to save the output JSON file."
    )
    args = parser.parse_args()

    gfa_file_to_process = args.gfa_file
    generate_dummy = False

    try:
        with open(args.gfa_file, 'r') as f_check: pass
    except Exception:
        print(f"Warning: Problem accessing GFA file '{args.gfa_file}'.")
        print("Generating a dummy 'example.gfa' for demonstration.\n")
        dummy_gfa_content = """H\tVN:Z:1.0
S\tseg1\tAAAA\tLN:i:4
S\tseg2\tTTTTT\tLN:i:5
S\tseg3\tGG\tLN:i:2
S\tseg4\tCCCCCCCCCC\tLN:i:10
S\tseg5\tNNN\tLN:i:3
S\tseg6\tACGT\tLN:i:4
S\tseg7\tGATTACA\tLN:i:7
S\tseg8\tTTCC\tLN:i:4
P\tGRCh38#0#chr1\tseg1+,seg2+,seg3-\t*
P\tGRCh38.chrX\tseg4+,seg5-\t*
P\tGRCh38_chr2_customName\tseg6+,seg7-\t*
P\tGRCh38.1alt\tseg1+,seg8+\t*
P\tMyAssembly_chr1\tseg1+,seg8-\t*
P\tother_contig_path\tseg4+,seg1-\t*
"""
        gfa_file_to_process = "example.gfa"
        with open(gfa_file_to_process, "w") as f: f.write(dummy_gfa_content)
        generate_dummy = True
        print(f"Dummy GFA file '{gfa_file_to_process}' created.\n")

    print(f"Processing GFA file: {gfa_file_to_process}")
    print(f"Looking for paths starting with '{args.prefix}'.")
    if args.chromosome:
        # Normalize and print target chromosome here for clarity before parsing
        temp_target_chr = str(args.chromosome).lower()
        if not temp_target_chr.startswith("chr"):
            normalized_filter_chr = f"chr{temp_target_chr}"
        else:
            normalized_filter_chr = temp_target_chr
        print(f"Filtering results for chromosome: {normalized_filter_chr}")


    node_data_for_json, _ = parse_gfa_to_grch38_positions( # We primarily need the first returned dict now
        gfa_file_to_process,
        args.prefix,
        args.chromosome
    )

    if node_data_for_json:
        if args.output_json:
            try:
                with open(args.output_json, 'w') as f_json:
                    json.dump(node_data_for_json, f_json, indent=4)
                print(f"\nOutput successfully saved to: {args.output_json}")
            except IOError as e:
                print(f"\nError: Could not write to JSON file {args.output_json}: {e}")
        else:
            print("\n--- Node Data (would be saved to JSON) ---")
            # Print a sample if no output file specified, for quick review
            count = 0
            for node_id, data in sorted(node_data_for_json.items()):
                print(f"Node {node_id}: {data}")
                count += 1
                if count >= 10 and len(node_data_for_json) > 15: # Print first 10 if many
                    print(f"... and {len(node_data_for_json) - count} more nodes.")
                    break
            if not node_data_for_json:
                 print("No node data generated based on criteria.")

    else:
        print("\nNo node data generated based on the criteria.")


    if generate_dummy:
        print(f"\nNote: A dummy GFA file '{gfa_file_to_process}' was created for this demonstration.")
        if args.chromosome is None and args.output_json is None:
            print(f"Try running with '--chromosome 1 --output_json output.json' or '--chr X --out out.json'.")