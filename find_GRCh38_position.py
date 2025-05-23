import re
import argparse
import json
import datetime # For timestamps

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
        dict: path_node_details for more granular path-specific info.
    """
    print("\n[LOG] Initializing GFA parsing function...")
    node_lengths = {}
    path_segments_temp = {}

    # --- Pass 1: Read S lines for lengths and P/W lines for path segments ---
    print("[LOG] Starting Pass 1: Reading GFA segments and path definitions...")
    found_segments_count = 0
    found_potential_paths_count = 0
    try:
        with open(gfa_filepath, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                if not parts: continue
                record_type = parts[0]

                if record_type == 'S':
                    found_segments_count +=1
                    if len(parts) < 3:
                        # This warning is already good
                        # print(f"Warning: Malformed S line at line {line_num}: {line.strip()}")
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
                        # print(f"Warning: Malformed {record_type} line at line {line_num}: {line.strip()}")
                        continue
                    path_name = parts[1]

                    if path_name.startswith(target_path_prefix) and "chr" in path_name.lower():
                        found_potential_paths_count +=1
                        segment_ids_orientations_str = parts[2]
                        if record_type == 'P' or (record_type == 'W' and ',' in segment_ids_orientations_str):
                            segment_definitions = segment_ids_orientations_str.split(',')
                            parsed_segments = []
                            for seg_orient in segment_definitions:
                                if not seg_orient: continue
                                node = seg_orient[:-1]
                                orientation = seg_orient[-1]
                                if orientation not in ['+', '-']:
                                    node = seg_orient
                                    orientation = '?'
                                parsed_segments.append({'id': node, 'orientation': orientation})
                            path_segments_temp[path_name] = parsed_segments
                        elif record_type == 'W':
                             print(f"Info: Single-step W line {path_name} at line {line_num} not processed by default path logic.")
        print(f"[LOG] Pass 1 complete. Found {len(node_lengths)} unique segments with lengths.")
        print(f"[LOG] Identified {len(path_segments_temp)} potential reference paths matching prefix and 'chr' pattern.")

    except FileNotFoundError:
        print(f"Error: GFA file not found at {gfa_filepath}")
        return {}, {}
    except Exception as e:
        print(f"An error occurred during GFA parsing (Pass 1): {e}")
        return {}, {}

    # --- Normalize target_chromosome_arg for comparison ---
    normalized_target_chromosome_filter = None
    if target_chromosome_arg:
        temp_target = str(target_chromosome_arg).lower()
        if not temp_target.startswith("chr"):
            normalized_target_chromosome_filter = f"chr{temp_target}"
        else:
            normalized_target_chromosome_filter = temp_target
        # Log about active filter is now more prominent here

    # --- Pass 2: Calculate coordinates and build output data ---
    print("[LOG] Starting Pass 2: Calculating GRCh38 positions and building output data...")
    if normalized_target_chromosome_filter:
        print(f"[LOG] Filtering for paths matching chromosome: '{normalized_target_chromosome_filter}'")
    else:
        print("[LOG] No specific chromosome filter applied; processing all matched reference paths.")

    output_node_data = {}
    path_node_details_output = {}
    processed_paths_count = 0

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
            # This warning is fine as is, if it occurs.
            # print(f"Warning: Could not reliably extract chromosome name from path '{path_name}'. Skipping this path.")
            continue

        if normalized_target_chromosome_filter and chromosome_name_from_path != normalized_target_chromosome_filter:
            continue # Skip if not the target chromosome

        processed_paths_count += 1
        # Optional: Log each path being processed if verbose mode is desired later
        # print(f"[LOG] Processing path: {path_name} (as {chromosome_name_from_path})")


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

            output_node_data[node_id] = {
                "strand": orientation,
                "length": node_len,
                "grch38_position": grch38_pos_str
            }
            path_specific_node_list.append({
                'id': node_id, 'len': node_len, 'strand': orientation,
                'chr': chromosome_name_from_path, 'start': start_pos, 'end': end_pos
            })
            current_chromosome_position = end_pos + 1

        if path_specific_node_list:
            path_node_details_output[path_name] = path_specific_node_list

    print(f"[LOG] Pass 2 complete. Processed {processed_paths_count} paths matching criteria.")
    print(f"[LOG] Generated mapping data for {len(output_node_data)} unique nodes.")
    print("[LOG] GFA parsing function finished.")
    return output_node_data, path_node_details_output

if __name__ == '__main__':
    start_time = datetime.datetime.now()
    print("===================================================")
    print(" GFA to GRCh38 Position Mapping Script ")
    print(f" Start Timestamp: {start_time.strftime('%Y-%m-%d %H:%M:%S')} ")
    print("===================================================")

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

    print("\n--- Input Parameters ---")
    print(f"GFA File: {gfa_file_to_process}")
    print(f"Path Prefix: '{args.prefix}'")
    if args.chromosome:
        # Normalize and print target chromosome here for clarity before parsing
        temp_target_chr_arg = str(args.chromosome).lower()
        if not temp_target_chr_arg.startswith("chr"):
            normalized_display_chr = f"chr{temp_target_chr_arg}"
        else:
            normalized_display_chr = temp_target_chr_arg
        print(f"Target Chromosome: '{args.chromosome}' (normalized to '{normalized_display_chr}' for filtering)")
    else:
        print("Target Chromosome: Not specified (will process all matching paths)")
    if args.output_json:
        print(f"Output JSON File: {args.output_json}")
    else:
        print("Output JSON File: Not specified (will print sample to console)")
    print("------------------------")


    try:
        # Check file existence before creating dummy, if user provided a path
        if not generate_dummy: # Avoid this check if we are about to create the dummy
            with open(args.gfa_file, 'r') as f_check: pass
    except Exception:
        print(f"\nWarning: Problem accessing GFA file '{args.gfa_file}'.")
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
        gfa_file_to_process = "example.gfa" # Ensure this is updated if dummy is created
        with open(gfa_file_to_process, "w") as f: f.write(dummy_gfa_content)
        generate_dummy = True
        print(f"Dummy GFA file '{gfa_file_to_process}' created for testing.\n")


    node_data_for_json, _ = parse_gfa_to_grch38_positions(
        gfa_file_to_process,
        args.prefix,
        args.chromosome
    )

    print("\n--- Output Summary ---")
    if node_data_for_json:
        if args.output_json:
            print(f"[LOG] Writing output to JSON file: {args.output_json}...")
            try:
                with open(args.output_json, 'w') as f_json:
                    json.dump(node_data_for_json, f_json, indent=4)
                print(f"[LOG] Output successfully saved to: {args.output_json}")
            except IOError as e:
                print(f"\nError: Could not write to JSON file {args.output_json}: {e}")
        else:
            print("[LOG] JSON output not requested. Displaying sample of results to console:")
            count = 0
            for node_id, data in sorted(node_data_for_json.items()):
                print(f"Node {node_id}: {data}")
                count += 1
                if count >= 5 and len(node_data_for_json) > 10: # Print first 5 if many
                    print(f"... and {len(node_data_for_json) - count} more nodes.")
                    break
            if not node_data_for_json: # Should be caught by outer if, but good for consistency
                 print("No node data was generated based on criteria.")
    else:
        print("No node data was generated based on the specified criteria.")

    if generate_dummy:
        print(f"\nNote: A dummy GFA file '{gfa_file_to_process}' was used for this demonstration.")

    end_time = datetime.datetime.now()
    print("---------------------------------------------------")
    print(f" Script Finished. Total execution time: {end_time - start_time} ")
    print(f" End Timestamp: {end_time.strftime('%Y-%m-%d %H:%M:%S')} ")
    print("===================================================")