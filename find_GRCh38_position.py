import json
import argparse
import sys


def extract_path_info_from_gfa(gfa_file_path, target_path_name_input):
    """
    Extracts nodes from a specified path in a GFA file, calculates their
    cumulative positions, and outputs the information in JSON format.

    Args:
        gfa_file_path (str): The path to the GFA file.
        target_path_name_input (str): The name of the path to extract.
                                     The script will try to match this directly,
                                     and if it contains tabs, it will also try
                                     parts of it (e.g., the last part like 'chr1',
                                     the first part like 'GRCh38') or an
                                     underscored version (e.g., 'GRCh38_0_chr1').

    Returns:
        str: A JSON string containing the node information for the target path,
             or an error message if the path is not found or an issue occurs.
    """
    segments = {}  # To store segment sequences and lengths {segment_id: {'seq': 'ACGT', 'len': 4}}
    path_segments_oriented = []
    found_path = False
    matched_path_name_in_gfa = ""  # To store the name of the path as found in GFA

    potential_path_names_to_try = [target_path_name_input]
    if '\t' in target_path_name_input:
        parts = target_path_name_input.split('\t')
        if len(parts) > 0:
            # Add last part (e.g., "chr1" from "GRCh38\t0\tchr1")
            if parts[-1] not in potential_path_names_to_try:
                potential_path_names_to_try.append(parts[-1])
            # Add first part (e.g., "GRCh38")
            if parts[0] not in potential_path_names_to_try:
                potential_path_names_to_try.append(parts[0])
        # Add version with underscores (e.g., "GRCh38_0_chr1")
        underscored_name = target_path_name_input.replace('\t', '_')
        if underscored_name not in potential_path_names_to_try:
            potential_path_names_to_try.append(underscored_name)

    # Ensure the original input is first if it wasn't already, then unique-ify while preserving order somewhat
    # (More robust way would be to use a list and check `not in` before appending)
    # Simple approach:
    unique_potential_paths = []
    for name in potential_path_names_to_try:
        if name not in unique_potential_paths:
            unique_potential_paths.append(name)
    potential_path_names_to_try = unique_potential_paths

    try:
        with open(gfa_file_path, 'r') as f_gfa:
            for line in f_gfa:
                clean_line = line.strip()
                if not clean_line:  # Skip empty lines
                    continue
                parts = clean_line.split('\t')
                if not parts:  # Skip lines that become empty after split
                    continue
                record_type = parts[0]

                if record_type == 'S':
                    if len(parts) < 3:
                        # print(f"Warning: Skipping malformed S record: {line.strip()}", file=sys.stderr)
                        continue
                    segment_id = parts[1]
                    sequence = parts[2]
                    segments[segment_id] = {'seq': sequence, 'len': len(sequence)}
                elif record_type == 'P':
                    if len(parts) < 3:
                        # print(f"Warning: Skipping malformed P record: {line.strip()}", file=sys.stderr)
                        continue
                    path_id_from_file = parts[1]

                    if path_id_from_file in potential_path_names_to_try:
                        matched_path_name_in_gfa = path_id_from_file
                        found_path = True
                        # Ensure parts[2] (segment names) exists and is not empty
                        if parts[2]:
                            segment_names_with_orientations = parts[2].split(',')
                            for seg_orient in segment_names_with_orientations:
                                if len(seg_orient) < 2:  # segment name + orientation
                                    # print(f"Warning: Skipping malformed segment entry '{seg_orient}' in path '{matched_path_name_in_gfa}'.", file=sys.stderr)
                                    continue
                                path_segments_oriented.append(
                                    {'id': seg_orient[:-1], 'strand': seg_orient[-1]}
                                )
                        else:  # Path found but has no segments listed in parts[2]
                            # print(f"Warning: Path '{matched_path_name_in_gfa}' found but has an empty segment list in GFA.", file=sys.stderr)
                            pass  # Path is found, but will result in empty node list later if this is the case.
                        break  # Found the path, stop reading P lines for this specific path

    except FileNotFoundError:
        return json.dumps({"error": f"File not found: {gfa_file_path}", "status": "error"})
    except Exception as e:
        return json.dumps({"error": f"Error reading GFA file '{gfa_file_path}': {str(e)}", "status": "error"})

    if not segments:
        return json.dumps({"error": "No segments (S records) found in the GFA file.", "status": "error"})

    if not found_path:
        error_message_main = f"Path matching '{target_path_name_input}' not found in the GFA file."

        # Create a list of other names that were tried, excluding the original input if it's the only one
        tried_additionally = [f"'{name}'" for name in potential_path_names_to_try if name != target_path_name_input]

        if tried_additionally:
            error_message_detail = f" Also tried interpreting it as: {', '.join(tried_additionally)}."
            full_error_message = error_message_main + error_message_detail
        else:
            full_error_message = error_message_main
        return json.dumps({"error": full_error_message, "status": "error"})

    if not path_segments_oriented and found_path:  # Path was declared but no segments in it
        # This can happen if P line is "P\tpath_name\t\t*" (empty segment list)
        return json.dumps({
            "message": f"Path '{matched_path_name_in_gfa}' was found but contains no segments.",
            "path_name": matched_path_name_in_gfa,
            "nodes": [],
            "status": "success_empty_path"
        })

    output_nodes = []
    current_cumulative_pos = 0

    for seg_info in path_segments_oriented:
        seg_id = seg_info['id']
        strand = seg_info['strand']

        if seg_id not in segments:
            return json.dumps({
                "error": f"Segment '{seg_id}' from path '{matched_path_name_in_gfa}' not found in S records.",
                "status": "error"
            })

        node_data = segments[seg_id]
        node_len = node_data['len']
        node_seq = node_data['seq']

        output_nodes.append({
            "node_id": seg_id,
            "grch38_position_start": current_cumulative_pos,  # Assuming this is generic chromosome name
            "strand_in_path": strand,
            "sequence": node_seq,
            "length": node_len
        })
        current_cumulative_pos += node_len

    # Successfully processed
    return json.dumps({
        "path_name_input": target_path_name_input,
        "path_name_found_in_gfa": matched_path_name_in_gfa,
        "nodes": output_nodes,
        "status": "success"
    }, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extracts node information from a specified path in a GFA file and outputs JSON. \n"
                    "The 'grch38_position_start' in the output is a cumulative position based on node lengths in the specified path.",
        formatter_class=argparse.RawTextHelpFormatter  # Allows for newlines in help text
    )
    parser.add_argument(
        "gfa_file",
        help="Path to the GFA file."
    )
    parser.add_argument(
        "path_name",
        help="Target path name to extract. Examples:\n"
             "- 'chr1' (if GFA P record is 'P\\tchr1\\t...')\n"
             "- 'GRCh38\\t0\\tchr1' (with actual tabs; the script will try to interpret this by looking for 'chr1', 'GRCh38', or 'GRCh38_0_chr1' as path names in GFA P records).\n"
             "  If your shell requires special syntax for tabs, use it (e.g., $'GRCh38\\t0\\tchr1' in bash)."
    )
    parser.add_argument(
        "--output_file", "-o",
        help="Optional. Path to save the JSON output. If not provided, prints to stdout.",
        required=False
    )

    if len(sys.argv) == 1:  # No arguments provided
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Call the main function
    json_result_str = extract_path_info_from_gfa(args.gfa_file, args.path_name)

    # Try to parse the JSON to check its status before printing/saving
    # This helps in providing better feedback if the function returned an error JSON
    try:
        result_data = json.loads(json_result_str)
    except json.JSONDecodeError:
        # Should not happen if extract_path_info_from_gfa always returns valid JSON string
        print(f"Error: The script produced invalid JSON output:\n{json_result_str}", file=sys.stderr)
        sys.exit(1)

    if args.output_file:
        try:
            with open(args.output_file, 'w') as f_out:
                f_out.write(json_result_str)
            print(f"Output successfully saved to {args.output_file}", file=sys.stderr)
            if result_data.get("status") != "success" and result_data.get("status") != "success_empty_path":
                print(f"Note: An issue was reported: {result_data.get('error') or result_data.get('message')}",
                      file=sys.stderr)
        except IOError as e:
            print(f"Error writing to output file '{args.output_file}': {e}", file=sys.stderr)
            print("\nJSON Result (stdout fallback):", file=sys.stderr)
            print(json_result_str)  # Print the JSON to stdout as a fallback
            sys.exit(1)
    else:
        # Print to stdout
        print(json_result_str)
        if result_data.get("status") != "success" and result_data.get("status") != "success_empty_path":
            # If printing to stdout, error messages are already part of the JSON
            # But a clear exit code can be helpful for scripting
            sys.exit(1)  # Exit with error code if not successful

    # Ensure a non-zero exit code if the process wasn't fully successful
    if result_data.get("status") != "success" and result_data.get("status") != "success_empty_path":
        sys.exit(1)