import json
import argparse
import sys
import re
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', stream=sys.stderr)


def parse_w_line_path(w_path_str):
    """
    Parses the <Path> string from a GFA W line (e.g., ">s1<s2>s3")
    into a list of oriented segments.
    Example: ">s1<s2" -> [{'id': 's1', 'strand': '+'}, {'id': 's2', 'strand': '-'}]
    """
    oriented_segments = []
    # Regex to find all occurrences of an orientation character followed by a segment ID.
    # Segment ID allows alphanumeric, underscore, dot, hyphen.
    matches = re.findall(r"([<>])([\w.-]+)", w_path_str)
    for orientation_char, seg_id in matches:
        strand = '+' if orientation_char == '>' else '-'
        oriented_segments.append({'id': seg_id, 'strand': strand})
    if not matches and w_path_str:  # Check if regex failed but path string was not empty
        logging.warning(
            f"Could not parse any segments from W-line path string: '{w_path_str}'. It might have an unexpected format.")
    return oriented_segments


def parse_p_line_path(p_segments_str):
    """
    Parses the SegmentNamesOrientations string from a GFA P line (e.g., "s1+,s2-,s3+")
    into a list of oriented segments.
    """
    oriented_segments = []
    if not p_segments_str:  # Handle empty segment list string
        return oriented_segments

    segment_entries = p_segments_str.split(',')
    for seg_orient in segment_entries:
        if len(seg_orient) < 2:
            logging.warning(f"Skipping malformed P-line segment entry '{seg_orient}'.")
            continue
        seg_id = seg_orient[:-1]
        strand = seg_orient[-1]
        if strand not in ['+', '-']:
            logging.warning(f"Skipping P-line segment entry '{seg_orient}' with invalid strand '{strand}'.")
            continue
        oriented_segments.append({'id': seg_id, 'strand': strand})
    return oriented_segments


def extract_path_info_from_gfa(gfa_file_path, target_path_name_input):
    segments = {}

    # For W line matching
    w_parts_target = None
    if '\t' in target_path_name_input:
        parts = target_path_name_input.split('\t')
        if len(parts) == 3:
            w_parts_target = tuple(parts)  # (SampleID, HaplotypeID, SeqName)
            logging.info(
                f"Input path '{target_path_name_input}' interpreted as potential W-line target: Sample='{parts[0]}', Haplotype='{parts[1]}', SeqName='{parts[2]}'")

    # For P line matching
    potential_p_names_to_try = [target_path_name_input]
    if '\t' in target_path_name_input:  # If it was a composite name
        parts = target_path_name_input.split('\t')
        if parts[-1] not in potential_p_names_to_try: potential_p_names_to_try.append(parts[-1])  # e.g., chr1
        if parts[0] not in potential_p_names_to_try: potential_p_names_to_try.append(parts[0])  # e.g., GRCh38
        underscored_name = target_path_name_input.replace('\t', '_')
        if underscored_name not in potential_p_names_to_try: potential_p_names_to_try.append(underscored_name)

    unique_potential_p_names = []
    for name in potential_p_names_to_try:
        if name not in unique_potential_p_names:
            unique_potential_p_names.append(name)
    potential_p_names_to_try = unique_potential_p_names

    # Store found path segments and type
    final_path_segments_oriented = []
    path_source_type = None  # "W" or "P"
    path_identifier_gfa = None  # The actual identifier matched in GFA

    # --- Pass 1: Collect all S records and identify W/P paths ---
    # We need all S records before processing paths, so a multi-pass approach is more robust.
    # Pass 1.A: Collect all S, W, P lines
    s_records = {}
    w_lines_raw = []
    p_lines_raw = []

    logging.info(f"Starting GFA file parsing for S, W, P records from: {gfa_file_path}")
    try:
        with open(gfa_file_path, 'r') as f_gfa:
            for i, line_content in enumerate(f_gfa):
                line_num = i + 1
                clean_line = line_content.strip()
                if not clean_line or clean_line.startswith('#'):
                    continue
                parts = clean_line.split('\t')
                if not parts:
                    continue
                record_type = parts[0]

                if record_type == 'S':
                    if len(parts) < 3:
                        logging.warning(f"L{line_num}: Malformed S record: {clean_line}")
                        continue
                    s_records[parts[1]] = {'seq': parts[2], 'len': len(parts[2])}
                elif record_type == 'W':
                    if len(parts) < 6:  # W <SampleID> <HaplotypeID> <SeqName> <SeqStart> <SeqEnd> <Path>
                        logging.warning(f"L{line_num}: Malformed W record: {clean_line}")
                        continue
                    w_lines_raw.append({'fields': parts, 'line_num': line_num})
                elif record_type == 'P':
                    if len(parts) < 3:  # P <PathName> <SegmentNames> <Orientations>
                        logging.warning(f"L{line_num}: Malformed P record: {clean_line}")
                        continue
                    p_lines_raw.append({'fields': parts, 'line_num': line_num})
        segments = s_records
        logging.info(f"Collected {len(segments)} S-records, {len(w_lines_raw)} W-lines, {len(p_lines_raw)} P-lines.")

    except FileNotFoundError:
        logging.error(f"File not found: {gfa_file_path}")
        return json.dumps({"error": f"File not found: {gfa_file_path}", "status": "error"})
    except Exception as e:
        logging.error(f"Error reading GFA file '{gfa_file_path}': {str(e)}")
        return json.dumps({"error": f"Error reading GFA file '{gfa_file_path}': {str(e)}", "status": "error"})

    if not segments:
        logging.warning("No S (segment) records found in the GFA file.")
        # It's possible a path exists but its segments are not defined, which is an error later.

    # Pass 1.B: Try to find and parse W path first
    if w_parts_target:
        for w_line_info in w_lines_raw:
            w_fields = w_line_info['fields']
            # W <SampleID> <HaplotypeID> <SeqName> <SeqStart> <SeqEnd> <Path>
            # parts[0]   parts[1]    parts[2]    parts[3]    parts[4]   parts[5]
            current_w_key = (w_fields[1], w_fields[2], w_fields[3])
            if current_w_key == w_parts_target:
                logging.info(f"L{w_line_info['line_num']}: Matched W-line with target: {w_parts_target}")
                w_path_str = w_fields[6]  # <Path> is the 7th field (index 6)
                parsed_segments = parse_w_line_path(w_path_str)
                if parsed_segments:  # Successfully parsed some segments
                    final_path_segments_oriented = parsed_segments
                    path_source_type = "W"
                    path_identifier_gfa = f"W:{current_w_key[0]}/{current_w_key[1]}/{current_w_key[2]}"
                    logging.info(
                        f"Successfully extracted {len(parsed_segments)} segments from W-line '{path_identifier_gfa}'.")
                    break  # Use the first matching W line
                else:
                    logging.warning(
                        f"L{w_line_info['line_num']}: W-line '{path_identifier_gfa}' matched but its path string '{w_path_str}' yielded no segments.")

    # Pass 1.C: If no W path found, try P path
    if not final_path_segments_oriented:  # if path_source_type is None
        logging.info(
            f"No W-line match found or W-line path was empty. Attempting to find P-line match for targets: {potential_p_names_to_try}")
        for p_line_info in p_lines_raw:
            p_fields = p_line_info['fields']
            p_name_in_file = p_fields[1]  # PathName field
            if p_name_in_file in potential_p_names_to_try:
                logging.info(f"L{p_line_info['line_num']}: Matched P-line with name: '{p_name_in_file}'")
                p_segments_str = p_fields[2]  # SegmentNames,Orientations
                parsed_segments = parse_p_line_path(p_segments_str)
                if parsed_segments:
                    final_path_segments_oriented = parsed_segments
                    path_source_type = "P"
                    path_identifier_gfa = f"P:{p_name_in_file}"
                    logging.info(
                        f"Successfully extracted {len(parsed_segments)} segments from P-line '{path_identifier_gfa}'.")
                    break  # Use the first matching P line
                else:
                    logging.warning(
                        f"L{p_line_info['line_num']}: P-line '{p_name_in_file}' matched but its segment string '{p_segments_str}' yielded no segments.")

    # --- Pass 2: Process the found path ---
    if not path_source_type:
        error_message_main = f"Path matching '{target_path_name_input}' not found as a W or P line."
        tried_additionally = [f"'{name}'" for name in potential_p_names_to_try if
                              name != target_path_name_input and name not in (w_parts_target if w_parts_target else [])]
        if tried_additionally:
            error_message_detail = f" For P-lines, also tried interpreting it as: {', '.join(tried_additionally)}."
            full_error_message = error_message_main + error_message_detail
        else:
            full_error_message = error_message_main
        logging.warning(full_error_message)
        return json.dumps({"error": full_error_message, "status": "error"})

    if not final_path_segments_oriented:  # Path source was identified but it had no segments
        msg = f"Path '{path_identifier_gfa}' (type {path_source_type}) was found but contains no actual segments in its definition."
        logging.warning(msg)
        return json.dumps({
            "message": msg,
            "path_name_input": target_path_name_input,
            "path_identifier_gfa": path_identifier_gfa,
            "path_source_type": path_source_type,
            "nodes": [],
            "status": "success_empty_path"
        })

    output_nodes = []
    current_cumulative_pos = 0
    for seg_info in final_path_segments_oriented:
        seg_id = seg_info['id']
        strand = seg_info['strand']

        if seg_id not in segments:
            err_msg = f"Segment '{seg_id}' from path '{path_identifier_gfa}' (type {path_source_type}) not found in S records."
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error"})

        node_data = segments[seg_id]
        output_nodes.append({
            "node_id": seg_id,
            "grch38_position_start": current_cumulative_pos,
            "strand_in_path": strand,
            "sequence": node_data['seq'],
            "length": node_data['len']
        })
        current_cumulative_pos += node_data['len']

    logging.info(f"Successfully processed {len(output_nodes)} nodes for path '{path_identifier_gfa}'.")
    return json.dumps({
        "path_name_input": target_path_name_input,
        "path_identifier_gfa": path_identifier_gfa,
        "path_source_type": path_source_type,
        "nodes": output_nodes,
        "status": "success"
    }, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Extracts node information from a specified path (W or P line) in a GFA file and outputs JSON. \n"
                    "W-lines are prioritized if the input path name matches the W-line key format (Sample\\tHaplotype\\tSeqName).\n"
                    "The 'grch38_position_start' in the output is a cumulative position based on node lengths in the specified path.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "gfa_file",
        help="Path to the GFA file."
    )
    parser.add_argument(
        "path_name",
        help="Target path name to extract. Examples:\n"
             "- For P-line: 'chr1' (if GFA P record is 'P\\tchr1\\t...')\n"
             "- For W-line: 'SampleA\\t0\\tchrX' (use actual tabs; e.g., $'SampleA\\t0\\tchrX' in bash)\n"
             "  If a 3-part tab-separated name is given, W-line matching is attempted first.\n"
             "  Otherwise, or if W-line fails, P-line matching is attempted based on this name or its parts."
    )
    parser.add_argument(
        "--output_file", "-o",
        help="Optional. Path to save the JSON output. If not provided, prints to stdout.",
        required=False
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug level logging."
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)  # Set root logger level
        logging.debug("Debug logging enabled.")

    logging.info(f"Processing GFA file: '{args.gfa_file}' for path: '{args.path_name}'")
    json_result_str = extract_path_info_from_gfa(args.gfa_file, args.path_name)

    try:
        result_data = json.loads(json_result_str)
    except json.JSONDecodeError:
        logging.critical(f"The script produced invalid JSON output:\n{json_result_str}")
        print(json_result_str, file=sys.stdout)  # print it anyway for inspection
        sys.exit(1)

    if args.output_file:
        try:
            with open(args.output_file, 'w') as f_out:
                f_out.write(json_result_str)
            logging.info(f"Output successfully saved to {args.output_file}")
        except IOError as e:
            logging.error(f"Error writing to output file '{args.output_file}': {e}")
            print(json_result_str, file=sys.stdout)  # Print to stdout as a fallback
            sys.exit(1)  # Exit with error code
    else:
        print(json_result_str, file=sys.stdout)

    if result_data.get("status") != "success" and result_data.get("status") != "success_empty_path":
        logging.warning(
            f"Process completed with status: {result_data.get('status')}. Error: {result_data.get('error') or result_data.get('message')}")
        sys.exit(1)
    else:
        logging.info(f"Process completed successfully with status: {result_data.get('status')}.")