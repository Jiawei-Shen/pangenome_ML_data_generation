import json
import argparse
import sys
import re
import logging
import subprocess

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', stream=sys.stderr)


def parse_w_line_path(w_path_str):
    oriented_segments = []
    matches = re.findall(r"([<>])([\w.-]+)", w_path_str)
    for orientation_char, seg_id in matches:
        strand = '+' if orientation_char == '>' else '-'
        oriented_segments.append({'id': seg_id, 'strand': strand})
    if not matches and w_path_str:
        logging.warning(f"Could not parse segments from W-line path string: '{w_path_str}'.")
    return oriented_segments


def parse_p_line_path(p_segments_str):
    oriented_segments = []
    if not p_segments_str: return oriented_segments
    segment_entries = p_segments_str.split(',')
    for seg_orient in segment_entries:
        if len(seg_orient) < 2:
            logging.warning(f"Skipping malformed P-line segment entry '{seg_orient}'.")
            continue
        seg_id, strand = seg_orient[:-1], seg_orient[-1]
        if strand not in ['+', '-']:
            logging.warning(f"Skipping P-line segment entry '{seg_orient}' with invalid strand '{strand}'.")
            continue
        oriented_segments.append({'id': seg_id, 'strand': strand})
    return oriented_segments


def find_path_line_with_grep_direct(gfa_file_path, user_key_input_raw):
    """
    Uses grep -P to find a line based on a user-provided key.
    The key can optionally start with '^', which will be stripped.
    It must then start with 'W\t' or 'P\t'.
    Returns (line_content, "W" or "P", identifier_for_json_output)
    """
    processed_key = user_key_input_raw
    if user_key_input_raw.startswith("^"):
        processed_key = user_key_input_raw[1:]
        logging.debug(
            f"User input '{user_key_input_raw}' started with '^'; stripped to '{processed_key}' for internal processing.")
    else:
        logging.debug(f"User input '{user_key_input_raw}' processed as is (no leading '^' found).")

    path_type = None
    grep_pattern_core = ""
    identifier_for_json_output = None

    # Enhanced debugging for string checks
    is_w_type = processed_key.startswith("W\t")
    is_p_type = processed_key.startswith("P\t")
    logging.debug(
        f"For processed_key='{processed_key}' (length {len(processed_key)}): is_w_type={is_w_type}, is_p_type={is_p_type}")
    # For further debugging, show byte representation if tabs are tricky
    # logging.debug(f"Processed key bytes: {processed_key.encode('utf-8')}")

    if is_w_type:
        parts = processed_key.split('\t', 3)
        logging.debug(
            f"W-type detected for '{processed_key}'. Parts after split('\\t', 3): {parts} (length: {len(parts)})")
        if len(parts) == 4:  # W, Sample, Haplotype, SeqName
            path_type = "W"
            sample_id, haplotype_id, seq_name = parts[1], parts[2], parts[3]
            grep_pattern_core = f"W\t{re.escape(sample_id)}\t{re.escape(haplotype_id)}\t{re.escape(seq_name)}"
            identifier_for_json_output = (sample_id, haplotype_id, seq_name)
            logging.info(
                f"Interpreted processed key '{processed_key}' as W-line key: Sample='{sample_id}', Haplotype='{haplotype_id}', SeqName='{seq_name}'")
        else:
            logging.error(
                f"Processed W-line key '{processed_key}' (from input '{user_key_input_raw}') not in expected format 'W\\tSampleID\\tHaplotypeID\\tSeqName'. Expected 4 parts after splitting on tab (W, S, H, N), got {len(parts)} parts: {parts}.")
            return None, None, None
    elif is_p_type:
        parts = processed_key.split('\t', 1)
        logging.debug(
            f"P-type detected for '{processed_key}'. Parts after split('\\t', 1): {parts} (length: {len(parts)})")
        if len(parts) == 2:  # P, PathName
            path_type = "P"
            path_name = parts[1]
            grep_pattern_core = f"P\t{re.escape(path_name)}"
            identifier_for_json_output = path_name
            logging.info(f"Interpreted processed key '{processed_key}' as P-line key: PathName='{path_name}'")
        else:
            logging.error(
                f"Processed P-line key '{processed_key}' (from input '{user_key_input_raw}') not in expected format 'P\\tPathName'. Expected 2 parts after splitting on tab (P, Name), got {len(parts)} parts: {parts}.")
            return None, None, None
    else:
        logging.error(
            f"Invalid path_name input '{user_key_input_raw}'. After optional '^' stripping (resulting in '{processed_key}'), key must start with 'W\\t' or 'P\\t' (using actual tabs).")
        return None, None, None

    # Final grep pattern always anchors to start of line, and requires a tab after the core key.
    pattern = f"^{grep_pattern_core}\t"
    cmd = ['grep', '-P', '-m', '1', pattern, gfa_file_path]
    logging.debug(f"Executing grep: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False, encoding='utf-8')
        if result.returncode == 0 and result.stdout.strip():
            line_content = result.stdout.strip()
            logging.info(f"{path_type}-line found via grep using key derived from '{user_key_input_raw}'.")
            return line_content, path_type, identifier_for_json_output
        elif result.returncode == 1:
            logging.info(
                f"No {path_type}-line found via grep for key derived from '{user_key_input_raw}' (pattern: '{pattern}').")
            return None, None, None
        else:
            logging.warning(
                f"grep for {path_type}-line with key derived from '{user_key_input_raw}' failed (RC={result.returncode}): {result.stderr.strip()}")
            return None, None, None
    except FileNotFoundError:
        logging.error("grep command not found. Please ensure grep is installed and in PATH.")
        return None, None, None
    except Exception as e:
        logging.error(f"Exception running grep for key derived from '{user_key_input_raw}': {e}")
        return None, None, None


# extract_path_info_from_gfa and if __name__ == '__main__' remain the same as the previous full script.
# ... (rest of the script from the previous full version) ...
def extract_path_info_from_gfa(gfa_file_path, target_path_key_input):
    # Step 1: Find the specific W or P line using grep based on user's direct key
    found_line_content, path_source_type, matched_identifier_parts = find_path_line_with_grep_direct(
        gfa_file_path, target_path_key_input
    )

    if not found_line_content:
        err_msg = f"Path for key '{target_path_key_input}' not found using grep."
        # Specific error/reason should have been logged by find_path_line_with_grep_direct
        return json.dumps({"error": err_msg, "status": "error_path_not_found_grep"})

    # Step 2: Parse the found path line
    final_path_segments_oriented = []
    path_identifier_gfa = ""  # For JSON output, e.g., "W:Sample/Hap/Seq" or "P:Name"

    line_parts = found_line_content.split('\t')

    if path_source_type == "W":
        if len(line_parts) < 7:  # W(0) Samp(1) Hap(2) SeqN(3) Start(4) End(5) Path(6)
            err_msg = f"Malformed W-line structure returned by grep: '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_w_line"})
        w_path_str = line_parts[6]
        final_path_segments_oriented = parse_w_line_path(w_path_str)
        # matched_identifier_parts is (Sample, Haplotype, SeqName)
        path_identifier_gfa = f"W:{matched_identifier_parts[0]}/{matched_identifier_parts[1]}/{matched_identifier_parts[2]}"
        logging.info(
            f"Parsed W-line '{path_identifier_gfa}'; {len(final_path_segments_oriented)} segment refs from path string.")
    elif path_source_type == "P":
        if len(line_parts) < 3:  # P(0) PathName(1) Segments(2)
            err_msg = f"Malformed P-line structure returned by grep: '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_p_line"})
        p_segments_str = line_parts[2]
        final_path_segments_oriented = parse_p_line_path(p_segments_str)
        # matched_identifier_parts is the PathName string
        path_identifier_gfa = f"P:{matched_identifier_parts}"
        logging.info(
            f"Parsed P-line '{path_identifier_gfa}'; {len(final_path_segments_oriented)} segment refs from segment string.")

    if not final_path_segments_oriented and path_source_type:
        msg = f"Path '{path_identifier_gfa}' (type {path_source_type}) definition found, but its path string yielded no segments."
        logging.warning(msg)
        return json.dumps({
            "message": msg, "path_name_input": target_path_key_input,
            "path_identifier_gfa": path_identifier_gfa, "path_source_type": path_source_type,
            "nodes": [], "status": "success_empty_path_definition"
        })

    # Step 3: Load S-records
    segments = {}
    logging.info(f"Loading S-records from GFA file: {gfa_file_path}")
    try:
        with open(gfa_file_path, 'r') as f_gfa:
            for line_content in f_gfa:
                clean_line = line_content.strip()
                if not clean_line or clean_line.startswith('#') or not clean_line.startswith('S\t'):
                    continue
                parts = clean_line.split('\t')
                if len(parts) < 3: continue
                segments[parts[1]] = {'seq': parts[2], 'len': len(parts[2])}
        logging.info(f"Loaded {len(segments)} S-records.")
        if not segments and final_path_segments_oriented:
            logging.warning("No S-records loaded, but path needs segments. Errors likely.")
    except FileNotFoundError:
        logging.error(f"File not found during S-record loading: {gfa_file_path}")
        return json.dumps({"error": f"File not found: {gfa_file_path}", "status": "error_file_not_found_s_load"})
    except Exception as e:
        logging.error(f"Error reading GFA for S-records: {str(e)}")
        return json.dumps({"error": f"Error reading GFA for S-records: {str(e)}", "status": "error_s_load_exception"})

    # Step 4: Process nodes
    output_nodes = []
    current_cumulative_pos = 0
    for seg_info in final_path_segments_oriented:
        seg_id, strand = seg_info['id'], seg_info['strand']
        if seg_id not in segments:
            err_msg = f"Segment '{seg_id}' from path '{path_identifier_gfa}' not found in S records."
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_segment_not_found_in_s_records"})
        node_data = segments[seg_id]
        output_nodes.append({
            "node_id": seg_id, "grch38_position_start": current_cumulative_pos,
            "strand_in_path": strand, "sequence": node_data['seq'], "length": node_data['len']
        })
        current_cumulative_pos += node_data['len']

    logging.info(f"Successfully processed {len(output_nodes)} nodes for path '{path_identifier_gfa}'.")
    return json.dumps({
        "path_name_input": target_path_key_input, "path_identifier_gfa": path_identifier_gfa,
        "path_source_type": path_source_type, "nodes": output_nodes, "status": "success"
    }, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Uses grep to find a W or P path in a GFA file based on a user-provided key, then extracts node information.\n"
                    "The 'grch38_position_start' in the output is the cumulative position along the specified path.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "gfa_file",
        help="Path to the GFA file."
    )
    parser.add_argument(
        "path_name",
        help="Target path key. Should effectively start with 'W\\t' or 'P\\t' (using actual tabs from your shell).\n"
             "A single leading '^' character will be automatically stripped by the script if present.\n"
             "Examples:\n"
             "- For P-line: $'P\\tpathName1' or $'^P\\tpathName1' (Bash syntax for tab)\n"
             "- For W-line: $'W\\tSampleID\\tHaplotypeID\\tSeqName' or $'^W\\tSampleID\\tHaplotypeID\\tSeqName'\n"
             "The parts of the key (like PathName, SampleID, etc.) will be regex-escaped by the script."
    )
    parser.add_argument(
        "--output_file", "-o",
        help="Optional. Path to save the JSON output. If not provided, prints to stdout.",
        required=False
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug level logging (includes grep commands executed)."
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr);
        sys.exit(1)
    args = parser.parse_args()

    if args.debug: logging.getLogger().setLevel(logging.DEBUG)

    try:
        subprocess.run(['grep', '--version'], capture_output=True, check=True, text=True, encoding='utf-8')
        logging.debug("grep command is available and functional.")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.critical(f"Critical: grep command not found or not functional. Error: {e}")
        sys.exit(1)

    logging.info(
        f"Starting GFA path extraction. File: '{args.gfa_file}', Target Path Key (raw input): '{args.path_name}'")
    json_result_str = extract_path_info_from_gfa(args.gfa_file, args.path_name)

    try:
        result_data = json.loads(json_result_str)
    except json.JSONDecodeError:
        logging.critical(f"Script produced invalid JSON: {json_result_str}")
        print(json_result_str, file=sys.stdout)
        sys.exit(1)

    if args.output_file:
        try:
            with open(args.output_file, 'w') as f_out:
                f_out.write(json_result_str)
            logging.info(f"Output saved to {args.output_file}")
        except IOError as e:
            logging.error(f"Error writing to '{args.output_file}': {e}")
            print("\nJSON Result (stdout fallback):", file=sys.stderr)
            print(json_result_str, file=sys.stdout)
            sys.exit(1)
    else:
        print(json_result_str, file=sys.stdout)

    status = result_data.get("status", "unknown_status")
    if status not in ["success", "success_empty_path_definition"]:
        logging.warning(
            f"Process completed with status '{status}'. Error: {result_data.get('error') or result_data.get('message', 'N/A')}")
        sys.exit(1)
    else:
        logging.info(f"Process completed successfully with status: '{status}'.")