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


def find_path_line_with_user_pattern(gfa_file_path, user_grep_pattern):
    """
    Uses grep -P -m 1 with the exact pattern provided by the user.
    Returns the found line content or None.
    """
    cmd = ['grep', '-P', '-m', '1', user_grep_pattern, gfa_file_path]
    logging.debug(f"Executing user-defined grep: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False, encoding='utf-8')
        if result.returncode == 0 and result.stdout.strip():
            line_content = result.stdout.strip()
            logging.info(f"Line found via grep using user pattern '{user_grep_pattern}'.")
            return line_content
        elif result.returncode == 1:  # grep found nothing
            logging.info(f"No line found via grep for user pattern '{user_grep_pattern}'.")
            return None
        else:  # grep command itself had an error
            logging.warning(
                f"User-defined grep pattern '{user_grep_pattern}' failed (RC={result.returncode}): {result.stderr.strip()}")
            return None
    except FileNotFoundError:
        logging.error("grep command not found. Please ensure grep is installed and in PATH.")
        return None
    except Exception as e:
        logging.error(f"Exception running user-defined grep with pattern '{user_grep_pattern}': {e}")
        return None


def extract_path_info_from_gfa(gfa_file_path, user_grep_pattern_input):
    # Step 1: Find the specific line using grep with the user's exact pattern
    found_line_content = find_path_line_with_user_pattern(gfa_file_path, user_grep_pattern_input)

    if not found_line_content:
        err_msg = f"Path not found using user-provided grep pattern: '{user_grep_pattern_input}'."
        return json.dumps({"error": err_msg, "status": "error_path_not_found_user_grep"})

    # Step 2: Parse the found path line to determine its type and content
    final_path_segments_oriented = []
    path_identifier_gfa = ""
    path_source_type = None

    line_parts = found_line_content.split('\t')
    actual_record_type = line_parts[0] if line_parts else None

    logging.info(
        f"Line found by grep ('{user_grep_pattern_input}') starts with record type: '{actual_record_type}'. Parsing accordingly.")

    if actual_record_type == "W":
        path_source_type = "W"
        if len(line_parts) < 7:  # W(0) Samp(1) Hap(2) SeqN(3) Start(4) End(5) Path(6)
            err_msg = f"Malformed W-line structure returned by grep: '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_w_line"})

        sample_id, haplotype_id, seq_name = line_parts[1], line_parts[2], line_parts[3]
        path_identifier_gfa = f"W:{sample_id}/{haplotype_id}/{seq_name}"  # For context in output
        w_path_str = line_parts[6]
        final_path_segments_oriented = parse_w_line_path(w_path_str)
        logging.info(
            f"Identified found line as W-line: '{path_identifier_gfa}'; {len(final_path_segments_oriented)} segment refs from path string.")
    elif actual_record_type == "P":
        path_source_type = "P"
        if len(line_parts) < 3:  # P(0) PathName(1) Segments(2)
            err_msg = f"Malformed P-line structure returned by grep: '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_p_line"})

        path_name_from_line = line_parts[1]
        path_identifier_gfa = f"P:{path_name_from_line}"  # For context in output
        p_segments_str = line_parts[2]
        final_path_segments_oriented = parse_p_line_path(p_segments_str)
        logging.info(
            f"Identified found line as P-line: '{path_identifier_gfa}'; {len(final_path_segments_oriented)} segment refs from segment string.")
    else:
        err_msg = f"Line found by grep pattern '{user_grep_pattern_input}' was not a W or P line. Found: '{found_line_content}'"
        logging.error(err_msg)
        return json.dumps({"error": err_msg, "status": "error_unexpected_line_type_from_grep"})

    if not final_path_segments_oriented and path_source_type:
        msg = f"Path '{path_identifier_gfa}' (type {path_source_type}) definition found by grep, but its path string yielded no segments."
        logging.warning(msg)
        return json.dumps({
            "message": msg, "path_name_input_pattern": user_grep_pattern_input,  # Changed to reflect it's a pattern
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
        "path_name_input_pattern": user_grep_pattern_input,  # Changed field name
        "path_identifier_gfa": path_identifier_gfa,
        "path_source_type": path_source_type, "nodes": output_nodes, "status": "success"
    }, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Uses a user-provided grep pattern to find a W or P path in a GFA file, then extracts node information.\n"
                    "The script executes 'grep -P -m 1 <user_pattern> <gfa_file>'.\n"
                    "The 'grch38_position_start' in the output is the cumulative position along the specified path.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "gfa_file",
        help="Path to the GFA file."
    )
    parser.add_argument(
        "path_grep_pattern",
        help="The exact grep pattern (Perl-compatible, PCRE) to find the desired W or P line.\n"
             "The script will run 'grep -P -m 1 \"YOUR_PATTERN\" GFA_FILE'.\n"
             "Example to find a specific W-line (ensure your shell passes tabs correctly, e.g., using $'...'):\n"
             "  $'^W\\tGRCh38\\t0\\tchr1\\t'\n"
             "Example for a P-line:\n"
             "  $'^P\\tmyPathName\\t'"
    )
    parser.add_argument(
        "--output_file", "-o",
        help="Optional. Path to save the JSON output. If not provided, prints to stdout.",
        required=False
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug level logging (includes the full grep command executed)."
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
        f"Starting GFA path extraction. File: '{args.gfa_file}', User Grep Pattern: '{args.path_grep_pattern}'")
    json_result_str = extract_path_info_from_gfa(args.gfa_file, args.path_grep_pattern)

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