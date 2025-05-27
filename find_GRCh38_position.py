import json
import argparse
import sys
import re
import logging
import subprocess

# Configure logging (can be overridden by --debug)
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
        if len(seg_orient) < 2:  # segment name + orientation
            logging.warning(f"Skipping malformed P-line segment entry '{seg_orient}'.")
            continue
        seg_id = seg_orient[:-1]
        strand = seg_orient[-1]
        if strand not in ['+', '-']:
            logging.warning(f"Skipping P-line segment entry '{seg_orient}' with invalid strand '{strand}'.")
            continue
        oriented_segments.append({'id': seg_id, 'strand': strand})
    return oriented_segments


def find_path_line_with_grep(gfa_file_path, w_parts_target, potential_p_names):
    """
    Uses grep -P to find a W or P line.
    Returns (line_content, "W" or "P", matched_identifier_parts_or_name)
    """
    # Try W-line first if w_parts_target is provided
    if w_parts_target:
        # W line fields: W(0) SampleID(1) HaplotypeID(2) SeqName(3) Start(4) End(5) Path(6)
        # Using re.escape for robustness, though for simple names it won't change them.
        pattern = f"^W\t{re.escape(w_parts_target[0])}\t{re.escape(w_parts_target[1])}\t{re.escape(w_parts_target[2])}\t"
        # Using -P for Perl-compatible regex to ensure \t is handled as tab
        cmd = ['grep', '-P', '-m', '1', pattern, gfa_file_path]
        logging.debug(f"Executing grep for W-line: {' '.join(cmd)}")
        try:
            # Using utf-8 encoding for output from grep
            result = subprocess.run(cmd, capture_output=True, text=True, check=False, encoding='utf-8')
            if result.returncode == 0 and result.stdout.strip():
                line_content = result.stdout.strip()
                # Log key identifying parts of the W-line found
                log_w_key = "\t".join(line_content.split('\t')[1:4])
                logging.info(f"W-line found via grep for key '{log_w_key}'.")
                return line_content, "W", w_parts_target  # Return the original target tuple
            elif result.returncode == 1:  # grep found nothing (normal case if no match)
                logging.info("No W-line found via grep matching the 3-part key.")
            else:  # grep command itself had an error
                logging.warning(f"grep command for W-line failed (RC={result.returncode}): {result.stderr.strip()}")
        except FileNotFoundError:  # grep command not found
            logging.error("grep command not found. Please ensure grep is installed and in PATH.")
            return None, None, None  # Critical error, cannot proceed with grep
        except Exception as e:
            logging.error(f"Exception running grep for W-line: {e}")
            return None, None, None

    # Try P-line if W not found or not targeted
    logging.debug(f"Attempting P-line grep with potential names: {potential_p_names}")
    for p_name in potential_p_names:
        pattern = f"^P\t{re.escape(p_name)}\t"
        cmd = ['grep', '-P', '-m', '1', pattern, gfa_file_path]
        logging.debug(f"Executing grep for P-line ('{p_name}'): {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False, encoding='utf-8')
            if result.returncode == 0 and result.stdout.strip():
                line_content = result.stdout.strip()
                logging.info(f"P-line found for name '{p_name}' via grep.")
                return line_content, "P", p_name  # Return the matched P-name
            elif result.returncode == 1:  # grep found nothing
                logging.info(f"No P-line found via grep for name '{p_name}'.")  # Expected if no match
            else:  # grep command error
                logging.warning(
                    f"grep command for P-line '{p_name}' failed (RC={result.returncode}): {result.stderr.strip()}")
        except FileNotFoundError:  # Should have been caught by W-line grep if grep is missing
            logging.error("grep command not found during P-line search. Please ensure grep is installed and in PATH.")
            return None, None, None  # Critical error
        except Exception as e:
            logging.error(f"Exception running grep for P-line '{p_name}': {e}")
            # Continue to try other P-names if one grep call fails for other reasons than grep missing

    return None, None, None  # No path line found by any grep attempt


def extract_path_info_from_gfa(gfa_file_path, target_path_name_input):
    # Determine target types for W and P lines
    w_parts_target = None
    if '\t' in target_path_name_input:
        parts = target_path_name_input.split('\t')
        if len(parts) == 3:
            w_parts_target = tuple(parts)
            logging.info(
                f"Input path '{target_path_name_input}' interpreted as potential W-line target: Sample='{parts[0]}', Haplotype='{parts[1]}', SeqName='{parts[2]}'")

    potential_p_names_to_try = [target_path_name_input]  # Always try the raw input for P-lines too
    if '\t' in target_path_name_input:  # If it was a composite name, derive other P-line candidates
        parts = target_path_name_input.split('\t')
        if parts[-1] not in potential_p_names_to_try: potential_p_names_to_try.append(parts[-1])  # e.g., chr1
        if parts[0] not in potential_p_names_to_try: potential_p_names_to_try.append(parts[0])  # e.g., GRCh38
        underscored_name = target_path_name_input.replace('\t', '_')
        if underscored_name not in potential_p_names_to_try: potential_p_names_to_try.append(underscored_name)

    unique_potential_p_names = []
    for name in potential_p_names_to_try:
        if name not in unique_potential_p_names:
            unique_potential_p_names.append(name)

    # Step 1: Find the specific W or P line using grep
    found_line_content, path_source_type, matched_identifier = find_path_line_with_grep(
        gfa_file_path, w_parts_target, unique_potential_p_names
    )

    if not found_line_content:
        # If grep itself was not found, find_path_line_with_grep would have logged an error
        # and returned None, None, None. This indicates no path matched.
        err_msg = f"Path matching '{target_path_name_input}' not found using grep for W or P lines."
        logging.warning(err_msg)
        return json.dumps({"error": err_msg, "status": "error_path_not_found_grep"})

    # Step 2: Parse the found path line (obtained from grep)
    final_path_segments_oriented = []
    path_identifier_gfa = ""  # This will be the string identifier like "W:S/H/N" or "P:Name"

    line_parts = found_line_content.split('\t')  # GFA is tab-delimited

    if path_source_type == "W":
        # W line fields: W(0) SampleID(1) HaplotypeID(2) SeqName(3) Start(4) End(5) Path(6)
        if len(line_parts) < 7:  # Minimum fields for a W line with a path
            err_msg = f"Malformed W-line structure returned by grep (expected at least 7 fields): '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_w_line"})
        w_path_str = line_parts[6]
        final_path_segments_oriented = parse_w_line_path(w_path_str)
        # matched_identifier is the tuple (Sample, Haplotype, SeqName) used for searching
        path_identifier_gfa = f"W:{matched_identifier[0]}/{matched_identifier[1]}/{matched_identifier[2]}"
        logging.info(
            f"Successfully parsed W-line '{path_identifier_gfa}'; found {len(final_path_segments_oriented)} segment references from its path string.")
    elif path_source_type == "P":
        # P line fields: P(0) PathName(1) SegmentNamesOrientations(2) CIGARs(3 - optional)
        if len(line_parts) < 3:  # Minimum fields for a P line with segments
            err_msg = f"Malformed P-line structure returned by grep (expected at least 3 fields): '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_p_line"})
        p_segments_str = line_parts[2]
        final_path_segments_oriented = parse_p_line_path(p_segments_str)
        path_identifier_gfa = f"P:{matched_identifier}"  # matched_identifier is the P-name string
        logging.info(
            f"Successfully parsed P-line '{path_identifier_gfa}'; found {len(final_path_segments_oriented)} segment references from its segment string.")

    if not final_path_segments_oriented and path_source_type:  # Path line found, but its path definition was empty
        msg = f"Path '{path_identifier_gfa}' (type {path_source_type}) was found and parsed, but its path definition string yielded no segments."
        logging.warning(msg)
        return json.dumps({
            "message": msg, "path_name_input": target_path_name_input,
            "path_identifier_gfa": path_identifier_gfa, "path_source_type": path_source_type,
            "nodes": [], "status": "success_empty_path_definition"
        })

    # Step 3: Load S-records (still requires reading the GFA file)
    segments = {}
    logging.info(f"Loading S-records from GFA file: {gfa_file_path}")
    try:
        with open(gfa_file_path, 'r') as f_gfa:
            for line_content in f_gfa:
                clean_line = line_content.strip()
                if not clean_line or clean_line.startswith('#') or not clean_line.startswith('S\t'):
                    continue  # Skip comments, empty lines, or non-S lines quickly
                parts = clean_line.split('\t')
                if len(parts) < 3: continue  # Malformed S-line
                segments[parts[1]] = {'seq': parts[2], 'len': len(parts[2])}
        logging.info(f"Loaded {len(segments)} S-records.")
        if not segments and final_path_segments_oriented:  # Path needs segments but none were loaded
            logging.warning(
                "No S-records loaded from GFA file, but path segments were identified. This will likely lead to errors.")

    except FileNotFoundError:  # Should have been caught by grep call if file truly missing
        logging.error(f"File not found during S-record loading: {gfa_file_path}")  # Should be rare here
        return json.dumps({"error": f"File not found: {gfa_file_path}", "status": "error_file_not_found_s_load"})
    except Exception as e:
        logging.error(f"Error reading GFA file for S-records: {str(e)}")
        return json.dumps(
            {"error": f"Error reading GFA file for S-records: {str(e)}", "status": "error_s_load_exception"})

    # Step 4: Process nodes using loaded S-records and parsed path segments
    output_nodes = []
    current_cumulative_pos = 0
    for seg_info in final_path_segments_oriented:
        seg_id = seg_info['id']
        strand = seg_info['strand']
        if seg_id not in segments:
            err_msg = f"Segment '{seg_id}' from path '{path_identifier_gfa}' (type {path_source_type}) was defined in the path, but not found among loaded S records."
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_segment_not_found_in_s_records"})

        node_data = segments[seg_id]
        output_nodes.append({
            "node_id": seg_id,
            "grch38_position_start": current_cumulative_pos,
            # Name 'grch38_position_start' is kept as per original request
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
        description="Uses grep to find a W or P path in a GFA file, then extracts node information and calculates cumulative positions.\n"
                    "W-lines are prioritized if the input path name matches the W-line key format (Sample\\tHaplotype\\tSeqName).\n"
                    "The 'grch38_position_start' in the output is generic for the cumulative position along the specified path.",
        formatter_class=argparse.RawTextHelpFormatter  # Allows for newlines in help text
    )
    parser.add_argument(
        "gfa_file",
        help="Path to the GFA file."
    )
    parser.add_argument(
        "path_name",
        help="Target path name to extract. Examples:\n"
             "- For P-line: 'chr1' (if GFA P record is 'P\\tchr1\\t...')\n"
             "- For W-line: 'SampleA\\t0\\tchrX' (use actual tabs; e.g., $'SampleA\\t0\\tchrX' in bash on Linux/macOS).\n"
             "  If a 3-part tab-separated name is given, W-line matching via grep is attempted first.\n"
             "  Otherwise, or if W-line match fails, P-line matching via grep is attempted based on this name or its derived parts."
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

    if len(sys.argv) == 1:  # No arguments provided
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)  # Set root logger level
        logging.debug("Debug logging enabled.")

    # Early check for grep availability
    try:
        # Corrected line: removed stderr=subprocess.DEVNULL as it conflicts with capture_output=True
        subprocess.run(['grep', '--version'], capture_output=True, check=True, text=True, encoding='utf-8')
        logging.debug("grep command is available and functional.")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.critical(
            f"Critical: grep command not found or not functional. This script relies on grep. Please ensure 'grep' is installed and in your system's PATH. Error: {e}")
        sys.exit(1)  # Exit if grep is not available

    logging.info(f"Starting GFA path extraction. File: '{args.gfa_file}', Target Path: '{args.path_name}'")
    json_result_str = extract_path_info_from_gfa(args.gfa_file, args.path_name)

    # Try to parse the JSON output to check its status before printing/saving
    # This helps in providing better feedback if the function returned an error JSON
    try:
        result_data = json.loads(json_result_str)
    except json.JSONDecodeError:
        # This should ideally not happen if extract_path_info_from_gfa always returns a valid JSON string
        logging.critical(
            f"The script internally produced invalid JSON output. This is a bug. Output was:\n{json_result_str}")
        print(json_result_str, file=sys.stdout)  # Print it anyway for inspection
        sys.exit(1)

    if args.output_file:
        try:
            with open(args.output_file, 'w') as f_out:
                f_out.write(json_result_str)
            logging.info(f"Output successfully saved to {args.output_file}")
        except IOError as e:
            logging.error(f"Error writing to output file '{args.output_file}': {e}")
            # If writing to file fails, also print the JSON to stdout as a fallback
            print("\nJSON Result (stdout fallback due to file write error):", file=sys.stderr)
            print(json_result_str, file=sys.stdout)
            sys.exit(1)  # Exit with error code
    else:
        # Print to stdout by default
        print(json_result_str, file=sys.stdout)

    # Provide a clear overall status message to stderr and set exit code
    status = result_data.get("status", "unknown_status")  # Default to unknown if status field is missing
    if status not in ["success", "success_empty_path_definition"]:
        # Error messages or warnings should have been logged by extract_path_info_from_gfa
        logging.warning(
            f"Process completed with status: '{status}'. Check logs for details. Reported error/message: {result_data.get('error') or result_data.get('message', 'N/A')}")
        sys.exit(1)  # Exit with error code if not fully successful
    else:
        logging.info(f"Process completed successfully with status: '{status}'.")