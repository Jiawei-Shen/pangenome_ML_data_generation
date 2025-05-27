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


def find_path_line_with_grep(gfa_file_path, w_parts_target, potential_p_names):
    """
    Uses grep to find a W or P line.
    Returns (line_content, "W" or "P", matched_identifier_parts_or_name)
    """
    # Try W-line first if w_parts_target is provided
    if w_parts_target:
        # W line fields: W(0) SampleID(1) HaplotypeID(2) SeqName(3) Start(4) End(5) Path(6)
        pattern = f"^W\t{re.escape(w_parts_target[0])}\t{re.escape(w_parts_target[1])}\t{re.escape(w_parts_target[2])}\t"
        cmd = ['grep', '-E', '-m', '1', pattern, gfa_file_path]
        logging.debug(f"Executing grep for W-line: {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False, encoding='utf-8')
            if result.returncode == 0 and result.stdout.strip():
                line_content = result.stdout.strip()
                logging.info(f"W-line found via grep: {line_content.split(chr(9))[0:4]}")  # Log key parts
                return line_content, "W", w_parts_target  # Return the original target tuple
            elif result.returncode == 1:
                logging.info("No W-line found via grep matching the 3-part key.")
            else:
                logging.warning(f"grep for W-line failed (RC={result.returncode}): {result.stderr.strip()}")
        except FileNotFoundError:
            logging.error("grep command not found. Please ensure grep is installed and in PATH.")
            return None, None, None  # Critical error, cannot proceed with grep
        except Exception as e:
            logging.error(f"Exception running grep for W-line: {e}")
            return None, None, None

    # Try P-line if W not found or not targeted
    logging.debug(f"Attempting P-line grep with names: {potential_p_names}")
    for p_name in potential_p_names:
        pattern = f"^P\t{re.escape(p_name)}\t"
        cmd = ['grep', '-E', '-m', '1', pattern, gfa_file_path]
        logging.debug(f"Executing grep for P-line ('{p_name}'): {' '.join(cmd)}")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False, encoding='utf-8')
            if result.returncode == 0 and result.stdout.strip():
                line_content = result.stdout.strip()
                logging.info(f"P-line found for '{p_name}' via grep: {line_content.split(chr(9))[0:2]}")
                return line_content, "P", p_name  # Return the matched P-name
            elif result.returncode == 1:
                logging.info(f"No P-line found via grep for '{p_name}'.")  # Expected if no match
            else:
                logging.warning(f"grep for P-line '{p_name}' failed (RC={result.returncode}): {result.stderr.strip()}")
        except FileNotFoundError:  # Should have been caught above, but good practice
            logging.error("grep command not found during P-line search.")
            return None, None, None
        except Exception as e:
            logging.error(f"Exception running grep for P-line '{p_name}': {e}")

    return None, None, None  # No path line found


def extract_path_info_from_gfa(gfa_file_path, target_path_name_input):
    # Determine target types
    w_parts_target = None
    if '\t' in target_path_name_input:
        parts = target_path_name_input.split('\t')
        if len(parts) == 3: w_parts_target = tuple(parts)

    potential_p_names_to_try = [target_path_name_input]
    if '\t' in target_path_name_input:
        parts = target_path_name_input.split('\t')
        if parts[-1] not in potential_p_names_to_try: potential_p_names_to_try.append(parts[-1])
        if parts[0] not in potential_p_names_to_try: potential_p_names_to_try.append(parts[0])
        underscored = target_path_name_input.replace('\t', '_')
        if underscored not in potential_p_names_to_try: potential_p_names_to_try.append(underscored)
    unique_potential_p_names = []
    for name in potential_p_names_to_try:
        if name not in unique_potential_p_names: unique_potential_p_names.append(name)

    # Step 1: Find the specific W or P line using grep
    found_line_content, path_source_type, matched_identifier = find_path_line_with_grep(
        gfa_file_path, w_parts_target, unique_potential_p_names
    )

    if not found_line_content:
        err_msg = f"Path matching '{target_path_name_input}' not found using grep for W or P lines."
        logging.warning(err_msg)
        return json.dumps({"error": err_msg, "status": "error_path_not_found_grep"})

    # Step 2: Parse the found path line
    final_path_segments_oriented = []
    path_identifier_gfa = ""
    line_parts = found_line_content.split('\t')

    if path_source_type == "W":
        if len(line_parts) < 7:  # W Samp Hap SeqName Start End Path
            err_msg = f"Malformed W-line returned by grep: '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_w_line"})
        w_path_str = line_parts[6]
        final_path_segments_oriented = parse_w_line_path(w_path_str)
        # matched_identifier is (Sample, Haplotype, SeqName)
        path_identifier_gfa = f"W:{matched_identifier[0]}/{matched_identifier[1]}/{matched_identifier[2]}"
        logging.info(
            f"Successfully extracted {len(final_path_segments_oriented)} segment refs from W-line path string.")
    elif path_source_type == "P":
        if len(line_parts) < 3:  # P PathName Segments CIGARs(opt)
            err_msg = f"Malformed P-line returned by grep: '{found_line_content}'"
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_malformed_grep_p_line"})
        p_segments_str = line_parts[2]
        final_path_segments_oriented = parse_p_line_path(p_segments_str)
        path_identifier_gfa = f"P:{matched_identifier}"  # matched_identifier is the P-name
        logging.info(
            f"Successfully extracted {len(final_path_segments_oriented)} segment refs from P-line segment string.")

    if not final_path_segments_oriented:
        msg = f"Path '{path_identifier_gfa}' (type {path_source_type}) was defined but its path string yielded no segments."
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
                    continue
                parts = clean_line.split('\t')
                if len(parts) < 3: continue
                segments[parts[1]] = {'seq': parts[2], 'len': len(parts[2])}
        logging.info(f"Loaded {len(segments)} S-records.")
        if not segments:
            logging.warning("No S-records loaded. Path processing will likely fail if segments are needed.")

    except FileNotFoundError:  # Should have been caught by grep call if file truly missing
        logging.error(f"File not found during S-record loading: {gfa_file_path}")
        return json.dumps({"error": f"File not found: {gfa_file_path}", "status": "error"})
    except Exception as e:
        logging.error(f"Error reading GFA for S-records: {str(e)}")
        return json.dumps({"error": f"Error reading GFA for S-records: {str(e)}", "status": "error"})

    # Step 4: Process nodes
    output_nodes = []
    current_cumulative_pos = 0
    for seg_info in final_path_segments_oriented:
        seg_id, strand = seg_info['id'], seg_info['strand']
        if seg_id not in segments:
            err_msg = f"Segment '{seg_id}' from path '{path_identifier_gfa}' not found in S records."
            logging.error(err_msg)
            return json.dumps({"error": err_msg, "status": "error_segment_not_found"})
        node_data = segments[seg_id]
        output_nodes.append({
            "node_id": seg_id, "grch38_position_start": current_cumulative_pos,
            "strand_in_path": strand, "sequence": node_data['seq'], "length": node_data['len']
        })
        current_cumulative_pos += node_data['len']

    logging.info(f"Successfully processed {len(output_nodes)} nodes for path '{path_identifier_gfa}'.")
    return json.dumps({
        "path_name_input": target_path_name_input, "path_identifier_gfa": path_identifier_gfa,
        "path_source_type": path_source_type, "nodes": output_nodes, "status": "success"
    }, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Uses grep to find a W or P path in a GFA, then extracts node info.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("gfa_file", help="Path to the GFA file.")
    parser.add_argument("path_name", help="Target path name (e.g., 'chr1' for P, or 'Sample\\tHap\\tSeq' for W).")
    parser.add_argument("--output_file", "-o", help="Optional output JSON file path.")
    parser.add_argument("--debug", action="store_true", help="Enable debug level logging.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr);
        sys.exit(1)
    args = parser.parse_args()

    if args.debug: logging.getLogger().setLevel(logging.DEBUG)

    # Check for grep availability
    try:
        subprocess.run(['grep', '--version'], capture_output=True, check=True, text=True, encoding='utf-8')
        logging.debug("grep command is available.")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        logging.critical(f"grep command not found or not functional. This script relies on grep. Error: {e}")
        sys.exit(1)

    logging.info(f"Processing GFA: '{args.gfa_file}', Target Path: '{args.path_name}'")
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
            print(json_result_str, file=sys.stdout)  # Fallback
            sys.exit(1)
    else:
        print(json_result_str, file=sys.stdout)

    status = result_data.get("status", "unknown_status")
    if status not in ["success", "success_empty_path_definition"]:
        logging.warning(
            f"Process ended with status '{status}'. Error/Message: {result_data.get('error') or result_data.get('message')}")
        sys.exit(1)
    else:
        logging.info(f"Process completed with status: {status}.")