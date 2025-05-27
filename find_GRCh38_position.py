import json


def extract_path_info_from_gfa(gfa_file_path, target_path_name_input):
    """
    Extracts nodes from a specified path in a GFA file, calculates their
    cumulative positions, and outputs the information in JSON format.

    Args:
        gfa_file_path (str): The path to the GFA file.
        target_path_name_input (str): The name of the path to extract
                                     (e.g., "GRCh38\t0\tchr1" or "chr1").

    Returns:
        str: A JSON string containing the node information for the target path,
             or an error message if the path is not found or an issue occurs.
    """
    segments = {}  # To store segment sequences and lengths {segment_id: {'seq': 'ACGT', 'len': 4}}
    path_segments_oriented = []
    found_path = False

    # Determine the actual path name to search for in the GFA P records
    # This logic tries to infer the actual GFA PathName from the potentially composite user input.
    # If target_path_name_input contains tabs, it's likely a composite identifier.
    # We might try matching parts of it, or the user should provide the exact GFA PathName.

    # For this version, we will first try to match target_path_name_input directly.
    # If that fails and it contains tabs, we'll try common interpretations.

    potential_path_names_to_try = [target_path_name_input]
    if '\t' in target_path_name_input:
        parts = target_path_name_input.split('\t')
        if len(parts) > 0 and parts[-1] not in potential_path_names_to_try:  # e.g., "chr1"
            potential_path_names_to_try.append(parts[-1])
        if len(parts) > 0 and parts[0] not in potential_path_names_to_try:  # e.g., "GRCh38"
            potential_path_names_to_try.append(parts[0])
        # Consider a version with underscores if tabs were replaced
        underscored_name = target_path_name_input.replace('\t', '_')
        if underscored_name not in potential_path_names_to_try:
            potential_path_names_to_try.append(underscored_name)

    try:
        with open(gfa_file_path, 'r') as f_gfa:
            for line in f_gfa:
                parts = line.strip().split('\t')
                record_type = parts[0]

                if record_type == 'S':
                    segment_id = parts[1]
                    sequence = parts[2]
                    # Optional fields like LN:i:<length> can also be parsed if present
                    # For now, relying on len(sequence)
                    segments[segment_id] = {'seq': sequence, 'len': len(sequence)}
                elif record_type == 'P':
                    path_id_from_file = parts[1]  # GFA PathName is parts[1]

                    if path_id_from_file in potential_path_names_to_try:
                        # Store the successfully matched name to use in messages
                        matched_path_name_in_gfa = path_id_from_file
                        found_path = True
                        segment_names_with_orientations = parts[2].split(',')
                        for seg_orient in segment_names_with_orientations:
                            path_segments_oriented.append(
                                {'id': seg_orient[:-1], 'strand': seg_orient[-1]}
                            )
                        break  # Found the path, stop reading P lines for this path

            # If path was found, ensure all its segments were loaded (if P appeared before all S)
            # This simple parser assumes S records are generally before or sufficiently populated
            # when P record is found. A multi-pass parser could be more robust.

    except FileNotFoundError:
        return json.dumps({"error": f"File not found: {gfa_file_path}"})
    except Exception as e:
        return json.dumps({"error": f"Error reading GFA file: {str(e)}"})

    if not segments:
        return json.dumps({"error": "No segments (S records) found in the GFA file."})

    if not found_path:
        # Constructing the error message carefully to avoid problematic f-string expressions
        error_message_main = f"Path matching '{target_path_name_input}' not found in the GFA file."
        tried_names_str = ", ".join(
            [f"'{name}'" for name in potential_path_names_to_try if name != target_path_name_input])
        if tried_names_str:
            error_message_detail = f" Also tried interpreting it as: {tried_names_str}."
            full_error_message = error_message_main + error_message_detail
        else:
            full_error_message = error_message_main
        return json.dumps({"error": full_error_message})

    if not path_segments_oriented:
        # Use matched_path_name_in_gfa if path was found but empty
        return json.dumps({"error": f"Path '{matched_path_name_in_gfa}' was found but contains no segments."})

    output_nodes = []
    current_cumulative_pos = 0

    for seg_info in path_segments_oriented:
        seg_id = seg_info['id']
        strand = seg_info['strand']

        if seg_id not in segments:
            return json.dumps(
                {"error": f"Segment '{seg_id}' from path '{matched_path_name_in_gfa}' not found in S records."})

        node_data = segments[seg_id]
        node_len = node_data['len']
        node_seq = node_data['seq']

        output_nodes.append({
            "node_id": seg_id,
            "grch38_position_start": current_cumulative_pos,
            "strand_in_path": strand,
            "sequence": node_seq,
            "length": node_len
        })
        current_cumulative_pos += node_len

    return json.dumps(output_nodes, indent=4)


if __name__ == '__main__':
    dummy_gfa_content = """H\tVN:Z:1.0
S\t1\tACGT\tLN:i:4
S\t2\tTGA\tLN:i:3
S\t3\tCCGGAA\tLN:i:6
S\t4\tN\tLN:i:1
S\tGRCh38_s1\tGATTACA
S\tGRCh38_s2\tTCAT
P\tGRCh38_path1\tGRCh38_s1+,GRCh38_s2-\t*
P\tchr1\t1+,2-,3+\t45M,3M,6M
P\tpathX\t1+,4-\t*
P\tGRCh38_0_chr1\t1+,2-\t*
"""
    file_path = "example.gfa"
    with open(file_path, "w") as f:
        f.write(dummy_gfa_content)

    # Scenario 1: User provides the exact path name as in GFA P record
    print("--- Scenario 1: Exact GFA Path Name 'chr1' ---")
    json_output_scenario1 = extract_path_info_from_gfa(file_path, "chr1")
    print("JSON Output:")
    print(json_output_scenario1)
    print("-" * 30)

    # Scenario 2: User provides a composite name "GRCh38\t0\tchr1"
    # The script will try "GRCh38\t0\tchr1", then "chr1", then "GRCh38", then "GRCh38_0_chr1"
    user_composite_path = "GRCh38\t0\tchr1"  # \t is a tab character
    # This will match "chr1" from the dummy GFA, or "GRCh38_0_chr1" if present.
    # In this dummy GFA, "chr1" and "GRCh38_0_chr1" both exist. "GRCh38_0_chr1" might be tried first if target_path_name_input is that.
    # The order in potential_path_names_to_try matters.

    print(f"--- Scenario 2: User composite path '{user_composite_path}' ---")  # Prints with tab
    # For a literal display showing \t:
    # print(f"--- Scenario 2: User composite path 'GRCh38\\t0\\tchr1' ---")
    json_output_scenario2 = extract_path_info_from_gfa(file_path, user_composite_path)
    print("JSON Output:")
    print(json_output_scenario2)  # Should find 'chr1' or 'GRCh38_0_chr1' based on current logic
    print("-" * 30)

    # Scenario 3: Path name with underscores that matches GFA
    print("--- Scenario 3: Path name 'GRCh38_0_chr1' ---")
    json_output_scenario3 = extract_path_info_from_gfa(file_path, "GRCh38_0_chr1")
    print("JSON Output:")
    print(json_output_scenario3)
    print("-" * 30)

    # Scenario 4: Path not found
    print("--- Scenario 4: Path not found 'nonexistent_path' ---")
    json_output_scenario4 = extract_path_info_from_gfa(file_path, "nonexistent_path")
    print("JSON Output:")
    print(json_output_scenario4)
    print("-" * 30)

    # Scenario 5: Path not found with composite name
    user_nonexistent_composite_path = "MyGenome\tX\tgene1"
    print(f"--- Scenario 5: Composite Path not found '{user_nonexistent_composite_path}' ---")
    json_output_scenario5 = extract_path_info_from_gfa(file_path, user_nonexistent_composite_path)
    print("JSON Output:")
    print(json_output_scenario5)
    print("-" * 30)