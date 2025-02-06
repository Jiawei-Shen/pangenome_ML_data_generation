import json
import argparse
import os
import subprocess

def extract_nodes(json_file_path):
    """Extracts node (segment) IDs for each genomic region from the given JSON file."""
    with open(json_file_path, "r") as f:
        data = json.load(f)

    node_dict = {}

    for region, details in data.items():
        if "segments" in details:
            node_ids = [segment["segment_id"] for segment in details["segments"]]
            node_dict[region] = node_ids

    return node_dict

def extract_reads_from_gam(gam_file, segment_ids):
    """Extracts reads from a GAM file that map to at least one of the given segment IDs."""
    result = subprocess.run(['vg', 'view', '-a', gam_file], stdout=subprocess.PIPE, text=True)

    if result.returncode != 0:
        print(f"Error running vg view on {gam_file}")
        return []

    matching_reads = set()

    for line in result.stdout.strip().split('\n'):
        read = json.loads(line)

        if 'path' in read:
            for mapping in read['path']['mapping']:
                if 'position' in mapping and 'node_id' in mapping['position']:
                    if str(mapping['position']['node_id']) in segment_ids:
                        matching_reads.add(read['name'])
                        break

    return matching_reads

def extract_reads_from_gaf(gaf_file, segment_ids):
    """Extracts reads from a GAF file that map to at least one of the given segment IDs."""
    matching_reads = set()

    with open(gaf_file, 'r') as f:
        for line in f:
            columns = line.strip().split("\t")
            if len(columns) > 5:
                path = columns[5]  # Segment/path information
                for node in segment_ids:
                    if node in path:
                        matching_reads.add(columns[0])  # Read name
                        break

    return matching_reads

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find reads mapped to nodes from a JSON file in a GAM or GAF file.")
    parser.add_argument("json_file", type=str, help="Input JSON file containing segment IDs.")
    parser.add_argument("alignment_file", type=str, help="Input GAM or GAF file for read extraction.")

    args = parser.parse_args()

    file_extension = os.path.splitext(args.alignment_file)[-1].lower()

    # Extract nodes from JSON
    nodes_dict = extract_nodes(args.json_file)
    all_segment_ids = {segment_id for segment_list in nodes_dict.values() for segment_id in segment_list}

    print(f"Extracted {len(all_segment_ids)} segment IDs from {args.json_file}")

    # Extract reads from GAM or GAF
    if file_extension == ".gam":
        reads = extract_reads_from_gam(args.alignment_file, all_segment_ids)
    elif file_extension == ".gaf":
        reads = extract_reads_from_gaf(args.alignment_file, all_segment_ids)
    else:
        raise ValueError("Unsupported file format. Please provide a GAM or GAF file.")

    # Print results
    if reads:
        print(f"\nReads mapped to extracted segment IDs ({len(reads)} total):")
        for read in reads:
            print(read)
    else:
        print("\nNo reads found for the extracted segment IDs.")
