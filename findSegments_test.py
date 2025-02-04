import json
import subprocess
import argparse
import sys

def extract_segments_from_json(gbz_file, path_range):
    # Run vg find and pipe to vg view to get JSON output
    cmd = f"vg find -x {gbz_file} -p {path_range} | vg view -v -j -"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    segment_info = {}

    for line in process.stdout:
        data = json.loads(line)

        # Extract node (segment) IDs
        if 'node' in data:
            for node in data['node']:
                node_id = node['id']

                # Check if this node is part of a path
                if 'path' in data:
                    for path in data['path']:
                        for mapping in path['mapping']:
                            if mapping['position']['node_id'] == node_id:
                                if node_id not in segment_info:
                                    segment_info[node_id] = {'segment_id': node_id, 'paths': []}
                                if path['name'] not in segment_info[node_id]['paths']:
                                    segment_info[node_id]['paths'].append(path['name'])

    process.stdout.close()
    process.wait()

    return list(segment_info.values())

def call_find_segments(gbz_file, chromosome, position, windowLength=100, ref="GRCh38#0"):
    """Call extract_segments_from_json with the given chromosome and position."""
    formatted_position = f"{ref}#{chromosome}:{position - windowLength}-{position + windowLength}"
    segments = extract_segments_from_json(gbz_file, formatted_position)
    print(f"Output for {formatted_position}: {segments}")
    return {formatted_position: segments}

def main(json_file, gbz_file):
    # Load JSON data
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)

    results = {}
    output_file = "segments_output.json"

    # Process each entry in the JSON
    for key, record in enumerate(data):
        chromosome = record.get('chromosome')
        position = record.get('position')

        if not chromosome or not position:
            print(f"Skipping record {key} due to missing chromosome or position.")
            continue

        print(f"\nProcessing {chromosome}:{position}")
        segment_data = call_find_segments(gbz_file, chromosome, position)
        results.update(segment_data)

        # Save results every 10 iterations
        if (key + 1) % 10 == 0:
            try:
                with open(output_file, 'w') as outfile:
                    json.dump(results, outfile, indent=4)
                print(f"Intermediate results saved to {output_file} after {key + 1} records.")
            except Exception as e:
                print(f"Error saving intermediate results: {e}")

    # Save final results
    try:
        with open(output_file, 'w') as outfile:
            json.dump(results, outfile, indent=4)
        print(f"Final results saved to {output_file}.")
    except Exception as e:
        print(f"Error saving final results: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process JSON file to extract segment IDs and paths.')
    parser.add_argument('-j', '--json_file', type=str, help='Path to the JSON file containing chromosome and position data.')
    parser.add_argument('-x', '--gbz_file', type=str, required=True, help='Path to the GBZ index or graph file.')

    args = parser.parse_args()
    main(args.json_file, args.gbz_file)
