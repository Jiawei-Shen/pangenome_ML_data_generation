import json
import subprocess
import sys


def call_find_segments(chromosome, position, ref="GRCh38"):
    """Call the findSegments.py script with the given chromosome and position."""
    formatted_position = f"{ref}#{chromosome}:{position}"
    cmd = ['python', 'findSegments.py', formatted_position]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Error processing {formatted_position}: {result.stderr}")
    else:
        print(f"Output for {formatted_position}:\n{result.stdout}")


def main(json_file):
    # Load JSON data
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)

    # # Check if the data is a dictionary of records
    # if not isinstance(data, dict):
    #     print("JSON file must contain a dictionary of records.")
    #     sys.exit(1)

    # Process each entry in the JSON
    for key, record in enumerate(data):
        chromosome = record.get('chromosome')
        position = record.get('position')

        if not chromosome or not position:
            print(f"Skipping record {key} due to missing chromosome or position.")
            continue

        print(f"\nProcessing {chromosome}:{position}")
        call_find_segments(chromosome, position)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python call_findSegments.py <data.json>")
        sys.exit(1)

    json_file = sys.argv[1]
    main(json_file)
