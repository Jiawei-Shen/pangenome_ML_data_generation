import json
import argparse
import subprocess
import threading
import queue

def load_segment_ids(json_file):
    """Loads segment IDs from the given JSON file and organizes them by region."""
    with open(json_file, "r") as f:
        data = json.load(f)

    segment_to_region = {}  # Maps segment_id to its corresponding region
    for region, details in data.items():
        if "segments" in details:
            for segment in details["segments"]:
                segment_id = segment["segment_id"]
                segment_to_region[segment_id] = region

    return segment_to_region

def process_gam_file(gam_file, segment_to_region, result_queue):
    """Iterates over a GAM file and assigns reads to the correct region based on segment alignment."""
    result_dict = {region: set() for region in set(segment_to_region.values())}

    # Stream the GAM file as JSON
    process = subprocess.Popen(['vg', 'view', '-a', gam_file], stdout=subprocess.PIPE, text=True)

    for line in process.stdout:
        read = json.loads(line.strip())

        if 'path' in read:
            for mapping in read['path']['mapping']:
                if 'position' in mapping and 'node_id' in mapping['position']:
                    segment_id = str(mapping['position']['node_id'])

                    if segment_id in segment_to_region:
                        region = segment_to_region[segment_id]
                        result_dict[region].add(read['name'])  # Store read names

    result_queue.put(result_dict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find reads mapped to nodes from a JSON file in a GAM file.")
    parser.add_argument("json_file", type=str, help="Input JSON file containing segment IDs.")
    parser.add_argument("gam_file", type=str, help="Input GAM file.")

    args = parser.parse_args()

    # Load segment IDs and their corresponding regions
    segment_to_region = load_segment_ids(args.json_file)
    print(f"Loaded {len(segment_to_region)} segment IDs from {args.json_file}")

    # Multi-threaded processing
    result_queue = queue.Queue()
    thread = threading.Thread(target=process_gam_file, args=(args.gam_file, segment_to_region, result_queue))
    thread.start()
    thread.join()

    # Collect results
    result_dict = result_queue.get()

    # Print results
    for region, reads in result_dict.items():
        print(f"\nRegion: {region}")
        if reads:
            print(f"  Reads mapped ({len(reads)} total):")
            for read in reads:
                print(f"    {read}")
        else:
            print("  No reads found for this region.")
