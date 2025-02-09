import subprocess
import json
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor


def extract_nodes(graph_xg, path_name, output_nodes="nodes_in_path.txt"):
    """Extract node IDs belonging to a given path using vg paths and vg view."""
    cmd = f"vg paths -x {graph_xg} -p {path_name} -V | vg view -v - | jq -r '.node[].id' > {output_nodes}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"[✔] Extracted node IDs for path '{path_name}' into {output_nodes}")


def process_read(line, node_ids):
    """Check if the read aligns to any of the extracted node IDs and group by node."""
    try:
        read = json.loads(line)
        read_info = {
            "read_name": read["name"],
            "sequence": read["sequence"],
            "mapping_quality": read.get("mapping_quality", 0)
        }

        mapped_nodes = set()  # Store unique node IDs the read maps to
        for mapping in read.get("path", {}).get("mapping", []):
            node_id = str(mapping["position"].get("node_id", ""))
            if node_id in node_ids:
                mapped_nodes.add(node_id)

        return [(node_id, read_info) for node_id in mapped_nodes] if mapped_nodes else None

    except json.JSONDecodeError:
        return None  # Return None if no match is found


def filter_reads(input_gam, nodes_file, output_json, threads=4):
    """Filter reads from GAM that align to the extracted nodes and group by node."""
    # Load node IDs into a set
    with open(nodes_file, "r") as f:
        node_ids = {line.strip() for line in f}

    # Process GAM file with multi-threading
    node_read_map = defaultdict(list)  # Dictionary to store reads grouped by node
    process = subprocess.Popen(["vg", "view", "-a", input_gam], stdout=subprocess.PIPE, text=True)

    processed_count = 0  # Counter for processed reads
    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = []
        for result in executor.map(lambda line: process_read(line, node_ids), process.stdout):
            if result:
                results.append(result)

            # Print progress every 1000 reads
            processed_count += 1
            if processed_count % 1000 == 0:
                print(f"[INFO] Processed {processed_count} reads...")

    # Collect results and group reads by node
    for result in results:
        for node_id, read_info in result:
            node_read_map[node_id].append(read_info)

    # Save to JSON
    with open(output_json, "w") as f:
        json.dump(node_read_map, f, indent=2)

    print(f"[✔] Filtered reads grouped by node saved to {output_json}")
    return output_json


def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(
        description="Extract reads from a GAM file that align to a specific path and group them by node.")

    parser.add_argument("-x", "--xg", required=True, help="Graph XG file")
    parser.add_argument("-g", "--gam", required=True, help="Input GAM file")
    parser.add_argument("-p", "--path", required=True, help="Path name to filter reads")
    parser.add_argument("-n", "--nodes", default="nodes_in_path.txt", help="Output file for extracted node IDs")
    parser.add_argument("-j", "--json", default="filtered_reads_by_node.json",
                        help="Output JSON file grouping reads by node")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for processing")

    args = parser.parse_args()

    # Run the extraction pipeline
    extract_nodes(args.xg, args.path, args.nodes)
    filter_reads(args.gam, args.nodes, args.json, args.threads)


if __name__ == "__main__":
    main()
