import subprocess
import json
import argparse
import re
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed


def extract_nodes_from_gfa(gfa_file, reference_name, chromosome, output_json="nodes.json", threads=4):
    """Extract nodes from GFA path and directly assign node information if available using multi-threading."""
    print("[INFO] Starting path extraction...")
    command = f"grep '^W' {gfa_file} | grep '{reference_name}' | grep '{chromosome}'"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    path_nodes = {}
    processed_count = 0
    for line in result.stdout.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) >= 7:
            sequence = parts[6].strip()
            nodes = re.findall(r'[><][^><]+', sequence)
            for n in nodes:
                strand = n[0]
                node_id = n[1:]
                path_nodes[node_id] = {"strand": strand}
                processed_count += 1
                if processed_count % 10000 == 0:
                    print(f"[INFO] Extracted {processed_count} nodes from paths...")

    print(f"[✔] Path extraction complete. {processed_count} nodes extracted.")

    path_node_ids = set(path_nodes.keys())
    lock = threading.Lock()
    node_data = {}
    processed_count = 0

    def process_chunk(chunk):
        local_data = {}
        nonlocal processed_count
        for line in chunk:
            if line.startswith('S'):
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    node_id = parts[1]
                    if node_id in path_node_ids:
                        sequence = parts[2]
                        start_offset = int(parts[5].split(':')[2])  # Extract SO from parts[5]
                        local_data[node_id] = {
                            "sequence": sequence,
                            "length": len(sequence),
                            "start_offset": start_offset
                        }
                        processed_count += 1
                        if processed_count % 10000 == 0:
                            print(f"[INFO] Processed {processed_count} nodes with sequence data...")

        with lock:
            node_data.update(local_data)

    print("[INFO] Starting segment processing...")
    with open(gfa_file, 'r') as gfa:
        lines = gfa.readlines()
        chunk_size = len(lines) // threads
        chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(process_chunk, chunk) for chunk in chunks]
        for future in as_completed(futures):
            future.result()

    print(f"[✔] Segment processing complete. {processed_count} nodes with sequences processed.")

    for node_id in path_nodes.keys():
        if node_id in node_data:
            path_nodes[node_id].update(node_data[node_id])

    with open(output_json, 'w') as f:
        json.dump({"nodes": path_nodes}, f, indent=2)

    print(f"[✔] Extracted nodes with IDs, strands, sequences, lengths, and start offsets using {threads} threads, and saved to {output_json}")



def load_nodes(nodes_json):
    """Load node IDs, strands, sequences, and lengths from JSON."""
    with open(nodes_json, "r") as f:
        data = json.load(f)
    return data["nodes"]


def process_read(line, node_info, node_read_map, lock, processed_count, output_json):
    """Check if the read aligns to any node and store it grouped by node."""
    try:
        read = json.loads(line)
        read_info = {
            "read_name": read["name"],
            "sequence": read["sequence"],
            "mapping_quality": read.get("mapping_quality", 0),
            "score": read.get("score", None),
            "quality": read.get("quality", ""),
            "path": read.get("path", {})
        }
        mapped_nodes = set()
        # print(read.get("path", {}).get("mapping", []))
        for mapping in read.get("path", {}).get("mapping", []):
            node_id = str(mapping["position"].get("node_id", ""))
            print("Hello: ", type(node_id), '\n', type(node_info[0]))
            if node_id in node_info:
                mapped_nodes.add(node_id)

        with lock:
            for node_id in mapped_nodes:
                if node_id not in node_read_map:
                    node_read_map[node_id] = {
                        "strand": node_info[node_id]["strand"],
                        "sequence": node_info[node_id]["sequence"],
                        "length": node_info[node_id]["length"],
                        "start_offset": node_info[node_id]["start_offset"],
                        "reads": []
                    }
                node_read_map[node_id]["reads"].append(read_info)

            processed_count[0] += 1
            if processed_count[0] % 5000000 == 0:
                print(f"[INFO] Processed {processed_count[0]} reads...")
                save_json(node_read_map, f"tmp/{output_json}_batch_{processed_count[0]}.json")
                node_read_map.clear()  # Clear memory after saving each batch

    except json.JSONDecodeError:
        return  # Skip invalid JSON reads


def filter_reads(input_gam, nodes_file, output_json, threads=4):
    """Filter reads from GAM that align to extracted nodes and group them by node in real-time using multi-threading."""
    print("[INFO] Starting read filtering...")
    node_info = load_nodes(nodes_file)
    node_read_map = {}
    lock = threading.Lock()
    processed_count = [0]

    process = subprocess.Popen(["vg", "view", "-a", input_gam], stdout=subprocess.PIPE, text=True)

    with ThreadPoolExecutor(max_workers=threads) as executor:
        for _ in executor.map(lambda line: process_read(line, node_info, node_read_map, lock, processed_count, output_json), process.stdout):
            pass

    print(f"[✔] Read filtering complete. {processed_count[0]} reads processed.")
    if node_read_map:
        save_json(node_read_map, f"tmp/{output_json}_batch_final.json")

    merge_json_files("tmp", output_json)
    print(f"[✔] Filtered reads grouped by node saved to {output_json}")


def save_json(data, output_file):
    """Saves the current JSON state to a file atomically."""
    temp_file = output_file + ".tmp"
    os.makedirs(os.path.dirname(temp_file), exist_ok=True)
    with open(temp_file, "w") as f:
        json.dump(data, f, indent=2)
    os.replace(temp_file, output_file)
    print(f"[✔] Saved progress to {output_file}")


def merge_json_files(temp_dir, output_file):
    """Merge all batch JSON files into one final JSON file."""
    print("[INFO] Merging all batch JSON files...")
    merged_data = {}

    batch_files = sorted(glob.glob(f"{temp_dir}/{output_file}_batch_*.json"))

    for batch_file in batch_files:
        print(f"  - Loading {batch_file}")
        with open(batch_file, 'r') as f:
            batch_data = json.load(f)
            for node_id, node_data in batch_data.items():
                if node_id not in merged_data:
                    merged_data[node_id] = node_data
                else:
                    merged_data[node_id]["reads"].extend(node_data["reads"])

    with open(output_file, 'w') as f:
        json.dump(merged_data, f, indent=2)
    print(f"[✔] Merged data saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Extract reads from a GAM file that align to specific nodes from a GFA file and group them by node.")
    parser.add_argument("-gfa", "--gfa_file", required=False, help="GFA graph file")
    parser.add_argument("-r", "--reference", required=False, help="Reference name to filter (e.g., GRCh38)")
    parser.add_argument("-c", "--chromosome", required=False, help="Chromosome name to filter (e.g., chr1)")
    parser.add_argument("-gam", "--gam", required=True, help="Input GAM file")
    parser.add_argument("-n", "--nodes", default="nodes.json", help="Output JSON file for extracted node data")
    parser.add_argument("-j", "--json", default="grouped_reads.json", help="Output JSON file grouping reads by node")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for processing")

    args = parser.parse_args()

    # extract_nodes_from_gfa(args.gfa_file, args.reference, args.chromosome, args.nodes, args.threads)
    filter_reads(args.gam, args.nodes, args.json, args.threads)


if __name__ == "__main__":
    main()
