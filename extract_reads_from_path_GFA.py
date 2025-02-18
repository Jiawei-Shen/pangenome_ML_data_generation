import subprocess
import json
import argparse
import threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import re

def extract_nodes_from_gfa(gfa_file, reference_name, chromosome, output_json="nodes.json"):
    """Extract nodes with ID, length, and strand symbol from the GFA file based on reference and chromosome."""
    command = f"grep '^W' {gfa_file} | grep '{reference_name}' | grep '{chromosome}'"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    node_info = []
    for line in result.stdout.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) >= 7:
            sequence = parts[6].strip()

            # Use regex to split and keep delimiters
            nodes = re.findall(r'[><][^><]+', sequence)

            for n in nodes:
                strand = n[0]  # '>' or '<'
                node_id = n[1:]  # Node ID after the symbol
                node_info.append({
                    "id": node_id,
                    "strand": strand,
                    "length": len(node_id)  # Adjust as needed
                })

    with open(output_json, 'w') as f:
        json.dump({"nodes": node_info}, f, indent=2)

    print(f"[\u2714] Extracted node IDs, strands, and lengths, and saved to {output_json}")

def load_nodes(nodes_json):
    with open(nodes_json, "r") as f:
        data = json.load(f)
    node_info = {}
    for node in data["nodes"]:
        node_id = str(node["id"])
        node_info[node_id] = {"strand": node["strand"], "length": node["length"]}
    return node_info

def process_read(line, node_info, node_read_map, lock):
    try:
        read = json.loads(line)
        read_info = {"read_name": read["name"], "sequence": read["sequence"], "path": read.get("path", {})}
        mapped_nodes = set()
        for mapping in read.get("path", {}).get("mapping", []):
            node_id = str(mapping["position"].get("node_id", ""))
            if node_id in node_info:
                mapped_nodes.add(node_id)
        with lock:
            for node_id in mapped_nodes:
                if node_id not in node_read_map:
                    node_read_map[node_id] = {"strand": node_info[node_id]["strand"], "length": node_info[node_id]["length"], "reads": []}
                node_read_map[node_id]["reads"].append(read_info)
    except json.JSONDecodeError:
        return

def filter_reads(input_gam, nodes_file, output_json, threads=4):
    node_info = load_nodes(nodes_file)
    node_read_map = {}
    lock = threading.Lock()
    process = subprocess.Popen(["vg", "view", "-a", input_gam], stdout=subprocess.PIPE, text=True)
    with ThreadPoolExecutor(max_workers=threads) as executor:
        for _ in executor.map(lambda line: process_read(line, node_info, node_read_map, lock), process.stdout):
            pass
    save_json(node_read_map, output_json)
    print(f"[\u2714] Filtered reads grouped by node saved to {output_json}")

def save_json(data, output_file):
    with open(output_file, "w") as f:
        json.dump(data, f, indent=2)
    print(f"[\u2714] Saved progress to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract reads from a GAM file that align to specific nodes from a GFA file and group them by node.")
    parser.add_argument("-gfa", "--gfa_file", required=True, help="GFA graph file")
    parser.add_argument("-r", "--reference", required=True, help="Reference name to filter (e.g., GRCh38)")
    parser.add_argument("-c", "--chromosome", required=True, help="Chromosome name to filter (e.g., chr1)")
    parser.add_argument("-gam", "--gam", required=True, help="Input GAM file")
    parser.add_argument("-n", "--nodes", default="nodes.json", help="Output JSON file for extracted node data")
    parser.add_argument("-j", "--json", default="grouped_reads.json", help="Output JSON file grouping reads by node")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for processing")
    args = parser.parse_args()

    extract_nodes_from_gfa(args.gfa_file, args.reference, args.chromosome, args.nodes)
    filter_reads(args.gam, args.nodes, args.json, args.threads)

if __name__ == "__main__":
    main()
