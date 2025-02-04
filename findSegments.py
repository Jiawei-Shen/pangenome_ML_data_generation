import subprocess
import json
import sys


def run_vg_find_position(position):
    """Find the node corresponding to a given position."""
    cmd = f'/home/jiawei/anaconda3/envs/pange/bin/vg find -x ./Pangenome_Graph/hprc-v1.1-mc-grch38.d9.gbz -p {position} | vg view -j -'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running vg find for position:", result.stderr)
        return None
    return json.loads(result.stdout)


def run_vg_find_neighbors(node_id, context=5):
    """Find neighboring nodes within the given context range."""
    cmd = f'/home/jiawei/anaconda3/envs/pange/bin/vg find -x ./Pangenome_Graph/hprc-v1.1-mc-grch38.d9.gbz -n {node_id} -c {context} | vg view -j -'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running vg find for neighbors:", result.stderr)
        return None
    return json.loads(result.stdout)


def main():
    if len(sys.argv) != 2:
        print("Usage: python findSegments.py <position>")
        sys.exit(1)

    position = sys.argv[1]
    position_data = run_vg_find_position(position)
    if not position_data or "node" not in position_data:
        print("No node found for given position.")
        sys.exit(1)

    node_id = position_data["node"][0]["id"]
    print(f"Found node ID: {node_id}")

    neighbors_data = run_vg_find_neighbors(node_id)
    if not neighbors_data:
        print("No neighboring nodes found.")
        sys.exit(1)

    print(json.dumps(neighbors_data, indent=4))


if __name__ == "__main__":
    main()
