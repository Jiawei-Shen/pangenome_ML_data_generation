import argparse
import gzip
import json
from google.protobuf.internal.decoder import _DecodeVarint32
from vg_pb2 import Alignment  # Import the compiled Protobuf message
import pysam


def read_gam_bgzf(filepath):
    with pysam.tabix_iterator(filepath) as f:
        data = f.read()
        print("Read", len(data), "bytes from BGZF-compressed GAM file")


def read_gam(file_path):
    """Reads a GAM file and parses alignments."""
    alignments = []

    with open(file_path, "rb") as f:  # GAM files are BGZF compressed
        buf = f.read()
        pos = 0
        print("Done reading the GAM file into the buffer.")

        while pos < len(buf):
            msg_len, new_pos = _DecodeVarint32(buf, pos)  # Read message length
            pos = new_pos
            print(msg_len, new_pos)
            message = buf[pos: pos + msg_len]  # Extract message
            pos += msg_len

            alignment = Alignment()
            alignment.ParseFromString(message)  # Parse Protobuf message
            alignments.append(alignment)

    return alignments


def gam_to_json(alignments):
    """Converts alignments to JSON format."""
    return [json.loads(aln.SerializeToJsonString()) for aln in alignments]


def main():
    parser = argparse.ArgumentParser(description="Parse a binary GAM file and output alignments in JSON format.")
    parser.add_argument("input_gam", default="filtered.gam", help="Path to the input GAM file")
    parser.add_argument("-o", "--output", help="Path to the output JSON file (optional)")
    parser.add_argument("-n", "--num", type=int, default=5, help="Number of alignments to display (default: 5)")

    args = parser.parse_args()

    print(f"Reading GAM file: {args.input_gam}")
    # read_gam_bgzf(args.input_gam)

    alignments = read_gam(args.input_gam)

    json_output = gam_to_json(alignments)

    if args.output:
        with open(args.output, "w") as out_file:
            json.dump(json_output, out_file, indent=2)
        print(f"Saved output to {args.output}")
    else:
        # Print the first N alignments
        print(json.dumps(json_output[:args.num], indent=2))


if __name__ == "__main__":
    main()
