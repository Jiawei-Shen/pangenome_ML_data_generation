import struct
import vg_pb2  # Make sure you have compiled vg.proto to vg_pb2.py
import gzip


def read_varint(file):
    """ Reads a varint from a binary file """
    shift = 0
    result = 0
    while True:
        byte = file.read(1)
        if not byte:
            raise EOFError("Unexpected end of file while reading varint")
        b = ord(byte)
        result |= (b & 0x7F) << shift
        if not (b & 0x80):
            break
        shift += 7
    return result


def parse_gam(file_path):
    """ Reads a GAM file in binary format and parses it into Alignment objects """
    alignments = []

    with gzip.open(file_path, "rb") as f:  # Read BGZF (GAM) file in binary mode
        while True:
            aln = vg_pb2.Alignment()
            try:
                length = read_varint(f)  # Read the length of the next alignment
                data = f.read(length)  # Read the protobuf-encoded data

                if len(data) != length:
                    print("Unexpected end of file while reading alignment.")
                    break

                aln.ParseFromString(data)  # Parse the binary data into a Protobuf message
                alignments.append(aln)

            except EOFError:
                break
            except Exception as e:
                print(f"Error parsing alignment: {e}")
                break

    return alignments


# Path to your binary GAM file
gam_file = "/scratch/jshen/Qichen_data/test_gam/LIB027514_223KKHLT4_S13_L006_001_hprc_hg38_v11.gam"

# Read alignments
alignments = parse_gam(gam_file)

# Print the first 5 parsed alignments
for i, aln in enumerate(alignments[:5]):
    print(f"Alignment {i + 1}: Read Name: {aln.name}, Sequence Length: {len(aln.sequence)}")
