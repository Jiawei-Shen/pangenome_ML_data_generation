import vg_pb2  # Import generated protobuf classes
import pysam  # For reading BGZF files


def parse_gam(file_path):
    alignments = []
    with pysam.BGZFile(file_path, "rb") as f:  # Open BGZF-compressed GAM file
        while True:
            aln = vg_pb2.Alignment()
            try:
                # Read the 4-byte little-endian length prefix
                length_bytes = f.read(4)
                if not length_bytes:
                    break  # End of file

                length = int.from_bytes(length_bytes, byteorder="little")
                data = f.read(length)  # Read the full protobuf message

                if len(data) != length:
                    raise ValueError("Unexpected end of file while reading alignment.")

                aln.ParseFromString(data)  # Deserialize the Protobuf message
                alignments.append(aln)

            except Exception as e:
                print(f"Error parsing alignment: {e}")
                break

    return alignments


# Example usage
gam_file = "/scratch/jshen/Qichen_data/test_gam/LIB027514_223KKHLT4_S13_L006_001_hprc_hg38_v11.gam"
alignments = parse_gam(gam_file)

# Print first 5 alignments
for aln in alignments[:5]:
    print(f"Read Name: {aln.name}, Sequence: {aln.sequence}")
