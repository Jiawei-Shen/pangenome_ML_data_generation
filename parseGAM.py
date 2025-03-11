import vg_pb2
import pysam
import struct


def parse_gam(file_path):
    alignments = []

    with pysam.BGZFile(file_path, "rb") as f:
        while True:
            aln = vg_pb2.Alignment()
            try:
                length_bytes = f.read(4)  # Read length prefix
                if not length_bytes:
                    break

                length = struct.unpack("<I", length_bytes)[0]
                data = f.read(length)

                if len(data) != length:
                    print("Unexpected end of file while reading alignment.")
                    break

                aln.ParseFromString(data)
                alignments.append(aln)

            except Exception as e:
                print(f"Error parsing alignment: {e}")
                break

    return alignments


# Use the original BGZF-compressed GAM file
gam_file = "/scratch/jshen/Qichen_data/test_gam/LIB027514_223KKHLT4_S13_L006_001_hprc_hg38_v11.gam"
alignments = parse_gam(gam_file)

# Print first 5 alignments
for i, aln in enumerate(alignments[:5]):
    print(f"Alignment {i + 1}: Read Name: {aln.name}, Sequence: {aln.sequence}")
