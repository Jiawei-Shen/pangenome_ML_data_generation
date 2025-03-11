import vg_pb2  # Ensure you have the protobuf definitions
import gzip

def parse_gam(gam_file):
    alignments = []
    with gzip.open(gam_file, "rb") as f:
        while True:
            size_bytes = f.read(4)
            if not size_bytes:
                break
            size = int.from_bytes(size_bytes, byteorder="little")
            data = f.read(size)
            alignment = vg_pb2.Alignment()
            alignment.ParseFromString(data)
            alignments.append(alignment)
    return alignments

if __name__ == "__main__":
    gam_file = "/scratch/jshen/Qichen_data/test_gam/LIB027514_223KKHLT4_S13_L006_001_hprc_hg38_v11.gam"  # Change this to your actual GAM file path
    alignments = parse_gam(gam_file)

    # Print some parsed alignments
    for aln in alignments[:3]:  # Print first 3 alignments
        print(aln)
