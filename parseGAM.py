import pysam
import vg_pb2

def parse_gam(file_path):
    alignments = []
    with pysam.AlignmentFile(file_path, "r") as gam_file:
        for alignment in gam_file:
            aln = vg_pb2.Alignment()
            # Assuming the alignment is in the format you expect, decode it into vg_pb2
            # This assumes your GAM file is already in Protobuf format
            aln.ParseFromString(alignment.query)
            alignments.append(aln)
    return alignments

# Example usage
gam_file = "/scratch/jshen/Qichen_data/test_gam/LIB027514_223KKHLT4_S13_L006_001_hprc_hg38_v11.gam"
alignments = parse_gam(gam_file)

# Print some details
for aln in alignments[:5]:  # Print first 5 alignments
    print(aln.name, aln.sequence)
