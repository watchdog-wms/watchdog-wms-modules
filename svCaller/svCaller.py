import argparse, sys, os # type: ignore

current_folder = os.path.dirname(os.path.abspath(__file__))
python_classes_path = os.path.join(current_folder, "python_classes")

parser = argparse.ArgumentParser(description="This script detects deletions given certain parameters.")

parser.add_argument("-bed_graph", "--bed", help="Path to bedgraph file")
parser.add_argument("-minimum_cluster_difference", "--min_cld", help="The mininum distance of two clusters, at which they still get combined", default=100)
parser.add_argument("-minimum_size", "--min_size", help="Path to file that contains variants, called by bcftools", default=2)
parser.add_argument("-maximum_z_score", "--max_z", help="Maximum z score threshold for coverage analysis", default=0.0)
parser.add_argument("-maximum_direct_z_score", "--max_direct", help="Maximum direct z score threshold for coverage analysis", default=-2.5)
parser.add_argument("-maximum_local_z_score", "--max_local", help="Maximum local z score threshold for coverage analysis", default=-6.0)
parser.add_argument("-local_range", "--range", help="Size of range/region before a certain position, used for the determination of local z Score parameters", default=500)
parser.add_argument("-pseudo_count", "--pc", help="Pseudo count for coverages over positions", default=1)
parser.add_argument("-tolerance", "--tol", help="Tolerance of insertion positions mapped to deletions", default=0.8)
parser.add_argument("-bamfile", "--bam", help="Path to bam file used for clipping patter analysis")
parser.add_argument("-outpath_deletions", "--out_del", help="Path to output txt file containing deletions")
parser.add_argument("-outpath_insertions", "--out_ins", help="Path to output txt file containing insertions")
parser.add_argument("-maximum_pattern_diff", "--max_patt_diff", help="Maximum distance of peaks of clipped reads to count them as insertion", default=10)
parser.add_argument("-minimum_surounding_z_score", "--min_sur_z", help="Minimum local z score for clipping pattern analysis", default=50)
parser.add_argument("-window_size", "--ws", help="Size of the window, whose position are controlled to be significantly low", default=20)
parser.add_argument("-min_z_score", "--min_z", help="Minimum z score for clipping pattern analysis", default=10)
parser.add_argument("-get_clipping_file", "--get_clp_file", help="Set this paramters as a path to get a file containing for each position the number of clipped reads")
parser.add_argument("-minimum_number_of_reads", "--min_reads", help="Minimum number of reads at which a position is permitted to be a peak", default=10)
parser.add_argument("-genome_propagate", "--gen_prop", help="Number of propagation to determine genome start/end", default=1000)
parser.add_argument("-gap_counter", "--gap", help="Maximum number of permitted consecutive gaps/0-coverage positions during the determination of the genome start/end", default=5)
parser.add_argument("-reference", "--ref", help="Path to reference genome")
parser.add_argument("-first_weight_score", "--fir_ws", help="Primary threshold for the score, which results from the verification of deletions with clipped sequences", default=0.0)
parser.add_argument("-second_weight_score", "--sec_ws", help="Secondary, more stringent threshold for the score, which results from the verification of deletions with clipped sequences", default=1.0)
parser.add_argument("-consensus_path", "--con_path", help="Path to the file containing the consensus sequences")
parser.add_argument("-matrix_pseudo_count", "--mpc", help="Small pseudo count for the log used for the computation of the PWMs", default=1)
parser.add_argument("-minimum_length", "--min_length", help="The minimum length of a consensus sequence", default=10)
parser.add_argument("-clipping_verfication_range", "--clp_ver_range", help="The range about clipped positons of deletions, where consesus sequences are tried to match on", default=100)


args = parser.parse_args()

sys.path.append(python_classes_path)

from classes import deletion_caller, insertion_caller

min_sur_z = float(args.min_sur_z)
window_size = int(args.ws)
tolerance = float(args.tol)
min_z = float(args.min_z)
bedgraph_path = args.bed
output_path_deletions = args.out_del
output_path_insertions = args.out_ins
min_cld = int(args.min_cld)
min_size = int(args.min_size)
max_z = float(args.max_z)
max_direct_z = float(args.max_direct)
max_local_z = float(args.max_local)
local_range = int(args.range)
pseudo_count = args.pc
max_patt_diff = int(args.max_patt_diff)
bam = args.bam
min_reads = int(args.min_reads)
gen_prop = int(args.gen_prop)
gap_counter = int(args.gap)
reference = str(args.ref)
fir_ws = float(args.fir_ws)
sec_ws = float(args.sec_ws)
con_path = args.con_path
mpc = float(args.mpc)
min_length = int(args.min_length)
clp_ver_range = int(args.clp_ver_range)



# Run deletionCaller
del_caller = deletion_caller(bedgraph_path, min_cld, min_size, max_z, max_direct_z, max_local_z, local_range, pseudo_count, gen_prop, gap_counter)
deletions = del_caller.get_deletion_clusters()
genome_start, genome_end = del_caller.get_genome_positions()
print("== Finished Coverage-Analysis ==")

# Run insertion caller
ins_caller = insertion_caller(bam, max_patt_diff, min_z, min_sur_z, window_size, tolerance, genome_start, genome_end, pseudo_count, min_reads, deletions, reference, fir_ws, sec_ws, con_path, mpc, min_length, clp_ver_range)
insertions, verified_deletions = ins_caller.get_insertions_and_verified_deletions()
ins_caller.write_insertions_to_file(output_path_insertions, insertions)
ins_caller.write_deletions_to_file(output_path_deletions, verified_deletions)

print("== Finished Clipping-Analysis and -Verification ==")

if args.get_clp_file is not None:
    ins_caller.write_clippings_to_file(args.get_clp_file)