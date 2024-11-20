import argparse, os
parser = argparse.ArgumentParser(description="This file computes replicate and variant caller consistent SNPs.")
parser.add_argument("-bcftools_replicates", "--bcf_rep", help="Path to VCF-file containing variants of first replicate called by bcftools.")
parser.add_argument("-varscan_replicates", "--var_rep", help="Path to VCF-file containing variants of first replicate called by Varscan.")
parser.add_argument("-outpath", "--out", help="Path to output file.")
args = parser.parse_args()

# Get vcf files of all replicates called by bcftools
bcftools_replicates = str(args.bcf_rep).split(",")
bcftools_replicates = [rep.strip() for rep in bcftools_replicates]

# Remove NA replicates
bcftools_replicates = [path for path in bcftools_replicates if os.path.exists(path)]

# Get vcf files of all replicates called by Varscan
varscan_replicates = str(args.var_rep).split(",")
varscan_replicates = [rep.strip() for rep in varscan_replicates]

# Remove NA replicates
varscan_replicates = [path for path in varscan_replicates if os.path.exists(path)]

# Get total number of predictions/vcf files
replicate_number = len(bcftools_replicates) + len(varscan_replicates)

if any(rep in varscan_replicates for rep in bcftools_replicates):
    raise Exception("Files created by bcftools and Varscan have identical names.")

snp_lists = {}

# Iterate over all vcf files called by bcftools and extract SNPs
for replicate in bcftools_replicates:
    with open(replicate, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("#") or line == "":
                continue
            
            # Get information of SNP
            line_splitted = line.split("\t")
            chr = line_splitted[0]
            position = line_splitted[1]
            reference = line_splitted[3]
            alternative = line_splitted[4]
            info = line_splitted[7]
            
            # Remove INDELs (since we only want to compute consistent SNPs)
            if info.startswith("INDEL"):
                continue
            # SNP list is separated first via chromosome (so if chr does not exist yet, introduce it for SNP list). Same story for replicate down below.
            if chr not in snp_lists.keys():
                snp_lists[chr] = {}
            if replicate not in snp_lists[chr].keys():
                snp_lists[chr][replicate] = {}
            
            # Add reference and alternative nucleotide to the position of SNP on the given chromosome in the underlying sample
            snp_lists[chr][replicate][position] = [reference, alternative]
    file.close()

# Do exact same stuff for replicates analyzed by Varscan
for replicate in varscan_replicates:
    with open(replicate, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith("Chrom") or line == "":
                continue
            
            line_splitted = line.split("\t")
            chr = line_splitted[0]
            position = line_splitted[1]
            reference = line_splitted[2]
            alternative = line_splitted[3]
            
            if chr not in snp_lists.keys():
                snp_lists[chr] = {}
            if replicate not in snp_lists[chr].keys():
                snp_lists[chr][replicate] = {}
                
            snp_lists[chr][replicate][position] = [reference, alternative]
    file.close()

consistent_snps = {}

# Iterating over one of these for dictionaries is enough, since we are only interested in their intersection (these are then the consistent SNPs)
for chr in snp_lists.keys():
    if len(snp_lists[chr].keys()) != replicate_number:
        continue
    first_replicate = next(iter(snp_lists[chr]))
    
    # Iterate over all SNPs/positions
    for pos in snp_lists[chr][first_replicate].keys():
        passed_all_reps = True
        # Iterate over all other replicates and check if SNP is present in those as well
        for replicate in snp_lists[chr].keys():
            if replicate == first_replicate:
                continue
            # Check if SNP is present in other replicate
            if pos in snp_lists[chr][replicate]:
                # Check also if alternative and reference nucleotide of the SNPs are identical as well
                if snp_lists[chr][first_replicate][pos][0] != snp_lists[chr][replicate][pos][0] or snp_lists[chr][first_replicate][pos][1] != snp_lists[chr][replicate][pos][1]:
                    # If the SNPs alternative or reference nucleotide is unequal to the given SNP from above, it is not consistent!
                    passed_all_reps = False
                    break
            # If SNP is not present in at least one of the other replicates, it is not consistent and thus has not passed all replicates -> set boolean parameter on false
            else:
                passed_all_reps = False
                break
        # If SNP/position is present and its alternative and reference nucleotide is identical in all replicates -> Consistent SNP
        if passed_all_reps:
            if chr not in consistent_snps.keys():
                consistent_snps[chr] = {}
            consistent_snps[chr][pos] = [snp_lists[chr][first_replicate][pos][0], snp_lists[chr][first_replicate][pos][1]]

# Write consistent SNPs in output file
with open(args.out, "w") as file:
    file.write(f"CHR\tPOS\tREF\tALT\n")
    for chr in consistent_snps.keys():
        for snp in consistent_snps[chr]:
            file.write(f"{chr}\t{snp}\t{consistent_snps[chr][snp][0]}\t{consistent_snps[chr][snp][1]}\n")
file.close()
