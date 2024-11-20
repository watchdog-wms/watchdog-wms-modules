import re, argparse

parser = argparse.ArgumentParser(description="This file matches variants to gene features.")
parser.add_argument("-gtf_path", "--gtf", help="Path to gtf file, which variants are matched on.")
parser.add_argument("-inpath", "--infile", help="Path to the input file, containing the variants (e.g. SNPs or deletions).")
parser.add_argument("-outpath", "--out", help="Path to output file.")
parser.add_argument("-mode", "--m", help="Select mode of variant matching. SNP for matching SNPs, DELETION for matching deletions and INSERTION for matching insertions.")
args = parser.parse_args()

gtf_path = args.gtf
out_path = args.out
input_file = args.infile
mode = args.m
inputs = {}

# Read in variants and save them in the inputs dictionary
with open(input_file, "r") as file:
    # If variants are deletions
    if mode == "DELETION":
        for line in file:
            line = line.strip()
            if line != "" and not line.startswith("#") and not line.startswith("CHR") and not line.startswith("="):
                line_splitted = line.split("\t")
                chr = line_splitted[0]
                start = int(line_splitted[1])
                end = int(line_splitted[2].strip())
                if chr not in inputs.keys():
                    inputs[chr] = [[start, end]]
                else:
                    inputs[chr].append([start, end])
    # If variants are SNPs
    elif mode == "SNP":
        for line in file:
            line = line.strip()
            if line != "" and not line.startswith("CHR"):
                line_splitted = line.strip().split("\t")
                chr = line_splitted[0]
                pos = int(line_splitted[1])
                ref = line_splitted[2]
                alt = line_splitted[3].strip()
                if chr not in inputs.keys():
                    inputs[chr] = [[pos, ref, alt]]
                else:
                    inputs[chr].append([pos, ref, alt])
    # If variants are insertions
    elif mode == "INSERTION":
        for line in file:
            line = line.strip()
            if line != "" and not line.startswith("#") and not line.startswith("CHR") and not line.startswith("="):
                line_splitted = line.split("\t")
                chr = line_splitted[0]
                position = int(line_splitted[1].strip())
                if chr not in inputs.keys():
                    inputs[chr] = [[position]]
                else:
                    inputs[chr].append([position])
    # If mode is not given or written incorrectly, raise an exception
    else:
        raise Exception("Mode not found. Mode has to be 'SNP', 'INSERTION' or 'DELETION' in upper case letters")

# If variant not in gene, skip all lines that correspond to this gene
skip_gene = False
# If variant not in transcript, skip all lines that correspond to this transcript
skip_transcript = False
# Regular expression for gene identifiers
pattern_gene = r'gene_id "([^"]+)"'
# Regular expression for transcript identifiers
pattern_transcript = r'transcript_id "([^"]+)"'

# Function that checks if a given variant is in a certain feature
def check_variant_hit(variant, feature):
    
    # Get variant start
    var_start = int(variant[0])
    # Get feature start and end
    feature_start = feature[0]
    feature_end = feature[1]
    
    # If variants are deletions also variant end is required
    if mode == "DELETION":
        var_end = variant[1]
    # If SNPs or insertions are matched on features, start and end of variants are equal
    else:
        var_end = var_start
    
    # Compute size of the feature
    feature_size = feature_end - feature_start
    
    # Size of the variant
    size = -1
    
    # If variant end is within the feature...
    if var_end >= feature_start and var_end <= feature_end:
        # ...Check if variant start is before the feature start -> Matching area/size = variant end - feature start
        if var_start <= feature_start:
            size = var_end - feature_start
        # ...Else if variant start is also within feature -> Matching area/size = variant end - variant start
        else:
            size = var_end - var_start
        return True, f"{size}bp ({round(size*100/feature_size, 1)}%)"
    # If variant start is within feature (variant end can not be in feature as well, since this case is checked above in the first if statement)
    elif var_start >= feature_start and var_start <= feature_end:
        # Matching area/size of this variant on this particular feature is obviously the feature end - variant start
        size = feature_end - var_start
        return True, f"{size}bp ({round(size*100/feature_size, 1)}%)"
    # If variant start is before feature start and variant end is behind feature end, the matchin are on the feature is simply the feature size
    elif var_start <= feature_start and var_end >= feature_end:
        size = feature_end - feature_start
        return True, f"{size}bp ({round(size*100/feature_size, 1)}%)"
    # Edge case, can be ignored
    else:
        return False, size
        
# Given exons, compute introns: Simply the sequences between the exons
def compute_introns(exons):
    introns = []
    for i in range(0, len(exons)):
        exon = exons[i]
        # Check if next exon and thus a sequence between two exons is available
        if i+1 < len(exons):
            next_exon = exons[i+1]
            intron_start = exon[1]+1
            intron_end = next_exon[0]-1
            introns.append([intron_start, intron_end])
    
    return introns
        


# Clear file, since information are written parallel to the analysis
with open(out_path, "w") as writer:
    writer.write(f"# This file contains annotated {mode}s\n")

# Iterate over all chromosomes
for chr in inputs.keys():
    
    # Initialize gene ID, transcript ID and exon list (empty)
    gene_id = ""
    transcript_id = ""
    exon_list = []
    
    # Iterate over all variants of the input (can only either be SNPs, insertions or deletions)
    for variant in inputs[chr]:
        # Open output file and write gained information instantly
        with open(out_path, "a") as writer:
            # Simply write type of variant and its corrsponding information/data in the output file
            var_start = int(variant[0])
            if mode == "DELETION":
                var_end = variant[1]
                writer.write(f"\n{chr}\t{mode}\t{var_start}-{var_end}\n")
            elif mode == "SNP":
                var_ref = variant[1]
                var_var = variant[2]
                writer.write(f"\n{chr}\t{mode}\t{var_start}\t{var_ref} -> {var_var}\n")
            else:
                writer.write(f"\n{chr}\t{mode}\t{var_start}\n")

            # Now GTF file is opened and all lines are iterated
            with open(gtf_path, "r") as file:
                for line in file:
                    
                    # Get information of line (chromosome, feature (identifier), feature start and end as well as general information about the feature)
                    line_split = line.split("\t")
                    chrom = line_split[0]
                    if chrom != chr:
                        continue
                    
                    ident = line_split[2]
                    start = int(line_split[3])
                    end = int(line_split[4])
                    infos = line_split[8]
                    
                    # If a certain transcript should be skipped (due to mismatch of the variant), skip all lines that don't correspond to a new transcript or gene
                    if ident != "transcript" and ident != "gene" and skip_transcript:
                        continue
                    # If a certain gene should be skipped (due to mismatch of the variant), skip all lines that don't correspond to a new gene
                    if ident != "gene" and skip_gene:
                        continue
                    
                    # If line contains information about gene
                    if ident == "gene":
                        # Set skipping parameter to false, until the match of the variant on this gene is examined
                        skip_gene = False
                        skip_transcript = False
                        
                        # If a gene was analyzed and the variant was matched on before, the exon list can contain exons -> Compute introns and check for matches
                        if len(exon_list) > 0:
                            intron_list = compute_introns(exon_list)
                            exon_list = []
                            for intron in intron_list:
                                is_hit, overlap = check_variant_hit(variant, intron)
                                if is_hit:
                                    intron_start = intron[0]
                                    intron_end = intron[1]
                                    if mode == "DELETION":
                                        writer.write(f"\t\tIntron\t{intron_start}\t{intron_end}\tOverlap: {overlap}\n")
                                    else:
                                        writer.write(f"\t\tIntron\t{intron_start}\t{intron_end}\n")
                        
                        # Check if variant matches and get overlap size
                        is_hit, overlap = check_variant_hit(variant, [start, end])
                        
                        # If gene does not match, skip all lines that correspond to the gene
                        if not is_hit:
                            skip_gene = True
                            continue
                        # If match of variant on feature was successful, get gene ID and write all information gathered in output file
                        match_gene = re.search(pattern_gene, infos)
                        if match_gene:
                            gene_id = match_gene.group(1)
                            if mode == "DELETION":
                                writer.write(f"Gene\t{gene_id}\t{start}\t{end}\tOverlap: {overlap}\n")
                            else:
                                writer.write(f"Gene\t{gene_id}\t{start}\t{end}\n")
                        else:
                            skip_gene = True
                            continue
                    # Exact same story for transcript    
                    elif ident == "transcript":
                        skip_transcript = False
                        if len(exon_list) > 0:
                            intron_list = compute_introns(exon_list)
                            exon_list = []
                            for intron in intron_list:
                                is_hit, overlap = check_variant_hit(variant, intron)
                                if is_hit:
                                    intron_start = intron[0]
                                    intron_end = intron[1]
                                    if mode == "DELETION":
                                        writer.write(f"\t\tIntron\t{intron_start}\t{intron_end}\tOverlap: {overlap}\n")
                                    else:
                                        writer.write(f"\t\tIntron\t{intron_start}\t{intron_end}\n")
                        is_hit, overlap = check_variant_hit(variant, [start, end])
                        if not is_hit:
                            skip_transcript = True
                            continue
                        
                        match_transcript = re.search(pattern_transcript, infos)
                        if match_transcript:
                            transcript_id = match_transcript.group(1)
                            if mode == "DELETION":
                                writer.write(f"\tTranscript\t{transcript_id}\t{start}\t{end}\t{infos.strip()}\tOverlap: {overlap}\n")
                            else:
                                writer.write(f"\tTranscript\t{transcript_id}\t{start}\t{end}\t{infos.strip()}\n")
                        else:
                            skip_transcript = True
                            continue
                    # If feature is an exon check if variant matches on it. Independent of the result, append exon to exon list
                    elif ident == "exon":
                        is_hit, overlap = check_variant_hit(variant, [start, end])
                        if is_hit:
                            if mode == "DELETION":
                                writer.write(f"\t\tExon\t{start}\t{end}\t{infos.strip()}\tOverlap: {overlap}\n")
                            else:
                                writer.write(f"\t\tExon\t{start}\t{end}\t{infos.strip()}\n")
                        exon_list.append([start, end])
                    # If feature is a coding sequence, check if variant matches
                    elif ident == "CDS":
                        is_hit, overlap = check_variant_hit(variant, [start, end])
                        if is_hit:
                            if mode == "DELETION":
                                writer.write(f"\t\tCDS\t{start}\t{end}\t{infos.strip()}\tOverlap: {overlap}\n")
                            else:
                                writer.write(f"\t\tCDS\t{start}\t{end}\t{infos.strip()}\n")
            # After all lines have been iterated, check the exon list of the final gene or transcript         
            if len(exon_list) > 0:
                intron_list = compute_introns(exon_list)
                exon_list = []
                for intron in intron_list:
                    is_hit, overlap = check_variant_hit(variant, intron)
                    if is_hit:
                        intron_start = intron[0]
                        intron_end = intron[1]
                        if mode == "DELETION":
                            writer.write(f"\t\tIntron\t{intron_start}\t{intron_end}\tOverlap: {overlap}\n")
                        else:
                            writer.write(f"\t\tIntron\t{intron_start}\t{intron_end}\n")
            file.close()
        writer.close()
                