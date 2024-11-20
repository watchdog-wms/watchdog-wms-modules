import argparse

parser = argparse.ArgumentParser(description="This script identifies the strain of the input data.")
parser.add_argument("-input_snps", "--input", help="Path to the file containing (replicate & variant caller) consistent SNPs of input data.")
parser.add_argument("-reference_snps", "--ref", help="Path to the txt-file containing the reference SNPs.", default = "")
parser.add_argument("-outpath", "--out", help="Path to the output txt-file containing the strain and further information.")
parser.add_argument("-config", "--con", help="Path to a config file containing the sample-strain identifications.")


args = parser.parse_args()

strains = {}
strain_distances = {}
predicted_strain = {}

# Read in config file that contains reference samples and their known afiliation to a certain virus strain
with open(args.con, "r") as file:
    for line in file:
        line = line.strip()
        if line != "" and not line.startswith("SAMPLE"):
            line_splitted = line.split("\t")
            sample = line_splitted[0]
            strain = line_splitted[1]
            if strain not in strains.keys():
                strains[strain] = [sample]
            else:
                strains[strain].append(sample)

# Dictionary containing the reference datasets as keys and the SNPs as values (in list)
data_snp = {}

consistent_snps = {}

# Open and extract reference SNPs file
with open(args.ref, "r") as file:
    for line in file:
        if line.strip() == "":
            continue
        
        line_splitted = line.split("\t")
        chr = line_splitted[0]
        dataset = line_splitted[1]
        position = line_splitted[2].strip()
        
        if chr not in data_snp.keys():
            data_snp[chr] = {}
        
        if dataset in data_snp[chr].keys():
            data_snp[chr][dataset].append(position)
        else:
            data_snp[chr][dataset] = [position]
file.close()

# Open input data/consistent SNPs
with open(args.input, "r") as file:
    for line in file:
        line = line.strip()
        if line == "":
            continue
        line_splitted = line.split("\t")
        chr = line_splitted[0]
        position = line_splitted[1]
        if chr not in consistent_snps.keys():
            consistent_snps[chr] = [position]
        else:
            consistent_snps[chr].append(position)  
file.close()

# Dictionary having the datasets as keys and their corresponding distance to the input data as value
distance = {}
data_info = {}

# Calculate distance of consistent SNPs to each dataset
for chr in data_snp.keys():
    
    if chr not in consistent_snps.keys():
        continue
    
    for dataset in data_snp[chr].keys():
        
        union = []
        intersection = []
        
        # Iterate over reference SNPs for each dataset
        for snp_ref in data_snp[chr][dataset]:
            # Append reference SNP to union anyway
            union.append(snp_ref)
            # If reference SNP is also in consistent SNPs append it to intersection as well
            if snp_ref in consistent_snps[chr]:
                intersection.append(snp_ref)
        # Append also all consistent SNPs to union
        for snp_in in consistent_snps[chr]:
            union.append(snp_in)
        
        # Unique to get sets
        union = list(set(union))
        intersection = list(set(intersection))
        # Compute distance between reference SNPs of a certain sample and the consistent SNPs
        tmp_distance = len(union) - len(intersection)
        # Safe this distance and also the size of union and intersection
        if chr not in data_info.keys():
            data_info[chr] = {}
            distance[chr] = {}
        data_info[chr][dataset] = [len(union), len(intersection)]
        distance[chr][dataset] = int(tmp_distance)
        print(dataset, tmp_distance)

# Clear output file
with open(args.out, "w") as file:
    pass
# Compute distance to strains
for chr in distance.keys():
    strain_distances[chr] = {}
    for strain in strains.keys():
        # Initialize all distance, union and intersection of a certain strain with 0
        strain_distances[chr][strain] =  [0, 0, 0]
        # Get the number of samples that correspond to this strain
        num_strain = len(strains[strain])
        
        # Iterate over all distances of the reference samples to the consistent SNPs
        for dataset in distance[chr].keys():
            # Only sum up the distances, unions and intersections of reference samples that correspond to the strain of interest that is currently iterated
            if dataset in strains[strain]:
                strain_distances[chr][strain][0] += distance[chr][dataset]
                strain_distances[chr][strain][1] += data_info[chr][dataset][0]
                strain_distances[chr][strain][2] += data_info[chr][dataset][1]
        # Average these values
        strain_distances[chr][strain][0] = round(strain_distances[chr][strain][0]/num_strain, 2)
        strain_distances[chr][strain][1] = round(strain_distances[chr][strain][1]/num_strain, 2)
        strain_distances[chr][strain][2] = round(strain_distances[chr][strain][2]/num_strain, 2)
    
    # Just a high number
    min_strain_dist = 100000000
    # Compute strain with lowest difference -> This is the prediction
    for strain in strain_distances[chr].keys():
        strain_dist = strain_distances[chr][strain][0]
        if strain_dist < min_strain_dist:
            min_strain_dist = strain_dist
            predicted_strain[chr] = [strain, min_strain_dist]



with open(args.out, "w") as file:
    for chr in strain_distances.keys():
        file.write(f"Predicted strain for {chr}: {predicted_strain[chr][0]} (Distance: {predicted_strain[chr][1]})\n\n")
        file.write("# All distances:\n")
        for strain in strain_distances[chr].keys():
            file.write(f"# {strain} (Distance: {strain_distances[chr][strain][0]}; Union: {strain_distances[chr][strain][1]}; Intersection: {strain_distances[chr][strain][2]})\n") 
        file.write("\n")
file.close()
