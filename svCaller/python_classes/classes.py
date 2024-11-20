import pandas as pd, numpy as np, re, pysam, warnings # type: ignore
warnings.simplefilter(action='ignore', category=FutureWarning)

# deletion_caller calls/identifies deletions based on a coverage analysis
class deletion_caller:
    
    # Ininitialize all global variables
    bedgraph_path = ""
    output_path = ""
    min_cld = -1
    min_size = -1
    max_z = -1.0
    max_direct_z = -1.0
    max_local_z = -1.0
    local_range = -1
    genome_start = -1
    genome_end = -1
    pseudo_count = -1
    
    deletions = {}
    genome_coverage = {}
    last_end_pos = -1
    
    global_cov_mean = -1.0
    global_cov_std = -1.0
    
    already_deletions = False
    
    local_cov_mean = -1.0
    local_cov_std = -1.0
    
    skip_position = False
    
    deletion_clusters = {}
    
    gen_prop = -1
    gaps = -1
    
    def __init__(self, bedgraph_path, min_cld, min_size, max_z, max_direct_z, max_local_z, local_range, pseudo_count, gen_prop, gap_counter):
        
        # Reset all parameters (Is needed, if deletion_caller is used several times in the same program)
        self.reset_global_parameters()
        
        # Set parameters
        self.bedgraph_path = bedgraph_path
        self.min_cld = int(min_cld)
        self.min_size = int(min_size)
        self.max_z = float(max_z)
        self.max_direct_z = float(max_direct_z)
        self.max_local_z = float(max_local_z)
        self.local_range = int(local_range)
        self.pseudo_count = int(pseudo_count)
        self.gen_prop = int(gen_prop)
        self.gaps = int(gap_counter)
        
        #### Call functions ####
        
        # First extract the coverage for all chromosomes
        self.extract_coverage()
        
        # Identify genome_starts and -ends for all chromosomes
        genome_starts, genome_ends = self.get_genome_positions()
        
        # Iterate over all chromosomes and do coverage analysis for each of them
        for chr in self.genome_coverage.keys():
            # Get chromosome specific genome_start and -end
            self.genome_start, self.genome_end = int(genome_starts[chr]), int(genome_ends[chr])
            # Assure that coverage counts are only in between the genome_start to genome_end range (remove others)
            self.genome_coverage[chr] = self.restrict_coverage(self.genome_coverage[chr], self.genome_start, self.genome_end)
            # Initialize empty deletions dict
            self.deletions[chr] = {}
            # Compute global mean and standard deviation of the whole coverage
            self.global_cov_mean, self.global_cov_std = self.compute_parameters(self.genome_coverage[chr])
            # Identify deletions by coverage analysis
            self.identify_deletions(chr)
            # Cluster neighbouring deletions together and display them as one
            self.cluster_deletions(chr)
        
    
    # This function is used to extract the coverage values (stored in a bedgraph file) and subsequently transfer them to a dictionary
    def extract_coverage(self):
        
        first_row = False
        last_end_pos = -1
        
        # Input bedgraph file
        bedgraph = pd.read_csv(self.bedgraph_path, header = None, sep="\t", names=["Chr", "Start", "End", "Cov"])
        
        for index, row in bedgraph.iterrows():
            # Extract all information from the bedgraph line
            chr = row.iloc[0]
            start_pos = row.iloc[1]
            end_pos = row.iloc[2]
            
            # Edge-case for first row to set last_end_pos
            if not first_row:
                last_end_pos = start_pos
                first_row = True
            
            # Initialize output-dictionary for corresponding chromosome
            if chr not in self.genome_coverage.keys():
                self.genome_coverage[chr] = {}
            
            # Extract coverage value
            coverage = row.iloc[3]
            
            # Positions with 0 coverage aren't listed in the bedgraph file. Include them in the output-dictionary with value = pseudo count!
            if start_pos != last_end_pos:
                for i in range(last_end_pos+1, start_pos+1):
                    self.genome_coverage[chr][i] = self.pseudo_count
            # If only "one" position is mentioned in the line
            if (end_pos - start_pos) == 1:
                self.genome_coverage[chr][end_pos] = coverage+self.pseudo_count
            # If multiple positions are bunched in one line (due to same consecutive coverage) extract them all individually
            else:
                for i in range(1, (end_pos+1)-start_pos):
                    self.genome_coverage[chr][start_pos+i] = coverage+self.pseudo_count
            last_end_pos = end_pos
    
    # Function to calculate z-Scores of a coverage value
    def compute_z_score(self, coverage, mean, std):
        return (coverage - mean)/std

    # Given the full or a subset of the coverage dictionary, compute the mean and std of the logarithmized coverage
    def compute_parameters(self, region_coverage):
        cov = np.array(list(region_coverage.values()))
        log_cov = np.log(cov)
        cov_mean = log_cov.mean()
        cov_std = log_cov.std()
        return cov_mean, cov_std

    # This function extracts a subset of the coverage dictionary and returns it as smaller coverage dictionary
    def get_local_coverage(self, position, max_distance, chr, detected_deletions):
        # Counts back propagated positions until max_distance is reached
        back_propagated_pos = 0
        # Is set to the previous position (start of backpropagation)
        tmp_pos = position-1
        # Subset/Output of coverage dictionary
        local_coverage = {}
        # Backpropagate either max_distance position or if we reach the genome start 
        while back_propagated_pos < max_distance and tmp_pos >= self.genome_start:
            # If current position isn't noted as deletion add this position and its coverage to the new dictionary and continue backpropagation
            if tmp_pos in detected_deletions.keys() and not detected_deletions[tmp_pos]:
                local_coverage[tmp_pos] = self.genome_coverage[chr][tmp_pos]
                back_propagated_pos += 1
                tmp_pos -= 1
            # If current position is noted as deletion don't add it to the dictionary and also don't count this position as backpropagated
            else:
                tmp_pos -= 1
        # Check if the function was able to backpropagate the desired length (can be restricted by the genome start)
        if back_propagated_pos == max_distance:
            return local_coverage
        # If not: Return empty dictionary (We want either the full or no information but not half stuff)
        else:
            return {}
    
    # Do coverage analysis and identify deletions with it
    def identify_deletions(self, chr):
        
        already_del = False
        local_mean = 0
        local_std = 0
        skip_position = False
        last_end = int(self.genome_start)
        detected_deletions = {}
                
        # Iterage over all genome positions
        # Idea: Update last_end position when no deletion, continue/skip updating when detecting a deletion
        # ==> When there is a deletion detected between two no-deletion regions, you should see the sudden jump of the last_end position
        for pos in self.genome_coverage[chr].keys():
        
            skip_position = False
            
            # Extract current coverage, log it and compute z-Score with global mean and std
            cov_current = self.genome_coverage[chr][pos]
            log_cov_current = np.log(cov_current)
            z = self.compute_z_score(log_cov_current, self.global_cov_mean, self.global_cov_std)
            # If z-Score is smaller (more significant) than the "direct"-threshold (more stringent threshold), directly mark this position as deletion
            if z < self.max_direct_z:
                detected_deletions[pos] = True
                continue
            # If z-Score is smaller (more significant) than the "overall"-threshold (less stringent) make further testings if position is deletion
            if z < self.max_z:
                # If there is no current deletion, extract a local coverage dict from the genome_coverage dict
                if not already_del:
                    local_cov = self.get_local_coverage(pos, self.local_range, chr, detected_deletions)
                    # If this local coverage dict was created successfully, compute a global mean and std of its logarithmized values
                    if len(local_cov) > 0:
                        local_mean, local_std = self.compute_parameters(local_cov)
                        # If local std is 0, z-Score can not be computed --> Skip this position
                        if local_std == 0:
                            skip_position = True
                    # If no coverage dict could be created, skip this position
                    else:
                        skip_position = True

                # If everything works as planned, compute local z-Score based on the local parameters and check if it is below the very stringent "local"-threshold
                # If a deletion is previous to this position and has not been interrupted yet (already_del = True), make use of the last defined local parameters, since new estimated parameters would include parts of the deletion and therefore skrew the analysis
                if not skip_position:
                    local_z = self.compute_z_score(log_cov_current, local_mean, local_std)
                    if local_z < self.max_local_z:
                        already_del = True
                        detected_deletions[pos] = True
                        continue
            # If no deletion was detected for the current position (coverage of this position passed no z-Score threshold)
            # Check if detected deletion is long enough   
            if (pos - last_end) >= self.min_size:
                start_del = last_end+1
                self.deletions[chr][start_del] = pos-1
            # If not unmark the corresponding position (so that they are not marke as deletions anymore)
            else:
                for i in range(last_end+1, pos):
                    detected_deletions[i] = False
            already_del = False
            last_end = pos
            detected_deletions[pos] = False
    
    # Cluster neighbouring deletions
    def cluster_deletions(self, chr):
        
        self.deletion_clusters[chr] = {}
        
        cluster_starts = []
        cluster_ends = []

        first_cluster = True

        # Iterate over detected deletions
        for del_pos in self.deletions[chr].keys():
            
            # Get start- and end position of deletion
            start_pos = int(del_pos)
            end_position = int(self.deletions[chr][del_pos])
            
            # Idea: Add start- and end positions of deletions in lists and check if there are any other deletions, whose start is closely to the end of the current cluster
            # Edge-case for first cluster
            if first_cluster:
                cluster_starts.append(start_pos)
                cluster_ends.append(end_position)
                first_cluster = False
            # Check if start position of new deletion is within a defined range (minimum_cluster_differenc/min_cld)
            # If so append its positions to the corresponding clusters
            elif start_pos <= (cluster_ends[-1] + self.min_cld):
                cluster_starts.append(start_pos)
                cluster_ends.append(end_position)
            # If deletion can not be added to current cluster, save current cluster and "open" a new cluster for this deletion
            else:
                self.deletion_clusters[chr][cluster_starts[0]] = cluster_ends[-1]
    
                # Reset clusters
                cluster_starts = [start_pos]
                cluster_ends = [end_position]

        # Adding last cluster
        if len(cluster_starts) > 0:
            self.deletion_clusters[chr][cluster_starts[0]] = cluster_ends[-1]
    
    #### Getter-functions and file-writing functions ####
    
    def get_deletion_clusters(self):
        return self.deletion_clusters

    def write_deletion_clusters_to_file(self, file_path):
        with open(file_path, "w") as file:
            file.write("CHR\tSTART\tEND\n")
            for chr in self.deletion_clusters.keys():
                for del_start in self.deletion_clusters[chr].keys():
                    file.write(f"{chr}\t{del_start}\t{self.deletion_clusters[chr][del_start]}\n")
        file.close()
    
    def write_deletions_to_file(self, file_path, deletions):
        with open(file_path, "w") as file:
            file.write("# DELETIONS\n")
            file.write("CHR\tPOSITION\tCLIPPING POSITIONS\n")
            for chr in deletions.keys():
                for deletion in deletions[chr]:
                    file.write(f"{chr}\t{deletion[0]}-{deletion[1]}\n")

    # This function computes the genome start and -end, since outliers can skrew these positions and with that the coverage analysis
    def get_genome_positions(self):
        start_positions = {}
        end_positions = {}
        # Compute start and end position for echa chromosome individually
        for chr in self.genome_coverage.keys():
            positions = list(self.genome_coverage[chr].keys())
            
            # Stores the current guess for the genome start position
            tmp_start = positions[0]
            # Stores the current position (for each iteration)
            tmp_pos_forward = tmp_start
            # Stores the number of propagated positions
            propagated_forward_pos = self.gen_prop
            # Counts gaps within our iteration
            gap_counter = 0
            # Is needed if start position is discarded
            no_new_start = False
            
            
            # Idea: Propagate over a predefinded number of positions
            while propagated_forward_pos > 0:
                
                # If there is a deletion/0-Coverage value
                if self.genome_coverage[chr][tmp_pos_forward] == self.pseudo_count:
                    # Check if our gap counter is still below a predefined number of max gaps
                    if gap_counter < self.gaps:
                        # Increase gap counter and propagate forward
                        gap_counter += 1
                        tmp_pos_forward += 1
                        propagated_forward_pos -= 1
                    # If there were already to many gaps, iterate forward (don't increase propagation counter) and discard current start
                    else:
                        tmp_pos_forward += 1
                        no_new_start = True
                    continue
                
                # If start was discarded because of too many gaps, define new start and reset the propagation- and gap counter
                # Also set the current position as new genome start prediction/guess
                if no_new_start:
                    tmp_start = tmp_pos_forward
                    tmp_pos_forward += 1
                    propagated_forward_pos = self.gen_prop-1
                    gap_counter = 0
                    no_new_start = False
                # If start wasn't discarded and the current position has a Coverage > 0, then reset gap counter (gaps are counted only consecutive) and propagate forward
                else:
                    gap_counter = 0
                    tmp_pos_forward += 1
                    propagated_forward_pos -= 1
            
            # Get last/best guess for genome start position and store it
            start_positions[chr] = tmp_start
            
            # Simply do the same stuff for the genome end (just reverse it and don't calculate with +1 but with -1)
            
            tmp_end = positions[-1]
            tmp_pos_reverse = tmp_end
            propagated_reverse_pos = self.gen_prop
            gap_counter = 0
            no_new_end = False
                        
            while propagated_reverse_pos > 0:
                if self.genome_coverage[chr][tmp_pos_reverse] == self.pseudo_count:
                    if gap_counter < self.gaps:
                        gap_counter += 1
                        tmp_pos_reverse -= 1
                        propagated_reverse_pos -= 1
                    else:
                        tmp_pos_reverse -= 1
                        no_new_end = True
                    continue
                    
                if no_new_end:
                    tmp_end = tmp_pos_reverse
                    tmp_pos_reverse -= 1
                    propagated_reverse_pos = self.gen_prop-1
                    gap_counter = 0
                    no_new_end = False
                else:
                    gap_counter = 0
                    tmp_pos_reverse -= 1
                    propagated_reverse_pos -= 1
            
            end_positions[chr] = tmp_end

        return start_positions, end_positions
    
    
    # Extract only the part of the genome coverage which lies inbetween the genome start and -end
    def restrict_coverage(self, coverage, genome_start, genome_end):
        new_coverage = {pos: cov for pos, cov in coverage.items() if genome_start <= pos <= genome_end}
        return new_coverage
        
    
    def reset_global_parameters(self):
        self.bedgraph_path = ""
        self.output_path = ""
        self.min_cld = -1
        self.min_size = -1
        self.max_z = -1.0
        self.max_direct_z = -1.0
        self.max_local_z = -1.0
        self.local_range = -1
        self.genome_start = -1
        self.genome_end = -1
        self.pseudo_count = -1

        self.deletions = {}
        self.genome_coverage = {}
        self.last_end_pos = -1

        self.global_cov_mean = -1.0
        self.global_cov_std = -1.0

        self.already_deletions = False

        self.local_cov_mean = -1.0
        self.local_cov_std = -1.0

        self.skip_position = False

        self.deletion_clusters = {}
        
        self.gen_prop = -1
        self.gaps = -1

# insertion_caller calls insertions and verifies deletions by clipping pattern analysis 
class insertion_caller:
    
    # Ininitialize all global variables
    bam_path = ""
    output_path = ""
    max_patt_diff = -1
    min_z = -1.0
    min_sur_z = -1.0
    max_hit = -1
    window_size = -1
    tolerance = -1.0
    genome_start = {}
    genome_end = {}
    pseudo_count = -1
    min_reads = -1
    
    genome_left = {}
    genome_right = {}
    
    z_left = {}
    z_right = {}
    
    chromosomes = []
    deletions = {}
    
    verified_deletions = {}
    insertions = {}
    
    deletions_consensus = {}
    insertions_consensus = {}
    
    reference = {}
    
    new_insertions = {}
    
    mpc = -1.0
    min_length = -1
    clp_ver_range = -1
    
    
    def __init__(self, bam_path, max_patt_diff, min_z, min_sur_z, window_size, tolerance, genome_start, genome_end, pseudo_count, min_reads, deletions, reference_path, fir_ws, sec_ws, con_path, mpc, min_length, clp_ver_range):
        
        # Reset all parameters (Is needed, if deletion_caller is used several times in the same program)
        self.reset_global_parameters()
        
        # Set parameters
        self.bam_path = bam_path
        self.max_patt_diff = int(max_patt_diff)
        self.min_z = float(min_z)
        self.min_sur_z = float(min_sur_z)
        self.window_size = int(window_size)
        self.tolerance = float(tolerance)
        self.genome_start = genome_start
        self.genome_end = genome_end
        self.pseudo_count = int(pseudo_count)
        self.min_reads = int(min_reads)
        # Deletions/Clustered deletions from the coverage analysis (-> deletion_caller)
        self.deletions = deletions
        self.mpc = float(mpc)
        self.min_length = int(min_length)
        self.clp_ver_range = int(clp_ver_range)
        
        self.extract_reference_seq(reference_path)
        
        # Count clipped reads
        self.compute_clippings()
        
        # Do clipping pattern analysis and verification for each chromosome individually
        for chr in self.chromosomes:
                                    
            # Skip if a chromosome was not detected during coverage analysis
            if chr not in self.deletions.keys():
                continue
            
            # Main outputs of insertion_caller
            self.verified_deletions[chr] = []
            self.insertions[chr] = []
            # Contains deletions that were not verified by clipped sequences and thus may contain inserted sequences
            self.new_insertions[chr] = []
            
            # Insertion calling and deletion verification
            self.run_genome_iteration(chr)
            
            # If at least one deletion was found:
            if len(self.verified_deletions[chr]) > 0:
                # Iterate over all deletions to verify them, if their clipped sequences match the oposite site of the reference genome
                for ver_del in self.verified_deletions[chr]:
                    # Get position weight matrices for the deletion (consensus sequences are not used here)
                    start_pwm, end_pwm, start_cons, end_cons = self.extract_consensus_seq(chr, ver_del[2]-1, ver_del[3]+1)

                    # Compute score for both clipped sequences/regions
                    score_start, score_end = self.compare_cons_to_ref(ver_del[2]-1, ver_del[3]+1, start_pwm, end_pwm, chr)
                    #print(score_start, score_end, start_cons, end_cons)
                    # If both clipped sequences fail to pass threshold append a warning suggesting an insertion or wrong prediction
                    if score_start < fir_ws and score_end < fir_ws:
                        ver_del.append("== Warning: The deletion might be wrong or an insertion was placed into this deletion ==")
                        self.new_insertions[chr].append([ver_del[2]-1, ver_del[3]+1])
                    # If one clipped sequence fails to pass threshold, check if the other score passes a more stringent threshold -> Append warning
                    elif score_start < fir_ws:
                        if score_end < sec_ws:
                            ver_del.append("== Warning: Likely Upstream-mismatch & weak Downstream-matches of clipped sequences ==")
                        else:
                            ver_del.append("== Warning: Likely Upstream-mismatch of clipped sequences ==")
                    elif score_end < fir_ws:
                        if score_start < sec_ws:
                            ver_del.append(ver_del.append("== Warning: Weak Upstream-matches & likely Downstream-mismatches of clipped sequences =="))
                        else:
                            ver_del.append("== Warning: Likely Downstream-mismatch of clipped sequences ==")
                    # If the scores of both clipped sequences pass threshold, append 0, suggesting no warning and all is fine
                    else:
                        ver_del.append(0)
            
            # Extract consensus sequences of the clipped sequence parts for insertions, which are used as input later on to identify the inserted sequence
            insertions_to_write = []
            if len(self.insertions[chr]) > 0:
                for insertion in self.insertions[chr]:
                    
                    # We have to look, what type of insertion/overlap we've got
                    second_ins_pos = -1
                    # If clipped reads "collide", the overlap is 0 but the second position is actually +1
                    if insertion[1] == 0:
                        second_ins_pos = insertion[0]+1
                        start_pwm, end_pwm, start_cons, end_cons = self.extract_consensus_seq(chr, insertion[0], second_ins_pos)
                    # If clipped reads overlap at least by one position, second position is +overlap-1
                    else:
                        second_ins_pos = insertion[0]+insertion[1]-1
                        # second position has to be at the first place, since it contains the right-clipped reads
                        start_pwm, end_pwm, start_cons, end_cons = self.extract_consensus_seq(chr, second_ins_pos, insertion[0])
                    
                    insertions_to_write.append([insertion[0], second_ins_pos, start_cons, end_cons])
                    
                # Get consensus sequences of deletions to verify if they may be insertions
                for new_insertion in self.new_insertions[chr]:
                    start_pwm, end_pwm, start_cons, end_cons = self.extract_consensus_seq(chr, new_insertion[0], new_insertion[1])
                    insertions_to_write.append([new_insertion[0], new_insertion[1], start_cons, end_cons])                
                
                self.write_consensus_to_file(con_path, insertions_to_write)
    
    def extract_reference_seq(self, ref_path):
        with open(ref_path, "r") as file:
            chr = ""
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    chr = line[1:]
                    self.reference[chr] = ""
                else:
                    self.reference[chr] += line
        file.close()
    
    # Calculate z-Score     
    def compute_z_score(self, input, mean, std):
        return (input - mean)/std

    # Identifies and counts clipped reads
    def compute_clippings(self):
                
        # Open bamfile and extract chromosomes
        bamfile = pysam.AlignmentFile(self.bam_path, "rb")
        self.chromosomes = list(bamfile.references)
        
        # Initialize for each chromosomes the dictionaries to count left clipped reads (genome_left) and right clipped reads (genome_right)
        for chr in self.chromosomes:
            self.genome_left[chr] = {}
            self.genome_right[chr] = {}
            for i in range(self.genome_start[chr], self.genome_end[chr]+1):
                self.genome_left[chr][i] = 0
                self.genome_right[chr][i] = 0

        # Iterate over all reads in bamfile
        for read in bamfile.fetch():
            
            # Extract important information
            read_start = int(read.reference_start)+1
            read_cigar = read.cigarstring
            read_chr = read.reference_name
            
            # Get number of matched positions (matched_positions)
            match = re.search(r'(\d+)M', read_cigar)
            if match:
                matched_positions = int(match.group(1))
            else:
                continue
            
            # Calculate (matched, so no clippings) read_end (not specified in bamfile) with matched_positions
            read_end = (read_start + matched_positions)-1
            
            # Ignore all reads not in the specified genome range
            if read_start < self.genome_start[chr] or read_end > self.genome_end[chr]:
                continue
            
            clipping = "NA"
            
            # Split cigar-String by matched positions (M) to identify clipping direction of read
            tmp_splitted = read_cigar.split("M")
            # If "S" is on the left, but not on the right site of "M": Left clipped read
            if "S" in tmp_splitted[0] and not "S" in tmp_splitted[1]:
                clipping = "L"
            # If "S" is noit on the left, but on the right site of "M": Right clipped read
            elif "S" in tmp_splitted[1] and not "S" in tmp_splitted[0]:
                clipping = "R"
            # If "S" is on the left and on the right site of "M": Left and Right clipped read ("B" := "Both")
            else:
                clipping = "B"
            
            # If read is left clipped, increase counter for left clipped reads at the position of the read start
            if clipping == "L":
                self.genome_left[read_chr][read_start] += 1
            # If read is right clipped, increase counter for right clipped reads at the position of the read end
            elif clipping == "R":
                self.genome_right[read_chr][read_end] += 1
                #if read_end == 120910:
                    #test.append(read)
            # Increase both
            elif clipping == "B":
                self.genome_left[read_chr][read_start] += 1
                self.genome_right[read_chr][read_end] += 1
                #if read_end == 120910:
                    #test.append(read)
            # If no clipping was recognized, just skip the read
            else:
                continue
        # Compute for each chromosomes the z-Scores
        for chr in self.chromosomes:
            
            # Logarithmize the total left and right reads speparately
            log_left_clips = np.array(list(self.genome_left[chr].values()))
            log_right_clips = np.array(list(self.genome_right[chr].values()))

            # Compute mean and standard deviation of logarithmized reads for left and right reads speparately
            mean_left_clips = np.mean(log_left_clips)
            std_left_clips = np.std(log_left_clips)
            mean_right_clips = np.mean(log_right_clips)
            std_right_clips = np.std(log_right_clips)
            
            self.z_left[chr] = {}
            self.z_right[chr] = {}
            
            # Compute z-Score for each position with the parameters from above
            for pos in range(self.genome_start[chr]+1, self.genome_end[chr]+1):
                self.z_left[chr][pos] = self.compute_z_score(self.genome_left[chr][pos], mean_left_clips, std_left_clips)
                self.z_right[chr][pos] = self.compute_z_score(self.genome_right[chr][pos], mean_right_clips, std_right_clips)
        
        
    # Function that can identify for a specific position if there is a significant peak
    def check_peak(self, pos, z_dict, chr):
        
        # If position has too few right- and left clipped reads, skip/return False
        if self.genome_left[chr][pos] < self.min_reads and self.genome_right[chr][pos] < self.min_reads:
            return False
        
        # Get z-Score for position (if left or right peak is determined by the z_dict parameter in the call of the function -> If you want e.g. check a left peak, provide z_left in function call)
        z_score = z_dict[pos]
        
        # If z-Score is larger (more significant) than a threshold, make further tests
        if z_score > self.min_z:
            
            tmp_left = []
            tmp_right = []
            local_z_left = -1.0
            local_z_right = -1.0
            
            # Check if peak is also locally significant
            # Iterate over a predefined range/window around position and add the peaks/number of clipped reads to the corresponding lists
            for j in range(pos-(self.window_size), pos+(self.window_size+1)):
                
                # Don't add peak of interest to list, otherwise it may skrew analysis
                if j == pos or j not in z_dict.keys():
                    continue
                
                tmp_left.append(self.genome_left[chr][j])
                tmp_right.append(self.genome_right[chr][j])
            
            # Compute local parameters of these lists
            local_mean_left = np.mean(tmp_left)
            local_std_left = np.mean(tmp_left)
            local_mean_right = np.mean(tmp_right)
            local_std_right = np.mean(tmp_right)
            
            # Compute both z-Scores
            if local_std_left != 0:
                local_z_left = self.compute_z_score(self.genome_left[chr][pos], local_mean_left, local_std_left)
            if local_std_right != 0:
                local_z_right = self.compute_z_score(self.genome_right[chr][pos], local_mean_right, local_std_right)
            
            # Check if both z-Scores fail a very stringent threshold
            if local_z_left < self.min_sur_z and local_z_right < self.min_sur_z:
                return False
            # If not, call the position a peak
            return True
        return False
    
    # Function that detects insertions and verifies deletions based on clipping pattern analysis
    def run_genome_iteration(self, chr):
        
        last_left_peak = -1
        last_right_peak = -1

        z_left, z_right = self.z_left[chr], self.z_right[chr]
        # Store insertions/deletions that were already verified -> Are used to compare similar insertions/deletions and decide which one is the best
        already_insertions = []
        already_verified = []
        
        # Idea: Iterate over genome, mark peaks and check the distance and directions (left clipped- or right clipped peak) to determine if insertion or deletion was found
        for i in range(self.genome_start[chr]+1, self.genome_end[chr]+1):
            # Check if we have left or/and right peak for a specific position
            is_left_peak = self.check_peak(i, z_left, chr)
            is_right_peak = self.check_peak(i, z_right, chr)
            
            # If we have two peaks at the same position, then this is a insertion with overlap = 1 --> Skip afterwards
            if is_left_peak and is_right_peak:
                # Set last detected left- and right peak
                last_left_peak = i
                last_right_peak = i
                # Add insertion
                self.insertions[chr].append([i, 1])
                already_insertions.append(i)
                continue
                
            # Deletions arise, when there is first a peak of right clipped reads followed by a peak of left clipped reads
            # ==> If we have at the current position a left clipped peak: Look back to find a right clipped peak --> Deletion
            if is_left_peak:
                # Set new last left peak
                last_left_peak = i
                
                # Start case
                if last_right_peak == -1:
                    continue
                
                # Actually an insertion can also be found with this clipping pattern, but only with overlap = 0 (when the reads "collide")
                if i-last_right_peak == 1:
                    self.insertions[chr].append([i-1, 0])
                    already_insertions.append(i-1)
                    continue
                
                # Define deletion start as the previously/last found right peak +1
                pot_deletion_start = last_right_peak+1
                # Define deletion end at the current detected left peak -1
                pot_deletion_end = i-1
                
                # Iterate over all detected deletions to verify the potential deletion:
                for del_start in self.deletions[chr].keys():
                    del_end = self.deletions[chr][del_start]
                    
                    # Calculate the tolerance (number of allowed deviating positions) based on the length of the deletions (from coverage analysis)
                    tolerance_value = int(((del_end-del_start)+1)*self.tolerance)
                    if tolerance_value == 0:
                        tolerance_value = 1
                    
                    # If deletion found by clipping pattern is within the tolerance window of a deletion found by coverage analysis
                    if abs(del_start-pot_deletion_start) <= tolerance_value and abs(del_end-pot_deletion_end) <= tolerance_value:
                        
                        # Now only check left to do is, if this deletion from the coverage analysis was already identified by another clipping pattern
                        # If that's the case:
                        if del_start in already_verified:
                            tmp_deletions = self.verified_deletions[chr]
                            for deletion_pair in tmp_deletions:
                                # Get the corresonding clipping pattern deletion
                                if deletion_pair[0] == del_start:
                                    # Compute for both (current clipping pattern and the previous detected clipping pattern) the sum of the z-Scores
                                    z_sum_current = z_left[i] + z_right[last_right_peak]
                                    z_sum_verified = z_left[deletion_pair[3]+1] + z_right[deletion_pair[2]-1]
                                    # Check if z-Score sum of current/new clipping pattern is better/greater than the previously found one
                                    if z_sum_current > z_sum_verified:
                                        # If so: remove previously found clipping pattern and add new one
                                        self.verified_deletions[chr].remove(deletion_pair)
                                        self.verified_deletions[chr].append([del_start, del_end, pot_deletion_start, pot_deletion_end])
                                    break
                        # If deletion from coverage analysis has not been matched to a clipping pattern deletion yet, then do it now
                        else:
                            self.verified_deletions[chr].append([del_start, del_end, pot_deletion_start, pot_deletion_end])
                            already_verified.append(del_start)
            
            # Insertions arise, when there is first a peak of left clipped reads followed by a peak of right clipped reads
            # ==> If we have at the current position a rigth clipped peak: Look back to find a left clipped peak --> Insertion/Overlap
            if is_right_peak:
                # Set current position as new last right peak
                last_right_peak = i
                
                # Start case
                if last_left_peak == -1:
                    continue
                
                # If the overlap is longer than a predefined threshold (too long overlaps make no sense) skip this peak
                if i-last_left_peak+1 > self.max_patt_diff:
                    continue
                
                # Do exactly the same stuff for insertions to check whether there is a better clipping pattern for the "same" insertion (start)
                if last_left_peak in already_insertions:
                    tmp_insertions = self.insertions[chr]
                    for insertion in tmp_insertions:
                        if insertion[0] == last_left_peak:
                            z_sum_current = z_left[last_left_peak] + z_right[i]
                            if insertion[1] == 0:
                                z_sum_insertion = z_right[insertion[0]] + z_left[insertion[0] + 1]
                            else:
                                z_sum_insertion = z_left[insertion[0]] + z_right[insertion[0] + (insertion[1]-1)]
                            if z_sum_current > z_sum_insertion:
                                self.insertions[chr].remove(insertion)
                                self.insertions[chr].append([last_left_peak, (i-last_left_peak)+1])
                            break
                # If no clipping patter was already found for this insertion, add it
                else:
                    self.insertions[chr].append([last_left_peak, (i-last_left_peak)+1])
                    already_insertions.append(last_left_peak)
    
    
    def extract_consensus_seq(self, chr, clp_start_pos, clp_end_pos):
        
        # Define the 4 nucleotides
        nucl = ["A", "C", "G", "T"]
                
        nucl_abund = self.get_nucleotide_abundancies(chr)
        
        # Function that computes the read end based on the cigar string and the read start
        def get_read_end(cigar, read_start):
            
            match = re.search(r'(\d+)M', cigar)

            if match:
                matched_positions = int(match.group(1))
            else:
                return -1, -1
            
            read_end = (read_start + matched_positions)-1
            return read_end, matched_positions
        
        # Get the length of the longest clipped sequence
        def get_longest_clip(seq_list):
            longest_clip = -1
            for clip in seq_list:
                length_clip = len(clip)
                if length_clip > longest_clip:
                    longest_clip = length_clip
            return longest_clip-1

        # Create Position-Weight-Matrix based on a list of clipped sequences
        def get_clipped_sequences(seq_list):

            longest_clip = get_longest_clip(seq_list)
            count_nucl = {}
            # Iterate over all position (from 0 to the end of the longest clipped sequence)
            for i in range(0, longest_clip):
                count_nucl[i] = {"A": self.mpc, "C": self.mpc, "G": self.mpc, "T": self.mpc, "Total": 4*self.mpc}
                # Iterate over all clipped sequences
                for clip in seq_list:
                    if i < len(clip) and clip[i] != "N":
                        # Increase total number of clipped sequences overlapping this position and of course the number of the corresponding nucleotide as well
                        count_nucl[i]["Total"] += 1
                        count_nucl[i][clip[i]] += 1
            return count_nucl
        
        def get_consensus_seq(counts):
            
            # Define nucleotides
            nucl = ["A", "C", "G", "T"]
            
            # Initialize output consensus sequence
            pred_consensus = ""
            # Initialize output PWM
            pwm = pd.DataFrame(columns=nucl)
            # Iterate over all positions of the PWM
            for pos in counts.keys():
                num_total = float(counts[pos]["Total"])
                # If total number of reads is below a threshold the consensus nucleotide can't be predicted with high confidence --> STOP!
                if num_total < 30:
                    break
                # Compute for each position the percentage of the nucleotides
                perc_A = float(counts[pos]["A"])/num_total
                perc_C = float(counts[pos]["C"])/num_total
                perc_G = float(counts[pos]["G"])/num_total
                perc_T = float(counts[pos]["T"])/num_total
                
                perc_A_norm = np.log2(perc_A) - np.log2(nucl_abund[0])
                perc_C_norm = np.log2(perc_C) - np.log2(nucl_abund[1])
                perc_G_norm = np.log2(perc_G) - np.log2(nucl_abund[2])
                perc_T_norm = np.log2(perc_T) - np.log2(nucl_abund[3])
                
                new_row = pd.DataFrame({"A": [perc_A_norm], "C": [perc_C_norm], "G": [perc_G_norm], "T": [perc_T_norm]})
                # Add new row to PWM
                pwm = pd.concat([pwm, new_row], ignore_index=True)
                
                # Initialize them in a list to calculate the nucleotide with the max percentage
                perc = [perc_A, perc_C, perc_G, perc_T]
                max_nucl_idx = perc.index(max(perc))
                max_nucl = nucl[max_nucl_idx]
                #max_perc = perc[max_nucl_idx]
                pred_consensus += max_nucl
                

            return pwm, pred_consensus
        

        # Get the position of peaks
        
        bamfile = pysam.AlignmentFile(self.bam_path, "rb")

        # Contains all clipped sequences of right clipped reads
        seq_list_start = []
        # Contains all clipped sequences of elft clipped reads
        seq_list_end = []
        
        # Iterate over all reads in the corresponding region (Now: right clipped reads)
        for read in bamfile.fetch(chr, clp_start_pos-1, clp_start_pos+1):
            
            # Get important parameters
            read_start = int(read.reference_start)+1
            read_cigar = read.cigarstring
            read_seq = read.query_sequence
            
            # Calculate read end and matched positions
            read_end, matched_positions = get_read_end(read_cigar, read_start)
            
            if read_end == -1:
                continue
            
            # Only take reads that are exactly at the peak position and have a clipped sequence
            if read_end != clp_start_pos or matched_positions == len(read_seq):
                continue
            # If a read has its end at the peak position, add its right clipped sequence to the corresponding list
            if read_end == clp_start_pos:
                seq_list_start.append(read_seq[matched_positions:])
        
        # Same stuff just for left clipped reads (apart of: see down below)
        for read in bamfile.fetch(chr, clp_end_pos-1, clp_end_pos+1):
            
            read_start = int(read.reference_start)+1
            read_cigar = read.cigarstring
            read_seq = read.query_sequence
            
            read_end, matched_positions = get_read_end(read_cigar, read_start)
            
            if read_end == -1:
                continue
            
            if read_start != clp_end_pos or matched_positions == len(read_seq):
                continue
            if read_start == clp_end_pos:
                # Get clipped sequence
                clipped_seq = read_seq[:len(read_seq)-(read_end-read_start)-1]
                # Reverse it, since computation is easier then
                rev_seq = clipped_seq[::-1]
                seq_list_end.append(rev_seq)

        # Get PWMs
        count_nucl_start = get_clipped_sequences(seq_list_start)
        count_nucl_end = get_clipped_sequences(seq_list_end)  

        # Get consensus sequences
        start_pwm, start_cons = get_consensus_seq(count_nucl_start)
        end_pwm, end_cons = get_consensus_seq(count_nucl_end)
        # Reverse consensus sequence
        end_cons = end_cons[::-1]
        # Clipped sequences were reversed for easier computation -> Reverse rows of PWM
        end_pwm = end_pwm.iloc[::-1].reset_index(drop=True)
        
        return start_pwm, end_pwm, start_cons, end_cons
            

    def compare_cons_to_ref(self, clp_start, clp_end, pwm_start, pwm_end, chr):
                
        # Extract sliding window arround clipped positions
        start_ref_string = self.reference[chr][clp_start-self.clp_ver_range:clp_start+self.clp_ver_range]
        end_ref_string = self.reference[chr][clp_end-self.clp_ver_range:clp_end+self.clp_ver_range]
        
        # Get length of consensus sequences
        len_cons_start = pwm_start.shape[0]
        len_cons_end = pwm_end.shape[0]
 
        best_kmer_score_start = -100000000000
        #best_kmer_start = -1
        
        # Iterate over all len_cons_end-mers in the window
        for i in range(0, len(start_ref_string)-len_cons_end+1):
            
            kmer = start_ref_string[i:i+len_cons_end]
            N_counter = 0
            kmer_score = 0
            
            # Calculate kmer-Score
            for index, nucl in enumerate(kmer):
                if nucl == "N":
                    N_counter += 1
                    continue
                kmer_score += pwm_end.loc[index, nucl]

            # Normalize kmer-Score by kmer length
            kmer_score_norm = kmer_score/(len(kmer)-N_counter)
            
            # Update best kmer-Score if necessary
            if kmer_score_norm > best_kmer_score_start:
                best_kmer_score_start = kmer_score_norm
                #best_kmer_start = clp_start-100+i
        
        # Same stuff, just for end of the deletion
        best_kmer_score_end = -100000000000
        
        # Iterate over all len_cons_end-mers in the window
        for i in range(0, len(end_ref_string)-len_cons_start+1):
            
            kmer = end_ref_string[i:i+len_cons_start]
            N_counter = 0
            kmer_score = 0
            
            # Calculate kmer-Score
            for index, nucl in enumerate(kmer):
                if nucl == "N":
                    N_counter += 1
                    continue
                kmer_score += pwm_start.loc[index, nucl]
            
            # Normalize kmer-Score by kmer length 
            kmer_score_norm = kmer_score/(len(kmer)-N_counter)
            
            # Update best kmer-Score if necessary
            if kmer_score_norm > best_kmer_score_end:
                best_kmer_score_end = kmer_score_norm
                #best_kmer_start = clp_start-100+i
        
        return best_kmer_score_start, best_kmer_score_end   
    
    # This function writes the consensus sequence to a file, so they can be used for further analysis
    def write_consensus_to_file(self, file_path, insertions):
        
        with open(file_path, "w") as file:
            for insertion in insertions:
                # Consensus will only be written to the output file, if they are long (and with that significant) enough
                if len(insertion[2]) > self.min_length:
                    file.write(f">{insertion[0]}_{insertion[1]}_START\n")
                    file.write(f"{insertion[2]}\n")
                if len(insertion[3]) > self.min_length:
                    file.write(f">{insertion[0]}_{insertion[1]}_END\n")
                    file.write(f"{insertion[3]}\n")
        file.close()
    
    # Get nucleotide abundancies of a chromosome
    def get_nucleotide_abundancies(self, chr):
        # Count the nucleotides
        num_A = self.reference[chr].count("A")
        num_C = self.reference[chr].count("C")
        num_G = self.reference[chr].count("G")
        num_T = self.reference[chr].count("T")
        num_N = self.reference[chr].count("N")
        
        # Get total number of nucleotides (substract number of N, so that the probability of A, C, G and T sums up to one)
        num_total = len(self.reference[chr]) - num_N
        
        # Return percentages and add up a small pseudo count, so that even 0 probabilities can be logged
        return [float(num_A)/float(num_total), float(num_C)/float(num_total), float(num_G)/float(num_total), float(num_T)/float(num_total)]
    
    #### Getter and writing functions ####
    
    def get_clippings(self, chr):
        return self.genome_left[chr], self.genome_right[chr]
    
    def get_z_scores(self, chr):
        return self.z_left[chr], self.z_right[chr]

    def get_chromosomes(self):
        return self.chromosomes
    
    def write_clippings_to_file(self, file_path):
        with open(file_path, "w") as file:
            file.write("Chr\tPosition\tLeft-clipped Reads\tRight-clipped Reads\n")
            for chr in self.genome_left.keys():
                for pos in self.genome_left[chr].keys():
                    file.write(f"{chr}\t{pos}\t{self.genome_left[chr][pos]}\t{self.genome_right[chr][pos]}\n")
        file.close()
    
    def write_z_scores_to_file(self, file_path, chr):
        with open(file_path, "w") as file:
            file.write("Position\tLeft z-Score\tRight z-Score\n")
            for pos in self.z_left[chr].keys():
                file.write(f"{pos}\t{self.z_left[chr][pos]}\t{self.z_right[chr][pos]}\n")
        file.close()
    
    def write_insertions_to_file(self, file_path, insertions):
        with open(file_path, "w") as file:
            file.write("# INSERTIONS\n")
            file.write("CHR\tPOSITION\tOVERLAP\n")
            for chr in insertions.keys():
                for insertion in insertions[chr]:
                    file.write(f"{chr}\t{insertion[0]}\t{insertion[1]}\n")
        file.close()
     
    def write_deletions_to_file(self, file_path, deletions):
        with open(file_path, "w") as file:
            file.write("# DELETIONS\n")
            file.write("CHR\tSTART\tEND\n")
            for chr in deletions.keys():
                for deletion in deletions[chr]:
                    # No warning for deletion
                    if deletion[4] == 0:
                        file.write(f"{chr}\t{deletion[0]}\t{deletion[1]}\n")   
                    # Warning for deletion
                    else:
                        file.write(f"{deletion[4]}\n")
                        file.write(f"{chr}\t{deletion[0]}\t{deletion[1]}\n")   

    def get_insertions_and_verified_deletions(self):
        return self.insertions, self.verified_deletions
    
    def reset_global_parameters(self):
        self.bam_path = ""
        self.output_path = ""
        self.max_patt_diff = -1
        self.min_z = -1.0
        self.min_sur_z = -1.0
        self.window_size = -1
        self.tolerance = -1.0
        self.genome_start = {}
        self.genome_end = {}
        self.pseudo_count = -1
        self.min_reads = -1

        self.genome_left = {}
        self.genome_right = {}

        self.z_left = {}
        self.z_right = {}

        self.chromosomes = []
        self.deletions = {}

        self.verified_deletions = {}
        self.insertions = {}
    
        self.deletions_consensus = {}
        self.insertions_consensus = {}

        self.reference = {}

        self.new_insertions = {}

        self.mpc = -1.0
        self.min_length = -1
        self.clp_ver_range = -1
    