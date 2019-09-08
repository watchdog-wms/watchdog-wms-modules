'''
Created on 22 Nov 2017

@author: friedl
'''

import subprocess
import os.path
import sys
import collections


def build_STAR_index(genome_file, target_folder, star_path ='STAR', star_thread_nr=1, overhang=100, annotation_file=None):
    '''
    calls STAR to generate a genome index
    'genome_file': multi-fasta file with the chromosome sequences
    'target_folder': output folder for writing the STAR index
    'star_path': path to the STAR executable (default: in PATH)
    'star_thread_nr': number of threads to use for STAR
    'annotation_file': include gene annotations as GTF file into the index
    'overhang': only required for the annotation file
    '''
    
    # indexing command without annotations
    index_command = [star_path, '--runThreadN', str(star_thread_nr), '--runMode', 'genomeGenerate', '--genomeDir', target_folder,
                     '--genomeFastaFiles', genome_file]
    
    # indexing command with gene annotations
    if(annotation_file is not None):
        index_command.append(['--sjdbGTFfile', annotation_file, '--sjdbOverhang', str(overhang)])
    
    # execute indexing command as child process
    print('Generate STAR index: '+' '.join(index_command))
    try:
        # set working directory to output directory -> log file is written to working directory
        subprocess.check_call(index_command,cwd=target_folder)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in STAR genomeGenerate call \n'+' '.join(index_command))
        raise

def map_with_STAR(reads1_file, reads2_file, index_folder, out_prefix, star_path='STAR', star_thread_nr=1):
    '''
    maps paired-end reads with star
    'reads1_file' & 'reads2_file': fastq files with the first and second reads of each pair (no singletons allowed) for paired end-data
        'reads2_file' = None for single-end data
    'index_folder': folder with the STAR index of the reference genome
    'out_prefix': prefix for all STAR output files
    'star_path': path to STAR executable, default in PATH
    'star_thread_nr': number of threads to use with STAR
    '''
    
    # parameters passed by the user
    mapping_command=[star_path, '--genomeDir', index_folder, '--runThreadN', str(star_thread_nr),
                     '--outFileNamePrefix', out_prefix, '--readFilesIn', reads1_file ]
    # single-end data: file2 = None, paired-end data: append file2 for to the '--readFilesIn' argument 
    if reads2_file is not None:
        mapping_command.append(reads2_file)
    
    # fixed parameters used for circular reads -> taken from circRNA_finder wrapper script
    mapping_command.extend(['--chimSegmentMin', '20', '--chimScoreMin', '1', '--alignIntronMax', '100000', '--outFilterMismatchNmax', '4',
                            '--alignTranscriptsPerReadNmax', '100000', '--outFilterMultimapNmax', '2'])
    
    # execute STAR as child process
    print('Mapping with STAR: '+' '.join(mapping_command))
    try:
        subprocess.check_call(mapping_command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in STAR genomeGenerate call \n'+' '.join(mapping_command))
        raise

# TODO: remove completely (modified awk script sufficient to generate final output of the module) if module calculations in agreement with circRNA_finder output!?
def run_circRNA_finder(star_mapping_prefix, out_prefix, circRNA_finder_path='postProcessStarAlignment.pl'):
    '''
    runs circRNA_finder to detect circular RNAs in the output of STAR
    'star_mapping_prefix': output directory + filename prefix for the mapping of STAR
    'out_prefix': output directory +filename prefix for the detected circRNAs
    'circRNA_finder_path': path to the script postProcessStarAlignment.pl (default current directory)
    '''
    
    # call circRNA_finder script working on STAR output: paths to the directories are required to end with "/"
    circrna_finder_command = ['perl', circRNA_finder_path, star_mapping_prefix, out_prefix]
    
    # execute circRNA_finder as child process
    print('Circular RNA detection with circRNA_finder: '+' '.join(circrna_finder_command))
    try:
        subprocess.check_call(circrna_finder_command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in circRNA_finder call \n'+' '.join(circrna_finder_command))
        raise
    
    # rename output files of circRNA_finder such that prefix of the STAR mapping is not appended to the out_prefix
    p1 = os.path.split(star_mapping_prefix)
    p2 = os.path.split(out_prefix)
    for ending in ['filteredJunctions.bed', 's_filteredJunctions.bed', 's_filteredJunctions_fw.bed']:
        os.rename(os.path.join(p2[0], p2[1]+p1[1]+ending), os.path.join(p2[0], p2[1]+ending))
        

def annotate_and_wirte_output(circrna_finder_prefix, star_prefix, final_out_file, library_type):
    '''
    creates a final output file in the circRNA format used by all watchdog modules
    the method identifies circular reads and annotates the circRNAs of circRNA_finder with the read ids
    furthermore, it adapts the strand of the identified circRNAs based on strandedness of the sequencing library and the strand of the canonical splice site
    input files:
        'circrna_finder_prefix'+s_filteredJunctions.bed: detected circRNAs of circRNA_finder
        'star_prefix'+Chimeric.out.junction: chimeric reads found by STAR
    output file:
        'final_out_file': annotated circRNA finder output with junction coordinates, junction read counts and matching junction read ids
    'library_type': 0=unstranded/unknown, 1=stranded, read1, 2=stranded, read2 (this type is the circRNA_finder default)
    '''
    
    # call awk on chimeric junction file of STAR
    awk_script = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'filterJunctionReads_modified.awk')
    junction_file = star_prefix+'Chimeric.out.junction'
    awk_command=['awk', '-f', awk_script, junction_file]
    awk_out = circrna_finder_prefix+'chimeric_circular.txt'
    print('Searching circular reads with '+' '.join(awk_command))
    with open(awk_out, 'wt') as awk_writer:
        subprocess.check_call(awk_command, stdout=awk_writer)
        
    # read awk result into a dictionary mapping circRNA coordinate -> read ids as list and coordinate -> splice signal
    circ_to_reads={}
    circ_to_signal={}
    with open(awk_out, 'rt') as awk_reader:   
        for line in awk_reader:
            junction_read = line.strip('\n').split(' ')
            # junctions with a splice signal 0=no signal, 1= signal (gtag), 2=reverse complement (ctac)
            splice_signal = int(junction_read[5])
            if(splice_signal>0):
                # add coordinate data
                coordinates = (junction_read[0], junction_read[1], junction_read[2], junction_read[3])
                if(coordinates not in circ_to_reads):
                    circ_to_reads[coordinates]=[junction_read[4]]
                    circ_to_signal[coordinates]=[splice_signal]
                else:
                    circ_to_reads[coordinates].append(junction_read[4])
                    circ_to_signal[coordinates].append(splice_signal)

    # check strandedness of the data
    analyze_library_type(circ_to_signal, library_type)
    
    # decide how to proceed: guess strand from splice site (unknown or unstranded) or infer from library type (stranded)
    final_circ_data={}
    # keep strand for all
    if(library_type==1):
        final_circ_data = circ_to_reads
    # change strand for all
    elif(library_type==2):
        for coordinate,reads in circ_to_reads.items():
            if(coordinate[3]=='+'):
                new_coordinate=(coordinate[0], coordinate[1], coordinate[2], '-')
            else:
                new_coordinate=(coordinate[0], coordinate[1], coordinate[2], '+')
            final_circ_data[new_coordinate]=reads
    # infer strand
    else:
        for coordinate in sorted(circ_to_reads.keys()):
            # keep strand
            if circ_to_signal[coordinate][0]==1:
                new_coordinate=coordinate
            # invert strand
            else:
                if(coordinate[3]=='+'):
                    new_coordinate=(coordinate[0], coordinate[1], coordinate[2], '-')
                else:
                    new_coordinate=(coordinate[0], coordinate[1], coordinate[2], '+')
            # add new coordinate to data -> merge data from splice site with same coordinates but on opposite strands
            if not new_coordinate in final_circ_data:
                final_circ_data[new_coordinate]=[]
            final_circ_data[new_coordinate].extend(circ_to_reads[coordinate])
    
    # sanity check: compare identified circ rna reads with the read counts reported by circrna finder
    error=check_with_cf_result(final_circ_data, circrna_finder_prefix+'s_filteredJunctions.bed', library_type)
    
    # write output with final predictions
    out_columns=['chr', 'start', 'end', 'strand', '#junction_reads', 'junction_reads_ID']
    with open(final_out_file, 'wt') as w:
        w.write('\t'.join(out_columns)+'\n')
        
        for coordinate, read_list in sorted(final_circ_data.items(), key=lambda x: (-len(x[1]), x[0])):
            read_id_str = ','.join(read_list)
            line_out = list(coordinate)+[str(len(read_list)), read_id_str]
            w.write('\t'.join(line_out)+'\n')
        
    return error



def analyze_library_type(splice_signal_data, defined_strandedness):
    '''
    analyzes strands of circular splice sites and strands of splice signals to learn about the strandedness of the library
    'splice_signal_data': mapping circrna coordinates -> list of splice signals (one signal for each read) 1=GTAG (same strand), 2=CTAC (opposite strand)
    'defined_strandedness': strandedness defined by user input
    this function is not necessary for the circrna identification, however it outputs a warning if the strandedness of the data looks different than the specified strandedness
    returns:
        boolean True = input strandedness matches data, False: data suggests different strandedness, None: insufficient circ reads to infer strandedness
        integer giving the inferred strandedness of the data
    '''
     
    # record strand and splice signal to learn about the strandedness of the library
    strand_data={('+',1):0, ('-',1):0, ('+',2):0, ('-',2):0}
    
    # inspect splice signals for all reads of a circular event
    for coord, signal_nr_list in splice_signal_data.items():
        min_sign = min(signal_nr_list)
        max_sign = max(signal_nr_list)
        # all reads at the same circular junction -> should have the same splice signal
        if not min_sign==max_sign:
            sys.stderr.write('ERROR: something is wrong with the splice signal annotation of the reads for '+str(coord)+'\n')
        else:
            strand_data[(coord[3],min_sign)]+=len(signal_nr_list)
        
    # print data about the strandedness of the data
    print('Library type')
    total_reads = sum(strand_data.values())
    print('number circularly spliced reads: '+str(total_reads))
    if total_reads>0:
        read_fractions = {k:v/total_reads for k,v in strand_data.items()}
        for k,v in sorted(read_fractions.items()):
            print(str(k)+': '+str(v))
    guessed_library_type=-1
    
    # do not guess library type if there are insufficient reads
    if(total_reads<40):
        print('Not sufficient circular reads to infer library type')
        return None, defined_strandedness
        
    # infer library type
    else:
        if(read_fractions[('+',1)] <0.05 and read_fractions[('-',1)]  <0.05):
            print('stranded with strand of read2')
            guessed_library_type=2
        elif(read_fractions[('-',2)] <0.05 and read_fractions[('+',2)]  <0.05):
            print('stranded with strand of read1')
            guessed_library_type=1
        else:
            print('unstranded library')
            guessed_library_type=0
            
        # compare inferred strandedness to strandedness passed as input parameter
        if defined_strandedness==guessed_library_type:
            print('in agreement with specified strandedness '+str(defined_strandedness))
            return True, guessed_library_type
        else:
            sys.stderr.write('WARNING: disagreement between inferred strandedness '+str(guessed_library_type)+' and specified strandedness '+str(defined_strandedness)+'\n'
                             'consider running with option --strandedLibrary '+str(guessed_library_type)+'\n')
            return False, guessed_library_type


# TODO: consider removing the call of circrna finder and use only the modified awk script
def check_with_cf_result(circ_data, cf_file, library_type):
    '''
    checks the computed circ_rnas against the raw output of circRNA_finder
    'circ_data': mapping of circrna coordinates to read ids derived from modified awk script
    'cf_file': file with the original circRNA_finder output
    'library_type': 0 (unstranded), 1 (stranded, read1), 2 (stranded, read2)
    '''
    
    mismatch=False
    
    # read data from circrna finder
    cf_coord_to_count=collections.defaultdict(int)
    with open(cf_file, 'rt') as cfr:
        for line in cfr:
            bed_cells = line.strip('\n').split('\t')
            # stranded based on read2: compare data is it is
            if library_type==2:
                cf_coord = (bed_cells[0], bed_cells[1], bed_cells[2], bed_cells[5])
            # stranded based on read1: compare circ rnas with inverted strand
            elif library_type==1:
                cf_coord = (bed_cells[0], bed_cells[1], bed_cells[2], '-') if bed_cells[5]=='+' else (bed_cells[0], bed_cells[1], bed_cells[2], '+')
            # unstranded: compare chromosome, start and end but not the strand
            else:
                cf_coord = (bed_cells[0], bed_cells[1], bed_cells[2])
            cf_abundance = int(bed_cells[4])
            cf_coord_to_count[cf_coord]+=cf_abundance
                
    # read processed data
    if library_type>0:
        calc_coord_to_count={k:len(v) for k,v in circ_data.items()}
    else:
        calc_coord_to_count={(k[0], k[1], k[2]):len(v) for k,v in circ_data.items()}
        
    # compare the data
    fp = set(calc_coord_to_count.keys())-set(cf_coord_to_count.keys())
    if(len(fp)>0):
        sys.stderr.write('ERROR: could not find some circ rnas in the circrna finder output\n'+'\n'.join(str(x) for x in fp)+'\n')
        mismatch=True
    fn = set(cf_coord_to_count.keys())-set(calc_coord_to_count.keys())
    if(len(fn)>0):
        sys.stderr.write('ERROR: missed some circ rnas from the circrna finder output\n'+'\n'.join([str(x) for x in fn])+'\n')
        mismatch=True
    for circ_rna in set(calc_coord_to_count.keys()) & set(cf_coord_to_count.keys()):
        count_exp = cf_coord_to_count[circ_rna]
        count_act = calc_coord_to_count[circ_rna]
        if not count_exp == count_act:
            sys.stderr.write('ERROR: read count mismatch for circ rna '+str(circ_rna)+' '+str(count_exp)+' vs. '+str(count_act)+'\n')
            mismatch=True
    
    return mismatch


if __name__ == '__main__':
    pass