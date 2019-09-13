'''
Created on 15 Nov 2017

@author: friedl
'''

import subprocess
import os.path
import sys


def buildIndex(genome_file, index_prefix, bwa_path='bwa'):
    ''' 
    builds a bwa index for a fasta file 
    'genome_file': genome sequence in fasta format for which the index is generated
    'index_prefix': output directory with file name prefix for the new bwa index
    'bwa_path': path to the BWA executable, default: bwa in $PATH
    '''
    
    bwa_index_command = [bwa_path, 'index', '-p', index_prefix, '-a', 'bwtsw', genome_file]
    print('Indexing reference\n'+' '.join(bwa_index_command))
    
    try:
        subprocess.check_call(bwa_index_command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in BWA call \n'+' '.join(bwa_index_command))
        raise


def runBWA(fastq_reads1, fastq_reads2, index_prefix, sam_out, bwa_threads=1, bwa_path='bwa', bwa_seed_size=19, bwa_score_threshold=30):
    '''
    runs bwa mem to map sequenced reads (paired-end data)
    'fastq_reads1' & 'fastq_reads_2': sequences of paired reads in fastq format
        for single-end data: fastq_reads2 is set to None
    'index_prefix': output directory with file name prefix for the bwa index (representing the reference genome)
    'sam_out': path for output file of BWA in SAM format
    'bwa_threads': number of threads to use for BWA
    'bwa_path': path to the BWA executable, default: bwa in $PATH
    'bwa_seed_size': minimum seed length (bwa parameter -k)
    'bwa_score_threshold': does not output alignments with alignment score below the threshodl (bwa parameter -T)
    '''
    
    # call for single-end data
    bwa_mem_command = [bwa_path, 'mem', '-T', str(bwa_score_threshold), '-k', str(bwa_seed_size), '-t', str(bwa_threads), index_prefix, fastq_reads1]
    # extend command for paired-end data
    if fastq_reads2 is not None:
        bwa_mem_command.append(fastq_reads2)
    print('Mapping with BWA\n'+' '.join(bwa_mem_command))
    
    try:
        subprocess.check_call(bwa_mem_command, stdout=open(sam_out, 'wt'))
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in BWA call \n'+' '.join(bwa_mem_command))
        raise
    


def runCIRI2(bwa_sam_file, ciri_out_file, genome_file, annotation_file=None, ciri_threads=1, ciri_script='CIRI2.pl', ciri_stringency='high', keep_tmp=False):
    '''
    runs CIRI2 to detect circRNAs in reads mapped by BWA
    'bwa_sam_file': output of BWA mem in SAM format
    'ciri_out_file': file for writing the output of CIRI2 -> log file written to ciri_out_file.log
    'genome_file': reference genome that was used for the mapping in 'bwa_sam_file'
    'annotation_file': gtf file with gene annotations for 'genome_file' (optional)
    'ciri_threads': number of threads to use for CIRI -> ATTENTION each thread requires >10GB RAM
    'ciri_script': path to the CIRI2 perl script'
    'ciri_stringency': influences filtering of detected circular RNAs according to circ reads, number of cigars and number of false positive reads,
        three possible values 'high', 'medium' or 'low'
    'keep_temp': if True, the temporary files of ciri are not deleted at the end of the run
    '''
    
    ciri_command = ['perl', ciri_script, '-T', str(ciri_threads), '-I', bwa_sam_file, '-O', ciri_out_file, '-F', genome_file]
    if annotation_file is not None:
        ciri_command+=['-A', annotation_file] 
    if keep_tmp is True:
        ciri_command+=['-D']
    if(ciri_stringency=='high'):
        ciri_command+=['-high']
    elif(ciri_stringency=='medium'):
        ciri_command+=['-low']
    elif(ciri_stringency=='low'):
        ciri_command+=['-0']
        
    print('Running CIRI2\n'+' '.join(ciri_command))
    
    # working directory for CIRIerror.log
    wkdir = os.path.dirname(ciri_out_file)
    
    # call CIRI and report if an error occurs
    try:
        subprocess.check_call(ciri_command, cwd=wkdir)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in CIRI call \n'+' '.join(ciri_command))
        raise


def writeOutput(raw_ciri_out_file, transformed_ciri_file):
    '''
    transforms CIRI output into a file format shared with circRNA finder module for easy annotation and merging of predictions
    'raw_ciri_out_file': output of CIRI2 script
    'transformed_ciri_file': file with 0-based coordinates and junction reads of all circRNAs reported in 'raw_ciri_out_file'
    '''
    
    #header of output file format
    out_columns=['chr', 'start', 'end', 'strand', '#junction_reads', 'junction_reads_ID']
    
    with open(raw_ciri_out_file, 'rt') as reader, open(transformed_ciri_file, 'wt') as writer:
        
        # read and write header line giving the table format
        first_line = reader.readline()
        if(first_line.startswith('circRNA_ID')):
            writer.write('\t'.join(out_columns)+'\n')
        else:
            raise ValueError('Unknown file format: Missing CIRI file header in '+str(raw_ciri_out_file)+'\n'+first_line)
        
        # each line corresponds to a detected circ rna
        for line in reader:
            # transform required columns
            old_row = line.strip('\n').split('\t')
            new_row = []
            if(len(old_row)<12):
                raise ValueError('Unknown file format: Expecting 12 tab-separated columns in '+str(raw_ciri_out_file)+'\n'+line)
            # keep chromosome name
            new_row.append(old_row[1])
            # transform start to 0-based coordinate
            new_row.append(str(int(old_row[2])-1))
            # 1-based inclusive end position = 0-based exclusive end position
            new_row.append(old_row[3])
            # keep strand
            new_row.append(old_row[10])
            # get junction reads, remove last ',' symbol from list of reads
            new_row.append(old_row[4])
            read_list=old_row[11].rstrip(',')
            new_row.append(read_list)
            writer.write('\t'.join(new_row)+'\n')
                
                
if __name__ == '__main__':
    pass