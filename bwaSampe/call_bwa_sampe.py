'''
Created on Jul 13, 2018

@author: friedl
'''

import subprocess
import sys

def run_bwa_sampe(genome_index, fastq_1, fastq_2, sai_1, sai_2, sam_out, bwa_path='bwa', index_to_ram=False):
    ''' calls bwa sampe to generate a sam file from 2 sai files for paired-end sequencing data'''
    
    # build bwa command
    command=[bwa_path, 'sampe', '-f', sam_out]
    if index_to_ram:
        command.append('-P')
    command+=[genome_index, sai_1, sai_2, fastq_1, fastq_2]
    
    # execute bwa command
    print('Running command:\n'+' '.join(command))
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in BWA call \n'+' '.join(command))
        raise