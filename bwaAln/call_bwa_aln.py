'''
Created on Jul 13, 2018

@author: friedl
'''

import subprocess
import sys

def run_bwa_aln(genome_index, fastq_in_file, sai_out_file, bwa_path='bwa', threads=1, R_param=None):
    ''' calls bwa aln '''
    
    # build bwa command
    command = [bwa_path, 'aln', '-t', str(threads), '-f', sai_out_file]
    if(R_param is not None):
        command+=['-R', str(R_param)]
    command+=[genome_index, fastq_in_file]
    
    # execute bwa command
    print('Running command:\n'+' '.join(command))
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in BWA call \n'+' '.join(command))
        raise