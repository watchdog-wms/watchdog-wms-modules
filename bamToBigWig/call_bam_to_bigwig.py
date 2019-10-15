'''
Created on Oct 8, 2018

@author: friedl
'''

import subprocess
import sys

def run_bam_coverage(bam_file, bigwig_file, binsize=1, process_number=1, bamcoverage_path='bamCoverage'):
    ''' builds and executes bamCoverage commmand to convert a bam file into a bigwig file'''
    
    command=[bamcoverage_path, '--bam', bam_file, '--outFileName', bigwig_file, '--outFileFormat', 'bigwig', '--binSize', str(binsize), '--numberOfProcessors', str(process_number)]
    
    # execute bamCoverage command
    print('Running command:\n'+' '.join(command))
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in bamCoverage call \n'+' '.join(command))
        raise