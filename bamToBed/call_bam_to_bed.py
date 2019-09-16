'''
Created on Oct 9, 2018

@author: friedl
'''

import sys
import subprocess

def run_bedtools_bamtobed(bam_file, bed_file, bedtools_path='bedtools', split=True):
    ''' calls external program bedtools bamtobed to convert a bam file into bed format '''
    
    command = [bedtools_path, 'bamtobed']
    if split:
        command.append('-split')
    command+=['-i', bam_file]
    
    with open(bed_file, 'wt') as bedw:
        try:
            subprocess.check_call(command, stdout=bedw)
        except subprocess.CalledProcessError:
            sys.stderr.write('Error in bedtools call \n'+' '.join(command))
            raise
