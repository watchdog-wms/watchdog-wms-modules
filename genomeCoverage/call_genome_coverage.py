'''
Created on Feb 5, 2019

@author: friedl
'''

import subprocess
import sys
import os

def bam_to_bedgraph(bam_path, bedgraph_path, bedtools_exec='bedtools'):
    ''' converts a bam file to bedgraph format using bedtools genome coverage '''
    
    command = [bedtools_exec, 'genomecov', '-bg', '-split', '-ibam', bam_path]
    print('Running command\n'+' '.join(command))
    
    with open(bedgraph_path, 'wt') as bgw:
        try:
            subprocess.check_call(command, stdout=bgw)
        except subprocess.CalledProcessError:
            sys.stderr.write('Error in bedtools call \n'+' '.join(command))
            raise


def bedgraph_to_tdf(bedgraph_path, genome_file, tdf_path, igvtools_exec='igvtools'):
    ''' converts a bedgraph file to tdf format using igvtools '''
    
    # set working directory for igvtools to write log file in the directory where the output file is located
    wdir = os.path.dirname(tdf_path)
    command = [igvtools_exec, 'toTDF', bedgraph_path, tdf_path, genome_file]
    print('Running command\n'+' '.join(command))
    
    try:
        subprocess.check_call(command, cwd=wdir)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in igvtools call \n'+' '.join(command))
        raise
    
    # remove log file of igv -> empty file
    os.remove(os.path.join(wdir, 'igv.log'))

