'''
Created on Oct 17, 2018

@author: friedl
'''

import subprocess
import sys

def call_phantom_peak_script(bamfile, summaryfile, plotfile, rdatafile, spp_script, rscript='Rscript', tmpdir=None, threads=1):
    
    # build command with input paths & output paths & threads
    command=[rscript, spp_script, '-c='+bamfile, '-p='+str(threads), '-savd='+rdatafile, '-savp='+plotfile, '-out='+summaryfile, '-rf']
    if tmpdir is not None:
        command.append('-tmpdir='+tmpdir)
        
    # run command
    try:
        subprocess.check_call(command)
    except subprocess.CalledProcessError:
        sys.stderr.write('Error in phantompeakqualtools call \n'+' '.join(command))
        raise
    