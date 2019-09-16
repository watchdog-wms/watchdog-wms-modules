'''
Created on Oct 17, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse
import watchdog_utils as wutils
from call_phantompeakqualtools import call_phantom_peak_script

def get_command_line_options():
    ''' commandline parser for the wrapper of phantompeakqualtools in watchdog'''
    
    parser = argparse.ArgumentParser(description='Wrapper for the phantompeakqualtools')
    
    # input & output files
    parser.add_argument('--inBam', required=True, default=None, type=wutils.valid_bam,
        metavar='reads.bam', help='file with mapped reads in bam format (ending *.bam required but no index required)')
    parser.add_argument('--outPrefix', required=True, default=None,
        metavar='opref', help='the module creates 3 files: opref.txt (summary), opref.pdf (phantom peak plot) & opref.rData (data used for plotting)')
    parser.add_argument('--tmpdir', required=False, default=None, type=wutils.valid_folder_path,
        metavar='tmpdir', help='folder for temporary files, the bam file will be copied there and renamed with a random suffix, default: tempdir() of R')
    
    #technical stuff: path to rscript, path to spp -> default in PATH
    parser.add_argument('--rscriptPath', required=False, default='Rscript', type=wutils.valid_exec,
        metavar='Rscript', help='path to Rscript executable, default: use Rscript from $PATH')
    # for upgrade to SLE 15: install spp from https://github.com/hms-dbmi/spp and 
    # manually add "library(caTools)" to the script run_spp.R from https://github.com/kundajelab/phantompeakqualtools
    parser.add_argument('--sppPath', required=True, default=None, type=wutils.valid_file_path,
        metavar='run_spp.R', help='path to run_spp.R script')
    
    # addtional settings
    parser.add_argument('--threads', required=False, default=1, type=wutils.positive_integer,
        metavar='threadNr', help='Number of threads, default: 1')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def create_outfiles(options):
    ''' creates all parent directories given in the output prefix and returns the resulting file paths for the summary, plot and rdata file '''
    
    wutils.make_parent_dirs(options.outPrefix)
    # summary file, plot file, data file
    return (options.outPrefix+'.txt', options.outPrefix+'.pdf', options.outPrefix+'.Rdata')


def main():
    ''' main method for running phantompeakqualtools as watchdog module '''
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')

    # check options and prepare output location
    _,o = get_command_line_options()
    paths = create_outfiles(o)
      
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the module for running phantompeakqualtools!\n')  
    
    # run the external program
    call_phantom_peak_script(o.inBam, paths[0], paths[1], paths[2], o.sppPath, o.rscriptPath, o.tmpdir, o.threads)

    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n-> check out the results:\n'+'\n'.join(paths))
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)
  
    
if __name__ == '__main__':
    main()