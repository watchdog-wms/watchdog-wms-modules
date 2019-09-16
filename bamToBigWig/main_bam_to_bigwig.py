'''
Created on Oct 8, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse

from call_bam_to_bigwig import run_bam_coverage
import watchdog_utils as wutils

## TODO: option for strand-specific data

def get_command_line_options():
    '''command line parser for the watchdog module that converting bam to bigwig format using deepTools bamCoverage'''
    
    parser = argparse.ArgumentParser(description='Wrapper for the deepTools bamCoverage')
    
    # input & output files
    parser.add_argument('--inBam', required=True, default=None, type=wutils.valid_indexed_bam,
        metavar='reads.bam', help='file with mapped reads in bam format with an index')
    parser.add_argument('--outBw', required=True, default=None,
        metavar='reads.bw', help='file for writing read data in bigWig format')
    parser.add_argument('--returnFilePath', required=False, default=None,
        metavar='watchdogdir/tabfile', help='flag that is used by watchdog only for handling the task output')
        
    #options for bamCoverage
    parser.add_argument('--bamCoveragePath', required=False, default='bamCoverage', type=wutils.valid_exec,
        metavar='bamCoverage', help='path to bamCoverage executable (default: use executable from PATH)')
    parser.add_argument('--binSize', required=False, default=1, type=wutils.positive_integer,
        metavar='binsize', help='resolution of bigwig file, setting binsize to a value >1 decreases the size of the bigwig file, but decreases the resolution of the genomic coverage (default:1)')
    parser.add_argument('--numberOfProcessors', required=False, default=1, type=wutils.positive_integer,
        metavar='prNr', help='number of parallel processes to use')
        
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def create_outfiles(options):
    ''' create all parent directories of the resulting bigwig file '''
     
    wutils.make_parent_dirs(options.outBw)


def main():
    ''' main method for converting bam -> bigwig using the watchdog module '''
    
    # check options and prepare output location
    _,o = get_command_line_options()
    create_outfiles(o)
      
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the module for converting bam to bigwig format !\n')
      
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # run the main program
    run_bam_coverage(o.inBam, o.outBw, o.binSize, o.numberOfProcessors, o.bamCoveragePath) 

    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n-> the bigwig file is located at \''+o.outBw+'\'')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)
  
    # write return file for watchdog
    if(o.returnFilePath is not None):
        returnVars = [('bigWigFile', o.outBw)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)
    
    
if __name__ == '__main__':
    main()

