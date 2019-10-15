'''
Created on Oct 9, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse
import watchdog_utils as wutils
from call_bam_to_bed import run_bedtools_bamtobed

def get_command_line_options():
    ''' commandline parser for the bamToBed watchdog module that for converting bam to bed format'''
    
    parser = argparse.ArgumentParser(description='Wrapper for the bedtools bamtobed')
    
    # input & output files
    parser.add_argument('--inBam', required=True, default=None, type=wutils.valid_file_path,
        metavar='reads.bam', help='file with mapped reads in bam format (no index required)')
    parser.add_argument('--outBed', required=True, default=None,
        metavar='reads.bed', help='file for writing read data in bed format')
    parser.add_argument('--returnFilePath', required=False, default=None,
        metavar='watchdogdir/tabfile', help='flag that is used by watchdog only for handling the task output')
        
    #options for bamtobed
    parser.add_argument('--bedtoolsPath', required=False, default='bedtools', type=wutils.valid_exec,
        metavar='bedtools', help='path to bedtools executable (default: use executable from PATH)')
    split_mux = parser.add_mutually_exclusive_group()
    split_mux.add_argument('--split', required=False, action='store_true', default=True,
                           help='if this flag is set, the aligned split reads (cigar contains N) will be written as separate intervals in the bed file, set by default')
    split_mux.add_argument('--nosplit', required=False, action='store_true', default=False,
                           help='if this flag is set, the aligned split reads (cigar contains N) will be written as a single interval in the bed file')
        
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def create_outfiles(options):
    ''' creates all parent directories of the resulting bed file that do not exist '''
    
    wutils.make_parent_dirs(options.outBed)


def main():
    ''' main method for converting bam -> bed using the watchdog module '''
    
    # check options and prepare output location
    _,o = get_command_line_options()
    create_outfiles(o)
      
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the module for converting bam to bed format !\n')
      
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # run the main program
    split_spliced = True
    if o.nosplit:
        split_spliced = False
    run_bedtools_bamtobed(bam_file=o.inBam, bed_file=o.outBed, bedtools_path=o.bedtoolsPath, split=split_spliced)

    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n-> the bed file is located at \''+o.outBed+'\'')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)
  
    # write return file for watchdog
    if(o.returnFilePath is not None):
        returnVars = [('bedFile', o.outBed)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)
    
    
if __name__ == '__main__':
    main()
