'''
Created on Jul 13, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse

from call_bwa_aln import run_bwa_aln
import watchdog_utils as wutils



def get_command_line_options():
    '''command line parser for running bwa aln as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Wrapper for the bwa aln tool')
    
    # input files
    parser.add_argument('--inReads', required=True, type=wutils.valid_uncompressed_fastq_path,
        metavar='dir/reads.fastq', help='fastq file with the sequenced reads')
    parser.add_argument('--bwaIndex', required=True, type=wutils.valid_bwa_index,
        metavar='dir/PrefixIndexFiles', help='Common prefix of bwa index files for the reference genome')

    # output files
    parser.add_argument('--outSai', required=True,
        metavar='dir/mapped_reads.sai', help='file for writing mapped reads in bwa format')
    parser.add_argument('--returnFilePath', required=False, default=None,
        metavar='watchdogdir/tabfile', help='flag that is used by watchdog only for handling the task output')
        
    #bwa options
    parser.add_argument('--bwaPath', required=False, default='bwa', type=wutils.valid_exec,
        metavar='bwa', help='path to BWA executable (default: use executable from PATH)')
    parser.add_argument('--threads', required=False, default=1, type=wutils.positive_integer, 
        metavar='NumOfThreads', help='number of threads to use for bwa aln (-t option of bwa), default: 1')
    parser.add_argument('--stopIfMoreThanBestHits', required=False, default=None, type=wutils.positive_integer,
        metavar='NumOfHits', help='stop searching when there are more than that many best hits (default: bwa default)')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def create_outfiles(options):
    ''' create all parent directories of the sai file with the output of bwa aln '''
    
    wutils.make_parent_dirs(options.outSai)
    
    
def main():
    ''' main method for running the bwa aln wrapper as watchdog module '''
    
    # check options and prepare output locations
    _,o = get_command_line_options()
    create_outfiles(o)
    
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the module for running bwa aln !\n')
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # run the main program
    run_bwa_aln(o.bwaIndex, o.inReads, o.outSai, o.bwaPath, o.threads, o.stopIfMoreThanBestHits)
    
    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n-> check out \''+o.outSai+'\' for the output of bwa aln')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)

    # write return file for watchdog
    if(o.returnFilePath is not None):
        returnVars = [('bwaSaiFile', o.outSai)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)



if __name__ == '__main__':
    main()
    
    
    
    
