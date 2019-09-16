'''
Created on Jul 13, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse

from call_bwa_sampe import run_bwa_sampe
import watchdog_utils as wutils


def get_command_line_options():
    '''command line parser for running bwa sampe as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Wrapper for the bwa sampe tool')
    
    # input files
    parser.add_argument('--inReads1', required=True, type=wutils.valid_uncompressed_fastq_path,
        metavar='dir/reads1.fastq', help='uncompressed fastq (.fq, .fastq) file with the sequenced reads')
    parser.add_argument('--inReads2', required=True, type=wutils.valid_uncompressed_fastq_path,
        metavar='dir/reads2.fastq', help='uncompressed fastq (.fq, .fastq) file with the sequenced reads (mates)')
    parser.add_argument('--inSai1', required=True, type=wutils.valid_file_path,
        metavar='dir/reads1.sai', help='output of bwa aln for the file given by inReads1')
    parser.add_argument('--inSai2', required=True, type=wutils.valid_file_path,
        metavar='dir/reads2.sai', help='output of bwa aln for the file given by inReads2')
    parser.add_argument('--bwaIndex', required=True, type=wutils.valid_bwa_index,
        metavar='dir/PrefixIndexFiles', help='Common prefix of bwa index files for the reference genome')

    # output files
    parser.add_argument('--outSam', required=True, default=None,
        metavar='dir/mapped_reads.sam', help='file for writing mapped reads in sam format')
    parser.add_argument('--returnFilePath', required=False, default=None,
        metavar='watchdogdir/tabfile', help='flag that is used by watchdog only for handling the task output')
        
    #bwa options
    parser.add_argument('--bwaPath', required=False, default='bwa', type=wutils.valid_exec,
        metavar='bwa', help='path to BWA executable (default: use executable from PATH)')
    pflag = parser.add_mutually_exclusive_group()
    pflag.add_argument('--indexInRam', required=False, action='store_true', default=False,
        help='option to load complete index into main memory (default: false)')
    pflag.add_argument('--noindexInRam', required=False, action='store_false', default=True,
        help='option to keep index partly on disk (default: true)')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def create_outfiles(options):
    ''' create all parent directories of the sam file with the aligned reads '''
     
    wutils.make_parent_dirs(options.outSam)
     
     
def main():
    ''' main method for running the bwa aln wrapper as watchdog module '''
     
    # check options and prepare output locations
    _,o = get_command_line_options()
    create_outfiles(o)
     
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the module for running bwa sampe !\n')
     
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
     
    # run the main program
    run_bwa_sampe(o.bwaIndex, o.inReads1, o.inReads2, o.inSai1, o.inSai2, o.outSam, o.bwaPath, o.indexInRam)
     
    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n-> check out \''+o.outSam+'\' for the output of bwa aln')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)
 
    # write return file for watchdog
    if(o.returnFilePath is not None):
        returnVars = [('bwaPairedSamFile', o.outSam)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)



if __name__ == '__main__':
    main()