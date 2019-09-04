'''
Created on 13 Nov 2017

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse
import os.path

import watchdog_utils as wutils
import algo_trimmedFastqPairFilter

def get_command_line_options():
    '''command line parser for running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Extracts paired reads from 2 fastq files.')
    
    # input 1 & 2 or input prefix
    parser.add_argument('--inPrefix', required=False, default=None, metavar='fastq_prefix',
                        help='reads in 2 fastq files: prefix1.[fastq|fq], prefix2.[fastq|fq]')
    parser.add_argument('--inReads1', required=False, default=None, type=wutils.valid_file_path, metavar='reads_1.fastq',
                        help='path to first fastq file with reads')
    parser.add_argument('--inReads2', required=False, default=None, type=wutils.valid_file_path, metavar='reads_2.fastq',
                        help='path to second fastq file with reads')
    
    # output 1,2,3 or output prefix
    parser.add_argument('--outPrefix', required=False, default=None, metavar='out_prefix',
                        help='writes output to three files: prefix1.fastq, prefix2.fastq, prefixsingleton.fastq')
    parser.add_argument('--outReads1', required=False, default=None, metavar='paired_reads_1.fastq',
                        help='output file for first reads of paired data')
    parser.add_argument('--outReads2', required=False, default=None, metavar='paired_reads_2.fastq',
                        help='output file for second reads of paired data')
    parser.add_argument('--outSingletons', required=False, default=None, metavar='singleton_reads.fastq',
                        help='output file for singleton reads without a mate')
    
    # watchdog return file
    parser.add_argument('--returnFilePath', required=False, default=None, metavar='watchdog_variables',
                        help='internal watchdog command line parameter')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options
    
    
def check_and_get_inputfiles(parser, cmdl_options):
    ''' check and extract paths of input files '''
    
    infiles=['', '']
    
    # 2 input file paths
    if(cmdl_options.inPrefix is None):
        for pos,(nr,option) in enumerate([('1', cmdl_options.inReads1), ('2', cmdl_options.inReads2)]):
            # check if option is used -> file existence handled by argument parser
            if(option is None):
                parser.error('Missing argument: inPrefix or inReads'+nr+' option is required!')
            #assign file name
            infiles[pos]=option
                
    # 1 prefix for both input files
    else:
        prefix = cmdl_options.inPrefix
        for pos,(nr,option) in enumerate([('1', cmdl_options.inReads1), ('2', cmdl_options.inReads2)]):
            if(option is not None):
                parser.error('Argument clash: use only inPrefix or inReads'+nr+'!')
            if(os.path.exists(prefix+nr+'.fastq')):
                if(os.path.exists(prefix+nr+'.fq')):
                    parser.error('Found 2 files matching the prefix: '+prefix+nr+'.fastq and '+prefix+nr+'.fq\n Use option inReads'+nr+' to specify your file!')
                else:
                    infiles[pos]=prefix+nr+'.fastq'
            elif(os.path.exists(prefix+nr+'.fq')):
                infiles[pos]=prefix+nr+'.fq'
            else:
                parser.error('Cannot find file matching the prefix: '+prefix+nr+'.[fastq|fq]')
    
    return infiles
                
    
def get_and_create_outfiles(parser, cmdl_options):   
    ''' check output options and create output files '''
    
    outfiles=['', '', '']
    
    # 3 separate paths for output file
    if(cmdl_options.outPrefix is None):
        for pos,(oname,option) in enumerate([('outReads1', cmdl_options.outReads1), ('outReads2', cmdl_options.outReads2), ('outSingletons', cmdl_options.outSingletons)]):
            # check if option is used
            if(option is None):
                parser.error('Missing argument: outPrefix or '+oname+' option is required!')
            outfiles[pos]=option
    
    # 1 shared prefix for the output 
    else:
        for pos,(oname,option) in enumerate([('outReads1', cmdl_options.outReads1), ('outReads2', cmdl_options.outReads2), ('outSingletons', cmdl_options.outSingletons)]):
            if(option is not None):
                parser.error('Argument clash: use only outPrefix or '+oname+'!')
        outfiles[0]=cmdl_options.outPrefix+'1.fastq'
        outfiles[1]=cmdl_options.outPrefix+'2.fastq'
        outfiles[2]=cmdl_options.outPrefix+'singleton.fastq'
    
    # make parent folders for every outfile
    for ofile in outfiles:
        wutils.make_parent_dirs(ofile)
    
    return outfiles


if __name__ == '__main__':
    
    print('Program call:')
    print(' '.join(sys.argv))
    
    p,o = get_command_line_options()
    inpaths = check_and_get_inputfiles(p, o)
    outpaths = get_and_create_outfiles(p, o)
    
    # start of execution
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Starting trimmedFastqPairFilter module...')
    
    print('Input:')
    print(inpaths)
    print('Output:')
    print(outpaths)
    
    algo_trimmedFastqPairFilter.filterFastqs(inpaths[0], inpaths[1], outpaths[0], outpaths[1], outpaths[2])
    
    # end of execution
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'trimmedFastqPairFilter module finished, check out \''+','.join(outpaths)+'\' for the results')
    
    # output resources consumed by the script
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=False)
    
    # write return variables
    if(o.returnFilePath is not None):
        returnVars = [('pairedReads1', outpaths[0]), ('pairedReads2', outpaths[1]), ('singletonReads', outpaths[2])]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)
    
    