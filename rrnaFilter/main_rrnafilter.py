'''
Created on Jun 29, 2018

@author: friedl
'''

import argparse
import os
import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')

import watchdog_utils as wutils
import rrna_filter_algo.runFilter

# TODO: option to throw exception if few non-rrna reads

def get_command_line_options():
    '''command line parser for the rrna filter running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='removes rrna reads from sequencing data')
        
    # input files
    groupInput = parser.add_argument_group('Input Options')
    groupInput.add_argument('--in1', required=True, type=wutils.valid_fastq_path,
        metavar='dir/reads1.fq.gz', help='first fastQ (.fq, .fastq) or gzipped fastQ (.fq.gz) file with the sequenced reads')
    groupInput.add_argument('--in2', required=False, type=wutils.valid_fastq_path, default = None, 
        metavar='dir/reads2.fq.gz', help='second fastQ (.fq, .fastq) or gzipped fastQ (.fq.gz) file with the sequenced reads (for paired-end data only)')
    groupInput.add_argument('--rrnaIndex', required=True, type=wutils.valid_bwa_index,
        metavar='dir/PrefixIndexFiles', help='Common prefix of bwa index files for rrna sequence')

    # output files
    groupOutput = parser.add_argument_group('Output Options')
    groupOutput.add_argument('--out1', required=True, type=wutils.valid_outfile_ending_fastq_or_fasta,
        metavar='dir/filtered_reads1.fq.gz', help='file (ending .fq.gz, .fq, .fastq or .fa) for writing non-rrna reads from in1')
    groupOutput.add_argument('--out2', required=False, type=wutils.valid_outfile_ending_fastq_or_fasta, default=None,
        metavar='dir/filtered_reads2.fq.gz', help='file (ending .fq.gz, .fq, .fastq or. fa) for writing non-rrna reads from in2 (for paired-end data)')
    groupOutput.add_argument('--sam', required=True, type=wutils.valid_outfile_ending_sam,
        metavar='dir/mapped_rrna.sam', help='sam file (ending .sam) for writing rrna reads from in1 and in2')
    groupOutput.add_argument('--workdir', required=False, default=os.getcwd(), type=wutils.valid_folder_path,
        metavar='dir', help='path to directory for writing large temporary files (content is deleted at the end of execution), default: current directory')
    tmp = groupOutput.add_mutually_exclusive_group()
    tmp.add_argument('--keepTmp', required=False, default=False, action='store_true',
        help='option to keep temporary files (default: false)')
    tmp.add_argument('--nokeepTmp', required=False, default=True, action='store_false',
        help='option to remove temporary files, used by watchdog boolean parameter (default:true)')
    groupOutput.add_argument('--returnFilePath', required=False, default=None,
        metavar='watchdogdir/tabfile', help='flag that is used by watchdog only for handling the task output')
        
    # filter options
    groupFilter = parser.add_argument_group('Filter Options', 'If an option (except for the pairFiltering) is unset, no filtering is performed')
    groupFilter.add_argument('--maxEditDistance', required=False, type=wutils.positive_integer_or_zero, default=None,
        metavar='max_ed', help='maximum allowed edit distance for a read alignment against rrna, default: not limited')
    groupFilter.add_argument('--maxMismatches', required=False, type=wutils.positive_integer_or_zero, default=None,
        metavar='max_mm', help='maximum allowed number of mismatches for a read alignment against rrna, default: not limited')
    groupFilter.add_argument('--maxIndels', required=False, type=wutils.positive_integer_or_zero, default=None,
        metavar='max_indel', help='maximum allowed number of indels for a read alignment against rrna, default: not limited')
    groupFilter.add_argument('--pairFiltering', required=False, type=int, choices=[1,2], default=2,
        help='Number of reads of a pair required to fulfil the options above, default: 2')
        
    #bwa options
    groupBWA = parser.add_argument_group('BWA Options')
    groupBWA.add_argument('--bwaPath', required=False, type=wutils.valid_exec, default='bwa',
        metavar='bwaPath', help='path to bwa executable')
    groupBWA.add_argument('--seedSize', required=False, type=wutils.positive_integer, default=25,
        metavar='seedsize', help='size of initial seed for bwa (-k option of bwa), default: 25')
    groupBWA.add_argument('--threads', required=False, type=wutils.positive_integer, default=1,
        metavar='NumOfThreads', help='number of threads to use for bwa (-t option of bwa), default: 1')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def check_input_files(parser, options):
    ''' check manually combination of input files -> for each input file should be an outputfile '''
    
    if(options.in2 is not None and options.out2 is None):
        parser.error('Missing option --out2 for paired-end data (flag in2 is used)!')
    if(options.out2 is not None and options.in2 is None):
        parser.error('Missing option --in2 for paired-end data (flag out2 is used)!')
    

def create_outfiles(options):
    ''' creates all parent folders for the output files if they do not exist'''
    
    wutils.make_parent_dirs(options.out1)
    wutils.make_parent_dirs(options.sam)
    if(options.out2 is not None):
        wutils.make_parent_dirs(options.out2)

        
def main():
    ''' main function calling the rrna filter '''
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # check options and prepare output locations
    p,o = get_command_line_options()
    check_input_files(p,o)
    create_outfiles(o)
    
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the rRNA filter module!\n')
    
    # run the main program
    rrna_filter_algo.runFilter.run(in1=o.in1, in2=o.in2, rrnaIndex=o.rrnaIndex, workdir=o.workdir, out1=o.out1, out2=o.out2, sam=o.sam,
                       seedSize=o.seedSize, threads = o.threads, maxEditDistance=o.maxEditDistance, maxMismatches=o.maxMismatches, maxIndels=o.maxIndels,
                       pairFiltering=o.pairFiltering, bwapath=o.bwaPath, keeptmp=o.keepTmp)
    
    # log
    end_timepoint=wutils.get_current_time()
    res_files = o.out1
    if(o.out2 is not None):
        res_files+= ', '+o.out2
    print(end_timepoint[1]+'Module finished successfully!\n-> check out \''+res_files+'\' for the non-rRNA reads')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)
    
    # write return file for watchdog
    if(o.returnFilePath is not None):
        returnVars = [('rrnaSAMFile', o.sam), ('filteredFQ1', o.out1)]
        if(o.out2 is None):
            returnVars.append(('filteredFQ2', 'not_defined_for_single_end'))
        else:
            returnVars.append(('filteredFQ2', o.out2))
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)
        

if __name__ == '__main__':
    main()


    
