'''
Created on Mar 27, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse

import watchdog_utils as wutils
import algo_mappingsum
import plot_mappingsum
import txttables.tablefunc as tf


def get_command_line_options():
    '''command line parser for the mappingsum running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Sophies first Watchdog Module "mappingSummary"'
                                    'to summarize read counts after different steps in a RNA-seq workflow')
    
    # options to summarize reads in FASTQ files (unmapped)
    qc = parser.add_argument_group('Summarize read counts before and after quality control with FASTQC')
    qc.add_argument('--basicStatsSummary', required=False, default=None, type=wutils.valid_file_path, metavar='Basic_Statistics.txt',
                    help='Output of the Watchdog Module mergeStatistics applied on the Basic Statistics reported by FASTQC')
    qc.add_argument('--rawRegex', required=False, default=None, type=wutils.valid_regex_with_one_group, metavar='regex',
                    help='regular expression (python re) to extract the sample name from a fastq file with untrimmed reads')
    qc.add_argument('--trimRegex', required=False, default=None, type=wutils.valid_regex_with_one_group, metavar='regex',
                    help='regular expression (python re) to extract the sample name from a fastq file with trimmed reads')
    
    # options to summarize reads in BAM files (mapped)
    mp = parser.add_argument_group('Summarize read counts before and after mapping with ContextMap')
    mp.add_argument('--idxstatsSummary', required=False, default=None, type=wutils.valid_file_path, metavar='Idxstats.txt',
                    help='Output of the Watchdog Module mergeStatistics applied on the Idxstatistics reported by the bamstats module')
    mp.add_argument('--bamRegex', required=False, default=None, type=wutils.valid_regex_with_one_group, metavar='regex',
                    help='regular expression (python re) to extract the sample name from a bam file with mapped reads')
    mp.add_argument('--chromosomeGroupingTable', required=False, default=None, type=wutils.valid_file_path, metavar='groups.tsv',
                    help='tab-separated table with a header with chromosome names in column 0 and groups in column 1')
    
    # options to modify output location
    out = parser.add_argument_group('Output files')
    out.add_argument('--countTable', required=True, default=None, type=str, metavar='counts.tsv',
                     help='path for writing a table with all extracted read counts')
    out.add_argument('--countPlot', required=False, default=None, type=str, metavar='count_plot.[svg|pdf|png]',
                     help='path for saving a summary plot of total, trimmed and mapped reads'+
                     'format is identified by file ending, all formats supported by pyplot are allowed')
    out.add_argument('--groupPlot', required=False, default=None, type=str, metavar='count_plot.[svg|pdf|png]',
                     help='path for saving a summary plot of the fraction of mapped reads for given groups of chromosomes'+
                     'format is identified by file ending, all formats supported by pyplot are allowed')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def check_input_files(parser, cmdl_options):
    '''
    check manually combination of input files -> 5 valid combinations, 2 invalid combinations: chromosome table without idxstats 
    ensure all required regexes and output paths are given
    '''
    
    #input files
    if(cmdl_options.idxstatsSummary is None and cmdl_options.chromosomeGroupingTable is not None):
        parser.error('Missing input file: --idxstatsSummary is required for --chromosomeGroupingTable')
    if(cmdl_options.basicStatsSummary is None and cmdl_options.idxstatsSummary is None):
        parser.error('Missing input files: either --basicStatsSummary or --idxstatsSummary is required')
        
    # regex without file or file without regex
    if(cmdl_options.idxstatsSummary is not None):
        if(cmdl_options.bamRegex is None):
            parser.error('Missing regex: --bamRegex is required for --idxstatsSummary')
    else:
        if(cmdl_options.bamRegex is not None):
            parser.error('--bamRegex is set without --idxstatsSummary')
    if(cmdl_options.basicStatsSummary is None):
        if(cmdl_options.rawRegex is not None):
            parser.error('Missing regex: --rawRegex is required for --basicStatsSummary')
        if(cmdl_options.trimRegex is not None):
            parser.error('Missing regex: --trimRegex is required for --basicStatsSummary')
    else:
        if(cmdl_options.rawRegex is None):
            parser.error('--rawRegex is set without --basicStatsSummary')
        if(cmdl_options.trimRegex is None):
            parser.error('--trimRegex is set without --basicStatsSummary')
            
    # plotting options
    if(cmdl_options.countPlot is not None):
        if(cmdl_options.idxstatsSummary is None):
            parser.error('Missing input file: --idxStats is required for --countPlot')
        if(cmdl_options.basicStatsSummary is None):
            parser.error('Missing input file: --basicStatsSummary is required for --countPlot')
    if(cmdl_options.groupPlot is not None):
        if(cmdl_options.idxstatsSummary is None):
            parser.error('Missing input file: --idxStats is required for --groupPlot')
        if(cmdl_options.chromosomeGroupingTable is None):
            parser.error('Missing input file: --chromosomeGroupingTable is required for --groupPlot')
    

def get_and_create_outfiles(cmdl_options):   
    ''' creates all parent folders for writing the output '''
    wutils.make_parent_dirs(cmdl_options.countTable)
    # parent directories for optional plots
    if(cmdl_options.countPlot is not None):
        wutils.make_parent_dirs(cmdl_options.countPlot)
    if(cmdl_options.groupPlot is not None):
        wutils.make_parent_dirs(cmdl_options.groupPlot)


if __name__ == '__main__':
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # check options and prepare output locations
    p,o = get_command_line_options()
    check_input_files(p,o)
    get_and_create_outfiles(o)
    
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to Sophies first watchdog module!')
    print('This module summarizes and visualizes number of reads of sequencing data.\n')
    
    # run everything
    tab = algo_mappingsum.calculate_read_numbers(o.basicStatsSummary, o.rawRegex, o.trimRegex, o.idxstatsSummary, o.bamRegex, o.chromosomeGroupingTable)
    tf.writeTable(tab, o.countTable, sep='\t', header=True)
    if(o.countPlot is not None):
        plot_mappingsum.barplot_read_numbers(tab, o.countPlot)
    if(o.groupPlot is not None):
        plot_mappingsum.barplot_chrom_groups(tab, o.groupPlot)
    
    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n ->check out \''+o.countTable+'\' for the results')

    