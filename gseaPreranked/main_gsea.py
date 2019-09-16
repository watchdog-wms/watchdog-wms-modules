'''
Created on Mar 8, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse
from os.path import join

import watchdog_utils as wutils
import utils.enrichment_analysis as ea


def get_command_line_options():
    '''command line parser for the GSEAPreranked wrapper running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Watchdog Module for running GSEAPreranked')
    # path to executable
    parser.add_argument('--gseaJar', required=True, default=None, type=wutils.valid_file_path, metavar='gsea.jar',
                        help='Path of the GSEA jar file')
    
    # name of the analysis
    parser.add_argument('--label', required=True, default=None, type=str, metavar='name',
                        help='name of the analysis, e.g. sample name')
    
    # output dir -> generate if necessary
    parser.add_argument('--outdir', required=True, default=None, type=str, metavar='outdir',
                        help='directory to store the results of GSEA')
    
    # options for the input table: location and format
    input_group = parser.add_argument_group('Input table')
    input_group.add_argument('--geneTab', required=True, default=None, type=wutils.valid_file_path, metavar='gene.tsv',
                             help='tab-separated table of genes with expression values/changes')
    header_gr = input_group.add_mutually_exclusive_group()
    header_gr.add_argument('--hasHeader', action='store_true', default=False,
                             help='indicates if the first line of the geneTab should be interpreted as header')
    header_gr.add_argument('--nohasHeader', action='store_true', default=False,
                             help='indicates if the first line of the geneTab should be interpreted as header')
    input_group.add_argument('--geneCol', required=False, default=0, type=wutils.positive_integer_or_zero, metavar='geneColPos',
                             help='0-based position of the column with gene names')
    input_group.add_argument('--rankCol', required=False, default=1, type=wutils.positive_integer_or_zero, metavar='rankColPos',
                             help='0-based position of the column with values to rank the genes, e.g. fold changes')
    
    # options for running gsea
    gsea_group = parser.add_argument_group('GSEA Prerank Options')
    gsea_group.add_argument('--geneset', required=False, default='hallmark', type=str,
                            choices=['go', 'hallmark', 'transcription_factor', 'oncogenic_signatures', 'immunologic_signatures'],
                            help='gene sets to test for enrichment, default: hallmark sets')
    gsea_group.add_argument('--genesetVersion', required=False, default='6.1', type=str, metavar='version',
                            help='version of MSigDB to use')
    gsea_group.add_argument('--scoring', required=False, default='unweighted', type=str, choices=['weighted', 'unweighted'],
                            help='unweighted: classic score based on ranks, weighted: score includes values used for ranking')
    gsea_group.add_argument('--plotNr', required=False, default=50, type=wutils.positive_integer, metavar='plot_nr',
                            help='create plots for "plot_nr" top scoring genes')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def get_and_create_outfiles(cmdl_options):   
    ''' 
    creates all parent folders for writing the output
    '''
        
    wutils.make_parent_dirs(join(cmdl_options.outdir,''))


if __name__ == '__main__':
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # process command-line options: generate output folder if it does not exist
    p,o = get_command_line_options()
    get_and_create_outfiles(o)
    
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Starting GSEA wrapper...\n')
    
    ### execution GSEA
    # create rank file
    ranktab = join(o.outdir, o.label+'.rnk')
    ea.create_gsea_from_tsv(o.geneTab, ranktab, o.hasHeader, o.geneCol, o.rankCol)
    
    # call gsea
    ea.run_gsea(ranktab, o.label, o.outdir, o.geneset, o.genesetVersion, o.scoring, o.gseaJar, o.plotNr)
    ### GSEA end
    
    # log
    end_timepoint=wutils.get_current_time()
    print('\n'+end_timepoint[1]+'GSEA finished, check out \''+o.outdir+'\' for the results')
    
    
