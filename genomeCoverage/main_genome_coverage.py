'''
Created on 10 Nov 2016

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse
import watchdog_utils as wutils

from call_genome_coverage import bam_to_bedgraph, bedgraph_to_tdf


def get_command_line_options():
    ''' commandline parser for the genomeCoverage watchdog module that for converting bam to bedgraph or tdf format'''

    parser = argparse.ArgumentParser(description='wrapper for bedtools genomecov (bam -> bedgraph) and igvtools toTDF (bedgraph->tdf)')
    
    # input
    parser.add_argument('--bam', required=True, default=None, type=wutils.valid_bam,
        metavar='mapping.bam', help='path to bam file whose genome coverage should be analyzed, sam format is not supported,'
        'an index is not required for the bam file, but the reads have to be sorted by genomic coordinates')
    parser.add_argument('--genome', required=False, default=None, type=wutils.valid_igv_genome_file,
        metavar='organism.genome', help='genome file created by the IGV for the genome that was used to create the bam file,'
        'required only for conversion to tdf format, file has to end with .genome or .chrom.sizes depending on the format')
    
    # output
    parser.add_argument('--outPrefix', required=True, default=None, type=str,
        metavar='folder/prefix', help='file name prefix for saving the bedgraph and the tdf file')
    
    # boolean for conversion to tdf format.
    # --tdf flag does not change anything but is required for watchdog, only the value of --notdf is accessed
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--tdf', required=False, action='store_true', default=True,
        help='transform bedgraph file into tdf format using igvtools (default behaviour)')
    group.add_argument('--notdf', required=False, action='store_false', default=True,
        help='do not transform bedgraph file into tdf')
    
    # paths to executables
    parser.add_argument('--bedtoolsPath', required=False, default='bedtools', type=wutils.valid_exec,
        metavar='bedtools', help='path to bedtools executable, use if bedtools is not in PATH')
    parser.add_argument('--igvtoolsPath', required=False, default='igvtools', type=wutils.valid_exec,
        metavar='igvtools', help='path to igvtools executable, use if igvtools is not in PATH')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def check_option_interactions(parser, options):
    ''' resolve boolean option for tdf creation and interaction of tdf and genome option'''
    
    tdf_flag = options.notdf
    if tdf_flag and options.genome is None:
        parser.error('Missing --genome option for creating a tdf file (flag --tdf is set)')
    return tdf_flag
    

def create_outfiles(options):
    ''' creates parent directories and full filenames for the given file name prefix'''
    
    # create parent directories
    wutils.make_parent_dirs(options.outPrefix)
    
    # full output paths
    bedgraph_file = options.outPrefix+'.bedgraph'
    tdf_file = options.outPrefix+'.bedgraph.tdf'
    
    return bedgraph_file, tdf_file



def main():
    ''' main method for converting bam to bedgraph and tdf '''
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # check options and prepare output location
    p,o = get_command_line_options()
    create_tdf = check_option_interactions(p,o)
    bedgraph_out, tdf_out = create_outfiles(o)
    
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the genomeCoverage Module!\n')
    
    # perform conversion bam -> bedgraph -> tdf 
    bam_to_bedgraph(o.bam, bedgraph_out, bedtools_exec=o.bedtoolsPath)
    if create_tdf:
        bedgraph_to_tdf(bedgraph_out, o.genome, tdf_out, igvtools_exec=o.igvtoolsPath)
    
    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!')
    if create_tdf:
        print('-> the tdf file is located at \''+tdf_out+'\'')
    else:
        print('-> the bedgraph file is located at \''+bedgraph_out+'\'')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)


if __name__ == '__main__':
    main()
