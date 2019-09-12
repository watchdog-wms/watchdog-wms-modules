'''
Created on 15 Nov 2017

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse
import os

import run_circrna_finder
import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for the circRNA_finder wrapper running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Runs circRNA_finder to detect circular RNAs in single-end or paired-end sequencing data.')
    
    # options for the input files
    input_group = parser.add_argument_group('Required Input')
    input_group.add_argument('--inReads1', required=False, default=None, type=wutils.valid_file_path, metavar='reads_1.fastq',
                        help='path to single-end fastq file or path to first fastq file with paired reads')
    input_group.add_argument('--inReads2', required=False, default=None, type=wutils.valid_file_path, metavar='reads_2.fastq',
                        help='path to second fastq file with paired reads (paired-end data only)')
    input_group.add_argument('--strandedLibrary', required=False, default=0, type=int, choices=[0,1,2],
                        help='indicates if the library is strand specific, 0 = unstranded/unknown, 1 = stranded (first read), 2 = stranded (second read), (default: 0),'
                        '\nif the library type is unstranded/unknown the strand is guessed from the strand of the AG-GT splice site')
    input_group.add_argument('--reference', required=False, default=None, type=wutils.valid_file_path, metavar='genome.fa',
                        help='path to (multi-)fasta file with the reference genome')
    input_group.add_argument('--inSTAR', required=False, default=None, type=wutils.valid_star_output, metavar='mapped_reads.sam',
                        help='output prefix of a STAR mapping that was created with STAR run with chimeric segment detection')
    
    # options for the output files
    output_group = parser.add_argument_group('Output Location')
    output_group.add_argument('--outPrefix', required=True, default=None, metavar='out/prefix',
                        help='path and file name prefix for all files produced by this module,'+
                        'the final file is named out/prefixcfCirc.txt')
    output_group.add_argument('--outCirc', required=False, default=None, metavar='cfCirc.txt',
                        help='final output of predicted CircRNAs')
    output_group.add_argument('--returnFilePath', required=False, default=None, metavar='watchdog_variables',
                        help='internal watchdog command line parameter')

    # options for STAR
    star_group = parser.add_argument_group('Optional Arguments for STAR')
    star_group.add_argument('--starPath', required=False, default='STAR', type=wutils.valid_exec, metavar='path/to/STAR',
                        help='specify a path to the STAR executable if STAR is not part of your PATH variable')
    star_group.add_argument('--starThreads', required=False, default=1, type=wutils.positive_integer, metavar='threadNr',
                        help='number of threads to use with STAR, default:1')
    star_group.add_argument('--starIndex', required=False, default=None, type=wutils.valid_star_index, metavar='index_folder',
                        help='STAR index for the reference genome, if no index is provided it is automatically created by the module using the file given by --reference')
 
    # options for circ_rna_finder
    circrna_finder_group = parser.add_argument_group('Optional Arguments for circRNA_finder')
    circrna_finder_group.add_argument('--cfPath', required=False, default='postProcessStarAlignment.pl', type=wutils.valid_file_path, metavar='postProcessStarAlignment.pl',
                        help='path to circRNA_finder perl script \'postProcessStarAlignment.pl\'')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def check_inputfiles(parser, cmdl_options):
    ''' check manually required input files: either one or two FASTQ files with an (indexed) reference genome or mapped reads with the chimeric junctions of STAR are required '''

    if(cmdl_options.inSTAR is None):
        if(cmdl_options.inReads1 is None):
            parser.error('Missing input file: --inReads1 is required!')
        if(cmdl_options.reference is None and cmdl_options.starIndex is None):
            parser.error('Missing reference genome: either --reference or --starIndex is required!')
        if(cmdl_options.reference is not None and cmdl_options.starIndex is not None):
            parser.error('Argument clash: --reference and --starIndex options are both set!')
    else:
        if(cmdl_options.inReads1 is not None or cmdl_options.inReads2 is not None):
            parser.error('Argument clash: --inSTAR and --inReads1/2 options are both set!')
        if(cmdl_options.reference is not None):
            parser.error('Unused argument --reference, --inSTAR option does not require a reference genome')
        if(cmdl_options.starIndex is not None):
            parser.error('Unused argument --starIndex, --inSTAR option does not require a star index')
            

def get_and_create_outfiles(cmdl_options):   
    ''' 
    creates all parent folders for writing the output
    '''
        
    wutils.make_parent_dirs(cmdl_options.outPrefix)
    if(cmdl_options.outCirc is not None):
        wutils.make_parent_dirs(cmdl_options.outCirc)

     
def main():
    ''' main method for the commandline interface of the circrna_finder module '''
    
    print('Program call:')
    print(' '.join(sys.argv))
    
    # handle command line options
    p,o = get_command_line_options()
    check_inputfiles(p, o)
    get_and_create_outfiles(o)
    
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Starting circRNA_finder wrapper')
    
    # call STAR
    if(o.inSTAR is None):
        
        # build STAR index
        if(o.starIndex is None):
            index_location = o.outPrefix+'star_index'
            os.mkdir(index_location)
            run_circrna_finder.build_STAR_index(o.reference, index_location, star_path=o.starPath, star_thread_nr=o.starThreads)
        else:
            index_location = o.starIndex
        
        # map reads with STAR
        star_file_prefix = o.outPrefix
        run_circrna_finder.map_with_STAR(o.inReads1, o.inReads2, index_location, star_file_prefix, star_path=o.starPath, star_thread_nr=o.starThreads)
    
    else:
        star_file_prefix = o.inSTAR
    
    # call circRNA_finder
    run_circrna_finder.run_circRNA_finder(star_file_prefix, o.outPrefix, circRNA_finder_path=o.cfPath)
    
    # create final output file
    if(o.outCirc is None):
        final_file = o.outPrefix+'cfCirc.txt'
    else:
        final_file = o.outCirc
    run_circrna_finder.annotate_and_wirte_output(o.outPrefix, star_file_prefix, final_file, o.strandedLibrary)
    
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'circRNA_finder finished, check out \''+final_file+'\' for the results')
    
    # output resources
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)
    
    # write watchdog return variables
    if(o.returnFilePath is not None):
        wutils.write_watchdog_return_file(o.returnFilePath, [('cfCircs', final_file)])



if __name__ == '__main__':
    main()

    #     sys.argv=['main_circrna_finder.py', '-h']
#     sys.argv=['main_circrna_finder.py',
#               '--inReads1', '/mnt/raidinput/tmp/friedl/test_single_end_circ.fastq',
#               '--starIndex', '/mnt/raidproj/proj/projekte/sequencing/Illumina/HSV_circ_rnas/hg19/star/',
#               '--outPrefix', '/mnt/raidinput/tmp/friedl/cf_test/se_test_',
#               '--cfPath', '/mnt/raidproj/proj/software/circRNA_finder/circRNA_finder-fbd522a8e8fadec27c01d03d839c2dda73b8a189/postProcessStarAlignment.pl']
# 
#     sys.argv=['main_circrna_finder.py',
#               '--inSTAR', '/mnt/raidinput/tmp/friedl/cf_test/se_test_',
#               '--outPrefix', '/mnt/raidinput/tmp/friedl/cf_test/se_test_',
#               '--cfPath', '/mnt/raidproj/proj/software/circRNA_finder/circRNA_finder-fbd522a8e8fadec27c01d03d839c2dda73b8a189/postProcessStarAlignment.pl',
#               '--strandedLibrary', '0']

#     sys.argv=['main_circrna_finder.py',
#               '--inReads1', '/mnt/raidinput/tmp/friedl/test_SRR2124299_1.fastq',
#               '--inReads2', '/mnt/raidinput/tmp/friedl/test_SRR2124299_2.fastq',
#               '--starIndex', '/mnt/raidproj/proj/projekte/sequencing/Illumina/HSV_circ_rnas/hg19/star/',
#               '--outPrefix', '/mnt/raidinput/tmp/friedl/cf_test_pe_r2/SRR2124299_',
#               '--cfPath', '/mnt/raidproj/proj/software/circRNA_finder/circRNA_finder-fbd522a8e8fadec27c01d03d839c2dda73b8a189/postProcessStarAlignment.pl',
#               '--strandedLibrary', '2']
