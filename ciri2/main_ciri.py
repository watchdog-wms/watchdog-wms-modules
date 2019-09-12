'''
Created on 15 Nov 2017

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse

import run_ciri as runner
import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for the CIRI2 wrapper running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Runs CIRI2 to detect circular RNAs in single-end or paired-end sequencing data.')
    
    # options for the input files
    input_group = parser.add_argument_group('Required Input')
    input_group.add_argument('--inReads1', required=False, default=None, type=wutils.valid_file_path, metavar='reads_1.fastq',
                        help='path to first fastq file with reads (for single-end or paired-end data first reads)')
    input_group.add_argument('--inReads2', required=False, default=None, type=wutils.valid_file_path, metavar='reads_2.fastq',
                        help='path to second fastq file with reads (for paired-end data second reads only)')
    input_group.add_argument('--inSAM', required=False, default=None, type=wutils.valid_file_path, metavar='mapped_reads.sam',
                        help='path to SAM file that was created with BWA Mem')
    input_group.add_argument('--reference', required=True, default=None, type=wutils.valid_file_path, metavar='genome.fa',
                        help='path to (multi-)fasta file with the reference genome')
    
    # options for the output files
    output_group = parser.add_argument_group('Output Location')
    output_group.add_argument('--outPrefix', required=True, default=None, metavar='out/prefix',
                        help='path and file name prefix for all files produced by this module,'+
                        'the final file is named out/prefixciriCirc.txt')
    output_group.add_argument('--outCirc', required=False, default=None, metavar='ciriCirc.txt',
                        help='final output of predicted CircRNAs')
    output_group.add_argument('--returnFilePath', required=False, default=None, metavar='watchdog_variables',
                        help='internal watchdog command line parameter')

    # options for BWA
    bwa_group = parser.add_argument_group('Optional Arguments for BWA')
    bwa_group.add_argument('--bwaPath', required=False, default='bwa', type=wutils.valid_exec, metavar='path/to/bwa',
                        help='specify a path to the BWA executable if bwa is not part of your PATH variable')
    bwa_group.add_argument('--bwaThreads', required=False, default=1, type=wutils.positive_integer, metavar='threadNr',
                        help='number of threads to use with BWA, default:1')
    bwa_group.add_argument('--bwaIndex', required=False, default=None, type=wutils.valid_bwa_index, metavar='index_prefix',
                        help='BWA index for the reference genome provided by the --reference option, if no index is provided it is automatically created by the module')
    bwa_group.add_argument('--bwaSeedSize', required=False, default=19, type=wutils.positive_integer, metavar='seed_length',
                        help='BWA -k parameter for the minimum seed length')
    bwa_group.add_argument('--bwaScoreThreshold', required=False, default=30, type=wutils.positive_integer, metavar='score',
                        help='BWA -T parameter for the minimum alignment score, default 30, but 19 recommended')
    
    # options for CIRI2
    ciri_group = parser.add_argument_group('Optional Arguments for CIRI2')
    ciri_group.add_argument('--ciriPath', required=False, default='CIRI2.pl', type=wutils.valid_file_path, metavar='CIRI2.pl',
                        help='path to CIRI2 perl script')
    ciri_group.add_argument('--ciriThreads', required=False, default=1, type=wutils.positive_integer, metavar='threadNr',
                        help='number of threads to use for CIRI2, default:1')
    ciri_group.add_argument('--ciriAnnotation', required=False, default=None, type=wutils.valid_file_path, metavar='annotation.gtf',
                        help='GTF file with gene annotations for the genome given in the --reference option,'
                        'if a GTF file is passed to this module, CIRI annotates all circRNAs with the corresponding gene')
    ciri_group.add_argument('--ciriStringency', required=False, default='high', type=str, choices=['high', 'medium', 'low'],
                        help='Controls how stringent CIRI filters the circRNAs based on circular reads, cigar strings and false positive reads')
    boolparam = ciri_group.add_mutually_exclusive_group()
    boolparam.add_argument('--ciriKeepTmpFiles', required=False, action='store_true', default=False,
                        help='if this flag is set, CIRI2 does not delete the temporary files at the end')
    boolparam.add_argument('--nociriKeepTmpFiles', required=False, action='store_false', default=True,
                        help='if this flag is set, CIRI2 does deletes the temporary files at the end (flag used only by watchdog)')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def check_inputfiles(parser, cmdl_options):
    ''' check manually required input files: either SAM file or FASTQ file (inReads1) are required '''

    if(cmdl_options.inSAM is None):
        if(cmdl_options.inReads1 is None):
            parser.error('Missing input file: --inReads1 is required!')
    else:
        if(cmdl_options.inReads1 is not None or cmdl_options.inReads2 is not None):
            parser.error('Argument clash: --inSAM and --inReads1/2 options are both set!')


def get_and_create_outfiles(cmdl_options):   
    ''' 
    creates all parent folders for writing the output
    '''
        
    wutils.make_parent_dirs(cmdl_options.outPrefix)
    if(cmdl_options.outCirc is not None):
        wutils.make_parent_dirs(cmdl_options.outCirc)


def main():
    ''' main method running the command line interface '''
    
    print('Program call:')
    print(' '.join(sys.argv))
    
    # parse command line
    p,o = get_command_line_options()
    check_inputfiles(p, o)
    get_and_create_outfiles(o)
    
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Starting CIRI2 wrapper...')
     
    # run BWA
    if(o.inSAM is None):
        # get index for BWA
        if(o.bwaIndex is None):
            index_prefix = o.outPrefix+'bwaIndex'
            runner.buildIndex(o.reference, index_prefix, o.bwaPath)
        else:
            index_prefix = o.bwaIndex
        # call bwa mem
        bwa_sam_file = o.outPrefix+'bwa.sam'
        runner.runBWA(o.inReads1, o.inReads2, index_prefix, bwa_sam_file, o.bwaThreads, o.bwaPath, o.bwaSeedSize, o.bwaScoreThreshold)
    else:
        bwa_sam_file = o.inSAM
    
    # run CIRI
    ciri_raw_output = o.outPrefix+'rawPred.txt'
    runner.runCIRI2(bwa_sam_file, ciri_raw_output, o.reference, o.ciriAnnotation, o.ciriThreads, o.ciriPath, o.ciriStringency, o.ciriKeepTmpFiles)
     
    # create the final output
    if(o.outCirc is None):
        final_output_file = o.outPrefix+'ciriCirc.txt'
    else:
        final_output_file = o.outCirc
    runner.writeOutput(ciri_raw_output, final_output_file)
    
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'CIRI2 finished, check out \''+final_output_file+'\' for the results')
     
    # output resources consumed by the script
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)
          
    # write watchdog return variables
    if(o.returnFilePath is not None):
        wutils.write_watchdog_return_file(o.returnFilePath, [('ciriCircs', final_output_file)])
    

if __name__ == '__main__':    
    main()
    
#     sys.argv=['main_ciri.py',
#               '--inReads1', '/mnt/raidinput/tmp/friedl/test_single_end_circ.fastq',
#               '--reference', '/mnt/raidproj/proj/projekte/sequencing/Illumina/HSV_circ_rnas/hg19/genome_hg19.fa',
#               '--outPrefix', '/mnt/raidinput/tmp/friedl/ciri_test/se_test_',
#               '--bwaIndex', '/mnt/raidproj/proj/projekte/sequencing/Illumina/HSV_circ_rnas/hg19/bwa/hg19',
#               '--bwaScoreThreshold', '19',
#               '--ciriAnnotation', '/mnt/raidproj/proj/projekte/sequencing/Illumina/HSV_circ_rnas/hg19/gene_annotation_hg19.gtf',
#               '--ciriPath', '/mnt/raidproj/proj/software/CIRI/CIRI_v2.0.6/CIRI2.pl']
      
    