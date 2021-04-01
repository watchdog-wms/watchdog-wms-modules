'''
Created on Jul 20, 2018

@author: friedl
'''

import argparse
import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')

import watchdog_utils as wutils
import algo_filter_bwa_sampe

def get_command_line_options():
    '''command line parser for filtering the output of bwa sampe (watchdog module)'''
    
    parser = argparse.ArgumentParser(description='Removes read pairs from sam/bam files created by bwa sampe')
    
    # input and output files
    parser.add_argument('--inSamBam', required=True, default=None, type=wutils.valid_mapping,
        metavar='in.sam/bam', help='path to mapped paired reads in sam or bam format (recognized by file ending)')
    parser.add_argument('--outSamBam', required=True, default=None, type=wutils.valid_outfile_ending_sam_bam,
        metavar='out.sam/bam', help='path to write remaining paired reads in sam or bam format (recognized by file ending)')
    # parameter for file paththat is only used for watchdog
    parser.add_argument('--returnFilePath', required=False, default=None, type=str,
        metavar='watchdogdir/tabfile', help='flag that is used by watchdog only for handling the task output')
    
    # boolean parameters
    unmapped_flag = parser.add_mutually_exclusive_group()
    unmapped_flag.add_argument('--removeUnmapped', required=False, default=True, action='store_true',
        help='use this flag to remove pairs with at least one unmapped read, default: true')
    unmapped_flag.add_argument('--noremoveUnmapped', required=False, default=False, action='store_true',
        help='use this flag to keep unmapped read pairs, default: false')
    properpair_flag = parser.add_mutually_exclusive_group()
    properpair_flag.add_argument('--removeImproperPairs', required=False, default=True, action='store_true',
        help='use this flag to remove pairs that are not properly paired according to bwa sampe, default: true')
    properpair_flag.add_argument('--noremoveImproperPairs', required=False, default=False, action='store_true',
        help='use this flag to keep read pairs that are not properly paired, default: false')
    isSingleEnd_flag = parser.add_mutually_exclusive_group()
    isSingleEnd_flag.add_argument('--isSingleEnd', required=False, default=False, action='store_true',
        help='use this flag to indicate that single end data should be filtered, default: false')
    isSingleEnd_flag.add_argument('--noisSingleEnd', required=False, default=True, action='store_true',
        help='use this flag to indicate that paired end data should be filtered, default: true')



    # integer values -> do not remove any reads based on qualities or hit number by default
    parser.add_argument('--removeMapqBelow', required=False, default=20, type=wutils.positive_integer_or_zero,
        metavar='minQuality', help='Remove all read pairs with at least one mate of mapping quality smaller than minQuality (taken from field "MAPQ" in SAM file), '
        'setting the option to 0 deactivates filtering based on mapping quality, default: 20')
    parser.add_argument('--removeMoreThanOptimalHits', required=False, default=1, type=wutils.positive_integer_or_zero,
        metavar='maxHits', help='Remove all read pairs with more than maxHits optimal alignment positions for at least one mate (based bwa aln specific tag "X0"), '
        'setting the option to 0 deactivates filtering based on hit number, default: 1 (=removes all multi-mapped reads)')

    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options

def process_option_interactions(options):
    ''' handle boolean watchdog parameters (modelled as mutually exclusive group) & filtering on hit count and quality '''
    
    removal_of_unmapped = not options.noremoveUnmapped
    removal_of_improper_pairs = not options.noremoveImproperPairs
    isSingleEndData = options.isSingleEnd
    
    if (not removal_of_unmapped) and removal_of_improper_pairs:
        print('WARNING: removal of improper pairs will also remove unmapped pairs\n')
    
    if options.removeMapqBelow>0 and not options.removeMoreThanOptimalHits==1 :
        print('WARNING: filtering on mapping quality >0 will also remove multi-mapped reads\n')
        
    if (options.removeMapqBelow>0 or options.removeMoreThanOptimalHits>0) and (not removal_of_unmapped):
        print('WARNING: filtering on mapping quality >0 or hit count >0 will also remove unmapped reads\n')
        
    quality_cutoff = options.removeMapqBelow
    if quality_cutoff==0 :
        quality_cutoff=None
        
    hit_cutoff = options.removeMoreThanOptimalHits
    if hit_cutoff==0:
        hit_cutoff=None
    
    return removal_of_unmapped, removal_of_improper_pairs, quality_cutoff, hit_cutoff, isSingleEndData


def check_input_file(options):
    ''' check if the module is able to process the input file, i.e. if it was created by bwa sampe '''
    
    passed, error_msg = algo_filter_bwa_sampe.check_sam_header(options.inSamBam)
    if error_msg is not None:
        raise ValueError('There was an error reading the input file -> there is an error in the sam or bam format\n'+error_msg)
    if not passed:
        print('WARNING: The input file was not created by bwa sampe, this module might not work properly on the input file\n')
    

def create_outfiles(options):
    ''' create all parent directories of the sam/bam file with the remaining reads '''
    
    wutils.make_parent_dirs(options.outSamBam)



def main():
    ''' main method for running the module'''
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # check options and prepare output locations
    _,o = get_command_line_options()
    re_unmapped, re_imppairs, qual_cut, hit_cut, single_end = process_option_interactions(o)
    check_input_file(o)
    create_outfiles(o)
    
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the module for fitlering the output of bwa sampe !\n')

    algo_filter_bwa_sampe.remove_reads(o.inSamBam, o.outSamBam, re_unmapped, re_imppairs, qual_cut, hit_cut, single_end)
    
    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n-> check out \''+o.outSamBam+'\' with the filtered read pairs')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=False)

    # write return file for watchdog
    if(o.returnFilePath is not None):
        returnVars = [('filteredPairs', o.outSamBam)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)

if __name__ == '__main__':
#     sys.argv=['', '--inSamBam', '/home/extproj/tmp/friedl/chip_1_bcl3_wt_2h.bam', '--outSamBam', '/home/extproj/tmp/friedl/test/filt_chip_1_bcl3_wt_2h.sam',
#               '--removeUnmapped', '--removeImproperPairs', '--removeMoreThanOptimalHits', '1', '--removeMapqBelow', '20']
    main()
    
