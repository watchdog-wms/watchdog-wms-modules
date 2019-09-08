'''
Created on 1 Dec 2017

@author: friedl
'''
import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse

import algo_circ_combination
import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for circRNA combination running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Combines the predictions of circularRNAs made with the modules for CIRI2 and circRNA_finder.')
    # paths for the 2 input files: two different predictions of circRNAs (CIRI2 & circRNA_finder)
    parser.add_argument('--inCircs1', required=True, default=None, type=wutils.valid_file_path, metavar='prediction1.circs.txt',
                        help='First prediction file with circRNAs and junction reads')
    parser.add_argument('--inCircs2', required=True, default=None, type=wutils.valid_file_path, metavar='prediction2.circs.txt',
                        help='Second prediction file with circRNAs and junction reads')
    # paths for the 3 output files: intersection and union of the outputs and intersection of coordinates with union of reads
    parser.add_argument('--outIntersection', required=True, default=None, metavar='intersection.circs.txt',
                        help='output path for the intersection of the predictions')
    parser.add_argument('--outUnion', required=True, default=None, metavar='union.circs.txt',
                        help='output path for the union of the predictions')
    parser.add_argument('--outIntersectedUnion', required=True, default=None, metavar='intersect_with_union_counts.circs.txt',
                        help='output path for the union of the predictions')
    # cut off for junction read count
    parser.add_argument('--minReads', required=False, default=2, type=wutils.positive_integer, metavar='cutoff',
                        help='Minimum number of predicted junction reads required for writing a circRNA into the output files,'+
                        'cutoff is applied independently to the intersection and the union of the predictions')
    
    # for watchdog compatibility
    parser.add_argument('--returnFilePath', required=False, default=None, metavar='watchdog_variables',
                        help='internal watchdog command line parameter')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def create_output_dirs(options):
    ''' create all directories for writing the combined circRNA predictions'''
    
    wutils.make_parent_dirs(options.outUnion)
    wutils.make_parent_dirs(options.outIntersection)
    wutils.make_parent_dirs(options.outIntersectedUnion)

def main():
    ''' main method for the commandline interface of the circ_combination module '''
    
    print('Program call:')
    print(' '.join(sys.argv))
    
    _,o = get_command_line_options()
    create_output_dirs(o)
    
    t_start = wutils.get_current_time()
    print(t_start[1]+'Start Combination of circRNA predictions')
    
    algo_circ_combination.combine_circular_rna(o.inCircs1, o.inCircs2, o.outIntersection, o.outUnion, o.outIntersectedUnion, o.minReads)
    
    t_end = wutils.get_current_time()
    print(t_end[1]+'Finished. Check out the results in '+o.outIntersection+', '+o.outIntersection+' and '+o.outIntersectedUnion)
    
    wutils.print_resources(t_end[0]-t_start[0], child_processes=False)
    
    if(o.returnFilePath is not None):
        returnVars = [('circIntersection', o.outIntersection), ('circUnion', o.outUnion), ('circIntersectedUnion', o.outIntersectedUnion)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)
    

if __name__ == '__main__':
    main()

        


    