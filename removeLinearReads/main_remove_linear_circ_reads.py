'''
Created on 5 Dec 2017

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse
import watchdog_utils as wutils
import algo_remove_linear_circ_reads as algo

def get_command_line_options():
    '''command line parser for the removeLinearReads watchdog module'''
    
    parser = argparse.ArgumentParser(description='Removes linearly mappable reads from a circRNA prediction.')
    
    # file paths
    parser.add_argument('--mapping', required=True, default=None, type=wutils.valid_mapping, metavar='mapping.[sam|bam]',
                        help='path to a SAM or BAM file with mapped reads from the sample for which circRNAs were predicted')
    parser.add_argument('--circRNAPrediction', required=True, default=None, type=wutils.valid_file_path, metavar='circs.txt',
                        help='predicted circRNAs from the CIRI2, circRNAfinder or the circCombination module')
    parser.add_argument('--circOut', required=True, default=None, metavar='filtered_circs.txt',
                        help='all circRNAs from the input file with at least minReads remaining circular reads '+
                        'after removing all linearly mappable reads from the lists circular junction reads')
    
    # required read pairing
    parser.add_argument('--paired', required=True, default=True, type=wutils.valid_string_boolean, metavar='yes|no',
                        help='indicates if SAM or BAM input file contains paired-end (yes) or single-end (no) data')
    
    # optional cutoff
    parser.add_argument('--minReads', required=False, default=2, type=wutils.positive_integer, metavar='cutoff',
                        help='Minimum number of predicted junction reads required for writing a circRNA to the outputfile, default:2')
    
    # for watchdog compatibility
    parser.add_argument('--returnFilePath', required=False, default=None, metavar='watchdog_variables',
                        help='internal watchdog command line parameter')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def main():
    ''' main method for calling removeLinearReads from commandline '''
    
    print('Program call:')
    print(' '.join(sys.argv))
    
    _,o = get_command_line_options()
    wutils.make_parent_dirs(o.circOut)
    
    t_start = wutils.get_current_time()
    print(t_start[1]+'Start Removal of Linearly Mappable Reads')
    
    algo.remove_linearly_mappable_reads(o.circRNAPrediction, o.mapping, o.circOut, o.minReads, o.paired)
    
    t_end = wutils.get_current_time()
    print(t_end[1]+'Finished. Check out the results in '+o.circOut)
    
    wutils.print_resources(t_end[0]-t_start[0], child_processes=False)
    
    if(o.returnFilePath is not None):
        returnVars = [('filteredCircs', o.circOut)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)


if __name__ == '__main__':
    main()

    
    
