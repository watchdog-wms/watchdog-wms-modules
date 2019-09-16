'''
Created on May 18, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import argparse

import algo_concatenate
import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for file concatenation running as watchdog module'''
    
    parser = argparse.ArgumentParser(description='Concatenates 2 or more files.')
    #input files passed as comma-separated list
    parser.add_argument('--inFile', required=True, default=None, metavar='file1.txt,file2.txt,...', type=wutils.valid_list_of_files,
                        help='comma-separated list of inputfiles')
    # paths for the output files
    parser.add_argument('--outFile', required=True, default=None, metavar='concatenated_file.txt',
                        help='path to save the concatenated files')
    # for watchdog compatibility
    parser.add_argument('--returnFilePath', required=False, default=None, metavar='watchdog_variables',
                        help='internal watchdog command line parameter')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def create_output_dirs(cmdl_options):
    ''' create parent directories of concatenated output file '''
    
    wutils.make_parent_dirs(cmdl_options.outFile)
    

if __name__ == '__main__':
    
    print('Program call:')
    print(' '.join(sys.argv))
    p,o = get_command_line_options()
    create_output_dirs(o)
    
    t_start = wutils.get_current_time()
    print(t_start[1]+'Welcome! Concatenation of files in progress...')
    
    algo_concatenate.concatenate_files(o.inFile, o.outFile)
    
    t_end = wutils.get_current_time()
    print(t_end[1]+'Finished. Check out the results in '+o.outFile)
    
    wutils.print_resources(t_end[0]-t_start[0], child_processes=False)
    
    if(o.returnFilePath is not None):
        returnVars = [('concatenatedFile', o.outFile)]
        wutils.write_watchdog_return_file(o.returnFilePath, returnVars)
        
        