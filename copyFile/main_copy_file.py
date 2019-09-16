'''
Created on Jul 3, 2018

@author: friedl
'''

import os.path
import shutil
import argparse
import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')

import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for the copyFile watchdog module'''
    
    parser = argparse.ArgumentParser(description='Copies a given file to a new location.')
    
    parser.add_argument('--sourcePath', required=True, default=None, type=wutils.valid_file_path, metavar='path/to/file',
                        help='path of the file to copy')
    parser.add_argument('--targetPath', required=True, default=None, type=str, metavar='new/path/to/file',
                        help='path of the new location of the file, all non-existing parent directories of the file are created')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def copy_file(source_path, target_path):
    ''' copies file of 'source_path' to 'target_path' (folder and filename) '''
    
    # check if file to copy really exists
    if(os.path.exists(source_path) and os.path.isfile(source_path)):
        # create directory of the target location
        wutils.make_parent_dirs(target_path)
        # copy file
        shutil.copy(source_path, target_path)
        

if __name__ == '__main__':
    
    # parse command line
    print('Program call:')
    print(' '.join(sys.argv))
    p,o = get_command_line_options()
    
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Copy file '+o.sourcePath+' ...')
    
    copy_file(o.sourcePath, o.targetPath)
    
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Successfully copied -> check out '+o.targetPath)
    