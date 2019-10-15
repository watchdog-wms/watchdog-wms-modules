'''
Created on Jul 3, 2018

@author: friedl
'''

import shutil
import os.path
import argparse
import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')

import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for the deleteFolder watchdog module'''
    
    parser = argparse.ArgumentParser(description='Deletes a folder and all its content.')
    
    parser.add_argument('--folder', required=True, default=None, type=wutils.valid_folder_path, metavar='path/to/folder',
                        help='path to the folder that will be deleted')
    
    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def delete_folder(folder_path):
    ''' deletes the path to the given folder if it exists'''
    
    if(os.path.exists(folder_path) and os.path.isdir(folder_path)):
        shutil.rmtree(folder_path)


if __name__ == '__main__':
    
    # parse command line
    p,o = get_command_line_options()
    print(' '.join(sys.argv))
    
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Deleting folder...')
    
    delete_folder(o.folder)
    
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Successfully deleted '+o.folder)
    
    
    