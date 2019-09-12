'''
Created on Jul 27, 2018

@author: friedl
'''

'''
Created on Jul 3, 2018

@author: friedl
'''

import argparse
import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')

import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for the createFolder watchdog module'''
    
    parser = argparse.ArgumentParser(description='Creates a folder and its parent directories.')
    
    parser.add_argument('--folderPath', required=True, default=None, type=str, metavar='path/to/folder',
                        help='folder that will be created')

    cmdl_options = parser.parse_args()
    return parser, cmdl_options


def main():
    ''' main method executing the watchdog module'''
    
    # parse command line
    print('Program call:')
    print(' '.join(sys.argv))
    _,o = get_command_line_options()

    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Create folder '+o.folderPath+' ...')
    
    wutils.create_folder(o.folderPath)
    
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Folder creates finished successfully!')
    

if __name__ == '__main__':
    main()

    