'''
Created on May 18, 2018

@author: friedl
'''

import gzip

def concatenate_files(infile_list, outfile):
    '''
    concatenates files from infile_list and saves the result in outfile
    for each file a new line is started (-> newline at file end not required)
    'infile_list': list of file paths, order of the files = order of concatenation
    'outfile': path to save the concatenated file
    '''
    
    with open(outfile, 'wt') as concatwr:
        for infile in infile_list:
            
            # decide if infile is compressed or not
            if(infile.endswith('.gz')):
                open_func = gzip.open
            else:
                open_func = open
            
            with open_func(infile, 'rt') as reader:
                for line in reader:
                    concatwr.write(line.strip('\n')+'\n')