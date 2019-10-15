'''
Created on 21 Nov 2016

@author: friedl
'''

import tempfile
import shutil

class WorkDirManager:
    '''
    class for handling working directory implemented as context manager
    '''
    
    def __init__(self, workdir, keepTmp=False):
        '''
        constructor, executed as usually
        '''
       
        self.path = workdir
        self.delete_on_close = not keepTmp 
    
    def __enter__(self):
        '''
        executed after constructor at beginning of with, creates directory
        '''
        
        self.pathToTemp= tempfile.mkdtemp(suffix='', prefix='rrnaFilter_', dir=self.path)
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        '''
        executed after end of with block, deletes directory
        '''
        if(self.delete_on_close):
            shutil.rmtree(self.pathToTemp)
            