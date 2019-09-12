'''
Created on 21 Nov 2016

@author: friedl
'''

import gzip
import subprocess
import os.path
import re

class BWAWrapper:
    '''
    wrapper class for executing bwa
    '''


    def __init__(self, readsIn1, readsIn2, index, pathToTemp, seedSize, threads, bwapath):
        '''
        Constructor passing all file paths for running bwa
        '''
        
        # pass input and output paths for bwa
        self.listInfiles =[readsIn1]
        if(readsIn2 is not None):
            self.listInfiles.append(readsIn2)
        self.indexPref = index
        self.outDir = pathToTemp
        
        # pass bwa options
        self.bwaPath = bwapath
        self.seed = seedSize
        self.numThreads = threads
        
        # lists of files created by this class
        self.readFiles=[]
        self.samFiles=[]
    
    
    def map(self):
        '''
        method for calling bwa on the reads of the files in the listInfiles attribute
        |listInfiles| = 1 -> single-end data, |listInfiles| = 2 -> paired-end data
        '''
        
        # consider the reads of each file separately
        for ifile in self.listInfiles:
            
            #unzip file with reads
            fqFile=self.unzip(ifile)
            self.readFiles.append(fqFile)
            
            # map fq file to rrna reference
            samFile = self.runBWA(fqFile)
            self.samFiles.append(samFile)
    
          
    def unzip(self, origFilePath):
        '''
        method for unzipping fq files, if necessary, returns the original path if file is uncompressed, else returns the path of the decompressed file
        '''
    
        # file is compressed -> unzip file
        if(origFilePath.endswith('.gz')):
            
            # get new file name + location
            print('Extracting file '+origFilePath+'\n')
            fqFile = os.path.join(self.outDir, re.sub('.gz$','', os.path.basename(origFilePath)) )
            
            # extract file
            with open(fqFile, 'wb') as g, gzip.open(origFilePath, 'rb') as f:
                for chunk in f:
                    g.write(chunk)
            return fqFile
        
        # file is uncompressed -> nothing to do
        else:
            return origFilePath

        
    def runBWA(self, fastqFile):
        '''
        method for calling bwa on the decompressed fastq files
        '''
        
        print('Mapping file '+fastqFile+'\n')
        
        # redirect stdout to file handle -> stout of bwa is written to a sam file
        samFile = os.path.join(self.outDir, re.sub('(.fastq$)|(.fq$)','', os.path.basename(fastqFile))+'.sam')
        with open(samFile, 'wt') as handleSam:
        
            # call bwa
            options=[self.bwaPath, 'mem', '-k', str(self.seed), '-t', str(self.numThreads), self.indexPref, fastqFile]
            print('BWA program call:\n'+' '.join(options)+'\n')
            subprocess.check_call(options, stdout=handleSam)
        
        # pass sam file on
        return samFile
        
            
            