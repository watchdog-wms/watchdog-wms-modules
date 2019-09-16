'''
Created on 30 Dec 2016

@author: friedl
'''

import re
import gzip

from contextlib import contextmanager


class FastaWriter:
    '''
    writes non-rrna reads in fasta format
    '''
    
    transition = str.maketrans('ATGC', 'TACG')
    
    def __init__(self, fastaFile, pair_nr=None):
        '''
        passes path of fasta file
        '''
        
        self.fasta_file = fastaFile
        self.pair_nr = pair_nr
        # valid options for pair_nr
        if(not(pair_nr is None or pair_nr==1 or pair_nr==2)):
            raise ValueError('Invalid pair number attribute '+str(pair_nr))
        # count written reads
        self.read_counter=0
    
    @contextmanager
    def open(self):
        ''' context manager for opening (before yield) and closing (after yield) the file handle'''
        
        self.fastaWriter=open(self.fasta_file, 'wt')
        yield self
        self.fastaWriter.close()
    
    def addRead(self, read):
        '''
        adds a read to the fasta file
        '''
        
        # read name and sequence
        readName = read.query_name
        nuclSeq = read.query_sequence
                        
        # reverse and complement read sequence if read is mapped to the reverse strand
        if(read.is_reverse):
            nuclSeq=nuclSeq[::-1]
            nuclSeq = nuclSeq.translate(self.transition)
                    
        # write read to fasta file: for paired reads append /1 or /2 to read names -> "ContextMap format"
        line ='>'+readName
        if(self.pair_nr is not None):
            line+='/'+str(self.pair_nr)
        line+='\n'+nuclSeq+'\n'
        self.fastaWriter.write(line)
        
        # update written reads
        self.read_counter+=1



class FastQWriter():
    '''
    writes non-rrna reads in fastq format
    '''
    
    def __init__(self, fastqFileOut, fastqFileIn):
        '''
        passes path to fastq files
        '''
        
        self.fastqOut = fastqFileOut
        self.readsToWrite = set()
        self.fastqIn = fastqFileIn
        self.openfnc = open
        # count written (pos 0) and input (pos 1) reads
        self.read_counter=[0,0]
        
        
    def addRead(self, read):
        '''
        adds a read to the fastq file
        '''
        
        self.readsToWrite.add(read.query_name)
        
    @contextmanager
    def open(self):
        '''
        contextmanger for writing the fastq file
        '''
        
        yield self
        with self.openfnc(self.fastqOut, 'wt') as fastQWriter, open(self.fastqIn, 'rt') as fastQReader:
        
            # iterate over lines in groups of 4 (name, sequence, +, qualities)
            for line in fastQReader:
                fq_record = [line.strip('\n')]+[fastQReader.readline().strip('\n') for _ in range(0,3)]
                # count reads in input
                self.read_counter[1]+=1 
                    
                # extract read name -> remove additional information: information after space or read number appended via /1 and /2
                readName = fq_record[0].split(' ')[0][1:]
                readName = re.sub('/[12]$', '', readName)
                
                # write read if it is non-rrna
                if(readName in self.readsToWrite):
                    fastQWriter.write('\n'.join(fq_record)+'\n')
                    self.read_counter[0]+=1
                    self.readsToWrite.remove(readName)
        
        # check if all non-rrna reads were found in the fastq file
        if(len(self.readsToWrite)>0):
            if(len(self.readsToWrite)>=10):
                raise ValueError('I\'m sorry, something went wrong & this might be a bug...\nCould not write reads for:\n'+','.join(list(self.readsToWrite)[:10])+',...')
            else:
                raise ValueError('I\'m sorry, something went wrong & this might be a bug...\nCould not write reads for:\n'+','.join(list(self.readsToWrite)))



class CompressedFastQWriter(FastQWriter):
    
    def __init__(self, fastqFileOut, fastqFileIn):
        '''
        passes path to fastq files
        '''
        
        self.fastqOut = fastqFileOut
        self.readsToWrite = set()
        self.fastqIn = fastqFileIn
        self.openfnc = gzip.open
        # count written and input reads
        self.read_counter=[0,0]
        