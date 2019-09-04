'''
Created on 4 Oct 2017

@author: friedl
'''

class FastqRecord(object):
    '''
    wrapper class for all information about a read given in a fastq file
    '''
    
    def __init__(self, readid, readname, sequence, qualities):
        self.readid = readid
        self.readname = readname
        self.sequence = sequence
        self.qualities = qualities
        
    def getReadID(self):
        ''' retrieves the read ID: the read ID is identical for reads of a pair'''
        return self.readid
    
    def getReadName(self):
        ''' retrieves the read name: the read name differs for reads of a pair and contains the read id '''
        return self.readname
    
    def getReadSequence(self):
        return self.sequence
    
    def getQualities(self):
        return self.qualities

  
    
def readRecord(fileHandle):
    ''' reads the next FastqRecord from a fastq file openend for reading (given as file handle 'fileHandle') '''
    # parse read name
    headerline = fileHandle.readline()
    # no more reads in the file
    if(headerline == ''):
        return None
    readname = headerline[1:].strip('\n')
    readid =readname.split(' ')[0]
    # parse read sequence
    sequence = fileHandle.readline().strip('\n')
    # skip line with + sign
    fileHandle.readline()
    # get per-base qualities
    qualities = fileHandle.readline().strip('\n')
    # generate a new record
    return FastqRecord(readid, readname, sequence, qualities)

def writeRecord(record, fileHandle):
    ''' writes a FastqRecord 'record' to a fastq file opened for writing (given as file handle 'fileHandle') '''
    fileHandle.write('@'+record.getReadName()+'\n')
    fileHandle.write(record.getReadSequence())
    fileHandle.write('\n')
    fileHandle.write('+\n')
    fileHandle.write(record.getQualities()+'\n')



def filterFastqs(in1, in2, out1, out2, remainder):
    '''
    algorithm to find read pairs in 2 fastq files giving first reads and second reads
    input files:
        'in1' fastq file with the first reads
        'in2' fastq file with the second reads
    output files:
        'out1' reads from 'in1' with a mate in 'in2'
        'out2' reads from 'in2' with a mate in 'in1'
        'remainder' reads from 'in1' and 'in2' without a matching mate
    '''
    
    with open(in1, 'rt') as reader1, open(in2, 'rt') as reader2, open(out1, 'wt') as writer1, open(out2, 'wt') as writer2:
        
        # data structures for the file handles
        readers = [reader1, reader2]
        writers = [writer1, writer2]
        
        # reads without mate from the two files: maps read id to read record object
        singletons = [{}, {}]
        
        # current read from file 1 and file 2
        cur_records = [readRecord(reader1), readRecord(reader2)]
        
        # variable indicating which record was updated -> is never set to True if the new read is None because of EOF (exception if both reads are None)
        record1_new = True
        
        # special case: first file is empty but second contains some reads
        if(cur_records[0] is None):
            record1_new = False
        
        # iterate until no reads are left
        while(cur_records[0] is not None or cur_records[1] is not None):
            
            # get the current read that was updated
            if(record1_new is True):
                new_pos=0
                old_pos=1
            else:
                new_pos=1
                old_pos=0
                
            # decide if the mate of the current read was already found
            
            # check new read against last read from the other file
            if(cur_records[old_pos] is not None and cur_records[old_pos].getReadID()==cur_records[new_pos].getReadID()):
                writeRecord(cur_records[old_pos], writers[old_pos])
                writeRecord(cur_records[new_pos], writers[new_pos])
                cur_records[old_pos] = readRecord(readers[old_pos])
                cur_records[new_pos] = readRecord(readers[new_pos])
                if(cur_records[0] is None):
                    record1_new = False
                else:
                    record1_new = True
            
            # check new read against the unpaired reads from the other file 
            else:
                if cur_records[new_pos].getReadID() in singletons[old_pos]:
                    writeRecord(singletons[old_pos][cur_records[new_pos].getReadID()], writers[old_pos])
                    writeRecord(cur_records[new_pos], writers[new_pos])
                    del singletons[old_pos][cur_records[new_pos].getReadID()]
                    cur_records[new_pos] = readRecord(readers[new_pos])
                    # end of file is reached
                    if(cur_records[new_pos] is None):
                        record1_new = not record1_new
                
                # mate of the current new read was not found -> add it as singleton
                else:
                    # update first read
                    record1_new = False
                    if(cur_records[0] is not None and (len(singletons[0])<=len(singletons[1]) or cur_records[1] is None)):
                        singletons[0][cur_records[0].getReadID()]=cur_records[0]
                        cur_records[0] = readRecord(readers[0])
                        record1_new = True
                        
                    # update second read if end of first file is reached or first read was not updated                        
                    if((cur_records[0] is None and cur_records[1] is not None) or record1_new is False):
                        singletons[1][cur_records[1].getReadID()]=cur_records[1]
                        cur_records[1] = readRecord(readers[1])
                        record1_new = False
                        if(cur_records[1] is None):
                            record1_new = not record1_new

    
    # write singleton reads to another fastq file
    with open(remainder, 'wt') as singletonwriter:
        for record in singletons[0].values():
            writeRecord(record, singletonwriter)
        for record in singletons[1].values():
            writeRecord(record, singletonwriter)

    