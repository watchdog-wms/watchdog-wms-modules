'''
Created on 19 Dec 2016

@author: friedl
'''

import sys
import pysam

from rrna_filter_algo.nonrrnaWriter import FastaWriter, FastQWriter, CompressedFastQWriter

class SingleEndWriter:
    '''
    class for writing result of filtering for single-end data
    '''

    def __init__(self, rrnaMapping, readsIn, samOut, readsOut, readFilter):
        '''
        Constructor passing file paths and filter
        '''
        
        # pass file paths for input and output
        self.rrnaMapping = rrnaMapping
        self.samOut = samOut
        self.readsOut = readsOut
        
        # actual rrna filter
        self.readfilter = readFilter
        
        # initialize writer for non-rrna reads
        if(readsOut.endswith('.fa')):
            self.nonRrnaWriter=FastaWriter(readsOut)
        elif(readsOut.endswith('.fq') or readsOut.endswith('.fastq')):
            self.nonRrnaWriter=FastQWriter(readsOut,readsIn)
        elif(readsOut.endswith('.fq.gz')):
            self.nonRrnaWriter=CompressedFastQWriter(readsOut,readsIn)
            
        # count rrna reads
        self.rrna_counter = 0
        
    
    def writeResults(self):
        '''
        writes sam file with rrna reads and fastq/fasta file with non-rrna reads
        '''
        
        print('Filtering content of \n'+self.rrnaMapping+'\nWriting to \n'+self.samOut+'\n'+self.readsOut)
        
        # read in sam file with rrna mapping
        with pysam.AlignmentFile(self.rrnaMapping, 'r') as samReader, self.nonRrnaWriter.open():
        
            # generate handle and header for writing sam file with final rrna reads  
            headerOut=samReader.header.to_dict()
            headerOut['PG'].append({'ID':'1','CL':' '.join(sys.argv)})
            with pysam.AlignmentFile(self.samOut, 'wh', header=headerOut) as samWriter:
            
                # iterate over all reads and apply filter
                for read in samReader:
                    # supplementary mapping -> read is already mapped elsewhere, consider each read only once
                    if(not read.is_supplementary):
        
                        # filter read -> pass it to writer for non-rrna reads
                        if(not self.readfilter.isRRNA(read)):
                            self.nonRrnaWriter.addRead(read)
                        # keep read -> write it to sam file for rrna reads
                        else:
                            samWriter.write(read)
                            self.rrna_counter+=1
        
        
class PairedEndWriter:
    '''
    class for writing result of filtering for single-end data
    '''
    
    def __init__(self, rrnaMapping, readsIn, samOut, readsOut, readFilter, pairFilter):
        
        # lists of file paths (length 2, one path for reads1 and reads2)  
        self.rrnaMapping = rrnaMapping
        self.readsOut = readsOut
        
        # path to sam file for writing rrna reads
        self.samOut = samOut
        
        #read filter objects
        self.readFilter = readFilter
        self.pairFilter = pairFilter
        
        # initialize writer for non-rrna reads
        self.nonRrnaWriter = []
        for p, (inf,outf) in enumerate(zip(readsIn, readsOut)):
            if(outf.endswith('.fa')):
                self.nonRrnaWriter.append(FastaWriter(outf,pair_nr=p+1))
            elif(outf.endswith('.fq') or outf.endswith('.fastq')):
                self.nonRrnaWriter.append(FastQWriter(outf,inf))
            elif(outf.endswith('.fq.gz')):
                self.nonRrnaWriter.append(CompressedFastQWriter(outf,inf))
        
        # maps for rrna and non-rrna reads with mapping:  readname->(is_reverse, has_mate, is_unmapped, mapping_start, chromosome_id)
        self.rrnaToStrand = {}
        self.nonrrnaToStrand = {}
        
        # count rrna reads of each strand
        self.rrna_counter=[0,0]
        
    def writeResults(self):
        '''
        writes sam file with rrna reads and fastq/fasta file with non-rrna reads
        '''
        
        ### set up dictionaries with rrna and non-rrna reads from reads2
        print('Reading content of \n'+self.rrnaMapping[1])
        with pysam.AlignmentFile(self.rrnaMapping[1], 'r') as samReader:
            
            for read in samReader:
                # supplementary mapping -> read is already mapped elsewhere, consider each read only once
                if(not read.is_supplementary):
                    
                    # non-rrna read
                    if(not self.readFilter.isRRNA(read)):
                        self.nonrrnaToStrand[read.query_name]=(read.is_reverse, False, read.is_unmapped, read.reference_start, read.reference_id)
                    # rrna read
                    else:
                        self.rrnaToStrand[read.query_name]=(read.is_reverse, False, read.is_unmapped, read.reference_start, read.reference_id)
            
        
        ### filter and write reads1
        print('Filtering content of \n'+self.rrnaMapping[0]+'\nWriting to \n'+self.samOut+'\n'+self.readsOut[0])
        with pysam.AlignmentFile(self.rrnaMapping[0], 'r') as samReader, self.nonRrnaWriter[0].open():
        
            # generate handle and header for writing sam file with final rrna reads  
            headerOut=samReader.header.to_dict()
            headerOut['PG'].append({'ID':'1','CL':' '.join(sys.argv)})
            samWriter = pysam.AlignmentFile(self.samOut, 'wh', header=headerOut)
            
            # iterate over reads1 and write them
            for read in samReader:
                # supplementary mapping -> read is already mapped elsewhere, consider each read only once
                if(not read.is_supplementary):
                    
                    # check for rrna
                    isRrnaDecision=self.readFilter.isRRNA(read)
                    # ensure consistent filtering of pairs -> change classification of mate (reads2) if necessary
                    self.pairFilter.updatePairs(self.rrnaToStrand, self.nonrrnaToStrand, read.query_name, isRrnaDecision)

                    # non-rrna read pair
                    if(read.query_name in self.nonrrnaToStrand):
                        self.nonRrnaWriter[0].addRead(read)
                        self.nonrrnaToStrand[read.query_name]=(read.is_reverse, True, read.is_unmapped, read.reference_start, read.reference_id)
                    
                    # rrna read pair
                    elif(read.query_name in self.rrnaToStrand):
                        read.is_paired=True
                        read.is_proper_pair=True
                        read.is_read1 = True
                        read.is_read2 = False
                        read.mate_is_unmapped =self.rrnaToStrand[read.query_name][2]
                        read.mate_is_reverse = self.rrnaToStrand[read.query_name][0]
                        read.next_reference_id=self.rrnaToStrand[read.query_name][4]
                        read.next_reference_start=self.rrnaToStrand[read.query_name][3]
                        samWriter.write(read)
                        self.rrna_counter[0]+=1
                        self.rrnaToStrand[read.query_name]=(read.is_reverse, True, read.is_unmapped, read.reference_start, read.reference_id)
                    
                    #read is unpaired: none of the maps has an entry for the read
                    else:
                        # unpaired rrna read
                        if(isRrnaDecision):
                            read.is_paired=True
                            read.is_proper_pair=False
                            read.is_read1 = True
                            read.is_read2 = False
                            read.mate_is_unmapped =True
                            samWriter.write(read)
                            self.rrna_counter[0]+=1
                        # unpaired non-rrna read
                        else:
                            self.nonRrnaWriter[0].addRead(read)
        
        
        ### filter and write reads2
        print('Filtering content of \n'+self.rrnaMapping[1]+'\nWriting to \n'+self.samOut+'\n'+self.readsOut[1])
        with pysam.AlignmentFile(self.rrnaMapping[1], 'r') as samReader, self.nonRrnaWriter[1].open():
            for read in samReader:
                
                # supplementary mapping -> read is already mapped elsewhere, consider each read only once
                if(not read.is_supplementary):
                    
                    # non-rrna read
                    if(read.query_name in self.nonrrnaToStrand):
                        self.nonRrnaWriter[1].addRead(read)
                    
                    # rrna read
                    elif(read.query_name in self.rrnaToStrand):
                        
                        #read is paired
                        if(self.rrnaToStrand[read.query_name][1]):
                            read.is_paired=True
                            read.is_proper_pair=True
                            read.is_read1 = False
                            read.is_read2 = True
                            read.mate_is_unmapped =self.rrnaToStrand[read.query_name][2]
                            read.mate_is_reverse = self.rrnaToStrand[read.query_name][0]
                            read.next_reference_id=self.rrnaToStrand[read.query_name][4]
                            read.next_reference_start=self.rrnaToStrand[read.query_name][3]
                        # singleton read
                        else:
                            read.is_paired=True
                            read.is_proper_pair=False
                            read.is_read1 = False
                            read.is_read2 = True
                            read.mate_is_unmapped =True
                        samWriter.write(read)
                        self.rrna_counter[1]+=1
                    
                    #else not necessary because of step 1
        
        samWriter.close()
        
        