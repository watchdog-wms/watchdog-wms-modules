'''
Created on 21 Nov 2016

@author: friedl
'''  

from rrna_filter_algo.workDirManager import WorkDirManager
from rrna_filter_algo.bwaWrapper import BWAWrapper
from rrna_filter_algo.readFilter import ReadFilter
from rrna_filter_algo.writeFilterRes import SingleEndWriter, PairedEndWriter
from rrna_filter_algo.pairFilter import *


def run(in1, in2, rrnaIndex, workdir, out1, out2, sam,
        keeptmp=False, seedSize=25, threads=1, bwapath='bwa', maxEditDistance=None, maxMismatches=None, maxIndels=None, pairFiltering=2):
    
    # create temporary directory in workdir
    with WorkDirManager(workdir, keeptmp) as workingDirectory:
        print('Using working directory:'+workingDirectory.pathToTemp+'\n')
        
        # map reads to rrna sequence with bwa
        mappingStep = BWAWrapper(in1, in2, rrnaIndex, workingDirectory.pathToTemp, seedSize, threads, bwapath)
        mappingStep.map()
        
        # object that decides if a mapped read is rrna or not
        filterSettings = ReadFilter(maxEditDistance, maxMismatches, maxIndels)
          
        # find non rrna reads in single-end data
        if(in2 is None):
            filterStep = SingleEndWriter(mappingStep.samFiles[0], mappingStep.readFiles[0], sam, out1, filterSettings)
            filterStep.writeResults()
            # check inputted & outputted reads
            check_read_sum(mappingStep.readFiles[0], filterStep.rrna_counter, filterStep.nonRrnaWriter.read_counter)
              
        # find non-rrna reads in paired-end data
        else:
            
            # policy for handling pairs
            pairSettings = PairFilterBoth()
            if(pairFiltering==1):
                pairSettings = PairFilterOne()
              
            # filtering step
            filterStep = PairedEndWriter(mappingStep.samFiles, mappingStep.readFiles, sam, [out1, out2], filterSettings, pairSettings)
            filterStep.writeResults()
            
            # check inputted & outputted reads
            check_read_sum(mappingStep.readFiles[0], filterStep.rrna_counter[0], filterStep.nonRrnaWriter[0].read_counter)
            check_read_sum(mappingStep.readFiles[1], filterStep.rrna_counter[1], filterStep.nonRrnaWriter[1].read_counter)
        
        

def check_read_sum(infile, rrna_reads, non_rrna_reads):
    
    # output in fasta format (writer does not count reads in infile)
    if isinstance(non_rrna_reads, int):
        with open(infile, 'rt') as wcl:
            line_count = sum(1 for _ in wcl)
        in_read_count=int(line_count/4)
        rrna_out_count=non_rrna_reads
    
    # output in fastq format (writer counts reads in infile and outfile
    else:
        in_read_count=non_rrna_reads[1]
        rrna_out_count=non_rrna_reads[0]
        
    
    if(not in_read_count==(rrna_reads+rrna_out_count)):
        print('Input '+str(in_read_count))
        print('Rrna '+str(rrna_out_count))
        print('Non rrna '+str(non_rrna_reads))
        raise ValueError('I\'m sorry, something went wrong & this might be a bug...\nFound different number of reads in input and output for '+str(infile))


if __name__ == '__main__':
    
    import os.path 
    folder = '/home/users/friedl/work/rrnaFilterTestData/mock_mock/'
    index = '/mnt/biostor1/Data/Databases/GENOMES/Homo/GRCh38/indices/bwa/rDNA_thomas'
    infq = ['/home/users/friedl/work/rrnaFilterTestData/mock_mock/mock-mock_R1.fastq',
          '/home/users/friedl/work/rrnaFilterTestData/mock_mock/mock-mock_R2.fastq']
    out = [os.path.join(folder,n) for n in ['1.fq', '2.fq', 'rrna.sam']]
    
    run(infq[0], infq[1], index, folder, out[0], out[1], out[2])
    
    