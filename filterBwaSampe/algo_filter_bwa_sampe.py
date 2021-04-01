'''
Created on Jul 20, 2018

@author: friedl
'''

from pysam import AlignmentFile #@UnresolvedImport
import numpy as np

def check_sam_header(input_file):
    ''' tries to parse the header of the sam or bam file used as input '''
    
    m = get_mode_string(input_file, write=False)
    try:
        with AlignmentFile(input_file, m) as af:
            # no tag for the program that generated the output
            if not 'PG' in af.header:
                return False, None
            # only one program tag -> processed only by bwa
            if not len(af.header['PG'])==1:
                return False, None
            if not 'CL' in af.header['PG'][0]:
                return False, None
            # check progrma call
            return 'bwa sampe' in af.header['PG'][0]['CL'], None
    
    # catch any errors caused by pysam being unable to read the file
    except Exception as e:
        return False, str(e)


def get_mode_string(file, write=False):
    ''' auxiliary function for handling reading and writing in sam or bam format '''
    
    #write file
    if write:
        if(file.endswith('.bam')):
            return 'wb'
        else:
            return 'wh'
    #read file
    else:
        if(file.endswith('.bam')):
            return 'rb'
        else:
            return 'r'


def remove_reads(infile, outfile, unmapped, improper_pairs, map_qual_cut, opt_hit_cut, single_end):
    
    # first scan of the mapping -> decision of removal based on both reads of a pair
    # the decision for each read is saved in the filter object and passed on to the final filter step, no files are written
    
    first_removal_step = MappingFilter(infile, unmapped, improper_pairs, map_qual_cut, opt_hit_cut, single_end)
    first_removal_step.run()
    
    # second scan of the mapping -> reads are written to the output file 
    # which reads are written and which are not written is decided by the object created in the first filtering step that returns a final decision for every read
    # the second scan of the file is necessary to preserve the order of reads in the file (i.e. sorting by coordinates is preserved)
    
    second_removal_step = WritingFilter(infile, outfile, first_removal_step)
    second_removal_step.run()


class WritingFilter():
    '''
    class for removing read pairs from a sam file:
    it finally writes the output file based on the decision made by the mapping filter
    '''
    
    def __init__(self, infile, outfile, mapping_filter):
        # path to sam or bam files
        self.infile = infile
        self.outfile = outfile
        # decisions made for every read in infile
        self.mapping_filter = mapping_filter
    
    def run(self):
        '''
        scans the input file, gets the decision for each read and writes it to the output file if it meets the criteria
        '''
        
        read_mode = get_mode_string(self.infile, write=False)
        write_mode = get_mode_string(self.outfile, write=True)
 
        with AlignmentFile(self.infile, read_mode) as reader, AlignmentFile(self.outfile, write_mode, header=reader.header) as writer:
        
            for rpos, read in enumerate(reader):
                
                # query for each read the MappingFilter object to find out if the read meets all criteria
                # if the mapping filter returns true, the read fulfills the criteria and it is written to the output file
                # otherwise is is not written

                final_decision = self.mapping_filter.get_decision_for_read_at_pos(rpos)
                if(final_decision):
                    writer.write(read)

        
class MappingFilter():
    
    def __init__(self, infile, mapped_check, proper_pair_check, mapq_threshold, x0_threshold, single_end):
        self.mapping_file = infile
        self.mapped_check = mapped_check
        self.proper_pair_check = proper_pair_check
        self.mapq_threshold = mapq_threshold
        self.x0_threshold = x0_threshold
        self.single_end_data = single_end
        
        # numpy array of booleans: position of boolean = position of read in mapping file, boolean value = True -> keep read, False -> remove read
        self._decision_array = None 
        
        
    def run(self):
        '''
        evaluates mapping status, proper pair status, mapping quality and hit number according to the given thresholds
        '''
        
        # filters read pairs based on mapping status, proper pair stauts, mapping quality and hit number
        # the decision for each read is saved as boolean in _decision_array, the position in the array corresponds to the position in the file
        # decision = True -> read is kept, decision = false -> read is removed
        
        # initialize the boolean array for saving all decision
        # length = number of reads in the mapping
        
        read_mode = get_mode_string(self.mapping_file, write=False)
        with AlignmentFile(self.mapping_file, read_mode) as samReader1:
            total_reads = sum([1 for read in samReader1])
        self._decision_array = np.empty(total_reads, dtype=np.bool_)
        
        
        # dictionary mapping the read name to its position in the sam/bam file
        # a read name is inserted the dictionary if the first read of the pair is found, the key is the readname and the value is the position of the first read
        # the read name is removed from the dictionary if the second read of the pair is found
        name_to_mate_pos = {}
        
        with AlignmentFile(self.mapping_file, read_mode) as samReader2:
            
            for rpos, read in enumerate(samReader2):
                
                # check mapping status, proper pair status, mapping quality and hit number of current read without considering the mate, leading to a decision for the read
                # each criterion is evaluated separately and the final decision is a logical AND of the individual decisions
                # true: read is kept, false: read is removed
                # then make a consistent choice for both reads of a pair: logical AND of the decisions of the 2 reads
                # the decision of the mate is retrieved via name_to_mate_pos giving the position of the decision in _decision_array
                
                # decision based on the mapping status
                # first case -> mapping status criterion is not applied, second case -> read is unmapped/mapped
                if not self.mapped_check:
                    mapped_decision = True
                else:
                    mapped_decision = (not read.is_unmapped)
                    
                # decision based on the proper pair property (as defined in the output of bwa sampe)
                # first case -> improper pair criterion is not applied, second case -> evaluate mapped and proper pair flag
                # special case: proper_pair_flag = True but also unmapped = True, this happend when a read is aligned across two chromosomes (concatenated in bwa algo)
                if not self.proper_pair_check:
                    properpair_decision = True
                else:
                    properpair_decision = (not read.is_unmapped) and read.is_proper_pair
                
                # decision based on the mapping quality
                # first case: mapping quality is not evaluated, second case: keep read if the mapping quality is sufficiently high
                # unmapped reads have a mapping score of 0
                if self.mapq_threshold is None :
                    mapq_decision = True
                else:
                    mapq_decision = read.mapping_quality >= self.mapq_threshold
                    
                # decision based on the hit number
                # first case: hit number is not evaluated, second case: read is unmapped (no hit, is removed), third case: evaluate x0 tag
                # special case: mapped reads without X0 tag -> read with XT:A:M = read has no proper alignment but was aligned near its mate with many mismatches, 
                # this special case is considered as hit count = 1
                if self.x0_threshold is None :
                    hit_decision = True
                elif read.is_unmapped :
                    hit_decision = False
                else:
                    hit_decision = (not read.has_tag('X0')) or read.get_tag('X0') <= self.x0_threshold
                    
                # final decision for the read without considering the mate
                read_decision = mapped_decision and properpair_decision and mapq_decision and hit_decision
                self._decision_array[rpos] = read_decision
                
                # look for the mate
                # first case: mate not seen yet -> insert read name into the dictionary
                # second case: mate was already decided -> remove read name from dictionary and create the final decision for the pair as logical AND
                # of the decisions of both reads
                if read.query_name not in name_to_mate_pos :
                    name_to_mate_pos[read.query_name]= rpos
                else:
                    mate_pos = name_to_mate_pos.pop(read.query_name)
                    pair_decision = self._decision_array[mate_pos] and self._decision_array[rpos]
                    self._decision_array[rpos] = pair_decision
                    self._decision_array[mate_pos] = pair_decision
        
        # error handling -> this should not happen if each read is paired and each read of a pair occurs exactly once in the file
        # if the error occurs, a programming error occurred or the input file was not generated with bwa sampe or was already filtered by another tool
        if(not self.single_end_data and len(name_to_mate_pos)>0):
            print('WARNING: Could not find mates for some reads:\n'+','.join(list(sorted(name_to_mate_pos.keys()))[:3] )+',... ('+str(len(name_to_mate_pos))+' reads)')
            print('These reads will be removed from the final output')
            for _,rpos in name_to_mate_pos.items():
                self._decision_array[rpos]=False
            del name_to_mate_pos
    
    
    def get_decision_for_read_at_pos(self, position):
        return self._decision_array[position]
        


# if __name__ == '__main__':
#     
#     import time
#     
#     # test space required for boolean array
#     x = np.empty(10000000000, dtype=np.bool_)
#     for p in range(0,10000000000):
#         x[p]=True
#     time.sleep(20)
    
    
    
