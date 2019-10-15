'''
Created on 28 Nov 2016

@author: friedl
'''

import re

class ReadFilter():
    '''
    class for deciding if a read is removed from or kept in the mapping
    '''


    def __init__(self, maxEditDistance, maxMismatches, maxIndel):
        '''
        Constructor
        '''
        
        self.editdistance = maxEditDistance
        self.mismatches = maxMismatches
        self.indels = maxIndel
    
    
    def isRRNA(self, read):
        '''
        method returns true if read is kept in mapping
        method returns false if read is removed from the mapping
        the decision is based on maximum allowed edit distance, mismatches and indels
        '''
        
        # unmapped reads are removed
        if(read.is_unmapped):
            return False
        
        # supplementary mapping: read was already mapped elsewhere -> this is a duplicate
        elif(read.is_supplementary):
            return False
        
        # mapped read: apply filtering criteria
        else:
            
            # remove read if it is only partially mapped
            if(re.search('[NSHP]', read.cigarstring)):
                return False
            
            # no threshold for editdistance: editdistance = None -> do nothing, else remove read if edit distance is too large
            if(self.editdistance is not None):
                if(getEditDistance(read)>self.editdistance):
                    return False
            
            # no threshold for mismatches: mismatches = None -> do nothing, else remove read if too many mismatches
            if(self.mismatches is not None):
                if(getMismatchCount(read)>self.mismatches):
                    return False
            
            # no threshold for indels: indels = None -> do nothing, else remove read if too many indels
            if(self.indels is not None):
                if(getIndelCount(read)>self.indels):
                    return False
            
            return True                       
        
        
        
def getEditDistance(read):
    '''
    calculates edit distance of an aligned read, returns none if information about edit distance is missing
    '''
    
    if(read.has_tag('NM')):
        return read.get_tag('NM')
    else:
        return None

def getMismatchCount(read):
    '''
    calculates mismatch count of an aligned read, returns none if information about mismatches is missing
    '''
    
    if(read.has_tag('MD')):
        mismatchString = read.get_tag('MD')
        return len(re.findall('\D+', mismatchString))-len(re.findall('\^', mismatchString))
    else:
        return None

def getIndelCount(read):
    '''
    calculates indel count (sum of indel lengths) of an aligned read
    '''
    
    if(read.cigartuples is not None):
        indelCount=0
        for el in read.cigartuples:
            # I or D part of cigar string
            if(el[0]==1 or el[0]==2):
                indelCount+=el[1]
        return indelCount
    else:
        return None