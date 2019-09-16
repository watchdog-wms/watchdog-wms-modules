'''
Created on 23 Dec 2016

@author: friedl
'''

class PairFilterBoth:
    '''
    read pair is classified as rrna if both read of the pair meet the filter criteria
    '''
    
    def __init__(self):
        pass
    
    def updatePairs(self, rrnaMap, nonrrnaMap, readname, filterdecision):
        '''
        readname: name of the current read
        filterdecision: true if read is rrna, false if read is not rrna based on the current read of the pair
        rrnaMap & nonrrnaMap: classification of the read pair based on the other mate of the pair
            -> rrnaMap contains readnames of rrna reads, nonrrnaMap contains read names of non-rrna reads
        maps are modified in place
        '''
        
        # should hopefully never happen
        if(readname in rrnaMap and readname in nonrrnaMap):
            raise ValueError('Duplicate classification for read '+str(readname))
        
        # current read is not rrna but mate is rrna -> pair is non-rrna
        if((not filterdecision) and readname in rrnaMap):
            nonrrnaMap[readname]=rrnaMap[readname]
            rrnaMap.pop(readname)

        
class PairFilterOne:
    '''
    read pair is classified as rrna if one read of the pair meets the filter criteria
    '''
    
    def __init__(self):
        pass
    
    def updatePairs(self, rrnaMap, nonrrnaMap, readname, filterdecision):
        '''
        readname: name of the current read
        filterdecision: true if read is rrna, false if read is not rrna based on the current read of the pair
        rrnaMap & nonrrnaMap: classification of the read pair based on the other mate of the pair
            -> rrnaMap contains readnames of rrna reads, nonrrnaMap contains read names of non-rrna reads
        maps are modified in place
        '''
        
        # should hopefully never happen
        if(readname in rrnaMap and readname in nonrrnaMap):
            raise ValueError('Duplicate classification for read '+str(readname))
        
        # current read is rrna but mate is not rrna
        if(filterdecision and readname in nonrrnaMap):
            rrnaMap[readname]=nonrrnaMap[readname]
            nonrrnaMap.pop(readname)
        
        