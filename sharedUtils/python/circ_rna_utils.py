'''
Created on 1 Dec 2017

@author: friedl
'''

import sys

class CircRNACoordinates(object):
    ''' hashable representation of circRNA coordinates '''
    
    def __init__(self, chromosome, start, end, strand):
        # chromosome name
        self.chromosome = chromosome
        # 0-based start included, 0-based end excluded [start; end) 
        self.start = start
        self.end = end
        # + or -
        self.strand = strand
    
    def __str__(self):
        return self.chromosome+':'+str(self.start)+'-'+str(self.end)+','+self.strand
    
    def __repr__(self):
        return self.__str__()
    
    def __hash__(self):
        return hash((self.chromosome, self.start, self.end, self.strand))
    
    def __eq__(self, other):
        if(not self.chromosome == other.chromosome):
            return False
        elif(not self.start == other.start):
            return False
        elif(not self.end == other.end):
            return False
        elif(not self.strand == other.strand):
            return False
        else:
            return True
        
    def copy(self):
        return CircRNACoordinates(self.chromosome, self.start, self.end, self.strand)
    
    
    @classmethod
    def line_to_circRNACoordinates(cls, line, coordinate_columns=[0,1,2,3], one_based=False):
        ''' transforms the first 4 columns of one line into a CircRNACoordinates object,
            coordinate_columns: parser assumes a tab-separated format, coordinate columns gives the positions of the chromosome, start, end and strand,
            one_based: option for handling one_based coordinate input
        '''
        
        # check input options -> coordinate columns = list of 4 positions
        if(not len(coordinate_columns)==4):
            raise ValueError('Invalid coordinate_columns '+str(coordinate_columns))
        
        # check tab-separated format -> all columns in coordinate columns should be covered
        data = line.split('\t')
        if(len(data)<=max(coordinate_columns)):
            sys.stderr.write('Cannot parse line: '+str(line))
            return None
        
        # chromosome name -> nothing to check
        chr_name = data[coordinate_columns[0]]
        
        # start position -> check and transform into 0-based coordinates if necessary
        str_start = data[coordinate_columns[1]]
        try:
            start = int(str_start)
        except ValueError:
            sys.stderr.write('Invalid start position '+str(str_start) +' in '+str(line)+'\n')
            return None
        if(one_based is True):
            start=start-1
        if(start<0):
            sys.stderr.write('Negative start position '+str(str_start)+' in '+str(line)+'\n')
            return None
        
        # end position -> check for valid integer and compare with start position, nothing to transform for one-based data
        str_end = data[coordinate_columns[2]]
        try:
            end = int(str_end)
        except ValueError:
            sys.stderr.write('Invalid start position '+str(str_end)+' in '+str(line)+'\n')
            return None
        if(end<=start):
            sys.stderr.write('End smaller than start in '+str(line)+'\n')
            return None
        
        # strand -> check for + or -
        strand = data[coordinate_columns[3]]
        if(not (strand=='+' or strand=='-')):
            sys.stderr.write('Invalid strand '+str(strand)+' in '+str(line)+'\n')
            return None
        
        coordinates = CircRNACoordinates(chr_name, start, end, strand)
        return coordinates
        
        
class CircRNAPrediction(object):
    ''' representation of the prediction results of a circRNA tool '''
    
    def __init__(self, circ_coordinates, circ_read_count, circ_id_list):
        self.circ_coordinates = circ_coordinates
        self.circ_read_count = circ_read_count
        self.circ_id_list = circ_id_list
    
    def combine_with(self, other_circ_prediction, mode):
        ''' create a new prediction object von two independent predicitons for the same circRNA on the same data set'''
        if(not self.circ_coordinates==other_circ_prediction.circ_coordinates):
            return None
        read_set1 = set(self.circ_id_list)
        read_set2 = set(other_circ_prediction.circ_id_list)
        if(mode=='union'):
            combined_ids = list(read_set1 | read_set2)
        elif(mode=='intersection'):
            combined_ids = list(read_set1 & read_set2)
        else:
            return None
        return CircRNAPrediction(self.circ_coordinates, len(combined_ids), combined_ids)
    
    def update_circ_reads(self, read_list):
        ''' replaces 'circ_id_list' with read_list and changes the circ_read_count according to the size of 'read_list' '''
        self.circ_id_list = list(set(read_list))
        if(not self.circ_read_count==len(self.circ_id_list)):
            self.circ_read_count = len(self.circ_id_list)
            
    def copy(self):
        ''' creates a copy of the current circRNA object '''
        return CircRNAPrediction(self.circ_coordinates.copy(), self.circ_read_count, self.circ_id_list.copy())
    
    
class CircRNAPredictionFileFormat(object):
    ''' methods for parsing and writing the file format for circ rna predictions'''

    header_format = ['chr', 'start', 'end', 'strand', '#junction_reads', 'junction_reads_ID']
    
    @classmethod
    def line_to_circRNAPrediction(cls, line):
        ''' transforms one line (without '\n' at the end) describing a predicted circRNA into a data structure '''
        
        data = line.split('\t')
        if(not len(data)==6):
            sys.stderr.write('Cannot parse line: '+str(line))
            return None
        
        # check and parse coordinates
        coordinates = CircRNACoordinates.line_to_circRNACoordinates(line)
        if coordinates is None:
            return None
        
        #check junction read count
        try:
            jcount = int(data[4])
        except ValueError:
            sys.stderr.write('Invalid junction read count '+str(data[4])+' in '+str(line)+'\n')
            return None
        if(jcount<1):
            sys.stderr.write('Invalid junction read count <=0 '+str(data[4])+' in '+str(line)+'\n')
            return None
        
        # check read id list
        jreads = data[5].split(',')
        # remove duplicates in the list
        jreads = list(set(jreads))
        if(not len(jreads)==jcount):
            sys.stderr.write('Junction read list does not match read count in '+str(line)+'\n')
            return None
        
        prediction = CircRNAPrediction(coordinates, jcount, jreads)
        
        return prediction

    @classmethod
    def check_file_format(cls, headerline, filename):
        ''' checks the expected header (without '\n' at the end) '''
        
        header = headerline.split('\t')
        if(not len(header)==6):
            raise ValueError('Invalid input file '+str(filename)+': Expecting 6 tab-separated columns but found '+str(len(header)))
        for col_pos, col_name in enumerate(cls.header_format):
            if(not header[col_pos]==col_name):
                raise ValueError('Invalid input file '+str(filename)+': Cannot find column '+str(col_name)+' at position '+str(col_pos))
    
    @classmethod
    def circprediction_to_line(cls, prediction):
        
        c = prediction.circ_coordinates
        output_coord ='\t'.join([c.chromosome, str(c.start), str(c.end), c.strand])
        quant = '\t'.join([str(prediction.circ_read_count), ','.join(prediction.circ_id_list)])
        return output_coord+'\t'+quant
    
    @classmethod
    def get_file_header(cls):
        return '\t'.join(cls.header_format)
    