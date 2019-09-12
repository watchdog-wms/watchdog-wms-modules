'''
Created on 5 Dec 2017

@author: friedl
'''

# for handling both BAM and SAM format
from pysam import AlignmentFile # @UnresolvedImport

import re
import sys

import circ_rna_utils as cutils

def remove_linearly_mappable_reads(circ_rna_prediction_file, mapping, updated_prediction_file, minReads=2, paired=True):
    ''' main steps of the algorithm '''
    
    prediction = read_circ_rna_prediction(circ_rna_prediction_file)
    circ_reads = build_circ_read_dictionary(prediction)
    identify_linear_reads(mapping, circ_reads)
    updated_prediction = update_predictions(prediction, circ_reads, minReads, paired)
    write_circ_rna_prediction(updated_prediction, updated_prediction_file)
    
    
def read_circ_rna_prediction(circ_rna_file):
    ''' 
    reads a file with predictions of circRNAs into a list of circRNE predictions
    'circ_rna_file': formatted file from the CIRI2, circRNA_finder or combine_circ module (CircRNAPredictionFileFormat)
    return: a list of CircRNPrediction objects
    '''
    
    # list of predictions
    res = []
    
    with open(circ_rna_file, 'rt') as r:
        # check the file header
        header = r.readline().rstrip('\n')
        cutils.CircRNAPredictionFileFormat.check_file_format(header, circ_rna_file)
        for line in r:
            line = line.rstrip('\n')
            circ_pred = cutils.CircRNAPredictionFileFormat.line_to_circRNAPrediction(line)
            # line is in valid format
            if(circ_pred is not None):
                res.append(circ_pred)       
        
    return res

def write_circ_rna_prediction(prediction_list, output_file):
    '''
    writes the predicted circRNAs back to a file
    'prediction_list': a list of CircRNPrediction objects
    'output_file': path for writing the predictions in the CircRNAPredictionFileFormat
    '''
    
    with open(output_file, 'wt') as w:
        # write the file header
        w.write(cutils.CircRNAPredictionFileFormat.get_file_header()+'\n')
        # write a line for each prediction
        for circ_pred in prediction_list:
            w.write(cutils.CircRNAPredictionFileFormat.circprediction_to_line(circ_pred)+'\n')

    
def build_circ_read_dictionary(prediction_list):
    '''
    initializes the data structure for recording linearly mappable circRNA reads
    'prediction_list': a list of CircRNPrediction objects
    return: a dictionary mapping read ids of the CircRNPrediction objects to a count (set to initial value 0)
    '''
    
    # mapping readID -> count of linearly mapped mates
    res = {}
    
    # iteration over all reads in all predicions
    for prediction in prediction_list:
        for read_id in prediction.circ_id_list:
            # add each id with count 0
            if(not read_id in res):
                res[read_id]=0
    
    return res

def identify_linear_reads(mapping, read_dictionary):
    '''
    extracts mappings of putative circular read pairs and counts the linearly mapped mates for each pair
    'mapping': a SAM or BAM file (produced by ContextMap)
    'read_dictionary': contains all read ids of putative circular pairs and a counter set to 0 for the linearly mappable mates of each pair
        -> read_dictionary is updated in place with the linear read counts from the mapping
    '''
    
    # check SAM or BAM format
    open_mode = 'r'
    if(re.search('.bam$', mapping, re.IGNORECASE)):
        open_mode+='b'
    
    # iterate over the mappings of puptative circular reads 
    with AlignmentFile(mapping, open_mode) as af:
        for mr in af:
            if(mr.query_name in read_dictionary):
                
                # designed for mappings with ContextMap -> each read aligned once, for other mappings: consider primary non-supplementary alignments only
                if (not mr.is_secondary and not mr.is_supplementary):
                
                    # linear read = mapped read without soft or hard clipped regions
                    if((not mr.is_unmapped) and re.search('^([0-9]+[MNID])+$', mr.cigarstring)):
                        read_dictionary[mr.query_name]+=1
                

def update_predictions(prediction_list, read_dictionary, minReads, paired):
    '''
    removes read pairs from circular RNA predictions with a linear mapping
    'prediction_list': a list of CircRNPrediction objects
    'read_dictionary': mapping of the read pair ids of the CircRNPrediction objects in 'prediction_list' to the number of linearly mappable mates of the pair
    'minReads': minimum number of supporting circular read pairs required for predicting a circRNA
    'paired': if reads are paired-end data
    return: an updated list of CircRNPrediction objects,
            read pairs with two linearly mappable mates are removed and
            circRNA predictions with less then minReads are removed
    '''
    
    # new list with the updated predictions fulfilling the minReads threshold
    updated_list = []
    
    # iterate over the predictions and update read id lists
    for prediction in prediction_list:
        new_read_id_list = []
        for read_id in prediction.circ_id_list:
            
            # paired-end data: at least one mate of the pair is not mapped linearly
            if paired:
                if(read_dictionary[read_id]<2):
                    new_read_id_list.append(read_id)
                # check for error
                elif(read_dictionary[read_id]>2):
                    sys.stderr.print('WARNING: More than 2 primary linear mappings for '+str(read_id)+'\n')
                    
            # single-end data: the read is not mapped linearly
            else:
                if(read_dictionary[read_id]<1):
                    new_read_id_list.append(read_id)
                elif(read_dictionary[read_id]>1):
                    sys.stderr.print('WARNING: More than 1 primary linear mapping for '+str(read_id)+'\n')
        
        # adapt read counts to new size of read id lists
        prediction.update_circ_reads(new_read_id_list)
        
        # check minReads threshold
        if(prediction.circ_read_count>=minReads):
            updated_list.append(prediction)
    
    return updated_list

