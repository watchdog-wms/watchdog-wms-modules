'''
Created on 30 Nov 2017

@author: friedl
'''

import module_utils.circ_rna_utils as cutil

def read_prediction_file(infile):
    ''' reads in a prediction file an returns a dictionary coordinates -> prediciton '''
    
    # mapping of circRNA coordinates to circRNA prediction results
    res = {}
    
    with open(infile, 'rt') as circ_reader:
        
        # parse and check the file header
        header = circ_reader.readline().strip('\n')
        cutil.CircRNAPredictionFileFormat.check_file_format(header, infile)

        for line in circ_reader:
            line = line.strip('\n')
            prediction = cutil.CircRNAPredictionFileFormat.line_to_circRNAPrediction(line)
            if(prediction is not None):
                if(not prediction.circ_coordinates in res):
                    res[prediction.circ_coordinates]=prediction
                else:
                    raise ValueError('Duplicate circRNA '+str(prediction.circ_coordinates))
        
        return res


def combine_circular_rna(infile1, infile2, intersectfile, unionfile, intersected_union_count_file, minreads=2):
    ''' main method of the algorithm, callable from main or other scripts '''
    
    # read input files: mapping coordinates->prediction
    tool1 = read_prediction_file(infile1)
    tool2 = read_prediction_file(infile2)
    
    # open output files
    with open(intersectfile, 'wt' ) as iw, open(unionfile, 'wt') as uw, open(intersected_union_count_file, 'wt') as iuw:
        iw.write('\t'.join(cutil.CircRNAPredictionFileFormat.header_format)+'\n')
        uw.write('\t'.join(cutil.CircRNAPredictionFileFormat.header_format)+'\n')
        iuw.write('\t'.join(cutil.CircRNAPredictionFileFormat.header_format)+'\n')
        
        #iterate over all coordinates
        for circ in set(tool1.keys()) | set(tool2.keys()):
            
            # union & intersection
            if(circ in tool1 and circ in tool2):
                union = tool1[circ].combine_with(tool2[circ], 'union')
                if(union.circ_read_count>=minreads):
                    uw.write(cutil.CircRNAPredictionFileFormat.circprediction_to_line(union)+'\n')
                    iuw.write(cutil.CircRNAPredictionFileFormat.circprediction_to_line(union)+'\n')
                intersect = tool1[circ].combine_with(tool2[circ], 'intersection')
                if(intersect.circ_read_count>=minreads):
                    iw.write(cutil.CircRNAPredictionFileFormat.circprediction_to_line(intersect)+'\n')
                
            # union file only
            elif(circ in tool1):
                if(tool1[circ].circ_read_count>=minreads):
                    uw.write(cutil.CircRNAPredictionFileFormat.circprediction_to_line(tool1[circ])+'\n')
            else:
                if(tool2[circ].circ_read_count>=minreads):
                    uw.write(cutil.CircRNAPredictionFileFormat.circprediction_to_line(tool2[circ])+'\n')



if __name__ == '__main__':
    pass