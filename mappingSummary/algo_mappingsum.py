'''
Created on Mar 27, 2018

@author: friedl
'''

import txttables.tablefunc as tf
from txttables.tableclass import Table

from collections import defaultdict
import re


def calculate_read_numbers(fastqc_summary, raw_regex, trim_regex, idx_summary, bam_regex, chrom_file):
    '''
    produces a table with read counts: columns = count types, rows = samples
    fastqc_summary: summary file of fastqc
    raw_regex: regular expression to match raw fastq files, group 1 -> samplename
    trim_regex: regular expression to match trimmed fastq files, group 1 -> samplename
    idx_summary: summary filw of idxstats on bamfiles
    bam_regex: regular expression to match bam files, group 1 -> samplename
    chrom_file: tab-separated table with chromosome name (col 0) and organism/group name (col 1)
    '''
    
    # table to return
    res_tab = None
    
    #process fastqc
    fastqc_tab =None
    if(fastqc_summary is not None):
        fastqc_tab = parse_fastqc(fastqc_summary, raw_regex, trim_regex)
    
    # process idxstats
    idx_tab =None
    if(idx_summary is not None):
        idx_tab = parse_idx(idx_summary, bam_regex)
    
    # join results from first 2 steps
    if(fastqc_tab is None):
        res_tab = idx_tab
    elif(idx_tab is None):
        res_tab = fastqc_tab
    else:
        res_tab = tf.joinTables(fastqc_tab, idx_tab, [(0,0)], joinType='inner')
        
    #process chromosome groupings
    if(chrom_file is not None):
        chrom_tab = parse_chrom_file(chrom_file, idx_summary, bam_regex)
        res_tab = tf.joinTables(res_tab, chrom_tab, [(0,0)], joinType='inner')
        
    return res_tab



def parse_fastqc(fastqc_summary, raw_regex, trim_regex):
    '''
    produces a table with read counts: columns = count types -> raw and trimmed, rows = samples
    fastqc_summary: summary file of fastqc
    raw_regex: regular expression to match raw fastq files, group 1 -> samplename
    trim_regex: regular expression to match trimmed fastq files, group 1 -> samplename
    '''
    
    # map sample name to list (raw, trimmed) counts
    sample_to_counts=defaultdict(lambda:[0,0])
    
    # read fastqc statistics and update dictionary
    qcTable = tf.readTable(fastqc_summary, sep='\t', header=True)
    for rowInd in range(0, qcTable.rowNum()):
        if(qcTable.get(rowInd, 0)=='Total Sequences'):
            # get filename of the row and compare it against the regular expression for trimmed and raw fastq files to find the sample name
            filename = qcTable.get(rowInd, 2)
            readcount = int(qcTable.get(rowInd, 1))
            raw_match = re.search(raw_regex,filename)
            trim_match = re.search(trim_regex,filename)
            # add readcount if samplename is found
            if(raw_match):
                sample_to_counts[raw_match.group(1)][0]+=readcount
            elif(trim_match):
                sample_to_counts[trim_match.group(1)][1]+=readcount
    
    # resulting table
    counts = Table()
    counts.addColumn(str, 'sample', None)
    counts.addColumn(int, 'raw', 0)
    counts.addColumn(int, 'trimmed', 0)
    
    # transform content of the dictionary into the resulting table
    for sample, readnrs in sorted(sample_to_counts.items()):
        counts.addRow([sample, readnrs[0], readnrs[1]])

    return counts



def parse_idx(idx_summary, bam_regex):
    '''
    produces a table with read counts: columns = count types -> 1 column for mapped reads, rows = samples
    idx_summary: summary filw of idxstats on bamfiles
    bam_regex: regular expression to match bam files, group 1 -> samplename
    '''
    
    # map sample name to list counts
    sample_to_mapped=defaultdict(lambda:0)
    
    #read idxstat table and update dictionary
    idxTable = tf.readTable(idx_summary, sep='\t', header=True)
    for row in range(0, idxTable.rowNum()):
        bam_match = re.search(bam_regex, idxTable.get(row, 4))
        if(bam_match):
            sample = bam_match.group(1)
            count = int(idxTable.get(row, 2))
            sample_to_mapped[sample]+=count
        
    # resulting table
    counts = Table()
    counts.addColumn(str, 'sample', None)
    counts.addColumn(int, 'mapped', 0)
    
    # transform content of the dictionary into the resulting table
    for sample, readnr in sorted(sample_to_mapped.items()):
        counts.addRow([sample, readnr])

    return counts


  
def parse_chrom_file(chrom_file, idx_summary, bam_regex):
    '''
    produces a table with read counts: columns = mapped reads for chromosome groups, rows = samples
    idx_summary: summary filw of idxstats on bamfiles
    bam_regex: regular expression to match bam files, group 1 -> samplename
    chrom_file: tab-separated table with chromosome name (col 0) and organism/group name (col 1)
    '''
    
    # read assignment of chromosomes to groups as dictionary
    chrom_to_group = {}
    chr_tab = tf.readTable(chrom_file, sep='\t', header=True, headerstart='#')
    for row in range(0, chr_tab.rowNum()):
        chrom_to_group[chr_tab.get(row, 0)]=chr_tab.get(row, 1)
        
    # map sample name to map group -> counts
    sample_to_mapped=defaultdict(lambda:defaultdict(int))
    
    #read idxstat table and update dictionary
    idxTable = tf.readTable(idx_summary, sep='\t', header=True)
    for row in range(0, idxTable.rowNum()):
        bam_match = re.search(bam_regex, idxTable.get(row, 4))
        chr_name = idxTable.get(row,0)
        if(bam_match and chr_name in chrom_to_group):
            sample = bam_match.group(1)
            count = int(idxTable.get(row, 2))
            sample_to_mapped[sample][chrom_to_group[chr_name]]+=count
            
    # resulting table
    counts = Table()
    counts.addColumn(str, 'sample', None)
    groups = sorted(list(set(chrom_to_group.values())))
    for val in groups:
        counts.addColumn(int, val, 0)
    
    # transform content of the dictionary into the resulting table
    for sample, readnr_dict in sorted(sample_to_mapped.items()):
        newRow=[sample]
        for val in groups:
            newRow.append(readnr_dict[val])
        counts.addRow(newRow)

    return counts

    