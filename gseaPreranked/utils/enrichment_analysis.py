'''
Created on 4 Aug 2017

@author: friedl
'''

# module to run external processes
import subprocess as sp

# python standard libraries for accessing the services of the operating system
import os
import sys

# python module for regular expressions
import re

# numpy library for numeric operations
import numpy as np

# python plotting libraries
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as mpatches

# own package for manipulating tables
import txttables.tablefunc as tf
import txttables.tableclass as tc



def create_gsea_from_edgeR(inFile, outFile):
    '''
    transforms the results of edgeR from DETest module into input for GSEA preranked
    'inFile': edgeR results with log2FoldChanges for all tested genes
    'outFile': file ending with *.rnk with the gene names and log2 fold changes given in 'inFile'
    '''
    
    # edgeR file = tsv files with gene name in column 5 and fold change in column 1
    create_gsea_from_tsv(inFile, outFile, True, 5, 1)


def create_gsea_from_tsv(inFile, outFile, hasHeader, gene_pos, rank_pos):
    '''
    transforms the results of a gene expression analysis into input for GSEA preranked
    'inFile': tab-separated table with fold changes for all genes
    'hasHeader': boolean that indicates if the first line of 'inFile' is a table header
    'gene_pos': 0-based position of the column with the gene names
    'rank_pos': 0-based position of the column with the fold change or other values to rank the genes
    'outFile': file ending with *.rnk with the gene names and fold changes changes given in 'inFile'
    '''
    
    # extract column with gene name and log fold change
    edgeTab = tf.readTable(inFile, sep='\t', header=hasHeader, colsToRead=[gene_pos,rank_pos])
    
    # convert gene names to upper case -> match GSEA database
    edgeTab.modifyColumn(0, lambda s: s.upper())
    
    # sort by log fold change
    edgeTab.sortRows([1], [lambda x:float(x)], [False])
    
    # write file back to disk
    tf.writeTable(edgeTab, outFile, sep='\t', header=False)
    
    
def run_gsea(rankFile, label, outDir, gene_set='go', set_version='5.2', scoring_scheme='weighted', gseajar='/home/proj/software/GSEA/gsea2-2.2.3.jar', plot_nr=50):
    '''
    runs gsea analysis named 'label' for genes given in 'rankFile' and writes all results to 'outDir'
    'gene_set' specifies the kind of gene sets used in the analysis, currently 'go' and 'hallmark' are supported
    'scoring_scheme' defines how the enrichment score is calculated
        weighted -> sum factor*log2 foldchange (ranking metric) per gene
        unweighted -> sum constant factor per gene (independent of ranking metric)
    '''
    
    # call gsea on preranked set of genes
    command =['java', '-cp', gseajar, '-Xmx2G', 'xtools.gsea.GseaPreranked']
    
    # parameters to pass file paths and label
    command.extend(['-rnk', rankFile, '-rpt_label', label, '-out', outDir])
    
    # use gene names for the analysis
    command.extend(['-collapse', 'false'])
    
    # select gene sets to test: GO terms, hallmark gene sets, transcription factor sets, immunologic signature or cancer signatures
    if(gene_set=='go'):
        command.extend(['-gmx', 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c5.all.v'+set_version+'.symbols.gmt'])
    elif(gene_set=='hallmark'):
        command.extend(['-gmx', 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/h.all.v'+set_version+'.symbols.gmt'])
    elif(gene_set=='transcription_factor'):
        command.extend(['-gmx', 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c3.tft.v'+set_version+'.symbols.gmt'])
    elif(gene_set=='oncogenic_signatures'):
        command.extend(['-gmx', 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c6.all.v'+set_version+'.symbols.gmt'])
    elif(gene_set=='immunologic_signatures'):
        command.extend(['-gmx', 'gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/c7.all.v'+set_version+'.symbols.gmt'])
        
    # scoring scheme for calculating the enrichment score
    if(scoring_scheme=='weighted'):
        command.extend(['-scoring_scheme', 'weighted',])
    elif(scoring_scheme=='unweighted'):
        command.extend(['-scoring_scheme', 'classic',])
    else:
        raise ValueError('Unknown scoring scheme '+str(scoring_scheme))

    # remaining parameters: algorithm and speed
    command.extend(['-set_max', '500', '-set_min', '15', '-nperm', '1000'])
    
    # remaining parameters: output size and format
    command.extend(['-plot_top_x', str(plot_nr), '-zip_report', 'false'])
    
    # repitition of default parameters
    command.extend(['-mode', 'Max_probe', '-norm', 'meandiv', '-include_only_symbols', 'true'])
    command.extend(['-make_sets', 'true', '-rnd_seed', 'timestamp', '-gui', 'false'])
    
    sp.check_call(command)


def generate_gsea_overview_table(basefolder, outfile, test_names=None, writeExcel=False, sort_by_score=False, col_sort_func=lambda x:x):
    '''
    summarizes results of GSEA for several samples
        takes all gsea results available in 'basefolder'
        selects all gene sets with fdr <0.05 in at least one run of gsea
        creates a summary table 'outfile' with the enrichment scores and fdrs for all gene sets in all runs
        (row = gene sets with fdr<0.05 for at least one sample, column = NES or FDR per gsea run)
    writeExcel = True: creates a summary table in xls format
    sort_by_score: flag to control the sorting order of gene sets
        True <-> sort by number of samples with significant FDR in descending order and in case of ties by sum of absolute enrichment scores in descending order
        False <-> sort lexicographically by gene set names
    '''
    
    # handle list of folders with gsea results from differently labeled tests
    if(isinstance(basefolder,str)):
        basefolder = [basefolder]
        
    # map: testname (compared conditions) -> table with merged NES and FDR value
    table_map = {}
    #read in all results tables given in the subdirectories of basefolder
    for bfolder in basefolder:
        for file in os.listdir(bfolder):
            content = file.split('.')
            if(content[1]=='GseaPreranked'):
                testname = content[0]
                if(test_names is None or testname in test_names):
                    tablist=[]
                    for direction in ['neg', 'pos']:
                        tabpath= os.path.join(bfolder, file, 'gsea_report_for_na_'+direction+'_'+content[2]+'.xls')
                        tab = tf.readTable(tabpath, sep='\t', header=True, colsToRead=['NAME', 'NES', 'FDR q-val'])
                        tab.changeColumnName('NES', 'NES_'+direction)
                        tab.changeColumnName('FDR q-val', 'FDR_'+direction)
                        tablist.append(tab)
                    # join results of negative and positive enrichment -> one column with NES values and one column with pvalue
                    final_nes_colname = 'NES_'+testname
                    final_fdr_colname = 'FDR_'+testname
                    merge_negpos = tf.joinTables(tablist[0], tablist[1], joinCols=[(0,0)], joinType='fullouter')
                    merge_negpos.addColumn(str, columnName=final_nes_colname, defaultValue='0')
                    merge_negpos.modifyColumn(final_nes_colname, _merge_nes_scores, wholeRow=True)
                    merge_negpos.addColumn(str, columnName=final_fdr_colname, defaultValue='1')
                    merge_negpos.modifyColumn(final_fdr_colname, _merge_fdr, wholeRow=True)
                    merge_negpos = tf.selectColumns(merge_negpos, colList=['NAME', final_nes_colname, final_fdr_colname])      
                    table_map[testname] = merge_negpos
    
    # get significant gene sets: q-value < 0.05
    sign_gene_sets=set()
    for testname, tab in table_map.items():
        for r in range(0, tab.rowNum()):
            if(float(tab.get(r, 'FDR_'+testname))<0.05):
                sign_gene_sets.add(tab.get(r, 'NAME'))

    # set up result table with significant gene sets
    res =tc.Table()
    res.addColumn(str, columnName='Gene Set')
    for gs in sorted(sign_gene_sets):
        res.addRow([gs])
    
    # add NES and FDR for each comparison: join tables
    for run in sorted(table_map.keys(), key=col_sort_func):
        res = tf.joinTables(res, table_map[run], joinCols=[(0,0)], joinType='leftouter')
    
    # replace None by pvalue 1 or enrichment score 0 -> facilitate plotting of the data
    for c in range(1, res.colNum()):
        if(c%2 == 0):
            res.modifyColumn(c ,modifying_function=lambda v: '1' if v=='None' else v)
        else:
            res.modifyColumn(c ,modifying_function=lambda v: '0' if v=='None' else v)
            
    # sort sets by number of significant samples and sum of absolute normalized enrichment score
    if(sort_by_score is True):
        res.addColumn(float, 'abs_sum_NES', 0)
        res.modifyColumn('abs_sum_NES', _sum_abs_nes, wholeRow=True)
        res.addColumn(int, 'count_FDR<0.05')
        res.modifyColumn('count_FDR<0.05', _count_sign_fdrs ,wholeRow=True)
        res.sortRows(['count_FDR<0.05', 'abs_sum_NES'], [lambda x:x, lambda x:x], [False, False])
    
    # write to file     
    tf.writeTable(res, outfile, sep='\t', header=True)
    if(writeExcel is True):
        components=outfile.split('.')
        excelout = '.'.join(components[0:len(components)-1])+'.xls'
        tf.writeExcelTable(res, excelout, header=True)
        
def _merge_nes_scores(table, row):
    '''
    auxiliary function to merge 2 enrichment scores from negative and positive enrichment
    '''
    
    # get positive and negative enrichment score for the current row = gene set
    neg_enr = table.get(row, 'NES_neg')
    pos_enr = table.get(row, 'NES_pos')
    
    # no enrichment score available
    if(neg_enr=='None' and pos_enr=='None'):
        return '0'
    
    # positive enrichment score only
    elif(neg_enr=='None' and pos_enr!='None'):
        return pos_enr
    
    # negative enrichment score only
    elif(neg_enr!='None' and pos_enr=='None'):
        return neg_enr
    
    # special case: enrichment for positive and negative fold changes
    else:
        sys.stderr.write('WARNING: Positive and negative enrichment score for '+table.get(row, 'NAME')+'\n')
        return neg_enr+','+pos_enr
       
def _merge_fdr(table, row):
    '''
    auxiliary function to merge 2 fdrs from negative and positive enrichment
    '''
    
    # get positive and negative enrichment score for the current row = gene set
    neg_fdr = table.get(row, 'FDR_neg')
    pos_fdr = table.get(row, 'FDR_pos')
    
    # no pvalue available
    if(neg_fdr=='None' and pos_fdr=='None'):
        return '1'
    
    # pvalue for positive enrichment only
    elif(neg_fdr=='None' and pos_fdr!='None'):
        return pos_fdr
    
    # pvalue for negative enrichment only
    elif(neg_fdr!='None' and pos_fdr=='None'):
        return neg_fdr
    
    # special case: pvalue for positive and negative enrichment
    else:
        sys.stderr.write('WARNING: Positive and negative FDR value for '+table.get(row, 'NAME')+'\n')
        return neg_fdr+','+pos_fdr
    
def _sum_abs_nes(table, row):
    '''
    auxiliary function for evaluating the sum of absolute values of all enrichment scores for a gene set
    '''
    
    # stores current sum of absolute enrichment score values
    cursum=0
    
    # iterate over all columns and consider only enrichment score columns (starting with NES*)
    for c in range(0, table.colNum()):
        cname = table.getColumnName(c)
        if(cname.startswith('NES')):
            
            # add absolute value of current enrichment score
            cursum+=abs(float(table.get(row,c)))
    
    # return final sum
    return cursum

def _count_sign_fdrs(table, row):
    '''
    auxiliary function for counting the samples with FDR<0.05 for a gene set
    '''
    
    # stores current number of samples with significant FDR
    curcount=0
    
    # iterate over all columns and consider only p value columns (starting with FDR*)
    for c in range(0, table.colNum()):
        cname = table.getColumnName(c)
        if(cname.startswith('FDR')):
            
            # increment count if FDR is below significance niveau of 0.05
            if(float(table.get(row,c))<0.05):
                curcount+=1
    
    # return final count
    return curcount
    

def plot_gsea_summary(summary_table, plotout, title, overview=False, nes_range=(-2,2), offset_left=4, offset_bottom=1.5, plot_top=None):
    '''
    visualization of nes and fdr values for a set of GSEA runs
    'summary_table': table produced by 'generate_overview_table' giving gene sets with their nes and fdr in different runs of GSEA
    'plotout': file for saving the plot
    'title': name of the plot
    'overview': flag for visualizing many gene sets, if set to True it does not label the gene sets and does not write pvalues into the cells
    'nes_range': range for colorcoding the normalized enrichment scores -> nes_range[0] = dark red, nes_range[1] = dark blue
    'left_offset': width of left part of the figure for the gene set names (given in inches), default: 4 inches, only used if overview=False
    plot structure: rows = gene sets, columns = GSEA runs, color of a cell = NES of the gene set in the GSEA run
    '''
    
    # read table with plotting data
    tab= tf.readTable(summary_table, sep='\t', header=True)
    if(plot_top is not None):
        tab = tf.selectRows(tab, lambda _,r: r<plot_top)
    
    # find column positions with NES scores
    nescols=[]
    for c in range(0, tab.colNum()):
        if(tab.getColumnName(c).startswith('NES')):
            nescols.append(c)
    
    # extract names GSEA runs = suffixes of NES columns
    names=[]
    for col in nescols:
        name_parts=tab.getColumnName(col).split('_')
        names.append('_'.join(name_parts[1:len(name_parts)]))
        
    # extract names of gene sets
    gene_sets = tab.getColumn(0)
    
    # extract plotting data from the table: normalized enrichment scores and adjusted pvalues
    plotarray = np.zeros((tab.rowNum(),len(nescols)))
    pvalarray = np.zeros((tab.rowNum(),len(nescols)))
    for rowInd in range(0, tab.rowNum()):
        # iterate over NES columns of all GSEA runs
        for arraypos, colInd in enumerate(nescols):
            # plotarray[0] -> row plotted at the bottom, plotarray[rowNum-1] -> row plotted at the top
            # normalized enrichment score
            plotarray[tab.rowNum()-1-rowInd][arraypos]=tab.get(rowInd,colInd)
            # get corresponding pvalue for the current sample
            pvalcol = re.sub('NES_', 'FDR_', tab.getColumnName(colInd))
            pvalarray[tab.rowNum()-1-rowInd][arraypos]=tab.get(rowInd,pvalcol)
    
    # set up figure size
    # space for the columns (fixed width per column)
    map_width=len(nescols)*1
    # space for sample names at the bottom
    #offset_bottom=1.5
    # space for figure title
    offset_top= 0.5
    # space for colorbar at the right
    offset_right = 1.2
    # annotated gene names -> height depends on the number of rows, create additional space at the left for gene set names
    if(overview is False):
        height=tab.rowNum()*0.25
    # do not annotated gene names -> height independent of row number, only small margin at left hand side
    else:
        height=10
        offset_left=0.2
    height = height+offset_bottom+offset_top
    width=map_width+offset_left+offset_right
    # define figure and plotting area
    f=plt.figure(figsize=(width, height))
    # position of plot area (without axis labels and color legend) left, bottom, width, heigth as fraction of total figure size
    f.add_axes([offset_left/width, offset_bottom/height, map_width/width, (height-offset_top-offset_bottom)/height])

    # add column names: samples
    plt.xticks(np.arange(0.5,len(plotarray),1), names, fontsize=10, rotation=90)
    plt.xlabel('Samples', fontsize=12)
    
    # add row names: gene sets (as sorted in summary table) in inverted order -> inverted oder=first element of table plotted as top of the colorplot
    if(overview is False):
        plt.yticks(np.arange(tab.rowNum()-0.5,0,-1), gene_sets, fontsize=10)
    else:
        plt.tick_params(axis='y', left='off', labelleft='off', which='both')
    
    # set up color map: transition blue (down reg) -> white -> red (up reg)
    colors = [(0, (5/255,113/255,176/255)), (0.375, (1, 1, 1)), (0.625, (1, 1, 1)), (1, (202/255,0/255,32/255))]
    cm = LinearSegmentedColormap.from_list('my_list',colors, N=200)
    
    # color plot of NES for gene sets vs. GSEA runs 
    plt.pcolor(plotarray, cmap=cm, vmin=nes_range[0], vmax=nes_range[1], edgecolors='black')
    
    # add pvalue information
    ax = plt.gca()
    # row of plot = y-coordinate, column of plot= x-coordinate
    for x in range(0, len(pvalarray[0])):
        for y in range(0, len(pvalarray)):
            pvaltext = '{:.2f}'.format(float(pvalarray[y][x]))
            if(overview is False):
                ax.text(x+0.5,y+0.5, pvaltext, color='black', horizontalalignment='center', verticalalignment='center', fontsize=10)
            # mark significant cells with a star
            if(float(pvalarray[y][x])<0.05):
                ax.plot(x+0.8, y+0.5, marker='*', color='gold', markersize=8)
    
    # add legend for NES color coding
    barax = f.add_axes([(width-offset_right+0.2)/width, offset_bottom/height, 0.3/width, (height-offset_top-offset_bottom)/height])
    b=plt.colorbar(cax=barax)
    b.set_label('Normalized Enrichment Score', fontsize=12)
    
    # add legend for p values
    if(overview is False):
        lax = f.add_axes([0,0,offset_left/width, offset_bottom/height])
        lax.set_axis_off()
        lax.plot(0.05, 0.1, marker='*', color='gold', markersize=8, transform=lax.transAxes)
        lax.text(0.07, 0.1, 'Significant at Level 0.05', transform=lax.transAxes, horizontalalignment='left', verticalalignment='center', fontsize=10)
        rect = mpatches.Rectangle(xy=(0.05,0.2), width=0.2, height=0.12, linewidth=1, edgecolor='black', facecolor='none', transform=lax.transAxes)
        lax.add_patch(rect)
        lax.text(0.27, 0.2, 'per Sample and Set', horizontalalignment='left', verticalalignment='bottom', fontsize=10, transform=lax.transAxes)
        lax.text(0.15, 0.26, 'Pvalue', transform=lax.transAxes, horizontalalignment='center', verticalalignment='center', fontsize=10)
        lax.text(0.05, 0.4, 'GSEA Adjusted Pvalue', transform=lax.transAxes, horizontalalignment='left', verticalalignment='bottom', fontsize=12)
    
    # set title
    plt.text(s=title, fontsize=14, x=(offset_left+map_width/2)/width, y=(height-offset_top+0.2)/height, transform=f.transFigure, horizontalalignment='center')

    plt.savefig(plotout)
    
    
def gsea_set_heatmap(geneset_file, edgeRcomplete_table, outfolder, genesetname=None, testnames=None, fc_range=(None,None)):
    '''
    'geneset_file': file downloaded from msigdb with the gene sets used for GSEA analysis
    'edgeRcomplete_table': table with log2 fold changes of all edgeR runs
    'outfolder': folder for storing selected fold change tables and heatmap plot
    'genesetname': if None, all gene sets are from 'geneset_file' are analyzed
        otherwise genesetname is the name of the set to analyze or a list of set names
    'testnames': if None, all fold changes columns of 'edgeRcomplete_table' are plotted
        otherwise testnames lists all fold change columns (DE test names) to plot
    'fc_range': tuple with minimal and maximal log2 fold change value for the colorbar, default: (None, None) = limits chosen automatically
    returns
        list of generated heatmap figures, order of figure corresponds to the order given in 'geneset_file'
    '''
    
    # import only used for this function (but not required watchdog module)
    import utils.differential_expression as de
    
    # list of generated figures
    fig_list=[]
    
    # read gene set file
    with open(geneset_file, 'rt') as gsreader:
        for line in gsreader:
            line = line.strip('\n')
            content = line.split('\t')
            cursetname = content[0]
            if(genesetname is None or cursetname==genesetname or cursetname in genesetname):
                targets = [ x.lower().capitalize() for x in content[2:]]
                print('Plotting '+cursetname)
                print('gene set: '+str(len(targets)))
                f = de.cluster_by_foldchange(edgeRcomplete_table, test_names=testnames, gene_names=targets,
                    plot_data_file = os.path.join(outfolder, cursetname+'.txt'), heatmapfile = os.path.join(outfolder, cursetname+'.svg'),
                    plot_range=fc_range, title=cursetname)
                fig_list.append(f)
    
    return fig_list

def gsea_set_scatter(geneset_file, edgeRcomplete_table, outfolder, genesetname=None, testpairs=None, alpha=0.01):
    '''
    'geneset_file': file downloaded from msigdb with the gene sets used for GSEA analysis
    'edgeRcomplete_table': table with log2 fold changes of all edgeR runs
    'outfolder': folder for storing the scatter plots
    'genesetname': if None, all gene sets are from 'geneset_file' are analyzed
        otherwise genesetname is the name of the set to analyze or a list of set names
    'testpairs': if None, all  pairs of fold changes columns of 'edgeRcomplete_table' are plotted
        otherwise testnames lists all fold change pairs (DE test names) to plot
    returns
        list of generated scatter figures, order of figure corresponds to the order given in 'geneset_file'
    '''
    
    # import only used for this function (but not required watchdog module)
    import utils.fc_scatter as fcs
    
    if(isinstance(edgeRcomplete_table, str)):
        edgeRcomplete_table = tf.readTable(edgeRcomplete_table, header=True, sep='\t')
    
    # list of generated figures
    fig_list=[]
    
    # read gene set file
    with open(geneset_file, 'rt') as gsreader:
        for line in gsreader:
            line = line.strip('\n')
            content = line.split('\t')
            cursetname = content[0]
            if(genesetname is None or cursetname==genesetname or cursetname in genesetname):
                targets = set([ x.lower().capitalize() for x in content[2:]])
                print('Plotting '+cursetname)
                print('gene set: '+str(len(targets)))
                
                # reduce fold changes to gene set
                plot_table = tf.selectRows(edgeRcomplete_table, lambda t,r: t.get(r, 'name') in targets)
                                
                # create scatter plot
                plotfile = os.path.join(outfolder, cursetname+'.svg')
                if(testpairs is not None):
                    f = fcs.multiple_fc_scatter_plots(plot_table, testpairs, plotfile,
                        colormode='all', alpha=alpha, adjustP=False, title=cursetname, add_dot_counts=True)
                else:
                    f = fcs.fc_scatter_overview(plot_table, plotfile, samples=None, allpairs=False, colormode='all', alpha=alpha, add_dot_counts=True)
                fig_list.append(f)
    
    return fig_list

def prep_go_enrichment_significant(infile, outfile, alpha=None, resultColumn='ID', fc_direction='both', generateBackground=False):
    '''
    generates one or more input files for gene enrichment analysis given the results of edgeR from DETest module
    'infile': results from edgeR (results for all genes or only significant genes)
    'outfile': prefix for input file for GO enrichment (with one gene per line) created by this function
    'resultColumn': one or more column names -> for each column name, a separate file is generated
    'alpha': significance level, if set to None all genes of 'infile' are used (e.g. input= edgeR.significant.csv)
            if a value is set for alpha, all genes with pvalue < alpha are selected
    'fc_direction': a single value or a list of values that can take 3 values:
            'both' (no filtering based on fold change),
            'up': take only genes with positive fold change
            'down': take only genes with negative fold change
            -> for each fc_direction, a separate file is generated
    '''
    
    # generate Backround -> whole edgeR output as input, requires to select significant genes
    if(generateBackground is True and alpha is None):
        raise ValueError('Generate Background is only possible if a cutoff is set for the pvalue (alpha not None)!')
    
    # read in the resutls from edgeR (either edgeR.significant.csv or edgeR.all.csv)  
    if not isinstance(resultColumn, list):
        resultColumn = [resultColumn] 
    tab=tf.readTable(infile, sep='\t', header=True, colsToRead=['log2FC', 'adj.PValue']+resultColumn)
    
    # generate background for enrichment analysis
    if (generateBackground is True):
        for wcol in resultColumn:
            tf.writeTable(tab, outfile+'_background_'+wcol+'.txt', sep='\t', header=False, colsToWrite=[wcol])
    
    # select significant genes if required
    if(alpha is not None):
        tab = tf.selectRows(tab, lambda t,r: float(t.get(r, 'adj.PValue'))<alpha)
    
    # select genes with the required fold change direction
    if not isinstance(fc_direction, list):
        fc_direction = [fc_direction]
    for fcdir in fc_direction:
        if(fcdir == 'up'):
            res = tf.selectRows(tab, lambda t,r: float(t.get(r, 'log2FC'))>0)
        elif(fcdir=='down'):
            res = tf.selectRows(tab, lambda t,r: float(t.get(r, 'log2FC'))<0)
        elif(fcdir=='both'):
            res = tab
        else:
            raise ValueError('Unknown fc_direction mode: '+str(fcdir))
        
        # iterate over list of desired gene identifiers and write an output file for each
        for wcol in resultColumn:
            tf.writeTable(res, outfile+'_'+fcdir+'_'+wcol+'.txt', sep='\t', header=False, colsToWrite=[wcol])
    
    