'''
Created on Mar 27, 2018

@author: friedl
'''

import matplotlib.pyplot as plt
import seaborn as sns


def barplot_read_numbers(sample_counts, plot_file):
    '''
    plots raw, trimmed and mapped reads as stacked bar plot
    sample_counts: table produced by algo_mappingsum (columns = count types, rows = samples)
    plot_file: file to save the barplot
    '''
    
    # data to plot
    samples = sample_counts.getColumn('sample')
    total = sample_counts.getColumn('raw')
    trim = sample_counts.getColumn('trimmed')
    mapped = sample_counts.getColumn('mapped')
    
    # style the plot
    sns.set()
    sns.set_context('notebook')
    sns.set_style('whitegrid')
    sns.set_palette('colorblind', desat=0.75)
    plt.figure()
    ax = plt.gca()
    ax.xaxis.grid(False)
    
    # do the actual plotting
    plt.bar(range(1,len(samples)+1), total, 0.8, align='center', label='total')
    plt.bar(range(1,len(samples)+1), trim, 0.8, align='center', label='trimmed')
    plt.bar(range(1,len(samples)+1), mapped, 0.8, align='center', label='mapped')
    
    # add labels for the plot
    plt.xticks(range(1,len(samples)+1), samples, rotation='90')
    plt.title('Reads remaining after QC and mapping step')
    plt.xlabel('Samples')
    plt.ylabel('Number of Reads')
    leg = plt.legend(fancybox=True, frameon=True)
    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_alpha(0.5)
    plt.tight_layout()
    
    # save the plot
    plt.savefig(plot_file)
    plt.close()



def barplot_chrom_groups(sample_counts, plot_file):
    '''
    plots % of mapped reads for a defined groups of chromosomes (e.g. organisms) as stacked bar plot
    sample_counts: table produced by algo_mappingsum (columns = count types, rows = samples)
    plot_file: file to save the barplot
    '''
    
    # data to plot : all columns left of 'mapped'
    mapped_pos = -1
    for c in range(0, sample_counts.colNum()):
        if(sample_counts.getColumnName(c)=='mapped'):
            mapped_pos = c
            break
    
    # data to plot
    samples = sample_counts.getColumn('sample')
    mapped = sample_counts.getColumn('mapped')
    bottom = [0 for _ in range(0, len(samples))]

    # style the plot
    sns.set()
    sns.set_context('notebook')
    sns.set_style('whitegrid')
    sns.set_palette('colorblind', desat=0.75)
    plt.figure()
    ax = plt.gca()
    ax.xaxis.grid(False)
    
    # do the actual plotting
    for group in range(mapped_pos+1, sample_counts.colNum()):
        # calculate percentage of mapped reads for current column
        values = [float(x)/float(y)*100 for x,y in zip(sample_counts.getColumn(group), mapped)]
        plt.bar(range(1,len(samples)+1), values, 0.8, bottom=bottom, align='center', label=sample_counts.getColumnName(group))
        # update position of the bar top
        bottom = [x+y for x,y in zip(values, bottom)]
        
    # add labels for the plot
    plt.xticks(range(1,len(samples)+1), samples, rotation='90')
    plt.title('Distribution of mapped reads')
    plt.xlabel('Samples')
    plt.ylabel('Percentage of mapped reads')
    leg = plt.legend(fancybox=True, frameon=True)
    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_alpha(0.5)
    plt.tight_layout()
    
    # save the plot
    plt.savefig(plot_file)
    plt.close()
    
    