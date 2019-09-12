'''
Created on 4 Dec 2017

@author: friedl
'''

import argparse
import os
import shutil
import errno
import time
import resource # @UnresolvedImport
import re
import datetime
import subprocess
import sys
import signal

# type functions for command line parsing

# check file ending for input file only if different formats are allowed, handled differently and recognized by the file ending (either by the watchdog module or external executables)

def valid_file_path(string):
    ''' checks if a input file exists'''
    
    if(not os.path.exists(string) or not os.path.isfile(string)):
        raise argparse.ArgumentTypeError('File does not exist: \''+string+'\'')
    return string

def valid_folder_path(string):
    ''' checks if a input folder exists'''
    
    if(not os.path.exists(string) or not os.path.isdir(string)):
        raise argparse.ArgumentTypeError('Folder does not exist: \''+string+'\'')
    return string

def valid_fastq_path(string, endings=['fq', 'fastq', 'fq.gz']):
    ''' checks if string is a valid path to a fastq file (based on file ending) '''
    
    valid_file_path(string)
    for ending in endings:
        if(string.endswith('.'+ending)):
            return string
    raise argparse.ArgumentTypeError('File is not in FASTQ format: \''+string+'\'')

def valid_uncompressed_fastq_path(string):
    ''' checks if string is a valid path to an uncompressed fastq file (based on file ending) '''
    
    valid_file_path(string)
    for ending in ['fq', 'fastq']:
        if(string.endswith('.'+ending)):
            return string
    
    # only reached if invalid file ending
    if(string.endswith('.gz')):
        raise argparse.ArgumentTypeError('Module does not support compressed FASTQ format: \''+string+'\'')
    raise argparse.ArgumentTypeError('File is not in FASTQ format: \''+string+'\'')

def valid_igv_genome_file(string):
    ''' igvtools recognizes file formats from file name endings -> genome files end with .chrom.sizes or .genome'''
    
    valid_file_path(string)
    valid_ending(string, ['.chrom.sizes', '.genome'])
    return string    

def valid_mapping(string):
    ''' checks if string is a valid path to a sam or bam file '''
    
    valid_file_path(string)
    if(not re.search('\.(sam|bam)$', string, re.IGNORECASE)):
        raise argparse.ArgumentTypeError('File is not in SAM or BAM format: \''+string+'\'')
    return string

def valid_bam(string):
    ''' checks if string is a valid path to a bam file (based on file ending) '''
    
    valid_file_path(string)
    if not string.endswith('.bam'):
        raise argparse.ArgumentTypeError('File is not in BAM format: \''+string+'\'')
    return string

def valid_indexed_bam(string):
    ''' checks if string is a valid path to a bam file and if an index exists for the file '''
    
    valid_file_path(string)
    if(not re.search('\.bam$', string, re.IGNORECASE)):
        raise argparse.ArgumentTypeError('File is not in BAM format: \''+string+'\'')
    if not (os.path.exists(string+'.bai') or os.path.exists(re.sub('bam$', 'bai', string))):
        raise argparse.ArgumentTypeError('BAM file is not indexed: \''+string+'\'')
    return string

def valid_bwa_index(string):
    ''' checks the files of a bwa index'''
    
    valid_file_path(string+'.amb')
    valid_file_path(string+'.ann')
    valid_file_path(string+'.bwt')
    valid_file_path(string+'.pac')
    valid_file_path(string+'.sa')
    return string

def valid_star_index(string):
    ''' checks the files of a STAR index'''
    
    for file_name in ['Genome', 'SA', 'SAindex', 'chrLength.txt', 'chrName.txt', 'chrNameLength.txt', 'chrStart.txt']:
        valid_file_path(os.path.join(string,file_name))
    return string


def valid_star_output(string):
    ''' check if all files required for circRNA finder exist'''
    
    for ending in ['Chimeric.out.junction', 'Chimeric.out.sam', 'SJ.out.tab']:
        valid_file_path(string+ending)
    return string

def valid_exec(string):
    ''' checks if the executable is in the PATH variable or if the path to the executable exists '''
    
    if(shutil.which(string)==None):
        valid_file_path(string)
    return string

def positive_integer(string):
    ''' checks if a input string represents a valid integer >=1 '''
    
    nr = int(string)
    if(nr<1):
        raise argparse.ArgumentTypeError('invalid positive integer value: \''+string+'\' -> integer >=1 is required!')
    return nr

def positive_integer_or_zero(string):
    ''' checks if a input string represents a valid integer >=0 '''
    
    nr = int(string)
    if(nr<0):
        raise argparse.ArgumentTypeError('invalid positive integer value: \''+string+'\' -> integer >=0 is required!')
    return nr

def valid_regex_with_one_group(string):
    
    #check regular expression by compiling it
    try:
        re.compile(string)
    except re.error as e:
        raise argparse.ArgumentTypeError('invalid python regular expression: \''+string+'\'-> '+str(e)) from None

    group_nr = re.compile(string).groups
    if(group_nr<1):
        raise argparse.ArgumentTypeError('regular expression without group: \''+string+'\'')
    
    return string

def valid_string_boolean(string):
    ''' checks if the string equals 'yes' or 'no' and translates it into True or False '''
    
    # check valid input
    if not string in ['yes', 'no']:
        raise argparse.ArgumentTypeError('only "yes" or "no" allow but not '+'"'+string+'"')
    
    # return boolean
    if string == 'yes':
        return True
    else:
        return False
    

def valid_list_of_files(string, separator=',', min_length=2):
    ''' checks if all files of the list are valid paths '''
    
    # get components of the list
    file_list = string.split(separator)
    
    # check number of files
    if(len(file_list)<min_length):
        raise argparse.ArgumentTypeError('file list with '+str(len(file_list))+' < '+str(min_length)+' elements')
    
    # check paths
    for file_path in file_list:
        valid_file_path(file_path)
    
    # return list of strings
    return file_list

# checks endings of output file paths

def valid_outfile_ending_sam(string):
    ''' checks if output file has sam ending'''
    
    valid_ending(string, ['sam'])
    return string

def valid_outfile_ending_sam_bam(string):
    ''' checks if output file ends with sam or bam '''
    
    valid_ending(string, ['sam', 'bam'])
    return string
    
def valid_outfile_ending_fastq_or_fasta(string):
    ''' checks if output file has ending for fastq (compressed and uncompressed) or fasta fq.gz, .fq, .fastq or. fa'''
   
    valid_ending(string, ['fq', 'fastq', 'fq.gz', 'fa'])
    return string
    
def valid_ending(string, endings):
    ''' checks if string has a valid ending from endings'''
    
    for ending in endings:
        if(string.endswith(ending)):
            return string
    raise argparse.ArgumentTypeError('Incorrect file ending for '+string+', expecting ending in '+str(endings))
    

# create directories for output files

def make_parent_dirs(folder_path):
    ''' folder_path: file name or file name prefix, the method creates all parent directories '''

    try:
        os.makedirs(os.path.dirname(folder_path))
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def create_folder(folder_path):
    ''' creates folder given by 'folder_path' and all non-existent parent directories, in contrast to make_parent_dirs 'folder_path' does not contain a filename '''
    
    try:
        os.makedirs(folder_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

# measure resources consumed by modules
        
def get_current_time():
    ''' returns the current time as formatted string and number'''
    
    time_point = time.time()
    time_string = time.strftime("[%d %b %Y %H:%M:%S]: ", time.localtime())
    return time_point, time_string

def print_resources(elapsed_time, child_processes):
    ''' writes consumed resources to stdout '''
    
    print('\nResources:')
    print('wall_clock_time='+str(elapsed_time)+'='+str(datetime.timedelta(seconds=elapsed_time)))
    
    rusage_vars = [(resource.RUSAGE_SELF, 'self')]
    if(child_processes is True):
        rusage_vars.append((resource.RUSAGE_CHILDREN, 'children'))
    for rusage_var, name in rusage_vars:
        r_data = resource.getrusage(rusage_var)
        r_data = re.sub('resource.struct_rusage\(|\)', '', str(r_data))
        for res_info in r_data.split(', '):
            print(name+'_'+res_info) 

            
# write watchdog tables

def write_watchdog_return_file(filename, variable_value_pairs):
    ''' writes file with the return variables as given tuple (varname, value) in 'variable_value_pairs' for watchdog'''
    
    with open(filename, 'wt') as watchdog_writer:
        for var_name, var_val in variable_value_pairs:
            watchdog_writer.write(var_name+'\t'+str(var_val)+'\n')
        watchdog_writer.write('?EOF!')

def write_watchdog_process_table(process_tab_file, table_rows):
    ''' writes a watchdog process table: 2d list with each element = rows, each row is presented as a list of strings'''
    
    with open(process_tab_file, 'wt') as pt_writer:
        # write rows
        for row in table_rows:
            pt_writer.write('\t'.join(row)+'\n')
            
# execute workflows

def execute_watchdog_workflow(workflow_file, watchdog_exec, mailfile, error_message, port=None, timeout=60):

    # execute watchdog with small test workflow
    os.setpgrp()
    command = [watchdog_exec, '-x', workflow_file, '-mailConfig', mailfile]
    if port is not None:
        command+=['-p', str(port)]
    print('Running: '+' '.join(command))
    p=subprocess.Popen(command)
    try:
        p.wait(timeout)
    except:
        sys.stderr.write(error_message)
        os.killpg(0, signal.SIGKILL)
        sys.exit(1)

# watchdog resume files

def remove_watchdog_resume_files(workflow_dir, workflow_name=None):
    ''' remvoes all watchdog resume files in a given directory for the workflow named 'workflow_name' (if set to None, all files are removed) '''

    for file in os.listdir(workflow_dir):
        if file.endswith('watchdog.status.log'):
            if workflow_name is None or file.startswith(workflow_name):
                os.remove(os.path.join(workflow_dir, file))
