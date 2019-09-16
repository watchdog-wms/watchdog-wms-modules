'''
Created on Dec 17, 2018

@author: friedl
'''

import sys, os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/../sharedUtils/python')
import os
import shutil
import argparse
import subprocess
import watchdog_utils as wutils

def get_command_line_options():
    '''command line parser for the watchdog module for analyzing readout in recount data '''
    
    parser = argparse.ArgumentParser(description='Calculates readout for every sample in a project from recount')
    
    # input
    project_opts = parser.add_mutually_exclusive_group(required=True)
    project_opts.add_argument('--projectID', required=False, default=None, type=str, metavar='SRP_ID',
                        help='project id of a sra project indexed in recount, it is possible to pass several project ids separated by ","')
    project_opts.add_argument('--projectFile', required=False, default=None, type=wutils.valid_file_path, metavar='project.file',
                        help='file with one line giving project ids (file content = all allowed values for projectID)')
    parser.add_argument('--geneTSV', required=True, default=None, type=wutils.valid_file_path, metavar='genes.tsv',
                        help='tab-separated file with genes, cooridantes, exonic basepairs and upstream and downstream regions')
    
    # output
    parser.add_argument('--outfolder', required=True, default=None, metavar='folder',
                        help='folder for saving final results, creates a subfolder for the project with a table of coverage values for every sample in the project')
    parser.add_argument('--tmpfolder', required=True, default=None, metavar='tmp',
                        help='folder for saving temporary data, creates a subfolder for the project (named projectID)')
    
    # options regulating memory on disk and speed
    parser.add_argument('--threads', required=False, default=1, type=wutils.positive_integer, metavar='threadNr',
                        help='number of threads to use, equivalent to number of samples processed in parallel')
    # controls if files *.bw (downloaded bigwig files) and *.sum.tsv are deleted (output of bwtool)
    boolparam1 = parser.add_mutually_exclusive_group()
    boolparam1.add_argument('--removeTmpSampleData', required=False, action='store_true', default=True,
                           help='if this flag is set, temporary files are deleted at the end (default behaviour)')
    boolparam1.add_argument('--noremoveTmpSampleData', required=False, action='store_true', default=False,
                           help='if this flag is set, temporary files are not deleted at the end (keep for reuse)')
    # controls if rse_gene.Rdata is deleted 
    boolparam2 = parser.add_mutually_exclusive_group()
    boolparam2.add_argument('--removeTmpProjectData', required=False, action='store_true', default=True,
                           help='if this flag is set, temporary files are deleted at the end (default behaviour)')
    boolparam2.add_argument('--noremoveTmpProjectData', required=False, action='store_true', default=False,
                           help='if this flag is set, temporary files are not deleted at the end (keep for reuse)')
    # run download from recount server in parallel -> download fails if too many queries in short time
    boolparam3 = parser.add_mutually_exclusive_group()
    boolparam3.add_argument('--downloadParallel', required=False, action='store_true', default=False,
                            help='downloads bigwig data in parallel, only recommended for big projects and requires a small number of threads (default not set)')
    boolparam3.add_argument('--nodownloadParallel', required=False, action='store_false', default=False,
                            help='first downloads all bigwig files with one thread and then processes all samples, recommended for long project lists with small projects (default set)')
    
    # pass path of R to use
    parser.add_argument('--Rscript', required=False, default='Rscript', type=wutils.valid_exec, metavar='Rscript',
                        help='path to Rscript executable (preferentially version 5.3)')
    
    cmdl_options=parser.parse_args()
    return parser, cmdl_options


def create_outfolders(options, projID):
    ''' creates outfolder and tmpfolder including a subdirectory for the current project if they do not exist, returns the path to the directory with all final files ''' 
    
    project_outfolder = os.path.join(options.outfolder, projID)
    wutils.create_folder(project_outfolder)
    project_tmpfolder = os.path.join(options.tmpfolder, projID)
    wutils.create_folder(project_tmpfolder)
    return project_outfolder, project_tmpfolder

def resolve_booleans(options):
    ''' transform removeTmpData mutually exclusive group into boolean '''
    
    
    # default behaviour -> set remove_tmp to true
    tmpSampleBool = True
    # --noremoveTmpData was added -> set remove_tmp to false
    if options.noremoveTmpSampleData:
        tmpSampleBool = False
        
    tmpProjBool = True
    if options.noremoveTmpProjectData:
        tmpProjBool = False
    
    parallelDownBool = False
    # --downloadPatallel was added-> set to true
    if options.downloadParallel:
        parallelDownBool = True
    
    return tmpSampleBool, tmpProjBool, parallelDownBool

def check_gene_tsv(options):
    ''' checks if all information for calculating the readout is given '''
    
    with open(options.geneTSV, 'rt') as r:
        header = r.readline().strip('\n').split('\t')
        
    for colname in ['chr', 'geneid', 'exonic_bps', 'upstream_start', 'upstream_end',  'downstream_start', 'downstream_end']:
        if not colname in header:
            raise ValueError('Missing column "'+colname+'" in file '+options.geneTSV)
        
def get_project_list(options):
    ''' resolves the mutually exclusive group for the project ids and returns a list of project ids '''
    
    if options.projectID:
        return options.projectID.split(',')
    else:
        with open(options.projectFile, 'rt') as pr:
            content = pr.readlines()
        if(len(content)>1):
            raise ValueError('Project file may only contain one line of project ids separated by ,')
        else:
            return content[0].strip('\n').split(',')
    
def remove_tmp_data(remove_sample_data, remove_proj_data, res_proj_dir, tmp_proj_dir):
    ''' removes temporary project folder or all temporary files if outdir and tmpdir are the same '''
    
    # remove some files
    if remove_sample_data or remove_proj_data:
        
        # tmp files in result folder or do not delete both data types
        if res_proj_dir == tmp_proj_dir or not(remove_sample_data and remove_proj_data):
            for filename in os.listdir(tmp_proj_dir):
                if remove_sample_data and (filename.endswith('.bw') or filename.endswith('.sum.tsv')):
                    os.remove(os.path.join(tmp_proj_dir, filename))
                if remove_proj_data and (filename.endswith('.Rdata') or filename.endswith('.bed')):
                    os.remove(os.path.join(tmp_proj_dir, filename))
        
        # separate tmp dir and remove everything from it
        else:
            shutil.rmtree(tmp_proj_dir)
            

def run_Rscript_recount(project_id, gene_file, out_folder, tmp_folder, threads, remove_tmp, parallel_download, rscript_exec):
    ''' call the Rscript for calculating coverages of the gene, the upstream and the downstream region '''
    
    location = os.path.dirname(os.path.realpath(__file__))
    command = [rscript_exec, os.path.join(location, 'get_readout_data_from_recount.R'), project_id, gene_file, out_folder, tmp_folder, str(threads), str(remove_tmp), str(parallel_download)]
    print('Running command:\n'+' '.join(command))
    subprocess.check_call(command)
    

def main():
    ''' main method for running the watchdog module for calculation of readout in recount data '''
    
    # print command
    print('Program call:')
    print(' '.join(sys.argv)+'\n')
    
    # check options and prepare output location
    _,o = get_command_line_options()
    remove_tmp_sample, remove_tmp_proj, parallel_down = resolve_booleans(o)
    check_gene_tsv(o)
    proj_list = get_project_list(o)
    
    # log
    start_timepoint=wutils.get_current_time()
    print(start_timepoint[1]+'Welcome to the module for calculating readout from recount !\n')
    
    # iterate over all projects
    for projID in proj_list:
        
        proj_timepoint = wutils.get_current_time()
        print(proj_timepoint[1]+'Processing project '+projID)
        
        # create output and tmp folder for the current project
        res_folder, tmp_folder = create_outfolders(o, projID)
    
        # run the R code if possible
        try:
            run_Rscript_recount(projID, o.geneTSV, res_folder, tmp_folder, o.threads, remove_tmp_sample, parallel_down, o.Rscript)
        except subprocess.CalledProcessError:
            sys.stderr.write('An error occurred in the R script \n')
            raise
        
        # clean up temporary data
        finally:
            remove_tmp_data(remove_tmp_sample, remove_tmp_proj, res_folder, tmp_folder)

    # log
    end_timepoint=wutils.get_current_time()
    print(end_timepoint[1]+'Module finished successfully!\n-> the readout data is located at \''+o.outfolder+'\'')
    wutils.print_resources(end_timepoint[0]-start_timepoint[0], child_processes=True)


if __name__ == '__main__':
#     sys.argv=['main_recount_readout.py', '--projectID', 'SRP056378', '--geneTSV', '/mnt/raidinput/tmp/friedl/recount/gene_data/selected_regions10000_5000.tsv', 
#               '--outfolder', '/mnt/raidinput/tmp/friedl/recount/out_test', '--tmpfolder', '/mnt/raidinput/tmp/friedl/recount/tmp_test',
#               '--Rscript', '/home/proj/software/R/Unix/R-3.5.1/bin/Rscript', '--noremoveTmpData']
    main()
    
    