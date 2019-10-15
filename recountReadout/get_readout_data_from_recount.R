library(recount)
library(recount.bwtool)

# readin and readout regions as granges object and optionally create a bed file with the GRanges
read_regions <- function(region_file, out_bed = NULL){
    
    # read table 
    genetab <- read.table(region_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # first return value -> gene ids
    geneids <- genetab$geneid
    
    # granges without strand for downstream regions
    downstream <- GRanges( seqnames=genetab$chr, ranges=IRanges(start=genetab$downstream_start, end=genetab$downstream_end),
                           gene_id=genetab$geneid, direction="d")
    # granges without strand for upstream regions
    upstream <- GRanges( seqnames=genetab$chr, ranges=IRanges(start=genetab$upstream_start, end=genetab$upstream_end),
                         gene_id=genetab$geneid, direction="u")
    # second return values -> regions to quantify
    all_regions <- c(downstream, upstream)
    
    # write GRanges object to a bed file
    if (!is.null(out_bed)) {
        rtracklayer::export(all_regions, con=out_bed, format = "BED")
    }
    
    # return geneids and granges in a named list
    return(list("geneids"=geneids, "regions"=all_regions, "bed_file"=out_bed))
}


# downloads rse_gene file for project 'projectid' if not available in folder 'folder/projectid'

download_gene_rse<-function(projectid, folder){
    
    # check if the rse_gene file was already downloaded
    gene_path <- file.path(folder, "rse_gene.Rdata")
    if(!file.exists(gene_path)){
        download_study(projectid, type="rse-gene", outdir=folder)
    }
    
    # return the file path
    return(gene_path)
}


# check if the sample was already analysed successfully (-> restart crashed task)

output_file_exists <- function(outfolder, sampleid, expected_size){
    
    outfile <- file.path(outfolder, paste(sampleid, ".tab", sep=""))
    
    # check if the file exsits
    if(!file.exists(outfile)){
        return(FALSE)
    }
    
    # try to read the table, reading fails if the task crashed during writing the file (incomplete line)
    tryCatch({
        res <- read.table(outfile, stringsAsFactors = FALSE, header=TRUE, sep="\t")
        # check if all rows (genes) were written
        if(nrow(res)==expected_size){
            return(TRUE)
        }
        else{
            return(FALSE)
        }
    },
    error = function(err){
        print("error")
        return(FALSE)
    })
    
}

# downloads bigWig file for sample 'sampleid' of project 'projectid' if not available in 'folder/projectid/'

download_bigwig <- function(projectid, sampleid, folder){
    
    # check if the bigwig file was already downloaded
    bw_path <- file.path(folder, paste(sampleid,".bw", sep=""))
    if(!file.exists(bw_path)){
        url <- recount_url[ recount_url$project==projectid & startsWith(basename(recount_url$url), sampleid) & endsWith(recount_url$url, "bw"), "url"]
        downloader::download(url, destfile=bw_path, mode="wb")
    }
    
    return(bw_path)
}


# calculates coverage

sample_coverage_bwtool <- function(sampleid, bed_path, bigwig_path, bw_sum_folder, bwtool_exec = 'bwtool'){
    
    # call bwtool sum
    script_path <- system.file("extdata", "jhpce", "sum.sh", package="recount.bwtool")
    bwsum_output <- file.path(bw_sum_folder, paste(sampleid, ".sum.tsv", sep=""))
    cmd <- paste("bash", script_path, bwtool_exec, bed_path, bigwig_path, bwsum_output)
    system(cmd)
    
    # read extracted counts
    res <- read.table(bwsum_output, header = FALSE, colClasses = list(NULL, NULL, NULL, "numeric"))
    colnames(res) <- sampleid
    return (as.matrix(res))
}

# write results to file 

write_coverage_results <- function(result, outfolder, sample_id, gene_data, gene_counts, sample_info){
    
    # scale the counts obtained for upstream and downstream regions
    rownames(result) <- paste(gene_data[["regions"]]$gene_id, gene_data[["regions"]]$direction, sep='_' )
    up_down_rse <- SummarizedExperiment(assays = list("counts" = result), colData = sample_info, rowRanges = gene_data[["regions"]])
    up_down<- assays(scale_counts(up_down_rse))$counts

    # create a table with scaled counts for the gene, upstream and downstream region
    res_tab <- cbind(rownames(gene_counts), gene_counts[, sample_id], up_down[paste(rownames(gene_counts), "u", sep="_"), sample_id], up_down[paste(rownames(gene_counts), "d", sep="_"), sample_id])
    colnames(res_tab) <- c("gene_id", "gene_scaled", "upstream_scaled", "downstream_scaled")
    # print(file.path(outfolder, paste(sample_id, ".tab", sep="")))
    write.table(res_tab, file.path(outfolder, paste(sample_id, ".tab", sep="")), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}


# function used to process several samples in parallel

process_sample<- function(sample_id, sample_gene_counts, sample_info, project_id, tmpfolder, outfolder, gene_data, removebw = TRUE){
    
    if( !output_file_exists(outfolder, sample_id, length(gene_data[["geneids"]])) ){
        
        message(paste(Sys.time(), "processing sample", sample_id))
        
        # download the bigwig file
        bw_file <- download_bigwig(project_id, sample_id, tmpfolder)
        # call bwtool summary
        res <- sample_coverage_bwtool(sample_id, gene_data[["bed_file"]], bw_file, tmpfolder)
        # write results table
        write_coverage_results(res, outfolder, sample_id, gene_data, sample_gene_counts, sample_info)
        
        # remove all temporary files: bigwig and tsv file
        if(removebw){
            unlink(bw_file)
        }
    }
    else{
        message(paste(Sys.time(), "skipping sample", sample_id))
    }
}


is_recount_project <- function(project_id){
    if (project_id %in% recount_url$project) {
        return(TRUE)
    }
    else {
        stop(paste("Unknown projectID", project_id, "-> Empty result"))
    }
}

# yields same results as method above but with parallel download of bigwig files
process_project<- function(project_id, gene_file, outfolder, tmpfolder, threads=1, remove_tmp_data = TRUE, download_parallel = FALSE){
    
    # check if project was already analyzed 
    if(file.exists(file.path(outfolder, paste(project_id, ".finished", sep="")))){
        return()
    }
    
    # check if there is a project corresponding to the given projectID
    if(!is_recount_project(project_id)){
        return()
    }

    # read gene data and create a bed file with readin and readout data
    bed_file <- file.path(tmpfolder, "readout_regions.bed")
    gene_data <- read_regions(gene_file, out_bed = bed_file)
    
    # load project data at gene level
    rse_file <- download_gene_rse(project_id, tmpfolder)
    load(rse_file)
    
    # save scaled counts for selected genes
    quant_gene <- assays(scale_counts(rse_gene))$counts
    #as.data.frame(quant_gene)
    quant_gene <- quant_gene[ gene_data[["geneids"]], , drop=FALSE]
    
    # save sample data
    sample_data <-colData(rse_gene)
    
    # clean up main memory
    remove(rse_gene)
    gc()
    
    # option to process samples in parallel
    if (threads==1) {
        bpparam = BiocParallel::SerialParam()
    }
    else {
        bpparam = BiocParallel::MulticoreParam(threads, timeout=3600)
    }
    
    # first download everything
    if(!download_parallel){
        for (sample_id in rownames(sample_data)){
            download_bigwig(project_id, sample_id, tmpfolder)
        }
    }
    
    res <- bpmapply(process_sample, rownames(sample_data), lapply(rownames(sample_data), function(x){quant_gene[,x,drop=FALSE]}), 
                   lapply(rownames(sample_data), function(x){sample_data[x,,drop=FALSE]}),
                   MoreArgs = list("project_id" = project_id, "tmpfolder" = tmpfolder, "outfolder"=outfolder, "gene_data"=gene_data, "removebw"=remove_tmp_data),
                   SIMPLIFY = FALSE, BPPARAM=bpparam)
    
    
    fileConnection<-file(file.path(outfolder, paste(project_id, ".finished", sep="")))
    writeLines(rownames(sample_data), fileConnection)
    close(fileConnection)
}


# # for debugging
# main_test <-function(){
# 
#     project_id <- "SRP056912" #"SRP056912" #SRP056378" #"SRP058633"
#     gene_datatab <- "/mnt/raidinput/tmp/friedl/recount/gene_data/selected_regions10000_5000.tsv"
#     outfolder <- "/mnt/raidinput/tmp/friedl/recount/out_test/sample2"
#     tmpfolder <- "/mnt/raidinput/tmp/friedl/recount/tmp_test/sample2"
#     threads <- 4
#     removeTmp <- FALSE #as.logical()
#     multi_thread_download <- FALSE
#     res <- process_project(project_id = project_id, gene_file = gene_datatab, outfolder = outfolder, tmpfolder = tmpfolder, threads = threads, remove_tmp_data = removeTmp, download_parallel = multi_thread_download)
# }
# 
# main_test()

main <-function(){

    args <- commandArgs(trailingOnly = TRUE)
    project <- args[1]
    gene_region_tsv <- args[2]
    res_folder <- args[3]
    tmp_folder <- args[4]
    threads <- as.numeric(args[5])
    remove_tmp <- as.logical(args[6])
    multi_thread_download <- as.logical(args[7])
    process_project(project_id = project, gene_file = gene_region_tsv, outfolder = res_folder, tmpfolder = tmp_folder, threads = threads, remove_tmp_data = remove_tmp, download_parallel = multi_thread_download)
}

main()
