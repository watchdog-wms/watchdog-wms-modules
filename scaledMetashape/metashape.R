# load lib
library(RColorBrewer)
library(gplots)
library(getopt)
library(plotrix)
library(Polychrome)

warnings()

opt <- NULL
# options to parse
spec <- matrix(c('coverageFiles',   'c', 1, "character", # dir of coverage files
                 'genelist',        'g', 2, "character", # optional genelist with clusters of genes
		 'experiment',      'x', 2, "character", # optional experiment type (PRO/ATAC)
                 'bedgraphTable',   't', 1, "character", # bedgraph table
                 'bedname',         'z', 1, "character", # bedfile name
                 'factor',          'f', 2, "character", # optional factor
                 'normLibSize',     'l', 0, "logical",   # norm or not
                 'normBinLength',   'n', 0, "logical",   # norm or not
                 'normShapeSum',    's', 0, "logical",   # norm or not
                 'aggregateFUN',    'a', 1, "character", # aggregate function mean/median
                 'metaFrame',       'm', 1, "numeric", # frame of metageneplot, +-TSS
                 'bins',            'b', 1, "numeric",   # number of bins
                 'plotname',        'p', 1, "character",   # plotname prefix
                 'config',          'i', 2, "character",   # config table for conditions and color
                 'confirmRun2EndFile', 'e', 1, "character"
), ncol=4, byrow=T)

# parse the parameters
opt <- getopt(spec)


factor_array <- c()
if (! is.null(opt$factor)) {
  factor_array <- strsplit(opt$factor, ",")[[1]]
}
genelist_array <- c()
if (! is.null(opt$genelist)) {
  genelist_array <- strsplit(opt$genelist, ",")[[1]]
}

if (!is.null(opt$experiment) && opt$experiment == "PRO") {
	ylabel <- "PRO-seq profile"
} else {
	ylabel <- "occupancy"
}

###################### PARAMETER CHECK ######################
if(is.null(opt$bedgraphTable)) {
  print("[ERROR] Path to bedgraph table not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$bedname)) {
  print("[ERROR] bedname not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$coverageFiles)) {
  print("[ERROR] Path to coverage files not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$bins)) {
  print("[ERROR] Number of bins not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$metaFrame)) {
  print("[ERROR] Length of metageneplot frame not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$aggregateFUN)) {
  print("[ERROR] Aggregation function not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$normLibSize)) {
  print("[ERROR] Normalizing lib size not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$normBinLength)) {
  print("[ERROR] Normalizing bin length not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$normShapeSum)) {
  print("[ERROR] Normalizing shape sum not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}



args <- commandArgs(trailingOnly = FALSE)
functionsFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../../modules/sharedUtils/R/functions.R', sep = '')
source(functionsFileName)
binFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../../modules/binGenome/R/binGenome.lib.R', sep = '')
source(binFileName)


config <- c()
if (! is.null(opt$config)) {
  config <- read.table(opt$config, header=TRUE, sep="\t", quote="", row.names=1, as.is=TRUE)
}


######## labels for x-axis
fixedLabelsStart <- paste0(seq(-opt$metaFrame, 1500, by=750), "bp")
fixedLabelsEnd <- paste(seq(-1500, 3000, by=750), "bp", sep="")
fixedLabelsStart[fixedLabelsStart == "0bp"] <- "TSS"
fixedLabelsEnd[fixedLabelsEnd == "0bp"] <- "TTS"
bodyLabels <- c("", paste(seq(20, 80, by=20), "%", sep=""), "")
fixedLabelsStartTotalBins <- 50 #50 bins
fixedLabelsEndTotalBins <- 50 #50, die fixed, nicht die bins aus dem parameter opt$bins
#linesOnTicks <- c(4,6,13,15)
linesOnTicks <- c(5,14)



transGeneMapping <- read.csv("/home/proj/projekte/sequencing/Illumina/PolIIPausing/TSS-profiling/workflow/genelists/gene2transcript.proseq_uninf.txt", sep="\t", head=T, as.is=T)
maxTranscriptIDs <- read.csv("/home/proj/projekte/sequencing/Illumina/PolIIPausing/TSS-profiling/workflow/genelists/transcripts.txt", sep="\t", head=T, as.is=T)
geneList <- list(transGeneMapping$geneID)
names(geneList) <- c("")


######### init
shapeStoreGenes <- list()
shapeStoreGeneIds <- list()
counterGenes <- 1
shapeOrderStoreIds <- list()
groupedDataStore <- list()
groupedStore <- list()
groupedCondStore <- list()
gene <- list()
names <- c()
shapeStoreReps <- list()
groupedStoreReps <- list()
num_genes <- 0
dfGeneListCollection <- list()


######### functions
lighten <- function(color, add=70) {
    col <- col2rgb(color)
    col <- col + add
    col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) x-add-add else x))), maxColorValue=255)
    return(col)
} 

color_distance <- function(cola, colb) {
	cola_rgb <- col2rgb(cola)
	colb_rgb <- col2rgb(colb)
	diff_red <- abs(cola_rgb[1]-colb_rgb[1]) / 255
	diff_green <- abs(cola_rgb[2]-colb_rgb[2]) / 255
	diff_blue <- abs(cola_rgb[3]-colb_rgb[3]) / 255
	col_diff_percent <- (diff_red + diff_green + diff_blue) / 3 * 100
	return(col_diff_percent)
}

getColPal <- function(col) {
  if (is.na(col) || is.null(col)) {
    #return(brewer.pal(5, "YlGnBu"))
    return(palette(rainbow(5)))
  }
  else if (col=="darkviolet") {
    return(brewer.pal(5, "Purples"))
  }
  else if (col=="magenta") {
    return(brewer.pal(5, "RdPu"))
  }
  else if (col=="blue") {
    return(brewer.pal(5, "Blues"))
  }
  else if (col=="cyan3") {
    return(brewer.pal(5, "GnBu"))
  }
  else if (col=="red") {
    return(brewer.pal(5, "Reds"))
  }
  else if (col=="orange") {
    return(brewer.pal(5, "Oranges"))
  }
  else if (col=="black") {
    return(brewer.pal(5, "Greys"))
  }
  else if (col=="green") {
    return(brewer.pal(5, "Greens"))
  }
  else if (col=="seagreen4") {
    return(brewer.pal(5, "YlGnBu"))
  }
  else if (col=="yellowgreen") {
    return(brewer.pal(5, "YlGn"))
  }
  else {
    return(brewer.pal(5, "Set1"))
  }
}

getMatrices <- function(cov) {
  binMatrixConstructor <- setClass("binMatrix", slots = c(cov="matrix", libsize="numeric", binLength="vector", attributes="list"))
  a <- list(file = NULL, "norm.libsize" = opt$normLibSize, "norm.binlength" = opt$normBinLength, "norm.sumShape" = opt$normShapeSum, "merged" = TRUE)
  return(binMatrixConstructor(cov = cov, libsize = 0, binLength = cov[, ncol(cov)], attributes = a))
}

readBinGene <- function(m, sumFile) {
  size <- read.csv(sumFile, sep="\t", head=T)[1,1]
  
  # make first col be the names
  rownames(m) <- m[, 1]
  m <- m[, -1]
  bLength <- m[, ncol(m)]
  m <- as.matrix(m[, -ncol(m)])
  
  # create the attributes
  a <- list("file" = file, "norm.libsize" = opt$normLibSize, "norm.binlength" = opt$normBinLength, "norm.sumShape" = opt$normShapeSum, "merged" = FALSE)
  return(binMatrixConstructor(cov = m, libsize = size, binLength = bLength, attributes = a))
}

                          

storeCoverageFiles <- function(bedgraphtable_subset, dirpath, conditions) {
  binFiles <- paste(dirpath, list.files(path=dirpath, pattern=paste0("*", opt$bins, ".coverage.csv")), sep="/")
  binFiles_length <<- length(readLines(binFiles[1]))
  
  bedgraphtable_subset$basename_pos <- trimws(gsub(".bedgraph", "", basename(bedgraphtable_subset$BEDGRAPH_POS)))
  content <- unique(bedgraphtable_subset$BEDGRAPH_NEG)
  if (length(content) == 1 && is.na(content)) {
    bedgraphtable_subset$basename_neg <- NA
  }
  else {
    bedgraphtable_subset$basename_neg <- trimws(gsub(".bedgraph", "", basename(bedgraphtable_subset$BEDGRAPH_NEG)))
  }
  bedgraphtable_subset$basename_pos <- ifelse(grepl("anti", bedgraphtable_subset$CONDITION), paste0(bedgraphtable_subset$basename_pos, "_anti"), bedgraphtable_subset$basename_pos)
  bedgraphtable_subset$basename_neg <- ifelse(grepl("anti", bedgraphtable_subset$CONDITION), paste0(bedgraphtable_subset$basename_neg, "_anti"), bedgraphtable_subset$basename_neg)
  
  shapeOrderStore <- list()
  
  # loop over all samples and store them, reorder it based on file names
  for(file in binFiles) {
    filename <- gsub(paste0("_", opt$bedname, ".*coverage.csv"), "", basename(file))
    print(paste0("loading '", filename, "'..."))
    tmp <- subset(bedgraphtable_subset, bedgraphtable_subset$basename_pos == filename | bedgraphtable_subset$basename_neg == filename)
    cond <- tmp$CONDITION[1]
    
    sumFile <- paste0(dirname(file), "/", filename, ".sum.csv")
    
    # create map, store it
    entries <- length(shapeOrderStore[[cond]])
    if(entries == 0) shapeOrderStore[[cond]] <- list()
    
    shapeOrderStore[[cond]][[entries+1]] <- readBinTable(file, sumFile)
    
    entries <- length(shapeStoreReps[[filename]])
    if(entries == 0) shapeStoreReps[[filename]] <- list()
    shapeStoreReps[[filename]][[entries+1]] <- readBinTable(file, sumFile)
  }
  print(paste0(names(shapeOrderStore), " shapeorderstor"))
  # group replicates
  for(cond in conditions) {
    replicates <- shapeOrderStore[[cond]]
    print(paste0("grouping: [cond: '", cond, "']; n=", length(replicates)))
    
    # create map
    entries <- length(groupedDataStore[[cond]])
    if(entries == 0) groupedDataStore[[cond]] <- list()
    
    # store it	
    groupedDataStore[[cond]][[entries+1]] <<- groupMatrices(replicates, mean)
  }
  
  reps <<- names(shapeStoreReps)
  for(r in reps) {
    replicates <- shapeStoreReps[[r]]
    print(paste("grouping: [rep: '", r, "']; n=", length(replicates), sep=""))
    groupedData <- groupMatrices(replicates, mean)
    
    # create map
    entries <- length(groupedStoreReps[[r]])
    if(entries == 0) groupedStoreReps[[r]] <- list()
    
    # store it	
    groupedStoreReps[[r]][[entries+1]] <<- groupedData
  }
}
resetLists <- function() {
  shapeStoreGenes <<- list()
  shapeStoreGeneIds <<- list()
  counterGenes <<- 1
  shapeOrderStoreIds <<- list()
  groupedDataStore <<- list()
  groupedStore <<- list()
  groupedCondStore <<- list()
  gene <<- list()
  names <<- c()
  
  shapeStoreReps <<- list()
  groupedStoreReps <<- list()
  dfGeneListCollection <<- list()
  entries <<- 0
}


metashapeplot <- function(plotfile, condition, exp) {
  metaDataList <- c()
  names_tmp <- c()
  pal <- c()
  num_uniq_conds <- length(unique(conditions))
  vec_col_uniq <- createPalette(num_uniq_conds+3, c("#ff0000", "#00ff00", "#0000ff"))
  tmp_col_vec <- c()
  for (cond in conditions) {
        configcol <- config[cond, "COLOUR"]
        if (! is.na(configcol)) {
                tmp_col_vec <- c(tmp_col_vec, configcol)
        }
  }
  counter <- 1
  metaDataList <- NULL
  print(paste0("going ", conditions, " through metadatalist"))
  print(length(conditions))
  for (cond in conditions) {
    configcol <- config[cond, "COLOUR"]
    if (is.na(configcol)) {
        for (stored in tmp_col_vec) {
                coldist <- color_distance(stored, vec_col_uniq[counter])
                print(paste0("dist of random to stored ", coldist))
                if (coldist < 20) {
                        counter <- counter + 1
                }
        }
        configcol <- vec_col_uniq[counter]
    }
    pal <- c(pal, configcol)

    metaDataList <- c(metaDataList, groupedDataStore[[cond]])
    print(paste0(length(groupedDataStore[[cond]]), " groupeddatastore ", cond))
    condname <- ifelse(is.na(config[cond, "NAME"]), paste0("'", cond, "'"), config[cond, "NAME"])
    names_tmp <- c(names_tmp, condname)
    print(paste0("length names ", length(names_tmp), "   length data ", length(metaDataList)))
    counter <- counter + 1
  }
  names(metaDataList) <- names_tmp
  

  pdf(plotfile, width=6, height=6)
  lc <- ifelse(length(conditions) == 2, 2, ifelse(length(conditions)==3, 3, 4))	
  
  plotGroup(metaDataList, performWilcoxTest=FALSE, palette=c("seagreen4", "blue"), normBylibSize = opt$normLibSize, normByShapeSum = opt$normShapeSum, aggregateFun = opt$aggregateFUN, normByBinlength = opt$normBinLength, title =paste0("(n=",binFiles_length-1, " genes)"), labels=bodyLabels, ylab=ylabel, legendPos = "top", verticalLinesOnTicks = linesOnTicks, fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, fixedLabelsEnd = fixedLabelsEnd, showLegend = TRUE, showXLab = TRUE, showYLab = TRUE, legendNCol = 2, cexLegend = 0.6, cexAxis = 1.0, shapeLineWidth = 2.5, cexLab = 1.0, cexTitle = 1.0, legendInset = 0.01, legendXIntersp = 1.1, legendYIntersp = 0.8, lineYlab =  3, lineTitle = 1, cexLegendSymbol = 1.1, legendSpacingBetween=40)
  dev.off()
}


print("loaded functions")


########### iterate directories ##########
bedgraphtable <- read.csv(opt$bedgraphTable, header=TRUE, sep="\t", quote="")
factors <- unique(bedgraphtable$FACTOR)
factors <- factors[!is.na(factors) & factors != ""]
experiments <- unique(bedgraphtable$EXPERIMENT)
experiments <- experiments[!is.na(experiments)]

clustertype <- "x"
vec_of_clusters_without_readentries <- c()


for (fac in factors) {
  for (exp in experiments) {
    if (is.null(opt$factor) || fac %in% factor_array) {
      datatable_subset <- subset(bedgraphtable, bedgraphtable$FACTOR == fac & bedgraphtable$EXPERIMENT == exp)
      datatable_subset <- subset(datatable_subset, datatable_subset$USE == "YES")
      if (nrow(datatable_subset) > 0) {
        print(paste0("processing experiment ", exp, ", factor ", fac, "..."))
      } else {
        next
      }
      
      dirpath <- paste0(opt$coverageFiles, "/", paste0("exp", exp, "-", fac))
      conditions <- unique(datatable_subset$CONDITION)
      replicates <- c(basename(datatable_subset$BEDGRAPH_POS), basename(datatable_subset$BEDGRAPH_POS))
      replicates <- replicates[!is.na(replicates)]
      replicates <- trimws(gsub(".bedgraph", "", replicates))
      
      resetLists()

      storeCoverageFiles(datatable_subset, dirpath, conditions)
      
      n <- ifelse(is.null(opt$plotname), paste0(dirpath, "/2-scaled_metagene_", fac, "_", exp, ".pdf"), paste0(dirpath, "/", opt$plotname, "_scaled_metagene_", fac, "_", exp, ".pdf"))
      metashapeplot(n, conditions, exp)
    }
  }
}

##########################################

# we are done!
endNice(opt$confirmRun2EndFile, list('scaledMetashapeOutputFolder' = opt$coverageFiles))
