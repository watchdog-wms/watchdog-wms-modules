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
                 'bedgraphTable',   't', 1, "character", # bedgraph table
		 'experiment',      'x', 2, "character", # experiment type (ATAC/PRO/H2AZ)
                 'bedname',         'z', 1, "character", # bedfile name
                 'factor',          'f', 2, "character", # optional factor
                 'normLibSize',     'l', 0, "logical",   # norm or not
                 'normBinLength',   'n', 0, "logical",   # norm or not
                 'normShapeSum',    's', 0, "logical",   # norm or not
		 'wilcox',          'w', 0, "logical",   # perform wilxoc test or not
                 'aggregateFUN',    'a', 1, "character", # aggregate function mean/median
                 'metaFrame',       'm', 1, "numeric", # frame of metageneplot, +-TSS
                 'bins',            'b', 1, "numeric",   # number of bins
                 'plotname',        'p', 1, "character",   # plotname prefix
		 'clusterPositions','q', 2, "character", # list with cluster2wt peak position for ablines
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
if (! is.null(opt$clusterPositions)) {
	clust2ablinepos <- read.csv(opt$clusterPositions, sep="\t", header=FALSE)
}
if (is.null(opt$experiment)) {
	ylabel <- "occupancy"
} else if (opt$experiment == "PRO") {
	ylabel <- "PRO-seq profile"
} else if (opt$experiment == "ATAC") {
	ylabel <- "ATAC-seq profile"
} else if (opt$experiment == "H2AZ") {
	ylabel <- "H2A.Z occupancy profile"
} else if (opt$experiment == "ATAC-PRO" || opt$experiment == "PRO-ATAC") {
        ylabel <- "ATAC-seq/PRO-seq profile"
} else if (opt$experiment == "ATAC-H2AZ" || opt$experiment == "H2AZ-ATAC") {
	ylabel <- "ATAC-seq/H2A.Z occupancy profile"	
} else if (opt$experiment == "PRO-H2AZ" || opt$experiment == "H2AZ-PRO") {
	ylabel <- "PRO-seq/H2A.Z occupancy profile"
} else if (opt$experiment == "cRNA") {
	ylabel <- "cRNA-seq profile"
} else if (opt$experiment == "dRNA") {
	ylabel <- "dRNA-seq profile"
} else if (opt$experiment == "CDK7") {
	ylabel <- "CDK7 profile"
} else if (opt$experiment == "NELF") {
	ylabel <- "PRO-seq profile"
} else {
	ylabel <- "occupancy"
}

if (is.null(opt$wilcox)) {
	performwilcox <- FALSE
} else if (opt$wilcox == TRUE || opt$wilcox == "true") {
	performwilcox <- TRUE
} else {
	performwilcox <- FALSE
}

###################### PARAMETER CHECK #####################
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
fixedLabelsStart <- paste0(seq(-opt$metaFrame, opt$metaFrame, by=750), "bp")
fixedLabelsStart[fixedLabelsStart == "0bp"] <- "TSS"
#fixedLabelsStart[fixedLabelsStart == "0bp"] <- "pause site"
fixedLabelsStartTotalBins <- opt$bins
fixedLabelsEndTotalBins <- opt$bins
middleline <- FALSE
verticalAblines <- ifelse(middleline, c(5), c(0))  #0 weil dann wird nichts gezeichnet außerhalb der scala


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
  size <- read.csv(sumFile, sep="\t", head=T)[1,1]   #dec=","
  
  # make first col be the names
  rownames(m) <- m[, 1]
  m <- m[, -1]
  bLength <- m[, ncol(m)]
  m <- as.matrix(m[, -ncol(m)])
  
  # create the attributes
  a <- list("file" = file, "norm.libsize" = opt$normLibSize, "norm.binlength" = opt$normBinLength, "norm.sumShape" = opt$normShapeSum, "merged" = FALSE)
  return(binMatrixConstructor(cov = m, libsize = size, binLength = bLength, attributes = a))
}



storeCoverageFilesSingle <- function(bedgraphtable_subset, dirpath, conditions) {
  binFiles <- paste(dirpath, list.files(path=dirpath, pattern=paste0("*", opt$bins, ".coverage.csv")), sep="/")
  num_genes <<- length(readLines(binFiles[1]))
  
  bedgraphtable_subset$basename_pos <- trimws(gsub(".bedgraph", "", basename(bedgraphtable_subset$BEDGRAPH_POS)))
  content <- unique(bedgraphtable_subset$BEDGRAPH_NEG)
  if (length(content) == 1 && is.na(content)) {
    bedgraphtable_subset$basename_neg <- NA
  }
  else {
    bedgraphtable_subset$basename_neg <- trimws(gsub(".bedgraph", "", basename(bedgraphtable_subset$BEDGRAPH_NEG)))
  }
  
  bedgraphtable_subset$basename_pos <- ifelse(grepl("anti", bedgraphtable_subset$CONDITION), 
                            paste0(bedgraphtable_subset$basename_pos, "_anti"), bedgraphtable_subset$basename_pos)
  bedgraphtable_subset$basename_neg <- ifelse(grepl("anti", bedgraphtable_subset$CONDITION), 
                            paste0(bedgraphtable_subset$basename_neg, "_anti"), bedgraphtable_subset$basename_neg)
  
  shapeOrderStore <- list()
  
  # loop over all samples and store them, reorder it based on file names
  for(file in binFiles) {
    filename <- gsub(paste0("_", opt$bedname, ".*coverage.csv"), "", basename(file))
    print(paste0("loading '", filename, "'..."))
    tmp <- subset(bedgraphtable_subset, bedgraphtable_subset$basename_pos == filename |
                    bedgraphtable_subset$basename_neg == filename)
    cond <- tmp$CONDITION[1]
    
    sumFile <- paste0(dirname(file), "/", filename, ".sum.csv")
    data <- read.csv(file, sep="\t", head=T)   #dec="," for if file has values xxx,yy anstatt xxx.yy
    for(i in 1:nrow(data)) {
      m <- as.data.frame(data[i,])
      geneId <- toString(m$ID)
      binDataGene <- readBinGene(m, sumFile)
      shapeStoreGenes[[counterGenes]] <- binDataGene
      shapeStoreGeneIds[[counterGenes]] <- geneId
      counterGenes <- counterGenes + 1
    }
    # create map
    entries <- length(shapeOrderStore[[cond]])
    if(entries == 0) shapeOrderStore[[cond]] <- list()
    entries <- length(shapeOrderStoreIds[[cond]])
    if(entries == 0) shapeOrderStoreIds[[cond]] <- list()
    # store it
    shapeOrderStore[[cond]][[entries+1]] <- shapeStoreGenes
    shapeOrderStoreIds[[cond]][[entries+1]] <- shapeStoreGeneIds
    counterGenes <- 1
  }
  print("bin files gelesen......")
  # group replicates
  for(cond in conditions) {
    replicates <- shapeOrderStore[[cond]]
    #print(replicates)
    ids <- shapeOrderStoreIds[[cond]]
    print(paste0("grouping: [cond: '", cond, "']; n=", length(replicates)))
    for(i in (1:length(replicates[[1]]))) {
      for(j in (1:length(replicates))) {
	gene[j] <- replicates[[j]][[i]]
      }
      names <- c(names, ids[[j]][[i]])
      groupedData <- groupMatrices(gene, mean)
      groupedStore[[i]] <- groupedData
      meta <- getSummedShape(groupedData, normBylibSize = opt$normLibSize, normByBinlength = opt$normBinLength, normLibsize = 10^9, normByShapeSum = opt$normShapeSum, aggregateFun = opt$aggregateFUN)
      meta <- meta@attributes[["meta"]]
      if(i == 1) {
        metaData <- as.data.frame(meta)
      } else {
        metaData <- cbind(metaData, meta)
      }
      gene <- list()
    }
    colnames(metaData) <- names
    names <- c()
    # store it	
    groupedDataStore[[cond]] <<- metaData
    groupedCondStore[[cond]] <<- groupedStore
  }
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
  bedgraphtable_subset$basename_pos <- ifelse(grepl("anti", bedgraphtable_subset$CONDITION), 
                          paste0(bedgraphtable_subset$basename_pos, "_anti"), bedgraphtable_subset$basename_pos)
  bedgraphtable_subset$basename_neg <- ifelse(grepl("anti", bedgraphtable_subset$CONDITION), 
                          paste0(bedgraphtable_subset$basename_neg, "_anti"), bedgraphtable_subset$basename_neg)
  
  shapeOrderStore <- list()
  
  # loop over all samples and store them, reorder it based on file names
  for(file in binFiles) {
    filename <- gsub(paste0("_", opt$bedname, ".*coverage.csv"), "", basename(file))
    print(paste0("loading '", filename, "'..."))
    tmp <- subset(bedgraphtable_subset, bedgraphtable_subset$basename_pos == filename | 
                    bedgraphtable_subset$basename_neg == filename)
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
  
  # group replicates
  for(cond in names(shapeOrderStore)) {
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


metagene_quantiles <- function(plotfile, conditions) {
  pdf(plotfile, width=9, height=6)
  for (cond in conditions) {
    
    #configname <- parse(text=config[cond, "NAME"])
    pal <- getColPal(config[cond, "COLOUR"])
    
    g <- ncol(groupedDataStore[[cond]])
    
    print(cond)
    data <- apply(groupedDataStore[[cond]], 1, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
    quantileList <- list()
    
    for(i in 1:nrow(data)) {
      quantileList[[rownames(data)[i]]] <- getMatrices(t(as.matrix(data[i,])))
    }
    

    tmp <- paste0(" (n = ", num_genes, " genes)'")
    condname <- ifelse(is.na(config[cond, "NAME"]), paste0("'", cond, "'"), config[cond, "NAME"])
    configname_added <- parse(text=gsub("'$", tmp, condname))
    plotGroup(quantileList, ylab="occupancy", asIsYlab=TRUE, labels="", palette=pal, title=configname_added, legendPos = "topleft", cexLegend=0.7, cexLegendSymbol = 1.5, legendInset = 0.06, legendXIntersp = 0.8, legendYIntersp = 1.0, cexTitle=0.9, normBylibSize = opt$normLibSize, normByShapeSum = opt$normShapeSum, aggregateFun = opt$aggregateFUN, normByBinlength = opt$normBinLength, fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, cexAxis=1.0, cexLab = 1.0)
  }
  dev.off()
}

metagene_conds <- function(plotfile, conditions, exp) {
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
    condname <- ifelse(is.na(config[cond, "NAME"]), paste0("'", cond, "'"), config[cond, "NAME"])
    condname <- gsub("X", " ", condname)
    names_tmp <- c(names_tmp, condname)
    counter <- counter + 1
  }
  names(metaDataList) <- names_tmp
  pdf(plotfile, width=6, height=6)
  
  lc <- ifelse(length(conditions) == 2, 2, ifelse(length(conditions)==3, 3, 4))
  plotGroup(metaDataList, performWilcoxTest=performwilcox, ylab=ylabel, asIsYlab=TRUE, palette=pal, labels="", title=paste0("(n = ", binFiles_length-1, " genes)"), cexTitle=1.0, normBylibSize = opt$normLibSize, normByShapeSum = opt$normShapeSum, aggregateFun = opt$aggregateFUN, normByBinlength = opt$normBinLength, legendPos = "top", fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, cexAxis=1.0, legendNCol = lc, legendSpacingBetween = 22, cexLegend = 0.5, shapeLineWidth = 2.5, cexLab = 1.0, legendInset = 0.01, legendXIntersp = 1.1, legendYIntersp = 0.8, cexLegendSymbol = 1.1)
  
  dev.off()
}

metagene_reps <- function(plotfile, bedgraphtable_subset) {
  bedgraphtable_subset <- bedgraphtable_subset[order(bedgraphtable_subset$CONDITION),]
  c <- 1
  temp <- bedgraphtable_subset$CONDITION[1]
  repnames <- c()
  metaDataList <- c()
  cv <- c()
  repcount <- 1
  counter <- 1
  num_uniq_conds <- length(unique(bedgraphtable_subset$CONDITION))
  vec_col_uniq <- createPalette(num_uniq_conds+5, c("#ff0000", "#00ff00", "#0000ff"))
  tmp_col_vec <- c()
  for (l in 1:nrow(bedgraphtable_subset)) {
	configcol <- config[bedgraphtable_subset$CONDITION[l], "COLOUR"]
    	if (! is.na(configcol)) {
        	tmp_col_vec <- c(tmp_col_vec, configcol)
    	}
  }
  for (l in 1:nrow(bedgraphtable_subset)) {
    if (bedgraphtable_subset$CONDITION[l] != temp) {
      c <- c + 1
      temp <- bedgraphtable_subset$CONDITION[l]
      repcount <- 1
    } 
    curr_rep <- trimws(gsub(".bedgraph", "", basename(bedgraphtable_subset$BEDGRAPH_POS[l])))
    if (grepl("anti", bedgraphtable_subset$CONDITION[l])) {
      curr_rep <- paste0(curr_rep, "_anti")
      curr_rep <- gsub("_pos_", "_neg_", curr_rep)
    }
    metaDataList <- c(metaDataList, groupedStoreReps[[curr_rep]])
    
    configcol <- config[bedgraphtable_subset$CONDITION[l], "COLOUR"]
    if (is.na(configcol)) {
	for (stored in tmp_col_vec) {
                coldist <- color_distance(stored, vec_col_uniq[counter])
                print(paste0("dist of random to stored ", coldist))
                if (coldist < 20) {
                        counter <- counter + 1
                }
        }
    	configcol <- vec_col_uniq[counter]
    	print(configcol)
    }
    if (repcount > 1) {
	colortemp <- configcol
    	for (i in 2:repcount) {
		configcol <- lighten(configcol)
	}
	if ((255-col2rgb(configcol)[1]) < 70 || (255-col2rgb(configcol)[2]) < 70 || (255-col2rgb(configcol)[3]) < 70) {
		print("color lightening reached white, starting over from random color")
		configcol <- rgb(sample(c(1:255),1), sample(c(1:255),1), sample(c(1:255),1), max=255)
	}
	cv <- c(cv, configcol)
    }
    else {
	cv <- c(cv, configcol)
    }
    
    tmp <- paste0(" rep. ", repcount, "'")
    condname <- ifelse(is.na(config[bedgraphtable_subset$CONDITION[l], "NAME"]),  paste0("'", bedgraphtable_subset$CONDITION[l], "'"), config[bedgraphtable_subset$CONDITION[l], "NAME"])
    configname_added <- gsub("'$", tmp, condname)
    configname_added <- gsub("X", " ", configname_added)
    repnames <- c(repnames, configname_added)
    
    repcount <- repcount + 1
    counter <- counter + 1
  }
  names(metaDataList) <- repnames
  
  lc <- ifelse(c == 2, 2, ifelse(c == 3, 3, 4))
  numconds <- length(unique(bedgraphtable_subset$CONDITION))
  if (length(cv)/numconds >= 4) {
	cv <- unname(createPalette(length(cv), "#ff0000"))
  } 
  
  pdf(plotfile, width=6, height=6)
  plotGroup(metaDataList, ylab=ylabel, asIsYlab=TRUE, performWilcoxTest=performwilcox, palette=cv, labels="", title=paste0("(n = ", binFiles_length-1, " genes)"), cexTitle=1.0, normBylibSize = opt$normLibSize, normByShapeSum = opt$normShapeSum, aggregateFun = "mean", normByBinlength = opt$normBinLength, legendPos = "top", fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, cexAxis=1.0, legendNCol = lc, legendSpacingBetween = 22, cexLegend = 0.5, shapeLineWidth = 2.5, cexLab = 1.0, legendInset = 0.01, legendXIntersp = 1.1, legendYIntersp = 0.8, cexLegendSymbol = 1.1)
  dev.off()
}


metagene_quantiles_clustered <- function(plotfile, condtions, genelist) {
  pdf(plotfile, width=9, height=6)
  for (cond in conditions) {
    configcol <- config[cond, "COLOUR"]
    pal <- getColPal(configcol)
    
    dfGeneListCollection <<- list()
    qlc <- getQuantileClusterPlots(cond, genelist, plotfile)
    c <- 1
    for (m in qlc) {
      tmp <- paste0(", ", clustertype, c, " (n = ", nrow(dfGeneListCollection[[c]]), " genes)'")
      condname <- ifelse(is.na(config[cond, "NAME"]), paste0("'", cond, "'"), config[cond, "NAME"])
      configname_added <- parse(text=gsub("'$", tmp, condname))
      
      if (is.null(m)) {
	print(paste0("no entries to plot in cluster ", c))
      	c <- c + 1
        next
      }

      plotGroup(m, labels="", ylab="occupancy", asIsYlab=TRUE, palette=pal, title=configname_added, legendPos = "topleft", cexLegend=0.5, cexLegendSymbol = 1.3, legendInset = 0.06, legendXIntersp = 0.8, legendYIntersp = 1.0, cexTitle=1.0, normBylibSize = opt$normLibSize, normByShapeSum = opt$normShapeSum, aggregateFun = opt$aggregateFUN, normByBinlength = opt$normBinLength, fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, cexAxis=1.0, cexLab = 1.0)
      c <- c + 1
    }
  }
  dev.off()
}
getQuantileClusterPlots <- function(dataName, genelist, plotfile) {
  df <- groupedDataStore[[dataName]]

  for(i in 1:ncol(df)) {
    for (j in 1:max(genelist[,2])) {
      if(toString(colnames(df)[i]) %in% subset(genelist, genelist[,2] == j)[,1]) {
        ifelse(length(dfGeneListCollection) < j, dfGeneListCollection[[j]] <<- (df[,i]),
               dfGeneListCollection[[j]] <<- rbind(dfGeneListCollection[[j]], (df[,i])))
      }
    }
  }

  dataGeneListCollection <- list()
  for (k in 1:max(genelist[,2])) {
    if (length(dfGeneListCollection[[k]]) == 0) {
	    print(paste0("empty gene list without reads for cluster ", k))
	    if (!(k %in% vec_of_clusters_without_readentries)) {
		vec_of_clusters_without_readentries <<- c(vec_of_clusters_without_readentries, k)
	    }
	    dataGeneListCollection[[k]] <- dfGeneListCollection[[k]]
	    cluster_genes <- subset(genelist, genelist[,2]==k)
	    cluster_genes <- as.data.frame(cluster_genes)
	    colnames(cluster_genes)[1] <- paste0("genes_cluster_", k)
	    write.table(cluster_genes[,1], paste0(dirname(plotfile), "/", strsplit(strsplit(basename(plotfile), "_")[[1]][length(strsplit(basename(plotfile), "_")[[1]])], "\\.")[[1]][1], "_NOPLOT_cluster_genes.txt"), sep="\t", col.names=TRUE, quote=FALSE, row.names=FALSE)
	    genes_in_cluster <- c()
	    genes_not_in_cluster <- c()
	    for (tempi in 1:ncol(df)) {
		if (toString(colnames(df)[tempi]) %in% cluster_genes[, 1]) {
			genes_in_cluster <- c(genes_in_cluster, toString(colnames(df)[tempi]))
		}
	    	else {
			genes_not_in_cluster <- c(genes_not_in_cluster, toString(colnames(df)[tempi]))
		}
	    }
	    write.table(genes_in_cluster, paste0(dirname(plotfile), "/", strsplit(strsplit(basename(plotfile), "_")[[1]][length(strsplit(basename(plotfile), "_")[[1]])], "\\.")[[1]][1], "_NOPLOT_in_cluster.txt"), sep="\t", quote=FALSE, row.names=FALSE)
	    write.table(genes_not_in_cluster, paste0(dirname(plotfile), "/", strsplit(strsplit(basename(plotfile), "_")[[1]][length(strsplit(basename(plotfile), "_")[[1]])], "\\.")[[1]][1], "_NOPLOT_not_in_cluster.txt"), sep="\t", quote=FALSE, row.names=FALSE)
    }
    else {
    	dataGeneList <- apply(dfGeneListCollection[[k]], 2, quantile, probs=c(0.1,0.25,0.5,0.75,0.9))
    	dataGeneListCollection[[k]] <- dataGeneList
    }
  }
  
  tmp2 <- function(datagenelist) {
    quantilelist <- list()
    for (i in 1:nrow(datagenelist)) {
      quantilelist[[rownames(datagenelist)[i]]] <- getMatrices(t(as.matrix(datagenelist[i,])))
    }
    return(quantilelist)
  }
  
  quantileListCollection <- list()
  for (m in 1:max(genelist[,2])) {
    if (length(dataGeneListCollection[[m]]) == 0) {
	quantileListCollection[[m]] <- dataGeneListCollection[[m]]
    }
    else {
	quantileListCollection[[m]] <- tmp2(dataGeneListCollection[[m]])
    }
  }
  return(quantileListCollection)
}

metagene_conds_clustered <- function(plotfile, conditions, exp, genelist, fac) {
  metaDataList <- c()
  names_tmp <- c()
  num_conds_uniq <- length(unique(conditions))
  vec_col_uniq <- createPalette(num_conds_uniq+3, c("#ff0000", "#00ff00", "#0000ff"))
  tmp_col_vec <- c()
  for (cond in conditions) {
	configcol <- config[cond, "COLOUR"]
    	if (! is.na(configcol)) {
        	tmp_col_vec <- c(tmp_col_vec, configcol)
    	}
  }
  pal <- c()
  counter <- 1
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
    
    condname <- ifelse(is.na(config[cond, "NAME"]), paste0("'", cond, "'"), config[cond, "NAME"])
    condname <- gsub("X", " ", condname)
    names_tmp <- c(names_tmp, condname)
    counter <- counter + 1
  }
  names(metaDataList) <- names_tmp
  
  lc <- ifelse(length(conditions) == 2, 2, ifelse(length(conditions)==3, 3, 4))

  pdf(plotfile, width=6, height=6)
  #clusternames <- c("increased PI", "slightly reduced PI", "strongly reduced PI")
  for (i in 1:max(genelist[,2])) {
    geneList <- list(subset(genelist, genelist[,2] == i)[,1])
    names(geneList) <- c(paste0(clustertype, i))
    #names(geneList) <- c(clusternames[i])
    if (i %in% vec_of_clusters_without_readentries) {
	print(paste0("skipping cluster ", i))
    	next
    }
    #labs <- fixedLabelsStart
    if (!is.null(opt$clusterPositions)) {
	fixedLabelsStart <- paste0(seq(-opt$metaFrame, opt$metaFrame, by=60), "bp")
	fixedLabelsStart[seq(2,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(3,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(4,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(5,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(6,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(7,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(8,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(9,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[seq(10,length(fixedLabelsStart), 10)] <- ""
	fixedLabelsStart[49] <- ""
	fixedLabelsStart[52] <- ""
	fixedLabelsStart[100] <- ""
	fixedLabelsStart[101] <- "3000bp"
	fixedLabelsStart[51] <- "TSS"
	verticalAblines <- c()
    	verticalAblines <- c(verticalAblines, (clust2ablinepos[which(clust2ablinepos == i), 2]))
	verticalAblines <- verticalAblines[!is.na(verticalAblines)]
	print(fixedLabelsStart)
    }

    if (!is.null(opt$clusterPositions)) {
	plotMergedShape(metaDataList, performWilcoxTest=performwilcox,  asIsYlab=TRUE, palette=pal, geneList, name=ylabel, legendPos = "top", bodyLabels = "", fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, showXLab = TRUE, showLegend = TRUE, showYLab = TRUE, legendNCol = lc, cexAxis = 1.0, legendSpacingBetween = 20, cexLegend = 0.6, shapeLineWidth = 2.5, cexLab = 1.0, legendInset = 0.01, legendXIntersp = 1.1, legendYIntersp = 0.8, cexLegendSymbol = 1.0, cexTitle = 1.0, verticalLinesOnTicks=verticalAblines)
} else {
	plotMergedShape(metaDataList, performWilcoxTest=performwilcox,  asIsYlab=TRUE, palette=pal, geneList, name=ylabel, legendPos = "top", bodyLabels = "", fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, showXLab = TRUE, showLegend = TRUE, showYLab = TRUE, legendNCol = lc, cexAxis = 1.0, legendSpacingBetween = 20, cexLegend = 0.6, shapeLineWidth = 2.5, cexLab = 1.0, legendInset = 0.01, legendXIntersp = 1.1, legendYIntersp = 0.8, cexLegendSymbol = 1.0, cexTitle = 1.0)
}
    
  }
  dev.off()
}

metagene_reps_clustered <- function(plotfile, bedgraphtable_subset, genelist) {
  bedgraphtable_subset <- bedgraphtable_subset[order(bedgraphtable_subset$CONDITION),]
  c <- 1
  temp <- bedgraphtable_subset$CONDITION[1]
  repnames <- c()
  cv <- c()
  metaDataList <- c()
  repcount <- 1
  counter <- 1
  num_conds_uniq <- length(unique(bedgraphtable_subset$CONDITION))
  vec_col_uniq <- createPalette(num_conds_uniq+5, c("#ff0000", "#00ff00", "#0000ff"))
  tmp_col_vec <- c()
  for (l in 1:nrow(bedgraphtable_subset)) {
	configcol <- config[bedgraphtable_subset$CONDITION[l], "COLOUR"]
    	if (! is.na(configcol)) {
        	tmp_col_vec <- c(tmp_col_vec, configcol)
    	}
  }
  for (l in 1:nrow(bedgraphtable_subset)) {
    if (bedgraphtable_subset$CONDITION[l] != temp) {
      c <- c + 1
      temp <- bedgraphtable_subset$CONDITION[l]
      repcount <- 1
    } 
    curr_rep <- trimws(gsub(".bedgraph", "", basename(bedgraphtable_subset$BEDGRAPH_POS[l])))
    if (grepl("anti", bedgraphtable_subset$CONDITION[l])) {
      curr_rep <- paste0(curr_rep, "_anti")
      curr_rep <- gsub("_pos_", "_neg_", curr_rep)
    }
    metaDataList <- c(metaDataList, groupedStoreReps[[curr_rep]])
    
    configcol <- config[bedgraphtable_subset$CONDITION[l], "COLOUR"]
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

    if (repcount > 1) {
        colortemp <- configcol
        for (i in 2:repcount) {
                configcol <- lighten(configcol)
        }
        cv <- c(cv, configcol)
    }
    else {
        cv <- c(cv, configcol)
    }
    
    tmp <- paste0(" rep. ", repcount, "'")
    condname <- ifelse(is.na(config[bedgraphtable_subset$CONDITION[l], "NAME"]), paste0("'", bedgraphtable_subset$CONDITION[l], "'"), config[bedgraphtable_subset$CONDITION[l], "NAME"])
    configname_added <- gsub("'$", tmp, condname)
    configname_added <- gsub("X", " ", configname_added)
    
    repnames <- c(repnames, configname_added)
    repcount <- repcount + 1
    counter <- counter + 1
  }
  names(metaDataList) <- repnames
  
  lc <- ifelse(c == 2, 2, ifelse(c == 3, 3, 4))
  numconds <- length(unique(bedgraphtable_subset$CONDITION))
  if (length(cv)/numconds >= 4) {
        cv <- unname(createPalette(length(cv), "#ff0000"))
  }
  
  pdf(plotfile, width=6, height=6)
  for (i in 1:max(genelist[,2])) {
    geneList <- list(subset(genelist, genelist[,2] == i)[,1])
    names(geneList) <- c(paste0(clustertype, i))
    if (i %in% vec_of_clusters_without_readentries) {
        print(paste0("skipping cluster ", i))
        next
    }

    plotMergedShape(metaDataList, asIsYlab=TRUE, geneList, performWilcoxTest=performwilcox, palette=cv, name=ylabel, legendPos = "top", bodyLabels = "", fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, showXLab = c(TRUE), showLegend = c(TRUE), showYLab = TRUE, legendNCol = lc, cexAxis = 1.0, legendSpacingBetween = 22, cexLegend = 0.5, shapeLineWidth = 2.5, cexLab = 1.0, legendInset = 0.01, legendXIntersp = 1.1, legendYIntersp = 0.8, cexLegendSymbol = 1.1, cexTitle = 1.0)
  }
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
      #datatable_subset$CONDITION <- gsub("_", " ", datatable_subset$CONDITION)
      #conditions <- unique(datatable_subset$CONDITION)
      #conditions <- c("mock", "1h", "2h", "4h", "8h")
      #conditions <- c("mock", "2h", "1h")
      #conditions <- c("XRN1-mock", "mock", "XRN1-8h", "8h")
      print(conditions)

      replicates <- c(basename(datatable_subset$BEDGRAPH_POS), basename(datatable_subset$BEDGRAPH_POS))
      replicates <- trimws(gsub(".bedgraph", "", replicates))
      
      storeCoverageFilesSingle(datatable_subset, dirpath, conditions)
      
      n <- ifelse(is.null(opt$plotname), paste0(dirpath, "/1-metashapeplot_quantiles_", fac, "_", exp, ".pdf"), 
                  paste0(dirpath, "/", opt$plotname, "_metashapeplot_quantiles_", fac, "_", exp, ".pdf"))
      metagene_quantiles(n, conditions)
      
      resetLists()
      storeCoverageFiles(datatable_subset, dirpath, conditions)
      
      n <- ifelse(is.null(opt$plotname), paste0(dirpath, "/2-metagene_per_condition_", fac, "_", exp, ".pdf"), 
                  paste0(dirpath, "/", opt$plotname, "_metagene_per_condition_", fac, "_", exp, ".pdf"))
      metagene_conds(n, conditions, exp)
      n <- ifelse(is.null(opt$plotname), paste0(dirpath, "/3-metagene_per_replicate_", fac, "_", exp, ".pdf"), 
                  paste0(dirpath, "/", opt$plotname, "_metagene_per_replicate_", fac, "_", exp, ".pdf"))
      metagene_reps(n, datatable_subset)
      
      if (! is.null(opt$genelist)) {
        for (l in genelist_array) {
          clustertype <- ifelse(grepl("random", l), "random ", ifelse(grepl("ward", l), "cluster ", "CPM top "))
      	  clustertype <- ifelse(grepl("random", l), "random ", ifelse(grepl("ward", l), "cluster ", ifelse(grepl("pi", l), "PI ", "CPM top ")))
          
	  
	  if (file.exists(l) || file.exists(paste0(dirpath, "/", l))) {
                print("reading own cluster file")
                gl <- c()
                if (file.exists(paste0(dirpath, "/", l))) {
                        gl <- read.csv(paste0(dirpath, "/", l), sep="\t", header=FALSE)
			print(paste0("reading own cluster file: ", dirpath, "/", l))
                }
                else {
                        gl <- read.csv(l, sep="\t", header=FALSE)
			print(paste0("reading own cluster file: ", l))
                }
                gl[,2] <- as.numeric(as.character(gl[,2]))
                glname <- gsub(".cluster.csv", "", basename(l))
                resetLists()
                storeCoverageFilesSingle(datatable_subset, dirpath, conditions)
                n <- ifelse(is.null(opt$plotname),
                        paste0(dirpath, "/1-metashapeplot_quantiles_cluster_", fac, "_", exp, "_", glname, ".pdf"),
                        paste0(dirpath, "/", opt$plotname, "_metashapeplot_quantiles_cluster_", fac, "_", exp, "_", glname, ".pdf"))
                metagene_quantiles_clustered(n, conditions, gl)
                resetLists()
                storeCoverageFiles(datatable_subset, dirpath, conditions)
                n <- ifelse(is.null(opt$plotname),
                        paste0(dirpath, "/2-metagene_per_condition_cluster_", fac, "_", exp, "_", glname, ".pdf"),
                        paste0(dirpath, "/", opt$plotname, "_metagene_per_condition_cluster_", fac, "_", exp, "_", glname, ".pdf"))
                metagene_conds_clustered(n, conditions, exp, gl, fac)
                n <- ifelse(is.null(opt$plotname),
                        paste0(dirpath, "/3-metagene_per_replicate_cluster_", fac, "_", exp, "_", glname, ".pdf"),
                        paste0(dirpath, "/", opt$plotname, "_metagene_per_replicate_cluster_", fac, "_", exp, "_", glname, ".pdf"))
                metagene_reps_clustered(n, datatable_subset, gl)
          }
	  else if (file.exists(paste0("/mnt/raidtmp/weiss/metageneplots/workflow/pro_all/exp4-PROseq_Pol2/", l))) {
            #nothing
          } else {
            print("gene list file does not exist...")
          }
        }
      }
    }
  }
}

##########################################

# we are done!
endNice(opt$confirmRun2EndFile, list('generateMetagenePlotsOutputFolder' = opt$coverageFiles))
