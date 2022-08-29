# load lib
library(RColorBrewer)
library(gplots)
library(getopt)
library(plotrix)

opt <- NULL
# options to parse
spec <- matrix(c('coverageFiles',   'c', 1, "character", # dir of coverage files
                 'cluster',         'k', 1, "numeric",   # number of k cluster
                 'bedgraphTable',   't', 1, "character", # bedgraph table
                 'bedname',         'z', 1, "character", # bedfile name
                 'factor',          'f', 2, "character", # optional factor
                 'normLibSize',     'l', 0, "logical",   # norm or not
                 'normBinLength',   'n', 0, "logical",   # norm or not
                 'normShapeSum',    's', 0, "logical",   # norm or not
                 'bins',            'b', 1, "numeric",   # number of bins
                 'aggregateFUN',    'a', 1, "character", # aggregate function mean/median
                 'cpm',             'm', 1, "character", # counts per million table
                 'plotname',        'p', 2, "character", # plotname prefix
                 'confirmRun2EndFile', 'e', 1, "character"
), ncol=4, byrow=T)

# parse the parameters
opt <- getopt(spec)


factor_array <- c()
if (! is.null(opt$factor)) {
  factor_array <- strsplit(opt$factor, ",")[[1]]
}

#factor_array <- strsplit(opt$factor, ",")[[1]]

###################### PARAMETER CHECK ######################
if(is.null(opt$bedgraphTable)) {
  print("[ERROR] Path to bedgraph table not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$bedname)) {
  print("[ERROR] bedname not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$cluster)) {
  print("[ERROR] Number of cluster not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$coverageFiles)) {
  print("[ERROR] Path to coverage files not set")
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


#########
current_cluster <- c()
groupedDataStore <- list()

clusterStore <- list()
dfGeneListCollection <- list()
binFiles_length <- 0

clustfiles <- ""

shapeStoreGenes <- list()
counterGenes <- 1
shapeStoreGeneIds <- list()
shapeOrderStore <- list()
shapeOrderStoreIds <- list()
gene <- list()
groupedStore <- list()
names <- c()
shapeStoreReps <- list()
groupedStoreReps <- list()

################################
storeCoverageFilesSingle <- function(bedgraphtable_subset, dirpath, conditions) {
  binFiles <- paste(dirpath, list.files(path=dirpath, pattern=paste0("*", opt$bins, ".coverage.csv")), sep="/")
  
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
  
  
  # loop over all samples and store them, reorder it based on file names
  for(file in binFiles) {
    filename <- gsub(paste0("_", opt$bedname, ".*coverage.csv"), "", basename(file))
    print(paste0("loading '", filename, "'..."))
    tmp <- subset(bedgraphtable_subset, bedgraphtable_subset$basename_pos == filename | bedgraphtable_subset$basename_neg == filename)
    cond <- tmp$CONDITION[1]
    
    sumFile <- paste0(dirname(file), "/", filename, ".sum.csv")
    data <- read.csv(file, sep="\t", head=T)
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
  # group replicates
  
  for(cond in conditions) {
    replicates <- shapeOrderStore[[cond]]
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
  }
}

cluster_wardD2 <- function(mergeMatrixNormed, outfile) {
  d <- dist(mergeMatrixNormed, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  cluster <- cutree(hc, k=opt$cluster)
  clusterSideBar <<- colRamp[cluster]
  
  h <- plotHeamap2(mergeMatrixNormed, clusterSideBar, c(101, 201, 301), c("black", "black", "black"), conditions)
  clusters_per_genes_ordered<-cluster[h$rowInd]
  clusters_in_order<-unique(clusters_per_genes_ordered) 
  clusters_in_order_rev <- rev(clusters_in_order)
  clusters_new<-rep(0, length(clusters_in_order))
  for(i in 1:length(clusters_in_order)) {
    clusters_new[clusters_in_order_rev[i]] <- i
  }
  geneNames <- names(cluster)
  cluster <- clusters_new[cluster]
  names(cluster) <- geneNames
  
  df <- c()
  for(i in 1:opt$cluster) {
    temp <- subset(cluster, cluster == i)
    dataFrame <- cbind(names(temp), temp)
    df <- rbind(df, dataFrame)
  }
  rownames(df) <- NULL
  current_cluster <<- df
  clustfiles <<- paste0(clustfiles, paste0(outfile, "/ward-D2.cluster.csv", sep=""), sep=",")
  write.table(df, paste(outfile, "/ward-D2.cluster.csv", sep=""), sep="\t", row.names = FALSE, col.names = FALSE, quote=F)
}

cluster_random <- function(mergeMatrixNormed, outfile) {
  genes <- rownames(mergeMatrixNormed)
  genes <- sample(genes)
  df <- c()
  chunksize <- floor(length(genes)/opt$cluster)
  for (i in 1:opt$cluster) {
    start <- (i-1)*chunksize + 1
    if (i == opt$cluster) {
      dat <- cbind(genes[start:length(genes)], i)
      df <- rbind(df, dat)
      break
    }
    dat <- cbind(genes[start:(i*chunksize)], i)
    df <- rbind(df, dat)
  }
  current_cluster <<- df
  clustfiles <<- paste0(clustfiles, paste0(outfile, "/random.cluster.csv", sep=""), sep=",")
  write.table(df, paste(outfile, "/random.cluster.csv", sep=""), sep="\t", row.names = FALSE, col.names = FALSE, quote=F)
}

cluster_featurecounts <- function(outfile) {
  fctable <- read.csv(opt$cpm, sep="\t", header=TRUE, quote="")
  fctable_mean <- as.data.frame(cbind(fctable$FeatureID, apply(fctable[, c(2:ncol(fctable))], 1, mean)))
  fctable_mean <- fctable_mean[order(fctable_mean$V2, decreasing=TRUE),]
  
  df <- c()
  chunksize <- floor(nrow(fctable_mean)/opt$cluster)
  for (i in 1:opt$cluster) {
    start <- (i-1)*chunksize + 1
    if (i == opt$cluster) {
      dat <- cbind(fctable_mean[start:nrow(fctable_mean), 1], i)
      df <- rbind(df, dat)
      break
    }
    dat <- cbind(fctable_mean[start:(i*chunksize), 1], i)
    df <- rbind(df, dat)
  }
  current_cluster <<- df
  clustfiles <<- paste0(clustfiles, paste0(outfile, "/fc_counts.cluster.csv", sep=""), sep=",")
  write.table(df, paste(outfile, "/fc_counts.cluster.csv", sep=""), sep="\t", row.names = FALSE, col.names = FALSE, quote=F)
}

heatmap_clustered <- function(plotfile, mergeMatrixNormed, conditions) {
  pdf(file=plotfile,  width = 6, height = 6)
  plotHeamap2(mergeMatrixNormed, clusterSideBar, c(101, 202, 303), c("black", "black", "black"), conditions)
  dev.off()
}

plotHeamap2 <- function(data, rowsidecol, colsep, sepcol, conditions) {
  breaks=seq(0,max(getMergedNormedMatrix(conditions)),0.01)
  colRamp <- colorRampPalette(c("white","yellow", "orange", "red", "darkred", "black"))
  heatmap.2(
    as.matrix(data), 
    Colv="false",
    dendrogram="row",
    trace="none",
    scale="none",
    hclustfun=function(x) hclust(x, method = "ward.D2"),
    labCol="",
    labRow="",
    xlab="",
    cexRow=0.7,
    density.info="none",
    key.xlab=NA,
    key.par=list(mar=c(3,2,3,2)),
    main="",
    breaks=breaks,
    colsep=colsep,
    sepcolor=sepcol,
    RowSideColors=rowsidecol,
    col=colRamp(length(breaks)-1)
  )
}

#ohne farben, mit mittlerer linie bei tss und ohne kante rechts
#plotHeamap2 <- function(data, rowsidecol, colsep, sepcol, conditions) {
#  breaks=seq(0,max(getMergedNormedMatrix(conditions)),0.01)
#  colRamp <- colorRampPalette(c("white","yellow", "orange", "red", "darkred", "black"))
#  heatmap.2(
#    as.matrix(data),
#    Colv="false",
#    dendrogram="row",
#    trace="none",
#    scale="none",
#    hclustfun=function(x) hclust(x, method = "ward.D2"),
#    labCol="",
#    labRow="",
#    xlab="",
#    cexRow=0.7,
#    density.info="none",
#    key.xlab=NA,
#    key.par=list(mar=c(3,2,3,2)),
#    main="",
#    breaks=breaks,
#    colsep=c(51),
#    sepcolor=c("magenta"),
#    sepwidth=c(0.3),
#    #RowSideColors=rowsidecol,
#    col=colRamp(length(breaks)-1)
#  )
#}

getMergedNormedMatrix <- function(conditions) {
  merged_mat <- c()
  for (cond in conditions) {
    tmp_mat <- as.matrix(groupedDataStore[[cond]])
    tmp_mat[is.na(tmp_mat) | is.nan(tmp_mat)] <- 0
    tmp_mat <- as.data.frame(t(tmp_mat))
    if (is.null(merged_mat)) {
      merged_mat <- tmp_mat
    }
    else {
      merged_mat <- merge(merged_mat, tmp_mat, by="row.names", all=TRUE)
      row.names(merged_mat) <- merged_mat$Row.names
      merged_mat$Row.names <- NULL
    }
  }
  merged_mat <- data.matrix(merged_mat)
  merged_mat[is.na(merged_mat) | is.nan(merged_mat)] <- 0
  mergeMatrixNormed <- sweep(merged_mat,1, apply(merged_mat, 1, max), "/")
  mergeMatrixNormed[is.nan(mergeMatrixNormed)] <- 0
  return(mergeMatrixNormed)
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



################################
bedgraphtable <- read.csv(opt$bedgraphTable, header=TRUE, sep="\t", quote="")
factors <- unique(bedgraphtable$FACTOR)
factors <- factors[!is.na(factors) & factors != ""]
experiments <- unique(bedgraphtable$EXPERIMENT)
experiments <- experiments[!is.na(experiments)]


colRamp <- (rainbow(opt$cluster))
library(Polychrome)
colRamp <- createPalette(opt$cluster,  c("#ff0000", "#00ff00", "#0000ff"))

tmp <- sapply(colRamp, color.id)
for (i in 1:length(tmp)) {
  colRamp[i] <- tmp[[i]][1]
}

for (fac in factors) {
  for (exp in experiments) {
    if (is.null(opt$factor) || fac %in% factor_array) {
      datatable_subset <- subset(bedgraphtable, bedgraphtable$FACTOR == fac & bedgraphtable$EXPERIMENT == exp)
      datatable_subset <- subset(datatable_subset, datatable_subset$USE == "YES")
      if (nrow(datatable_subset) > 0 && datatable_subset$CLUSTERING[1] == "YES") {
        print(paste0("processing experiment ", exp, ", factor ", fac, "..."))
      } else {
        next
      }
      
      dirpath <- paste0(opt$coverageFiles, "/", paste0("exp", exp, "-", fac))
      conditions <- unique(datatable_subset$CONDITION)
      replicates <- c(basename(datatable_subset$BEDGRAPH_POS), basename(datatable_subset$BEDGRAPH_POS))
      replicates <- trimws(gsub(".bedgraph", "", replicates))
      
      storeCoverageFilesSingle(datatable_subset, dirpath, conditions)
      
      normedMatrix <- getMergedNormedMatrix(conditions)
      cluster_wardD2(normedMatrix, dirpath)
      cluster_random(normedMatrix, dirpath)
      cluster_featurecounts(dirpath)
      
      n <- ifelse(is.null(opt$plotname), paste0(dirpath, "/4-ward-D2Clustering_", opt$cluster, "_cluster.pdf"), 
                  paste0(dirpath, "/", opt$plotname, "_-ward-D2Clustering_", opt$cluster, "_cluster.pdf"))
      
      heatmap_clustered(n, normedMatrix, conditions)
    }
  }
}

# we are done!
endNice(opt$confirmRun2EndFile, list('coverageFiles' = opt$coverageFiles, 'bedname' = opt$bedname, 'clusterfiles' = clustfiles))


