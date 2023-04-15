
options(warn=-1)

# load lib
library(getopt)
library(plyr)
library(tidyr)
library(gtools)
#library(PTXQC)

opt <- NULL
# options to parse
spec <- matrix(c('out',              'o', 1, "character",   # output dir
		              'testCondition',       't', 1, "character",   # test condition
                 'controlCondition',     'c', 1, "character",   # control condition
		              'sampleAnnotation',    's', 1, "character",  #sample annotation file
                 'confirmRun2EndFile', 'e', 1, "character"
), ncol=4, byrow=T)

# parse the parameters
opt <- getopt(spec)


###################### PARAMETER CHECK ######################
if(is.null(opt$testCondition)) {
  print("[ERROR] test condition not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$out)) {
  print("[ERROR] Path to input file not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$controlCondition)) {
  print("[ERROR] control condition not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$sampleAnnotation)) {
  print("[ERROR] sample annotation file not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}


args <- commandArgs(trailingOnly = FALSE)
functionsFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../../modules/sharedUtils/R/functions.R', sep = '')
source(functionsFileName)


#########

c1_c2 <- paste0(opt$controlCondition, "_", opt$testCondition)
c2_c1 <- paste0(opt$testCondition, "_", opt$controlCondition)


samples <- read.csv(opt$sampleAnnotation, sep="\t", header=TRUE)

dexseq <- read.csv(paste0(opt$out, "/", c2_c1, "/DEXSeq.csv"), sep="\t", header=TRUE)
lfc <- colnames(dexseq)[grepl("log2fold", colnames(dexseq))]
wantedcols <- c("groupID", "padj", lfc, "genomicData.start", "genomicData.end", "genomicData.strand", "transcripts")
dexseq <- dexseq[, colnames(dexseq) %in% wantedcols]

dexseq <- subset(dexseq, !grepl("filled", dexseq$transcripts))



amssregs <- read.csv(paste0(opt$out, "/windows/amss.regions"), sep="\t", header=TRUE)
amssregs <- subset(amssregs, !is.na(amssregs$start_seq_idx))
res <- c()
for (line in 1:nrow(amssregs)) {
  print(amssregs$window[line])
  base <- paste0(opt$out, "/windows/", amssregs$window[line])
  counts <- read.csv(paste0(base, "/normalized_aggregated_counts.tsv"), sep="\t", header=TRUE)
  t <- subset(counts, counts$start>=amssregs$start_seq_idx[line] & counts$start<=amssregs$end_seq_idx[line])
  if (amssregs$dir[line]==c1_c2) {
    i <- t[which(t[,4]==max(t[,4])),2]
    i2 <- t[which(t[,4]==max(t[,4])),4][1]
  } else {
    i <- t[which(t[,6]==max(t[,6])),2]
    i2 <- t[which(t[,6]==max(t[,6])),6][1]
  }
  maxx <- i[ceiling(length(i)/2)]
  maxy <- round(i2, 5)

  
  s <- amssregs$start_seq_idx[line]
  e <- amssregs$end_seq_idx[line]
  d <- subset(dexseq, dexseq$genomicData.start==s & dexseq$genomicData.end==e &
                 dexseq$groupID==amssregs$window[line])[1,c(2,3,6)]

  d <- unname(d)
  rownames(d) <- NULL
  if (nrow(subset(dexseq, dexseq$genomicData.start==s & dexseq$genomicData.end==e &
                  dexseq$groupID==amssregs$window[line])) == 0) {
    d <- t(c("", "", ""))
  }
  z <- cbind(unname(amssregs[line,]), maxx, maxy, d)
  res <- rbind(res, z)

}
res <- as.data.frame(res)
colnames(res) <- c("start_seq_idx", "end_seq_idx", "score", "AUC", "length",
                   "direction", "window", "max_x", "max_y", "padj", colnames(dexseq)[3], "strand")
res <- res[,c(7,6,1,2,3,5,8,9,10,11,12)]
res$direction <- gsub("_", ">", res$direction)

print("first for loop done")



p <- paste0(opt$out, "windows/")
windows <- list.dirs(path=p, recursive=FALSE)
fin <- c()
for (gene in windows) {
  reg <- strsplit(gene, "/")[[1]]
  reg <- reg[length(reg)]
  subscores <- subset(res, res$window==reg)
  if (nrow(subscores)==0) {
    next
  }
  print(gene)
  reps <- list.files(path=paste0(gene, "/counts/"), pattern="counts")
  reps <- mixedsort(reps)
  strangenametempfix <- FALSE
  if (grepl("filter", reps)) {
    reps <- gsub("filtered\\.", "", reps)
    strangenametempfix <- TRUE
  }
  
  #SOLVE IT BY sampleAnnotation must be ordered (cond1 a-c, cond2 a-c)
  conds <- subset(samples, samples$condition==opt$controlCondition)
  numhalf <- nrow(conds)
  
  direction <- subscores$direction[1]
  repdf <- c()
  for (line in 1:(nrow(samples)/2)) {
  	rep <- samples$sample[line]
  	pair <- samples$sample[line+numhalf]
    if (strangenametempfix) {
      	c1 <- read.csv(paste0(gene, "/counts/", gsub("\\.counts", "\\.filtered\\.counts", rep, ".counts")), sep="\t", header=FALSE)
      	c2 <- read.csv(paste0(gene, "/counts/", gsub("\\.counts", "\\.filtered\\.counts", pair, ".counts")), sep="\t", header=FALSE)
    } else {
     		c1 <- read.csv(paste0(gene, "/counts/", rep, ".counts"), sep="\t", header=FALSE)
      	c2 <- read.csv(paste0(gene, "/counts/", pair, ".counts"), sep="\t", header=FALSE)
    }
    c1$V2 <- (c1$V2/sum(c1$V2))*100
    c2$V2 <- (c2$V2/sum(c2$V2))*100
    c <- merge(c1, c2, by="V1")
    colnames(c) <- c("idx", "c1", "c2")
    idx <- c$idx
    c$c1 <- as.numeric(c$c1)
    c$c2 <- as.numeric(c$c2)
    if (gsub(">", "_", direction)==c1_c2) {
      		c$diff <- c$c1 - c$c2
    } else {
      		c$diff <- c$c2 - c$c1
    }
    repdf <- cbind(repdf, c[,4])
  }
  repdf <- as.data.frame(repdf)
  repdf$idx <- idx
  
  for (line in 1:nrow(subscores)) {
    s <- subscores$start_seq_idx[line]
    e <- subscores$end_seq_idx[line]
    repdfsub <- subset(repdf, repdf$idx>=s & repdf$idx<=e)
    nreps <- ncol(repdfsub)-1
    results <- c()
    for (i in 1:nreps) {
      sumdiff <- sum((repdfsub[,i]))
      results <- c(results, sumdiff)
    }
    l <- cbind(reg, s, e, t(results))
    fin <- rbind(fin, l)
  }
}
fin <- as.data.frame(fin)
colnames(fin)[4:ncol(fin)] <- paste0("rep", 1:(ncol(fin)-3))
colnames(fin)[1:3] <- c("window", "start_seq_idx", "end_seq_idx")

s <- merge(res, fin, by=c("window", "start_seq_idx", "end_seq_idx"))
s <- s[,c(1,4,2,3,5:ncol(s))]
write.table(s, paste0(opt$out, "/scores.tsv"), sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

print("done writing scores table")

################################

# we are done!
endNice(opt$confirmRun2EndFile, list('out' = opt$out))


