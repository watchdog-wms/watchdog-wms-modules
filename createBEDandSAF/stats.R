args = commandArgs(trailingOnly=TRUE)

dir <- args[1]
prefix <- args[2]

numPos <- read.csv(paste0(dir, "mapping_numPositions.stats"), sep="\t", header=TRUE)
numTranscripts <- read.csv(paste0(dir, "mapping_numTranscripts.stats"), sep="\t", header=TRUE)
numGenes <- read.csv(paste0(dir, "mapping_numGenes.stats"), sep="\t", header=TRUE)



library(ggplot2)
require(data.table)

print(file.exists(paste0(dir, "mapping.info")))
if (file.exists(paste0(dir, "mapping.info"))) {
  mapping <- read.csv(paste0(dir, "mapping.info"), sep="\t", header=FALSE)
  entries <- nrow(mapping)
  mapping$V7 <- as.numeric(as.character(mapping$V7))
  mapping <- subset(mapping, !is.na(mapping$V7))

  ggplot(mapping, aes(`V7`)) + geom_histogram(bins=50) + ylab("# occurrences") + 
    xlab("distance to nearest annotated TSS") +
    ggtitle(paste0(entries, " entries")) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE))
  ggsave(paste0(dir, prefix, "_distance2TSS_hist.png"))
  

  ggplot(mapping, aes(`V7`)) + stat_ecdf(geom="step") + 
    ylab("# occurrences") + 
    xlab("distance to nearest annotated TSS") +
    ggtitle(paste0(entries, " entries")) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE))
  ggsave(paste0(dir, prefix, "_distance2TSS_cum_orig.png"))

  ggplot(mapping, aes(`V7`)) + stat_ecdf(geom="step") + 
    ylab("# occurrences") + 
    xlab("distance to nearest annotated TSS") +
    ggtitle(paste0(entries, " entries")) +
    coord_cartesian(xlim=c(0,10000))
  ggsave(paste0(dir, prefix, "_distance2TSS_cum.png"))
  
  
  #b <- list.files(path=dir, pattern=".bed")[1]
  #print(b)
  #bed <- read.csv(paste0(dir, "/", b), sep="\t", header=FALSE)
  #bed$peak <- bed$V2 + 3000
  #colnames(bed)[1] <- "chr"
  #colnames(bed)[6] <- "strand"
  #colnames(bed)[7] <- "peak"
  #colnames(mapping)[1] <- "chr"
  #colnames(mapping)[2] <- "strand"
  #colnames(mapping)[3] <- "peak"
  #m <- merge(bed, mapping, by=c("chr", "strand", "peak"))
  #m2 <- subset(m, m$V7 <= 100)
  #ggplot(m2, aes(`V7`)) + stat_ecdf(geom="step") + 
   # ylab("# occurrences") + 
   # xlab("distance to nearest annotated TSS") +
   # ggtitle(paste0(nrow(m2), " entries")) +
    #scale_x_continuous(labels = function(x) format(x, scientific = FALSE))
  #ggsave(paste0(dir, prefix, "_distance2TSS_cum_limiteddist.png"))
}



png(paste0(dir, prefix, "_numPosPerRegion.png"), width = 500, height = 500)
ggplot(numPos, aes(num_pos)) + geom_bar() + ylab("# occurrences") + xlab("# unique positions / region") +
  ggtitle(paste0(nrow(numPos), " entries"))
dev.off()


png(paste0(dir, prefix, "_numTranscriptsPerRegion.png"), width = 500, height = 500)
clear_mapping <- subset(numTranscripts, numTranscripts$num_transcripts == 1)
percent_unambigious <- nrow(clear_mapping)/nrow(numTranscripts)

ggplot(numTranscripts, aes(num_transcripts)) + geom_bar() + ylab("# occurrences") + 
  xlab("# unique transcripts / region") +
  ggtitle(paste0(nrow(numTranscripts), " entries")) + 
  annotate("text", label=paste0(round(percent_unambigious*100, 2), "% clear mapping"), x=4, y=3000)
dev.off()



png(paste0(dir, prefix, "_numGenesPerRegion.png"), width = 500, height = 500)
clear_mapping <- subset(numGenes, numGenes$num_genes == 1)
percent_unambigious <- nrow(clear_mapping)/nrow(numGenes)

ggplot(numGenes, aes(num_genes)) + geom_bar() + ylab("# occurrences") + xlab("# unique genes / region") +
  ggtitle(paste0(nrow(numGenes), " entries")) +
  annotate("text", label=paste0(round(percent_unambigious*100, 2), "% mapped"), x=2, y=3000)
dev.off()
