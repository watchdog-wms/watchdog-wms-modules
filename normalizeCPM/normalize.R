# load lib
library(getopt)

opt <- NULL
# options to parse
spec <- matrix(c('sums',   's', 1, "character", # dir of coverage files
                 'outputFile',         'o', 1, "character",   # number of k cluster
                 'counts',         'c', 1, "character",   # number of k cluster
                 'confirmRun2EndFile', 'e', 1, "character"
), ncol=4, byrow=T)

# parse the parameters
opt <- getopt(spec)


###################### PARAMETER CHECK ######################
if(is.null(opt$sums)) {
  print("[ERROR] Path to input file not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$counts)) {
  print("[ERROR] Path to input file not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$outputFile)) {
  print("[ERROR] Path to output file not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}



args <- commandArgs(trailingOnly = FALSE)
functionsFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../../modules/sharedUtils/R/functions.R', sep = '')
source(functionsFileName)
binFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../../modules/binGenome/R/binGenome.lib.R', sep = '')
source(binFileName)


#########
counts <- read.csv(opt$counts, sep="\t", header=TRUE)
sums <- read.csv(opt$sums, sep="\t", header=FALSE)


normedTable <- counts[,1]
for (col in 2:(ncol(counts))) {
  sample <- colnames(counts)[col]
  entry <- subset(sums, sums[,1] == sample)
  normedTable <- cbind(normedTable, counts[, col]/entry[1,2]*1000000)
}
colnames(normedTable) <- colnames(counts)

write.table(normedTable, opt$outputFile, sep="\t", quote=FALSE, row.names=FALSE)


################################

# we are done!
endNice(opt$confirmRun2EndFile, list('normedCounts' = opt$outputFile))


