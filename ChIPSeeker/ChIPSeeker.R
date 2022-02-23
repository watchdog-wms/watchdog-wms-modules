library('getopt')
# include functions
args <- commandArgs(trailingOnly = FALSE)
functionsFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../sharedUtils/R/functions.R', sep = '')
source(functionsFileName)

# names to keep of annotated peaks
keepCols <- c("seqnames", "start", "end", "width", "strand", "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "geneId", "distanceToTSS", "ENSEMBL", "SYMBOL", "GENENAME") 
renameCols <- c("chr", "peakStart", "peakEnd", "peakWidth", "peakStrand", "annotation", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", "geneId", "distanceToTSS", "ENSEMBL", "SYMBOL", "GENENAME") 


# options to parse
spec <- matrix(c('bedFiles', 'b', 1, 'character',
		 'annoDb', 'a', 1, 'character',
		 'txdb', 't', 1, 'character',
		 'outputDir', 'o', 1, 'character',
		 'promotorUpstream', 'u', 2, 'integer',
		 'promotorDownstream', 'd', 2, 'integer',
		 'resample', 'r', 2, 'integer',
		 'conf', 'c', 2, 'integer',
		 'version', 'v', 2, 'logical',
		 'confirmRun2EndFile', 'z', 1, 'character'), ncol = 4, byrow = T)


# parse the parameters
opt <- getopt(spec)

if(!is.null(opt$version) && opt$version) {
	cat(paste(as.character(packageVersion('ChIPseeker')), '\n', sep = ''))
	endNice(opt$confirmRun2EndFile)
}

# define required libs
requiredPackages <- c('GenomicFeatures', 'ChIPseeker')

# load some libs
endErrorMissingPackage(requiredPackages)
loadPackages(requiredPackages)

# set default values
if(is.null(opt$promotorUpstream)) { opt$promotorUpstream <- 3000 }
if(is.null(opt$promotorDownstream)) { opt$promotorDownstream <- 3000 }
if(is.null(opt$resample)) { opt$resample <- 1000 }
if(is.null(opt$conf)) { opt$conf <- 0.95 }

# verify that input parameters are given
verifyInputNotEmpty(opt, 'bedFiles', TRUE)
verifyInputNotEmpty(opt, 'annoDb', TRUE)
verifyInputNotEmpty(opt, 'txdb', TRUE)
verifyInputNotEmpty(opt, 'outputDir', TRUE)

# make new txdb friom gtf or gff3
if(file.exists(opt$txdb)) {
	print(paste('[INFO] Converting gff/gtf3 annotation to  \'', opt$txdb, '\' to txdb.', sep = ''))
	txdb <- makeTxDbFromGFF(opt$txdb)
} else { 
	# try to load it as library
	print(paste('[INFO] Loading txdb library \'', opt$txdb, '\'...', sep = ''))
	endErrorMissingPackage(opt$txdb)
	library(opt$txdb, character.only = TRUE)
	txdb <- get(ls(paste('package:', opt$txdb, sep = ''))[[1]])
}

# check if input files are there
peakFiles <- unlist(strsplit(opt$bedFiles, ','))
for(peakFile in peakFiles) {
	verifyFileExistence(peakFile, TRUE)
}

# try to install packages if missing
if(!testPackage(opt$annoDb)) {
	source('https://bioconductor.org/biocLite.R')
	biocLite(opt$annoDb)
}
# load organism database
library(opt$annoDb, character.only = TRUE)

# create output folder
if(!dir.exists(opt$outputDir)) {
	dir.create(opt$outputDir, recursive = TRUE, showWarnings = FALSE)
	if(!dir.exists(opt$outputDir)) {
		print(paste('[ERROR] Failed to create output folder \'', opt$outputDir, '\'.', sep = ''))
		quit('no', 17) # writing failed
	}
}

# get promotors from txdb annotation
promoter <- getPromoters(TxDb = txdb, upstream = opt$promotorUpstream, downstream = opt$promotorDownstream)
promoterRange <- c(-opt$promotorUpstream, opt$promotorDownstream)

# make plots for each sample
bnl <- c()
for(peakFile in peakFiles) {
	bn <- basename(peakFile)
	print(paste('[INFO] processing \'', bn ,'\'...', sep = ''))
	bn <- gsub('\\.[^\\.]+$', '', bn)
	bn <- gsub('\\.GPS_events$', '', bn)
	bnl <- c(bnl, bn)

	# read data and annotate them
	peak <- readPeakFile(peakFile, head=F)
	if(length(peak) > 25) {
		pdf(getOutputFile(opt$outputDir, paste(bn, '.pdf', sep = '')))

		tagMatrix <- getTagMatrix(peak, windows = promoter)
		peakAnno <- annotatePeak(peak, tssRegion = promoterRange, TxDb = txdb, annoDb = opt$annoDb)

		# make plots
		print(covplot(peak))
		if(is.matrix(tagMatrix) && sum(tagMatrix) > 0) {
			tagHeatmap(tagMatrix, xlim = promoterRange, color = 'red', title = 'tags at promotors')
			print(plotAvgProf(tagMatrix, xlim = promoterRange, xlab ='genomic region (5\'->3\')', ylab = 'read count frequency', conf = opt$conf, resample = opt$resample))
		}
		print(plotAnnoBar(peakAnno))
		upsetplot(peakAnno)
		print(plotDistToTSS(peakAnno, title = 'Distribution of peaks relative to TSS'))

		# write peaks to disk
		peakData <- as.data.frame(peakAnno)
		peakData <- peakData[, keepCols]
		colnames(peakData) <- renameCols
		write.table(peakData, file = getOutputFile(opt$outputDir, paste(bn, '.csv', sep = '')), quote=F, sep="\t", row.names=F)
		dev.off() 
	}
}

# plots for all samples
if(length(peakFiles) > 1) {
	print('[INFO] processing all samples together...')
	pdf(getOutputFile(opt$outputDir, 'allSamples.pdf'))
	# read samples and annotate them
	tagMatrixList <- lapply(peakFiles, getTagMatrix, windows = promoter)
	peakAnnoList <- lapply(peakFiles, annotatePeak, TxDb = txdb, tssRegion = promoterRange, verbose = FALSE)
	names(tagMatrixList) <- bnl
	names(peakAnnoList) <- bnl

	# make plots
	tagHeatmap(tagMatrixList, xlim = promoterRange, color=NULL)
	print(plotAvgProf(tagMatrixList, xlim = promoterRange, conf = opt$conf, resample = opt$resample, facet = 'row'))
	print(plotAnnoBar(peakAnnoList))
	print(plotDistToTSS(peakAnnoList))
	
	dev.off()
}

# we are done!
endNice(opt$confirmRun2EndFile, list('ChIPSeekerOutputFolder' = opt$outputDir))

