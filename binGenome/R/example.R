library(getopt)

opt <- NULL
# options to parse
spec <- matrix(c('outputFile', 'o', 1, "character",	# path to output PDF
		 'testDataDir', 't', 1, "character",	# path to test data folder
		 'installDir', 'i', 1, "character"	# path to folder with R scripts
		 ), ncol=4, byrow=T)

# parse the parameters
opt <- getopt(spec)

#opt$outputFile <- "/tmp/test.pdf"
#opt$installDir <- "/home/proj/software/watchdog/modules/binGenome/R/"
#opt$testDataDir <- "/home/proj/software/watchdog/modules/binGenome/test_data/"

###################### PARAMETER CHECK ######################
if(is.null(opt$outputFile)) {
	print(paste("[ERROR] Path to output file not set: '", opt$outputFile,"'", sep=""))
	quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$installDir)) {
	print(paste("[ERROR] Path to install dir of genome binner R scripts not set: '", opt$installDir,"'", sep=""))
	quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$testDataDir)) {
	print(paste("[ERROR] Path to test data not set: '", opt$testDataDir,"'", sep=""))
	quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}

if(!is.null(opt$installDir) && !file.exists(opt$installDir)) {
	print(paste("[ERROR] Install dir of genome binner R scripts does not exist: '", opt$installDir,"'", sep=""))
	quit(save = "no", status = 15, runLast = FALSE) # status = 15 <--> missing input files
}
if(!is.null(opt$testDataDir) && !dir.exists(opt$testDataDir)) {
	print(paste("[ERROR] Test data dir does not exist: '", opt$testDataDir,"'", sep=""))
	quit(save = "no", status = 15, runLast = FALSE) # status = 15 <--> missing input files
}
# create output folder if required
basedir <- dirname(opt$outputFile)
if(!dir.exists(basedir)) {
	dir.create(basedir, showWarnings = F, mode = "0755")
	if(!dir.exists(basedir)) {
		print(paste("[ERROR] Failed to create output dir '", basedir,"'.", sep=""))
		quit(save = "no", status = 17, runLast = FALSE) # status = 17 <--> writing failed
	}
}


# load lib
source(paste(opt$installDir, "/binGenome.lib.R", sep="")) 

# config for plots
inputFolder <- opt$testDataDir
removeEndings <- TRUE	# remove Gencode version endings of gene names
legendPos <- "top"	# legend position
nameYAxis <- "show case" 	# name on y-axis
metaDataNorm <- NULL	# do not norm the data
isGeneList <- TRUE	# gene lists must be mapped to maximal expressed transcripts

# labels for x-axis
fixedLabelsStart <- paste(seq(-3000, 1500, by=750), "bp", sep="")
fixedLabelsEnd <- paste(seq(-1500, 3000, by=750), "bp", sep="")
fixedLabelsStart[fixedLabelsStart == 0] <- "TSS"
fixedLabelsEnd[fixedLabelsEnd == 0] <- "TTS"
bodyLabels <- c("", paste(seq(20, 80, by=20), "%", sep=""), "")
linesOnTicks <- c(5,7,12,14)
fixedLabelsStartTotalBins <- 90
fixedLabelsEndTotalBins <- 90

# read some pre-processed input files
transGeneMapping <- read.csv(paste(inputFolder, "geneTranscriptMapping.csv", sep="/"), sep="\t", head=T, as.is=T)
maxTranscriptIDs <- scan(paste(inputFolder, "maxExpressedTranscripts.csv", sep="/"), what="character")
transGeneMapping <- transGeneMapping[transGeneMapping$transcript_id %in% maxTranscriptIDs, ]

# load binned chipSeq data
covPattern <- "*.coverage.csv"
binFiles <- paste(inputFolder, list.files(path=inputFolder, pattern=covPattern), sep="/")

# loop over all samples and store them
shapeStore <- list()
counter <- 1
for(file in binFiles) {
	name <- gsub(".coverage.+", "", basename(file))
	print(paste("loading '", name, "'...", sep=""))
	dir <- dirname(file)
	sumFile <- paste(dir, "/", name, ".sum.csv", sep="")
	binData <- readBinTable(file, sumFile)

	# get max transcripts
	binData <- filterByIDs(binData, maxTranscriptIDs)

	# store it for later use
	shapeStore[[counter]] <- binData
	counter <- counter + 1
}

# reorder it based on file names
counter <- 1
shapeOrderStore <- list()
for(file in binFiles) {
	name <- gsub("chip_", "", basename(file))
	cond <- gsub("_.+$", "", name)
	exp <- gsub(paste("^", cond, "_", sep=""), "", name)
	exp <- gsub("_.+$", "", exp)
	
	# create map
	entries <- length(shapeOrderStore[[exp]])
	if(entries == 0) shapeOrderStore[[exp]] <- list()
	entries <- length(shapeOrderStore[[exp]][[cond]])
	if(entries == 0) shapeOrderStore[[exp]][[cond]] <- list()

	# store it
	shapeOrderStore[[exp]][[cond]][[entries+1]] <- shapeStore[[counter]]
	counter <- counter + 1
}

# group replicates
groupedDataStore <- list()
expTypes <- names(shapeOrderStore)
for(exp in expTypes) {
	typeStore <- shapeOrderStore[[exp]]
	conditions <- names(typeStore)
	for(cond in conditions) {
		replicates <- typeStore[[cond]]
		print(paste("grouping: [cond: '", cond, "',  antibody: '", exp, "']; n=", length(replicates), sep=""))
		groupedData <- groupMatrices(replicates, mean)

		# create map
		entries <- length(groupedDataStore[[exp]])
		if(entries == 0) groupedDataStore[[exp]] <- list()
		entries <- length(groupedDataStore[[exp]][[cond]])
		if(entries == 0) groupedDataStore[[exp]][[cond]] <- list()
		
		# store it	
		groupedDataStore[[exp]][[cond]][[entries+1]] <- groupedData
	}
}


# create shape plots
metaDataList <- c(groupedDataStore[["PolII"]][["Ctl"]], groupedDataStore[["Ser2P"]][["Ctl"]])
names(metaDataList) <- c("PolII", "Ser2P")
geneList <- list(transGeneMapping[1:1000, "gene_id"], transGeneMapping[1001:1500, "gene_id"])
names(geneList) <- c("random set 1", "random set 2")

# with gene list; without wilcoxon test;
pdf(opt$outputFile, width=9, height=6)
plotMergedShape(metaDataList, normByShapeSum = F, aggregateFun = median, geneList, name = nameYAxis, isListGeneList = isGeneList, transGeneMapping = transGeneMapping, maxTranscripts = maxTranscriptIDs, metaDataNorm = NULL, legendPos = legendPos, bodyLabels = bodyLabels, verticalLinesOnTicks = linesOnTicks, fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, fixedLabelsEnd = fixedLabelsEnd, showLegend = c(TRUE, FALSE), showXLab = c(FALSE, TRUE), showYLab = TRUE, legendNCol = 2, cexLegend = 1.3, shapeLineWidth = 2.5, cexLab = 1.4, legendInset = 0.05, legendXIntersp = 0.8, legendYIntersp = 0.3, lineYlab =  3, lineTitle = 1, cexLegendSymbol = 2.5)
plotMergedShape(metaDataList, normByShapeSum = T, aggregateFun = median, geneList, name = nameYAxis, isListGeneList = isGeneList, transGeneMapping = transGeneMapping, maxTranscripts = maxTranscriptIDs, metaDataNorm = NULL, legendPos = legendPos, bodyLabels = bodyLabels, verticalLinesOnTicks = linesOnTicks, fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, fixedLabelsEnd = fixedLabelsEnd, showLegend = c(TRUE, FALSE), showXLab = c(FALSE, TRUE), showYLab = TRUE, legendNCol = 2, cexLegend = 1.3, shapeLineWidth = 2.5, cexLab = 1.4, legendInset = 0.05, legendXIntersp = 0.8, legendYIntersp = 0.3, lineYlab =  3, lineTitle = 1, cexLegendSymbol = 2.5)

# with wilcoxon test
plotMergedShape(metaDataList,  aggregateFun = mean, geneList, name = nameYAxis, isListGeneList = isGeneList, transGeneMapping = transGeneMapping, maxTranscripts = maxTranscriptIDs, metaDataNorm = NULL, legendPos = legendPos, bodyLabels = bodyLabels, verticalLinesOnTicks = linesOnTicks, fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, fixedLabelsEnd = fixedLabelsEnd, showLegend = c(TRUE, FALSE), showXLab = c(FALSE, TRUE), showYLab = TRUE, legendNCol = 2, cexLegend = 1.3, shapeLineWidth = 2.5, cexLab = 1.4, legendInset = 0.05, legendXIntersp = 0.8, legendYIntersp = 0.3, lineYlab =  3, lineTitle = 1, cexLegendSymbol = 2.5, performWilcoxTest = T)

# create 4 clusters based on Ser2P and drop biggest one afterwards
dataToCluster <- groupedDataStore[["Ser2P"]][["Ctl"]][[1]]
cl <- clustKmeans(dataToCluster, 4)
cll <- unlist(lapply(cl, length))
cll <- cll == max(cll)
cl <- cl[!cll]

# plot the clustering
plotMergedShape(metaDataList, cl, name = "clustering test", isListGeneList = FALSE, metaDataNorm = NULL, legendPos = legendPos, bodyLabels = bodyLabels, verticalLinesOnTicks = linesOnTicks, fixedLabelsStartTotalBins = fixedLabelsStartTotalBins, fixedLabelsEndTotalBins = fixedLabelsEndTotalBins, fixedLabelsStart = fixedLabelsStart, fixedLabelsEnd = fixedLabelsEnd, showLegend = c(TRUE, FALSE, FALSE), showXLab = c(FALSE, FALSE, TRUE), showYLab = c(FALSE, TRUE, FALSE), legendNCol = 2, cexLegend = 1.6, shapeLineWidth = 2.5, cexLab = 1.5, legendInset = 0.05, legendXIntersp = 0.8, legendYIntersp = 0.3, lineYlab =  3, lineTitle = 1, cexLegendSymbol = 3.5, performWilcoxTest = T, asIsYlab = T, cexAxis= 1.6, oma=c(6.5,1,0.25,0.25), cexTitle = 2.2, numberOfYTicks = c(0, 0.01, 0.02), minYA = 0, maxYA = 0.02, mar = c(1.5,4.5,2.5,0.5), wilcoxTestPVTransformFUN = function(x) { defaultPvalueColorTransformer(x, defaultColor = "green", pvThresholds = c(10^-20, 10^-30, 10^-40))})
dev.off()
