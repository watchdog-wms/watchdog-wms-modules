# load some libs
require(getopt)
require(DEXSeq)
require(GenomicFeatures)
require(GenomicRanges)
require(BiocParallel)

# options to parse
spec <- matrix(c('controlCondition', 'a', 1, "character",
		  'testCondition', 'b', 1, "character",
		  'countFile', 'c', 1, "character",
		  'flattedGTFAnnotation', 'd', 1, "character",
		  'sampleAnnotation', 'e', 1, "character",
		  'featureAnnotation', 'f', 1, "character",
		  'featureAnnotationID', 'g', 1, "character",
		  'featureAnnotationName', 'h', 1, "character",
		  'excludeSamples', 'i', 1, "character",
		  'pValueCutoff', 'j', 1, "numeric",
		  'minKeepReads', 'k', 1, "numeric",
		  'output', 'm', 1, "character",
		  'threads', 'n', 1, "numeric",
		  'confirmRun2EndFile', 'p', 1, "character"), ncol=4, byrow=T)

# parse the parameters
opt = getopt(spec);
# we do no more checking for arguments because we expect that all checking is done before!

vs <- paste(opt$testCondition, opt$controlCondition, sep="_")
# excludeSamples some of the samples.
if(!is.null(opt$excludeSamples)) 
{
	opt$excludeSamples <- sort(unlist(strsplit(opt$excludeSamples, ",")))
}

# build expressionSet
exprs <- as.matrix(read.table(opt$countFile, header=TRUE, sep="\t", row.names=1, as.is=TRUE, check.names = FALSE, stringsAsFactors = FALSE))
pData <- read.table(opt$sampleAnnotation, row.names=1, header=TRUE, sep="\t", stringsAsFactors = FALSE)

if(!is.null(opt$featureAnnotation)) {
	features <- read.table(opt$featureAnnotation, header=TRUE, sep="\t", as.is=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
}

# remove samples, if wished
if(!is.null(opt$excludeSamples)) 
{
	opt$excludeSamples <- opt$excludeSamples[pData[opt$excludeSamples, ] == opt$controlCondition | pData[opt$excludeSamples, ] == opt$testCondition]
	exprs <- exprs[, ! colnames(exprs) %in% opt$excludeSamples]
	vs <- paste(vs, "_excludeSamples_", paste(opt$excludeSamples, collapse="_"), sep="")
}

# later for DE gene test testCondition / controlCondition is computed though controlCondition should be control
opt$output <- paste(opt$output, "/", vs, sep="")
dir.create(opt$output, showWarnings = F, mode = "0755")
# removed not used rows
pDataTmp <- pData
pDataTmp$sample <- rownames(pDataTmp)
pDataTmp <- pDataTmp[pDataTmp$condition == opt$controlCondition | pDataTmp$condition == opt$testCondition, ]

pData <- as.data.frame(pDataTmp[, "condition"])
rownames(pData) <- pDataTmp$sample

S1 <- pDataTmp[pDataTmp$condition == opt$controlCondition, "sample"]
S2 <- pDataTmp[pDataTmp$condition == opt$testCondition, "sample"]
exprs <- exprs[, colnames(exprs) %in% c(S1, S2)]
colnames(pData) <- "condition"
if(class(pData$condition) == "factor") {
	pData$condition <- droplevels(pData$condition)
}

pData <- as.data.frame(pData[colnames(exprs), ])
rownames(pData) <- colnames(exprs)
colnames(pData)[1] <- "condition"

if(nrow(exprs) == 0) {
	print(paste("[ERROR] No gene passed the low expression filtering!", sep=""))
	quit(save = "no", status = 14, runLast = FALSE) # status = 14 <--> invalid arguments
}

# prepare the stuff for use with featureCount input
aggregates <- read.delim(opt$flattedGTFAnnotation, stringsAsFactors = FALSE, header = FALSE)
colnames(aggregates) <- c("chr", "source", "class", "start", "end", "ex", "strand", "ex2", "attr")
aggregates$strand <- gsub("\\.", "*", aggregates$strand)
aggregates <- aggregates[which(aggregates$class == "exonic_part"), ]
aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
transcripts <- strsplit(transcripts, "\\+")
exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr)
exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, end = aggregates$end), strand = aggregates$strand)
names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":") # changed, to work with newer version: removed E after :
names(transcripts) <- names(exoninfo)

# build DEXSeqDataSetgff
groupID <- sapply(strsplit(rownames(exprs), ":"), "[", 1)

tmp <- exprs
tmp <- as.data.frame(unlist(apply(tmp, 1, sum)))
tmp$names <- groupID
colnames(tmp)[1] <- "count"
agg <- aggregate(cbind(count) ~ names, tmp, FUN = "sum")
exprs <- exprs[groupID %in% agg[agg$count >= opt$minKeepReads*ncol(exprs), "names"], ]

groupID <- sapply(strsplit(rownames(exprs), ":"), "[", 1)
featureID <- sapply(strsplit(rownames(exprs), ":"), "[", 2)

x <- sapply(strsplit(names(exoninfo), ":"), "[", 1)

# remove features with too long names
found <- groupID %in% x & !is.na(featureID)
groupID <- groupID[found]
featureID <- featureID[found]
exoninfo <- exoninfo[x %in% groupID]
transcripts <- transcripts[x %in% groupID]
exprs <- exprs[rownames(exprs) %in% paste(groupID, featureID, sep=":"), ]

# ensure correct sorting of all of the stuff
exprs <- exprs[sort(rownames(exprs)), ]
groupID <- sapply(strsplit(rownames(exprs), ":"), "[", 1)
featureID <- sapply(strsplit(rownames(exprs), ":"), "[", 2)
exoninfo <- exoninfo[sort(names(exoninfo)), ]
transcripts <- transcripts[sort(names(transcripts))]

pData$condition <- factor(pData$condition)
dxd <- DEXSeqDataSet(exprs, pData, design= ~ sample + exon + condition:exon, featureID, groupID, featureRanges=exoninfo, transcripts=transcripts)

# go!!!
BPPARAM=MulticoreParam(workers=opt$threads)
dxd<-estimateSizeFactors(dxd)
dxd<-estimateDispersions(dxd,BPPARAM=BPPARAM)
plotDispEsts(dxd)

dxd<-testForDEU(dxd,BPPARAM=BPPARAM)
dxd<-estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
dxr1<-DEXSeqResults(dxd)
dframeTest <- is(dxr1, "DFrame")

# change name for better reading
if(!is.null(opt$featureAnnotation)) {
	features$gene_id <- features[, opt$featureAnnotationID]
	mapID <- function(x, features) {
		# test, if we must split it at ':'
		id_exon <- unlist(strsplit(x, ":"))
		r <- features[features$gene_id == id_exon[1], opt$featureAnnotationName]
		if(length(r) == 0) {
			return(x)
		}
		else {
			if(length(id_exon) == 2)
			{
				r <- paste(r, ":", id_exon[2], sep="")
			}
			return(r)
		}
	}
	mapIDList <- function(x, features) {
		return(paste(unlist(lapply(x, mapID, features)), collapse="+"))
	}
	if(dframeTest) {
		dxr1@listData$groupID <- unlist(lapply(strsplit(dxr1$groupID, "\\+"), mapIDList, features))
		names(dxr1@listData$genomicData) <- unlist(lapply(strsplit(names(dxr1$genomicData), "\\+"), mapIDList, features))
	} else {
		dxr1$groupID <- unlist(lapply(strsplit(dxr1$groupID, "\\+"), mapIDList, features))
		names(dxr1$genomicData) <- unlist(lapply(strsplit(names(dxr1$genomicData), "\\+"), mapIDList, features))
	}	
	rownames(dxr1) <- unlist(lapply(strsplit(rownames(dxr1), "\\+"), mapIDList, features))
}

# write it to the disk
DEXSeqHTML(dxr1, FDR=opt$pValueCutoff,color=c("#FF000080","#0000FF80"), path=opt$output, BPPARAM=BPPARAM)
if(dframeTest) {
	dxr1@listData$transcripts <- unlist(lapply(dxr1$transcripts, paste, collapse=";")) 
} else {
	dxr1$transcripts <- lapply(dxr1$transcripts, paste, collapse=";")
}

if(dframeTest) {
write.table(dxr1@listData, file = paste(opt$output, "DEXSeq.csv", sep="/"), append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE);
} else{
write.table(dxr1, file = paste(opt$output, "DEXSeq.csv", sep="/"), append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}

# write a file that we know, the script run to its end
if(!is.null(opt$confirmRun2EndFile)) {
	file.create(opt$confirmRun2EndFile, showWarnings = FALSE)
}

