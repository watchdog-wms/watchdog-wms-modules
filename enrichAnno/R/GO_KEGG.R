# load some libs
library(getopt)
library(clusterProfiler)
library(pathview)
library(KEGGREST)

endNice <- function(file) {
	# write a file that we know, the script run to its end
	if(!is.null(file)) {
		file.create(file, showWarnings = FALSE)
	}
	quit("no")
}

# options to parse
spec <- matrix(c('testFiles', 't', 1, "character",
		  'backgroundFile', 'b', 1, "character",
		  'orgDB', 'd', 1, "character",
		  'keggDBName', 'k', 1, "character",
		  'pValueCutoff', 'p', 1, "numeric",
		  'output', 'o', 1, "character",
		  'plotKegg', 'n', 0, "logical",
                  'foldchangeCol', 'f', 0, "character",
		  'confirmRun2EndFile', 'c', 1, "character"), ncol=4, byrow=T)

# parse the parameters
opt = getopt(spec);

# try to install packages if missing
if(!is.element(opt$orgDB, installed.packages()[,1])) {
	source("http://bioconductor.org/biocLite.R")
	biocLite(opt$orgDB)
}
if(!is.element("KEGGREST", installed.packages()[,1])) {
	source("http://bioconductor.org/btryCatchiocLite.R")
	biocLite("KEGGREST")
}

# test if KEGG name is avail
testLibResult <- tryCatch({
    keggInfo(opt$keggDBName)
}, warning = function(w) {}, error = function(e) {
	return(0)
}, finally = {})
if(is.numeric(testLibResult)) {
	print(paste("[WARNING] Organism", opt$keggDBName, "not supported by KEGGREST"))
	opt$keggDBName <- NULL
	opt$plotKegg <- FALSE
}

# create output folder
if(!dir.exists(dirname(opt$output))) {
	dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
}

backgroundNames <- as.character(read.csv(opt$backgroundFile, sep="\t", head=T)[, 1])
library(opt$orgDB, character.only = TRUE)
orgDB <- opt$orgDB
keggDBName <- opt$keggDBName
pValueCutoff <- opt$pValueCutoff
output <- paste(opt$output, ".csv", sep="")

# read in all test names
genes.entrez <- list()
files <- unlist(strsplit(opt$testFiles, ";"))
names <- c()
for(i in seq(1:length(files))) {
	testNames <- as.character(read.csv(files[i], sep="\t", head=T)[, 1])
	if(length(testNames) > 0) {
		testNames = gsub("\\.[0-9]+", "", testNames)
		entrez <- bitr(testNames, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgDB)
		name <- gsub("(/.+)*/", "", dirname(files[i]))
		genes.entrez[[i]] <- entrez$ENTREZID	
		names <- c(names, name)
	}
}
# no file contained any gene names
if(length(names) == 0) {
	print("no test file contained any lines")
	endNice(opt$confirmRun2EndFile)
}

names(genes.entrez) <- names

# GO/KEGG enrichment analysis
backgroundNames = gsub("\\.[0-9]+", "", backgroundNames)
background.entrez <- bitr(backgroundNames, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgDB)
background.entrez <- background.entrez$ENTREZID

# init the stuff
compGOMF <- NULL
compGOBP <- NULL
compGOCC <- NULL
compKEGG <- NULL

# make tests
if(!is.null(opt$keggDBName)) {
	tryCatch({compKEGG <- compareCluster(geneCluster = genes.entrez, fun = "enrichKEGG", universe = background.entrez, pvalueCutoff = pValueCutoff, pAdjustMethod = "BH", organism=keggDBName)}, error = function(e) {})
}
tryCatch({compGOMF <- compareCluster(geneCluster = genes.entrez, fun = "enrichGO", universe = background.entrez, pvalueCutoff = pValueCutoff, pAdjustMethod = "BH", OrgDb = orgDB, ont="MF")}, error = function(e) {})
tryCatch({compGOBP <- compareCluster(geneCluster = genes.entrez, fun = "enrichGO", universe = background.entrez, pvalueCutoff = pValueCutoff, pAdjustMethod = "BH", OrgDb = orgDB, ont="BP")}, error = function(e) {})
tryCatch({compGOCC <- compareCluster(geneCluster = genes.entrez, fun = "enrichGO", universe = background.entrez, pvalueCutoff = pValueCutoff, pAdjustMethod = "BH", OrgDb = orgDB, ont="CC")}, error = function(e) {})

if(!is.null(compGOMF) || !is.null(compGOBP) || !is.null(compGOCC) || !is.null(compKEGG)) {
	pdff <- gsub(".csv", ".enrich.pdf", output)
	pdf(file=pdff)
	if(!is.null(compGOMF)) {
		print(dotplot(compGOMF, showCategory = 40, title = "GO Enrichment Analysis: MF", font.size=7))
	}
	if(!is.null(compGOBP)) {
		print(dotplot(compGOBP, showCategory = 40, title = "GO Enrichment Analysis: BP", font.size=7))
	}
	if(!is.null(compGOCC)) {
		print(dotplot(compGOCC, showCategory = 40, title = "GO Enrichment Analysis: CC", font.size=7))
	}
	if(!is.null(compKEGG)) {
		print(dotplot(compKEGG, showCategory = 40, title = "KEGG Pathway Enrichment Analysis", font.size=7))
	}
	dev.off()

	# write the stuff in files
	if(!is.null(compGOMF)) { write.table(compGOMF@compareClusterResult, gsub(".csv", ".GOMF.csv", output), append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE) }
	if(!is.null(compGOBP)) { write.table(compGOBP@compareClusterResult, gsub(".csv", ".GOBP.csv", output), append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE) }
	if(!is.null(compGOCC)) { write.table(compGOCC@compareClusterResult, gsub(".csv", ".GOCC.csv", output), append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE) }
	if(!is.null(compKEGG)) { write.table(compKEGG@compareClusterResult, gsub(".csv", ".KEGG.csv", output), append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE) }
}

# plot kegg stuff
if(!is.null(opt$keggDBName) && !is.null(opt$plotKegg) && !is.null(compKEGG)) {
	if(!dir.exists("/tmp/KEGG/")) {
		dir.create("/tmp/KEGG/", recursive = TRUE, showWarnings = FALSE)
	}
	for(i in seq(1:length(files))) {
		tryCatch({
			name <- gsub("(/.+)*/", "", dirname(files[i]))
			selectIDs <- compKEGG@compareClusterResult[compKEGG@compareClusterResult$Cluster == name, c("ID", "geneID")]

			expression <- read.csv(gsub("significant-", "significant", gsub("[0-9]+-fold", "significant", files[i])), sep="\t", head=T)[c("ID", opt$foldchangeCol)]
			colnames(expression) <- c("ID", opt$foldchangeCol)
			expression$ID <- gsub("\\.[0-9]+", "", expression$ID)
			tmp <- bitr(expression$ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = orgDB)
			expression <- merge(expression, tmp, by.x="ID", by.y="ENSEMBL")
			converted <- c(expression$log2FC)
	 		names(converted) <-expression$ENTREZID
			dir.create(gsub(".csv", "_KEGG/", output))
			setwd(gsub(".csv", "_KEGG/", output))
			ii <- 1
			while(ii <= nrow(selectIDs)) {
				n <- selectIDs[ii, "ID"]
				genesOfInterest <- unlist(strsplit(selectIDs[ii, "geneID"], "/"))
				expressionInterest <- converted[names(converted) %in% genesOfInterest]
				maxFC <- max(abs(c(min(expressionInterest, na.rm=T), max(expressionInterest, na.rm=T))))
				pathview(gene.data  = converted, pathway.id = n, species = keggDBName, limit = list(gene=maxFC, cpd=1), kegg.dir="/tmp/KEGG/", node.sum="median", bins=50, low = list(gene = "red"), high = list(gene = "green"), kegg.native = TRUE)
				ii<-ii+1
			}
		}, error = function(e) {})
	}
}
endNice(opt$confirmRun2EndFile)



