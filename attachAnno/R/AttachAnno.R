library(getopt)
library(stringi)

endNice <- function(file) {
	# write a file that we know, the script run to its end
	if(!is.null(file)) {
		file.create(file, showWarnings = FALSE)
	}
	quit("no")
}

fixEscape <- function(x) {
	return(stri_unescape_unicode(x))
}

# options to parse
spec <- matrix(c('targetFile', 't', 1, "character",
		  'outputFile', 'o', 1, "character",
		  'targetIDcolumn', 'i', 1, "character",
		  'annotationFile', 'a', 1, "character",
		  'annotationIDcolumn', 'c', 1, "character",
		  'annotationSep', 'y', 1, "character",
		  'targetSep', 's', 1, "character",
                  'confirmRun2EndFile', 'z', 1, "character"), ncol=4, byrow=T)

# parse the parameters
opt = getopt(spec)
# create output folder
if(!dir.exists(dirname(opt$output))) {
	dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)
}

# read target file
targetFile <- read.csv(opt$targetFile, sep=fixEscape(opt$targetSep), head=T)
# check, if ID column is there
if(length(which(colnames(targetFile) == opt$targetIDcolumn)) != 1) {
	print(paste("[ERROR] Column with ID '", opt$targetIDcolumn, "' was not found in target file '", opt$targetFile, "'.", sep=""))
	quit(save = "no", status = 14, runLast = FALSE) # status = 14 <--> invalid arguments
}

# read in all annotation files
annoFiles <- unlist(strsplit(opt$annotationFile, ","))
annoFilesSep <- unlist(strsplit(opt$annotationSep, ","))
annoFilesID <- unlist(strsplit(opt$annotationIDcolumn, ","))

# ensure that sep and id array has the correct length
if(length(annoFilesSep) == 1) {
	annoFilesSep <- rep(annoFilesSep, length(annoFiles))
}
if(length(annoFilesID) == 1) {
	annoFilesID <- rep(annoFilesID, length(annoFiles))
}

m <- targetFile
# merge all the annotation files
for(i in 1:length(annoFiles)) {
	filename <- annoFiles[i]
	sepA <- annoFilesSep[i]
	id <- annoFilesID[i]
	print(paste("reading '", filename, "'...", sep=""))
	anno <- read.csv(filename, sep=fixEscape(sepA), head=T)

	# ensure that mapping ID is there
	if(length(which(colnames(anno) == id)) != 1) {
		print(paste("[ERROR] Column with ID '", id, "' was not found in annotation file '", filename, "'.", sep=""))
		quit(save = "no", status = 14, runLast = FALSE) # status = 14 <--> invalid arguments
	}

	# merge target and annotation
	m <- merge(m, anno, by.x=opt$targetIDcolumn, by.y=id, all.x=T, all.y=F)
}

# write results
write.table(m, file=opt$outputFile, quote=F, sep=fixEscape(opt$targetSep), row.names=F)

# we are done!
endNice(opt$confirmRun2EndFile)



