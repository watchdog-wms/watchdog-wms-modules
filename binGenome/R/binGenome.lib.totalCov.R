# load required R packages
if(exists("TARGETS") && TARGETS == 1) {
	appendRequirePackages(c("rtracklayer", "cluster", "pdist", "RColorBrewer"))
} else {
	require(rtracklayer)
	require(cluster)
	require(pdist)
	require(RColorBrewer)
}
# TODO: shorter functions + more comments ?

# writes GRanges to a BED file
# ranges: instance of GRanges
# filename: path to file to store the GRanges
# IDCol: name of the column part of the meta data that is used as ID
writeGRangesToBed <- function(ranges, filename, IDCol) {
	if(!is.null(ranges$score)) {
		scores  <- ranges$score
		scores[is.na(scores)] <- "."
	}
	else {
		scores <- "."
	}
	df <- cbind(as.character(seqnames(ranges)), as.character(start(ranges)), as.character(end(ranges)), as.character(unlist(ranges@elementMetadata[IDCol])), scores, as.character(strand(ranges)))
	write.table(df, file = filename, sep="\t", quote=F, row.names=F, col.names=F)
}

# writes GRanges to a SAF file
# ranges: instance of GRanges
# filename: path to file to store the GRanges
# IDCol: name of the column part of the meta data that is used as ID
writeGRangesToSAF <- function(ranges, filename, IDCol) {
	df <- cbind(as.character(unlist(ranges@elementMetadata[IDCol])), as.character(seqnames(ranges)), as.character(start(ranges)), as.character(end(ranges)), as.character(strand(ranges)))
	colnames(df) <- c("GeneID", "Chr", "Start", "End", "Strand")
	write.table(df, file = filename, sep="\t", quote=F, row.names=F)
}

# create a GRanges instance from a BED file
# filename: path to the BED file with columns "seqnames", "start", "end", "width", "strand"
# dropMetaData: [optional] if TRUE only columns with names "seqnames", "start", "end", "width", "strand" are kept
readGRangesFromBed <- function(filename, dropMetaData = FALSE) {
	data <- read.csv(filename, sep="\t", head=T)
	meta <- NULL
	if(!dropMetaData) {
		notMeta <- c("seqnames", "start", "end", "width", "strand")
		meta <- data[, !colnames(data) %in% notMeta]

	}
	gr <- GRanges(seqnames = data$seqnames, strand = data$strand, ranges = IRanges(start = data$start, end = data$end), meta)
	return(gr)
}

# dist function
dist.eucl.lin <- function(x1, x2) { pdist(x1, x2)@dist[1] }
applyLinEuclDistOnMatrix <- function(i, d1, d2) { a <- d1[i,]; b <- d2[i,]; a <- a/sum(a); b <- b/sum(b); dist.eucl.lin(a, b) }


################################################################################
####### functions to process gene bin tables created by genomeBinner.jar #######
################################################################################
# class definition with functions: 'binMatrix'
binMatrixConstructor <- setClass("binMatrix", slots = c(cov="matrix", libsize="numeric", totalCov="numeric", binLength="vector", attributes="list"))
setGeneric("getAttribute", function(matrix, name) standardGeneric("getAttribute"))
setGeneric("normByLibsizeFun", function(matrix, ...) standardGeneric("normByLibsizeFun")) # normLibsize is optional
setGeneric("normByBinlengthFun", function(matrix) standardGeneric("normByBinlengthFun"))
setGeneric("normByShapeSumFun", function(matrix) standardGeneric("normByShapeSumFun"))
setGeneric("normByShapeMaxFun", function(matrix) standardGeneric("normByShapeMaxFun"))
setGeneric("getSummedShape", function(matrix, ...) standardGeneric("getSummedShape")) # normBylibSize, normByBinlength, normLibsize is optional
setGeneric("filterByIDs", function(data, filterIDs, ...) standardGeneric("filterByIDs")) # removeEndings is optional
setGeneric("getClusterCountKmeans", function(matrix, ...) standardGeneric("getClusterCountKmeans")) # cluster.max, iter.max is optional
setGeneric("clustKmeans", function(matrix, k, ...) standardGeneric("clustKmeans"))
setGeneric("divideByIDLists", function(matrix, filterIDLists, ...) standardGeneric("divideByIDLists")) # removeEndings is optional
setGeneric("distance", function(matrix1, matrix2, ...) standardGeneric("distance")) # FUN is optional
setGeneric("shapeCount", function(matrix) standardGeneric("shapeCount"))
setGeneric("getTopCoveredElements", function(matrix, n) standardGeneric("getTopCoveredElements"))
setGeneric("reduce2TopCoveredElements", function(data, n, ...) standardGeneric("reduce2TopCoveredElements"))

# returns the number of regions part of the bin matrix
# matrix: instance of binMatrix class
setMethod("shapeCount", signature("binMatrix"), function(matrix) {
	return(nrow(matrix@cov))
})


# gets an attribute stored together with the bin matrix or NULL if attribute does not exist
# matrix: instance of binMatrix class
# name: attribute name
setMethod("getAttribute", signature("binMatrix", "character"), function(matrix, name) {
	attr <- matrix@attributes
	index <- which(names(attr) == name)
	if(length(index) == 0) {
		print(paste("[ERROR] attribute with name '", name, "' does not exist!", sep=""))
		return(NULL)
	}
	return(attr[[name]])
})


# norms the bin matrix by the library size stored in attribute 'norm.libsize' if not done yet
# matrix: instance of binMatrix class
# normLibsize: [optional] number of reads to which the matrix should be normed
setMethod("normByLibsizeFun", signature("binMatrix"), function(matrix, ..., normLibsize = 10^9) {
	if(getAttribute(matrix, "norm.libsize") == FALSE) {
		new <- matrix
 		new@cov <- new@cov * (normLibsize / new@libsize)
		new@attributes[["norm.libsize"]] <- TRUE
		return(new)
	}
	return(matrix)
})


# norms the bin matrix by the length of the bins (each column is normed by its corresponding bin size) if not done yet
# matrix: instance of binMatrix class
setMethod("normByBinlengthFun", signature("binMatrix"), function(matrix) {
	if(getAttribute(matrix, "norm.binlength") == FALSE) {
		new <- matrix

		bLength <- lapply(new@binLength, function(b) { unlist(strsplit(b, split=',')) })
		lMatrix <- unlist(lapply(bLength, function(x) { x <- unlist(strsplit(x, split=':')); mod <- 1:length(x) %% 2 == 1; rep(x[mod], x[!mod]) }))
		lMatrix <- as.numeric(matrix(lMatrix, byrow=T, nrow=nrow(new@cov)))

	 	new@cov <- new@cov / lMatrix
		new@attributes[["norm.binlength"]] <- TRUE
		return(new)
	}
	return(matrix)
})


# norms each shape to its sum --> sum of all bins within a single shape = 1
# matrix: instance of binMatrix class
setMethod("normByShapeSumFun", signature("binMatrix"), function(matrix) {
	if(getAttribute(matrix, "norm.sumShape") == FALSE) {
		new <- matrix
		new@cov <- t(apply(matrix@cov, 1, function(x) { s <- sum(x); x/max(1, s) }))
		new@attributes[["norm.sumShape"]] <- TRUE
		return(new)
	}
	return(matrix)
})

# norms each shape to its maximal value
# matrix: instance of binMatrix class
setMethod("normByShapeMaxFun", signature("binMatrix"), function(matrix) {
	if(getAttribute(matrix, "norm.maxShape") == FALSE) {
		new <- matrix
		new@cov <- t(apply(matrix@cov, 1, function(x) { m <- max(x); x/max(1, m) }))
		new@attributes[["norm.maxShape"]] <- TRUE
		return(new)
	}
	return(matrix)
})

# aggregates the data by columns (= multiple regions) in order to get a meta-shape
# data is normalized in several steps in that order: 1) library size 2) bin length 3) intensity
# if 3) is enabled and aggregate is 'mean': sum of bins of resulting meta-shapes is 1
# matrix: instance of binMatrix class
# normLibsize: [optional] number of reads to which the matrix should be normed
# normBylibSize: [optional] enables library size normalization
# normByBinlength: [optional] enables bin length normalization
# normByShapeSum: [optional] enables shape intensity normalization by sum of shape
# aggregateFun: [optional] function that is used to aggregate data of one bin
# normByShapeMax: [optional] enables shape intensity normalization by max of shape
setMethod("getSummedShape", signature("binMatrix"), function(matrix, ..., normBylibSize = FALSE, normByBinlength = FALSE, normLibsize = 10^9, normByShapeSum = FALSE, normByFeatureCount = NULL, aggregateFun = mean, normByShapeMax = FALSE) {
	if(!is.null(normByFeatureCount)) {
		print(paste("[ERROR] Attribute 'normByFeatureCount' of function 'getSummedShape' is deprecated!", sep=""))
		quit("no", 1)
	}
	new <- matrix

	if(normBylibSize == TRUE) {
		new <- normByLibsizeFun(matrix, normLibsize = normLibsize)
	}
	if(normByBinlength == TRUE) {
		new <- normByBinlengthFun(new)
	}
	if(normByShapeMax == TRUE) {
		new <- normByShapeMaxFun(new)
	}
	if(normByShapeSum == TRUE) {
		new <- normByShapeSumFun(new)
	}
	cov <- new@cov
	cs <- apply(new@cov, 2, FUN = aggregateFun)

	new@attributes[["meta"]] <- cs
	new@attributes[["cov"]] <- cov
	return(new)
})

# filter the entries of a binMatrix by a vector of IDs
# data: instance of binMatrix class
# filterIDs: IDs that should be part of the resulting binMatrix instance
# removeEndings: [optional] removes Gencode version numbers from names which are part of the suffix (.[0-9]+$)
setMethod("filterByIDs", signature("binMatrix", "vector"), function(data, filterIDs, ..., removeEndings = FALSE) {
	new <- data
	namesToTest <- rownames(new@cov)
	if(removeEndings) {
		namesToTest <- gsub("\\.[0-9]+$", "", namesToTest)
	}
	keep <- namesToTest %in% filterIDs
	res <- data@cov[keep, ]
	if(!is.matrix(res)) {
		res <- t(as.matrix(res, byrow=F, ncol=length(which(keep))))
		rownames(res) <- rownames(data@cov)[keep]
	}
	new@cov <- res
	new@binLength <- data@binLength[keep]
	new@totalCov <- data@totalCov[keep]
	return(new)
})

# filter the entries of a dataframe by a vector of IDs
# data: instance of data.frame with rownames are IDs that are filtered
# filterIDs: IDs that should be part of the resulting data.frame instance
# removeEndings: [optional] removes Gencode version numbers from names which are part of the suffix (.[0-9]+$)
setMethod("filterByIDs", signature("data.frame", "vector"), function(data, filterIDs, ..., removeEndings = FALSE) {
	namesToTest <- rownames(data)
	if(removeEndings) {
		namesToTest <- gsub("\\.[0-9]+$", "", namesToTest)
	}
	keep <- namesToTest %in% filterIDs
	new <- data[keep, ]
	return(new)
})


# find optimal number of cluster for k-means clustering using clusGap() function
# matrix: instance of binMatrix class
# cluster.max: [optional] maximal number of clusters
# iter.max: [optional] maximal number of iterations
setMethod("getClusterCountKmeans", signature("binMatrix"), function(matrix, ..., cluster.max = 10, iter.max = 500) {
	stat <- clusGap(matrix@cov, FUN = kmeans, nstart = 25, K.max = cluster.max, B = iter.max)
	gap <- as.data.frame(stat$Tab)$gap
	m <- max(gap)
	return(which(gap == m)[1])
})

# cluster shapes using k-means
# matrix: instance of binMatrix class
# k: nuber of clusters 
# iter.max: [optional] maximal number of iterations until convergence
setMethod("clustKmeans", signature("binMatrix", "numeric"), function(matrix, k, ..., iter.max = 1000) {
	# ensure that it is normalized as otherwise it does not make much sense
	new <- normByBinlengthFun(normByLibsizeFun(matrix))
	# cluster it using k-means
	clusters <- kmeans(new@cov, k, iter.max = iter.max)
	ids <- lapply(1:k, function(i, c) { names(c$cluster[c$cluster == i]) }, clusters)
	# reorder it by size
	o <- order(unlist(lapply(ids, length)), decreasing = T)
	ids <- lapply(o, function(i) { ids[[i]] })
	names(ids) <- paste("cluster", 1:k, sep="")
	return(ids)
})

# divides a binMatrix object into multiple objects based on ID lists (e.g. created by clustKmeans)
# matrix: instance of binMatrix class
# filterIDLists: list of lists with IDs that should be part of the resulting (sub)-binMatrix instances
# removeEndings: [optional] removes Gencode version numbers from names which are part of the suffix (.[0-9]+$)
setMethod("divideByIDLists", signature("binMatrix", "list"), function(matrix, filterIDLists, ..., removeEndings = FALSE) {
	l <- lapply(filterIDLists, function(filterIDs) { filterByIDs(matrix, filterIDs, removeEndings = removeEndings) })
	names(l) <- paste(names(l), " (n=", lapply(l, shapeCount), ")", sep="")
	return(l)
})

# divides a data.frame object into multiple df based on ID lists (e.g. created by clustKmeans)
# matrix: instance of data.frame with rownames are IDs that are filtered
# filterIDLists: list of lists with IDs that should be part of the resulting (sub)-binMatrix instances
# removeEndings: [optional] removes Gencode version numbers from names which are part of the suffix (.[0-9]+$)
setMethod("divideByIDLists", signature("data.frame", "list"), function(matrix, filterIDLists, ..., removeEndings = FALSE) {
	l <- lapply(filterIDLists, function(filterIDs) { filterByIDs(matrix, filterIDs, removeEndings = removeEndings) })
	names(l) <- paste(names(l), " (n=", lapply(l, shapeCount), ")", sep="")
	return(l)
})

# calculates pairwise distances for all equally named entries in two binMatrix objects
# matrix1: instance of binMatrix class
# matrix2: instance of binMatrix class
# FUN: [optional] function to calculate distances (default: norm shapes to 1; calculate euclidian distance between resulting vectors)
setMethod("distance", signature("binMatrix", "binMatrix"), function(matrix1, matrix2, ..., FUN = applyLinEuclDistOnMatrix) {
	targetIds <- intersect(rownames(matrix1@cov), rownames(matrix2@cov))
	matrix1 <- filterByIDs(matrix1, targetIds)@cov
	matrix2 <- filterByIDs(matrix2, targetIds)@cov
	matrix1 <- matrix1[targetIds, ]
	matrix2 <- matrix2[targetIds, ]
	rho <- unlist(lapply(1:length(targetIds), FUN = FUN, matrix1, matrix2))
})

# aggregates more than one binMatrix instances to one using the FUN function
# ensures, that samples are normed by lib size and bin length before function is applied
# matrix: list of instance of binMatrix class
# FUN: [optional] function to aggregate samples (default: median)
# normLibsize: [optional] number of reads to which the matrix should be normed (default: 10^9)
groupMatrices <- function(binMatrixList, FUN = median, normLibsize = 10^9) {
	# ensure that all matrices are normed
	normedList <- list()
	i <- 1
	for(mat in binMatrixList) {
		if(i == 1) {
			totalCov <- mat@totalCov
		} else {
			totalCov <- totalCov+mat@totalCov
		}
		new <- normByLibsizeFun(mat, normLibsize = normLibsize)
		new <- normByBinlengthFun(new)
		normedList[[i]] <- new@cov
		i <- i + 1
	}
	# group values
	covData <- apply(simplify2array(normedList), c(1,2), FUN)
	return(createGroupedMatrix(covData, binMatrixList[[1]]@binLength, totalCov))
}

# creates a object of class 'binMatrix' from a file written by genomeBinner.jar
# file:
# sumFile:
# removeEndings:
readBinTable <- function(file, sumFile, removeEndings = FALSE) {
	size <- read.csv(sumFile, sep="\t", head=T)[1,1]
	m <- read.csv(file, sep="\t", head=T)

	# remove the endings
	if(removeEndings) m[, 1] <- gsub("\\.[0-9]+$", "", m[, 1])

	# make first col be the names
	rownames(m) <- m[, 1]
	m <- m[, -1]
	bLength <- m[, ncol(m)]
	m <- m[, -ncol(m)]
	m <- as.matrix(m)
	total <- apply(m, 1, sum)

	# create the attributes
	a <- list("file" = file, "norm.libsize" = FALSE, "norm.binlength" = FALSE, "norm.sumShape" = FALSE, "merged" = FALSE, "norm.maxShape" = FALSE, "ratio" = FALSE)
	
	obj <- binMatrixConstructor(cov = m, libsize = size, binLength = bLength, totalCov = total, attributes = a)
	return(obj)
}

# function to create a aggregated 'binMatrix' instance
# cov: matrix containing the coverage of the shapes (and IDs as colnames)
# bLength: vector containing the length of the bins
createGroupedMatrix <- function(cov, bLength, total) {
	# create the attributes
	a <- list(file = NULL, "norm.libsize" = TRUE, "norm.binlength" = TRUE, "norm.sumShape" = FALSE, "merged" = TRUE, "norm.maxShape" = FALSE, "ratio" = FALSE)

	obj <- binMatrixConstructor(cov = cov, libsize = 0, totalCov = total, binLength = bLength, attributes = a)
	return(obj)
}

# get the IDs of the top n covered elements
# matrix: instance of binMatrix class
# n: number of elements to return
setMethod("getTopCoveredElements", signature("binMatrix"), function(matrix, n) {
	ids <- rownames(matrix@cov)
	if(length(ids) < n) n <- length(ids)
	cov <- matrix@totalCov
	or <- order(cov, decreasing = T)[1:n]
	return(ids[or])
})

# filters the matrix based on the top n covered elements
# matrix: instance of binMatrix class
# n: number of elements to return
# ...: is passed to filterByIDs
setMethod("reduce2TopCoveredElements", signature("binMatrix"), function(data, n, ...) {
	ids <- getTopCoveredElements(data, n)
	return(filterByIDs(data, ids, ...))
})

modifyTitleForParse <- function(title, addSpace = FALSE) {
	if(length(title) > 1) {
		unlist(lapply(title, FUN = modifyTitleForParse, addSpace = addSpace))
	} else {
		# replace >= and <= in title 
		modTitle <- title
		modTitle <- gsub(">=", "'*''>=''*'" , modTitle) 
		modTitle <- gsub("<=", "'*''<=''*'" , modTitle)

		# replace some other values [A-Za-z%] that are enclosed in *...* with '*...*' for parse
		if(grepl("\\*([A-Za-z%]+?)\\*", modTitle, perl=T) && !grepl("'\\*([A-Za-z%]+?)\\*'", modTitle, perl=T)) {
			modTitle <- gsub("\\*([A-Za-z%]+?)\\*","'*\\1*'", modTitle, perl=T)
		}
		# ensure that title string starts and ends with "'"
		if(!grepl("^'.*'$", modTitle)) {
			modTitle <- paste("'", modTitle, "'", sep="")
		}
		modTitle <- gsub("(\\*%)|(%\\*)", "%" , modTitle)

		# add a spacer
		if(addSpace) {
			modTitle <- gsub("^'", "'     ", modTitle)
		}
		return(modTitle)
	}
}

# function that transforms a list of pvalues into a list of colors
# pvs: pvalues
# defaultColor: [optional] default color if pvalue is not smaller or equal than any threshold
# pvThresholds: [optional] thresholds that are required to assign a specific color (p <= t)
# colors: [optional] colors associated with the thresholds; must have the same length as pvThresholds
defaultPvalueColorTransformer <- function(pvs, defaultColor = "lightgrey", pvThresholds = c(10^-3, 10^-10, 10^-15), colors = c("yellow", "orange", "red")) {
	if(length(pvThresholds) != length(colors)) {
		print(paste("[ERROR] Number of pValue thresholds must be equal to number of colors. (", length(pvThresholds), " vs. ", length(colors), ")", sep="")) 
		return(-1)
	}
	# add default value and name the values 
	pvThresholds <- c(pvThresholds, Inf)
	colors <- c(colors, defaultColor)
	names(pvThresholds) <- colors
	
	# get the smallest thresholds that is bigger than the pvalue 
	return(unlist(lapply(pvs, function(x, t) { names(t[min(t[x <= t]) == t]) }, pvThresholds)))
}

# creates a shape plot for a list of binMatrix instances (all shapes in one plot window)
# binMatrixList: list of binMatrix instances for which the shapes should be plotted
# labels: body labels of (scaled) X axis; assumes that labels are equally spaced over all bins and that first and last position is labeled;
# title: [optional] title of the shape panel; default: ""
# name: [optional] name added as prefix to the auto-generated Y-label; default: ""
# smooth: [optional] use function to smooth the shapes (e.g. smooth()); default: NULL
# showYLab: [optional] display Y label at the left side; default: TRUE
# palette: [optional] colors used for the shapes in the same order as in binMatrixList; default: NULL
# normBylibSize: [optional]; see getSummedShape(...); default: TRUE
# normByBinlength: [optional]; see getSummedShape(...); default: TRUE
# normLibsize: [optional]; see getSummedShape(...); default: TRUE
# normByShapeSum: [optional]; see getSummedShape(...); default: TRUE
# aggregateFun: [optional]; see getSummedShape(...); default: mean
# verticalLinesOnTicks: [optional] draw vertical lines at these bin positions; vector with bin indices (1:N); default: NULL
# fixedLabelsStart: [optional] X labels for fixed-length start bins; see labels; default: NULL
# fixedLabelsEnd: [optional] X labels for fixed-length end bins; see labels; default: NULL
# fixedLabelsStartTotalBins: [optional] number of bins defined as fixed-length start bins; default: 0
# fixedLabelsEndTotalBins: [optional] number of bins defined as fixed-length end bins; default: 0
# legendNCol: [optional] number of columns in legend; default: 1
# ylimScale: [optional] maximal y-value of all shapes is multiplied with that factor (e.g. allows zoom); default: 1
# ylab: [optional] custom y-axis label; if not set it is determined based on the applied normalization; default: NULL
# showXLab: [optional] show or hide labels on x-axis tick marks; default: TRUE
# showLegend: [optional] show or hide legend; default: TRUE
# minYA: [optional] custom minimal y value; only in combination with maxYA; default: NULL
# maxYA: [optional] custom maximal y value; only in combination with minYA; default: NULL
# legendPos: [optional] position of the legend; parameter of legend(...); e.g. top, topright; default: "top"
# performWilcoxTest: [optional] perform position-wise, paired wilcox-test and indicates corrected p-values with colors if two matrices are given; default: FALSE
# wilcoxTestPVTransformFUN: [optional] function that transforms a list of pvalues into a list of colors 
# cexAxis: [optional] cex factor for axis; default: 1.3
# cexLab: [optional] cex factor for axis labels; default: 1.2
# cexTitle: [optional] cex factor for title; default: 1.6
# shapeLineWidth: [optional] line width for the shapes; default: 1.75
# ablineLineType: [optional] line type of vertical lines; default: 3
# cexLegend: [optional] cex factor for legend labels; default: 1.5
# legendInset: [optional] inset distance(s) from the margins as a fraction of the plot region when legend is placed by keyword; parameter of legend(...); default: .006
# legendXIntersp: [optional] character interspacing factor for horizontal spacing of legend; parameter of legend(...); default: 0.6
# legendYIntersp: [optional] character interspacing factor for vertical spacing of legend; parameter of legend(...); default: 0.6
# cexLegendSymbol: [optional] cex factor for legend symbols; default: 2
# legendSpacingBetween: : [optional] horizontal spacing between two legend items; default: 25
# lineYlab: [optional] controls horizontal spacing of ylab; default: 3
# lineTitle: [optional] controls vertical spacing of title; default: 1
# numberOfYTicks: [optional] number of ticks on the y-axis; automatically determined by default; default: NULL
# asIsYlab: [optional] if set to TRUE, ylab is taken as it is without modification; default: FALSE
# normByShapeMax: [optional]; see getSummedShape(...); default: FALSE
plotGroup <- function(binMatrixList, labels, title = "", name = "", smooth = NULL, showYLab = TRUE, palette = NULL, normBylibSize = TRUE, normByBinlength = TRUE, normLibsize = 10^9, normByShapeSum = FALSE, aggregateFun = mean, verticalLinesOnTicks = NULL, fixedLabelsStart = NULL, fixedLabelsEnd = NULL, fixedLabelsStartTotalBins = 0, fixedLabelsEndTotalBins = 0, legendNCol = 1, ylimScale = 1, ylab = NULL, showXLab = TRUE, showLegend = TRUE, minYA = NULL, maxYA = NULL, legendPos = "top", performWilcoxTest = FALSE, wilcoxTestPVTransformFUN = defaultPvalueColorTransformer, cexAxis=1.3, cexLab=1.2, cexTitle=1.6, shapeLineWidth=1.75, ablineLineType=3, cexLegend=1.5, legendInset=.006, legendXIntersp=0.6, legendYIntersp = 0.6, cexLegendSymbol = 2, legendSpacingBetween = 25, lineYlab = 3, lineTitle = 1, numberOfYTicks = NULL, asIsYlab = FALSE, normByShapeMax = FALSE) {
	n <- length(binMatrixList) # get number of shapes to plot

	# get automatic y-label if none is set
	if(!asIsYlab) {
		if(is.null(ylab)) {
			if(showYLab) {
				normByFeatureCount <- identical(aggregateFun, mean)
				ylab <- "coverage per bin"
				if(normByBinlength) ylab <- "avg. coverage per bp"
				if(normBylibSize) ylab <- paste(ylab, " / library size * ", normLibsize, sep="")
				if(normByShapeMax) ylab <- "coverage normed by shape max"
				if(normByShapeSum) ylab <- "coverage normed by shape sum"
				if(normByFeatureCount) ylab <- paste(ylab, " / # shapes", sep="")
				if(normByShapeSum && normByFeatureCount) ylab <- "occupancy"
				if(name != "") ylab <- paste(name, ylab, sep=" ")
				if(!normByFeatureCount) {
					print(paste("[WARN] Custom aggregate function ('aggregateFun') for bins is used in plotGroup() but no custom 'ylab' label is set!", sep=""))
				}
			} else ylab <- ""
		}
	}
	# get default colors if none are set
	if(is.null(palette)) {
		nA <- max(3, n)
		set <- "Set1"
		if(nA > 9) {
			set <- "Set3"
		}
		if(nA > 12) {
			print("Not more than 12 colors can be set automatically!")
			return(FALSE)
		}
		palette <- brewer.pal(nA, set)
		# switch colors 1 and 2
		tmpc <- palette[1]
		palette[1] <- palette[2]
		palette[2] <- tmpc
	}

	# find max Y value
	maxY <- -10^16
	minY <- -maxY
	metaStore <- list()
	covStore <- list()
	# loop over all matrices
	for(i in 1:n) {
		binMatrix <- binMatrixList[[i]]
		# calculate the shape
		meta <- getSummedShape(binMatrix, normBylibSize = normBylibSize, normByBinlength = normByBinlength, normLibsize = normLibsize, normByShapeSum = normByShapeSum, aggregateFun = aggregateFun, normByShapeMax = normByShapeMax)

		# store the coverage and meta data of that shape
		covStore[[i]] <- meta@attributes[["cov"]]
		meta <- meta@attributes[["meta"]]
		metaStore[[i]] <- meta

		# get max and min Y value
		minY <- min(minY, meta)
		maxY <- max(maxY, meta)
	}
	# scale maximal Y value
	maxY <- maxY * ylimScale
	# overwrite min and max Y with custom value
	if(!is.null(minYA) && !is.null(maxYA)) {
		minY <- minYA
		maxY <- maxYA
	}
	# calulcate Y range
	YRange <- maxY-minY
	# set some range, if none is given (sample normed with itself)
	if(maxY == 1 && minY == 1) {
		minY <- 1 - 0.001
		maxY <- 1 + 0.001
	}
	# set minimal YRange to display
	YRange <- max(YRange, 0.001)

	# get some space for p-values when a statistic test if performed
	pvs <- NULL
	if(performWilcoxTest) {
		minY <- minY - 0.07*(YRange)
		pvPos <- minY + 0.005*(YRange)

		# perform position-wise, paired wilcox-test if two matrices are given
		if(n == 2) {
			a <- covStore[[1]]
			b <- covStore[[2]]
			# TODO: replace with better test 
			if(nrow(a) >= 10) {
				pvs <- c()
				for(p in 1:ncol(a)) {
					pv <- wilcox.test(a[,p], b[,p], exact=F, paired=T)$p.val
					pvs <- c(pvs, pv)
				}
				# correct p-values using bonferroni to be conservative
				pvs <- p.adjust(pvs, method = "bonferroni", n=length(pvs))
				pvs[is.nan(pvs)] <- 1
				# colors indicate result of tests
				pvsCol <- wilcoxTestPVTransformFUN(pvs)
			}
		}
	}

	# plot the shape lines
	range <- NULL
	for(i in 1:n) {
		meta <- metaStore[[i]]
		# smooth the shapes if required
		if(!is.null(smooth)) meta <- as.numeric(do.call(smooth, list(meta)))

		# add plot window, title and more...
		if(i == 1) {
			N <- length(meta) # number of bins
			range <- (1:N) - round(N/2) # get pseudo x values

			# prepare title for parse
			modTitle <- modifyTitleForParse(title, addSpace = TRUE)

			# create a new but empty plot window with the correct ranges
			plot(-1, -1, type="l", xaxt="n", yaxt="n", ylim=c(minY, maxY), xlim=c(min(range), max(range)), xlab="", ylab="", main=NULL)
			title(parse(text=modTitle), line=lineTitle, cex.main=cexTitle)
			# add Y label if labels should be dispayed
  		if(showYLab) {
				mtext(side = 2, text=ylab, cex=cexLab, line=lineYlab)
			}
		}

		# plot the actual shape
		lines(range, meta, type="l", col=palette[i], lwd=shapeLineWidth, xpd = FALSE)
	}

	# create x-axis ticks
	# create main interval
	labelsRange <- (N-fixedLabelsStartTotalBins-fixedLabelsEndTotalBins)
	if(length(labels >= 2)) {
		interval <- round(seq(1+fixedLabelsStartTotalBins, N+1-fixedLabelsEndTotalBins, by=(labelsRange/(length(labels)-1))))
	}
	else {
		interval <- round(seq(1+fixedLabelsStartTotalBins, N+1-fixedLabelsEndTotalBins, by = labelsRange/15))
		labels <- rep("", length(interval))
	}
	# add fixed bin ticks
	if(!is.null(fixedLabelsStartTotalBins) && length(fixedLabelsStart) > 0) {
		intervalStart <- round(seq(1, fixedLabelsStartTotalBins+1, by=fixedLabelsStartTotalBins/(length(fixedLabelsStart)-1)))
		interval <- c(intervalStart, interval[-1])
		labels <- c(fixedLabelsStart, labels[-1])
	}
	if(!is.null(fixedLabelsEndTotalBins) && length(fixedLabelsEnd) > 0) {
		intervalEnd <- round(seq(1, fixedLabelsEndTotalBins+1, by=fixedLabelsEndTotalBins/(length(fixedLabelsEnd)-1)))
		intervalEnd <- intervalEnd + fixedLabelsStartTotalBins + labelsRange
		interval <- c(interval[-length(interval)], intervalEnd)
		labels <- c(labels[-length(labels)], fixedLabelsEnd)
	}

	# axis for ticks
	interval[1] <- min(1, interval[1])
	interval[length(interval)] <- N
	interval <- interval - round(N/2)

	if(showXLab) {
		axis(1, at=interval, labels=labels, las=2, cex.axis=cexAxis)
	}
	else {
		axis(1, at=interval, labels=rep("", length(labels)), las=2, cex.axis=cexAxis)
	}
	if(!is.null(numberOfYTicks)) {
		if(length(numberOfYTicks) == 1) {
			axis(2, cex.axis=cexAxis, at = seq(signif(minY, 2), signif(maxY, 2), length.out=numberOfYTicks))
		} else {
			axis(2, cex.axis=cexAxis, at = numberOfYTicks)
		}
	} else {
		axis(2, cex.axis=cexAxis)
	}

	# plot ablines if wished
	if(!is.null(verticalLinesOnTicks)) {
		for(i in verticalLinesOnTicks) {
			abline(v=interval[i], lty=ablineLineType, xpd=FALSE)
		}
	}

	# add a legend
	if(showLegend) {
		legend(x=legendPos, parse(text=modifyTitleForParse(names(binMatrixList))), pt.cex = cexLegendSymbol, pt.bg=palette, ncol=legendNCol, bty='n', cex=cexLegend, inset = legendInset, x.intersp = legendXIntersp, y.intersp = legendYIntersp, pch=22, text.width=legendSpacingBetween)
	}
	# plot p-values if test was performed
	if(!is.null(pvs)) {
		points(range, rep(pvPos, length(range)), col=pvsCol, lwd=2, pch=15)
	}
}

# creates a shape plot with multiple panels based on a list of lists of IDs
# all 4 normalizations are applied
# binMatrixList: list of different datasets which will be plotted into one panel e.g. conditions (instances of binMatrix)
# interestLists: list of gene lists used to define serveral groups (panels)
# name:  name for the y-axis; default: "TODO"
# metaDataNorm: [optional] list of instances of binMatrix that are used to norm binMatrixList; must have the same length as binMatrixList; default: NULL
# isListGeneList: [optional] if interestLists are gene IDs this flag enables mapping to transcript IDs; default: FALSE
# transGeneMapping: [optional] dataframe containing columns 'gene_id' and 'transcript_id'; default: NULL
# maxTranscripts: [optional] names of maximal expressed transcripts; default: NULL
# removeEndings: [optional] removes Gencode version numbers from IDs which are part of the suffix (.[0-9]+$); default: FALSE
# psCount:  [optional] pseudocount that is added when a ratio is calculated; default: 0.025
# ylimScale:  [optional] factor that is multiplied with the maximal observed Y value; default: 1
# fixedLabelsStartTotalBins: [optional] number of bins defined as fixed-length start bins; default: 0
# fixedLabelsEndTotalBins: [optional] number of bins defined as fixed-length end bins; default: 0
# legendPos: [optional] position of the legend; parameter of legend(...); e.g. top, topright; default: NULL
# bodyLabels:  [optional] body labels of (scaled) X axis; assumes that labels are equally spaced over all bins and that first and last position is labeled; default: NULL
# fixedLabelsStart: [optional] X labels for fixed-length start bins; see labels; default: NULL
# fixedLabelsEnd: [optional] X labels for fixed-length end bins; see labels; default: NULL
# smoothingFun: [optional] use function to smooth the shapes (e.g. smooth()); default: NULL
# verticalLinesOnTicks: [optional] draw vertical lines at these bin positions; vector with bin indices (1:N); default: NULL
# showLegend: [optional] vector of the same length as interestLists or 1 value; TRUE if legend should be displayed in the panel; default: NULL
# showXLab: [optional] vector of the same length as interestLists or 1 value; TRUE if x-axis labels should be displayed in the panel; default: NULL
# showYLab: [optional] vector of the same length as interestLists or 1 value; TRUE if y-axis labels should be displayed in the panel; default: NULL
# mar: [optional] inner margin of layout; vector with 4 values (bottom, left, top, right)
# oma: [optional] outer margin of layout; vector with 4 values (bottom, left, top, right)
# numberOfYTicks: [optional] number of ticks on the y-axis; or positions of ticks; automatically determined by default; default: NULL
# aggregateFun: [optional] function that is used to aggregate data of one bin; default: mean
# normBylibSize: [optional]; see getSummedShape(...); default: TRUE
# normByBinlength: [optional]; see getSummedShape(...); default: TRUE
# normByShapeSum: [optional]; see getSummedShape(...); default: TRUE
# normByShapeMax: [optional]; see getSummedShape(...); default: FALSE
# normName: [optional]; name of the norm data for the y-axis
# applyAggregatedNorm: [optional]; aggregate the meta shape first and then calculate the ratio instead of doing it for all regions separately; default: FALSE
# ... see parameters of plotGroup()
plotMergedShape <- function(binMatrixList, interestLists, name = "TODO", metaDataNorm = NULL, isListGeneList = FALSE, transGeneMapping = NULL, maxTranscripts = NULL, removeEndings = FALSE, psCount = 0.025, ylimScale = 1, fixedLabelsEndTotalBins = 0, fixedLabelsStartTotalBins = 0, legendPos = "top", bodyLabels = NULL, fixedLabelsStart = NULL, fixedLabelsEnd = NULL, smoothingFun = NULL, verticalLinesOnTicks = NULL, showYLab = NULL, palette = NULL, legendNCol = 1, showXLab = NULL, showLegend = NULL, minYA = NULL, maxYA = NULL, performWilcoxTest = FALSE, wilcoxTestPVTransformFUN = defaultPvalueColorTransformer, cexAxis=1.3, cexLab=1.2, cexTitle=1.6, shapeLineWidth=1.75, ablineLineType=3, cexLegend=1.5, legendInset=.006, legendXIntersp=0.6, legendYIntersp = 0.6, cexLegendSymbol = 2, legendSpacingBetween = 25, mar = c(0.5,4.5,2.5,0.5), oma=c(5.25,0.25,0.25,0.25), lineYlab = 3, lineTitle = 1, numberOfYTicks = NULL, asIsYlab = FALSE, aggregateFun = mean, normBylibSize = TRUE, normByBinlength = TRUE, normByShapeSum = TRUE, normByShapeMax = FALSE, normName = NULL, applyAggregatedNorm = FALSE) {
	minY <- 10^16
	maxY <- -minY
	nameOrg <- name

	# convert gene list to transcript list
	if(isListGeneList) {
		if(is.null(transGeneMapping) || is.null(maxTranscripts)) {
			print("[ERROR] binGenome.lib.R: transGeneMapping and maxTranscripts can not be null if it is a gene list!") 
			return(-1)
		}

		# get maximal expressed transcripts
		interestLists <- lapply(interestLists, function(x, tgm, maxT, rme) {
			tgm <- as.data.frame(tgm)
			if(rme) { tgm$gene_id <- gsub("\\.[0-9]+$", "", tgm$gene_id) }
			r <- tgm[tgm$gene_id %in% x, "transcript_id"]
			return(as.character(r[r %in% maxT]))
			}, tgm = transGeneMapping, maxT = maxTranscripts, rme = removeEndings)
	}

	for(i in 1:length(interestLists)) {
		int <- interestLists[[i]]
		if(length(int) == 0) {
			print(paste("[ERROR] binGenome.lib.R: list ", i, " does not contain any valid ID!", sep=""))
			return(-1)
		}
	}

	# norm if given 
	if(!is.null(metaDataNorm)) {
		# test if both lists have the same length
		if(length(binMatrixList) != length(metaDataNorm)) {
			print("binGenome.lib.R: binMatrixList and metaDataNorm must have the same length!") 
			return(-1)
		}
		# adjust name
	  	if(is.null(normName)) {
			# adjust name
	  		normNameVar <- names(metaDataNorm)[1]
		}
		else {
			normNameVar <- normName
			# ensure that conditions have the same names
			for(i in 1:length(binMatrixList)) {
				if(names(binMatrixList)[i] != names(metaDataNorm)[i]) {
					print(paste("binGenome.lib.R: names of data and norm conditions vary: ", names(binMatrixList)[i], " vs. ", names(metaDataNorm)[i], sep="")) 
					return(-1)
				}
			}
		}
		name <- paste(name, " / ", normNameVar, " ratio", sep="")
		normByShapeSum <- FALSE
		normByShapeMax <- FALSE
		print(paste("[WARN] Parameter 'normByShapeSum' and 'normByShapeMax' was set to FALSE as norm matrix is used in plotMergedShape().", sep=""))
		ylab <- name
	}
	else {
		ylab <- NULL
	}
	if(asIsYlab) { ylab = name }

	# divide lists 
	metaDataDivided <- lapply(binMatrixList, function(x, interestLists, removeEndings) { divideByIDLists(x, interestLists, removeEndings = removeEndings)}, interestLists, removeEndings)

	n <- length(interestLists) # number of panels
	N <- length(metaDataDivided) # number of shapes within panel
	data2PlotStore <- list()
	# find Y-scale factor per panel
	for(i in seq(1, n)) {
		data2Plot <- lapply(1:N, function(x, i, data) { data[[x]][[i]] }, i, metaDataDivided)
		names(data2Plot) <- names(binMatrixList)
		data2debug <- list()
		print(name)

		# loop over all shapes
		L <- length(data2Plot)
		for(ii in seq(1, L)) {
			# norm if given
			if(!is.null(metaDataNorm)) {
				# get IDs part of both
				data2norm <- data2Plot[[ii]]
				norm <- metaDataNorm[[ii]]
				targetIds <- intersect(rownames(data2norm@cov), rownames(norm@cov))

				# restrict to genes part of both
				data2norm <- filterByIDs(data2norm, targetIds)
				norm <- filterByIDs(norm, targetIds)

				if(!applyAggregatedNorm) {				
					# normalize each shape to the sume of 1 
					data2norm <- normByShapeSumFun(data2norm)
					norm <- normByShapeSumFun(norm)

					# calculate ratio with pseudo count
					nenner <- norm@cov+psCount
					zaehler <- data2norm@cov+psCount

					# ignore positions with nenner == 0 (use NA instead)
					calcPos <- nenner > 0
					ratio <- matrix(data = NA, nrow = nrow(nenner), ncol = ncol(nenner))
					
					ratio[calcPos] <- zaehler[calcPos] / nenner[calcPos]
					data2Plot[[ii]]@cov <- ratio
				} else {
					# calculate meta shape for both matrices
					meta_data2norm <- getSummedShape(data2norm, normBylibSize = FALSE, normByBinlength = normByBinlength, normByShapeSum = TRUE, aggregateFun = aggregateFun, normByShapeMax = FALSE)
					meta_norm <- getSummedShape(norm, normBylibSize = FALSE, normByBinlength = normByBinlength, normByShapeSum = TRUE, aggregateFun = aggregateFun, normByShapeMax = FALSE)	
	
					# calculate the ratio between these shapes		
					ratio <- (meta_data2norm@attributes[["meta"]]) / (meta_norm@attributes[["meta"]])
					data2Plot[[ii]]@cov <- as.matrix(t(ratio))
				}

				# now it's a ratio plot
				data2Plot[[ii]]@attributes[["ratio"]] <- TRUE
			}
			# get the shape to plot
			meta <- getSummedShape(data2Plot[[ii]], normBylibSize = normBylibSize, normByBinlength = normByBinlength, normByShapeSum = normByShapeSum, aggregateFun = aggregateFun, normByShapeMax = normByShapeMax)
			meta <- meta@attributes[["meta"]]

			minY <- min(minY, meta, na.rm = T)
			maxY <- max(maxY, meta, na.rm = T)
		}
		data2PlotStore[[i]] <- data2Plot
	}
	
	# overwrite min and max Y with custom value
	if(!is.null(minYA) && !is.null(maxYA)) {
		minY <- minYA
		maxY <- maxYA
	}

	if(is.null(showLegend)) { showLegend <- rep(FALSE, n) } else { if(length(showLegend) == 1) { showLegend <- rep(showLegend, n) }}
	if(is.null(showXLab)) { showXLab <- rep(FALSE, n) } else { if(length(showXLab) == 1) { showXLab <- rep(showXLab, n) }}
	if(is.null(showYLab)) { showYLab <- rep(FALSE, n) } else { if(length(showYLab) == 1) { showYLab <- rep(showYLab, n) }}

	# apply ylim scale to maxY 
	maxY <- maxY * ylimScale

	# plot it
	graphics::layout(mat = matrix(1:n, nrow = n, ncol = 1), heights = rep(c(1)/n, n), widths = 1)
	par(mar = mar, oma=oma, xpd = NA, bg=NA)
	for(i in seq(1, n)) {
		data2Plot <- data2PlotStore[[i]]
		mainName <- paste(names(interestLists)[i], " (n=", nrow(data2Plot[[1]]@cov), ")", sep="")
		showLegendPanel <- showLegend[i]
		showXLabPanel <- showXLab[i]
		showYLabPanel <- showYLab[i]

		# plot the panel
		plotGroup(data2Plot, bodyLabels, showXLab = showXLabPanel, showLegend = showLegendPanel, title = mainName, name = name, normBylibSize = normBylibSize, normByShapeSum = normByShapeSum, aggregateFun = aggregateFun, normByBinlength = normByBinlength, smooth=smoothingFun, fixedLabelsStart=fixedLabelsStart, fixedLabelsEnd=fixedLabelsEnd, fixedLabelsStartTotalBins=fixedLabelsEndTotalBins, fixedLabelsEndTotalBins=fixedLabelsStartTotalBins, verticalLinesOnTicks=verticalLinesOnTicks, legendNCol = legendNCol, showYLab = showYLabPanel, ylab = ylab, minYA = minY, maxYA = maxY, legendPos = legendPos, palette = palette, performWilcoxTest = performWilcoxTest, wilcoxTestPVTransformFUN = wilcoxTestPVTransformFUN, cexAxis=cexAxis, cexLab=cexLab, cexTitle=cexTitle, shapeLineWidth=shapeLineWidth, ablineLineType=ablineLineType, cexLegend=cexLegend, legendInset=legendInset, legendXIntersp=legendXIntersp, legendYIntersp=legendYIntersp, cexLegendSymbol = cexLegendSymbol, legendSpacingBetween = legendSpacingBetween, lineYlab = lineYlab, lineTitle = lineTitle, numberOfYTicks = numberOfYTicks, asIsYlab = asIsYlab, normByShapeMax = normByShapeMax)	
	}
}
