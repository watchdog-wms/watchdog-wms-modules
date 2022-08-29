
library(getopt)
#library(akmedoids)  #lässt sich nicht mehr installieren

warnings()

opt <- NULL
# options to parse
spec <- matrix(c('coveragefiles',   'c', 1, "character", # dir of coverage files
                 'genelist',         'g', 1, "character",   # genelist of cluster
                 'outdir',   'o', 1, "character", # dir to output
		 'exp', 't', 1, "character", #exp is atac or pro chip etc, when atac - peaks can be near together
                 'confirmRun2EndFile', 'e', 1, "character"
), ncol=4, byrow=T)

# parse the parameters
opt <- getopt(spec)


if(is.null(opt$coveragefiles)) {
  print("[ERROR] Path to coveragefiles not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$outdir)) {
  print("[ERROR] Path to output dir not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}
if(is.null(opt$genelist)) {
  print("[ERROR] Path to genelist not set")
  quit(save = "no", status = 11, runLast = FALSE) # status = 11 <--> missing arguments
}

#source("/home/proj/software/watchdog/modules/binGenome/R/binGenome.lib.R")
#source("/home/proj/software/watchdog/modules/sharedUtils/R/functions.R")
args <- commandArgs(trailingOnly = FALSE)
functionsFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../../modules/sharedUtils/R/functions.R', sep = '')
source(functionsFileName)
binFileName <- paste(dirname(sub('--file=', '', args[grep('--file=', args)])), '/../../modules/binGenome/R/binGenome.lib.R', sep = '')
source(binFileName)



all<-FALSE


covfiledir <- opt$coveragefiles
out <- opt$outdir
clustergenes <- read.csv(opt$genelist, sep="\t", header=FALSE)


clustername <- gsub(".txt", "", basename(opt$genelist))
print(paste0("cluster: ", clustername))




### reading coverage files in
binFiles <- paste(covfiledir, list.files(path=covfiledir, pattern="*101.coverage.csv"), sep="/")
mocks <- list()
wts <- list()
for(file in binFiles) {
	cov <- file
	filename <- gsub("_pro_all.101.coverage.csv", "", basename(file))
	sum <- paste0(dirname(cov), "/", filename, ".sum.csv")
	binmatrixconstruct <- readBinTable(cov, sum)
	cond <- ifelse(grepl("[mock|Mock]", filename), "mock", "WT")
	print(paste0("reading ", filename, " -> ", cond))

	if (cond == "mock") {
		entries <- length(mocks)
        	mocks[[entries+1]] <- binmatrixconstruct
	} else {
		entries <- length(wts)
                wts[[entries+1]] <- binmatrixconstruct
	}
}


mockdata <- groupMatrices(mocks, mean)
wtdata <- groupMatrices(wts, mean)
#print(rownames(mockdata@cov)) #genes
#print(colnames(mockdata@cov)) #bins


mockdata <- normByLibsizeFun(mockdata, normLibsize=10^9)
wtdata <- normByLibsizeFun(wtdata, normLibsize=10^9)

mockdata <- normByBinlengthFun(mockdata)
wtdata <- normByBinlengthFun(wtdata)

mockdata@cov <- apply(mockdata@cov, c(1,2), median)
wtdata@cov <- apply(wtdata@cov, c(1,2), median)



### get genes of cluster
mockdata@cov <- subset(mockdata@cov, rownames(mockdata@cov) %in% clustergenes[,1])
wtdata@cov <- subset(wtdata@cov, rownames(wtdata@cov) %in% clustergenes[,1])



mockmeta <- getSummedShape(mockdata,
                         normBylibSize = TRUE, normByBinlength = TRUE, normLibsize = 10^9,
                         normByShapeSum = TRUE, aggregateFun = mean, normByShapeMax = FALSE)
wtmeta <- getSummedShape(wtdata,
                         normBylibSize = TRUE, normByBinlength = TRUE, normLibsize = 10^9,
                         normByShapeSum = TRUE, aggregateFun = mean, normByShapeMax = FALSE)

mockmeta <- mockmeta@attributes[["meta"]]
wtmeta <- wtmeta@attributes[["meta"]]

num_total_bins <- 101
range <- c(1:num_total_bins)
labs <- seq(-3000, 3000, 750)
labs[5] <- "TSS"


### metagene plots for checking shape
pdf(paste0(out, "/", clustername, "-extrema.pdf"), paper="a4r", width=18)



normedtable <- as.data.frame(cbind(c(1:num_total_bins), mockmeta, wtmeta))
colnames(normedtable) <- c("bin", "mock", "WT")



### start collecting extrema
max_peak_mock <- which.max(normedtable$mock)
max_peak_wt <- which.max(normedtable$WT)

print(paste0(max_peak_mock, "  ", max_peak_wt))


x <- normedtable[,1]
y_mock <- normedtable[,2]
# Maxima
max_x_mock <- x[ggpmisc:::find_peaks(y_mock)]
max_y_mock <- y_mock[ggpmisc:::find_peaks(y_mock)]
# Minima
min_x_mock <- x[ggpmisc:::find_peaks(-y_mock)]
min_y_mock <- y_mock[ggpmisc:::find_peaks(-y_mock)]


max_mock <- as.data.frame(cbind(max_x_mock, max_y_mock))
max_mock$type <- "max"
colnames(max_mock)[1] <- "bin"
colnames(max_mock)[2] <- "mean"
min_mock <- as.data.frame(cbind(min_x_mock, min_y_mock))
min_mock$type <- "min"
colnames(min_mock)[1] <- "bin"
colnames(min_mock)[2] <- "mean"

extrema_mock <- as.data.frame(rbind(min_mock, max_mock))
extrema_mock <- extrema_mock[order(extrema_mock[,1]),]
colnames(extrema_mock) <- c("bin", "mean", "type")
write.table(extrema_mock, paste0(out, "/", clustername, "-extrema_mock.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)




y_wt <- normedtable[,3]
# Maxima
max_x_wt <- x[ggpmisc:::find_peaks(y_wt)]
max_y_wt <- y_wt[ggpmisc:::find_peaks(y_wt)]
# Minima
min_x_wt <- x[ggpmisc:::find_peaks(-y_wt)]
min_y_wt <- y_wt[ggpmisc:::find_peaks(-y_wt)]


max_wt <- as.data.frame(cbind(max_x_wt, max_y_wt))
max_wt$type <- "max"
colnames(max_wt)[1] <- "bin"
colnames(max_wt)[2] <- "mean"
min_wt <- as.data.frame(cbind(min_x_wt, min_y_wt))
min_wt$type <- "min"
colnames(min_wt)[1] <- "bin"
colnames(min_wt)[2] <- "mean"

extrema_wt <- as.data.frame(rbind(min_wt, max_wt))
extrema_wt <- extrema_wt[order(extrema_wt[,1]),]
colnames(extrema_wt) <- c("bin", "mean", "type")
write.table(extrema_wt, paste0(out, "/", clustername, "-extrema_wt.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



### search three peaks for mock and wt
print("WT")
wt_middle_peak <- extrema_wt[which.max(extrema_wt$mean),1]
wt_middle_peak_idx <- which(extrema_wt$bin==wt_middle_peak)
t2 <- extrema_wt[(wt_middle_peak_idx+1):nrow(extrema_wt),]
wt_right_peak <- t2[which.max(t2$mean), 1]
wt_right_peak_idx <- which(t2$bin==wt_right_peak)
wt_right_peak_idx <- wt_middle_peak_idx + wt_right_peak_idx
t3 <- extrema_wt[1:(wt_middle_peak_idx-1),]
wt_left_peak <- t3[which.max(t3$mean), 1]
wt_left_peak_idx <- which(t3$bin==wt_left_peak)
wt_peaks <- c(wt_middle_peak, wt_right_peak, wt_left_peak)
t4 <- extrema_wt[wt_left_peak_idx:wt_right_peak_idx,]
wt_min_middle <- t4[which.min(t4$mean),1]
print(paste0("middle min peak wt ", wt_min_middle))
print(wt_peaks)



print("mock")
mock_middle_peak <- extrema_mock[which.max(extrema_mock$mean),1]
mock_middle_peak_idx <- which(extrema_mock$bin==mock_middle_peak)
t2 <- extrema_mock[(mock_middle_peak_idx+1):nrow(extrema_mock),]
mock_right_peak <- t2[which.max(t2$mean), 1]
mock_right_peak_idx <- which(t2$bin==mock_right_peak)
mock_right_peak_idx <- mock_middle_peak_idx + mock_right_peak_idx
t3 <- extrema_mock[1:(mock_middle_peak_idx-1),]
mock_left_peak <- t3[which.max(t3$mean), 1]
mock_left_peak_idx <- which(t3$bin==mock_left_peak)
mock_peaks <- c(mock_middle_peak, mock_right_peak, mock_left_peak)
t4 <- extrema_mock[mock_left_peak_idx:mock_right_peak_idx,]
mock_min_middle <- t4[which.min(t4$mean),1]
print(paste0("middle min peak mock ", mock_min_middle))
print(mock_peaks)



mock_peaks_orig_copy <- mock_peaks
wt_peaks_orig_copy <- wt_peaks



### compute always 3 peaks middle+left/right and then look if left/right have sufficient distance of the middle peaks height? if not -> delete this peak
check_peak_height_diff <- function(extrematable, peakidx) {
	if (extrematable[peakidx,1] >= 90) {
		return(FALSE)
	}
	if (extrematable[peakidx,1] <= 10) {
		return(FALSE)
	}
	
	peaky <- extrematable[peakidx, 2]

	leftmaxy <- extrematable[peakidx-2, 2]
	ydiffleft <- round(abs(peaky-leftmaxy),2)

	rightmaxy <- extrematable[peakidx+2, 2]
        ydiffright <- round(abs(peaky-rightmaxy), 2)

	leftminy <- extrematable[peakidx-1, 2]
        ydiffleftmin <- round(abs(peaky-leftminy),2)

	rightminy <- extrematable[peakidx+1, 2]
        ydiffrightmin <- round(abs(peaky-rightminy), 2)

	one_ten <- 0.10*max(extrematable$mean)
	one_ten <- ifelse(one_ten >= 0.01, round(one_ten, 2), round(one_ten, 3))
	twentyth <- 0.2*max(extrematable$mean)
	twentyth <- ifelse(twentyth >= 0.01, round(twentyth, 2), round(twentyth, 3))
	third <- 1/3*max(extrematable$mean)
	third <- ifelse(third >= 0.01, round(third, 2), round(third, 3))

				
	if (extrematable[peakidx-2,2] >= (0.95*extrematable[peakidx,2]) && abs(extrematable[peakidx-2,1]-extrematable[peakidx,1]) <= 3) {	
		x1 <- peakidx-2
		print("links genauso hoch")
	}
	else {
		x1 <- peakidx
	}
	if (extrematable[peakidx+2,2] >= (0.95*extrematable[peakidx,2]) && abs(extrematable[peakidx+2,1]-extrematable[peakidx,1]) <= 3) {
		x2 <- peakidx+2
		print("rechts genauso hoch")
	}
	else {
		x2 <-peakidx
	}


	if (round(abs(extrematable[x2+1, 2]-extrematable[x2, 2]),2) >= one_ten && peaky >= twentyth && (extrematable[x1-1, 2]+0.001) < extrematable[x1, 2]) {
                print("rechts genug tief")
		return(TRUE)
        }
	else if (round(abs(extrematable[x1-1, 2]-extrematable[x1, 2]),2) >= one_ten && peaky >= twentyth && (extrematable[x2+1,2]+0.001) < extrematable[x2,2]) {
                print("links genug tief")
		return(TRUE)
        }
	else {
		return(FALSE)
	}
					
}

print("")
print("mock")
print("left")
print(mock_peaks)
if (check_peak_height_diff(extrema_mock, mock_left_peak_idx) == FALSE) {
	mock_peaks <- mock_peaks[!mock_peaks == mock_left_peak]
}
print("right")
if (check_peak_height_diff(extrema_mock, mock_right_peak_idx) == FALSE) {
        mock_peaks <- mock_peaks[!mock_peaks == mock_right_peak]
}
for (p in mock_peaks) {
	if (p > 80) {
		mock_peaks <- mock_peaks[!mock_peaks == p]
	}
	if (p < 30) {
                mock_peaks <- mock_peaks[!mock_peaks == p]
        }
}
for (p in wt_peaks) {
        if (p > 80) {
                wt_peaks <- wt_peaks[!wt_peaks == p]
        }
        if (p < 30) {
                wt_peaks <- wt_peaks[!wt_peaks == p]
        }
}



if (length(which(extrema_mock$mean >= 0.8*max(extrema_mock$mean))) > 5) {
	print("!!!!!! achtung: Kurve ist sehr zackig -> sehr noisy")
	mock_peaks <- c()
}

if (opt$exp != "atac") {
	if ((mock_right_peak - mock_middle_peak) <= 4) {
		mock_peaks <- mock_peaks[!mock_peaks == mock_right_peak]
	}
	if ((mock_middle_peak - mock_left_peak) <= 4) {
	        mock_peaks <- mock_peaks[!mock_peaks == mock_left_peak]
	}
}



print("")
print("wt")
print("left")
if (check_peak_height_diff(extrema_wt, wt_left_peak_idx) == FALSE) {
        wt_peaks <- wt_peaks[!wt_peaks == wt_left_peak]
}
print("right")
if (check_peak_height_diff(extrema_wt, wt_right_peak_idx) == FALSE) {
        wt_peaks <- wt_peaks[!wt_peaks == wt_right_peak]
}
for (p in wt_peaks) {
        if (p > 80) {
                wt_peaks <- wt_peaks[!wt_peaks == p]
        }
}


if (opt$exp != "atac") {
	if ((wt_right_peak - wt_middle_peak) <= 4) {
	        wt_peaks <- wt_peaks[!wt_peaks == wt_right_peak]
	}
	if ((wt_middle_peak - wt_left_peak) <= 4) {
	        wt_peaks <- wt_peaks[!wt_peaks == wt_left_peak]
	}
}

if (length(which(extrema_wt$mean >= 0.85*max(extrema_wt$mean))) > 5) {
        print("!!!!!! achtung: Kurve ist sehr zackig -> sehr noisy")
        wt_peaks <- c()
}


print("writing identified peaks")
tl <- length(mock_peaks_orig_copy)+length(wt_peaks_orig_copy)
df <- cbind(rep(clustername, tl), c(rep("mock", length(mock_peaks_orig_copy)), rep("wt", length(wt_peaks_orig_copy))), c(mock_peaks_orig_copy, wt_peaks_orig_copy))
colnames(df) <- c("cluster", "cond", "peak")
write.table(df, paste0(out, "/", clustername, ".extrema"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


### classify cluster to group/transition
#group 1 = one peak --> 1
#group 2 = 2 peaks, equal height --> 3
#group 3 = 2 peaks, first higher --> 2
#group 4 = 2 peaks, second higher --> 4
#group 5 = noisy, more peaks --> 5


mock_group <- ""
wt_group <- ""
elbow_mock <- NA
elbow_wt <- NA

if (length(mock_peaks) == 1) {
	mock_group <- 1
} else if (length(mock_peaks) == 0) {
	mock_group <- 5
} else if (round(abs(extrema_mock[extrema_mock$bin==mock_peaks[1],2]-extrema_mock[extrema_mock$bin==mock_peaks[2],2]),2) < 0.01 || (length(mock_peaks)==3 && round(abs(extrema_mock[extrema_mock$bin==mock_peaks[2],2]-extrema_mock[extrema_mock$bin==mock_peaks[3],2]),2) < 0.01)) {
	mock_group <- 3
} else if (extrema_mock[mock_middle_peak_idx, 2] > extrema_mock[mock_right_peak_idx, 2] && mock_right_peak %in% mock_peaks) {
        mock_group <- 2
} else if (extrema_mock[mock_middle_peak_idx, 2] > extrema_mock[mock_left_peak_idx, 2] && mock_right_peak %in% mock_peaks) {
        mock_group <- 4
} else {
	mock_group <- 5
}

if (length(wt_peaks) == 1) {
        wt_group <- 1
} else if (length(wt_peaks) == 0) {
	wt_group <- 5
} else if (round(abs(extrema_wt[extrema_wt$bin==wt_peaks[1],2]-extrema_wt[extrema_wt$bin==wt_peaks[2],2]),2) < 0.01 || (length(wt_peaks)==3 && round(abs(extrema_wt[extrema_wt$bin==wt_peaks[2],2]-extrema_wt[extrema_wt$bin==wt_peaks[3],2]),2) < 0.01)) {
        wt_group <- 3
} else if (extrema_wt[wt_middle_peak_idx, 2] > extrema_wt[wt_right_peak_idx, 2] && wt_right_peak %in% wt_peaks) {
        wt_group <- 2
} else if (extrema_wt[wt_middle_peak_idx, 2] > extrema_wt[wt_left_peak_idx, 2] && wt_left_peak %in% wt_peaks) {
        wt_group <- 4
} else {
        wt_group <- 5
}

print("gruppe kurve")


plot(-1, -1, type="l", ylim=c(min(mockmeta, wtmeta), max(mockmeta, wtmeta)), xlim=c(1,num_total_bins), xlab="", ylab="", line.main=1, font.main=1, main=paste0(gsub("_", " ", clustername), ", ", nrow(clustergenes), " genes"), xaxt="n")
lines(range, mockmeta, type="l", col="green", lwd=1.75)
lines(range, wtmeta, type="l", col="blue", lwd=1.75)
axis(1, at=seq(1,101,12.5), labels=labs)
abline(v=mock_peaks, col="red", lty=2)
abline(v=wt_peaks, col="red", lty=2)
mtext("Corrected extrema", outer=TRUE,  cex=2, line=-1.5, font=2)
mtext(paste0("mock group = ", mock_group, ", wt group = ", wt_group), side=1, line=3, adj=0.5)





convert_to_bases <- function(number) {
	CONSTANT_BINSIZE <- 60 #58.823529
	return((number*CONSTANT_BINSIZE-3000)-60) #middle point of bases/bin
}




dev.off()




#### quantify difference between curves around TSS +- 750bp
max_peak_mock
mock_min_middle
diff_right <- 0
diff_left <- 0
if (opt$exp == "h2az") {
	for (i in (mock_min_middle+1):(mock_min_middle+1+13)) {
        	diff_right <- diff_right + (normedtable[i,3]-normedtable[i,2])
	}
	for (i in (mock_min_middle-1):(mock_min_middle-13)) {
        	diff_left <- diff_left + (normedtable[i,3]-normedtable[i,2])
	}
	diff_right <- round(diff_right, 3)
	diff_left <- round(diff_left, 3)
} else {
	for (i in (max_peak_mock+1):(max_peak_mock+13)) {
		diff_right <- diff_right + (normedtable[i,3]-normedtable[i,2])
	}
	for (i in (max_peak_mock-1):(max_peak_mock-13)) {
		diff_left <- diff_left + (normedtable[i,3]-normedtable[i,2])
	}
	diff_right <- round(diff_right, 3)
	diff_left <- round(diff_left, 3)
	if (max_peak_mock < 30) {
		diff_right <- 0
		for (i in (max_peak_wt+1):(max_peak_wt+13)) {
        		diff_right <- diff_right + (normedtable[i,3]-normedtable[i,2])
		}
		diff_left <- 0
		for (i in (max_peak_wt-1):(max_peak_wt-13)) {
        		diff_left <- diff_left + (normedtable[i,3]-normedtable[i,2])
		}
		diff_right <- round(diff_right, 3)
		diff_left <- round(diff_left, 3)
	}
}
print(paste0("rechts: ", diff_right, " und links: ", diff_left))


print("### RESULT ###")
print(paste0("cluster ", clustername, ": ", mock_group, " --> ", wt_group))



TSS <- 51

mh <- "("
md <- "("
for (p in mock_peaks) {
	mh <- paste0(mh, round(extrema_mock[which(extrema_mock$bin == p), 2],4), ",")
	md <- paste0(md, (p-TSS), ",")
}
mh <- gsub(",$", "", mh)
mh <- paste0(mh, ")")
md <- gsub(",$", "", md)
md <- paste0(md, ")")

wh <- "("
wd <- "("
for (p in wt_peaks) {
        wh <- paste0(wh, round(extrema_wt[which(extrema_wt$bin == p), 2],4), ",")
        wd <- paste0(wd, (p-TSS), ",")
}
wh <- gsub(",$", "", wh)
wh <- paste0(wh, ")")
wd <- gsub(",$", "", wd)
wd <- paste0(wd, ")")





res <- cbind(clustername, mock_group, wt_group, nrow(clustergenes), mh, wh, md, wd, (mock_middle_peak-TSS), (wt_middle_peak-TSS), diff_right, diff_left)
res <- as.data.frame(res)
colnames(res) <- c("cluster", "mock_group", "wt_group", "num_genes", "mock_heights", "wt_heights", "mock_dist", "wt_dist", "highestmock_dist", "highestwt_dist", "curve_diff_right", "curve_diff_left")

if (file.exists(paste0(out, "/cluster.classification"))) {
	write.table(res, paste0(out, "/cluster.classification"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
} else {
	write.table(res, paste0(out, "/cluster.classification"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)
}


### converted to bases
md <- paste((convert_to_bases(mock_peaks)-convert_to_bases(TSS)), collapse=",")
md <- paste0("(", md, ")")
wd <- paste((convert_to_bases(wt_peaks)-convert_to_bases(TSS)), collapse=",")
wd <- paste0("(", wd, ")")
res <- cbind(clustername, mock_group, wt_group, nrow(clustergenes), mh, wh, md, wd, convert_to_bases((mock_middle_peak-TSS)), convert_to_bases((wt_middle_peak-TSS)), diff_right, diff_left)
res <- as.data.frame(res)
if (file.exists(paste0(out, "/cluster_bases.classification"))) {
        write.table(res, paste0(out, "/cluster_bases.classification"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
} else {
        write.table(res, paste0(out, "/cluster_bases.classification"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, append=TRUE)
}




# we are done!
endNice(opt$confirmRun2EndFile, list('classifyPeaksOutputFolder' = opt$outdir))

