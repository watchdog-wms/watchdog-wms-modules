
options(warn=-1)

library(tidyr)
library(stringr)
library(gtools)

##############
### files for DEXSeq
##############

args = commandArgs(trailingOnly=TRUE)


dir <- args[1]
#dir <- "/mnt/raidtmp/weiss/quantCurves/H3K27ac/diffwindows/chr6-xxxxx-yyyy/"

sampleanno <- read.csv(as.character(args[2]), sep="\t", header=TRUE)
condname1 <- unique(sampleanno$condition)[1]
condname2 <- unique(sampleanno$condition)[2]
c1_c2 <- paste0(condname1, "_", condname2)
c2_c1 <- paste0(condname2, "_", condname1)

givenstrand <- as.character(args[3])


############################
### pooled amss from both dirs and filled up
###########################


annot_fc <- paste0(dir, "/annotation_fc_filled_up.single")
annot <- paste0(dir, "/annotation_filled_up.single")


files <- list.files(path=dir, pattern="*final_amss*")
files <- paste0(dir, "/", files)

file_c1 <- files[grepl(c1_c2, files)]
file_c2 <- gsub(c1_c2, c2_c1, file_c1)


combinedataframes <- function(f1, f2) {
  d1 <- read.csv(f1, sep="\t", header=TRUE)
  d1 <- subset(d1, !is.na(d1$start_seq_idx))
  d2 <- read.csv(f2, sep="\t", header=TRUE)
  d2 <- subset(d2, !is.na(d2$start_seq_idx))
  if (nrow(d1)>0) {
    if (nrow(d2)>0) {
      #beide
      arr <- strsplit(f1, "/")[[1]]
      arr <- arr[arr!=""]
      dir <- arr[grepl("tsv", arr)]
      dir <- ifelse(grepl(c1_c2, dir), c1_c2, c2_c1)
      diff <- arr[length(arr)-2]
      diff <- ifelse(diff=="windows", "diff", "nondiff")
      d1$dir <- dir
      d1$diff <- diff
      
      arr2 <- strsplit(f2, "/")[[1]]
      arr2 <- arr2[arr2!=""]
      dir <- arr2[grepl("tsv", arr2)]
      dir <- ifelse(grepl(c1_c2, dir), c1_c2, c2_c1)
      diff <- arr2[length(arr2)-2]
      diff <- ifelse(diff=="windows", "diff", "nondiff")
      d2$dir <- dir
      d2$diff <- diff
      
      d <- as.data.frame(rbind(d1, d2))
      d <- d[order(d$start_seq_idx),]
      return(d)
    } else {
      #d1
      arr <- strsplit(f1, "/")[[1]]
      arr <- arr[arr!=""]
      diff <- arr[length(arr)-2]
      diff <- ifelse(diff=="windows", "diff", "nondiff")
      d1$dir <- c1_c2
      d1$diff <- diff
      d1 <- d1[order(d1$start_seq_idx),]
      return(d1)
    }
  } else {
    if (nrow(d2)>0) {
      #d2
      arr2 <- strsplit(f2, "/")[[1]]
      arr2 <- arr2[arr2!=""]
      diff <- arr2[length(arr2)-2]
      diff <- ifelse(diff=="windows", "diff", "nondiff")
      d2$dir <- c2_c1
      d2$diff <- diff
      d2 <- d2[order(d2$start_seq_idx),]
      return(d2)
    } else {
      #none
      print("no lines")
      return(NA)
    }
  }
}

fillup_from_start2end <- function(s, e, df) {
  allnums <- seq(s, e, 1)
  allnums <- allnums[!(allnums %in% df$V4)]
  incl <- c()
  df$V4 <- as.numeric(as.character(df$V4))
  df$V5 <- as.numeric(as.character(df$V5))
  incl <- c(incl, unlist(apply(df[,c(4,5)], 1, function(x) {seq(x[1], x[2],1)})))
  allnums <- allnums[!(allnums %in% incl)]
  missingseqs <- split(allnums, cumsum(c(TRUE, diff(allnums)!=1)))
  missingidx <- c()
  missingidx <- c(missingidx, unlist(sapply(missingseqs, function(x) range(x))))
  missingidx <- as.data.frame(matrix(missingidx, ncol=2, byrow=TRUE))
  missingidx$c <- df$V1[1]
  missingidx$n <- "dexseq"
  missingidx$d <- "exonic_part"
  missingidx$p <- "."
  missingidx$s <- givenstrand
  missingidx$pp <- "."
  return(missingidx)
}

together <- function(file, file2, fc) {
  tmp <- combinedataframes(file, file2)
  if (!is.na(tmp) && nrow(tmp)>0) {
    arr <- strsplit(file, "/")[[1]]
    arr <- arr[arr!=""]
    pos <- arr[grepl("chr", arr)]
    chr <- strsplit(pos, "-")[[1]][1]
    start <- as.numeric(as.character(strsplit(pos, "-")[[1]][2]))
    end <- as.numeric(as.character(strsplit(pos, "-")[[1]][3]))
    geneid <- paste0(chr, "-", start, "-", end)
    
    entry <- c(chr, "dexseq", "aggregate_gene", start, end, ".", givenstrand, ".", 
               paste0("gene_id \"", geneid, "\""))
    gff <- entry
    #print(tmp)
    for (line in 1:nrow(tmp)) {
      entry <- c(chr, "dexseq", "exonic_part", tmp$start_seq_idx[line], tmp$end_seq_idx[line], ".", givenstrand, ".", 
                 paste0("transcripts \"", 
                        geneid, "-", tmp$start_seq_idx[line], "-", tmp$end_seq_idx[line], "-", 
                        tmp$dir[line], "-", tmp$diff[line],
                        "\"; exonic_part_number \"", 
                        line, "\"; gene_id \"", geneid, "\""))
      if (fc) {
        entry <- str_replace(entry,"transcripts \"", paste0("htname \"",geneid,":",line,"\"; transcripts \""))
      }
      gff <- rbind(gff, entry)
    }
    gff <- as.data.frame(gff)
    rownames(gff) <- NULL
    gff$V4 <- as.numeric(as.character(gff$V4))
    gff <- gff[order(gff$V4),] #sort start positions asc
   

    missing_regs <- fillup_from_start2end(start,end,gff[2:nrow(gff),])
    fixed_cols <- c(chr, "dexseq", "exonic_part", ".", gff$V7[1], ".", NA)
    missing_regs <- cbind(missing_regs, rbind(fixed_cols))
    missing_regs <- missing_regs[,c(3,4,5,1,2,6,7,8,9)]
    colnames(missing_regs) <- colnames(gff)
    missing_regs$V9 <- NA
    filled <- rbind(gff, missing_regs)
    filled <- filled[order(filled$V4),]
    rownames(filled) <- NULL

    exline <- filled[!is.na(filled$V9),9][2]
    filled$V9 <- apply(filled, 1, function(x) {
      if (is.na(x[9])) {
	newcor <- paste0(x[4], "-", x[5])
	if (!grepl(c1_c2, exline)) {
		exline <- str_replace(exline, c2_c1, c1_c2)
	}
        news <- str_replace(exline, paste0("\\d+-\\d+-", c1_c2, "-diff"),
                           paste0(newcor, "-filled-nondiff"))
	return(news)
      } else {
        return(x[9])
      }
    })
    c <- 1
    filled$V9 <- apply(filled, 1, function(x) {
      idx <- ifelse(c<10, paste0(0, c), c)
      if (grepl(":", x[9])) {
        x[9] <- str_replace(x[9], ":\\d+", paste0(":", idx))
      }
      if (grepl("exonic", x[9])) {
        x[9] <- str_replace(x[9], "number \"\\d+", paste0("number \"", idx))
        c <<- c+1
      }
      return(x[9])
    })
    
    if (grepl("Inf", filled$V4)) {
      print("ERROR should not have Inf as start positions")
      print("")
    }
    #print(filled)
    return(filled)
  } else {
	  return(NA)
  }
}



starttime <- Sys.time()
sorted_filledup <- together(file_c1, file_c2, TRUE)
if (!is.na(sorted_filledup)) {
  write.table(sorted_filledup, annot_fc, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
print(paste0("TIME annotation FC ", (Sys.time()-starttime)))



starttime2 <- Sys.time()
sorted_filledup <- together(file_c1, file_c2, FALSE)
if (!is.na(sorted_filledup)) {
  write.table(sorted_filledup, annot, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
print(paste0("TIME annotation ", (Sys.time()-starttime2)))

