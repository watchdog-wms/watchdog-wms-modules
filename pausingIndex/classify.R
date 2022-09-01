#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

df <- read.csv(args[1], sep="\t", header=TRUE)

df[,4] <- as.numeric(as.character(df[,4]))
df[,6] <- as.numeric(as.character(df[,6]))


df$group <- ifelse(df$PI >= 1, ifelse(df$body_rpkm >= 2, "paused_expressed", "paused_notExpressed"), ifelse(df$body_rpkm >= 2, "notPaused_expressed", "notPaused_notExpressed"))
df$group[is.na(df$group)] <- "too_low"


df <- df[order(df$PI, decreasing=TRUE, na.last=TRUE), ]

write.table(df, paste0(args[2], "/result_sense.table"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



