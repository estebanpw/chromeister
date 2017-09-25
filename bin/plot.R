#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("USE: Rscript --vanilla plot.R <matrix>")
}


path_mat = args[1]


data <- read.csv(path_mat, header = FALSE, sep = " ", skip=1)


png(paste(path_mat, ".png", sep=""), width = length(data[,1]), height = length(data[,1]))
image(t(as.matrix(data)), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n')
dev.off()
