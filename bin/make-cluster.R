#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("USE: Rscript --vanilla make-cluster.R <csv>")
}


path = args[1]

inputcsv <- read.csv(path, sep ="-")

# Make indexing table

indexing_table <- table(inputcsv[,1])
n_species <- length(indexing_table)

for(i in 1:n_species){
  indexing_table[i] <- i
}


distmat <- matrix(1, nrow=n_species, ncol=n_species)
rownames(distmat) <- rownames(indexing_table)
colnames(distmat) <- rownames(indexing_table)


for(i in 1:length(inputcsv[,1])){
  redirX <- (as.character(inputcsv[i,1]))
  redirY <- (as.character(inputcsv[i,2]))
  
  distmat[indexing_table[redirX],indexing_table[redirY]] <- as.numeric(inputcsv[i,3])
  distmat[indexing_table[redirY],indexing_table[redirX]] <- as.numeric(inputcsv[i,3])
}

cluster <- hclust(dist(distmat))
plot(cluster, main = "Clustering based on automatic scoring", xlab = "Organisms")

