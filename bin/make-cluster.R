#!/usr/bin/env Rscript
library(ape)

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("USE: Rscript --vanilla make-cluster.R <csv>")
}


path = args[1]

inputcsv <- read.csv(path, sep ="-", header=FALSE)

# Make indexing table

indexing_table <- table(inputcsv[,1])
n_species <- length(indexing_table)

for(i in 1:n_species){
  indexing_table[i] <- i
}


distmat <- matrix(NA, nrow=n_species, ncol=n_species)
rownames(distmat) <- rownames(indexing_table)
colnames(distmat) <- rownames(indexing_table)

for(i in 1:n_species){
  distmat[i,i] <- 0
}


for(i in 1:length(inputcsv[,1])){
  redirX <- (as.character(inputcsv[i,1]))
  redirY <- (as.character(inputcsv[i,2]))
  
  distmat[indexing_table[redirX],indexing_table[redirY]] <- as.numeric(inputcsv[i,3])
  distmat[indexing_table[redirY],indexing_table[redirX]] <- as.numeric(inputcsv[i,3])
}

cluster <- hclust(dist(distmat), method = "average")

# Rooted hierarchical cluster
#plot(cluster, main = "Clustering based on automatic scoring", xlab = "Organisms")

# Unrooted
plot(as.phylo(cluster), type = "unrooted", cex = 0.6, main = "Unrooted clustering based on automatic")

# Fan cluster
#plot(as.phylo(cluster), type = "fan")



