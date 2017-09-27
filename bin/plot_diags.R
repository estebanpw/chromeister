#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("USE: Rscript --vanilla plot.R <matrix>")
}


path_mat = args[1]




data <- as.matrix(read.csv(path_mat, sep = " "))
percentage <- 0.1




d <- col(data) - row(data)
groups <- split(data, d)
acum <- c()
for(i in 1:length(groups)){
  acum <- c(acum, sum(unlist(groups[i])))
}

indexes <- matrix(0, nrow=length(acum), ncol = 3)
indexes[,1] <- c(1:length(acum))
indexes[,2] <- acum

#sort
indexes[order(indexes[,2], decreasing=TRUE),]

#crop by percentage
n_percent <- percentage * length(acum)

for(i in 1:n_percent){
  indexes[i,3] <- 1
}

# Now resort based on first column to have them again in order
indexes[order(indexes[,1], decreasing=FALSE),]

lastacum <- 0

for(i in length(data[1,]):1){
  x <- i
  y <- 1

  if(indexes[i,3] == 0){
    while(x < length(data[1,]) && y <= length(data[,1])){
      
      
      if(x < length(data[1,]) && y <= length(data[,1])){
        data[x,y] <- 0  
      }
      
      y <- y + 1
      x <- x + 1
      
      
    }  
  }
  lastacum <- i
  
}

for(i in 2:length(data[,1])){
  x <- 1
  y <- i
  
  if(indexes[i,3] == 0){
    
    while(x < length(data[1,]) && y <= length(data[,1])){
      
      if(x < length(data[1,]) && y <= length(data[,1])){
        data[x,y] <- 0
      }
      y <- y + 1
      x <- x + 1
      
    }  
  }
}

png(paste(path_mat, ".png", sep=""), width = length(data[,1]), height = length(data[,1]))
image((data), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n')
dev.off()



