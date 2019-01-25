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


# reverse
data_r <- data
for(i in 1:length(data_r[,1])){
  data_r[i,] <- rev(data_r[i,])
}
d <- col(data_r) - row(data_r)
groups <- split(data_r, d)
acum_r <- c()
for(i in 1:length(groups)){
  acum_r <- c(acum, sum(unlist(groups[i]))/length(unlist(groups[i])) )
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


# same in reverse
#put into indexes to sort
indexes_r <- matrix(0, nrow=length(acum_r), ncol = 3)
indexes_r[,1] <- c(1:length(acum_r))
indexes_r[,2] <- acum_r

#sort
indexes_r[order(indexes_r[,2], decreasing=TRUE),]

#crop by percentage
n_percent <- percentage * length(acum_r)

for(i in 1:n_percent){
  indexes_r[i,3] <- 1
}

# Now resort based on first column to have them again in order
indexes_r[order(indexes_r[,1], decreasing=FALSE),]




finalmat <- matrix(0, nrow=length(data[,1]), ncol=length(data[1,]))


for(i in length(data[1,]):1){
  x <- i
  y <- 1

  if(indexes[i,3] == 1){
    while(x < length(data[1,]) && y <= length(data[,1])){
      
      
      if(x < length(data[1,]) && y <= length(data[,1])){
        finalmat[x,y] <- 1  
      }
      
      y <- y + 1
      x <- x + 1
      
      
    }  
  }
  
}

for(i in 2:length(data[,1])){
  x <- 1
  y <- i
  
  if(indexes[i,3] == 1){
    
    while(x < length(data[1,]) && y <= length(data[,1])){
      
      if(x < length(data[1,]) && y <= length(data[,1])){
        finalmat[x,y] <- 1
      }
      y <- y + 1
      x <- x + 1
      
    }  
  }
}


#same for reverse
for(i in length(data_r[1,]):1){
  x <- i
  y <- 1
  
  if(indexes[i,3] == 1){
    while(x < length(data_r[1,]) && y <= length(data_r[,1])){
      
      
      if(x < length(data_r[1,]) && y <= length(data_r[,1])){
        finalmat[x,length(data_r[,1])-y] <- 1  
      }
      
      y <- y + 1
      x <- x + 1
      
      
    }  
  }
  
}

for(i in 2:length(data_r[,1])){
  x <- 1
  y <- i
  
  if(indexes[i,3] == 1){
    
    while(x < length(data_r[1,]) && y <= length(data_r[,1])){
      
      if(x < length(data_r[1,]) && y <= length(data_r[,1])){
        finalmat[x,length(data_r[,1])-y] <- 1
      }
      y <- y + 1
      x <- x + 1
      
    }  
  }
}

# fully write last

for(i in 1:length(data[,1])){
  for(j in 1:length(data[1,])){
    if(finalmat[i,j] == 0){
      
      data[i,j] <- 0
    }
  }
}

png(paste(path_mat, ".png", sep=""), width = length(data[,1]), height = length(data[,1]))
image((data), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n')
dev.off()



