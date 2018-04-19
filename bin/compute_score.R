
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("USE: Rscript --vanilla plot.R <matrix>")
}



path_mat = args[1]


fancy_name <- strsplit(path_mat, "/")
fancy_name <- (fancy_name[[1]])[length(fancy_name[[1]])]





data <- as.matrix(read.csv(path_mat, sep = " "))


# Max of columns

len_i <- 1000
len_j <- 1000


score_density <- data
aux_density <- data


pmax_pos <- which.max(aux_density[,1])
for(i in 1:len_i){
  
  cmax_pos <- which.max(aux_density[i,]) # get max of row
  
  if((aux_density[i,cmax_pos]) > 0){ # if it has value
    #row_from_col <- which.max(aux_density[,cmax_pos]) # get max of the column pointed by maximum of row
    score_density[i,] <- 0 # put all others in row to 0 i.e. only use this max, EXCEPT for the maximum in the column
    #score_density[row_from_col, cmax_pos] <- 1
    score_density[i,cmax_pos] <- 1
    pmax_pos <- cmax_pos
  }
}


pmax_pos <- which.max(aux_density[1,])


for(i in 1:len_i){

  cmax_pos <- which.max(aux_density[,i]) # get max of column

  if((aux_density[cmax_pos,i]) > 0){
    score_density[,i] <- 0
    score_density[cmax_pos,i] <- 1
    pmax_pos <- cmax_pos
  }
}









score_copy <- score_density
diag_len <- 4

for(i in 6:(len_i-5)){
  for(j in 6:(len_j-5)){
    
    value <- 0
    for(w in (-diag_len/2):(diag_len/2)){
      if(score_density[i+w,j+w] > 0){
        value <- value + 1
      }
    }
    
    if(value >= diag_len){
      for(k in 1:5){
        score_copy[i+k,j+k] <- 1
      }
      for(k in 1:5){
        score_copy[i-k,j-k] <- 1
      }
    }else if(score_copy[i,j]==0){
      score_copy[i,j] <- 0
    }
  }
}

for(i in 6:(len_i-5)){
  for(j in 6:(len_j-5)){
    
    value <- 0
    for(w in (-diag_len/2):(diag_len/2)){
      if(score_density[i-w,j+w] > 0){
        value <- value + 1
      }
    }
    
    if(value >= diag_len){
      for(k in 1:5){
        score_copy[i-k,j+k] <- 1
      }
      for(k in 1:5){
        score_copy[i+k,j-k] <- 1
      }
    }else if(score_copy[i,j]==0){
      score_copy[i,j] <- 0
    }
  }
}


# Kernel to remove single points

for(i in 1:(length(score_copy[,1]))){
  for(j in 1:(length(score_copy[1,]))){
    
    value <- 0
    
    min_i <- max(1, i-1)
    max_i <- min(len_i, i+1)
    min_j <- max(1, j-1)
    max_j <- min(len_j, j+1)
    
    value <- sum(score_copy[min_i:max_i, min_j:max_j])
    
    if(value < 2) score_copy[i,j] <- 0
  }
}

# To compute the score
score <- 0
pmax_pos <- which.max(score_copy[,1])
dist_th <- 1.5
besti <- 1
bestj <- pmax_pos
dvec1 <- abs(which.max(score_copy[,2]) - which.max(score_copy[,1]))
dvec2 <- abs(which.max(score_copy[,3]) - which.max(score_copy[,2]))
dvec3 <- abs(which.max(score_copy[,4]) - which.max(score_copy[,3]))
for(i in 5:len_i){
  
  #print(paste(paste(paste(dvec1, dvec2), dvec3), mean(c(dvec1, dvec2, dvec3))))
  distance <- mean(c(dvec1, dvec2, dvec3))

  
  
  dvec1 <- dvec2
  dvec2 <- dvec3
  dvec3 <- abs(which.max(score_copy[,i]) - which.max(score_copy[,i-1]))
  
  # If there is a 0 or we are too far away just add max distance!
  if(distance > dist_th || distance == 0){
    score <- score + len_i
  }
  
}


score <- (score/(len_i^2))
print(score)


png(paste(path_mat, ".filt.png", sep=""), width = length(data[,1]), height = length(data[,1]))
image(t(score_copy), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n', main = paste(fancy_name, paste("filt. score=", score)))
dev.off()


# Output pixel coordinates of highly conserved signals

# To clear it in case it existed
write(c(paste("X", "Y")), file = paste("hits-XY-", paste(fancy_name, ".hits", sep=""), sep=""), append = FALSE, sep = "\n")

for(i in 1:(length(score_copy[,1]))){
  for(j in 1:(length(score_copy[1,]))){
    if(score_copy[i,j] != 0){
      write(c(paste(i, j)), file = paste("hits-XY-", paste(fancy_name, ".hits", sep=""), sep=""), append = TRUE, sep = "\n")
    }
  }
}