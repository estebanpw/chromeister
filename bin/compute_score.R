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
	  
	  cmax_pos <- which.max(aux_density[i,])
  
  if((aux_density[i,cmax_pos]) > 0){
	      score_density[i,] <- 0
      score_density[i,cmax_pos] <- 1
          pmax_pos <- cmax_pos
        }
}


pmax_pos <- which.max(aux_density[1,])


for(i in 1:len_i){
	  
	  cmax_pos <- which.max(aux_density[,i])
  
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
dist_th <- 3
for(i in 2:len_i){
	  
	  cmax_pos <- which.max(score_copy[,i])
  # taxicab distance
  distance <- abs(cmax_pos - pmax_pos)

    if(max(score_copy[,i]) == 0){
	        distance <- len_i
    }
    
    if(distance > dist_th){
	        score <- score + distance
      }

      pmax_pos <- cmax_pos
}


score <- (score/(len_i^2))



print(score)

png(paste(path_mat, ".filt.png", sep=""), width = length(data[,1]), height = length(data[,1]))
image(t(score_copy), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n', main = paste(fancy_name, paste("filt. score=", score)))
dev.off()


# Output pixel coordinates of highly conserved signals

# To clear it in case it existed
write(c("X", "Y"), file = paste("hits-XY-", fancy_name, sep=""), append = FALSE, sep = "\n")

for(i in 1:(length(score_copy[,1]))){
	for(j in 1:(length(score_copy[1,]))){
		if(score_copy[i,j] != 0){
			write(c(paste(i, j)), file = paste("hits-XY-", fancy_name, sep=""), append = TRUE, sep = "\n")
		}
	}
}



