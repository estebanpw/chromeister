#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
	  stop("USE: Rscript --vanilla plot.R <matrix>")
}


path_mat = args[1]


fancy_name <- strsplit(path_mat, "/")
fancy_name <- (fancy_name[[1]])[length(fancy_name[[1]])]





data <- as.matrix(read.csv(path_mat, sep = " "))


# png(paste(path_mat, "_raw.png", sep=""), width = length(data[,1]), height = length(data[,1]))
# image((data), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n', main = paste(fancy_name, " before filtering"))
# dev.off()


# Max of columns

len_i <- length(data[,1])
len_j <- length(data[1,])

pmax_pos <- which.max(data[,1])


for(i in 1:len_i){
	  
	  cmax_pos <- which.max(data[,i])
  
  if((data[cmax_pos,i]) > 0){
	      data[,i] <- 0
      data[cmax_pos,i] <- 1
          
          
          # if(i > 1 && i < len_i && cmax_pos > 1 && cmax_pos < len_j){
          #   data[cmax_pos, i-1] <- 1
          #   data[cmax_pos, i+1] <- 1
          #   data[cmax_pos-1, i] <- 1
          #   data[cmax_pos+1, i] <- 1
          # }
          
          pmax_pos <- cmax_pos
        }
}


pmax_pos <- which.max(data[1,])
for(i in 1:len_i){
	  
	  cmax_pos <- which.max(data[i,])
  
  if((data[i, cmax_pos]) > 0){
	      data[i,] <- 0
      data[i, cmax_pos] <- 1
          
          
          # if(i > 1 && i < len_i && cmax_pos > 1 && cmax_pos < len_j){
          #   data[cmax_pos, i-1] <- 1
          #   data[cmax_pos, i+1] <- 1
          #   data[cmax_pos-1, i] <- 1
          #   data[cmax_pos+1, i] <- 1
          # }
          
          pmax_pos <- cmax_pos
        }
}



# it was here before!!

#png(paste("max-metric/", paste(path_mat, "_post.bmp", sep="")), width = length(data[,1]), height = length(data[,1]))
#image(t(data), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n', main = paste(path_mat, paste("after. score=", "Unknown")))
#dev.off()

# Convolutional kernel (use 1 for a 3 by 3 kernel)
dim_kernel <- 8

# To compute the score
score <- 0
dist_th <- 10
i_prev <- 0
for(i in 1:len_i){
	  
	  cmax_pos <- which.max(data[,i])
  
  if(i > 2*dim_kernel+1 && i < len_i-(2*dim_kernel+1) && cmax_pos > 2*dim_kernel+1 && cmax_pos < len_j - (2*dim_kernel+1)){
	      #data[cmax_pos, i] <- 0 # temporary set to zero
	      # kernel
	      v <- matrix(0, nrow = dim_kernel+1+dim_kernel, dim_kernel+1+dim_kernel)
      
      idx1 <- 1
          for(k in (cmax_pos-dim_kernel):(cmax_pos+dim_kernel)){
		        idx2 <- 1
            for(w in (i-dim_kernel):(i+dim_kernel)){
		            v[idx1, idx2] <- data[k,w]
	            idx2 <- idx2 + 1
		          }
	          idx1 <- idx1 + 1
	        }
          #print(v)
          conv <- sum(v * 1/64)
	      if(conv >= 2/64){
		            data[cmax_pos, i] <- 1
	        i_prev <- i
		    }else{
			          data[cmax_pos, i] <- 0
		    }
	      #print(conv)
	      
	    }
    
    # taxicab distance
    distance <- abs(cmax_pos - pmax_pos) + abs(i_prev - i)
    
    if(distance > dist_th){
	        if(max(data[,i]) == 0){
			      distance <- len_i
        }
        score <- score + distance
	  }
      
      pmax_pos <- cmax_pos
}

# Growing pixels
copydata <- data
for(i in 1:len_i){
	  
	  cmax_pos <- which.max(data[,i])
  
  
  if(data[cmax_pos, i] > 0 && i > 1 && i < len_i && cmax_pos > 1 && cmax_pos < len_j){
	      copydata[cmax_pos, i-1] <- 1
      copydata[cmax_pos, i+1] <- 1
          copydata[cmax_pos-1, i] <- 1
          copydata[cmax_pos+1, i] <- 1
	    }
    
}

# Print score
#print(paste("The score is ", score/(len_i^2)))

score <- (score/(len_i^2))
# -0.0185714 + 1.17143 x + 3.57143 x^2

# Quadratic fit
#score <- 1.17143*score + 3.57143*(score^2)
#score <- max(0, min(1, score))
print(score)


png(paste(path_mat, ".filt.png", sep=""), width = length(data[,1]), height = length(data[,1]))
image(t(copydata), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n', main = paste(fancy_name, paste("filt. score=", score)))
dev.off()



