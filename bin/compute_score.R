#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 1){
  stop("USE: Rscript --vanilla plot.R <matrix>")
}


path_mat = args[1]




data <- as.matrix(read.csv(path_mat, sep = " "))


derivative <- function(vector){
  acum <- 0
  for(i in 2:length(vector)){
    acum <- acum + (vector[i] - vector[i-1])
  }
  
  return(acum)
}

deviation_keeper <- function(vector, th_sd){
  m <- mean(vector)
  sdev <- sd(vector)
  kept <- c()
  last_valid <- m
  
  for(i in 1:length(vector)){
    if((vector[i] - m)/sdev < th_sd){
      kept <- c(kept, vector[i])
      last_valid <- vector[i]
    }else{
      kept <- c(kept, last_valid)
      
    }
  }
  return (kept)
}

# accumulate hits
x <- 0
y <- 0



#adulterado

# for(i in 1:length(data[,1])){
#   data[i,i] <- 1000
# }
acum <- c()
acum_r <- c()

d <- col(data) - row(data)
groups <- split(data, d)
acum <- c()
for(i in 1:length(groups)){
  acum <- c(acum, sum(unlist(groups[i]))/length(unlist(groups[i])) )
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
  acum_r <- c(acum_r, sum(unlist(groups[i]))/length(unlist(groups[i])) )
}


# adulterado
acum <- acum[250:1750]
acum_r <- acum_r[250:1750]

#acum <- acum/(max(acum))
#acum_r <- acum_r/(max(acum_r))

acum_total <- (acum + acum_r)/2

#std_acum <- rnorm(length(acum_total), mean = mean(acum_total), sd = sd(acum_total))

# squared difference

smooth <- deviation_keeper(acum_total, 3)

sqd_diff <- sqrt(abs(acum_total^2 - smooth^2))



#plot(ecdf(acum_total), main = "TOTAL", xlim = c(0, max(acum_total, smooth)))
#plot(ecdf(smooth), main = "STANDARD", xlim = c(0, max(acum_total, smooth)))


# plot(acum, type="n", xlim = c(0, length(acum)), ylim = c(0, max(acum)), ylab= "Hits", xlab = "Length in bins", main = paste("Forward Histogram of hits ", name))
# for(i in 2:length(acum)){
# 
#   lines(c(i-1,i), c(acum[i-1],acum[i]), type="l")
# }
# 
# plot(acum_r, type="n", xlim = c(0, length(acum_r)), ylim = c(0, max(acum_r)), ylab= "Hits", xlab = "Length in bins", main = paste("Reverse Histogram of hits ", name))
# for(i in 2:length(acum_r)){
# 
#   lines(c(i-1,i), c(acum_r[i-1],acum_r[i]), type="l")
# }

#plot(acum_total, type="n", xlim = c(0, length(acum_total)), ylim = c(0, max(acum_total, smooth)), ylab= "Hits", xlab = "Length in bins", main = paste("TOTAL Histogram of hits ", name))
#for(i in 2:length(acum_total)){
  
#  lines(c(i-1,i), c(acum_total[i-1],acum_total[i]), type="l")
#}


#plot(smooth, type="n", xlim = c(0, length(smooth)), ylim = c(0, max(acum_total, smooth)), ylab= "Hits", xlab = "Length in bins", main = paste("TOTAL Histogram of hits ", name))
#for(i in 2:length(smooth)){
  
#  lines(c(i-1,i), c(smooth[i-1],smooth[i]), type="l")
#}


# plot(std_acum, type="n", xlim = c(0, length(std_acum)), ylim = c(0, max(acum_total, std_acum)), ylab= "Hits", xlab = "Length in bins", main = paste("STANDARD Histogram of hits ", name))
# for(i in 2:length(std_acum)){
# 
#   lines(c(i-1,i), c(std_acum[i-1],std_acum[i]), type="l")
# }






#lines(c(0, length(acum)/2), c(0, max(acum)), type ="l", col = "red", lwd=3)
#lines(c(length(acum)/2, length(acum)), c(max(acum),0), type ="l", col = "red", lwd=3)


















































# 
# #put into indexes to sort
# indexes <- matrix(0, nrow=length(acum), ncol = 3)
# indexes[,1] <- c(1:length(acum))
# indexes[,2] <- acum
# indexes[,3] <- 0
# 
# #sort
# indexes <- indexes[order(indexes[,2], decreasing=TRUE),]
# 
# #crop by percentage
# n_percent <- percentage * length(acum)
# 
# for(i in 1:n_percent){
#   indexes[i,3] <- 1
# }
# 
# # Now resort based on first column to have them again in order
# indexes <- indexes[order(indexes[,1], decreasing=FALSE),]
# 
# # same in reverse
# #put into indexes to sort
# indexes_r <- matrix(0, nrow=length(acum_r), ncol = 3)
# indexes_r[,1] <- c(1:length(acum_r))
# indexes_r[,2] <- acum_r
# 
# #sort
# indexes_r <- indexes_r[order(indexes_r[,2], decreasing=TRUE),]
# 
# #crop by percentage
# n_percent <- percentage * length(acum_r)
# 
# for(i in 1:n_percent){
#   indexes_r[i,3] <- 1
# }
# 
# # Now resort based on first column to have them again in order
# indexes_r <- indexes_r[order(indexes_r[,1], decreasing=FALSE),]
# 
# 
# 
# 
# 
# finalmat <- matrix(0, nrow=length(data[,1]), ncol=length(data[1,]))
# 
# curr_diag <- 1
# for(i in length(data[1,]):1){
#   x <- i
#   y <- 1
# 
#   if(indexes[curr_diag,3] == 1){
#     
#     while(x < length(data[1,]) && y <= length(data[,1])){
#       
#       
#       if(x < length(data[1,]) && y <= length(data[,1])){
#         finalmat[x,y] <- 1  
#       }
#       
#       y <- y + 1
#       x <- x + 1
#       
#       
#     }  
#   }
#   curr_diag <- curr_diag + 1
# }
# 
# 
# 
# for(i in 2:length(data[,1])){
#   x <- 1
#   y <- i
#   
#   if(indexes[i,3] == 1){
#     
#     while(x < length(data[1,]) && y <= length(data[,1])){
#       
#       if(x < length(data[1,]) && y <= length(data[,1])){
#         finalmat[x,y] <- 1
#       }
#       y <- y + 1
#       x <- x + 1
#       
#     }  
#   }
#   curr_diag <- curr_diag + 1
# }
# 
# 
# 
# 
# 
# curr_diag <- length(indexes_r[,1])
# #same for reverse
# for(i in length(data_r[1,]):1){
#   x <- i
#   y <- 1
# 
#   if(indexes_r[curr_diag,3] == 1){
#     while(x < length(data_r[1,]) && y <= length(data_r[,1])){
# 
# 
#       if(x < length(data_r[1,]) && y <= length(data_r[,1])){
#         finalmat[x,length(data_r[,1])-y] <- 1
#       }
# 
#       y <- y + 1
#       x <- x + 1
# 
# 
#     }
#   }
#   curr_diag <- curr_diag - 1
# 
# }
# for(i in 2:length(data_r[,1])){
#   x <- 1
#   y <- i
# 
#   if(indexes_r[curr_diag,3] == 1){
# 
#     while(x < length(data_r[1,]) && y <= length(data_r[,1])){
# 
#       if(x < length(data_r[1,]) && y <= length(data_r[,1])){
#         finalmat[x,length(data_r[,1])-y] <- 1
#       }
#       y <- y + 1
#       x <- x + 1
# 
#     }
#   }
#   curr_diag <- curr_diag - 1
# }
# 
# 
# # fully write last
# 
# for(i in 1:length(data[,1])){
#   for(j in 1:length(data[1,])){
#     if(finalmat[i,j] == 0){
#       
#       data[i,j] <- 0
#     }else{
#       data[i,j] <- 1
#     }
#   }
# }
# 














# png(paste(name, ".bmp", sep=""), width = length(data[,1]), height = length(data[,1]))
# 
# image((data), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n')
# dev.off()


#print(ks.test(acum_total, smooth))
# The highest score is length(acum_total) i.e. the case when all is 1 - 0, however this is not possible.
print(sum(sqd_diff))




