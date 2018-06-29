
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
  stop("USE: Rscript --vanilla plot.R <matrix> <matsize>")
}






growing_regions <- function(mat, reward = 6, penalty = 5, sidePenalty = 3, MAXHSPS = 500, TH = 10, WSIZE = 7){

  #write(mat, file ="data.txt", ncolumns = 200)
  
  len <- min(length(mat[1,]), length(mat[,1]))
  HSPS <- matrix(0, nrow=MAXHSPS, ncol=5)
  idx <- 1
  lH <- round(WSIZE/2) - 1
  rH <- round(WSIZE/2) + 1
  
  
  #print(paste("lH: ", lH, "rH: ", rH))
  
  if(WSIZE %% 2 == 0){
    print("WSIZE MUST BE ODD")
    stop()
  }
  
  
  
  i <- 1
  #readline(prompt="Press [enter] to continue")
  while(i < len){

    value <- max(mat[i,]) * reward
    if(value == 0) i <- i + 1
    pos <- which.max(mat[i,])
    # these two hold ending frag
    endfrag <- pos
    j <- i
    count_penalties <- 1

    #print("-----------------")

    while(value > 0 && j < len){

      #print(paste("took values ", j, endfrag, "which have score of ", mat[j, endfrag], "current value is: ", value))
      
      
      

      # Reset position used
      
      mat[max(1,j-1), max(1,endfrag-2)] <- 0
      mat[max(1,j-1), max(1,endfrag-1)] <- 0
      mat[max(1,j-1), endfrag] <- 0
      mat[max(1,j-1), min(len, endfrag+1)] <- 0
      mat[max(1,j-1), min(len, endfrag+2)] <- 0
      
      mat[j, max(1,endfrag-2)] <- 0
      mat[j, max(1,endfrag-1)] <- 0
      mat[j, endfrag] <- 0
      mat[j, min(len, endfrag+1)] <- 0
      mat[j, min(len, endfrag+2)] <- 0
      
      #print(paste("Erasing: (",max(1,j-1), max(1,endfrag-2),")(",max(1,j-1), max(1,endfrag-1),")(",max(1,j-1), endfrag,
      #            ")(",max(1,j-1), min(len, endfrag+1),")(",max(1,j-1), min(len, endfrag+2),")(",j, max(1,endfrag-2),
      #            ")(",j, max(1,endfrag-1),")(",j, endfrag,")(",j, min(len, endfrag+1),")(",j, min(len, endfrag+2),")"))
      
      
      # Go for next
      j <- j + 1
      # Check next, could be reverse or forward
      mleft <- max(1, endfrag-lH)
      mright <- min(len, endfrag+lH)
      window <- mat[j, mleft:mright]
      

      #print(paste("mleft: ", mleft, "mright", mright, "j is: ", j))
      #print(window)

      v <- max(window)
      selected <- which.max(window)
      # Make it rather go diagonally
      #print(paste("WIndow len: ", length(window)))
      chose_diagonal <- FALSE
      if(length(window) == WSIZE && v == window[lH]){
        selected <- lH
        chose_diagonal <- TRUE
      }
      
      if(length(window) == WSIZE && v == window[rH]){
        selected <- rH
        chose_diagonal <- TRUE
      }
      
      #print(paste("Selected value be like ", selected))
      

      if(v != 0){
        
        
        
        endfrag <- (mleft + selected ) # To make the indexing
        if(length(window) == WSIZE) endfrag <- endfrag - 1
        endfrag <- max(1, min(len, endfrag))
        #print(paste("\t", " endfragnew = ", endfrag, "max of window: ", max(window), "on position", which.max(window)))
      }

      #print(paste("Chose diagonal is ", chose_diagonal))
      
      # If no similarity is found
      if(v == 0){
        value <- value - count_penalties * penalty
        count_penalties <- count_penalties + 1
        #print("Got penalty #########")
        
      }else{
        
        # Similarity is found
        
        if(!chose_diagonal){
          value <- value + count_penalties * (-sidePenalty)
          count_penalties <- count_penalties + 1
          #print("Got SIDE penalty @@@ #########")
        }else{
          count_penalties <- 1
          value <- value + reward
          #print("Got reward thou #########")
        }
        
      }
      #readline(prompt="Press [enter] to continue")

    }
    # Check len of frag
    if(j-i > TH){
      # HSPS[idx, 1] <- i
      # HSPS[idx, 2] <- pos
      # HSPS[idx, 3] <- j
      # HSPS[idx, 4] <- endfrag
      # HSPS[idx, 5] <- abs(i-j)
      
      HSPS[idx, 1] <- pos
      HSPS[idx, 2] <- i
      HSPS[idx, 3] <- endfrag
      HSPS[idx, 4] <- j
      HSPS[idx, 5] <- abs(i-j)
      
      
      
      #print(paste("ACCEPT", i, pos, j, endfrag, "x-y", j-i, abs(endfrag-pos), sep = " "))
      idx <- idx + 1
    }else{
      #print("REJECT")
    }
    # This will prevent overlappign lines, I think
    #i <- j
    
    if(idx == MAXHSPS) break
  }
  
  
  
  
  return (HSPS)
}

detect_events <- function(HSPS, sampling){
  
  DIAG_SEPARATION <- 10
  # same as HSPS but adding the event
  output <- matrix(0, nrow=length(HSPS[,1]), ncol=1+length(HSPS[1,]))
  colnames(output) <- c("x1", "y1", "x2", "y2", "len", "event")
  j <- 0
  for(i in 1:(length(HSPS[,1]))){
    
    if(sum(HSPS[i,]) > 0){
      j <- j + 1
      is_inverted = FALSE
      is_diagonal = TRUE
      
      if(HSPS[i,1] > HSPS[i,3]) is_inverted = TRUE
      if(abs(HSPS[i,1] - HSPS[i,2]) > DIAG_SEPARATION && abs(HSPS[i,3] - HSPS[i,4]) > DIAG_SEPARATION) is_diagonal = FALSE
      
      output[i,1] <- HSPS[i,1] * sampling
      output[i,2] <- HSPS[i,2] * sampling
      output[i,3] <- HSPS[i,3] * sampling
      output[i,4] <- HSPS[i,4] * sampling
      output[i,5] <- HSPS[i,5] * sampling
      
      if(is_diagonal) output[i,6] <- "synteny block"
      if(is_diagonal && is_inverted) output[i,6] <- "inversion"
      if(!is_diagonal && !is_inverted) output[i,6] <- "transposition"
      if(!is_diagonal && is_inverted) output[i,6] <-"inverted transposition"
    }
    
    
  }
  
  return (output[1:j,])
}



paint_frags <- function(HSPS, l, sampling){
  
  plot(c(HSPS[1,1]*sampling, HSPS[1,3]*sampling), c(HSPS[1,2]*sampling, HSPS[1,4]*sampling), xlim = c(1,l*sampling), ylim = c(1,l*sampling), type="l", xlab="X-genome",ylab="Y-genome")
  for(i in 2:length(HSPS[,1])){
    if(sum(HSPS[i,]) > 0){
      lines(c(HSPS[i,1]*sampling, HSPS[i,3]*sampling), c(HSPS[i,2]*sampling, HSPS[i,4]*sampling))
    }
  }
}

supersample <- function(mat, upscale){
  l <- min(length(mat[1,]), length(mat[,1]))
  size <- round(l*upscale)
  m <- matrix(0, nrow=size, ncol=size)
  
  for(i in 1:size){
    for(j in 1:size){
      mleft <- max(1, i-1)
      mright <- min(l, i+1)
      mup <- max(1, j-1)
      mdown <- min(l, j+1)
      
      ri <- max(1, i/upscale)
      rj <- max(1, j/upscale)
      
      if(mat[ri, rj] > 0){
        m[i, j] <- 1
      }
    }
  }
  return (m)
}


downsample <- function(mat, downscale){
  l <- min(length(mat[1,]), length(mat[,1]))
  size <- round(l/downscale)
  m <- matrix(0, nrow=size, ncol=size)
  
  for(i in 1:l){
    for(j in 1:l){
      mup <- max(1, i-1)
      mdown <- min(l, i+1)
      mleft <- max(1, j-1)
      mright <- min(l, j+1)
      
      #print(paste("Matrix at ", i, j))
      #print(mat[mup:mdown, mleft:mright])
      
      if(sum(mat[mup:mdown, mleft:mright]) > 0){
        #print(paste("goes to ", round(i/downscale), round(j/downscale)))
        m[round(i/downscale), round(j/downscale)] <- 1
      }
    }
  }
  
  return (m)
}







path_mat = args[1]


fancy_name <- strsplit(path_mat, "/")
fancy_name <- (fancy_name[[1]])[length(fancy_name[[1]])]

# Read sequence lenghts
con <- file(path_mat,"r")
seq_lengths <- readLines(con, n=2)
seq_x_len <- as.numeric(seq_lengths[1])
seq_y_len <- as.numeric(seq_lengths[2])
close(con)


data <- as.matrix(read.csv(path_mat, sep = " ", skip=2))

# Get sequence names
name_x <- strsplit(fancy_name, "-")[[1]][1]
name_y <- strsplit(fancy_name, "-")[[1]][2]

# Max of columns

len_i <- as.numeric(args[2])
len_j <- as.numeric(args[2])


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


sampling_value <- 5
submat <- downsample(score_copy, sampling_value)
m <- growing_regions((submat), WSIZE = 7, TH = 5, penalty = 15)
events <- detect_events(m, sampling_value)
events <- rbind(events, c(0,0,0,0,0,"Null event"))
write(as.character(c(seq_x_len, seq_y_len)), file=paste(path_mat,".events.txt", sep=""), append = TRUE, sep =",", ncolumns=2)
write(as.character(c("x1", "y1", "x2", "y2", "len", "event")), file=paste(path_mat,".events.txt", sep=""), append = TRUE, sep =",", ncolumns=6)
for(i in 1:length(events[,1])){
  write(as.character(events[i,]), file=paste(path_mat,".events.txt", sep=""), append = TRUE, sep =",", ncolumns=6)  
}



coords1 <- round(seq(from=0, to=1, by=0.2)*seq_x_len)
coords2 <- round(seq(from=0, to=1, by=0.2)*seq_y_len)

final_image <- apply((t(score_copy)), 2, rev)

png(paste(path_mat, ".filt.png", sep=""), width = length(data[,1]), height = length(data[,1]))
image(t(final_image), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n', main = paste(fancy_name, paste("filt. score=", score)), xlab = name_x, ylab = name_y, axes = FALSE)
axis(1, tick = TRUE, labels = (coords1), at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))
axis(2, tick = TRUE, labels = rev(coords2), at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))
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
