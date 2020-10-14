#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 2){
  stop("USE: Rscript --vanilla plot.R <matrix> <matsize>")
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



#data <- as.matrix(read.csv(path_mat, sep = " ", skip=2))
data <- as.matrix(read.table(path_mat, skip=2, sep=" "))

# Get sequence names
name_x <- strsplit(fancy_name, "-")[[1]][1]
name_y <- strsplit(fancy_name, "-")[[1]][2]

# Max of columns

len_i <- as.numeric(args[2])
len_j <- as.numeric(args[2])


score_density <- data
aux_density <- data


# To compute the score
score <- 1
#print(score)
cat(score, "\n")


sampling_value <- 5
events <- c()
events <- rbind(events, c(0,0,0,0,0,"Null event"))
write(as.character(c(seq_x_len, seq_y_len)), file=paste(path_mat,".events.txt", sep=""), append = TRUE, sep =",", ncolumns=2)
write(as.character(c("x1", "y1", "x2", "y2", "len", "event")), file=paste(path_mat,".events.txt", sep=""), append = TRUE, sep =",", ncolumns=6)
for(i in 1:length(events[,1])){
  write(as.character(events[i,]), file=paste(path_mat,".events.txt", sep=""), append = TRUE, sep =",", ncolumns=6)  
}



coords1 <- round(seq(from=0, to=1, by=0.2)*seq_x_len)
coords2 <- round(seq(from=0, to=1, by=0.2)*seq_y_len)


score_copy <- data

#print(score_copy)


themax <- max(score_copy)

for(i in 1:(length(score_copy[,1]))){
  for(j in 1:(length(score_copy[1,]))){
	myv <- score_copy[i,j]

	for(k in seq(0.1,0.9,0.1)){
		realk <- 1-k
		if(myv > realk * themax){
			#print(paste(myv, "is bigger than", realk*themax, "so i assign", realk*10))
			score_copy[i,j] <- realk*10
			break
		}
	}

	#if(myv > 0.9 * themax) score_copy[i,j] <- 6
	#if(myv > 0.9) score_copy[i,j] <- 4
	#if(myv > 30) score_copy[i,j] <- 3
	#if(myv > 20) score_copy[i,j] <- 2
	#if(myv > 0) score_copy[i,j] <- 1

  }
}


final_image <- apply((t(score_copy)), 2, rev)

png(paste(path_mat, ".filt.png", sep=""), width = length(data[,1]), height = length(data[,1]))
image(t(final_image), col = grey(seq(1, 0, length = 256)), xaxt='n', yaxt='n', main = paste(fancy_name, paste("filt. score=", "unknown")), xlab = name_x, ylab = name_y, axes = FALSE)

# col = grey(seq(1, 0, length = 256))

axis(1, tick = TRUE, labels = (coords1), at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))
axis(2, tick = TRUE, labels = rev(coords2), at = c(0.0, 0.2, 0.4, 0.6, 0.8, 1))
tmp <- dev.off()


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
