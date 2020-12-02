import cv2
import numpy as np
import matplotlib
import sys
from matplotlib.pyplot import imshow
from matplotlib import pyplot as plt

if(len(sys.argv) != 2):
  print("Error, use: ", sys.argv[0], " <raw matrix>")
  exit()

# Get length of the fastas
readheader = open(sys.argv[1], "r")
xlen = int(readheader.readline())
ylen = int(readheader.readline())
readheader.close()

# Load input matrix which must be generated with the R script compute-score (w/o grid)
matrix = np.loadtxt(sys.argv[1], skiprows = 2)

# Get dimension of the matrix
dim = matrix.shape[0]

# Name output file as *.mat.events.txt
name = sys.argv[1].replace(".raw.txt", ".events.txt")

# Open output file to store events
events_file = open(name, "w")

# Write header
events_file.write("xStart,yStart,xEnd,yEnd,length,description\n")

# Map {0,1} to {0,255} for cv2
img = matrix * 255

# Hough transform parameters
rho = 1  # distance resolution in pixels of the Hough grid
theta = np.pi / 270  # angular resolution in radians of the Hough grid
threshold = 10  # minimum number of votes (intersections in Hough grid cell)
min_line_length = 5  # minimum number of pixels making up a line
max_line_gap = 3  # maximum gap in pixels between connectable line segments

# Create a color-image from the input one
line_image = np.uint8(np.zeros((matrix.shape[0], matrix.shape[1], 3)))
img_color = np.copy(line_image) * 0
img_color[(img > 0)] = (0,255,0)

# Plot the original image with color
#cv2.imshow("Original input with color", img_color)
#cv2.waitKey(0)

# Shortcut all smoothing and edge detection. The output matrix does not require so much post-processing
lines = cv2.HoughLinesP(np.uint8(img), rho, theta, threshold, np.array([]), min_line_length, max_line_gap)

# Create lines with color for each detected line to add them later to the plot
for line in lines:
  for x1,y1,x2,y2 in line:
    cv2.line(line_image,(x1,y1),(x2,y2),(0,0,255), 1)
    # Classify events according to their coordinates
    x1t = int((x1 / dim) * xlen)
    x2t = int((x2 / dim) * xlen)
    y1t = int((y1 / dim) * ylen)
    y2t = int((y2 / dim) * ylen)
    events_file.write(str(x1t)+","+str(y1t)+","+str(x2t)+","+str(y2t)+"\n")

# Close events since we are not writing to it anymore
events_file.close()

# Draw the lines on the  image
lines_edges = cv2.addWeighted(img_color, 0.5, line_image, 1, 0)

# Plot the detected lines separately
#cv2.imshow("Detected lines", line_image)
#cv2.waitKey(0)

# Plot the blending between original image and detected lines
cv2.imshow("Blend between original and detection", lines_edges)
cv2.waitKey(0)




'''
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

      output[i,1] <- as.numeric(round(((HSPS[i,1] * sampling) / dim) * xlen, 0))
      output[i,2] <- as.numeric(round(((HSPS[i,2] * sampling) / dim) * ylen, 0))
      output[i,3] <- as.numeric(round(((HSPS[i,3] * sampling) / dim) * xlen, 0))
      output[i,4] <- as.numeric(round(((HSPS[i,4] * sampling) / dim) * ylen, 0))


      max_len_x <- as.numeric(as.numeric(output[i,3]) - as.numeric(output[i,1]))
      max_len_y <- as.numeric(as.numeric(output[i,4]) - as.numeric(output[i,2]))
      if(HSPS[i,1] > HSPS[i,3]) {
         max_len_x <- as.numeric(as.numeric(output[i,1]) - as.numeric(output[i,3]))
      }

      output[i,5] <- max(max_len_x, max_len_y)

      if(is_diagonal) output[i,6] <- "synteny block"
      if(is_diagonal && is_inverted) output[i,6] <- "inversion"
      if(!is_diagonal && !is_inverted) output[i,6] <- "transposition"
      if(!is_diagonal && is_inverted) output[i,6] <-"inverted transposition"
    }


  }

  return (output[1:j,])
}
''' 

