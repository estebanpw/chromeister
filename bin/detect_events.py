import cv2
import numpy as np
import matplotlib
import sys
from matplotlib.pyplot import imshow
from matplotlib import pyplot as plt

if(len(sys.argv) < 2):
  print("Error, use: ", sys.argv[0], " <raw matrix> [optional: plot]")
  exit()
elif(len(sys.argv) == 3 and (sys.argv[2] != "plot" and sys.argv[2] != "png")):
  print("Error, unrecognized parameter :", sys.argv[2])
  exit()


# Get length of the fastas
readheader = open(sys.argv[1], "r")
xlen = int(readheader.readline())
ylen = int(readheader.readline())
readheader.close()

# Load input matrix which must be generated with the R script compute-score (w/o grid)
matrix = np.loadtxt(sys.argv[1], skiprows = 2)
matrix = np.flip(matrix, axis=0)

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
if(len(sys.argv) == 3 and sys.argv[2] == "plot"):
  cv2.imshow("Blend between original and detection", lines_edges)
  cv2.waitKey(0)

if(len(sys.argv) == 3 and sys.argv[2] == "png"):
  namepng = sys.argv[1].replace(".raw.txt", ".events.png")
  cv2.imwrite(namepng, lines_edges) 



