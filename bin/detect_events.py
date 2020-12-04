#!/usr/bin/env python
import cv2
import numpy as np
import sys

if(len(sys.argv) < 2):
  print("Error, use: ", sys.argv[0], " <raw matrix> [plot|png]")
  exit(0)
elif(len(sys.argv) == 3 and (sys.argv[2] != "plot" and sys.argv[2] != "png")):
  print("Error, unrecognized parameter :", sys.argv[2])
  exit(0)


# Get length of the fastas
readheader = open(sys.argv[1], "r")
xlen = int(readheader.readline())
ylen = int(readheader.readline())
readheader.close()

# Load input matrix which must be generated with the R script compute-score (w/o grid)
matrix = np.loadtxt(sys.argv[1], skiprows = 2)
# Flip alongside the x axis to have same orientation as output plot
matrix = np.flip(matrix, axis=0)

# Get dimension of the matrix
dim = matrix.shape[0]

# Name output file as *.mat.events.txt
name = sys.argv[1].replace(".raw.txt", ".events.txt")

# Open output file to store events
events_file = open(name, "w")

# Write header
events_file.write("xStart,yStart,xEnd,yEnd,strand,approximate length,description\n")

# Map {0,1} to {0,255} for cv2
img = matrix * 255

# Hough transform parameters
rho = 1  # distance resolution in pixels of the Hough grid
theta = np.pi / 270  # angular resolution in radians of the Hough grid
threshold = 10  # minimum number of votes (intersections in Hough grid cell)
min_line_length = 5  # minimum number of pixels making up a line
max_line_gap = 3  # maximum gap in pixels between connectable line segments

# Classification of events parameters
DIAG_SEPARATION = int(dim * 0.015)

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
    # Classify events according to their coordinates
    x1t = int((x1 / dim) * xlen)
    x2t = int((x2 / dim) * xlen)
    # Switch y1 by y2 (because of how Hough's transform reports the direction of the lines
    y1t = int((y2 / dim) * ylen)
    y2t = int((y1 / dim) * ylen)
    # Now compare direction and classify it
    is_inverted = False
    is_diagonal = True
    if(y2t > y1t): is_inverted = True
    if(is_inverted == False and abs(x1 - y1) > DIAG_SEPARATION): is_diagonal = False # Note: diagonal separation must be calculated over the DIMxDIM grid and not over the real coordinates
    if(is_inverted == True  and abs(x1 - y2) > DIAG_SEPARATION): is_diagonal = False
    event_type = "conserved block"
    if(is_diagonal  and is_inverted):        event_type = "inversion"
    if(not is_diagonal and not is_inverted): event_type = "transposition"
    if(not is_diagonal and is_inverted):     event_type = "inverted transposition"
    # Calculate length
    mxlen = x2t - x1t
    mylen = y2t - y1t
    if(is_inverted): mylen = y1t - y2t
    # Calculate strand
    strand = "f"
    if(is_inverted): strand = "r"

    # Create a line for the plot with different color
    '''
    if(is_inverted): cv2.line(line_image,(x1,y1),(x2,y2),(0,0,255), 1)
    else:
      if(is_diagonal): cv2.line(line_image,(x1,y1),(x2,y2),(255,0,0), 1)
      else: cv2.line(line_image,(x1,y1),(x2,y2),(0,255,0), 1)
    '''
    # B, G, R (its inverted)
    # conserved blocks  = Red
    # inverted blocks   = Green
    # transposed blocks = Blue
    if(is_inverted): cv2.line(line_image,(x1,y1),(x2,y2),(0,255,0), 1)
    else:
      if(is_diagonal): cv2.line(line_image,(x1,y1),(x2,y2),(0,0,255), 1) 
      else: cv2.line(line_image,(x1,y1),(x2,y2),(255,0,0), 1)
    

    # Write to file
    events_file.write(str(x1t)+","+str(y1t)+","+str(x2t)+","+str(y2t)+","+strand+","+str(max(mxlen, mylen))+","+event_type+"\n")

# Close events since we are not writing to it anymore
events_file.close()

# Draw the lines on the  image
#lines_edges = cv2.addWeighted(img_color, 0.5, line_image, 1, 0)

# Plot the detected lines separately
#cv2.imshow("Detected lines", line_image)
#cv2.waitKey(0)

# Plot the blending between original image and detected lines
if(len(sys.argv) == 3 and sys.argv[2] == "plot"):
  #cv2.imshow("Blend between original and detection", lines_edges)
  cv2.imshow("Blend between original and detection", line_image)
  cv2.waitKey(0)

if(len(sys.argv) == 3 and sys.argv[2] == "png"):
  namepng = sys.argv[1].replace(".raw.txt", ".events.png")
  #cv2.imwrite(namepng, lines_edges) 
  cv2.imwrite(namepng, line_image) 



