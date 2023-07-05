# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 11:25:41 2022

@author: Sommer Lab
"""

import matplotlib.pyplot as plt
#from google.colab import files
#import csv
#import itertools
import numpy as np 
#import io
#from matplotlib import pyplot as plt
import matplotlib.patches as patches
#import imageio
#import os
import configparser

#Input .dat file
filename = r"C:\Users\Sommer Lab\Documents\Data\2022\02 Feb 2022\8 Mar\400C.dat"
#Region of interest for taking the counts
#define new borders by shrinking the region of interest inwards by the amounts given 
BordersX_pre_bin = np.array((900,470)) 
BordersY_pre_bin = np.array((280,260))


# BordersX_pre_bin[0] += -1
# BordersY_pre_bin[0] += -1
#Input the corresponding .cfg file
filename_cfg = filename.split('.dat')[0]+".cfg"
#filename_cfg = "scan0.cfg"
config = configparser.ConfigParser()
config.read(filename_cfg)

#Input parameters for pictures
Pictures_Taken = int(config['Acquisition']['NumberinKineticSeries'])
num_frames = Pictures_Taken
# num_frames = 100  is this variable the same as Pictures_Taken?
data_type = np.uint16


#before binning
height1 = int(config['FullImage']['VerticalEnd']) - int(config['FullImage']['VerticalStart']) + 1
width1 = int(config['FullImage']['HorizontalEnd']) - int(config['FullImage']['HorizontalStart']) + 1
ymin1 = int(0) #origin placed at zero by python
ymax1 = height1 
xmin1 = int(0)
xmax1 = width1  

#bin values
binx = int(config['FullImage']['HorizontalBin'])
biny = int(config['FullImage']['VerticalBin'])

#redefine x and y with borders
xmin = int(int(xmin1 + BordersX_pre_bin[0])/binx)
xmax = int(int(xmax1 - BordersX_pre_bin[1])/binx)
ymin = int(int(ymin1 + BordersY_pre_bin[0])/biny)
ymax = int(int(ymax1 - BordersY_pre_bin[1])/biny)

#height and width after binning
height = int(height1/biny)
width = int(width1/binx)


yrange = ymax - ymin
xrange = xmax - xmin


#OPEN the file and load the data
file = open(filename,"rb")
content = file.read()
data_array = np.frombuffer(content,dtype=data_type)
images = np.reshape(data_array,(num_frames,height,width))

#Preview
PreviewIndex = 2
plt.figure()
# plt.imshow(images[PreviewIndex,:,:],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=5000)
#old settings: vmin=150,vmax=35000, vmin and vmax set the apparent brightness of the displayed image, but do not affect the counts file/plot
vmin = 1500 
vmax = 10000
plt.imshow(images[PreviewIndex,:,:],cmap="gray", origin="lower",interpolation="nearest",vmin=vmin,vmax=vmax)
# Create a Rectangle patch, use xmin-.5, ymin-.5 because origin = "lower" sets the left border at x = -.5 and the bottom border at y = -.5
rect = patches.Rectangle((xmin-.5, ymin-.5), xrange, yrange, linewidth=1, edgecolor='r', facecolor='none')
# Add the patch to the Axes
plt.gca().add_patch(rect)
plt.show()
plt.figure()
# plt.imshow(images[PreviewIndex,ymin:ymax,xmin:xmax],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=500)
plt.imshow(images[PreviewIndex,ymin:ymax,xmin:xmax],cmap="gray", origin="lower",interpolation="nearest",vmin=vmin,vmax=vmax)
plt.show()


#Make the plot:
total_counts = np.zeros(num_frames)

for i in range(num_frames):
    # im_temp = images[i,ymin:ymax,xmin:xmax] this is the old version
    im_temp = images[i,ymin:ymax, xmin:xmax]
    intensity = np.sum(im_temp)
    total_counts[i] = intensity
 
#save the csv file    
# filename = filename.split('.dat')[0]+"_counts.csv"
# np.savetxt(filename, total_counts, delimiter = "\n")

plt.figure()
plt.xlabel("Photo Number")
plt.ylabel("Number of Counts")
plt.plot(total_counts)

plt.show()
Picture_Number = list(range(1,Pictures_Taken))

