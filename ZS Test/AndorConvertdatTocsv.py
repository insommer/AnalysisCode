# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 11:25:41 2022

@author: Sommer Lab
"""

import matplotlib.pyplot as plt
#from google.colab import files
import csv
import itertools
import numpy as np 
import io
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import imageio
import os

#Input .dat file
filename = r"C:\Users\Sommer Lab\Documents\Data\2022\01 Jan 2022\25 Jan\scan1.dat"


#Input parameters for pictures
Pictures_Taken = 100
num_frames = 100
data_type = np.uint16

roi_mode = "full"
if roi_mode == "full":
    height = 2160
    width = 2560
    #set parameters below for region of interest, represented by red rectangle in the figure that follows
    ymax=1200
    ymin=800
    xmax = 2000
    xmin = 1400
# elif roi_mode == "cropped":
#     #TODO: put in numbers
#     height = 500
#     width = 500
#     ymax=400
#     ymin=100
#     xmax = 400
#     xmin = 100
    
elif roi_mode == "cropped":
    #TODO: put in numbers
    height = 708/biny
    width = 1148/binx
    ymax=1316/biny
    ymin=608/biny
    xmax = 2214/binx
    xmin = 1066
       
    
    
    
    
    
#ROI:
yrange = ymax-ymin
xrange = xmax - xmin

#OPEN the file and load the data
file = open(filename,"rb")
content = file.read()
data_array = np.frombuffer(content,dtype=data_type)
images = np.reshape(data_array,(num_frames,height,width))

#Preview
PreviewIndex = 2
plt.figure()
plt.imshow(images[PreviewIndex,:,:],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=5000)
# Create a Rectangle patch
rect = patches.Rectangle((xmin, ymin), xrange, yrange, linewidth=1, edgecolor='r', facecolor='none')
# Add the patch to the Axes
plt.gca().add_patch(rect)
plt.show()
plt.figure()
plt.imshow(images[PreviewIndex,ymin:ymax,xmin:xmax],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=500)
plt.show()

#Make the plot:
total_counts = np.zeros(num_frames)

for i in range(num_frames):
    im_temp = images[i,ymin:ymax,xmin:xmax]
    intensity = np.sum(im_temp)
    total_counts[i] = intensity
    
# #save the csv file    
# filename = filename.split('.dat')[0]+"_counts.csv"
# np.savetxt(filename, total_counts, delimiter = "\n")

plt.figure()
plt.xlabel("Photo Number")
plt.ylabel("Number of Counts")
plt.plot(total_counts)
# plt.switch_backend('automatic')
# xvalues = np.arange(1,101,1)
# print("xvalues = "+str(xvalues))
# print("total_counts = "+str(total_counts))
# plt.scatter(xvalues, total_counts)

plt.show()
Picture_Number = list(range(1,Pictures_Taken))

