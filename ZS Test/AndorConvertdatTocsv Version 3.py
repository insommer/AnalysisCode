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
filename = r"C:\Users\Sommer Lab\Documents\Data\2022\01 Jan 2022\21 Jan\scan1.dat"
#Region of interest for taking the counts
BordersX_pre_bin = np.array((0,0)) #define new borders by shrinking the region of interest inwards by the amounts given 
BordersY_pre_bin = np.array((0,0))


#Input the corresponding .cfg file
filename_cfg = filename.split('.dat')[0]+".cfg"
config = configparser.ConfigParser()
config.read(filename_cfg)

#Input parameters for pictures
Pictures_Taken = int(config['Acquisition']['NumberinKineticSeries'])
num_frames = Pictures_Taken
# num_frames = 100  is this variable the same as Pictures_Taken?
data_type = np.uint16


# height = int(int(config['Acquisition']['CustomAOIHeight'])/int(config['FullImage']['VerticalBin'])) #divide by bin bc it groups pixels together
# width = int(int(config['Acquisition']['CustomAOIWidth'])/int(config['FullImage']['HorizontalBin']))
ymin1 = int(config['FullImage']['VerticalStart'])
ymax1 = int(config['FullImage']['VerticalEnd'])
xmin1 = int(config['FullImage']['HorizontalStart'])
xmax1 = int(config['FullImage']['HorizontalEnd'])

binx = int(config['FullImage']['HorizontalBin'])
biny = int(config['FullImage']['VerticalBin'])

# height = int((ymax1-ymin1 +1)/int(config['FullImage']['VerticalBin']))
# width = int((xmax1-xmin1 +1)/int(config['FullImage']['HorizontalBin']))

#update borders after binning
BordersX = BordersX_pre_bin/binx
BordersY = BordersY_pre_bin/biny


xmin = int(xmin1 + BordersX_pre_bin[0]/binx)
xmax = int(xmax1 - BordersX_pre_bin[1]/binx)
ymin = int(ymin1 + BordersY_pre_bin[0]/biny)
ymax = int(ymax1 - BordersY_pre_bin[1]/biny)
height = int((ymax1-ymin1 +1-BordersY_pre_bin[0]-BordersY_pre_bin[1])/biny)
width = int((xmax1-xmin1 +1-BordersX_pre_bin[0]-BordersX_pre_bin[1])/binx)

print("height = "+str(height))
print("width = "+str(width))
print("ymin1 = "+str(ymin1))
print("ymax1 = "+str(ymax1))
print("xmin1 = "+str(xmin1))
print("xmax1 = "+str(xmax1))

print("xmin = "+str(xmin))
print("xmax = "+str(xmax))
print("ymin = "+str(ymin))
print("ymax = "+str(ymax))

# #old values
# height = 2160
# width = 2560

#these values can be adjusted to change the region of interest in the pictures for photon count collection
# ymax=1200
# ymin=800
# xmax = 2000
# xmin = 1400

#ROI:
#xmin = BordersX[0]
#xmax = 


yrange = int((ymax1-ymin1 +1)/biny)
xrange = int((xmax1-xmin1 +1)/binx)

#OPEN the file and load the data
file = open(filename,"rb")
content = file.read()
data_array = np.frombuffer(content,dtype=data_type)



# images = np.reshape(data_array,(num_frames,height,width))
# below doesn't work
# images = np.reshape(data_array,(height, num_frames, width))
# images = np.reshape(data_array,(num_frames,width, height))
# images = np.reshape(data_array,(width,height,num_frames))
# images = np.reshape(data_array,(height, width, num_frames))
images = np.reshape(data_array,(num_frames,height,width))

# print("images = "+str(images))
print("len(images)[0] = "+str(len(images[0])))
print("len(images)[0][0] = "+str(len(images[0][0])))


# This is the old code for the preview including the red rectangle over the region of interest. 
# Not necessary if we already crop the imaging region in Andor?

#Preview
PreviewIndex = 2
plt.figure()
# plt.imshow(images[PreviewIndex,:,:],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=5000)
plt.imshow(images[PreviewIndex,:,:],cmap="gray", origin="lower",interpolation="nearest")
# Create a Rectangle patch
rect = patches.Rectangle((xmin, ymin), xrange, yrange, linewidth=1, edgecolor='r', facecolor='none')
# Add the patch to the Axes
plt.gca().add_patch(rect)
plt.show()
plt.figure()
# plt.imshow(images[PreviewIndex,ymin:ymax,xmin:xmax],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=500)
plt.imshow(images[PreviewIndex,ymin:ymax,xmin:xmax],cmap="gray", origin="lower",interpolation="nearest")
plt.show()



#Preview
# plt.figure()
# # plt.imshow(images[PreviewIndex,ymin:ymax,xmin:xmax],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=500)
# plt.imshow(images[0,:,:],cmap="gray",origin = "lower", interpolation = "nearest") 
# plt.show()




#Make the plot:
total_counts = np.zeros(num_frames)

for i in range(num_frames):
    # im_temp = images[i,ymin:ymax,xmin:xmax] this is the old version
    im_temp = images[i,ymin:ymax, xmin:xmax]
    intensity = np.sum(im_temp)
    total_counts[i] = intensity
print("total_counts = "+str(total_counts))
print("max - min = "+str(max(total_counts)-min(total_counts)))    
#save the csv file    
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

