# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 10:28:06 2022

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
import scipy.optimize
import math

print("Input scan number.")
scan_num = str(input())
#Input .dat file
filename = r"C:\Users\Sommer Lab\Documents\Data\2022\04 Apr 2022\14 Apr\scan"+scan_num+".dat"
#Region of interest for taking the counts
#define new borders by shrinking the region of interest inwards by the amounts given 
BordersX_pre_bin = np.array((0,0)) 
BordersY_pre_bin = np.array((0,0))


# BordersX_pre_bin[0] += -1
# BordersY_pre_bin[0] += -1
#Input the corresponding .cfg file
filename_cfg = filename.split('scan'+scan_num+'.dat')[0]+"config.cfg"
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



#Preview
# plt.figure()
# plt.imshow(images[PreviewIndex,ymin:ymax,xmin:xmax],cmap="gray", origin="lower",interpolation="nearest",vmin=150,vmax=500)
# plt.imshow(images[0,:,:],cmap="gray",origin = "lower", interpolation = "nearest") 
# plt.show()




#Make the plot:
total_counts = np.zeros(num_frames)

for i in range(num_frames):
    # im_temp = images[i,ymin:ymax,xmin:xmax] this is the old version
    im_temp = images[i,ymin:ymax, xmin:xmax]
    intensity = np.sum(im_temp)
    total_counts[i] = intensity
 
#save the csv file    
#filename = filename.split('scan'+scan_num+'.dat')[0]+"Atom_count.csv"
#np.savetxt(filename, total_counts, delimiter = "\n")

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
Picture_Number = np.array(list(range(0,Pictures_Taken)))




#Counts to Atoms#
N_counts = total_counts #total counts in MOT
P = 29 #MOT beam power in mW
t_exp = float(config['Acquisition']['ExposureTime']) #s
Picture_Time = Picture_Number*float(config['Acquisition']['KineticCycleTime'])

#dont change parameters here
quantum_eff = .5
r_beam = .082
I1 = 2*P/(np.pi*r_beam**2)
I = 6*I1
Isat = 25.4 #W/m^2
s = I/Isat
gamma = 36.898e6
delta = 26e6*2*np.pi
R_scat = gamma*.5*s/(1+s+(2*delta/gamma)**2)
cos_theta = 150/np.sqrt(3.175**2+150**2)
solid_angle = 2*np.pi*(1-cos_theta)
N_atoms = (4*np.pi*N_counts)/(quantum_eff*R_scat*t_exp*solid_angle)

plt.figure()
plt.xlabel("Time elapsed (s)")
plt.ylabel("Number of Atoms")
plt.plot(Picture_Time, N_atoms, 'o')

atoms_vs_time = np.stack((Picture_Time, N_atoms))
#print(atoms_vs_time)

ready_to_save = 'false'
if ready_to_save == 'true':
    csvfilename = filename.split('scan'+scan_num+'.dat')[0]+"Atom_count.csv"
    np.savetxt(csvfilename, atoms_vs_time, delimiter = ",") 
    


print("Is this a load or decay curve? Type 'Load' or 'Decay'.")
Load_or_Decay = str(input())

#Fit parameters
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

atom_max = max(N_atoms)
array = N_atoms
value = atom_max*np.exp(-1)
#print("atom_max: ",atom_max)
#print("atom_max/e:",atom_max/np.exp(1))
#print("value: ",value)
emin1 = find_nearest(array, value)
#print("find nearest: ",find_nearest(array, value))
#print("emin1: ", emin1)
finder = np.where(N_atoms == emin1)
#print("where: ", finder)
#print(finder[0][0])
array_number = int(finder[0])
#print("array_number: ", array_number)
#######################################This is the time for the function to reach e**-1 of max value
emin1_time = Picture_Time[array_number]
#print("emin1_time: ",emin1_time)

#Fit of data
if Load_or_Decay == 'Decay':
    def Load_Decay(x, m, t, b):
         return m * np.exp(-t * x) + b
    p0 = (atom_max, 1/emin1_time, 0) # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(Load_Decay, Picture_Time, N_atoms, p0)
    m, t, b = params
    #Quality of fit
    squaredDiffs = np.square(N_atoms - Load_Decay(Picture_Time, m, t, b))
    squaredDiffsFromMean = np.square(N_atoms - np.mean(N_atoms))
    rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
    print(f"R² = {rSquared}")
    print(f"Y = {m} * e^(-{t} * x) + {b}")
    # plot the results
    plt.plot(Picture_Time, N_atoms, '.', label="data")
    plt.plot(Picture_Time, Load_Decay(Picture_Time, m, t, b), '--', label="fitted")
    plt.title("Fitted Exponential Curve")
    
if Load_or_Decay == 'Load':
    def Load_Decay(x, m, t, b):
         return -m * np.exp(-t * x) + b
    p0 = (atom_max, (1-math.log(math.e-1))/emin1_time, atom_max) # start with values near those we expect
    params, cv = scipy.optimize.curve_fit(Load_Decay, Picture_Time, N_atoms, p0)
    m, t, b = params
    #Quality of fit
    squaredDiffs = np.square(N_atoms - Load_Decay(Picture_Time, m, t, b))
    squaredDiffsFromMean = np.square(N_atoms - np.mean(N_atoms))
    rSquared = 1 - np.sum(squaredDiffs) / np.sum(squaredDiffsFromMean)
    print(f"R² = {rSquared}")
    print(f"Y = {m} * e^(-{t} * x) + {b}")
    # plot the results
    plt.plot(Picture_Time, N_atoms, '.', label="data")
    plt.plot(Picture_Time, Load_Decay(Picture_Time, m, t, b), '--', label="fitted")
    plt.title("Fitted Exponential Curve")

