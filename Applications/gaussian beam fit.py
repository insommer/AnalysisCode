# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 15:17:30 20x22

@author: Sommer Lab
"""
#import rawpy
#import imageio
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ImageAnalysis import ImageAnalysisCode

####################################
#Set the date, the folder name, and the file name
####################################

data_path =r"D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data"
date = '5/2/2024'
data_folder = r'\FLIR\Beam Waist'
filename = r'\ODT 3400.raw'


dataLocation = ImageAnalysisCode.GetDataLocation(date,DataPath=data_path)
path = dataLocation + data_folder + filename



# date = r'/2023/12-2023/12 Dec 2023'
# data_folder = r'/'
# filename = r'/Laser focus profile between two paths.raw'

# path = data_location + date + data_folder + filename

####################################
# Choose the camera
# 1 for the old Point Grey Chameleon, 
# 1.2 for the new Point Grey Chameleon, 
# 2 for Basler dart
####################################
camera = 1.2

# path = r'./FLIR/Camera Position Check.raw'
#aw = rawpy.imread(path)

if camera == 1:
    img = np.fromfile(path, dtype = np.uint16)
    width = 1288
    height = 964
    pixelsize_um = 3.75#microns
    
if camera == 1.2:
    img = np.fromfile(path, dtype = np.uint16)
    width = 2048
    height = 1536
    pixelsize_um = 3.45#microns

elif camera == 2:
    img = np.fromfile(path, dtype = np.uint8)
    width = 3840
    height = 2160
    pixelsize_um = 2


img_array = np.reshape(img, (height, width))
rowstart = 0
rowend = -1
columnstart = 0
columnend = -1
height = height + rowend - rowstart
width = width + columnend - columnstart
img_array = img_array[rowstart:rowend, columnstart:columnend]

#plt.imshow(img_array)
#plt.show()
#plt.figure()

sum_vs_x = np.sum(img_array,0)
sum_vs_y = np.sum(img_array,1)

x0 = np.argmax(sum_vs_x)
y0 = np.argmax(sum_vs_y)

# print(x0,y0)

#slices
slice_vs_x = img_array[y0+10,:]
slice_vs_y= img_array[:,x0]

#plt.plot(slice_vs_x)

def gaussianBeam(x, amp, center, w, offset):
     return offset + amp*np.exp(-2*(x-center)**2/w**2)

def fitgaussian(xdata, ydata,xc,do_plot=True, plot_title = False):
    popt, pcov = curve_fit(gaussianBeam, xdata, ydata,p0=[np.max(ydata), xc, 100, 5e3])
    #print(popt)
    if (do_plot):
        if plot_title:
            plt.title(plot_title)
        plt.plot(xdata,ydata)
        plt.plot(xdata,gaussianBeam(xdata, *popt))
    return popt,pcov

xvalues = pixelsize_um*np.arange(width)
yvalues = pixelsize_um*np.arange(height)


# popt,pcov =fitgaussian(xvalues, sum_vs_x, (x0)*pixelsize_um, plot_title = "sum vs. x")
# plt.figure()
# popt2,pcov2 = fitgaussian(yvalues, sum_vs_y, y0*pixelsize_um, plot_title = "sum vs. y")

popt,pcov =fitgaussian(xvalues, slice_vs_x, (x0)*pixelsize_um, plot_title = "slice vs. x")
plt.figure()
popt2,pcov2 = fitgaussian(yvalues, slice_vs_y, y0*pixelsize_um, plot_title = "slice vs. y")

print("X radius = {} um".format(round(popt[2],2)))
print("Y radius = {} um".format(round(popt2[2],2)))

print("X center = {} um".format(popt[1]))
print("Y center = {} um".format(popt2[1]))

# plt.show()
plt.figure()
plt.imshow(img_array)
plt.show()

#x0,y0 = 461,488

