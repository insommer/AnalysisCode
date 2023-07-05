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

path = r'./90%.raw'
#aw = rawpy.imread(path)
img = np.fromfile(path, dtype = np.uint16)
width = 1288
height = 964
pixelsize_um=3.75#microns
img_array = np.reshape(img, (height, width))

#plt.imshow(img_array)
#plt.show()
#plt.figure()

sum_vs_x = np.sum(img_array,0)
sum_vs_y = np.sum(img_array,1)

x0 = np.argmax(sum_vs_x)
y0 = np.argmax(sum_vs_y)

print(x0,y0)

#slices
slice_vs_x = img_array[y0,:]
slice_vs_y= img_array[:,x0]

#plt.plot(slice_vs_x)

def gaussianBeam(x, amp, center, w, offset):
    return offset + amp*np.exp(-2*(x-center)**2/w**2)

def fitgaussian(xdata, ydata,xc,do_plot=True):
    popt, pcov = curve_fit(gaussianBeam, xdata, ydata,p0=[np.max(ydata), xc, 20, 0])
    #print(popt)
    if (do_plot):
        plt.plot(xdata,ydata)
        plt.plot(xdata,gaussianBeam(xdata, *popt))
    return popt,pcov

xvalues = pixelsize_um*np.arange(width)
yvalues = pixelsize_um*np.arange(height)

popt,pcov =fitgaussian(xvalues, slice_vs_x, (x0)*pixelsize_um)
plt.figure()
popt2,pcov2 = fitgaussian(yvalues, slice_vs_y, y0*pixelsize_um)

print("X radius =",popt[2])
print("Y radius =",popt2[2])

plt.show()
plt.figure()
plt.imshow(img_array)

#x0,y0 = 461,488

