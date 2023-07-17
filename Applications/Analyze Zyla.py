# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:34:22 2023

@author: Sommer Lab
"""
from ImageAnalysis import ImageAnalysisCode
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate


data_location = r'C:/Users/Sommer Lab/Documents/Data/'

####################################
#Set the date and the folder name
####################################
date = r'/2023/07-2023/17 Jul 2023'
data_folder = r'/Andor/ODT Align'

data_folder = data_location + date + data_folder

####################################
#Parameter Setting
####################################
examNum = 5 #Only look at the last several results.
examFrom = None
do_plot = True

if examFrom is None:
    examFrom = -examNum
    
examUntil = examFrom + examNum
if examUntil == 0:
    examUntil = None

t_exp = 10e-6
picturesPerIteration = 3
ms = 1e-3

class SIUnits:
    m = 1.0
    um = 1e-6*m
units=SIUnits()

# tof_array = np.loadtxt('./ODT/Andor/TOF/TOF_list_ms.txt')*ms


rowstart = 0
rowend =-1
columnstart = 0
columnend = -1


rowstart = 100
rowend =220
columnstart = 150
columnend = 350
# rowstart = 200
# rowend =-300
# columnstart = 0
# columnend = -1

params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "zyla")      
images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)
images_array = images_array[examFrom: examUntil]

# ImageAnalysisCode.ShowImagesTranspose(images_array)

Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=1.0,  
                subtract_burntin=0, preventNAN_and_INF=True)
# plt.figure()
# plt.imshow(np.array(images_array[0][0]-images_array[0][2],dtype=np.float64)/(images_array[0][1]-images_array[0][2]),vmin=0,vmax=1.1)
# plt.imshow(images_array[0][0]-images_array[0][1])

imgNo = len(columnDensities)
angle_deg= 2 #rotates ccw

AtomNumberList=[]
widths_x = []
widths_y = []


if do_plot == True:
    fig, axs = plt.subplots(imgNo,3, figsize=(3.2*3, 2*imgNo), squeeze = False)
    plt.subplots_adjust(hspace=0.14, wspace=0.12)

for ind in range(imgNo):
    rotated_ = rotate(columnDensities[ind], angle_deg, reshape = False)[rowstart:rowend,columnstart:columnend]
    # rotated_=columnDensities[ind]
    if ind==0: #first time
        rotated_columnDensities =np.zeros((imgNo, *np.shape(rotated_)))
    rotated_columnDensities[ind] = rotated_

    #preview:
    dx=params.camera.pixelsize_meters/params.magnification    
    
    popt0, popt1 = ImageAnalysisCode.fitgaussian2D(rotated_columnDensities[ind], dx=dx, 
                                                  do_plot = do_plot, ax=axs[ind], Ind=ind, imgNo=imgNo,
                                                  title="1D density",
                                                  ylabel1D="1d density (atoms/$\mu$m)", xlabel1D="position ($\mu$m)",
                                                  title2D="column density",
                                                  xscale_factor=1/units.um, yscale_factor=units.um)
    
    # popt0, popt1 = ImageAnalysisCode.fitgaussian2(rotated_columnDensities[ind],dx=dx, do_plot = True, title="1D density",
    #                                               ylabel1D="1d density (atoms/$\mu$m)", xlabel1D="position ($\mu$m)",
    #                                               title2D="column density",
    #                                               xscale_factor=1/units.um, yscale_factor=units.um)
    
    
    wy = abs(popt1[2])
    AtomNumberY = popt1[0]* wy*(2*np.pi)**0.5 
    AtomNumberList.append(AtomNumberY)
    print("\n{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
    # print(popt1)
    
    if popt0 is not None:
        width_x = popt0[2]/units.um
        width_y = popt1[2]/units.um
        print("RMS cloud size x: {:.2f} um".format(width_x))
        print("RMS cloud size y: {:.2f} um".format(width_y))
    
        widths_x.append(width_x)
        widths_y.append(width_y)


fig.tight_layout()
print('\nThe average number of atoms:{:.2e}'.format(np.mean(AtomNumberList)))
    
print("Mean RMS width x: {:.2f} +/- {:.2f} um".format(np.mean(widths_x), np.std(widths_x)))
print("Mean RMS width y: {:.2f} +/- {:.2f} um".format(np.mean(widths_y), np.std(widths_y)))

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  

ax1.plot(widths_y, 'tab:orange')
ax1.set_ylabel('Y Widths (µm)', color='tab:orange')
ax1.tick_params(axis="y", labelcolor='tab:orange')

ax2.plot(AtomNumberList, 'tab:green')
ax2.set_ylabel('Atom Number', color='tab:green')
ax2.tick_params(axis="y", labelcolor='tab:green')


#Temperature fit
# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
#                                              do_plot = True, save_folder = data_folder)

# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
#                                              do_plot = True, save_folder = data_folder)