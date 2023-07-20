# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:34:22 2023

@author: Sommer Lab
"""
from ImageAnalysis import ImageAnalysisCode
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import os

data_location = r'C:/Users/Sommer Lab/Documents/Data/'

####################################
#Set the date and the folder name
####################################
date = r'/2023/07-2023/20 Jul 2023'
data_folder = r'/Andor/Test TOF'

data_folder = data_location + date + data_folder

# data_folder = './lens 54mm wait 200ms TOF' 
t_exp = 10e-6
picturesPerIteration = 3
sec=1
ms = 1e-3*sec
# tof_array = np.loadtxt(data_folder+"/TOF_ms.txt")*ms


variable="TOF (ms)"

list_file_name = variable+".txt"
list_file = os.path.join(data_folder, list_file_name)

if os.path.exists(list_file):
    tof_array = np.loadtxt(list_file)*ms
    print("Loaded")
else:
    List = '''
0
0.5
1
1.5
2
2.5
3
3.5
4
    '''
    tof_array = np.array(List.split('\n')[1:-1], dtype='float')
    np.savetxt(list_file, tof_array)
    tof_array = tof_array * ms


rowstart = 10
rowend =-10
columnstart = 10
columnend = -50

params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "zyla")      
images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)

#ImageAnalysisCode.ShowImagesTranspose(images_array)

Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=None,  
                subtract_burntin=0, preventNAN_and_INF=False)

imgNo = len(columnDensities)
angle_deg= 1 #rotates ccw

for ind in range(imgNo):
    rotated_ = rotate(columnDensities[ind][rowstart:rowend,columnstart:columnend], angle_deg, reshape = False)
    if ind==0: #first time
        rotated_columnDensities =np.zeros((imgNo, *np.shape(rotated_)))
    rotated_columnDensities[ind] = rotated_

#preview:
dx=params.camera.pixelsize_meters/params.magnification
# popt0, popt1 = ImageAnalysisCode.fitgaussian2(rotated_columnDensities[3],dx=dx, do_plot = True, title="",
#                                               ylabel1D="1d density (atoms/m)", xlabel1D="distance (m)")

#Temperature fit
# plt.figure(figsize=(8,3))
# plt.subplot(1,2,1)
popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
                                             do_plot = True, save_folder = data_folder, newfig=True)
# plt.subplot(1,2,2)
popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
                                             do_plot = True, save_folder = data_folder, newfig=True)
plt.tight_layout()