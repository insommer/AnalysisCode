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
date = r'/2023/09-2023/07 Sep 2023'
data_foldername = 'Evap Time (ms)_High PD 2.45 V Unchanged Evap Time 2 50 ms'

data_folder = './Andor/'+ data_foldername
data_folder = data_location + date + data_folder


t_exp = 10e-6
picturesPerIteration = 3
sec=1
ms = 1e-3*sec

variable = data_foldername.split('_')[0]
# variable = 'wait (ms)'

list_file_name = variable+".txt"
list_file = os.path.join(data_folder, list_file_name)

if os.path.exists(list_file):
    variable_array = np.loadtxt(list_file)
else:
    List = '''
50
100
350
300
450
100
200
350
250
300
450
300
50
250
400
400
450
100
50
250
150
150
500
400
500
150
200
500
200
350
'''
    variable_array = np.array(List.split('\n')[1:-1], dtype='float')
    np.savetxt(list_file, variable_array)

rowstart = 0
rowend = -1
columnstart = 0
columnend = -1#300

rowstart =100
rowend = -130
columnstart = 100
columnend = -100

params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "zyla")      
images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)

#ImageAnalysisCode.ShowImagesTranspose(images_array)

Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=None,  
                subtract_burntin=0, preventNAN_and_INF=False)

imgNo = len(columnDensities)
angle_deg= 1 #rotates ccw
atom_numbers=[]
sizes = []
ycenter = []
for ind in range(imgNo):
    rotated_ = rotate(columnDensities[ind][rowstart:rowend,columnstart:columnend], angle_deg, reshape = False)
    if ind==0: #first time
        rotated_columnDensities =np.zeros((imgNo, *np.shape(rotated_)))
    rotated_columnDensities[ind] = rotated_
    if ind < 10:
        do_plot = True
    else:
        do_plot = True
    #preview:
    dx=params.camera.pixelsize_meters/params.magnification
    popt0, popt1 = ImageAnalysisCode.fitgaussian2(rotated_columnDensities[ind],dx=dx, do_plot = 1, title="",
                                                  ylabel1D="1d density (atoms/m)", xlabel1D="distance (m)")
    
    if popt1 is None:
        sizes.append(np.nan)
        atom_numbers.append(np.nan)
        ycenter.append(np.nan)        
    else:
        wy = abs(popt1[2])
        AtomNumberY = abs(popt1[0]* wy*(2*np.pi)**0.5)
        # AtomNumberY = popt1[0]
        print("{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
        atom_numbers.append(AtomNumberY)
        sizes.append(wy)
        ycenter.append(popt1[1])

#Calculate the average and std for the results.        
variable_unique = np.unique(variable_array)

sizes_array = np.array(sizes)
sizes_avg = [ sizes_array[variable_array==ii].mean() for ii in variable_unique ]
sizes_std = [ sizes_array[variable_array==ii].std() for ii in variable_unique ]

atom_numbers_array = np.array(atom_numbers)
atom_numbers_avg = [ atom_numbers_array[variable_array==ii].mean() for ii in variable_unique ]
atom_numbers_std = [ atom_numbers_array[variable_array==ii].std() for ii in variable_unique ]

ycenter_array = np.array(ycenter)
ycenter_avg = [ ycenter_array[variable_array==ii].mean() for ii in variable_unique ]
ycenter_std = [ ycenter_array[variable_array==ii].std() for ii in variable_unique ]

#plotting
plt.figure(figsize=(15,4))
plt.subplot(1,3,1)
plt.errorbar(variable_unique, atom_numbers_avg, atom_numbers_std)
# plt.plot(variable_unique, atom_numbers_avg)
# plt.plot(variable_array, atom_numbers, 'x')

plt.xlabel(variable)
plt.ylabel("Atom Number")

plt.subplot(1,3,2)
plt.errorbar(variable_unique, sizes_avg, sizes_std)
# plt.plot(variable_unique, sizes_avg)
# plt.plot(variable_array, sizes, 'x')
plt.xlabel(variable)
plt.ylabel("cloud size (m)")

plt.subplot(1,3,3)
plt.errorbar(variable_unique, ycenter_avg, ycenter_std)
# plt.plot(variable_unique, ycenter_avg)
# plt.plot(variable_array, ycenter, 'x')
plt.xlabel(variable)
plt.ylabel("ycenter (m)")

plt.tight_layout()
plt.savefig(data_folder+"/Amplitude vs "+variable+ ".png")
plt.show()

#Temperature fit
# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
#                                              do_plot = True, save_folder = data_folder)

# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
#                                              do_plot = True, save_folder = data_folder)