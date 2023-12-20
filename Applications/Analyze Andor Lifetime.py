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

#data_location = r'C:/Users/Sommer Lab/Documents/Data/'
data_location = r'Z:/ats317group/Data/'


####################################
#Set the date and the folder name
####################################
date = r'/2023/12-2023/12 Dec 2023'
data_folder = r'/Andor/B Field On Lens at 102 Aligned'

data_folder = data_location + date + data_folder


# data_folder = './Andor/lifetime' 
t_exp = 10e-6
subtract_bg = 1

subtract_burntin = 1
picturesPerIteration = 4 if subtract_burntin else 3

sec=1
ms = 1e-3*sec

variable="wait (s)"

list_file_name = variable+".txt"
list_file = os.path.join(data_folder, list_file_name)

if os.path.exists(list_file):
    times = np.loadtxt(list_file)
else:
    List = '''
10
30
50
100
200
300
400
500
600
10
30
50
100
200
300
400
500
600
    '''
    times = np.array(List.split('\n')[1:-1], dtype='float')
    np.savetxt(list_file, times)
    

# times = np.loadtxt(data_folder+'/wait_ms.txt')*ms



rowstart = 300
rowend = -300
columnstart=300
columnend=-300

params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "zyla")      
images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)

#ImageAnalysisCode.ShowImagesTranspose(images_array)

Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=None,  
                subtract_burntin=subtract_burntin, preventNAN_and_INF=False)

imgNo = len(columnDensities)
print(imgNo)
angle_deg= 2 #rotates ccw
#ROTATE
for ind in range(imgNo):
    rotated_ = rotate(columnDensities[ind], angle_deg, reshape = False)[rowstart:rowend,columnstart:columnend]
    if ind==0: #first time
        rotated_columnDensities =np.zeros((imgNo, *np.shape(rotated_)))
    rotated_columnDensities[ind] = rotated_

#Gauss fit for atom number:
dx=params.camera.pixelsize_meters/params.magnification
AtomNumbers = []
for ind in range(imgNo):
    poptx, popty = ImageAnalysisCode.fitgaussian2D(rotated_columnDensities[ind],dx=dx, do_plot = True, title=str(ind),
                                              ylabel1D="1d density (atoms/m)", xlabel1D="distance (microns)",
                                              subtract_bg = subtract_bg,
                                              xscale_factor=1e6)

    wy = abs(popty[2])
    AtomNumberY = popty[0]* wy*(2*np.pi)**0.5 
    AtomNumbers.append(AtomNumberY)
    print("{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
    
    
print(times)  
print(AtomNumbers)  

# AN = np.array(AtomNumbers)
# AN = np.mean(AN.reshape(-1,4), axis=1)

# WT = np.array(times)
# WT = np.mean(WT.reshape(-1,4), axis=1)

#print(WT)

plt.figure(figsize=(4,3))
plt.plot(times, AtomNumbers,'o')
plt.xlabel(variable)
plt.ylabel("Atom number")
plt.tight_layout()


ImageAnalysisCode.fit_exponential(times, AtomNumbers, offset=0, xlabel=variable, ylabel="Atom Number",legend=True)
