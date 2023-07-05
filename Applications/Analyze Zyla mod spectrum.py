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

data_folder = './Andor/trap mod 100ms wait tmod 50ms amplitude 0.5V_1' 
t_exp = 10e-6
picturesPerIteration = 3
sec=1
ms = 1e-3*sec

if os.path.exists(data_folder+"/freq_kHz.txt"):
    freq_kHz = np.loadtxt(data_folder+"/freq_kHz.txt")
else:
    List = '''
2
2.4
2.8
3.2
3.6
4
4.4
4.8
5.2
5.6
6
6.4
6.8
7.2
7.6
8
8.4
8.8
9.2
9.6
10
10.4
10.8
11.2
11.6
12
12.4
12.8
13.2
13.6
14
14.4
14.8
15.2
15.6
16
    '''
    freq_kHz = np.array(List.split('\n')[1:-1], dtype='float')
    np.savetxt(data_folder+"/freq_kHz.txt", freq_kHz, fmt='%.2f')

rowstart = 100
rowend =220
columnstart = 175
columnend = 375


# rowstart = 0
# rowend = -1
# columnstart = 0
# columnend = -1#300

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
    popt0, popt1 = ImageAnalysisCode.fitgaussian2(rotated_columnDensities[ind],dx=dx, do_plot = do_plot, title="",
                                                  ylabel1D="1d density (atoms/m)", xlabel1D="distance (m)")
    
    if popt1 is None:
        sizes.append(np.nan)
        atom_numbers.append(np.nan)
    else:
        wy = abs(popt1[2])
        # AtomNumberY = popt1[0]* wy*(2*np.pi)**0.5 
        AtomNumberY = popt1[0]
        print("{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
        atom_numbers.append(AtomNumberY)
        sizes.append(wy)
plt.figure(figsize=(7,3))
plt.subplot(1,2,1)
plt.plot(freq_kHz, atom_numbers, 'x')
plt.xlabel("freq (kHz)")
plt.ylabel("Amplitude")
plt.subplot(1,2,2)
plt.plot(freq_kHz, sizes, 'x')
plt.xlabel("freq (kHz)")
plt.ylabel("cloud size (m)")
plt.tight_layout()
plt.savefig(data_folder+"/spectrum.png")
plt.show()

#Temperature fit
# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
#                                              do_plot = True, save_folder = data_folder)

# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
#                                              do_plot = True, save_folder = data_folder)