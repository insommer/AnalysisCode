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
date = r'/2023/07-2023/18 Jul 2023'
data_folder = r'/Andor/trap mod 100ms wait tmod 50ms amp 0.5V_2'

data_folder = data_location + date + data_folder

t_exp = 10e-6
picturesPerIteration = 3
sec=1
ms = 1e-3*sec

if os.path.exists(data_folder+"/freq_kHz.txt"):
    freq_kHz = np.loadtxt(data_folder+"/freq_kHz.txt")
else:
    List = '''
15.2
15.6
10.8
11.2
10.4
4
6.4
8.4
11.6
14.8
2
10
3.2
7.6
2.4
12
9.2
8
6
11.6
6.4
6.8
10.4
12.8
8.8
5.6
14
10
10.8
14.8
2.8
8.4
8.8
12.4
15.2
4.4
7.2
9.2
14
12.4
12.8
13.6
14.4
16
9.6
13.2
4.8
8
14.4
9.6
2
4.4
2.8
2.4
12
6.8
5.2
16
3.6
4
5.6
11.2
4.8
15.6
13.6
5.2
3.6
7.6
6
13.2
3.2
7.2
    '''
    freq_kHz = np.array(List.split('\n')[1:-1], dtype='float')
    np.savetxt(data_folder+"/freq_kHz.txt", freq_kHz)

rowstart = 150
rowend = 320
columnstart = 160
columnend = 420


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
        do_plot = False
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


freq_kz_unique = np.unique(freq_kHz)

atom_numbers_array = np.array(atom_numbers)
atom_numbers_avg = [ atom_numbers_array[freq_kHz==ii].mean() for ii in freq_kz_unique ]
atom_numbers_std = [ atom_numbers_array[freq_kHz==ii].std() for ii in freq_kz_unique ]

cloudsize_array = np.array(sizes)
cloudsize_avg = [ np.nanmean(cloudsize_array[freq_kHz==ii]) for ii in freq_kz_unique ]
cloudsize_std = [ np.nanstd(cloudsize_array[freq_kHz==ii]) for ii in freq_kz_unique ]

fig, axs = plt.subplots(1, 2, figsize=(7,3))
axs[0].errorbar(freq_kz_unique, atom_numbers_avg, atom_numbers_std)
axs[0].set_xlabel("freq (kHz)")
axs[0].set_ylabel("Amplitude")

axs[1].errorbar(freq_kz_unique, cloudsize_avg, cloudsize_std)
axs[1].set_xlabel("freq (kHz)")
axs[1].set_xlabel("cloud size (m)")

plt.tight_layout()
plt.savefig(data_folder+"/spectrum.png")
plt.show()



#Temperature fit
# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
#                                              do_plot = True, save_folder = data_folder)

# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
#                                              do_plot = True, save_folder = data_folder)