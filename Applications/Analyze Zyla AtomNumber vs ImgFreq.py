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

data_folder = './cMOT freq scan_2' 
t_exp = 10e-6
picturesPerIteration = 3
sec=1
ms = 1e-3*sec

variable="img freq (MHz)"

list_file_name = variable+".txt"
list_file = os.path.join(data_folder, list_file_name)

if os.path.exists(list_file):
    freq_kHz = np.loadtxt(list_file)
else:
    List = '''
210
214
218
222
226
230
234
238
242
246
250
    '''
    freq_kHz = np.array(List.split('\n')[1:-1], dtype='float')
    np.savetxt(list_file, freq_kHz)

rowstart =0#550
rowend =-1#710
columnstart = 0#550#550
columnend = -1#800#800

params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "zyla")      
images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)

#ImageAnalysisCode.ShowImagesTranspose(images_array)

Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=None,  
                subtract_burntin=0, preventNAN_and_INF=False)

print("Number_of_atoms = "+str(Number_of_atoms/1e6))

imgNo = len(columnDensities)
angle_deg= 0 #rotates ccw
atom_numbers=[]
sizes = []
for ind in range(imgNo):
    rotated_ = rotate(columnDensities[ind][rowstart:rowend,columnstart:columnend], angle_deg, reshape = False)
    if ind==0: #first time
        rotated_columnDensities =np.zeros((imgNo, *np.shape(rotated_)))
    rotated_columnDensities[ind] = rotated_
    do_plot = True
    # if ind < 10:
    #     do_plot = False
    # else:
    #     do_plot = True
    #preview:
    dx=params.camera.pixelsize_meters/params.magnification
    popt0, popt1 = ImageAnalysisCode.fitgaussian2(rotated_columnDensities[ind],dx=dx, do_plot = do_plot, title="column density",
                                                  ylabel1D="1d density (atoms/m)", xlabel1D="distance (m)")
    if popt1 is None:
        sizes.append(np.nan)
        atom_numbers.append(np.nan)
    else:
        wy = abs(popt1[2])
        AtomNumberY = popt1[0]* wy*(2*np.pi)**0.5 
        # AtomNumberY = popt1[0]
        print("{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
        atom_numbers.append(AtomNumberY)
        sizes.append(wy)

for ind in range(imgNo):
    plt.figure()
    plt.imshow(ratio_array[ind])
    # plt.imshow(images_array[ind][0])
    plt.colorbar()
    
        
plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
plt.plot(freq_kHz, atom_numbers, 'x')
plt.xlabel(variable)
plt.ylabel("Atom Number")
plt.subplot(1,2,2)
plt.plot(freq_kHz, sizes, 'x')
plt.xlabel(variable)
plt.ylabel("cloud size (m)")
plt.tight_layout()
plt.savefig(data_folder+"/Atom number vs "+variable+ ".png")
plt.show()

#Temperature fit
# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
#                                              do_plot = True, save_folder = data_folder)

# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
#                                              do_plot = True, save_folder = data_folder)