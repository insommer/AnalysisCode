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
date = r'/2023/07-2023/31 Jul 2023'
data_folder = r'/Andor/RF_spec_225.5-228.3'

data_folder = data_location + date + data_folder

t_exp = 10e-6
picturesPerIteration = 3
sec=1
ms = 1e-3*sec

variable="freq (MHz)"

list_file_name = variable+".txt"
list_file = os.path.join(data_folder, list_file_name)

if os.path.exists(list_file):
    cicero_values = np.loadtxt(list_file)
else:
    List = '''
225.5
225.6
225.7
225.8
225.9
226
226.1
226.2
226.3
226.4
226.5
226.6
226.7
226.8
226.9
227
227.1
227.2
227.3
227.4
227.5
227.6
227.7
227.8
227.9
228
228.1
228.2
228.3
    '''
    cicero_values = np.array(List.split('\n')[1:-1], dtype='float')
    np.savetxt(list_file, cicero_values)

rowstart =0#550
rowend =-200#710
columnstart = 0#550#550
columnend = -500#800#800

params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "zyla")      
images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)

#ImageAnalysisCode.ShowImagesTranspose(images_array)

Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=None,  
                subtract_burntin=0, preventNAN_and_INF=False)

print("Number_of_atoms = "+str(Number_of_atoms/1e6))

imgNo = len(columnDensities)
angle_deg= 0 #rotates ccw
atom_numbers_x=[]
atom_numbers_y=[]
amplitudes_x=[]
amplitudes_y=[]
atom_counts =[]
sizes = []
NUM_ROWS=4 #rows of plots per figure
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
    row = ind % NUM_ROWS #current row within the figure
    new_figure = (row == 0) 
    try:
        popt0, popt1 = ImageAnalysisCode.fitgaussian2(rotated_columnDensities[ind],dx=dx, do_plot = do_plot, title="column density",
                                                  ylabel1D="1d density (atoms/m)", xlabel1D="distance (m)",new_figure=new_figure,
                                                  num_rows=NUM_ROWS, row=row)
    except:
        popt0, popt1 = None, None
        
    if popt1 is None or popt0 is None or popt0[0]>1e10 or popt1[0]>1e10:
        sizes.append(np.nan)
        atom_numbers_x.append(np.nan)
        atom_numbers_y.append(np.nan)
        amplitudes_x.append(np.nan)
        amplitudes_y.append(np.nan)
    else:
        #popt0 is x fit
        #popt1 is y fit
        fit_w_x = abs(popt0[2])
        fit_w_y = abs(popt1[2])
        AtomNumberX = popt0[0]* fit_w_x*(2*np.pi)**0.5 
        AtomNumberY = popt1[0]* fit_w_y*(2*np.pi)**0.5 
        # AtomNumberY = popt1[0]
        print("{}. Atom Number from gauss fit (X) = {:.2e}".format(ind, AtomNumberX))
        atom_numbers_x.append(AtomNumberX)
        atom_numbers_y.append(AtomNumberY)
        amplitudes_x.append(popt0[0])
        amplitudes_y.append(popt1[0])
        sizes.append(fit_w_x)
    count = np.sum(rotated_)*dx*dx
    atom_counts.append(count)

# for ind in range(imgNo):
#     plt.figure()
#     plt.imshow(ratio_array[ind])
#     # plt.imshow(images_array[ind][0])
#     plt.colorbar()
    
        
plt.figure(figsize=(8,4))
plt.subplot(1,2,1)
plt.plot(cicero_values, atom_numbers_x, 'x',label="x")
plt.plot(cicero_values, atom_numbers_y, 'o',label="y")
# plt.plot(cicero_values, atom_counts, '.-',label="count")
plt.xlabel(variable)
plt.ylabel("Atom Number")
plt.subplot(1,2,2)
plt.plot(cicero_values, amplitudes_x, 'x',label="x")
plt.plot(cicero_values, amplitudes_y, 'o',label="y")
plt.xlabel(variable)
plt.ylabel("fit amplitude")
plt.legend()
plt.tight_layout()
plt.savefig(data_folder+"/Atom number vs "+variable+ ".png")
plt.show()

#Temperature fit
# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
#                                              do_plot = True, save_folder = data_folder)

# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
#                                              do_plot = True, save_folder = data_folder)