import sys
import os
from scipy.optimize import curve_fit
from scipy.ndimage import rotate
from scipy.integrate import simpson
 
# # getting the name of the directory
# # where the this file is present.
# current = os.path.dirname(os.path.realpath(__file__))
 
# # Getting the parent directory name
# # where the current directory is present.
# parent = os.path.dirname(current)
 
# # adding the parent directory to
# # the sys.path.
# sys.path.append(parent)
 
# # now we can import the module in the parent
# # directory.
from ImageAnalysis import ImageAnalysisCode
import numpy as np
import matplotlib.pyplot as plt

####################################
#Set the date and the folder name
####################################
date = '9/12/2023'
data_folder = r'/FLIR/odt align'

data_folder = ImageAnalysisCode.GetDataLocation(date) + data_folder

####################################
#Parameter Setting
####################################
examNum = 6 #The number of runs to exam.
examFrom = None#Set to None if you want to check the last several runs. 
do_plot = True

####################################
####################################

examFrom, examUntil = ImageAnalysisCode.GetExamRange(examNum, examFrom)
    
# data_folder =  './FLIR/odt align'
t_exp = 10e-6
picturesPerIteration = 3
# t0 = 40e-6


rowstart = 80
rowend = 225
columnstart = 130
columnend = 320

# rowstart = 0
# rowend = -1
# columnstart = 0
# columnend = -1

#config = ImageAnalysisCode.LoadConfigFile(dataFolder = data_folder)
binsize=4

params = ImageAnalysisCode.ExperimentParams(t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "chameleon")      
# images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)
# images_array = ImageAnalysisCode.loadSeriesRAW(params = params, picturesPerIteration=picturesPerIteration, data_folder = data_folder)
images_array = ImageAnalysisCode.loadSeriesPGM(picturesPerIteration=picturesPerIteration, data_folder = data_folder, 
                                               binsize=binsize, file_encoding = 'binary')
images_array = images_array[examFrom: examUntil]



# Number_of_atoms, columnDensities = ImageAnalysisCode.flsImaging(images_array, params = params, firstFrame=0, rowstart = 0, rowend = -1, 
#                                                                 columnstart =0, columnend = -1, subtract_burntin = 1)


ImageAnalysisCode.ShowImagesTranspose(images_array)


# print(data_folder+r'\\image-02242023185744-1.Raw')
# print(ImageAnalysisCode.loadRAW(data_folder+'//image-02242023185744-1.Raw'))

Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=1, rowstart = rowstart, rowend = rowend, columnstart = columnstart, columnend = columnend, 
                subtract_burntin=0, preventNAN_and_INF=True)


centerx = 920
centery = 600
radius = 250
columnDensities = ImageAnalysisCode.CircularMask(columnDensities, centerx=centerx/binsize, centery=centery/binsize,
                                                  radius=radius/binsize)

plt.figure()
plt.imshow(ratio_array[0], vmin = 0, cmap = 'gray') #vmax = 1.5
plt.colorbar()
plt.show()


# for count, x in enumerate(columnDensities):       
#     widthx, center_x, widthy, center_y = ImageAnalysisCode.fitgaussian(columnDensities[count], do_plot = 1)
#     print("Center x:",center_x)
#     print("Center y:",center_y)
# tof_array = np.loadtxt(data_folder + '/tof_list.txt')*1e-3+t0
# poptx, pcovx, popty, pcovy = ImageAnalysisCode.thermometry(params, columnDensities, tof_array, do_plot = True, data_folder = data_folder)


center_x_array = np.zeros(len(images_array))
center_y_array = np.zeros(len(images_array))

for count, x in enumerate(columnDensities):
    widthx, center_x, widthy, center_y = ImageAnalysisCode.fitgaussian(columnDensities[count],title = "Vertical Column Density",\
                                                                        vmax = None, do_plot = 1, save_column_density=0,\
                                                                            column_density_xylim=(columnstart, columnend, rowstart, rowend))
    center_x_array[count] = center_x
    center_y_array[count] = center_y
    print("Center x:",center_x)
    print("Center y:",center_y)
    print("Number of atoms:{}e6".format(round(Number_of_atoms[count]/(1e6))))
    if (count+1)%2 == 0:
        print("difference in center x:",center_x_array[count]-center_x_array[count - 1])
        print("difference in center y:",center_y_array[count]-center_y_array[count - 1])

'''
ImageAnalysisCode.imageFreqOptimization(np.loadtxt(data_folder+"/imgfreq.txt"), Number_of_atoms, ratio_array)
plt.imshow(ratio_array[0][rowstart:rowend,columnstart:columnend],vmin=0,vmax=1.2,cmap="gray")
densityvsrow = np.sum(n2d[0][rowstart:rowend,columnstart:columnend], 1)
print("densityvsrow = "+str(np.shape(densityvsrow)))
plt.figure(figsize=(4,3))
plt.plot(densityvsrow)
'''

def gaussianBeam(x, amp, center, w, offset, slope):
    return offset + amp*np.exp(-2*(x-center)**2/w**2) + slope*x


def fitgaussian(xdata, ydata, do_plot=True):
    popt, pcov = curve_fit(gaussianBeam, xdata, ydata,p0=[3e9, 670, 10, 5e9, 3e7])
    #print(popt)
    if (do_plot):
        plt.plot(xdata,ydata)
        plt.plot(xdata,gaussianBeam(xdata, *popt))
        plt.title("amplitude = {:.2e}".format(popt[0]))
        plt.ylabel("1D atomic density arb units")
        plt.xlabel("vertical row index")
        plt.tight_layout()
        # plt.savefig("atom count plot.png", dpi = 600)
    return popt,pcov

'''
#odt-specific code below
angle_deg= 0 #rotates ccw
rotated_columnDensities = rotate(columnDensities[0][rowstart:rowend,columnstart:columnend], angle_deg, reshape = False)
plt.imshow(rotated_columnDensities)
plt.colorbar()
densityvsrow = np.sum(rotated_columnDensities, 1)*deltaY
print("densityvsrow = "+str(np.shape(densityvsrow)))
plt.figure(figsize=(4,3))
plt.plot(densityvsrow)

nstart = 659
nstop = 685
xvalues = np.arange(nstart, nstop)
popt, pcov = fitgaussian(xvalues, densityvsrow[nstart:nstop], do_plot = 1)
odt_fit = gaussianBeam(xvalues, popt[0], popt[1], popt[2], popt[3], popt[4])-popt[4]*xvalues
num_atoms = simpson(odt_fit, xvalues)*deltaX
print("Number of atoms in ODT: {}e6".format(num_atoms/(1e6)))
plt.show()  
'''


