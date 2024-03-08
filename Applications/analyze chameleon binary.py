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
date = '3/7/2024'
data_path =r"D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data"
data_folder = r'/FLIR/Test'


####################################|
#Parameter Setting
####################################
examNum = 5 #The number of runs to exam.
examFrom = None #Set to None if you want to check the last several runs. 
do_plot = True

showTimestamp = True
variablesToDisplay = None
variablesToDisplay = [
    'wait',
    'ODT Position',
    'ZSBiasCurrent',
    'VerticalBiasCurrent',
    'CamBiasCurrent'
    ]
# variablesToDisplay = ['wait','cMOT coil', 'ZSBiasCurrent', 'VerticalBiasCurrent', 'CamBiasCurrent']

variableFilterList = None
variableFilterList = [
    # 'wait==3', 
    # 'VerticalBiasCurrent==9.5',
    # 'ZSBiasCurrent==5.5'
    ] # NO SPACE around the operator!


rowstart = 80
rowend = 225
columnstart = 130
columnend = 320

rowstart = 0
rowend = -1
columnstart = 0
columnend = -1

binsize=4

centerx = 165 * binsize
centery = 80 * binsize
radius = 150

####################################
####################################

dataLocation = ImageAnalysisCode.GetDataLocation(date, DataPath=data_path)
data_folder = dataLocation + data_folder
variableLog_folder = dataLocation + r'/Variable Logs'
examFrom, examUntil = ImageAnalysisCode.GetExamRange(examNum, examFrom)
    
# data_folder =  './FLIR/odt align'
t_exp = 10e-6
picturesPerIteration = 3
# t0 = 40e-6

#config = ImageAnalysisCode.LoadConfigFile(dataFolder = data_folder)


params = ImageAnalysisCode.ExperimentParams(date, t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "chameleon")      
# images_array = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder=data_folder)
# images_array = ImageAnalysisCode.loadSeriesRAW(params = params, picturesPerIteration=picturesPerIteration, data_folder = data_folder)
images_array, fileTime = ImageAnalysisCode.loadSeriesPGM(picturesPerIteration=picturesPerIteration, data_folder = data_folder, 
                                               binsize=binsize, file_encoding = 'binary', 
                                               examFrom=examFrom, examUntil=examUntil, return_fileTime=1)

variableLog = ImageAnalysisCode.LoadVariableLog(variableLog_folder)
logTime = ImageAnalysisCode.Filetime2Logtime(fileTime, variableLog)

# Number_of_atoms, columnDensities = ImageAnalysisCode.flsImaging(images_array, params = params, firstFrame=0, rowstart = 0, rowend = -1, 
#                                                                 columnstart =0, columnend = -1, subtract_burntin = 1)

if variableFilterList:        
    filterList = ImageAnalysisCode.VariableFilter(logTime, variableLog, variableFilterList)
    images_array = np.delete(images_array, filterList, 0)
    logTime = list(np.delete(logTime, filterList, 0))

ImageAnalysisCode.ShowImagesTranspose(images_array, logTime, variableLog, 
                                      variablesToDisplay, showTimestamp=showTimestamp)



Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                firstFrame=0, correctionFactorInput=1, rowstart = rowstart, rowend = rowend, columnstart = columnstart, columnend = columnend, 
                subtract_burntin=0, preventNAN_and_INF=True)


# columnDensities = ImageAnalysisCode.CircularMask(columnDensities, centerx=centerx/binsize, centery=centery/binsize,
#                                                   radius=radius/binsize)

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

for count, img in enumerate(columnDensities):
    
    # Xdistribution = images_array[count, 1].sum(axis=0)
    # Ydistribution = images_array[count, 1].sum(axis=1)
    
    sutter_widthx, sutter_center_x, sutter_widthy, sutter_center_y = ImageAnalysisCode.fitgaussian(images_array[count, 1], 
                                                                                                   title='shutter fitting', do_plot=0)
    _, vmax = ImageAnalysisCode.CircularMask(columnDensities[count], centerx=sutter_center_x, centery=sutter_center_y,
                                                      radius=radius/binsize)
    
    widthx, center_x, widthy, center_y = ImageAnalysisCode.fitgaussian(columnDensities[count],title = "Vertical Column Density",
                                                                        vmax = vmax, do_plot = 1, save_column_density=0,
                                                                        column_density_xylim=(columnstart, columnend, rowstart, rowend),
                                                                        count=count, logTime=logTime, variableLog=variableLog, 
                                                                        variablesToDisplay=variablesToDisplay, showTimestamp=False)
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