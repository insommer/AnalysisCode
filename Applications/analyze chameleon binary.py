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
date = '11/11/2024'
data_path =r"D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data"
data_path = '../Test Data'

# data_folder = r'/FLIR/Cam bias scan og position'
data_folder = r'/FLIR/Test top'

# plt.rcParams['image.interpolation'] = 'nearest'


####################################
#Parameter Setting
####################################
examNum = None #The number of runs to exam.
examFrom = None #Set to None if you want to check the last several runs. 
do_plot = True

showTimestamp = True
variablesToDisplay = None
variablesToDisplay = [
    # 'wait',
    'Lens_Position'
    # 'ODT Position',
    # 'ZSBiasCurrent',
    # 'VerticalBiasCurrent',
    # 'CamBiasCurrent'
    ]

variableFilterList = None
variableFilterList = [
    ] # NO SPACE around the operator!

rowstart = 10
rowend = -10
columnstart = 10
columnend = -10

rowstart = 300
rowend = 750
columnstart = 450
columnend = 900

# rowstart = 10
# rowend = 1524
# columnstart = 10
# columnend = 2038


binsize=1

centerx = 165 * binsize
centery = 80 * binsize
radius = 300

####################################
####################################

dataLocation = ImageAnalysisCode.GetDataLocation(date, DataPath=data_path)
data_folder = dataLocation + data_folder
variableLog_folder = dataLocation + r'/Variable Logs'
examFrom, examUntil = ImageAnalysisCode.GetExamRange(examNum, examFrom)
    
t_exp = 10e-6
picturesPerIteration = 3
# t0 = 40e-6

#config = ImageAnalysisCode.LoadConfigFile(dataFolder = data_folder)

params = ImageAnalysisCode.ExperimentParams(date, t_exp = t_exp, picturesPerIteration= picturesPerIteration, cam_type = "chameleon")      
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
                                                                                                          params=params,
                firstFrame=0, correctionFactorInput=1, rowstart=rowstart, rowend=rowend, columnstart=columnstart,
                columnend=columnend, subtract_burntin=0, preventNAN_and_INF=True)



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
    masked, vmax = ImageAnalysisCode.CircularMask(columnDensities[count], centerx=sutter_center_x, centery=sutter_center_y,
                                                      radius=radius/binsize)
    
    widthx, center_x, widthy, center_y = ImageAnalysisCode.fitgaussian(columnDensities[count], title = "Vertical Column Density",
                                                                        vmax = vmax, do_plot = 1, save_column_density=0,
                                                                        column_density_xylim=(columnstart, columnend, rowstart, rowend),
                                                                        count=count, logTime=logTime, variableLog=variableLog, 
                                                                        variablesToDisplay=variablesToDisplay, showTimestamp=True)
    center_x_array[count] = center_x
    center_y_array[count] = center_y
    print("Center x:",center_x)
    print("Center y:",center_y)
    print("Number of atoms:{}e6".format(round(Number_of_atoms[count]/(1e6))))
    if (count+1)%2 == 0:
        print("difference in center x:",center_x_array[count]-center_x_array[count - 1])
        print("difference in center y:",center_y_array[count]-center_y_array[count - 1])
