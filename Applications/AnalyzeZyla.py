# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:34:22 2023

@author: Sommer Lab
"""
DEBUG=True
def log(msg,indent=0):
    if DEBUG:
        print("#########  "+'\t'*indent+ msg)
log("Starting script Analyze Zyla...")
from ImageAnalysis import ImageAnalysisCode
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import pandas as pd
import os

####################################
#Set the date and the folder name
####################################
data_path =r"Z:\ats317group\Data"
date = '10/10/2023'
# data_folder = [r'/Andor/MOT and ODT Movement', r'/Andor/MOT and ODT Movement_3',
#                r'/Andor/MOT and ODT Movement data_1', r'/Andor/MOT and ODT Movement data', 
#                r'/Andor/MOT and ODT Movement data_2', r'/Andor/MOT and ODT Movement data_4',
#                r'/Andor/ODT Position 5']
data_folder1 = [r'/Andor/Position 1 Bias Scan']
# data_folder1 = [r'/spool test']
####################################
#Parameter Setting
####################################
# repetition = 1 #The number of identical runs to be averaged. 
# examNum = 99 #The number of runs to exam.
# examFrom1 = None #Set to None if you want to check the last several runs. 
# plotPWindow = 3
# do_plot = True
# uniformscale = 0

variablesToDisplay = ['wait','cMOT coil', 'VerticalBiasCurrent', 'ZSBiasCurrent'] 
# showTimestamp = False

# variableFilterList = None
# variableFilterList = ['wait==30']#, 'VerticalBiasCurrent==3'] # NO SPACE around the operator!

# pictureToHide = None
# # pictureToHide = [8] # list(range(0,10,2))

# subtract_bg = True
# signal_feature = 'narrow' 


####################################
####################################
class SIUnits:
    m = 1.0
    um = 1e-6*m
units=SIUnits()

def AnalyzeZyla(date, 
                data_folders, 
                data_path='C:/Users/Sommer Lab/Documents/Data/', 
                repetition=1, 
                examNum=3, 
                examFrom=None, 
                plotPWindow=5, do_plot=True, uniformscale=0, 
                variablesToDisplay= variablesToDisplay,
                showTimestamp=False, variableFilterList=None, 
                pictureToHide=None,
                subtract_bg=True, 
                signal_feature='narrow', 
                rowstart=0, rowend=-1, 
                columnstart=0, columnend=-1,
                picturesPerIteration = 3):
    
    dataLocation = ImageAnalysisCode.GetDataLocation(date,DataPath=data_path)
    data_folders = [ dataLocation + f for f in data_folders ]
    variableLog_folder = dataLocation + r'/Variable Logs'
    examFrom, examUntil = ImageAnalysisCode.GetExamRange(examNum, examFrom, repetition)
    
    
    params = ImageAnalysisCode.ExperimentParams(t_exp = 10e-6, picturesPerIteration= picturesPerIteration, cam_type = "zyla")
    images_array = None
    log("Loading spooling series...")
    for ff in data_folders:
        log("loading "+os.path.split(ff)[-1],1)
        if images_array is None:
            images_array, fileTime = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder = ff, 
                                                                       return_fileTime=1)
        else:
            _images_array, _fileTime = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder = ff, 
                                                                           return_fileTime=1)
            images_array = np.concatenate([images_array, _images_array], axis=0)
            fileTime = fileTime + _fileTime
    
    images_array = images_array[examFrom: examUntil]
    fileTime = fileTime[examFrom: examUntil]
    
    variableLog = ImageAnalysisCode.LoadVariableLog(variableLog_folder)
    logTime = ImageAnalysisCode.Filetime2Logtime(fileTime, variableLog)
        
    if variableFilterList is not None and variableLog is not None:    
        filteredList = ImageAnalysisCode.VariableFilter(logTime, variableLog, variableFilterList)
        images_array = np.delete(images_array, filteredList, 0)
        logTime = np.delete(logTime, filteredList, 0)
    
    if pictureToHide is not None:
        images_array = np.delete(images_array, pictureToHide, 0)
        if logTime is not None:
            logTime = np.delete(logTime, pictureToHide, 0)
    
    # ImageAnalysisCode.ShowImagesTranspose(images_array)
    log("Calculating column densities...")
    Number_of_atoms, N_abs, ratio_array, columnDensities, deltaX, deltaY = ImageAnalysisCode.absImagingSimple(images_array, 
                    firstFrame=0, correctionFactorInput=1.0,  
                    subtract_burntin=0, preventNAN_and_INF=True)
    # plt.figure()
    # plt.imshow(np.array(images_array[0][0]-images_array[0][2],dtype=np.float64)/(images_array[0][1]-images_array[0][2]),vmin=0,vmax=1.1)
    # plt.imshow(images_array[0][0]-images_array[0][1])
    
    
    
    imgNo = len(columnDensities)
    angle_deg= 2 #rotates ccw
    
    AtomNumbers=[]
    widths_x = []
    widths_y = []
    centers_x = []
    centers_y = []
    variableDicts = variableLog.to_dict('index')
    dataDicts = []
    
    if uniformscale:
        vmax = columnDensities.max()
        vmin = columnDensities.min()
    else:
        vmax = None
        vmin = 0
    log("Running gaussian fits...")
    for ind in range(imgNo):
        log("gauss fit for index {}".format(ind))
        
        currentDataDict = variableDicts[logTime[ind]]
        
        plotInd = ind % plotPWindow
        if do_plot == True and plotInd == 0:
            # if ind//plotPWindow>0:
            #     fig.tight_layout()
            plotNo = min(plotPWindow, imgNo-ind)
            fig, axs = plt.subplots(plotNo , 3, figsize=(3*3, 1.8*plotNo), squeeze = False)
            plt.subplots_adjust(hspace=0.14, wspace=0.12)
            
        rotated_ = rotate(columnDensities[ind], angle_deg, reshape = False)[rowstart:rowend,columnstart:columnend]
        # rotated_=columnDensities[ind]
        if ind==0: #first time
            rotated_columnDensities =np.zeros((imgNo, *np.shape(rotated_)))
        rotated_columnDensities[ind] = rotated_
    
        #preview:
        dx=params.camera.pixelsize_meters/params.magnification
        
        popt0, popt1 = ImageAnalysisCode.fitgaussian2D(rotated_columnDensities[ind], dx=dx, 
                                                      do_plot = do_plot, ax=axs[plotInd], Ind=plotInd, imgNo=plotNo,
                                                      subtract_bg = subtract_bg, signal_feature = signal_feature, 
                                                      vmax = vmax, vmin = vmin,
                                                      title="1D density", title2D="column density",
                                                      xlabel1D="position ($\mu$m)", ylabel1D="1d density (atoms/$\mu$m)",                                                  
                                                      xscale_factor=1/units.um, yscale_factor=units.um)
        
        if variablesToDisplay is not None and variableLog is not None:
            variablesToDisplay = [ii.replace(' ','_') for ii in variablesToDisplay]
            axs[plotInd,0].text(0,1, 
                            variableLog.loc[logTime[ind]][variablesToDisplay].to_string(name=showTimestamp).replace('Name','Time'), 
                            fontsize=5, ha='left', va='top', transform=axs[plotInd,0].transAxes, 
                            bbox=dict(boxstyle="square", ec=(0,0,0), fc=(1,1,1), alpha=0.7))
            
            
        if popt0 is not None and popt1 is not None:
                
            amp_x, center_x, width_x, _ = popt0/units.um
            amp_y, center_y, width_y, _ = popt1/units.um
            
            # guess = [amp_g, center_g, w_g, offset_g]
            
            
            # wx = abs(popt0[2])
            AtomNumberX = amp_x * width_x * (2*np.pi)**0.5 * units.um * units.um
            
            # wy = abs(popt1[2])
            AtomNumberY = amp_y * width_y * (2*np.pi)**0.5 * units.um * units.um
            
            AtomNumbers.append(AtomNumberY)
            print("\n{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
            # width_x = popt0[2]/units.um
            
            print("RMS cloud size x: {:.2f} um".format(width_x))
            print("RMS cloud size y: {:.2f} um".format(width_y))
            print("x center: {:.2f} um".format(center_x))
            print("y center: {:.2f} um".format(center_y))
            centers_x.append(center_x)
            centers_y.append(center_y)
            widths_x.append(width_x)
            widths_y.append(width_y)
            currentDataDict['amplitude_xfit']=amp_x
            currentDataDict['amplitude_yfit']=amp_y
            currentDataDict['AtomNumber_yfit']=AtomNumberY
            currentDataDict['center_x']=center_x
            currentDataDict['center_y']=center_y
            currentDataDict['width_x']=width_x
            currentDataDict['width_y']=width_y
        else:
            log("fit failed for index {}".format(ind),2)
        dataDicts.append(currentDataDict)
    
    # runDataFrame = pd.DataFrame()
    
    # runDataFrame['AtomNumber_yfit']=AtomNumbers
    # runDataFrame['centers_x']=centers_x
    # runDataFrame['widths_x']=widths_x
    # runDataFrame['centers_y']=centers_y
    # runDataFrame['widths_y']=widths_y
    
    fig.tight_layout()
    
    print('\nThe average number of atoms:{:.2e}'.format(np.mean(AtomNumbers)))
        
    print("Mean RMS width x: {:.2f} +/- {:.2f} um".format(np.mean(widths_x), np.std(widths_x)))
    print("Mean RMS width y: {:.2f} +/- {:.2f} um".format(np.mean(widths_y), np.std(widths_y)))
    
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    
    if repetition > 1:
        widths_y = np.array(widths_y).reshape(-1, repetition)
        AtomNumbers = np.array(AtomNumbers).reshape(-1, repetition)
        
        widths_y_std = widths_y.std(axis=1)
        AtomNumbers_std = AtomNumbers.std(axis=1)
        
        widths_y = widths_y.mean(axis=1)
        AtomNumbers = AtomNumbers.mean(axis=1)
    else:
        widths_y_std = None
        AtomNumbers_std = None
    
    xx = np.arange(len(widths_y))
    
    ax1.errorbar(xx, widths_y, widths_y_std, capsize=8, color='tab:orange')
    ax1.plot(xx, widths_y, '.', color='tab:orange')
    ax1.set_ylabel('Y Widths (µm)', color='tab:orange')
    ax1.tick_params(axis="y", labelcolor='tab:orange')
    
    ax2.errorbar(xx, AtomNumbers, AtomNumbers_std, capsize=5, color='tab:green')
    ax2.plot(xx, AtomNumbers, '.-', color='tab:green')
    ax2.set_ylabel('Atom Number', color='tab:green')
    ax2.tick_params(axis="y", labelcolor='tab:green')
    
    return dataDicts



# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()

# ax1.plot(xx, widths_y, '.-', color='tab:orange')
# ax1.set_ylabel('Y Widths (µm)', color='tab:orange')
# ax1.tick_params(axis="y", labelcolor='tab:orange')

# ax2.plot(xx, AtomNumbers, '.-', color='tab:green')
# ax2.set_ylabel('Atom Number', color='tab:green')
# ax2.tick_params(axis="y", labelcolor='tab:green')

# fig.tight_layout()
# plt.rcParams.update({'font.size': 12})
# plt.figure(figsize=(5,4))
# plt.plot(centers_y, AtomNumbers,'o')
# plt.xlabel(r"Trap Position ($\mu$m)")
# plt.ylabel("Atom Number")
# plt.tight_layout()
# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
# #Save the numbers:
# np.savetxt(data_folder+"/AtomNumbers.txt",AtomNumbers)
# np.savetxt(data_folder+"/centers_y.txt", centers_y)
# plt.show()


#Temperature fit
# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="y", 
#                                              do_plot = True, save_folder = data_folder)

# popt, pcov = ImageAnalysisCode.thermometry1D(params, rotated_columnDensities, tof_array, thermometry_axis="x", 
#                                              do_plot = True, save_folder = data_folder)


if __name__ == "__main__":
    AnalyzeZyla(date, data_folder1,
                data_path)