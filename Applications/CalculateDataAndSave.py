# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 17:11:59 2023

@author: Sommer Lab
"""

from ImageAnalysis import ImageAnalysisCode
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import pandas as pd
import os





DEBUG=True
def log(msg,indent=0):
    if DEBUG:
        print("#########  "+'\t'*indent+ msg)
        





def CalculateFromZyla(dayFolderPath, 
                dataFolders, 
                variableLog=None,                
                repetition=1, 
                examNum=None, 
                examFrom=None, 
                plotPWindow=3, do_plot=False, uniformscale=0, 
                variablesToDisplay= None, variableFilterList=None, showTimestamp=False, pictureToHide=None,
                subtract_bg=True, signal_feature='narrow', 
                rowstart=10, rowend=-10, 
                columnstart=10, columnend=-10,
                angle_deg= 2, #rotates ccw
                picturesPerIteration = 3, 
                lengthFactor=1e-6):
    
    
    dataFolderPaths = [ os.join.path(dayFolderPath, f) for f in dataFolders ]
    examFrom, examUntil = ImageAnalysisCode.GetExamRange(examNum, examFrom, repetition)
    
    params = ImageAnalysisCode.ExperimentParams(t_exp = 10e-6, picturesPerIteration= picturesPerIteration, cam_type = "zyla")
    images_array = None
    NoOfRuns = []
    
    log("Loading spooling series...")
    for ff in dataFolderPaths:
        log("loading "+os.path.split(ff)[-1],1)
        if images_array is None:
            images_array, fileTime = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder = ff, 
                                                                       return_fileTime=1)
            NoOfRuns.append(len(fileTime))
        else:
            _images_array, _fileTime = ImageAnalysisCode.LoadSpooledSeries(params = params, data_folder = ff, 
                                                                           return_fileTime=1)
            images_array = np.concatenate([images_array, _images_array], axis=0)
            fileTime = fileTime + _fileTime
            NoOfRuns.append(len(fileTime))
    
    images_array = images_array[examFrom: examUntil]
    fileTime = fileTime[examFrom: examUntil]
    
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
    
    
    imgNo = len(columnDensities)
    
    
    AtomNumberX=[]
    AtomNumberY=[]
    widths_x = []
    widths_y = []
    centers_x = []
    centers_y = []
    
    if uniformscale:
        vmax = columnDensities.max()
        vmin = columnDensities.min()
    else:
        vmax = None
        vmin = 0
        
    log("Running gaussian fits...")
    for ind in range(imgNo):        
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
                                                      do_plot = do_plot, ax=None, Ind=None, imgNo=None,
                                                      subtract_bg = subtract_bg, signal_feature = signal_feature, 
                                                      vmax = vmax, vmin = vmin,
                                                      title="1D density", title2D="column density",
                                                      xlabel1D="position ($\mu$m)", ylabel1D="1d density (atoms/$\mu$m)",                                                  
                                                      xscale_factor=1/lengthFactor, yscale_factor=lengthFactor)
        
        if variablesToDisplay is not None and variableLog is not None:
            variablesToDisplay = [ii.replace(' ','_') for ii in variablesToDisplay]
            axs[plotInd,0].text(0,1, 
                            variableLog.loc[logTime[ind]][variablesToDisplay].to_string(name=showTimestamp).replace('Name','Time'), 
                            fontsize=5, ha='left', va='top', transform=axs[plotInd,0].transAxes, 
                            bbox=dict(boxstyle="square", ec=(0,0,0), fc=(1,1,1), alpha=0.7))
        
        
        if popt0 is None:
            center_x, width_x, atomNumberX = [np.nan] * 3
        else:            
            amp_x, center_x, width_x, _ = popt0
            atomNumberX = amp_x * width_x * (2*np.pi)**0.5
            
        if popt1 is None:
            center_y, width_y, atomNumberY = [np.nan] * 3
        else:                    
            amp_y, center_y, width_y, _ = popt1
            atomNumberY = amp_y * width_y * (2*np.pi)**0.5 
            
            # print("\n{}. Atom Number from gauss fit = {:.2e}".format(ind, AtomNumberY))
            
            # print("RMS cloud size x: {:.2f} um".format(width_x))
            # print("RMS cloud size y: {:.2f} um".format(width_y))
            # print("x center: {:.2f} um".format(center_x))
            # print("y center: {:.2f} um".format(center_y))
            centers_x.append(center_x)
            centers_y.append(center_y)
            widths_x.append(width_x)
            widths_y.append(width_y)
            AtomNumberX.append(atomNumberX)
            AtomNumberY.append(atomNumberY)
            
    return pd.DataFrame( np.array( [centers_y, widths_y, AtomNumberY, centers_x, widths_x, AtomNumberX] ).T,
                        index=logTime, columns=['Ycenter', 'Ywidth', 'AtomNumber', 'Xcenter', 'Xwidth', 'AtomNumberX'])
        

    
    # print('\nThe average number of atoms:{:.2e}'.format(np.mean(AtomNumbers)))
        
    # print("Mean RMS width x: {:.2f} +/- {:.2f} um".format(np.mean(widths_x), np.std(widths_x)))
    # print("Mean RMS width y: {:.2f} +/- {:.2f} um".format(np.mean(widths_y), np.std(widths_y)))
    



























totalDataPath =r"C:\Users\Sommer Lab\Documents\Data"
date = '9/26/2023'

dataCSV_filename = 'data_10_3.csv'
dayFolderPath = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)
dataFolders = [r'/Andor/Test_1']

dataCSV_filePath = os.path.join(dayFolderPath, dataCSV_filename)


if os.path.exists(dataCSV_filePath):
    dataCSV = pd.read_csv(dataCSV_filePath)
    dataCSV.time = pd.to_datetime(dataCSV.time)
    dataCSV.set_index('time', inplace=True)

else:
    variableLog_folder = os.path.join(dayFolderPath, 'Variable Logs')
    dataCSV = ImageAnalysisCode.LoadVariableLog(variableLog_folder)
    
    
    
results = AnalyzeZyla(dayFolderPath, 
                dataFolders, 
                variableLog = dataCSV,                
                repetition=1)

combined = results.join( dataCSV, how='outer', lsuffix='_l')

combined.to_csv( os.path.join(dayFolderPath, dataCSV_filename) )

    
