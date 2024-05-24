# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:34:22 2023

@author: Sommer Lab
"""
from ImageAnalysis import ImageAnalysisCode
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import pandas as pd
import os
from scipy import constants

####################################
#Set the date and the folder name
####################################
dataRootFolder = r"D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data"
date = '5/24/2024'
data_folder = [
    # r'Andor/ODT 400 Modulation 0.1 V 10-50 kHz Variable tmod_1',
    # r'Andor/D1 bias scan Negative Polarity', 
    # r'Andor/D1 bias scan Positive Polarity'
    # 'Andor\Misaligned ODT vs wait',
    # 'Andor\Med Field Wait_2',
    # 'Andor\Misaligned ODT vs wait',
    # 'Andor\ODT 400 Misalign',
    'Andor\Test',
    'Andor\Test ODT 400',
    # 'Andor\Test ODT 400',
    # 'Andor\Test 3 PPI_1'
    ]
####################################
#Parameter Setting
####################################
reanalyze = 1
saveresults = 1
overwriteOldResults = 0

repetition = 1 #The number of identical runs to be averaged.
subtract_burntin = 0
skipFirstImg = 'auto'
# skipFirstImg = 0
examNum = None #The number of runs to exam.
examFrom = None #Set to None if you want to check the last several runs. 




variableFilterList = [
    # [# 'wait==50', 
    # # 'VerticalBiasCurrent==0'
    # 'fmod_kHz==0',
    # # 'Evap_Tau==0.1',
    # # 'Evap_Time_1==2'
    # ], 
    # [
    # 'TOF==0',
    # 'Evap_Tau==0.1',
    # 'Evap_Time_1==2']
    ] # NO SPACE around the operator!

pictureToHide = None
# pictureToHide = [0,1,2,3] # list(range(0,10,2))

subtract_bg = 1
signal_feature = 'wide' 
signal_width = 10 #The narrower the signal, the bigger the number.
fitbgDeg = 5
rotateAngle = 0.5 #rotates ccw


rowstart = 10
rowend = -10
columnstart = 10
columnend = -10

columnstart=400
columnend=1150 

# columnstart=500
# columnend=800

# columnstart=300
# columnend=1000 

# rowstart = 1010
# rowend = 1110

# rowstart = 250
# rowend = 500

# rowstart = 1000	#ODT 400
# rowend = 1125	

# rowstart -= 200
# rowend += 200

# rowstart += 20
# rowend -= 20

# columnstart -= 300
# columnend += 300

####################################
####################################
dayfolder = ImageAnalysisCode.GetDataLocation(date, DataPath=dataRootFolder)
dataPath = [ os.path.join(dayfolder, f) for f in data_folder]

# dataPath = [
#     # r'D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data\2024\05-2024\07 May 2024\Andor\ODT 1900 Misalign',
#     # r'D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data\2024\05-2024\23 May 2024\Andor\Test',
# #     # r'D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data\2024\05-2024\15 May 2024\Andor\ODT 1900 Modulation 0.1 V 10-50 kHz Variable tmod',
# #     # r'D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data\2024\05-2024\15 May 2024\Andor\ODT 2650 Modulation 0.1 V 10-50 kHz Variable tmod',
# #     r'D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data\2024\05-2024\15 May 2024\Andor\ODT 3400 Modulation 0.1 V 10-50 kHz Variable tmod_1',
# #     r'D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data\2024\05-2024\16 May 2024\Andor\ODT 4150 Modulation 0.1 V 10-50 kHz Variable tmod'
# ]
    



# variableLog_folder = dayFolder + r'/Variable Logs'
examFrom, examUntil = ImageAnalysisCode.GetExamRange(examNum, examFrom, repetition)

params = ImageAnalysisCode.ExperimentParams(date, t_exp = 10e-6, picturesPerIteration=None, cam_type = "zyla")
dxMicron = params.camera.pixelsize_microns/params.magnification   #The length in micron that 1 pixel correspond to. 
dxMeter = params.camera.pixelsize_meters/params.magnification    #The length in meter that 1 pixel correspond to. 



#%%

columnDensities, variableLog = ImageAnalysisCode.PreprocessZylaImg(*dataPath, examRange=[examFrom, examUntil], 
                                                                   rotateAngle=rotateAngle, 
                                                                   rowstart=rowstart, rowend=rowend, 
                                                                   columnstart=columnstart, columnend=columnend,
                                                                   subtract_burntin=subtract_burntin, 
                                                                   skipFirstImg=skipFirstImg, 
                                                                   showRawImgs=0, rebuildCatalogue=0)
#%%
        
popts, bgs = ImageAnalysisCode.FitColumnDensity(columnDensities, dx = dxMicron, mode='both', yFitMode='single',
                                                subtract_bg=0, Xsignal_feature='wide', Ysignal_feature='narrow')

results = ImageAnalysisCode.AnalyseFittingResults(popts, logTime=variableLog.index) 

if variableLog is not None:
    results = results.join(variableLog)

if saveresults:
    ImageAnalysisCode.SaveResultsDftoEachFolder(results, overwrite=overwriteOldResults)    
#%%
# results = results[ results.fmod_kHz <12 ]

# %%
# ImageAnalysisCode.PlotFromDataCSV(results, 'fmod_kHz', 'YatomNumber', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   groupbyX=1, threeD=0,
#                                   figSize = 0.5
#                                   )

# ImageAnalysisCode.PlotFromDataCSV(results, 'fmod_kHz', 'Ywidth', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   # groupby='ODT_Position', 
#                                    groupbyX=1, 
#                                   threeD=0,
#                                   figSize = 0.5
#                                   )


# ImageAnalysisCode.PlotFromDataCSV(results, 'IterationNum', 'YatomNumber', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   # groupby='ODT_Position', 
#                                    groupbyX=1, 
#                                   threeD=0,
#                                   figSize = 0.5
#                                   )
fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
results.YatomNumber.plot(title='Atom Number', linestyle='', marker='.')

fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
results.Ycenter.plot(title='y Position', linestyle='', marker='.')

fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
results.Xcenter.plot(title='x Position', linestyle='', marker='.')

# %%

intermediatePlot = 1
plotPWindow = 3
plotRate = 1
uniformscale = 0
rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9}

variablesToDisplay = [
                    # # 'Coil_medB', 
                        'wait',
                        # 'ODT_Position',
                        # 'fmod_kHz',
                        # 'tmod_ms',
                        # 'Evap_Tau',
                        'VerticalBiasCurrent',
                        # 'YatomNumber',
                        # 'CamBiasCurrent'
                        'Folder',
                        'ODT Misalign'
                        
                      ]
showTimestamp = False
# variablesToDisplay=None
textY = 1
textVA = 'bottom'

if intermediatePlot:
    # ImageAnalysisCode.ShowImagesTranspose(images_array, uniformscale=False)
    ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, 
                                          dx=dxMicron, 
                                            # filterLists=[['fmod_kHz<12']],
                                           plotRate=plotRate, plotPWindow=plotPWindow,
                                            variablesToDisplay = variablesToDisplay,
                                           showTimestamp=showTimestamp,
                                          variableLog=results, 
                                          logTime=variableLog.index,
                                          textLocationY=1, rcParams=rcParams)
    
    # ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, dx=dx, 
    #                                       plotRate=1, plotPWindow=plotPWindow,
    #                                       variablesToDisplay = variablesToDisplay,
    #                                       variableLog=variableLog, logTime=variableLog.index,
    #                                       textLocationY=0.8, rcParams=rcParams)

    # xx = np.arange(len(imgs_oneD[0]))
    # fig, axes = plt.subplots(fileNo, 1, sharex=True, layout='constrained')
    # for ii in range(fileNo):        
    #     axes[ii].plot(imgs_oneD[ii], '.')
    #     axes[ii].plot(xx, ImageAnalysisCode.Gaussian(xx, *popt_Basler[ii]))
    #     axes[ii].text(0.9,0.8, files[ii], transform=axes[ii].transAxes)

    # c, w = np.array(popt_Basler).mean(axis=0)[1:-1]
    # axes[-1].set(xlim=[c-15*w, c+15*w])
#%%
