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
date = '5/7/2024'
data_folder = [
    # r'Andor/Potential Modulation_5',
    r'Andor/Potential Modulation ODT 400 Low Frequency 500 ms',


    # r'Andor/ODT 2650 ZS bias 5_2',
    # r'Andor/ODT 2650 ZS bias 5_3',
    # r'Andor/ODT 2650 ZS bias 5_1',
    # r'Andor/ODT 2650 ZS bias 5_4',
    # r'Andor/ODT 2650 ZS bias Negative Vertical',
    # r'Andor/ODT 2650 ZS bias Negative Vertical_1',
    # r'Andor/ODT 3400 Bias Scan',
    # r'Andor/ODT 3400 Bias Scan_1'
    # r'Andor/ODT 3400 Bias Fine Scan',
    # r'Andor/ODT 3400 Bias Fine Scan_2',
    # r'Andor/ODT 400 ZS bias 5.1_1',
    # r'Andor/ODT 400 ZS bias 5.1',
    # r'Andor/ODT 400 ZS bias 4.8',
    # r'Andor/ODT 400 ZS bias 4.8_1',
    # r'Andor/ODT 1150 Bias Scan_1',
    # r'Andor/ODT 1150 Bias Scan_2',
    # r'Andor/ODT 1150 Bias Scan_3',
    # r'Andor/ODT 3400 Coarse Scan_4',
    # r'Andor/ODT 3400 Coarse Scan_5',
    # r'Andor/ODT 3400 Coarse Scan_6',
    # r'Andor/ODT 1150 Scan',
    # r'Andor/Vary EvapTime and Tau',
    # r'Andor/D1 bias scan Negative Polarity', 
    # r'Andor/D1 bias scan Positive Polarity'
    ]
####################################
#Parameter Setting
####################################
repetition = 1 #The number of identical runs to be averaged.
subtract_burntin = 0
examNum = 5 #The number of runs to exam.
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


# rowstart = 10
# rowend = -10
# columnstart = 10
# columnend = -10

columnstart=800
columnend=1200

# # ODT 400
rowstart = 1000
rowend = 1125
# columnstart=750
# columnend=1150

# # ODT 2560
# rowstart = 350
# rowend = 450
# columnstart=800
# columnend=1150

# # ODT 3400
# rowstart = 120
# rowend = 230
# columnstart=800
# columnend=1150

# rowstart = 660
# rowstart = 500
# rowend = -10
# columnstart = 600
# columnend = -200

# rowstart =750 #ODT 2675
# rowend = 830
# rowstart =800 #ODT1150
# rowend = 900
# # rowstart =616 #ODT1675
# # rowend = 651
# rowstart =570 #ODT1900
# rowend = 670
# # rowstart = 800 #ODT990
# # rowend = 835

# rowstart = 888 #ODT700
# rowend = 923
# rowstart = 1078 #ODT50
# rowend = 1113

# rowstart = 443 #ODT3800
# rowend = 478

rowstart -= 50
rowend += 50

####################################
####################################
dayfolder = ImageAnalysisCode.GetDataLocation(date, DataPath=dataRootFolder)
dataPath = [ os.path.join(dayfolder, f) for f in data_folder]
# variableLog_folder = dayFolder + r'/Variable Logs'
examFrom, examUntil = ImageAnalysisCode.GetExamRange(examNum, examFrom, repetition)

params = ImageAnalysisCode.ExperimentParams(date, t_exp = 10e-6, picturesPerIteration=None, cam_type = "zyla")
dxMicron = params.camera.pixelsize_microns/params.magnification    #The length in micron that 1 pixel correspond to. 
dxMeter = params.camera.pixelsize_meters/params.magnification    #The length in meter that 1 pixel correspond to. 

columnDensities, variableLog = ImageAnalysisCode.PreprocessZylaImg(*dataPath, examFrom=examFrom, examUntil=examUntil, 
                                                                   rotateAngle=rotateAngle, 
                                                                   rowstart=rowstart, rowend=rowend, 
                                                                   columnstart=columnstart, columnend=columnend,
                                                                   subtract_burntin=subtract_burntin, 
                                                                   showRawImgs=1)
#%%
        
popts, bgs = ImageAnalysisCode.FitColumnDensity(columnDensities, dx = dxMicron, mode='both', yFitMode='single',
                                                subtract_bg=0, Xsignal_feature='wide', Ysignal_feature='narrow',
                                                rowstart=rowstart, rowend=rowend, columnstart=columnstart, columnend=columnend)

results = ImageAnalysisCode.AnalyseFittingResults(popts, logTime=variableLog.index) 

if variableLog is not None:
    results = results.join(variableLog)
# results.to_csv('Test.csv')
#%%
# results = results[ results.YatomNumber < 1e7 ]

# %%
ImageAnalysisCode.PlotFromDataCSV(results, 'fmod_kHz', 'YatomNumber', 
                                  # iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                  groupbyX=1, threeD=0,
                                  figSize = 0.5
                                  )

ImageAnalysisCode.PlotFromDataCSV(results, 'fmod_kHz', 'Ywidth', 
                                  # iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                  groupbyX=1, threeD=0,
                                  figSize = 0.5
                                  )

# ImageAnalysisCode.PlotFromDataCSV(results, 
#                                   'ZSBiasCurrent',
#                                   'YatomNumber', 
#                                   ['VerticalBiasCurrent>2.4', 'VerticalBiasCurrent<2.8'],
#                                   iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   groupbyX=1, threeD=0,
#                                   figSize = 0.5
#                                   )

# %%

intermediatePlot = 1
plotPWindow = 6
plotRate = 1
uniformscale = 0
rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9}

variablesToDisplay = [
                    # # 'Coil_medB', 
                        'TOF',
                        'fmod_kHz',
                        'tmod_ms',
                        # 'Evap_Tau',
                       # 'ZSBiasCurrent',
                       # 'VerticalBiasCurrent',
                        # 'CamBiasCurrent'
                      ]
showTimestamp = False
# variablesToDisplay=None
textY = 1
textVA = 'bottom'

if intermediatePlot:
    # ImageAnalysisCode.ShowImagesTranspose(images_array, uniformscale=False)
    ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, 
                                          dx=dxMicron, 
                                            # filterLists=[['fmod_kHz>32']],
                                           plotRate=plotRate, plotPWindow=plotPWindow,
                                          variablesToDisplay = variablesToDisplay,
                                          variableLog=variableLog, 
                                          logTime=variableLog.index,
                                          textLocationY=0.8, rcParams=rcParams)
    
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
    
# %%
# import matplotlib.dates as mdates

# fig, ax = plt.subplots(1,1, figsize=(8,6), layout='constrained')
# # plt.plot(results.wait, results.YatomNumber)


# ax.plot(results.index, results.Ycenter.values, '.')
# ax.set(xlabel='time (Day HH:MM)', ylabel='y center (µm)')
# ax.set_xticks(ax.get_xticks()[::2])
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))

# %%

# ImageAnalysisCode.PlotFromDataCSV(results, 
#                                   xVariable='ZSBiasCurrent',
#                                   yVariable='YatomNumber', 
#                                   iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   groupbyX=1, threeD=0,
#                                   figSize = 0.5
#                                   )

# # %%

# atomNumber = results[['YatomNumber', 'ZSBiasCurrent', 'VerticalBiasCurrent']].groupby(['ZSBiasCurrent', 'VerticalBiasCurrent']).mean().YatomNumber.unstack()
# fig, ax = plt.subplots(1,1, figsize=(8,6), layout='constrained')

# plt.pcolormesh(atomNumber.columns.values, atomNumber.index.values, atomNumber.values)