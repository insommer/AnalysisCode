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

####################################
#Set the date and the folder name
####################################
dataRootFolder = r"D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data"
date = '4/12/2024'
data_folder = [
    r'Andor/ODT 1150 Bias Scan',
    r'Andor/ODT 1150 Bias Scan_1',

    # r'Andor/ODT 1900',
    # r'Andor/Test'
    ]
####################################
#Parameter Setting
####################################
repetition = 1 #The number of identical runs to be averaged.
subtract_burntin = 0
skipFirstImg = 1
examNum = None #The number of runs to exam.
examFrom = None #Set to None if you want to check the last several runs. 
plotPWindow = 6
intermediatePlot = 1
uniformscale = 0
rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9}

variablesToDisplay = [
                    # # 'Coil_medB', 
                        'wait',
                        # 'ODT Misalign',
                        'ODT Position',
                      'ZSBiasCurrent',
                      'VerticalBiasCurrent',
                        # 'CamBiasCurrent'
                      ]
showTimestamp = False
# variablesToDisplay=None
textY = 1
textVA = 'bottom'

variableFilterList = None
variableFilterList = [
    # 'wait==50', 
    ] # NO SPACE around the operator!

pictureToHide = None
# pictureToHide = [0,1,2,3] # list(range(0,10,2))

subtract_bg = 0
signal_feature = 'narrow' 
signal_width = 10 #The narrower the signal, the bigger the number.
fitbgDeg = 5
rotateAngle = 0.5 #rotates ccw

rowstart = 10
rowend = -10
columnstart = 10
columnend = -10

# rowstart = 660
# rowstart = 500
# rowend = -10
# columnstart = 600
# columnend = -200

columnstart = 770
columnend = 1100


# rowstart =250
# rowend = -500
# columnstart = 600
# columnend = -550

# rowstart =200
# # rowend = -150
# columnstart = 500
# columnend = -300

# rowstart =750 #ODT 2675
# rowend = 830
# # rowstart =616 #ODT1675
# # rowend = 651
# rowstart =570 #ODT1900
# rowend = 670
rowstart =790 #ODT1150
rowend = 890
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

pPI = 4 if subtract_burntin else 3
params = ImageAnalysisCode.ExperimentParams(date, t_exp = 10e-6, picturesPerIteration=pPI, cam_type = "zyla")

columnDensities, variableLog = ImageAnalysisCode.PreprocessZylaImg(*dataPath, examFrom=examFrom, examUntil=examUntil, 
                                                                   rotateAngle=rotateAngle, subtract_burntin=subtract_burntin, 
                                                                   skipFirstImg=skipFirstImg)
columnDensities = columnDensities[:, rowstart:rowend, columnstart:columnend]

dx = params.camera.pixelsize_microns/params.magnification #The length in micron that 1 pixel correspond to. 
YcolumnDensities = columnDensities.sum(axis=2) * dx / 1e6**2

popts = []
bgs = []
for ydata in YcolumnDensities:
    popt, bg = ImageAnalysisCode.fitMultiGaussian(ydata, dx=dx, 
                                                  subtract_bg=subtract_bg, signal_feature=signal_feature, 
                                                  fitbgDeg=3, amp=1, width=3, denoise=0)
    popts.append(popt)
    bgs.append(bg)
    
# XcolumnDensities = columnDensities.sum(axis=1) * dx / 1e6**2
# poptsX = []
# for xdata in XcolumnDensities:
#     poptx,_ = ImageAnalysisCode.fitSingleGaussian(xdata, dx=dx,
#                                                   subtract_bg=1, signal_feature='wide')
#     poptsX.append(poptx)

results = ImageAnalysisCode.AnalyseFittingResults(popts, logTime=variableLog.index)

if variableLog is not None:
    results = results.join(variableLog)
# results.to_csv('0305.csv')

# %%
ImageAnalysisCode.PlotFromDataCSV(results, 'ZSBiasCurrent', 'YatomNumber', 
                                  iterateVariable='VerticalBiasCurrent', 
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
if intermediatePlot:
    # ImageAnalysisCode.ShowImagesTranspose(images_array, uniformscale=False)
    ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, dx=dx, 
                                          plotPWindow=plotPWindow,
                                          variablesToDisplay = variablesToDisplay,
                                          variableLog=variableLog, logTime=variableLog.index,
                                          textLocationY=0.8, rcParams=rcParams)

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