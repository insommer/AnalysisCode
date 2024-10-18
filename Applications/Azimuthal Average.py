# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 13:41:52 2024

@author: insommer
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
dataRootFolder = r"C:\Lehigh University Dropbox\Ariel Sommer\Sommer Lab Shared\Data"
date = '8/30/2024'
# date = '11/07/2023'

data_folder = [
    # r'Andor/Xcenter vs wait'
    r'Andor/Probe evap',
    # r'Andor/Probe evap',
    # r'Andor/ODT -1 Align',
    ]
####################################
#Parameter Setting
####################################
reanalyze = 1
saveresults = 0
overwriteOldResults = 0

repetition = 1 #The number of identical runs to be averaged.
subtract_burntin = 0
skipFirstImg = 'auto'
# skipFirstImg = 0

rotateAngle = 1.5 #rotates ccw
# rotateAngle = 0 #rotates ccw

examNum = 6 #The number of runs to exam.
examFrom = 0 #Set to None if you want to check the last several runs. 
showRawImgs = 0

# Set filters for the data, NO SPACE around the operator.
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
    ] 

pictureToHide = None
# pictureToHide = [0,1,2,3] # list(range(0,10,2))

subtract_bg = 0
signal_feature = 'wide' 
signal_width = 10 #The narrower the signal, the bigger the number.
fitbgDeg = 5

rowstart = 10
rowend = -10
columnstart = 10
columnend = -10


# columnstart=600
# columnend= 1100

# rowstart = 800
# rowend = -900

# rowstart = 730
# rowend = 1070

# # second pass
# rowstart = 1330
# rowend = 1375


####################################
####################################
dayfolder = ImageAnalysisCode.GetDataLocation(date, DataPath=dataRootFolder)
dataPath = [ os.path.join(dayfolder, f) for f in data_folder]

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
                                                                   showRawImgs=showRawImgs, rebuildCatalogue=0,
                                                                   # filterLists=[['Evap_timestep>0.2']]
                                                                   )
#%%
autoCrop = 1
if autoCrop:
    columnDensities_croped = ImageAnalysisCode.AutoCrop(columnDensities, sizes=[200, 100])
    print('ColumnDensities auto cropped.')
# %% Plot the first 6 images.
fig, axes = plt.subplots(1, 6, figsize=(10, 2), sharex=True, sharey=True, layout='constrained')

for i in range(6):
    axes[i].imshow(columnDensities_croped[i], cmap='gray')
    
#%% Show the pictures
arange, sizes = ImageAnalysisCode.PlotArangeAndSize(len(columnDensities_croped), 
                                                    sizes_ratio=(3, 2))
fig, axes = plt.subplots(*arange, figsize=sizes, 
                         squeeze = False, sharex=1, sharey=1, layout='constrained')
for ii, ax in enumerate(axes.flatten()):
    ax.imshow(columnDensities_croped[ii], cmap='gray')
    ax.text(0.05, 0.05, '{}'.format(ii), ha='left', va='bottom', transform=ax.transAxes,
            bbox=dict(boxstyle="square", ec=(0,0,0), fc=(1,1,1), alpha=0.7))
#%% Calculate azimuthal average
result = []
for img in columnDensities_croped:
    result.append(ImageAnalysisCode.AzimuthalAverage(img, sigma=1, 
                                                     do_plot=0))

#%% Plot the results
arange, sizes = ImageAnalysisCode.PlotArangeAndSize(len(result))
fig, axes = plt.subplots(*arange, figsize=sizes,
                          layout='constrained')

for ii, ax in enumerate(axes.flatten()):
    ImageAnalysisCode.plotRadialAtomDensity(*result[ii], dx=dxMicron, ax=ax, 
                                            linestyle='-', addAxieLabel=1)
    ax.text(0.05, 0.05, '{}'.format(ii), ha='left', va='bottom', transform=ax.transAxes,
            bbox=dict(boxstyle="square", ec=(0,0,0), fc=(1,1,1), alpha=0.7))
    
