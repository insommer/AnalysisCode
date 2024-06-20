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
date = '6/20/2024'
data_folder = [
    # r'Andor/High Field Round-Trip Survival'
    # r'Andor/B Field Scan_2'
    # r'Andor/ODT 400 Modulation 0.1 V 10-50 kHz Variable tmod_1',
    # r'Andor/D1 bias scan Negative Polarity', 
    # r'Andor/D1 bias scan Positive Polarity'
    # 'Andor\Misaligned ODT vs wait',
    # 'Andor\Med Field Wait_2',
    # 'Andor\Misaligned ODT vs wait',
    'Andor\Imaging Frequency Scan Lens at 195',
    'Andor\Imaging Frequency Scan Lens at 200',
    'Andor\Imaging Frequency Scan Lens at 198',
    'Andor\Imaging Frequency Scan Lens at 197',
    'Andor\Imaging Frequency Scan Lens at 197.5',


    # 'Andor\GM Temperature',
    # 'Andor\With Low Servo Late ODT LowServo1 0.3',
    # 'Andor\Test'
    ]
####################################
#Parameter Setting
####################################
reanalyze = 1
saveresults = 1
overwriteOldResults = 0

repetition = 1 #The number of identical runs to be averaged.
subtract_burntin = 0
rotateAngle = 1 #rotates ccw
skipFirstImg = 'auto'
# skipFirstImg = 0
examNum = None #The number of runs to exam.
examFrom = None #Set to None if you want to check the last several runs. 
showRawImgs = 0

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

rowstart = 10
rowend = -10
columnstart = 10
columnend = -10


columnstart=1050
columnend=1350 

rowstart = 970
rowend = 1050

rowstart -= 0
rowend += 0

# # rowstart += 100
# # rowend -= 100

# columnstart -= 200
# columnend += 200

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
# if not reanalyze:
#     resultsList = []
#     for pp in dataPath:
#         resutlsPath = os.path.join(pp, 'results.pkl')        
#         if os.path.exists(resutlsPath):
#             with open(resutlsPath, 'rb') as f:
#                 resultsList.append( pickle.load(f) )
            

#%%
columnDensities, variableLog = ImageAnalysisCode.PreprocessZylaImg(*dataPath, examRange=[examFrom, examUntil], 
                                                                   rotateAngle=rotateAngle, 
                                                                   rowstart=rowstart, rowend=rowend, 
                                                                   columnstart=columnstart, columnend=columnend,
                                                                   subtract_burntin=subtract_burntin, 
                                                                   skipFirstImg=skipFirstImg, 
                                                                   showRawImgs=showRawImgs, rebuildCatalogue=0)

# variableLog.Lens_Position = 1.85

#%%
        
popts, bgs = ImageAnalysisCode.FitColumnDensity(columnDensities, dx = dxMicron, mode='both', yFitMode='single',
                                                subtract_bg=0, Xsignal_feature='wide', Ysignal_feature='narrow')

results = ImageAnalysisCode.AnalyseFittingResults(popts, logTime=variableLog.index) 

if variableLog is not None:
    results = results.join(variableLog)

if saveresults:
    ImageAnalysisCode.SaveResultsDftoEachFolder(results, overwrite=overwriteOldResults)    
#%%
results = results[ results.Ywidth < 10 ]

# %%
ImageAnalysisCode.PlotFromDataCSV(results, 'HF_AOM_Freq', 'Ywidth', 
                                   iterateVariable='Lens_Position', 
                                   # filterLists=[['Ywidth<10']],
                                  groupbyX=1, threeD=0,
                                  figSize = 0.5
                                  )

ImageAnalysisCode.PlotFromDataCSV(results, 'HF_AOM_Freq', 'YatomNumber', 
                                  # iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                  groupbyX=1, threeD=0,
                                  figSize = 0.5
                                  )

# ImageAnalysisCode.PlotFromDataCSV(results, 'fmod_kHz', 'Ywidth', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   # groupby='ODT_Position', 
#                                    groupbyX=1, 
#                                   threeD=0,
#                                   figSize = 0.5
#                                   )


# ImageAnalysisCode.PlotFromDataCSV(results, 'B_Field', 'YatomNumber', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   # groupby='ODT_Position', 
#                                     groupbyX=1, 
#                                   threeD=0,
#                                   figSize = 0.5
#                                   )
fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
results.YatomNumber.plot(title='Atom Number', linestyle='', marker='.')

# fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
# results.Ycenter.plot(title='y Position', linestyle='', marker='.')

# fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
# results.Xcenter.plot(title='x Position', linestyle='', marker='.')

# %%

intermediatePlot = 1
plotPWindow = 3
plotRate = 1
uniformscale = 0
rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9}

variablesToDisplay = [
                    # # 'Coil_medB', 
                        'wait',
                        'TOF',
                        'HF_AOM_Freq',
                        # 'FB Voltage',
                        # 'B_Field',
                        
                        # 'ODT_Position',
                        # 'fmod_kHz',
                        # 'tmod_ms',
                        # 'Evap_Tau',
                        # 'VerticalBiasCurrent',
                        # 'YatomNumber',
                        # 'HF_AOM_Freq',
                        # 'CamBiasCurrent',
                        'Lens Position'
                        
                      ]
showTimestamp = False
# variablesToDisplay=None
textY = 1
textVA = 'bottom'

if intermediatePlot:
    # ImageAnalysisCode.ShowImagesTranspose(images_array, uniformscale=False)
    ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, 
                                          dx=dxMicron, 
                                            # filterLists=[['YatomNumber>12']],
                                           plotRate=plotRate, plotPWindow=plotPWindow,
                                            variablesToDisplay = variablesToDisplay,
                                           showTimestamp=showTimestamp,
                                          variableLog=results, 
                                          logTime=variableLog.index,
                                          uniformscale=uniformscale,
                                          textLocationY=0.9, rcParams=rcParams,
                                          figSizeRate=1)
    
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


# dataFrame = dataFrame[ dataFrame.width_y < 600 ]
# # dataFrame = dataFrame[ dataFrame.Evap_Tau > 0.05 ]


# %%
# var2 = 'Evap_Time_1'
# var1 = 'D1_Re_Attn'
# df = ImageAnalysisCode.multiVariableThermometry(
#                                                 # results[results.Ywidth<22], 
#                                                 results, 
                                                
#                                                 var1, var2, 
#                                                 fitXVar='TOF',  fitYVar='Ywidth',do_plot=1, add_Text=0)

# %%
# cols = ['PSD', 'T (K)', 'AtomNum']

# for r in cols:
#     df1 = df[r].unstack()
#     fig, ax = plt.subplots(1,1, figsize=[4,3], layout='constrained')
#     cax = ax.pcolormesh(df1.columns, df1.index, df1.values, cmap='viridis')
#     ax.set(xlabel=var2, ylabel=var1, title=r)
#     fig.colorbar(cax, ax=ax)

# fig, ax = plt.subplots(1,1, layout='constrained')
# ax.plot(df['T (K)'], df.PSD, '.')
# ax.set(xlabel='T (K)', ylabel='PSD')

# fig, ax = plt.subplots(1,1, layout='constrained')
# ax.plot(df['AtomNum'], df.PSD, '.')
# ax.set(xlabel='AtomNum', ylabel='PSD')

# %%
# print('======================')
# print('The phase space density is:\n{:.5e}'.format(df.PSD.iloc[0]))

