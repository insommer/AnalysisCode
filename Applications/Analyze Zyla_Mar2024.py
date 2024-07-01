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
date = '7/1/2024'
# date = '11/07/2023'

data_folder = [
    # r'Andor/High Field Round-Trip Survival'
    # r'Andor/B Field Scan_2'
    # r'Andor/ODT 400 Modulation 0.1 V 10-50 kHz Variable tmod_1',
    # r'Andor/D1 bias scan Negative Polarity', 
    # r'Andor/D1 bias scan Positive Polarity'
    # r'Andor\Misaligned ODT vs wait',
    # r'Andor\Med Field Wait_2',
    # r'Andor\Misaligned ODT vs wait',
    # r'Andor\Vary Evap Time 2_1'
    # r'Andor\gray molasses more D1 power_1'
    # r'Andor\With Low Servo Late ODT LowServo1 0.29 Changed Lens_1',
    # r'Andor\before evap thermometry low field Digital Mag Off',
    r'Andor\GM temp longer tof'
    # r'Andor\PSD Before Evap',
    ]
####################################
#Parameter Setting
####################################
reanalyze = 1
saveresults = 1
overwriteOldResults = 1

repetition = 1 #The number of identical runs to be averaged.
subtract_burntin = 1
rotateAngle = 0 #rotates ccw
skipFirstImg = 'auto'
# skipFirstImg = 0
examNum = None #The number of runs to exam.
examFrom = None #Set to None if you want to check the last several runs. 
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


# # columnstart=1000
# columnend=1100

# rowstart = 990
# rowend = 1000

# # # rowstart -= -5
# # # rowend += 5

# rowstart -= 0
# rowend += 80

# rowstart -= 300
# rowend += 300



# columnstart -= 700
# # columnend += 200

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

# if autoCrop:
#     columnDensities = ImageAnalysisCode.AutoCrop(columnDensities, xsize=200, ysize=50)
#%%
        
popts, bgs = ImageAnalysisCode.FitColumnDensity(columnDensities, dx = dxMicron, mode='both', yFitMode='single',
                                                subtract_bg=subtract_bg, Xsignal_feature='wide', Ysignal_feature='wide')

results = ImageAnalysisCode.AnalyseFittingResults(popts, logTime=variableLog.index) 

if variableLog is not None:
    results = results.join(variableLog)

if saveresults:
    ImageAnalysisCode.SaveResultsDftoEachFolder(results, overwrite=overwriteOldResults)    
#%%
# mask = (results.Ywidth < 5.5)  & (results.Ywidth > 3.9) & (results.HF_AOM_Freq<=317)
# results = results[ mask ]
# columnDensities = columnDensities[mask]
# popts[0] = np.array(popts[0])[mask]
# popts[1] = np.array(popts[1])[mask]
# bgs[0] = np.array(bgs[0])[mask]
# bgs[1] = np.array(bgs[1])[mask]


# %%
# ImageAnalysisCode.PlotFromDataCSV(results, 'HF_AOM_Freq', 'Ywidth', 
#                                    iterateVariable='Lens_Position', 
#                                    # filterLists=[['Ywidth<10']],
#                                   groupbyX=1, threeD=0,
#                                   figSize = 0.5
#                                   )

# ImageAnalysisCode.PlotFromDataCSV(results, 'HF_AOM_Freq', 'YatomNumber', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   groupbyX=1, threeD=0,
#                                   figSize = 0.5
#                                   )
ImageAnalysisCode.PlotFromDataCSV(results, 'TOF', 'YatomNumber', 
                                  # iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                   groupbyX=0, threeD=0,
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
                        'Evap_Time_2',
                        'TOF'
                        #'HF_AOM_Freq',
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
                        #'Lens Position'
                        
                      ]
showTimestamp = False
# variablesToDisplay=None
textY = 1
textVA = 'bottom'

if intermediatePlot:
    # ImageAnalysisCode.ShowImagesTranspose(images_array, uniformscale=False)
    ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, 
                                          dx=dxMicron, 
                                            # filterLists=[['Evap_Time_2==1.5']],
                                           plotRate=plotRate, plotPWindow=plotPWindow,
                                            variablesToDisplay = variablesToDisplay,
                                           showTimestamp=showTimestamp,
                                          variableLog=variableLog, 
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
var2 = 'Evap_Time_2'
var1 = 'LowServo1'
df = ImageAnalysisCode.multiVariableThermometry(
                                                # results[results.Ywidth<22], 
                                                results, 
                                                
                                                var1, var2, 
                                                fitXVar='TOF',  fitYVar='Ywidth',do_plot=1, add_Text=0)

# %%
# cols = ['PSD', 'T (K)', 'AtomNum']

# for r in cols:
#     df1 = df[r].unstack()
#     fig, ax = plt.subplots(1,1, figsize=[4,3], layout='constrained')
#     cax = ax.pcolormesh(df1.columns, df1.index, df1.values, cmap='viridis')
#     ax.set(xlabel=var2, ylabel=var1, title=r)
#     fig.colorbar(cax, ax=ax)

# %%
# plt.rcParams['figure.figsize'] = [4, 3]

colnames = ['PSD', 'AtomNum', 'T (K)']

fig, axes = plt.subplots(1,3, figsize=[8, 2.5], layout='constrained')
for ii, ax in enumerate(axes): 
    df[colnames[ii]].plot(ls='',marker='x', ax=ax)
    ax.set(title=colnames[ii])
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))

fig, ax = plt.subplots(1,1, figsize=[4, 3], layout='constrained')

Amean = results[results.TOF==0].groupby([var1, var2]).mean().YatomNumber.values
Astd = results[results.TOF==0].groupby([var1, var2]).std().YatomNumber.values
ax.errorbar(np.arange(len(Amean)), Amean, Astd)
ax.set(title=colnames[1])
ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))



# fig, ax = plt.subplots(1,1, layout='constrained')
# ax.plot(df['T (K)'], df.PSD, '.')
# ax.set(xlabel='T (K)', ylabel='PSD')

# fig, ax = plt.subplots(1,1, layout='constrained')
# ax.plot(df['AtomNum'], df.PSD, '.')
# ax.set(xlabel='AtomNum', ylabel='PSD')

# %%
print('======================')
print('The phase space density is:\n{:.5e}'.format(df.PSD.iloc[0]))
