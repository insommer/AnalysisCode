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
date = '9/3/2024'
# date = '11/07/2023'

data_folder = [
    # r'Andor/Imaging Frequency Scan_6',
    # r'Andor/Imaging Frequency Scan_4',

    # r'Andor/Xcenter vs wait'
    r'Andor/MF probe evap',
    # r'Andor/Probe evap',
    # r'Andor/ODT -1 Align',
    # r'Andor\gray molasses more D1 power_1'
    # r'Andor\With Low Servo Late ODT LowServo1 0.29 Changed Lens_1',
    # r'Andor\before evap thermometry low field Digital Mag Off',
    # r'Andor\GM Temp Vary D1 Re Attn_3',
    # r'Andor\RF Freq Scan 225 - 235 RF for 5 ms'
    # r'Andor\RF test with function output',
    # r'Andor\RF Freq Scan ZSBiasCurrent 0.7_1',
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


columnstart=600
columnend= 1100

# rowstart = 800
# rowend = -900

rowstart = 730
rowend = 1070

# # second pass
# rowstart = 1330
# rowend = 1375

# # first pass
# rowstart = 1365
# rowend = 1400

# rowstart = 500
# rowend = -500

# rowstart -= 40
# rowstart += 80
# rowend += 80

rowstart += 150
rowend -= 110   

columnstart += 140
columnend -= 200

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
                                                                   showRawImgs=showRawImgs, rebuildCatalogue=0,
                                                                   # filterLists=[['Evap_timestep>0.2']]
                                                                   )

autoCrop = 0
if autoCrop:
    columnDensities = ImageAnalysisCode.AutoCrop(columnDensities, sizes=[120, 80])
    print('ColumnDensities auto cropped.')
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

# ImageAnalysisCode.PlotFromDataCSV(results, 'RF_FRQ_MHz', 'YatomNumber', 
# #                                   # iterateVariable='VerticalBiasCurrent', 
# #                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   groupbyX=1, threeD=0,
#                                   figSize = 0.5
#                                   )
# ImageAnalysisCode.PlotFromDataCSV(results, 'LF_AOM_freq', 'Ywidth', 
#                                    iterateVariable='Lens_Position', 
#                                    # filterLists=[['wait>=10', 'wait<=100']],
#                                    groupbyX=1, threeD=0,
#                                    figSize = 0.5, do_fit=0
#                                    )

# ImageAnalysisCode.PlotFromDataCSV(results, 'wait', 'Xcenter', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   # groupby='ODT_Position', 
#                                     groupbyX=1, 
#                                   threeD=0,
#                                   figSize = 0.5
#                                   )


# ImageAnalysisCode.PlotFromDataCSV(results, 'CamBiasCurrent', 'YatomNumber', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   # groupby='ODT_Position', 
#                                     groupbyX=1, 
#                                   threeD=0,
#                                   figSize = 0.5
#                                   )

# ImageAnalysisCode.PlotFromDataCSV(results, 'wait', 'Xcenter', 
#                                   # iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   # groupby='ODT_Position', 
#                                     groupbyX=1, 
#                                   threeD=0,
#                                   figSize = 0.5
#                                   )

# fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
# results.YatomNumber.plot(title='Atom Number', linestyle='', marker='.')

# fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
# results.Ycenter.plot(title='y Position', linestyle='', marker='.')

# fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
# results.Xcenter.plot(title='x Position', linestyle='', marker='.')

# %%

intermediatePlot = 1
plotPWindow = 4
plotRate = 1
uniformscale = 0
rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9}

variablesToDisplay = [
                    # # 'Coil_medB', 
                        'TOF',
                        # 'wait',
                        # 'LF_AOM_freq',
                        # 'Lens_Position',
                        # 'FB Voltage',
                        # 'B_Field',
                        # 'ODT_Position',
                        # 'fmod_kHz',
                        # 'tmod_ms',
                        # 'Evap_Tau',
                        # 'VerticalBiasCurrent',
                        # 'B_spikeTime',
                        # 'HF_AOM_Freq',
                        # 'CamBiasCurrent',
                        #'Lens Position',
                        
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
                                          figSizeRate=1, sharey='col')
    
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


#%% LIFETIME MEASUREMENT

# popt,_ = ImageAnalysisCode.fit_exponential(results['wait'], results['YatomNumber'],
#     dx=1, doplot = True, label="", title="Trap Lifetime", newfig=True, xlabel="Wait Time (ms)", ylabel="Y Atom Number", 
#     offset = None, 
#     legend=True)

#%% OSCILLATION OF CLOUD

# dfmean = results.groupby('wait')[['Xcenter', 'Ycenter']].mean().reset_index()
# dfstd = results.groupby('wait')[['Xcenter', 'Ycenter']].std().reset_index()

# plt.figure()
# plt.errorbar(dfmean['wait'], dfmean['Ycenter'], yerr=dfstd['Ycenter'], fmt='-o')
# plt.ylabel('Ycenter')
# plt.xlabel('wait')
# plt.tight_layout()

# plt.figure()
# plt.errorbar(dfmean['wait'], dfmean['Xcenter'], yerr=dfstd['Xcenter'], fmt='-o')
# plt.ylabel('Xcenter')
# plt.xlabel('wait')
# plt.tight_layout()

# %% THERMOMETRY

# # var2 = 'D1_Re_Attn'
var2 = 'Evap_timestep'
# var2 = 'LowServo1'
df = ImageAnalysisCode.multiVariableThermometry(
                                                results, 
                                                # var1, 
                                                var2, 
                                                fitXVar='TOF',  fitYVar='Ywidth',do_plot=1, add_Text=1)


#%% ASPECT RATIO CALCULATION
# width_mean = results.groupby('TOF')[['Xwidth', 'Ywidth']].mean().reset_index()
# width_std = results.groupby('TOF')[['Xwidth', 'Ywidth']].std().reset_index()


# aspectRatio = width_mean['Ywidth'] / width_mean['Xwidth']
# # aspectRatio_error = aspectRatio * np.sqrt( (width_std['Xwidth']/width_std['Xwidth'])**2 + (width_std['Ywidth']/width_std['Ywidth'])**2)

# plt.figure()
# plt.plot(width_mean['TOF'], aspectRatio, '-o')
# plt.xlabel('Time-of-flight (ms)')
# plt.ylabel('Aspect Ratio')

# %% 2-D plot when have two variable parameters
# cols = ['PSD', 'T (K)', 'AtomNum']

# for r in cols:
#     df1 = df[r].unstack()
#     fig, ax = plt.subplots(1,1, figsize=[4,3], layout='constrained')
#     cax = ax.pcolormesh(df1.columns, df1.index, df1.values, cmap='viridis')
#     ax.set(xlabel=var2, ylabel=var1, title=r)
#     fig.colorbar(cax, ax=ax)

# %% 1-D plot when only vary one parameter
plt.rcParams['figure.figsize'] = [4, 3]

colnames = ['PSD', 'AtomNum', 'T (K)']

fig, axes = plt.subplots(1,3, figsize=[8, 2.5], layout='constrained')
for ii, ax in enumerate(axes): 
    df[colnames[ii]].plot(ls='',marker='x', ax=ax)
    ax.set(title=colnames[ii], )
    # ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    ax.set_yscale('log')

# fig, ax = plt.subplots(1,1, figsize=[4, 3], layout='constrained')

# Amean = results[results.TOF==0].groupby([var1, var2]).mean().YatomNumber.values
# Astd = results[results.TOF==0].groupby([var1, var2]).std().YatomNumber.values
# ax.errorbar(np.arange(len(Amean)), Amean, Astd, marker='.')
# ax.set(title=colnames[1])
# ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))



# fig, ax = plt.subplots(1,1, layout='constrained')
# ax.plot(df.AtomNum, df.PSD, '.')
# ax.set(xlabel='AtomNum', ylabel='PSD')
# ax.set_yscale('log')
# ax.set_xscale('log')



# # fig, ax = plt.subplots(1,1, layout='constrained')
# # ax.plot(df['AtomNum'], df.PSD, '.')
# # ax.set(xlabel='AtomNum', ylabel='PSD')

# # %%
# print('======================')
# print('The phase space density is:\n{}'.format(df[['AtomNum', 'T (K)']]))
