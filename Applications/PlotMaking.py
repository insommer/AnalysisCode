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
date = '4/17/2024'
data_folder = [
    r'Andor/ODTs',
    # r'Andor/D1 bias scan Negative Polarity', 
    # r'Andor/D1 bias scan Positive Polarity'
    ]

#%% Load the Zyla pictures and convert to columDensities
dayfolder = ImageAnalysisCode.GetDataLocation(date, DataPath=dataRootFolder)
dataPath = [ os.path.join(dayfolder, f) for f in data_folder]
params = ImageAnalysisCode.ExperimentParams(date, t_exp = 10e-6, picturesPerIteration=None, cam_type = "zyla")
dxMicron = params.camera.pixelsize_microns/params.magnification    #The length in micron that 1 pixel correspond to. 

columnDensities, variableLog = ImageAnalysisCode.PreprocessZylaImg(*dataPath, 
                                                                   examRange = ImageAnalysisCode.GetExamRange(examNum=None, examFrom=None),                                                                   
                                                                   rotateAngle=0.5, 
                                                                   rowstart = 10, 
                                                                   rowend = -10, 
                                                                   columnstart = 800,
                                                                   columnend = 1100,
                                                                   subtract_burntin=0, 
                                                                   showRawImgs=0)
#%% Fit the columdensities and analyze the results
popts, bgs = ImageAnalysisCode.FitColumnDensity(columnDensities, dx = dxMicron, mode='both', yFitMode='single',
                                                subtract_bg=0, Xsignal_feature='wide', Ysignal_feature='narrow')

resultsODT = ImageAnalysisCode.AnalyseFittingResults(popts, logTime=variableLog.index) 

resultsODT = resultsODT.join(variableLog)
# results.to_csv('Test.csv')

ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, 
                                      dx=dxMicron, 
                                        # filterLists=[['fmod_kHz>32']],
                                       plotRate=1, plotPWindow=6,
                                      # variablesToDisplay = variablesToDisplay,
                                      variableLog=variableLog, 
                                      logTime=variableLog.index,
                                      textLocationY=0.8, rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9})
#%%
####################################
#Set the date and the folder name
####################################
dataRootFolder = r"D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data"
date = '4/17/2024'
data_folder = [
    # r'Andor/ODTs',
    r'Andor/D1 bias scan Negative Polarity', 
    r'Andor/D1 bias scan Positive Polarity'
    ]

#%% Load the Zyla pictures and convert to columDensities
dayfolder = ImageAnalysisCode.GetDataLocation(date, DataPath=dataRootFolder)
dataPath = [ os.path.join(dayfolder, f) for f in data_folder]
params = ImageAnalysisCode.ExperimentParams(date, t_exp = 10e-6, picturesPerIteration=None, cam_type = "zyla")
dxMicron = params.camera.pixelsize_microns/params.magnification    #The length in micron that 1 pixel correspond to. 

columnDensities, variableLog = ImageAnalysisCode.PreprocessZylaImg(*dataPath, 
                                                                   examRange = ImageAnalysisCode.GetExamRange(examNum=None, examFrom=None),                                                                   
                                                                   rotateAngle=0.5, 
                                                                   rowstart = 10, 
                                                                   rowend = -10, 
                                                                   columnstart = 10,
                                                                   columnend = -10,
                                                                   subtract_burntin=1, 
                                                                   showRawImgs=0)
#%% Fit the columdensities and analyze the results
popts, bgs = ImageAnalysisCode.FitColumnDensity(columnDensities, dx = dxMicron, mode='both', yFitMode='single',
                                                subtract_bg=0, Xsignal_feature='wide', Ysignal_feature='narrow')

resultsGM = ImageAnalysisCode.AnalyseFittingResults(popts, logTime=variableLog.index) 

if variableLog is not None:
    resultsGM = resultsGM.join(variableLog)
# results.to_csv('Test.csv')

ImageAnalysisCode.plotImgAndFitResult(columnDensities, popts, bgs=bgs, 
                                      dx=dxMicron, 
                                        # filterLists=[['fmod_kHz>32']],
                                       plotRate=1, plotPWindow=6,
                                      # variablesToDisplay = variablesToDisplay,
                                      variableLog=variableLog, 
                                      logTime=variableLog.index,
                                      textLocationY=0.8, rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9})

resultsGM.YatomNumber = resultsGM.YatomNumber  / 1.4

# %%
plt.rcParams['font.size'] = 15
plt.rcParams['xtick.labelsize'] = 13
plt.rcParams['ytick.labelsize'] = 13

fig1, ax1 = ImageAnalysisCode.PlotFromDataCSV(resultsGM, 'Ycenter', 'YatomNumber', 
                                  # iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                  groupby='ODT_Position', 
                                  threeD=0,
                                  figSize = 0.5
                                  )

fig2, ax2 = ImageAnalysisCode.PlotFromDataCSV(resultsODT, 'Ycenter', 'YatomNumber', 
                                  # iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                  groupby='ODT_Position', 
                                  threeD=0,
                                  figSize = 0.5
                                  )

# %%
for ax in [ax1, ax2]:
    ax.set(xlabel = 'Vertical Position (µm)', ylabel = 'Atom Number')
ax1.set_title('Gray Molasses')
ax2.set_title('ODT')
# %%


ODTmean = resultsODT.groupby('ODT_Position').mean()
ODTsem = resultsODT.groupby('ODT_Position').sem()
GMmean = resultsGM.groupby('ODT_Position').mean()
GMsem = resultsGM.groupby('ODT_Position').sem()

x = ODTmean.Ycenter
x1 = GMmean.Ycenter
A = ODTmean.YatomNumber
B = GMmean.YatomNumber
errorA = ODTsem.YatomNumber
errorB = GMsem.YatomNumber

ratio = A / B
error = ratio * ( (errorA/A)**2 + (errorB/B)**2 )**0.5

fig1, ax1 = plt.subplots(1,1, figsize=(5,4), layout='constrained')
ax1.errorbar(x, A, errorA)

fig2, ax2 = plt.subplots(1,1, figsize=(5,4), layout='constrained')
ax2.errorbar(x1, B, errorB)

for ax in [ax1, ax2]:
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    ax.set(xlabel = 'Vertical Position (µm)', ylabel = 'Atom Number')
    
ax1.set_title('ODT')
ax2.set_title('Gray Molasses')

fig, ax = plt.subplots(1,1, figsize=(5,4), layout='constrained')
ax.errorbar(x, ratio, error)
ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
ax.set(xlabel = 'Vertical Position (µm)', ylabel = 'Atom Number Ratio', ylim=(0,None))


# %%

ax2.set_xticks( ax2.get_xticks()[1:-1:2] )

# import matplotlib.dates as mdates

# fig, ax = plt.subplots(1,1, figsize=(8,6), layout='constrained')
# # plt.plot(results.wait, results.YatomNumber)


# ax.plot(results.index, results.Ycenter.values, '.')
# ax.set(xlabel='time (Day HH:MM)', ylabel='y center (µm)')
# ax.set_xticks(ax.get_xticks()[::2])
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))
