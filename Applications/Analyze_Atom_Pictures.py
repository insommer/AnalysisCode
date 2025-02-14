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
dataRootFolder = '../Test Data'

date = '11/11/2024'

data_folder = [
    r'Test top',
    ]
####################################
#Parameter Setting'
####################################
cameras = [
    'zyla',
    'chameleon'
]

reanalyze = 1
saveresults = 0
overwriteOldResults = 1

examNum = 3 #The number of runs to exam.
examFrom = None #Set to None if you want to check the last several runs. 
autoCrop = 0
showRawImgs = 0

# in the format of [zyla, chameleon]
runParams = {
    'subtract_burntin': [0, 0],
    'skip_first_img': ['auto', 0],
    'rotate_angle': [0, 0], #rotates ccw
    'ROI': [
        # [10, -10, 10, -10],
        [300, 430, 650, 1000],
        [10, -10, 10, -10]
    ], # rowStart, rowEnd, colStart, colEnd
    
    'subtract_bg': [1, 0], 
    'y_feature': ['wide', 'wide'], 
    'x_feature': ['wide', 'wide'], 
    'y_peak_width': [10, 10], # The narrower the signal, the bigger the number.
    'x_peak_width': [10, 10], # The narrower the signal, the bigger the number.
    'fitbgDeg': [5, 5],
    
    'optical_path': ['side', 'top']
}

# Set filters for the data, NO SPACE around the operator.
variableFilterList = [
    # [# 'wait==50', 
    # # 'VerticalBiasCurrent==0'
    # 'fmod_kHz==0',
    # # 'Evap_Tau==0.1',
    # # 'Evap_Time_1==2'
    # ], 
    ] 

####################################
dayfolder = ImageAnalysisCode.GetDayFolder(date, root=dataRootFolder)
paths_zyl = [ os.path.join(dayfolder, 'Andor', f) for f in data_folder]
paths_cha = [ os.path.join(dayfolder, 'FLIR', f) for f in data_folder]
runParams['paths'] = [paths_zyl, paths_cha]

runParams['expmntParams'] = np.vectorize(ImageAnalysisCode.ExperimentParams)(
    date, axis=runParams['optical_path'], cam_type=cameras)

runParams['dx_micron'] = np.vectorize(lambda a: a.camera.pixelsize_microns / a.magnification)(runParams['expmntParams'])

runParams = pd.DataFrame.from_dict(runParams, orient='index', columns=['zyla', 'chameleon'])
examRange = ImageAnalysisCode.GetExamRange(examNum, examFrom)

####################################
####################################


#%%
# if not reanalyze:
#     resultsList = []
#     for pp in dataPath:
#         resutlsPath = os.path.join(pp, 'results.pkl')        
#         if os.path.exists(resutlsPath):
#             with open(resutlsPath, 'rb') as f:
#                 resultsList.append( pickle.load(f) )



#%%
OD = {}
varLog = {}
fits = {}
results = {}

for cam in cameras:
    params = runParams[cam]

    OD[cam], varLog[cam] = ImageAnalysisCode.PreprocessBinImgs(*params.paths, camera=cam, examRange=examRange,
                                                     rotateAngle=params.rotate_angle, 
                                                               ROI=params.ROI,
                                                      subtract_burntin=params.subtract_burntin, 
                                                      skipFirstImg=params.skip_first_img,
                                                      showRawImgs=showRawImgs, 
                                                      #!!!!!!!!!!!!!!!!!
                                                      #! Keep rebuildCatalogue = 0 unless necessary!
                                                      rebuildCatalogue=0,
                                                      ##################
                                                      # filterLists=[['TOF<1']]
                                                     )

    if autoCrop:
        OD[cam] = ImageAnalysisCode.AutoCrop(OD[cam], sizes=[120, 70])
        print('opticalDensity auto cropped.')

    # columnDensities[cam] = OD[cam] / params.expmntParams.cross_section
    # popts[cam], bgs[cam]
    fits[cam] = ImageAnalysisCode.FitColumnDensity(OD[cam]/params.expmntParams.cross_section, 
                                                    dx = params.dx_micron, mode='both', yFitMode='single',
                                                    subtract_bg=params.subtract_bg, Xsignal_feature=params.x_feature, 
                                                              Ysignal_feature=params.y_feature)

    results[cam] = ImageAnalysisCode.AnalyseFittingResults(fits[cam][0], logTime=varLog[cam].index)
    results[cam] = results[cam].join(varLog[cam])

    if saveresults:
        ImageAnalysisCode.SaveResultsDftoEachFolder(results[cam], overwrite=overwriteOldResults)    


# %%
# ImageAnalysisCode.PlotResults(results, 'HF_AOM_Freq', 'Ywidth', 
#                                    iterateVariable='Lens_Position', 
#                                    # filterLists=[['Ywidth<10']],
#                                   groupbyX=1, threeD=0,
#                                   figSize = 0.5
#                                   )

for cam in cameras:
    ImageAnalysisCode.PlotResults(results[cam], 'fmod_kHz', 'Ywidth', 
                                  # iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                  # groupby='ODT_Position', 
                                    groupbyX=1, 
                                  threeD=0,
                                  figSize = 0.5
                                  )

# fig, ax = plt.subplots(figsize=(5,4), layout='constrained') 
# results[cam].Xcenter.plot(title='x Position', linestyle='', marker='.')

# %%

    intermediatePlot = 1
    plotPWindow = 4
    plotRate = 1
    uniformscale = 0
    rcParams = {'font.size': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
                # 'image.interpolation': 'nearest'
                }

    variablesToDisplay = [
                        # # 'Coil_medB', 
                            # 'TOF',
                            # # 'ODT_Misalign',
                            'wait',
                            'fmod_kHz',
                            # 'Lens_Position'

                          ]
    showTimestamp = True
    textY = 1
    textVA = 'bottom'

    if intermediatePlot:
        ImageAnalysisCode.plotImgAndFitResult(OD[cam]/runParams[cam].expmntParams.cross_section, 
                                              fits[cam][0], bgs=fits[cam][1], 
                                              dx=runParams[cam].dx_micron, 
                                              imgs2=OD[cam], 
                                                # filterLists=[['LowServo1==0.5']],
                                               plotRate=plotRate, plotPWindow=plotPWindow,
                                                variablesToDisplay = variablesToDisplay,
                                               showTimestamp=showTimestamp,
                                              variableLog=varLog[cam], 
                                              logTime=varLog[cam].index,
                                              uniformscale=uniformscale,
                                              textLocationY=0.9, rcParams=rcParams,
                                              figSizeRate=1, 
                                              sharey=False
                                             )



