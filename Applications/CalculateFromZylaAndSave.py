# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 17:11:59 2023

@author: Sommer Lab
"""

from ImageAnalysis import ImageAnalysisCode
import numpy as np
# import matplotlib.pyplot as plt
from scipy.ndimage import rotate
import pandas as pd
import os
import datetime

totalDataPath = r"D:\Dropbox (Lehigh University)\Sommer Lab Shared\Data"
# totalDataPath =r"C:\Users\Sommer Lab\Documents\Data"
date = '1/8/2024'
dayFolder = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)

dataFolders = [
    # r'Andor/Bias Scan ODT at 1012',
    r'Andor/Bias Scan ODT at 1012_1'
    ]

saveToCSV = 1
writeToExistingFile = 1
Calculate = 1

targetFileName = 'Bias Scan ODT at 1012'
targetFolder = r'Z:\ats317group\Data\ODT Move Data Set 3'
targetFolder = dayFolder

print('###2')

if Calculate:    
    variableLogFolder = os.path.join(dayFolder, 'Variable Logs')
    variableLog = ImageAnalysisCode.LoadVariableLog(variableLogFolder)
    print('VariableLog loaded.')

    results = ImageAnalysisCode.CalculateFromZyla(dayFolder, dataFolders, 
                                                  variableLog = variableLog,
                                                  
                                                   rowstart = 600,
                                                   rowend = -300,
                                                   columnstart = 600,
                                                   columnend = -200,
                                                  # rowstart =560,
                                                  # rowend = -500,
                                                  # columnstart = 780,
                                                  # columnend = -550,
                                                  angle_deg= 1,
                                                  
                                                  subtract_bg=1, 
                                                  signal_feature='wide',
                                                  # signal_width=40, 
                                                  subtract_burntin=0,
                                                  plotRate=1,
                                                  plotPWindow=5,
                                                  variablesToDisplay = [
                                                       'ZSBiasCurrent',
                                                       'VerticalBiasCurrent',
                                                      # 'Coil_medB', 
                                                      # 'wait'
                                                      ]
                                                  )
    print('###4')


if Calculate and saveToCSV:
    if not os.path.exists(targetFolder):  
        os.mkdir(targetFolder)
    targetFilePath = os.path.join(targetFolder, targetFileName) + datetime.datetime.strptime(date, '%m/%d/%Y').strftime('_%b%d.csv')
    
    if os.path.exists(targetFilePath):
        if writeToExistingFile:
            dataCSV = pd.read_csv(targetFilePath)
            dataCSV.time = pd.to_datetime(dataCSV.time)
            dataCSV.set_index('time', inplace=True)
            
            intersection = dataCSV.index.intersection(results.index)
            results = pd.concat( [dataCSV.drop(intersection), results] )
        else:
            ii = 1
            targetFilePath = targetFilePath.replace('.csv', '_' + str(ii) + '.csv')
            while os.path.exists(targetFilePath):
                targetFilePath = targetFilePath.replace('_' + str(ii) + '.csv', '_' + str(ii+1) + '.csv')
                ii += 1                
                    
    results.to_csv( targetFilePath )
    print('Results saved to\n{}.'.format(targetFilePath))
    
#%%

if not Calculate:
    fileName = 'Variable Wait_1_Jan08'
    filePath = os.path.join(targetFolder, fileName) + '.csv'
    filePath = os.path.join(targetFolder, targetFileName) + datetime.datetime.strptime(date, '%m/%d/%Y').strftime('_%b%d.csv')
    results = pd.read_csv(filePath)

ImageAnalysisCode.PlotFromDataCSV(results, 
                                  xVariable='ZSBiasCurrent',
                                  yVariable='AtomNumber', 
                                  iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
                                  groupbyX=1, threeD=0,
                                  figSize = 0.5
                                  )

# fig, ax = ImageAnalysisCode.PlotFromDataCSV(results, 
#                                             xVariable='wait', 
#                                             yVariable='Xcenter', 
#                                     filterByAnd=['wait>10'], 
#                                     # filterByAnd=["Folder==r'Andor/No ODT Bias Scan_3'"],
#                                     # iterateVariable='VerticalBiasCurrent',
#                                     # groupbyX=1, threeD=0
#                                     figSize = 0.5
#                                     )
