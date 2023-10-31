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

totalDataPath =r"Z:\ats317group\Data"
date = '10/24/2023'
dataFolders = [r'Andor/No ODT Bias Scan_1']


saveToCSV = 1
writeToExistingFile = 1
Calculate = 1

targetFileName = 'No ODT Bias Scan_CalculatedOnOCT31'
targetFolder = r'Z:\ats317group\Data\Analysis Resutls in csv'


dayFolder = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)

if Calculate:    
    variableLogFolder = os.path.join(dayFolder, 'Variable Logs')
    variableLog = ImageAnalysisCode.LoadVariableLog(variableLogFolder)
    
    results = ImageAnalysisCode.CalculateFromZyla(dayFolder, dataFolders, 
                                                  variableLog = variableLog,
                                                  # rowstart = 400,
                                                  # rowend = 650,
                                                  # columnstart = 400,
                                                  # columnend = 700,
                                                  subtract_bg=False, 
                                                  signal_width=40
                                                  )
    
targetFilePath = os.path.join(targetFolder, targetFileName) + datetime.datetime.strptime(date, '%m/%d/%Y').strftime('_%b%d.csv')

if Calculate and saveToCSV:
    
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
    print('Results saved.')
    
#%%

if not Calculate:
    results = pd.read_csv(targetFilePath)

# ImageAnalysisCode.PlotFromDataCSV(results, 
#                                   'Ycenter', 'AtomNumber', 
#                                   iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   groupbyX=1, threeD=1)

fig, ax = ImageAnalysisCode.PlotFromDataCSV(results, 'ZSBiasCurrent', 'AtomNumber', 
                                    # filterByAnd=['wait==30', 'AtomNumber>1e4'], 
                                    # filterByAnd=["Folder==r'Andor/No ODT Bias Scan_3'"],
                                    iterateVariable='VerticalBiasCurrent',
                                    groupbyX=1, threeD=0
                                    )
