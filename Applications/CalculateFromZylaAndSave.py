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
# totalDataPath =r"C:\Users\Sommer Lab\Documents\Data"
date = '11/30/2023'

dataFolders = [
    r'Andor/ODT Position 1200 Bias Scan CMOT',
    ]

saveToCSV = 1
writeToExistingFile = 1
Calculate = 1

targetFileName = 'ODT Position 1200 Bias Scan CMOT'
targetFolder = r'Z:\ats317group\Data\ODT Move Data Set 3'

print('###2')

dayFolder = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)

if Calculate:    
    variableLogFolder = os.path.join(dayFolder, 'Variable Logs')
    variableLog = ImageAnalysisCode.LoadVariableLog(variableLogFolder)
    print('VariableLog loaded.')

    results = ImageAnalysisCode.CalculateFromZyla(dayFolder, dataFolders, 
                                                  variableLog = variableLog,
                                                  
                                                  # rowstart = 400,
                                                  # rowend = 650,
                                                  # columnstart = 400,
                                                  # columnend = 700,
                                                  
                                                    # rowstart = 400,
                                                    # rowend = -350,
                                                    # columnstart = 600,
                                                    # columnend = -670,
                                                  
                                                  subtract_bg=0, 
                                                  signal_width=40, 
                                                  subtract_burntin=1,
                                                  plotRate=0.2,
                                                  plotPWindow=4,
                                                  variablesToDisplay = ['ZSBiasCurrent',
                                                                        'VerticalBiasCurrent']                                                  
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
    fileName = 'No ODT Bias Scan_Oct24'
    filePath = os.path.join(targetFolder, fileName) + '.csv'
    results = pd.read_csv(filePath)

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
