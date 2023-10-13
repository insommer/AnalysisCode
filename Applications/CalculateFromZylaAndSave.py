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

totalDataPath =r"Z:\ats317group\Data"
date = '10/12/2023'
dataCSV_filename = 'ODT Position 5_1 and _2 .csv'

dataFolders = []
dataFolders = [r'Andor/Position 5 Bias Scan_1', r'Andor/Position 5 Bias Scan_2']
saveToCSV = 1
writeToExistingFile = 1
Calculate = 1


dayFolder = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)
dataCSV_filePath = os.path.join(dayFolder, dataCSV_filename)
fileExist = os.path.exists(dataCSV_filePath)

if Calculate:    
    variableLogFolder = os.path.join(dayFolder, 'Variable Logs')
    variableLog = ImageAnalysisCode.LoadVariableLog(variableLogFolder)
    
    results = ImageAnalysisCode.CalculateFromZyla(dayFolder, dataFolders, 
                                                  variableLog = variableLog,
                                                  rowstart = 400,
                                                  rowend = 580,
                                                  columnstart = 400,
                                                  columnend = 700)
    
if saveToCSV:
    if fileExist:
        if writeToExistingFile:
            dataCSV = pd.read_csv(dataCSV_filePath)
            dataCSV.time = pd.to_datetime(dataCSV.time)
            dataCSV.set_index('time', inplace=True)
            
            intersection = dataCSV.index.intersection(results.index)
            results = pd.concat( [dataCSV.drop(intersection), results] )
        else:
            dataCSV_filename = dataCSV_filename.replace('.csv', '_1.csv')
                    
    results.to_csv( os.path.join(dayFolder, dataCSV_filename) )
    print('Results saved.')
    
#%%

ImageAnalysisCode.PlotFromDataCSV(dataCSV_filePath, 
                                  'ZSBiasCurrent', 'AtomNumber', 
                                  iterateVariable='VerticalBiasCurrent', groupbyX=1)
