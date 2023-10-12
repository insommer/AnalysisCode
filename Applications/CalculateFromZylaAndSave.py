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
date = '10/10/2023'
dataCSV_filename = 'Position 1 Bias Scan_no redundant index.csv'
dataFolders = []
dataFolders = [r'Andor/Position 1 Bias Scan_2']
saveToCSV = 1
writeToExistingFile = 1
Calculate = 1


dayFolderPath = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)
dataCSV_filePath = os.path.join(dayFolderPath, dataCSV_filename)
fileExist = os.path.exists(dataCSV_filePath)

if Calculate:    
    variableLog_folder = os.path.join(dayFolderPath, 'Variable Logs')
    variableLog = ImageAnalysisCode.LoadVariableLog(variableLog_folder)
    
    results = ImageAnalysisCode.CalculateFromZyla(dayFolderPath, dataFolders, 
                                                  variableLog = variableLog)
    
    
if fileExist:
    if saveToCSV:
        if writeToExistingFile:
            dataCSV = pd.read_csv(dataCSV_filePath)
            dataCSV.time = pd.to_datetime(dataCSV.time)
            dataCSV.set_index('time', inplace=True)
            
            intersection = dataCSV.index.intersection(results.index)
            dataCSV = pd.concat( [dataCSV.drop(intersection), results] )
            
        else:
            dataCSV_filename = dataCSV_filename.replace('.csv', '_1.csv')
            dataCSV = results
        
        dataCSV.to_csv( os.path.join(dayFolderPath, dataCSV_filename) )



#%% Plot

ImageAnalysisCode.PlotFromDataCSV(dataCSV_filePath, 
                                  'ZSBiasCurrent', 'AtomNumber', 
                                  iterateVariable='VerticalBiasCurrent', groupbyX=1,
                                  filterlist=['wait==30'])

    
