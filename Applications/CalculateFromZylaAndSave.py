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
date = '9/29/2023'
dataCSV_filename = 'Data Folder 5-9_Calculated on 1019.csv'

dataFolders = []
dataFolders = [r'Andor/ODT position 9']
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
                                                  rowend = 650,
                                                  columnstart = 400,
                                                  columnend = 700,
                                                  subtract_bg=True, 
                                                  signal_width=40)
    
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

if not Calculate:
    result = pd.read_csv(dataCSV_filePath)

# ImageAnalysisCode.PlotFromDataCSV(dataCSV_filePath, 
#                                   'Ycenter', 'AtomNumber', 
#                                   iterateVariable='VerticalBiasCurrent', 
#                                   # filterByAnd=['VerticalBiasCurrent>7.6', 'VerticalBiasCurrent<8'],
#                                   groupbyX=1, threeD=1)

fig, ax = ImageAnalysisCode.PlotFromDataCSV(results, 'Ycenter', 'AtomNumber', 
                                   filterByAnd=['wait==30', 'AtomNumber>1e4'], 
                iterateVariable='VerticalBiasCurrent', groupby='ODT_Position')
