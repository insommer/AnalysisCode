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
dataCSV_filename = 'Data Folder 5-9.csv'
dataFolders = [r'Andor/ODT position 5', r'Andor/ODT position 6',
               r'Andor/ODT position 7', r'Andor/ODT position 8', 
               r'Andor/ODT position 9']
saveToCSV = 1
writeToExistingFile = 0

dayFolderPath = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)
dataCSV_filePath = os.path.join(dayFolderPath, dataCSV_filename)
fileExist = os.path.exists(dataCSV_filePath)

if fileExist:
    if writeToExistingFile:
        dataCSV = pd.read_csv(dataCSV_filePath)
        dataCSV.time = pd.to_datetime(dataCSV.time)
        dataCSV.set_index('time', inplace=True)
    else:
        raise FileExistsError("The file already exists. Change the filename if you don't want to overite the older file.")

else:
    variableLog_folder = os.path.join(dayFolderPath, 'Variable Logs')
    dataCSV = ImageAnalysisCode.LoadVariableLog(variableLog_folder)

print('Variable Log Loaded.')
    
results = ImageAnalysisCode.CalculateFromZyla(dayFolderPath, dataFolders, 
                                              variableLog = dataCSV, repetition=1)

# if fileExist:
    
# else:
dataCSV = results.join( dataCSV, how='outer', lsuffix='_L' )
    

if saveToCSV:
    dataCSV.to_csv( os.path.join(dayFolderPath, dataCSV_filename) )

    
