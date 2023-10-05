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
date = '9/26/2023'

dataCSV_filename = 'data_10_5.csv'
dayFolderPath = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)
dataFolders = [r'Andor/Test', r'Andor/Test_1']

dataCSV_filePath = os.path.join(dayFolderPath, dataCSV_filename)


if os.path.exists(dataCSV_filePath):
    dataCSV = pd.read_csv(dataCSV_filePath)
    dataCSV.time = pd.to_datetime(dataCSV.time)
    dataCSV.set_index('time', inplace=True)
else:
    variableLog_folder = os.path.join(dayFolderPath, 'Variable Logs')
    dataCSV = ImageAnalysisCode.LoadVariableLog(variableLog_folder)
    
    
results = ImageAnalysisCode.CalculateFromZyla(dayFolderPath, 
                dataFolders, 
                variableLog = dataCSV,                
                repetition=1)

combined = results.join( dataCSV, how='outer', lsuffix='_l' )

combined.to_csv( os.path.join(dayFolderPath, dataCSV_filename) )

    
