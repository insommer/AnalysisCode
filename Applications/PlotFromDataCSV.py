# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 16:23:02 2023

@author: Sommer Lab
"""

from ImageAnalysis import ImageAnalysisCode
import os
import pandas as pd

folderName = r'Z:\ats317group\Data\Analysis Resutls in csv'
fileNames = os.listdir(folderName)
# fileNames = ['ODT Position 1 Bias Scan_Oct10.csv',
#             'ODT Position 1 Bias Scan_Oct12.csv']

df = []
for f in fileNames:
    dataCSV_filePath = os.path.join(folderName, f)
    df.append(pd.read_csv(dataCSV_filePath))
    
df = pd.concat(df)



ImageAnalysisCode.PlotFromDataCSV(df, 'ZSBiasCurrent', 'AtomNumber', 
                                  iterateVariable='VerticalBiasCurrent', 
                                  filterByAnd=[
                                      'AtomNumber<1000000', 
#                                       'VerticalBiasCurrent>=7.5',
                                      'ODT_Position==5'
                                  ],
                                  groupbyX=1, threeD=0,
                                 viewElev=20, viewAzim=-20
                                 )



# df.groupby()