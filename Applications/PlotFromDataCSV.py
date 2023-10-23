# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 16:23:02 2023

@author: Sommer Lab
"""

from ImageAnalysis import ImageAnalysisCode
import os
import pandas as pd

totalDataPath =r"Z:\ats317group\Data"
date = '10/17/2023'
dataCSV_filename = 'No ODT Bias Scan.csv'

dayFolder = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)
dataCSV_filePath = os.path.join(dayFolder, dataCSV_filename)
df = pd.read_csv(dataCSV_filePath)

ImageAnalysisCode.PlotFromDataCSV(df, 'ZSBiasCurrent', 'AtomNumber', 
                                  iterateVariable='VerticalBiasCurrent', 
                                  # filterByAnd=['AtomNumber<1000000', 'VerticalBiasCurrent>=7.5'],
                                  groupbyX=1, threeD=1)