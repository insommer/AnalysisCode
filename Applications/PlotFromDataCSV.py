# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:14:44 2023

@author: insommer
"""

import sys
# sys.path.append(r'C:\Users\Sommer Lab\Documents\Analysis Code')
import os

# from ImageAnalysis import ImageAnalysisCode
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime
from ImageAnalysis import ImageAnalysisCode


totalDataPath =r"Z:\ats317group\Data"
date = '9/29/2023'
dataCSV_filename = 'Data Folder 5-9.csv'
iterateVariableName = 'VerticalBiasCurrent'
xVariable = 'Ycenter'
yVariable = 'AtomNumber'

dayFolderPath = ImageAnalysisCode.GetDataLocation(date, DataPath=totalDataPath)
dataCSV_filePath = os.path.join(dayFolderPath, dataCSV_filename)

if not os.path.exists(dataCSV_filePath):
    raise FileNotFoundError("The file does not exist!")

df = pd.read_csv(dataCSV_filePath)

#Filter the data as needed. 
df = df[ (df.wait==30) & (df.Ywidth < 0.00006) ]

iterateVariable = df[iterateVariableName]
iterateVariable = iterateVariable.unique()

fig, ax = plt.subplots(figsize=(8,5))
for ii in (iterateVariable):
    dfSelect = df[ (df.VerticalBiasCurrent==ii) ] 
    plt.plot( dfSelect[xVariable], dfSelect[yVariable], '.', label = '{} = {}'.format(iterateVariableName, ii))
plt.legend()
ax.set(xlabel=xVariable, ylabel=yVariable)
fig.tight_layout()
plt.show()