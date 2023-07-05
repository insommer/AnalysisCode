# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 10:14:33 2022

@author: Sommer Lab
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import csv

#df = pd.read_csv(r'C:\Users\Sommer Lab\Documents\Data\2022\01 Jan 2022\31 Jan\scan1_counts.csv')
#print(df)



#filename = r'C:\Users\Sommer Lab\Documents\Data\2022\01 Jan 2022\31 Jan\scan1_counts.csv'


#content = np.loadtxt(filename)
#print(content)


#file = open(r'C:\Users\Sommer Lab\Documents\Data\2022\01 Jan 2022\31 Jan\scan1_counts.csv',"rb")
#content = file.read()
#print(content)


num_frames = 100

start_index = 23
end_index = 29
num_csvfiles = end_index-start_index +1
location = r"C:\Users\Sommer Lab\Documents\Data\2022\01 Jan 2022\31 Jan\scan"
matrix = np.zeros((num_csvfiles,num_frames))
matrix[:,:] = np.NaN


for matrix_index,i in enumerate(range(start_index, end_index + 1)):
    file = location+str(i)+"_counts.csv"
    content = np.loadtxt(file)
    content = content - np.min(content)
    content = content / np.max(content)
    matrix[matrix_index,:] = content
 
#matrix[0,:] = 0    
 
print(matrix)         


plt.imshow(matrix, cmap='jet', vmax=np.max(matrix), aspect = 10)
plt.colorbar()
plt.show()


# matrix = np.zeros((num_frames,num_csvfiles))




# #scans = ["scan1", "scan2", "scan3", "scan4", "scan5", "scan6", "scan7", "scan8", "scan9", "scan10", "scan11"]



# #filename = r"C:\Users\Sommer Lab\Documents\Data\2022\01 Jan 2022\25 Jan\scan1_counts.csv"


# matrix_tmp = np.array([])

# for i in range(num_csvfiles):
#         file = location+str(i+1)+"_counts.csv"
#         content = pd.read_csv(file)
#         len_tmp = len(content)
#         print("there are:", len_tmp, "entries in this csv file")
#         matrix_tmp.append(content)
#         #np.frombuffer(content,dtype=int)


# print(matrix_tmp)




#file_names = [directory_name.split(".")[0]+str(x)+"_counts.csv" for x in scans]
