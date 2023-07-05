# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 

@author: Sommer Lab
"""

import csv
import numpy as np 
from matplotlib import pyplot as plt
import configparser

#inputs
data_folder = r"C:\Users\Sommer Lab\Documents\Data\2022\03 Mar 2022\10 Mar"
scan_number = 1


filename = data_folder+"scan"+scan_number+"_counts.csv"
filename_cfg = filename.split('_counts.csv')[0]+".cfg"
config = configparser.ConfigParser()
config.read(filename_cfg)
print("filename = "+str(filename))


#read counts from the _counts.csv file
Pictures_Taken = int(config['Acquisition']['NumberinKineticSeries'])
time_per_picture = float(config['Acquisition']['KineticCycleTime']) 
time_axis = time_per_picture*np.linspace(1, Pictures_Taken, Pictures_Taken)

#Make the plot:
total_counts = np.zeros(Pictures_Taken)
with open(filename, "r") as f:
    reader = csv.reader(f, delimiter="\n")
    for i, row in enumerate(reader):
       total_counts[i] = row[0]

#convert number of counts to number of atoms
quantum_eff = .5
P = 27.5e-3
r_beam = .01
I1 = 2*P/(np.pi*r_beam**2)
I = 6*I1
Isat = 25 #W/m^2
s = I/Isat
gamma = 36.898e6
delta = 26e6*2*np.pi
R_scat = gamma*.5*s/(1+s+(2*delta/gamma)**2)
t_exp = .005 #s
cos_theta = 150/np.sqrt(3.175**2+150**2)
solid_angle = 2*np.pi*(1-cos_theta)
num_atoms = (4*np.pi*total_counts)/(quantum_eff*R_scat*t_exp*solid_angle)


#make the plot
plt.figure()
plt.rcParams.update({'font.size':12})
plt.title('scan'+str(scan_number))
plt.xlabel("Time (s)")
plt.ylabel("Number of atoms in MOT")
plt.plot(time_axis, num_atoms)
plt.show()
# plt.savefig(location.split('scan')[0]+'\Graphs\scan'+str(scan)+'.png', dpi = 300)





