#!/usr/bin/env python

import os
import sys

import pandas as pd
import numpy as np
import scipy.optimize as optim

def linear_fit(x,a,b):
    return a*x+b

data_directory = '/User/sns9/Research/Transcription/Microscopy_FISH_mRNA-counts_Rep1/landscape'

xlist = [1,2,5,10,20]
sets = list(range(1,6))

os.chdir(data_directory)

dataframes = []
data_size = 0

x = []

for i in xlist:
    x.append(i)
    for s in sets:
        filename = 'trial'+str(int(i))+'_'+str(int(s))+'_MI_landscape.csv'
        dataframes.append(pd.read_csv(filename))
        data_size += 1

x = np.array(x)

data_shape = dataframes[-1].shape
index_set = dataframes[-1].index
columns_set = list(dataframes[-1].columns)
#print(dataframes[-1].index[1:])

output_frame = dataframes[-1].copy()

for i in index_set[1:]:
    for j in columns_set[1:]:
        y = []

        for k in range(0,data_size):
            y.append(dataframes[k].at[i,j])

        popt, pcov = optim.curve_fit(linear_fit,x,np.array(y))

        output_frame.at[i,j] = popt[1]


output_frame.to_csv('LimitingMI_surface.csv')
