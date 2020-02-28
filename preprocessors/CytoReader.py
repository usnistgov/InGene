"""Conversion of flow cytometry data into conditional output distributions for SEMIL
Author: Swarnavo Sarkar
Email: swarnavo.sarkar@nist.gov

If you are using SEMIL and any of the pre and postprocessing code, please cite 'Mutual Information Landscapes as a Performance Metric for Biochemical Reaction Networks'
"""


import glob #filenames and pathnames utility
import os   #operating sytem utility
import numpy as np
import pandas as pd
import pickle
import random as rand
import sys
import math

data_directory = '/Users/sns9/Research/IMS_project/FeedbackExpDec18/WTA'

#data_directory = '/Volumes/Shared_Data/GSF-IMS/Mutual information paper/2019-02-06-0942_IPTG_IPTG_gradient/Cytometry data'

current_dir = os.getcwd()
os.chdir(data_directory)

# Plate and duplicate label
plate_label = ['B','E']
rep_label = 'lacX'
tag = 'lacXWT_ml'
filter_string = 'lacX-'
conc_separator = '-'
plate_separator = '_'
data_fraction = 1
data_fractions = [1,2,5,10]
n_samples = list(range(1,6))

def extract_data(data_set,index_set):
    extracted_data = []

    for i in index_set:
        try:
            extracted_data.append(data_set.values[i])
        except IndexError:
            print(i)
            sys.exit()

    return np.array(extracted_data)

def compute_mean_variance(data):
    min_offset = min(data)
    shifted = data - min_offset

    mean_response = np.mean(shifted)
    var_response = np.var(shifted)

    percentiles = np.percentile(shifted,[5.0,95.0])

    return mean_response, var_response, percentiles

coli_files = glob.glob('*.frame_pkl')

filenames = [file.rsplit('.',1)[0] for file in coli_files]

coli_frame = [ pickle.load(open(file, 'rb')) for file in coli_files ]

gated_data = [frame.loc[frame['is_cell']] for frame in coli_frame]
singlet_data = [frame.loc[frame['is_singlet']] for frame in coli_frame]
all_data = [frame for frame in coli_frame]

fl_channel = 'BL1-A-MEF'
glob_min = 1000000
glob_max = 0

data_covered = []

location_string = {}
wt_string = {}

means = {}
vars = {}
percents = {}

index_set = None
data_size = 0

conclist = []
datas = {}

for i, data, gated, singlet, file in zip(range(len(all_data)), all_data, gated_data, singlet_data, coli_files):
    index_set = None
    for j in range(1):
        label, plate_no = filenames[i].split(plate_separator)
        this_plate = plate_no[0]

        if (plate_label[0] in plate_no or plate_label[1] in plate_no) and rep_label in label: # or plate_label[1] in plate_no
            print(label)

            conc_v = float(label.lstrip(filter_string))#conc_separator)[1])
            if conc_v!=0.0:
                expo = math.log(conc_v)/math.log(2.0)
                if abs(expo-int(expo))<1e-16:
                    conc_value = str(conc_v)
                else:
                    conc_value = str(conc_v*1000)

                if conc_value not in data_covered:
                    data_covered.append(conc_value)
                    conclist.append(float(conc_value))

                    datas[conc_value] = singlet[fl_channel]
                    print(len(singlet[fl_channel]))

                    glob_max = max(glob_max,max(singlet[fl_channel]))
                    glob_min = min(glob_min,min(singlet[fl_channel]))
conclist.sort()
print(conclist)

print(glob_max,glob_min)

bin_edge = np.linspace(0.0, glob_max-glob_min,500)

for k in range(0,len(bin_edge)-1):
    if k==0:
        locstring = str(0.5*(bin_edge[k]+bin_edge[k+1]))
    else:
        locstring += ','+str(0.5*(bin_edge[k]+bin_edge[k+1]))

dir_tag = '4_'+plate_label[0]+plate_label[1]

os.chdir('..')
#os.chdir(current_dir)
try:
    os.mkdir(dir_tag+tag)
except OSError:
    pass

os.chdir(dir_tag+tag)

f = open('samples.txt','w')

for c in conclist:
    if conclist.index(c)%4==0:
        print(c,file=f)

f.close()

f = open('response.csv','w')
print('i,g,+,-',file=f)

for c in conclist:
    if conclist.index(c)%4==0:
        cs = str(c)
        darray = datas[cs].values

        darray = darray - glob_min

        m = np.mean(darray)
        pc = np.percentile(darray,[5.0,95.0])

        print(cs+','+str(m)+','+str(pc[1]-m)+','+str(m-pc[0]),file=f)

f.close()

f = open('expressions.csv','w')

for c in conclist:
    if conclist.index(c)%4==0:
        cs = str(c)
        darray = datas[cs].values
        darray_list = list(darray)

        hist, b_edges = np.histogram(darray_list,bin_edge)

        total_wt = sum(list(hist))

        pdfstring = ''

        for h in hist:
            pdfstring += str(float(h)/float(total_wt))+','

        pdfstring.rstrip(',')

        print(locstring,file=f)
        print(pdfstring,file=f)

f.close()

for df in data_fractions:
    for k in n_samples:
        if df>=1:
            df_str = str(int(df))
        else:
            df_s = str(df).split('.')
            df_str = df_s[0]+'p'+df_s[1]

        f = open('expressions'+df_str+'_'+str(k)+'.csv','w')

        for c in conclist:
            if conclist.index(c)%4==0:
                cs = str(c)
                darray = datas[cs].values
                darray_list = list(darray)
                sample_size = int(len(darray_list)/df)

                d_sampled = rand.choices(darray_list,k=sample_size)

                hist, b_edges = np.histogram(np.array(d_sampled),bin_edge)

                total_wt = sum(list(hist))

                pdfstring = ''

                for h in hist:
                    pdfstring += str(float(h)/float(total_wt))+','

                pdfstring.rstrip(',')

                print(locstring,file=f)
                print(pdfstring,file=f)

        f.close()

m = 0.5*math.log10(conclist[0]) + 0.5*math.log10(conclist[-1])
v = 0.5*math.log10(conclist[0])**2 + 0.5*math.log10(conclist[-1])**2
s = v - m**2
print('Max std: ',math.sqrt(s))
