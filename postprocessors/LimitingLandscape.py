#!/usr/bin/env python

import os
import sys

import pandas as pd
import numpy as np
import scipy.optimize as optim
from BetaDist import BetaDist
import math
import scipy.special as sp
from scipy.stats import beta
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# Enter upper and lower bounds for the probability distribution
#ub, lb = 6.0,-3.0
ub, lb = 0.4,0.0
plot_title = 'e = 500'

# Enter axes label strings
y_name = r"$\sigma(\mu)$"
x_name = r"$\langle \mu \rangle$"

# Enter tick locations
x_ticks = [0.05,0.15,0.25,0.35]
y_ticks = [0.025,0.05,0.075,0.1,0.125]

for hd in range(0,6):
    data_directory = '/Users/sns9/Research/SaraWalker_collaboration/data/mh_500/Hamming Distance  '+str(hd)+'/6landscape'

    #data_directory = '/Users/sns9/Research/IMS_project/partial_expts/2_BElacXweak2'
    #data_directory = '/Users/sns9/Research/Transcription/SED/Flow_FISH_mRNA_Rep3/landscape'
    #data_directory = '/Volumes/Shared_Data/GSF-IMS/E-Coli/pLMSF-lacI/2019-05-02-1156_IPTG_IPTG_gradient/GlacI'

    data_fractions = [1,2,5]
    n_samples = list(range(1,6))

    def create_radial_law(MI_matrix,H_matrix):
        max_H = math.log10(np.max(H_matrix)*1.05)
        min_H = -2#math.log10(np.min(H_matrix)*0.95)
        bin_size = 10#int(self.resolution/2)
        H_bins = np.linspace(min_H,max_H,bin_size+1)
        d_H = H_bins[1]-H_bins[0]
        MI_samples = {}

        for k in range(0,bin_size):
            MI_samples[k] = []

        rs, cs = MI_matrix.shape[0],MI_matrix.shape[1]

        for i in range(0,rs):
            for j in range(0,cs):
                H = H_matrix[i,j]
                MI = MI_matrix[i,j]

                if H>0.0:
                    bin_loc = max(int((math.log10(H)-min_H)/d_H),0)
                    MI_samples[bin_loc].append(MI_matrix[i,j])

        f = open('MI_rate.csv','w')
        print('H,I,+,-',file=f)
        f.close()

        for k in range(0,bin_size):
            if len(MI_samples[k])>0:
                outstring = str(10**(0.5*(H_bins[k]+H_bins[k+1])))
                this_mean = np.mean(MI_samples[k])
                outstring += ','+str(this_mean)
                outstring += ','+str(max(MI_samples[k])-this_mean)
                outstring += ','+str(this_mean-min(MI_samples[k]))

                print(outstring,file=open('MI_rate.csv','a'))
            else:
                print('0,0,0,0',file=open('MI_rate.csv','a'))

    def linear_fit(x,a,b):
        return a*x+b

    os.chdir(data_directory)

    dataframes = []
    xlist = []

    data_size = 0

    for df in data_fractions:
        if df>=1:
            df_str = str(int(df))
        else:
            df_s = str(df).split('.')
            df_str = df_s[0]+'p'+df_s[1]

        for n in n_samples:
            xlist.append(float(df))
            filename = 'trial'+df_str+'_'+str(n)+'_MI_landscape.csv'
            dataframes.append(pd.read_csv(filename))
            data_size += 1

    x = np.array(xlist)

    # data_size = 0
    # for filename in file_list:
    #     dataframes.append(pd.read_csv(filename))
    #     data_size += 1
    dataframes[-1].rename(columns={'Unnamed: 0': ''},inplace=True)

    data_shape = dataframes[-1].shape
    index_set = dataframes[-1].index
    columns_set = list(dataframes[-1].columns)
    rows_set = list(dataframes[-1][columns_set[0]])

    column_range = range(0,len(columns_set))

    output_frame = dataframes[-1].copy()
    output_pcov = dataframes[-1].copy()

    max_MI = 0
    MI_var = 0
    max_coord = [0,0]

    for i in index_set[1:]:
        for jj in column_range[1:]:
            y = []

            for k in range(0,data_size):
                j = dataframes[k].columns[jj]
                y.append(dataframes[k].at[i,j])

            popt, pcov = optim.curve_fit(linear_fit,x,np.array(y))

            if popt[1]>max_MI:
                max_MI = popt[1]
                max_coord = [float(rows_set[i]),float(j)]
                m_coord = [i,j]
                MI_var = pcov[1,1]
                #print(popt[0],popt[1])

            j = output_frame.columns[jj]
            output_frame.at[i,j] = popt[1]
            output_pcov.at[i,j] = pcov[1,1]

    #output_transpose = output_frame.T

    # of = open('c_values.csv','w')
    # kk = 0
    # for df in data_fractions:
    #     for n in n_samples:
    #         print(str(df)+','+str(dataframes[kk].at[m_coord[0],m_coord[1]]),file=of)
    #         kk += 1
    # of.close()

    # Substitute columns_set
    # for k in output_frame.index.values:
    #     output_frame.at[k,''] = 10**output_frame.at[k,'']

    output_pcov.to_csv('LimitingPcov.csv',header=False,index=False)

    #output_frame.T.to_csv('LimitingMI_surface.csv',header=False)
    output_frame.to_csv('LimitingMI_surface.csv',index=False)

    f = open('C-coord.txt','w')
    print('C = ',max_MI,file=f)
    print('C_var = ',MI_var,file=f)
    print('Mean = ',max_coord[1],file=f)
    print('Std = ',max_coord[0],file=f)

    # #Creating relative entropy landscape
    # rel_ent_frame = dataframes[-1].copy()
    # #C_beta = BetaDist(lb,ub,math.log10(max_coord[0]),max_coord[1])
    # C_beta = BetaDist(lb,ub,math.log10(max_coord[1]),max_coord[0])
    #
    # B_C = sp.gamma(C_beta.p_beta)*sp.gamma(C_beta.q_beta)/sp.gamma(C_beta.p_beta + C_beta.q_beta)
    #
    # for i in index_set[1:]:
    #     for j in columns_set[1:]:
    #         #beta_obj = BetaDist(lb,ub,math.log10(float(rows_set[i])),float(j))
    #         beta_obj = BetaDist(lb,ub,math.log10(float(j)),float(rows_set[i]))
    #
    #         this_B = sp.gamma(beta_obj.p_beta)*sp.gamma(beta_obj.q_beta)/sp.gamma(beta_obj.p_beta + beta_obj.q_beta)
    #         this_di_p = sp.digamma(beta_obj.p_beta)
    #         this_di_q = sp.digamma(beta_obj.q_beta)
    #         this_di_pq = sp.digamma(beta_obj.p_beta + beta_obj.q_beta)
    #
    #         try:
    #             v = math.log(B_C) - math.log(this_B) + (beta_obj.p_beta - C_beta.p_beta)*this_di_p + (beta_obj.q_beta - C_beta.q_beta)*this_di_q
    #             v += (C_beta.p_beta - beta_obj.p_beta + C_beta.q_beta - beta_obj.q_beta)*this_di_pq
    #             v *= 1.0/math.log(2.0)
    #         except ValueError:
    #             print(math.log10(float(rows_set[i])),float(j))
    #             sys.stdout.flush()
    #             sys.exit()
    #
    #         rel_ent_frame.at[i,j] = v
    #
    # #rel_ent_frame.T.to_csv('Hdrop.csv',header=False)
    # rel_ent_frame.to_csv('Hdrop.csv',index=False)
    #
    MI_matrix = output_frame.to_numpy()[:,1:]
    #H_matrix = rel_ent_frame.to_numpy()[:,1:]
    #create_radial_law(MI_matrix,H_matrix)

    # Write CC pdf
    m, s = math.log10(max_coord[1]), math.log10(max_coord[0])

    s_m = (m-lb)/(ub-lb)
    s_v = (s**2)/((ub-lb)**2)

    p = s_m*(s_m*(1-s_m)/s_v - 1)
    q = (s_m*(1-s_m)/s_v - 1) - p

    xx = np.linspace(0.0,1.0,1000)
    sx = xx*(ub-lb) + lb
    xpf = beta.cdf(xx,p,q,0,1)

    ff = open('cc_pdf.csv','w')

    for k in range(0,1000):
        print(str(10**sx[k])+','+str(xpf[k]),file=ff)

    ff.close()

    c_pdf = beta.pdf(xx,p,q,0,1)/(ub-lb)
    ff = open('p_c.csv','w')

    for k in range(0,1000):
        print(str(10**sx[k])+','+str(c_pdf[k]),file=ff)

    ff.close()

    # plot landscape
    out_mat = output_frame.to_numpy()
    MI_matrix = out_mat[:,1:]

    std_range = out_mat[:,0]

    mean_l = []

    for l in list(output_frame)[1:]:
        mean_l.append(float(l))

    mean_range = np.array(mean_l)

    # print(mean_range)
    # print(std_range)
    # print(MI_matrix)

    # Set up for contour
    X = np.zeros(MI_matrix.shape)
    Y = np.zeros(MI_matrix.shape)

    for i in range(0,int(len(list(mean_range)))):
        for j in range(0,int(len(list(std_range)))):
            X[i,j] = mean_range[j]
            Y[i,j] = std_range[i]

    MImax = np.max(MI_matrix)
    MImin = np.min(MI_matrix)

    levels = list(np.linspace(0,2,5))

    fig = plt.figure()
    plt.axes().set_aspect(3)

    hot_reversed = plt.cm.get_cmap('hot_r')

    C = plt.pcolor(mean_range,std_range,MI_matrix,cmap=hot_reversed,vmin=0.0,vmax=3.5)
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.title(plot_title,fontsize=14)

    #plt.ylabel(r"$\sigma(\log_{10} \mathrm{I_{ex}})$",fontsize=16)
    #plt.xlabel(r"$10^{\langle \log_{10} \mathrm{I_{ex}} \rangle}$",fontsize=16)


    plt.yticks(y_ticks)
    plt.xticks(x_ticks)

    plt.ylabel(y_name,fontsize=16)
    plt.xlabel(x_name,fontsize=16)

    #cs = plt.contour(X,Y,MI_matrix,np.array(levels[1:4]),colors='#00008b')
    #plt.clabel(cs,fontsize='large')
    plt.tick_params(which='major',length=6, width=1)
    plt.tick_params(which='minor',length=4, width=1)


    #plt.plot(max_coord[1],max_coord[0],color='white',marker='.',markersize=10)

    #MI_text = str("%.2f"%round(max_MI,2))+'('+str("%.1f"%round(max_coord[1],1))+','+str("%.2f"%round(max_coord[0],2))+')'
    #plt.text(0.5*max_coord[1],1.1*max_coord[0],MI_text,fontsize=14,color='white')

    #plt.xscale('log')
    cbar = fig.colorbar(C)

    cbar.ax.tick_params(labelsize=10)

    #cbar_fit = ax.imshow(MI_matrix,cm)

    fig.tight_layout()

    #plt.show()

    fig.savefig('Limit_landscape.pdf', bbox_inches='tight')
    fig.savefig('Limit_landscape.png', bbox_inches='tight', dpi=300)
