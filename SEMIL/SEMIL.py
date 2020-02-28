"""Code to compute Mutual Information landscape from gene expression data
"""
import numpy as np
import random as rand
import os
import math
import sys
import time
from MakeInputDict import MakeInputDict
from MIcalculator import MIcalculator
from SROMgenerator import SROMgenerator
from scipy.optimize import minimize
from scipy import optimize
import matplotlib.pyplot as plt
import scipy.special as sp
from BetaDist import BetaDist
from scipy.stats import beta
import matplotlib.colors as colors
import copy

class SEMIL(MakeInputDict):
    def __init__(self,data_path,input_lines,_flag):
        self.SROM_obj = SROMgenerator(data_path,input_lines)

        response_file = self.SROM_obj.input_folder+'/'+self.SROM_obj.input_dict['Response file']
        self.MI_instance = MIcalculator(response_file,len(self.SROM_obj.sample_locations))
        #self.MI_instance.setloc(self.SROM_obj.sample_locations)
        #self.MI_instance.read_responses(response_file)
        self.flag = _flag
        self.casename = self.SROM_obj.input_dict['Casename'].rstrip('\r\n')
        self.provided_wts = False

        #os.system("cp land_template.vsz "+self.SROM_obj.input_folder)

        self.read_search_range()
        self.setup_bounds()

        os.chdir(self.SROM_obj.input_folder)

        #self.SROM_obj.input_folder = os.getcwd()

    def set_SROM_wts(self,_wt_dicts,_wt_locs):
        self.provided_wts = True
        self.wt_dict = _wt_dicts
        self.wt_locs = _wt_locs

    def read_search_range(self):
        m_range = self.SROM_obj.input_dict['Mean range']
        s_range = self.SROM_obj.input_dict['Std range']
        resolve = self.SROM_obj.input_dict['Landscape resolution'].rstrip('\r\n').split(',')
        self.resolution = [int(float(resolve[0])),int(float(resolve[1]))]

        m_range = m_range.split(',')
        self.mean_range = [float(m_range[0]),float(m_range[1])]

        s_range = s_range.split(',')
        self.std_range = [float(s_range[0]),float(s_range[1])]

    def setup_bounds(self):
        self.all_bounds = []

        self.all_bounds.append((self.mean_range[0],self.mean_range[1]))
        self.all_bounds.append((self.std_range[0],self.std_range[1]))

        self.x0 = np.array([0.5*(self.mean_range[0]+self.mean_range[1]),0.5*(self.std_range[0]+self.std_range[1])])

        self.MI_cutoff = math.log(len(self.SROM_obj.sample_locations))/math.log(2.0)
        
        #self.all_bounds.append((0.5,2.0))
        #self.all_bounds.append((0.1,1.5))

        #for k in range(0,len(self.SROM_obj.sample_locations)):
        #    self.all_bounds.append((0.0,1.0))

    def compute_MI(self,x0):
        if self.provided_wts==False:
            self.SROM_obj.set_dist(x0[0],x0[1])
            self.SROM_obj.create_beta_distribution()
            self.SROM_obj.initialize()
            self.SROM_obj.computeSROM()
        else:
            self.SROM_obj.sample_locations = copy.deepcopy(self.wt_locs)
            self.SROM_obj.p_SROM = copy.deepcopy(self.wt_dict[x0])

        self.MI_instance.get_input_pdf(self.SROM_obj.sample_locations,self.SROM_obj.p_SROM)
        self.MI_instance.compute_conditional_entropy()
        #self.MI_instance.read_input_wts(self.SROM_obj.p_SROM)
        #self.MI_instance.compute_conditional_entropy()

        MI = self.MI_instance.compute_mutual_information()

        if self.flag=='l':
            return_value = MI
        else:
            return_value = 0.5*(self.MI_cutoff-MI)**2
            print('MI track: ', MI,x0,return_value)	#self.SROM_obj.p_SROM

        return return_value

    def consProb(self,p_set):
        consVal = -1.0 + sum(p_set)

        return consVal

    def compute_CC(self):
        x0 = np.array([1.0,1.0])

        #cons = [{'type':'eq','fun':self.consProb}]
        #cons = [{'type':'ineq','fun': lambda x: x[1] - x[0]*2}]
        #result = minimize(self.compute_MI,self.x0,(),'SLSQP',None,None,None,self.all_bounds,options={'ftol':1e-10,'eps':0.0002,'maxiter':1000})
        result = minimize(self.compute_MI,self.x0,(),'TNC',None,None,None,self.all_bounds,None,options={'ftol':1e-6,'eps':0.05})

        #result = optimize.differential_evolution(self.compute_MI,self.all_bounds,tol=0.01)

        self.sol = result.x

        self.flag = 'l'
        self.CC = self.compute_MI(self.sol)

        print(self.sol)

        self.final_mu = result.x[0]
        self.final_std = result.x[1]

        f = open(self.casename+'_C.csv','w')
        print('C,Mean,Std',file=f)
        print(str(self.CC)+','+str(self.final_mu)+','+str(self.final_std),file=f)
        f.close()

        self.write_optimal_input_distribution(self.final_mu,self.final_std)

    def display_result(self):
        #print self.final_mu, self.final_std
        self.compute_MI(self.sol)
        print('Final conditional entropy: ', self.MI_instance.conditional_entropy)
        print('Final mean entropy: ', self.MI_instance.mean_entropy)
        print('MI: ', self.MI_instance.compute_mutual_information(), sum(self.MI_instance.mean_response_pdf.values()))
        print('Result: ', self.sol)

        mean = 0
        std = 0

        for i in range(0,len(self.SROM_obj.sample_locations)):
            mean += self.SROM_obj.sample_locations[i]*self.SROM_obj.p_SROM[i]
            std += (self.SROM_obj.sample_locations[i]**2)*self.SROM_obj.p_SROM[i]

        std += -(mean**2)

        std = math.sqrt(std)

        print(mean,std)

    def write_results(self):
        os.chdir(self.SROM_obj.input_folder)
        filename = self.casename+'_maxpdf.csv'
        outputfile = open(filename,'w')
        outputfile.close()

        print('Locations,Wts,cdf',file=open(filename,'a'))

        cdf = 0.0
        low_value = 10**float(self.SROM_obj.input_dict['Lower limit'])
        outstring = str(low_value)+','+str(0.0)+','+str(cdf)
        print(outstring, file=open(filename,'a'))

        for i in range(0,self.SROM_obj.sample_size):
            print(str(self.SROM_obj.sample_locations[i])+','+str(self.SROM_obj.p_SROM[i])+','+str(cdf), file=open(filename,'a'))
            cdf += self.SROM_obj.p_SROM[i]
            print(str(self.SROM_obj.sample_locations[i])+','+str(self.SROM_obj.p_SROM[i])+','+str(cdf), file=open(filename,'a'))

        up_value = 10**float(self.SROM_obj.input_dict['Upper limit'])
        print(str(up_value)+','+str(0.0)+','+str(cdf), file=open(filename,'a'))

        self.SROM_obj.write_true_CDF(self.casename)

    def plot_landscape(self,mean_range,std_range,MI_matrix):
        X = np.zeros(shape=(self.resolution[0],self.resolution[1]))
        Y = np.zeros(shape=(self.resolution[0],self.resolution[1]))

        m_range = list(10**np.array(mean_range))

        if hasattr(self,'CC'):
            c = "%.2f"%self.CC
            m = "%.2f"%self.final_mu
            st = "%.2f"%self.final_std
            CC_text = c+' ('+st+','+m+')'

        MImax = np.max(MI_matrix)
        MImin = np.min(MI_matrix)

        levels = list(np.linspace(MImin,MImax,6))
        print(levels[2:5])

        for i in range(0,int(self.resolution[0])):
            for j in range(0,int(self.resolution[1])):
                X[i,j] = std_range[j]
                Y[i,j] = m_range[i]

        fig = plt.figure()

        hot_reversed = plt.cm.get_cmap('hot_r')

        C = plt.pcolor(std_range,m_range,MI_matrix,cmap=hot_reversed)
        plt.rc('xtick',labelsize=14)
        plt.rc('ytick',labelsize=14)
        plt.xlabel(r"$\sigma(\log_{10} \mathrm{I_{ex}})$",fontsize=16)
        plt.ylabel(r"$10^{\langle \log_{10} \mathrm{I_{ex}} \rangle}$",fontsize=16)
        cs = plt.contour(X,Y,MI_matrix,np.array(levels[2:5]),cmap='Blues_r')
        plt.clabel(cs,fontsize='large')

        if hasattr(self,'CC'):
            plt.plot(self.final_std,10**self.final_mu,color='white',marker='.',markersize=10)
            plt.text(1.1*self.final_std,10**(0.9*self.final_mu),CC_text,fontsize=14,color='white')

        plt.yscale('log')
        fig.colorbar(C)

        fig.tight_layout()

        #plt.show()

        fig.savefig(self.casename+'.pdf', bbox_inches='tight')

    def write_optimal_input_distribution(self,m,s):
        ub = self.SROM_obj.upper_lim
        lb = self.SROM_obj.lower_lim

        s_m = (m-lb)/(ub-lb)
        s_v = (s**2)/((ub-lb)**2)

        p = s_m*(s_m*(1-s_m)/s_v - 1)
        q = (s_m*(1-s_m)/s_v - 1) - p

        xx = np.linspace(0.0,1.0,1000)
        sx = xx*(ub-lb) + lb
        true_x = 10**sx
        xcdf = beta.cdf(xx,p,q,0,1)
        #xpdf = beta.pdf(xx,p,q,0,1)
        xpdf = beta.pdf(sx,p,q,lb,(ub-lb))

        opti_cdf = np.zeros(shape=(1000,3))
        opti_cdf[:,0] = true_x
        opti_cdf[:,1] = xcdf
        opti_cdf[:,2] = xpdf

        np.savetxt('cc_pdf.csv',opti_cdf,delimiter=',')

    def get_MI_landscape(self):
        if self.provided_wts==False:
            self.wt_dict = {}

        mean_range = list(np.linspace(self.mean_range[0],self.mean_range[1],self.resolution[0]))
        std_range = list(np.linspace(self.std_range[0],self.std_range[1],self.resolution[1]))

        MI_matrix = np.zeros(shape=(len(mean_range),len(std_range)))

        for i in range(0,len(mean_range)):
            #this_line = str(10**mean_range[i])

            for j in range(0,len(std_range)):
                input_set = np.array([mean_range[i],std_range[j]])

                if self.provided_wts==False:
                    mi_value = self.compute_MI(input_set)

                    self.wt_dict[str(i)+','+str(j)] = copy.deepcopy(self.SROM_obj.p_SROM)
                else:
                    mi_value = self.compute_MI(str(i)+','+str(j))

                #print self.SROM_obj.p_SROM

                MI_matrix[i,j] = mi_value

            #print(i,j, 'completed')
            #sys.stdout.flush()

        max_MI = np.max(MI_matrix)
        coord = np.where(MI_matrix == max_MI)
        m_max, s_max = mean_range[coord[0][0]], std_range[coord[1][0]]

        os.chdir(self.SROM_obj.input_folder)

        self.write_optimal_input_distribution(m_max,s_max)

        of = open('C-coord.txt','w')
        print('C = ',max_MI,file=of)
        print('Mean = '+str(10**m_max),file=of)
        print('Std = '+str(10**s_max),file=of)
        of.close()

        filename = 'MI_landscape-'+self.casename+'.csv'
        wfile = open(filename,'w')
        wfile.close()

        m_values = ''

        for ms in mean_range:
            if self.SROM_obj.sample_space=='log':
                m_values += ','+str(10**ms)
            else:
                m_values += ','+str(ms)

        print(m_values,file=open(filename,'a'))

        for j in range(0,len(std_range)):
            if self.SROM_obj.sample_space=='log':
                this_line = str(10**std_range[j])
            else:
                this_line = str(std_range[j])

            for i in range(0,len(mean_range)):
                this_line += ','+str(MI_matrix[i,j])

            print(this_line,file=open(filename,'a'))

        #self.plot_landscape(mean_range,std_range,MI_matrix)
        #self.write_veusz_file()

        # cf = open('C-coord.txt')
        # print('C = '+str(max(max(MI_matrix))),file=cf)
        # cf.close()

    def write_veusz_file(self):
        with open(self.casename+'_landscape.vsz','w') as f:
            for line in all_lines:
                if 'AddImportPath' in line:
                    this_line = line.replace('/testpath',os.getcwd())
                    print(this_line,file=f)
                elif 'ImportFile2D' in line:
                    split_set = line.split()
                    this_line = line.replace('test_landscape.csv',self.casename+'_MI_landscape.csv')
                    print(this_line,file=f)
                else:
                    print(line,file=f)

    def compute_relative_entropy(self):
        mean_range = list(np.linspace(self.mean_range[0],self.mean_range[1],self.resolution[0]))
        std_range = list(np.linspace(self.std_range[0],self.std_range[1],self.resolution[0]))

        self.rel_ent = np.zeros(shape=(len(mean_range),len(std_range)))

        C_beta = BetaDist(self.SROM_obj.lower_lim,self.SROM_obj.upper_lim,self.sol[0],self.sol[1])
        print('p,q',C_beta.p_beta,C_beta.q_beta)

        self.B_C = sp.gamma(C_beta.p_beta)*sp.gamma(C_beta.q_beta)/sp.gamma(C_beta.p_beta + C_beta.q_beta)

        for m in range(0,len(mean_range)):
            for s in range(0,len(std_range)):
                beta_obj = BetaDist(self.SROM_obj.lower_lim,self.SROM_obj.upper_lim,mean_range[m],std_range[s])
                this_B = sp.gamma(beta_obj.p_beta)*sp.gamma(beta_obj.q_beta)/sp.gamma(beta_obj.p_beta + beta_obj.q_beta)
                this_di_p = sp.digamma(beta_obj.p_beta)
                this_di_q = sp.digamma(beta_obj.q_beta)
                this_di_pq = sp.digamma(beta_obj.p_beta + beta_obj.q_beta)

                self.rel_ent[m,s] = math.log(self.B_C/this_B) + (beta_obj.p_beta - C_beta.p_beta)*this_di_p + (beta_obj.q_beta - C_beta.q_beta)*this_di_q
                self.rel_ent[m,s] += (C_beta.p_beta - beta_obj.p_beta + C_beta.q_beta - beta_obj.q_beta)*this_di_pq

        os.chdir(self.SROM_obj.input_folder)

        filename = self.casename+'_relent.csv'
        wfile = open(filename,'w')
        wfile.close()

        std_values = ''

        for std in std_range:
            std_values += ','+str(std)

        print(std_values,file=open(filename,'a'))

        for i in range(0,len(mean_range)):
            this_line = str(10**mean_range[i])

            for j in range(0,len(std_range)):
                this_line += ','+str(self.rel_ent[i,j])

            print(this_line,file=open(filename,'a'))

    def write_land_and_gradient(self,MI_matrix,H_matrix,means,stds):
        MI_file = self.casename+'_MI_landscape.csv'
        wfile = open(MI_file,'w')
        wfile.close()

        H_file = self.casename+'_H.csv'
        wfile = open(H_file,'w')
        wfile.close()

        m_values = ''

        for ms in means:
            m_values += ','+str(10**ms)

        print(m_values,file=open(MI_file,'a'))
        print(m_values,file=open(H_file,'a'))

        for j in range(0,len(stds)):
            MI_line = str(stds[j])
            H_line = str(stds[j])

            for i in range(0,len(means)):

                MI_line += ','+str(MI_matrix[i,j])
                H_line += ','+str(H_matrix[i,j])

            print(MI_line,file=open(MI_file,'a'))
            print(H_line,file=open(H_file,'a'))

    def create_radial_law(self,MI_matrix,H_matrix):
        max_H = math.log10(np.max(H_matrix)*1.05)
        min_H = -2#math.log10(np.min(H_matrix)*0.95)
        bin_size = 10#int(self.resolution/2)
        H_bins = np.linspace(min_H,max_H,bin_size+1)
        d_H = H_bins[1]-H_bins[0]
        MI_samples = {}

        for k in range(0,bin_size):
            MI_samples[k] = []

        for i in range(0,int(self.resolution[0])):
            for j in range(0,int(self.resolution[1])):
                H = H_matrix[i,j]
                MI = MI_matrix[i,j]

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

    def landscape_and_gradient(self):
        self.compute_CC()
        self.write_results()
        #self.flag = 'l'

        mean_range = list(np.linspace(self.mean_range[0],self.mean_range[1],self.resolution[0]))
        std_range = list(np.linspace(self.std_range[0],self.std_range[1],self.resolution[1]))
        MI_matrix = np.zeros(shape=(len(mean_range),len(std_range)))
        H_matrix = np.zeros(shape=(len(mean_range),len(std_range)))

        C_beta = BetaDist(self.SROM_obj.lower_lim,self.SROM_obj.upper_lim,self.sol[0],self.sol[1])
        self.B_C = sp.gamma(C_beta.p_beta)*sp.gamma(C_beta.q_beta)/sp.gamma(C_beta.p_beta + C_beta.q_beta)

        for i in range(0,len(mean_range)):
            this_line = str(10**mean_range[i])

            for j in range(0,len(std_range)):
                input_set = np.array([mean_range[i],std_range[j]])
                mi_value = self.compute_MI(input_set)
                #print self.SROM_obj.p_SROM

                MI_matrix[i,j] = mi_value

                beta_obj = BetaDist(self.SROM_obj.lower_lim,self.SROM_obj.upper_lim,mean_range[i],std_range[j])
                this_B = sp.gamma(beta_obj.p_beta)*sp.gamma(beta_obj.q_beta)/sp.gamma(beta_obj.p_beta + beta_obj.q_beta)
                this_di_p = sp.digamma(beta_obj.p_beta)
                this_di_q = sp.digamma(beta_obj.q_beta)
                this_di_pq = sp.digamma(beta_obj.p_beta + beta_obj.q_beta)

                try:
                    H_matrix[i,j] = math.log(self.B_C) - math.log(this_B) + (beta_obj.p_beta - C_beta.p_beta)*this_di_p + (beta_obj.q_beta - C_beta.q_beta)*this_di_q
                    H_matrix[i,j] += (C_beta.p_beta - beta_obj.p_beta + C_beta.q_beta - beta_obj.q_beta)*this_di_pq
                    H_matrix[i,j] *= 1.0/math.log(2)
                except ValueError:
                    #print(self.B_C,this_B)
                    #(mean_range[i],std_range[j])
                    H_matrix[i,j] = 100.0

        self.plot_landscape(mean_range,std_range,MI_matrix)
        self.write_land_and_gradient(MI_matrix,H_matrix,mean_range,std_range)
        self.create_radial_law(MI_matrix,H_matrix)
