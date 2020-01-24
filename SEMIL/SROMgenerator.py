"""Class for computing Stochastic Reduced Order Model for log normal distribution for a fixed set of points.
"""

import os, sys
import numpy as np
import math
#import matplotlib.pyplot as plt
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import fmin_slsqp
from scipy.optimize import root
from MakeInputDict import MakeInputDict
import string
#from LogNormalDist import LogNormalDist
from scipy.optimize import minimize
from BetaDist import BetaDist

class SROMgenerator(MakeInputDict):
    def __init__(self,data_path,input_lines):
        """ This function reads the lines from the input file and creates the dictionary of input quantities.
        Args:
            input_lines (list): The list of input data
        """
        self.make_dict(input_lines)

        self.no_rand_variables = 1

        self.input_folder = data_path+'/'+self.input_dict['Input folder']
        self.sample_file = self.input_dict['Sample file']                       # No. of samples for SROM
        self.total_moments = int(self.input_dict['Moment order'])                    # No. of moments
        #self.int_size = int(input_dict['Integration intervals'])              # intervals for numerical integration
        #self.optim_method = self.input_dict['Optimization method']                   # Method for optimization

        # Parameters for log normal distribution
        #self.log_mu = float(self.input_dict['Log mean'])
        #self.log_std = float(self.input_dict['Log standard deviation'])

        # Parameters for beta distribution
        self.upper_lim = float(self.input_dict['Upper limit'])
        self.lower_lim = float(self.input_dict['Lower limit'])

        # weights for error functionals
        self.cdf_wt = float(self.input_dict['CDF error weight'])
        self.mmt_wt = float(self.input_dict['Moment error weight'])

        self.dist_type = self.input_dict['Distribution type']
        self.sample_space = self.input_dict['Sample space']

        # List of random variables whose SROM we want to construct
        self.target_variables = []

        # Sigma for smooth approximation of SROM CDF
        self.smooth_sigma = 0.01

        self.read_sample_locations()

        # if self.dist_type=='beta':
        #     self.create_beta_distribution()
        # elif self.dist_type=='normal':
        #     self.create_normal_distribution()
        # else:
        #     sys.exit()

    def read_sample_locations(self):
        self.sample_locations = []

        old_dir = os.getcwd()

        os.chdir(self.input_folder)

        ifile = open(self.sample_file)
        all_lines = ifile.readlines()
        ifile.close()

        for line in all_lines:
            self.sample_locations.append(float(line.rstrip('\r\n')))

        self.sample_size = len(self.sample_locations)
        self.sample_wts = (1.0/float(self.sample_size))*np.ones(shape=(len(self.sample_locations)))
        self.p_SROM = self.sample_wts

        os.chdir(old_dir)

    def set_dist(self,_mu,_std):
    	self.log_mu = _mu
    	self.log_std = _std

    def create_beta_distribution(self):
        self.distribution = BetaDist(self.lower_lim,self.upper_lim,self.log_mu,self.log_std)
        #print self.distribution.p_beta, self.distribution.q_beta
        self.distribution.get_sample_locations(self.sample_locations,self.sample_space)
        self.distribution.create_CDF_range()
        self.distribution.create_moments(self.total_moments)

        # Set Lagrange multiplier
        self.lagMult = 1.0

        self.setupBounds()

    def create_normal_distribution(self):
        self.distribution = NormalDist(self.lower_lim,self.upper_lim,self.log_mu,self.log_std)
        #print self.distribution.p_beta, self.distribution.q_beta
        self.distribution.get_sample_locations(self.sample_locations,self.sample_space)
        self.distribution.create_CDF_range()
        self.distribution.create_moments(self.total_moments)

        # Set Lagrange multiplier
        self.lagMult = 1.0

        self.setupBounds()

    def create_distribution(self):
        self.distribution = LogNormalDist(self.log_mu,self.log_std)
        self.distribution.create_CDF_range(self.sample_locations)
        self.distribution.create_moments(self.total_moments)

        # Set Lagrange multiplier
        self.lagMult = 1.0

        self.setupBounds()

    def compute_true_moments(self):
    	self.true_moments = []

    	for i in xrange(1,self.total_moments+1):
    	    self.true_moments.append(math.exp(i*self.log_mu + 0.5*(i*self.log_std)**2))

    def initialize(self):
        # Get bounds for the SROM initialguess
        #self.initial_guesses()

        # Set Lagrange multiplier
        self.lagMult = 1.0

        self.setupBounds()
        #self.normalize_correlations()

    def setupBounds(self):
        pBounds = []

        for sample_no in range(0,self.sample_size):
            pBounds.append((0.0,1.0))

            #self.allBounds = xBounds + pBounds
        self.allBounds = pBounds

    def computeSROM(self):
        self.optimSLSQP()

    def optimSLSQP(self):
        cons = [{'type':'eq','fun':self.consProb}]

        #res = minimize(self.objfunSLSQP,totalArray,(),'SLSQP',None,None,None,self.allBounds,cons,options={'ftol':1e-8,'maxiter':2000,'eps':1e-5})
        #x, fn, its, imode, smode = fmin_slsqp(self.objfunSLSQP,self.p_SROM,[],cons,[],None,self.allBounds,None,None,None,(),10000,1e-12,0,None,1)

        result = minimize(self.objfunSLSQP,self.p_SROM,(),'SLSQP',None,None,None,self.allBounds,cons,options={'ftol':1e-16})
        #result = minimize(self.objfunSLSQP,self.p_SROM,(),'SLSQP',self.objfunSLSQPprime,None,None,self.allBounds,cons,options={'ftol':1e-16})

        self.extract_solutions(result.x)

    def objfunSLSQP(self,totalArray):
        total_error = 0.0
        # Extract the probability values
        self.p_SROM = totalArray[len(totalArray)-self.sample_size:len(totalArray)]

        error_1 = self.cdf_wt*self.distribution.compute_CDF_error(self.p_SROM)
        error_2 = self.mmt_wt*self.distribution.compute_moment_error(self.p_SROM)

        total_error = error_1 + error_2

        return total_error

    def consProb(self,p_set):
        consVal = -1.0 + sum(p_set)

        return consVal

    def objfunSLSQPprime(self,totalArray):
        total_error = 0.0
        #fprimeSamp = np.zeros(shape=(self.sample_size*self.no_rand_variables))
        fprimeProb = np.zeros(shape=(self.sample_size))

        # Extract the probability values
        self.p_SROM = totalArray[len(totalArray)-self.sample_size:len(totalArray)]

        # Extract the data values
        #for variable_no in range(0,self.no_rand_variables):
        #    low_idx = variable_no*self.sample_size
        #    up_idx = (variable_no+1)*self.sample_size
        #
        #    x_set = totalArray[low_idx:up_idx]
        #    for sample_no in range(0,self.sample_size):
        #        self.x_SROM[variable_no,sample_no] = x_set[sample_no]

        #total_prime = np.zeros(shape=((self.no_rand_variables+1)*self.sample_size))
        idx = 0

        for variable_no in range(0,self.no_rand_variables):
            #variable_derivative = self.target_variables[variable_no].sample_derivative(self.x_SROM[variable_no,:],self.p_SROM)
            prob_deriv_CDF, prob_deriv_moment = self.distribution.prob_derivative(self.p_SROM)

            for sample_no in range(0,self.sample_size):
                idx = variable_no*self.sample_size + sample_no
                #fprimeSamp[idx] = variable_derivative[sample_no]

                fprimeProb[sample_no] += self.cdf_wt*prob_deriv_CDF[sample_no] + self.mmt_wt*prob_deriv_moment[sample_no]

        #fprimeSamp = self.sampleDeriv(xSet,pSet)
        #fprimeProb = self.probDeriv(xSet,pSet)

        #totalPrime = np.hstack((fprimeSamp,fprimeProb))
    	#totalPrime = fprimeProb

        return fprimeProb

    def getCDFvalue(self,xInt,xSet,pSet):
        sum = 0.0
        for sampNo in range(0,len(xSet)):
            #sum += 0.5*pSet[sampNo]*(1.0 + math.erf((xInt-xSet[sampNo])/np.sqrt(self.var2)))
            if xInt >= xSet[sampNo]:
                sum += pSet[sampNo]
        return sum

    def write_output(self):
        outputfile = open('SROM.csv','w')
        outputfile.close()

        print('SROMpts,SROMwts', file=open('SROM.csv','a'))

        outstring = str((self.sample_locations[0]/10.0))+','+'0.0'

        print(outstring,file=open('SROM.csv','a'))

        cdf = 0.0

        for i in xrange(0,self.sample_size):
            outstring = str(self.sample_locations[i])+','+str(cdf)+','+str(self.p_SROM[i])
            print(outstring, file=open('SROM.csv','a'))

            cdf += self.p_SROM[i]

            outstring = str(self.sample_locations[i])+','+str(cdf)+','+str(self.p_SROM[i])
            print(outstring, file=open('SROM.csv','a'))

        outstring = str(self.sample_locations[-1]*10)+','+str(cdf)+','+str(self.p_SROM[i])

        print(outstring,file=open('SROM.csv','a'))

        ofile = open('moments.csv','w')
        ofile.close()

        for i in xrange(0,self.total_moments):
            outputstring = str(self.distribution.true_moments[i])+','+str(self.distribution.moments[i])
            print(outputstring,file=open('moments.csv','a'))

    def write_summary(self):
        print('Mean: ', self.distribution.true_moments[0])
        print('Median: ', 10**self.distribution.log_mu)
        #print('Standard deviation: ', math.sqrt(self.distribution.true_moments[1]-self.distribution.true_moments[0]**2)

    def extract_solutions(self,x):
        # Extract the probability values
        self.p_SROM = x

        # Extract the data values
        #for variable_no in range(0,self.no_rand_variables):
        #    low_idx = variable_no*self.sample_size
        #    up_idx = (variable_no+1)*self.sample_size
        #
        #    x_set = x[low_idx:up_idx]
        #    rescaled_x = self.target_variables[variable_no].rescale_x_values(x_set)
        #    self.x_SROM[variable_no,:] = rescaled_x

    def plot_result(self):
        data_keys = self.target_variables[0].CDF_dict.keys()
        data_size =len(data_keys)
        data_keys.sort()

        trueCDF = np.zeros(shape=(data_size))
        sromCDF = np.zeros(shape=(data_size))

        for xNo in range(0,data_size):
            x = data_keys[xNo]
            trueCDF[xNo] = self.target_variables[0].CDF_dict[x]
            sromCDF[xNo] = self.get_SROM_CDF(0,x)

        plt.plot(data_keys,trueCDF,c='red')
        plt.plot(data_keys,sromCDF,c='blue')
        plt.show()

    def get_SROM_CDF(self,variable_no,x):
        srom_cdf = 0.0

        for srom_idx in range(0,self.sample_size):
            if x >= self.x_SROM[variable_no,srom_idx]:
                srom_cdf += self.p_SROM[srom_idx]

        return srom_cdf

    def write_true_CDF(self,casename):
        self.distribution.write_true_CDF(casename)
