#!/usr/bin/env python

from scipy.stats import lognorm
from scipy.stats import norm
import string
import numpy as np
import sys
import math
import numpy.random as rand

class LogNormalDist(object):
    def __init__(self,_log_mu,_log_std):
        self.log_mu = _log_mu
        self.log_std = _log_std
        
        self.mu = 10**(self.log_mu + 0.5*self.log_std**2)
        self.var = (10**(self.log_std**2)-1)*(10**(2*self.log_mu+self.log_std**2))
        
        self.s = self.log_std
        self.scale = 10**self.log_mu
        
#    def compute_moments(self,_n_moments):
#        self.n_moments = _n_moments
#        self.moments = []
#
#	for i in xrange(1,_n_moments+1):
#	    self.moments.append(10.0**(i*self.log_mu + 0.5*(i*self.log_std)**2))

    def create_CDF_range(self,_sample_locs):
        self.sample_locations = []
        self.sample_loc = []
        
        for sample in _sample_locs:
            self.sample_locations.append(sample)
            self.sample_loc.append(sample/500.0)
        
        self.cdf_x = []
        self.cdf_wts = []
        
        ll = 0.1*min(self.sample_locations)/500.0
        ul = 10.0*max(self.sample_locations)/500.0
        
        size = 250
        
        self.log_x = list(np.linspace(math.log10(ll),math.log10(ul),size+1))
        self.cdf_wts = []
        
        ofile = open('true_cdf.csv','w')
        
        print >> ofile, 'Truepts,Truewts'
        
        for sample in self.log_x:
            self.cdf_x.append(10**sample)
            #self.cdf_wts.append(lognorm.cdf(self.cdf_x[-1],self.s,self.log_mu,self.scale))
            self.cdf_wts.append(norm.cdf(sample,self.log_mu,scale=self.log_std))
            
            print >> ofile, str(self.cdf_x[-1])+','+str(self.cdf_wts[-1])
            
        ofile.close()
        
    def create_moments(self,n_moments):
        self.true_moments = []
        
        for k in xrange(0,n_moments):
            exponent = (k+1)*self.log_mu + 0.5*((k+1)*self.log_std)**2
            self.true_moments.append(10**exponent)
            
        #print self.true_moments
        
    def compute_CDF_error(self,_trial_wts):
        SROM_error = 0.0
        SROM_CDF = 0.0
        SROM_idx = 0
        
        location = self.sample_loc[0]
        
        for i in xrange(0,len(self.cdf_x)-1):
            if self.cdf_x[i]>location:
                SROM_CDF += _trial_wts[SROM_idx]
                SROM_idx += 1
                
                if SROM_idx<len(_trial_wts):
                    location = self.sample_loc[SROM_idx]
                else:
                    location = self.cdf_x[-1]+0.1
                
            SROM_error += (self.cdf_x[i+1]-self.cdf_x[i])*(SROM_CDF-self.cdf_wts[i])**2
            
        return (0.5*SROM_error)
    
    def compute_moment_error(self,_trial_wts):
        self.moments = []
        
        for i in xrange(1,len(self.true_moments)+1):
            loc_moment = 0.0
            for j in xrange(0,len(self.sample_loc)):
                loc_moment += _trial_wts[j]*self.sample_loc[j]**i
                
            self.moments.append(loc_moment)
            
        error = 0
        
        for i in xrange(0,len(self.true_moments)):
            error += ((math.log10(self.true_moments[i])-math.log10(self.moments[i])))**2#/self.true_moments[i])**2
            
        return (0.5*error)
        
    def scale_parameters(self):
        self.scaled_mean = (self.mean - self.lower_lim)/(self.upper_lim - self.lower_lim)
        self.scaled_variance = self.variance/((self.upper_lim - self.lower_lim)**2)
    
    def get_variable_bounds(self):
        return [self.lower_lim, self.upper_lim]
            
    def initialize_dist_structs(self,total_intervals,total_moments):
        self.get_variable_bounds()
        self.scale_parameters()
        self.compute_shape_params()
        self.create_CDF_dict(total_intervals)
        self.get_moments(total_moments)
    
    def sample_derivative(self,x_set,p_set):
        """This function computes the derivative of the SROM error functional
        with respect to the random variables
        """
        eachSampError = np.zeros(shape=(len(x_set)))
        pFactor = np.zeros(shape=(len(x_set)))
        
        for sampNo in range(0,len(x_set)):
            pFactor[sampNo] = -p_set[sampNo]/np.sqrt(2*math.pi*self.delta_sigma**2)
            
        # Sample derivative error from CDF error function
        for x_key in self.CDF_dict.keys():
            x_index = int(x_key/self.interval_size)
            true_cdf = self.CDF_dict[x_key]
            
            SROM_cdf = 0.0
            for srom_idx in range(0,len(x_set)):
                if x_set[srom_idx]<=x_key:
                    SROM_cdf += p_set[srom_idx]
                    
            commonIntegrand = (true_cdf-SROM_cdf)*self.interval_size
            
            for sampNo in range(0,len(x_set)):
                eachSampError[sampNo] += pFactor[sampNo]*commonIntegrand*math.exp(-self.inverse_variance*(x_key-x_set[sampNo])**2)

        # Sample derivative error from moment error function
        SROMmnts = np.zeros(shape=(self.total_moments,1))
        for sampNo in range(0,len(x_set)):
            for mo in range(0,len(self.scaled_moments)):
                SROMmnts[mo] += (x_set[sampNo]**(mo+1))*p_set[sampNo]
        
        for sampNo in range(0,len(x_set)):
            for mo in range(0,len(self.scaled_moments)):
                mFac = ((SROMmnts[mo]-self.scaled_moments[mo])/self.scaled_moments[mo]**2)
                eachSampError[sampNo] += mFac*(mo+1)*p_set[mo]*(x_set[sampNo]**mo)
            
        return eachSampError
    
    def prob_derivative(self,x_set,p_set):
        """This function computes the derivative of the SROM functional
        with respect to the probability weights
        """
        eachProbError = np.zeros(shape=(len(p_set)))

        # Prob derivative error from CDF error function
        for x_key in self.CDF_dict.keys():
            x_index = int(x_key/self.interval_size)
            true_cdf = self.CDF_dict[x_key]
            
            SROM_cdf = 0.0
            for srom_idx in range(0,len(x_set)):
                if x_set[srom_idx]<=x_key:
                    SROM_cdf += p_set[srom_idx]
                    
            commonIntegrand = (true_cdf-SROM_cdf)*self.interval_size
            
            for probNo in range(0,len(p_set)):
                eachProbError[probNo] += 0.5*commonIntegrand*(1+math.erf((x_key-x_set[probNo])*np.sqrt(self.inverse_variance)))
            
        # Prob derivative error from moment error function
        SROMmnts = np.zeros(shape=(self.total_moments,1))
        for sampNo in range(0,len(x_set)):
            for mo in range(0,self.total_moments):
                SROMmnts[mo] += (x_set[sampNo]**(mo+1))*p_set[sampNo]
        
        for probNo in range(0,len(p_set)):
            for mo in range(0,self.total_moments):
                mFac = ((SROMmnts[mo]-self.scaled_moments[mo])/self.scaled_moments[mo]**2)
                eachProbError[probNo] += mFac*(x_set[probNo]**(mo+1))
            
        return eachProbError
    
    def scale_x_values(self,x_set):
        scaled_x = np.zeros(shape=(len(x_set)))
        
        for sample_no in range(0,len(x_set)):
            scaled_x[sample_no] = (x_set[sample_no] - self.lower_lim)/(self.upper_lim - self.lower_lim)
            
        return scaled_x
    
    def rescale_x_values(self,x_set):
        rescaled_x = np.zeros(shape=(len(x_set)))
        
        for sample_no in range(0,len(x_set)):
            rescaled_x[sample_no] = x_set[sample_no]*(self.upper_lim - self.lower_lim) + self.lower_lim
            
        return rescaled_x