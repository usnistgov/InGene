from scipy.stats import norm
import string
import numpy as np
import sys
import math
import numpy.random as rand

class NormalDist(object):
    def __init__(self,_ll,_ul,_mu,_std):
        # Beta distribution limits, mean, and standard deviation
        self.lower_lim = _ll
        self.upper_lim = _ul
        self.mean = _mu
        self.variance = _std**2
        self.scale = _std
        
        self.scale_parameters()
        #self.compute_shape_params()

    def scale_parameters(self):
        self.scaled_mean = (self.mean - self.lower_lim)/(self.upper_lim - self.lower_lim)
        self.scaled_variance = self.variance/((self.upper_lim - self.lower_lim)**2)
        
    def create_moments(self,_n_moments):
        self.n_moments = _n_moments
        
        self.true_moments = []
        
        for i in xrange(0,self.n_moments):
            self.true_moments.append(0.0)
        
        for i in range(0,len(self.cdf_x)-1):         
            for j in xrange(0,self.n_moments):
                local_moments = ((self.cdf_x[i]+self.cdf_x[i+1])**(j+1))
                self.true_moments[j] += 0.5*local_moments*(self.pdf_wts[i]+self.pdf_wts[i+1])*(self.cdf_x[i+1]-self.cdf_x[i])
                #self.true_moments.append(norm.moment(i+1,loc=self.mean,scale=self.scale))
            
        #print self.true_moments
            
    def get_sample_locations(self,_sample_locs):
        self.sample_locations = []
        
        for i in xrange(0,len(_sample_locs)):
            this_sample = math.log10(_sample_locs[i]/500)
            self.sample_locations.append(this_sample)

    def create_CDF_range(self):
        """Create a truncated normal distribution within the range
        """
        size = 500
        
        self.cdf_x = list(np.linspace(self.lower_lim,self.upper_lim,size+1))
        self.cdf_wts = []
        self.pdf_wts = []
        
        for sample in self.cdf_x:
            self.pdf_wts.append(norm.pdf(sample,loc=self.mean,scale=self.scale))
            
        total_weight = 0.0
        
        for i in range(0,size-1):
            total_weight += 0.5*(self.pdf_wts[i]+self.pdf_wts[i+1])*(self.cdf_x[i+1]-self.cdf_x[i])
        
        ofile = open('true_cdf.csv','w')
        
        print >> ofile, 'Truepts,Truecdfs,Truewts'
        
        cdf_value = 0
        self.cdf_wts.append(0.0)
        
        for i in range(0,size):
            print >> ofile, str(self.cdf_x[i])+','+str(self.cdf_wts[-1])+','+str(self.pdf_wts[i]/total_weight)
            
            if i<size-1:
                cdf_value += 0.5*(self.pdf_wts[i]+self.pdf_wts[i+1])*(self.cdf_x[i+1]-self.cdf_x[i])
                self.cdf_wts.append(cdf_value/total_weight)
            
            #self.cdf_wts.append(norm.cdf(sample,loc=self.mean,scale=self.scale))
            #print >> ofile, str(sample)+','+str(self.cdf_wts[-1])+','+str(norm.pdf(sample,loc=self.mean,scale=self.scale))
            
        ofile.close()
        
    def compute_CDF_error(self,_trial_wts):
        SROM_error = 0.0
        SROM_CDF = 0.0
        SROM_idx = 0
        
        location = self.sample_locations[0]
        
        for i in xrange(0,len(self.cdf_x)-1):
            if self.cdf_x[i]>location:
                SROM_CDF += _trial_wts[SROM_idx]
                SROM_idx += 1
                
                if SROM_idx<len(_trial_wts):
                    location = self.sample_locations[SROM_idx]
                else:
                    location = self.cdf_x[-1]+0.1
                    
            this_cdf_wt = 0.5*(self.cdf_wts[i]+self.cdf_wts[i+1])
                
            SROM_error += (self.cdf_x[i+1]-self.cdf_x[i])*((SROM_CDF-this_cdf_wt)**2)
            
        return (0.5*SROM_error)
    
    def compute_moment_error(self,_trial_wts):
        self.moments = []
        
        for i in xrange(1,self.n_moments+1):
            loc_moment = 0.0
            for j in xrange(0,len(self.sample_locations)):
                loc_moment += _trial_wts[j]*self.sample_locations[j]**i
                
            self.moments.append(loc_moment)
            
        error = 0
        
        for i in xrange(0,len(self.true_moments)):
            if abs(self.true_moments[i])>0.0:
                error += ((self.true_moments[i]-self.moments[i])/self.true_moments[i])**2
            else:
                error += (self.true_moments[i]-self.moments[i])**2
            
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
    
    def prob_derivative(self,p_set):
        """This function computes the derivative of the SROM functional
        with respect to the probability weights
        """
        eachProbError = np.zeros(shape=(len(p_SROM)))

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