"""Class for computing error functinals of a discrete random variable against beta distributions.
"""

from scipy.stats import beta
import string
import numpy as np
import sys
import math
import numpy.random as rand

class BetaDist(object):
    def __init__(self,_ll,_ul,_mu,_std):
        # Beta distribution limits, mean, and standard deviation
        self.lower_lim = _ll
        self.upper_lim = _ul
        self.mean = _mu
        self.variance = _std**2

        self.scale_parameters()
        self.compute_shape_params()

        self.approx_sigma = 0.01

    def scale_parameters(self):
        self.scaled_mean = (self.mean - self.lower_lim)/(self.upper_lim - self.lower_lim)
        self.scaled_variance = self.variance/((self.upper_lim - self.lower_lim)**2)

    def compute_shape_params(self):
        self.p_beta = self.scaled_mean*(self.scaled_mean*(1.0-self.scaled_mean)/self.scaled_variance - 1.0)
        self.q_beta = (self.scaled_mean*(1.0-self.scaled_mean)/self.scaled_variance - 1.0) - self.p_beta
        #print self.p_beta, self.q_beta, self.scaled_mean, self.scaled_variance
        #sys.stdout.flush()

    def create_moments(self,_n_moments):
        self.n_moments = _n_moments

        self.true_moments = []

        for i in range(0,self.n_moments):
            self.true_moments.append(beta.moment(i+1,self.p_beta,self.q_beta))

        #print self.true_moments

    def get_sample_locations(self,_sample_locs,_sample_space):
        self.sample_locations = []
        self.sample_space = _sample_space

        for i in range(0,len(_sample_locs)):
            if self.sample_space=='log':
                this_sample = math.log10(_sample_locs[i])
            else:
                this_sample = _sample_locs[i]

            self.sample_locations.append((this_sample-self.lower_lim)/(self.upper_lim-self.lower_lim))

    def create_CDF_range(self):
        size = 250

        self.cdf_x = list(np.linspace(0.0,1.0,size+1))
        self.cdf_wts = []

        for sample in self.cdf_x:
            self.cdf_wts.append(beta.cdf(sample,self.p_beta,self.q_beta))

    def write_true_CDF(self,casename):
        filename = casename+'true_cdf.csv'
        ofile = open(filename,'w')
        ofile.close()

        print('Truepts,Truecdfs',file=open('true_cdf.csv','a'))

        print(len(self.cdf_x))

        for sample in self.cdf_x:
            self.cdf_wts.append(beta.cdf(sample,self.p_beta,self.q_beta))

            if self.sample_space=='log':
                scaled_back = 10**(sample*(self.upper_lim - self.lower_lim) + self.lower_lim)
            else:
                scaled_back = (sample*(self.upper_lim - self.lower_lim) + self.lower_lim)

            print(str(scaled_back)+','+str(self.cdf_wts[-1]),file=open(filename,'a'))

    def compute_CDF_error(self,_trial_wts):
        SROM_error = 0.0
        SROM_CDF = 0.0
        SROM_idx = 0

        location = self.sample_locations[0]

        for i in range(0,len(self.cdf_x)-1):
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

        for i in range(1,self.n_moments+1):
            loc_moment = 0.0
            for j in range(0,len(self.sample_locations)):
                loc_moment += _trial_wts[j]*self.sample_locations[j]**i

            self.moments.append(loc_moment)

        error = 0

        for i in range(0,len(self.true_moments)):
            error += ((self.true_moments[i]-self.moments[i])/self.true_moments[i])**2

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

    def prob_derivative(self,_trial_wts):
        """This function computes the derivative of the SROM functional
        with respect to the probability weights
        """
        eachProb_CDF = np.zeros(shape=(len(_trial_wts)))
        eachProb_moments = np.zeros(shape=(len(_trial_wts)))

        location = self.sample_locations[0]

        SROM_CDF = 0.0
        SROM_idx = 0

        # Prob derivative error from CDF error function
        #for x_key in self.CDF_dict.keys():
        for i in range(0,len(self.cdf_x)-1):
            # if self.cdf_x[i]>location:
            #     SROM_CDF += _trial_wts[SROM_idx]
            #     SROM_idx += 1

            SROM_CDF = 0.0

            for sampNo in range(0,len(_trial_wts)):
                SROM_CDF += 0.5*_trial_wts[sampNo]*(1+math.erf((self.cdf_x[i]-self.sample_locations[sampNo])/(math.sqrt(2)*self.approx_sigma)))

            true_cdf = 0.5*(self.cdf_wts[i]+self.cdf_wts[i+1])

            commonIntegrand = (SROM_CDF-true_cdf)*(self.cdf_x[i+1]-self.cdf_x[i])

            for probNo in range(0,len(_trial_wts)):
                eachProb_CDF[probNo] += 0.5*commonIntegrand*(1+math.erf((self.cdf_x[i]-self.sample_locations[probNo])/(math.sqrt(2)*self.approx_sigma)))

        # Prob derivative error from moment error function
        SROMmnts = np.zeros(shape=(self.n_moments,1))
        for sampNo in range(0,len(_trial_wts)):
            for mo in range(0,self.n_moments):
                SROMmnts[mo] += (self.sample_locations[sampNo]**(mo+1))*_trial_wts[sampNo]

        for probNo in range(0,len(_trial_wts)):
            for mo in range(0,self.n_moments):
                mFac = ((SROMmnts[mo]-self.true_moments[mo])/self.true_moments[mo]**2)
                eachProb_moments[probNo] += mFac*(self.sample_locations[probNo]**(mo+1))

        return eachProb_CDF, eachProb_moments

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
