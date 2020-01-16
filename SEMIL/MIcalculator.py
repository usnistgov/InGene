"""Class for computing mutual information for an alphabet and discrete input distribution
"""
import numpy as np
import random as rand
import os
import math
from MakeInputDict import MakeInputDict

class MIcalculator(MakeInputDict):
    def __init__(self,response_file,_input_size):
        #self.read_input_pdf(input_file)
        self.input_size = _input_size
        self.read_responses(response_file)
        #self.compute_mean_response()

    def read_input_pdf(self,input_file):
        ifile = open(input_file,'r')
        all_lines = ifile.readlines()
        ifile.close()

        self.input_loc = []
        self.input_wts = []

        for line in all_lines:
            this_set = line.rstrip('\r\n').split(',')
            self.input_loc.append(float(this_set[0]))
            self.input_wts.append(float(this_set[1]))

        self.input_size = len(self.input_loc)

        print(sum(self.input_wts))

    def get_input_pdf(self,_input_loc,_input_wts):
        self.input_loc = _input_loc
        self.input_wts = _input_wts

    def read_responses(self,response_file):
        ifile = open(response_file,'r')
        all_lines = ifile.readlines()
        ifile.close()

        self.conditional_entropy = 0.0

        #print all_lines

        self.response_pdfs = {}
        self.response_wts = {}
        self.partial_entropy = []

        #response_size = len(self.input_loc)
        response_size = self.input_size

        for i in range(0,response_size):
            this_entropy = 0.0
            loc_set = all_lines[2*i].rstrip('\r\n').split(',')
            wt_set = all_lines[2*i+1].rstrip('\r\n').split(',')

            self.response_pdfs[i] = {}

            set_size = len(loc_set)

            for k in range(0,set_size):
                if loc_set[k]!='':
                    #self.response_pdfs[i][int(loc_set[k])] = float(wt_set[k])
                    self.response_pdfs[i][loc_set[k]] = float(wt_set[k])
                    try:
                        this_entropy += float(wt_set[k])*math.log(float(wt_set[k]))
                    except ValueError:
                        pass

            self.partial_entropy.append(this_entropy)

            #self.conditional_entropy -= self.input_wts[i]*this_entropy
            #print sum(self.response_pdfs[i].values())

    def compute_conditional_entropy(self):
        self.conditional_entropy = 0.0

        for i in range(0,self.input_size):
            self.conditional_entropy -= self.input_wts[i]*self.partial_entropy[i]

        #print self.conditional_entropy

    def compute_mutual_information(self):
        key_set = []
        self.mean_response_pdf = {}
        self.mean_entropy = 0.0

        for k in range(0,self.input_size):
            for key in self.response_pdfs[k].keys():
                if key not in key_set:
                    key_set.append(key)
                    self.mean_response_pdf[key] = 0.0

        for key in key_set:
            for response_no in range(0,self.input_size):
                try:
                    self.mean_response_pdf[key] += self.input_wts[response_no]*self.response_pdfs[response_no][key]
                except KeyError:
                    pass

        for value in self.mean_response_pdf.values():
            try:
                self.mean_entropy -= value*math.log(value)
            except ValueError:
                pass

        #print self.mean_entropy
        #print 'Mutual information: ', self.mean_entropy - self.conditional_entropy

        self.mutual_information = (self.mean_entropy - self.conditional_entropy)/math.log(2.0)

        return max(0.0,self.mutual_information)
