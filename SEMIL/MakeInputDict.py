# Code to obtain shear modulus distribution from coarse graininig
import numpy as np
import os

class MakeInputDict(object):
    def __init__(self):
        pass
    
    def make_dict(self,_inputs):
        self.input_dict = {}
        
        for this_input in _inputs:
            if "=" in this_input:
                this_set = this_input.split("=")
                this_set[0] = this_set[0].rstrip(' ')
                this_set[1] = this_set[1].rstrip('\r\n').lstrip(' ')
                
                self.input_dict[this_set[0]] = this_set[1]