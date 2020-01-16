# Run gel modeling code
import sys
from optparse import OptionParser
from MIcalculator import MIcalculator

def processOptions():
    parser = OptionParser()
    parser.add_option("-r", dest="response", help="Name of the file containing the inputs.", default="")
    #parser.add_option("-f", dest="inputFormat", help="Use 'show' to see the input format. Or specify the name to create an empty inputfile with the placeholders.", default="")
    parser.add_option("-i", dest="input", help="Number of sample trajectories to be created.", default="")
    
    [options, args] = parser.parse_args()
    
    return options.response, options.input
    
if __name__ == '__main__':
    response_file, input_file = processOptions()
    
    MI_instance = MIcalculator(response_file,input_file)
