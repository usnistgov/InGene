import sys, os
from optparse import OptionParser
from SEMIL import SEMIL

def processOptions():
    parser = OptionParser()
    #parser.add_option("-r", dest="response", help="Name of the file containing the inputs.", default="")
    parser.add_option("-i", dest="inputfile", help="Number of sample trajectories to be created.", default="")
    parser.add_option("-f", dest="flag", help="Flag for landscapce(l) or channel capacity(c).", default="b")

    [options, args] = parser.parse_args()

    InputFile = open(options.inputfile, 'r')
    InputLines = InputFile.readlines()

    return InputLines, options.flag

if __name__ == '__main__':
    input_lines, flag = processOptions()

    current_dir = os.getcwd()
    print(current_dir)

    MIL_instance = SEMIL(current_dir,input_lines,flag)

    if flag=='l':
        MIL_instance.get_MI_landscape()
    elif flag=='c':
        MIL_instance.compute_CC()
        MIL_instance.display_result()
        MIL_instance.write_results()
    else:
        MIL_instance.landscape_and_gradient()
