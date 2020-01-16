import sys, os
from optparse import OptionParser
from SROMgenerator import SROMgenerator

def processOptions():
    totalList = []
    parser = OptionParser()
    parser.add_option("-i", dest="inputfile", help="Name of the file containing the inputs for a photon transport simulation. Use option -f to find the format.", default="")
    parser.add_option("-f", dest="inputFormat", help="Use 'show' to see the input format. Or specify the name to create an empty inputfile with the placeholders.", default="")
    
    [options, args] = parser.parse_args()

    if len(options.inputfile) != 0:
        try:
            InputFile = open(options.inputfile, 'r')
        except IOError:
            print >> sys.stderr , "ERROR : Cannot open inputfile. Check inputfile name."
            
        InputLines = InputFile.readlines()
    else:
        print >> sys.stderr , "ERROR : Please try 'python getSROM.py -h' to learn about correct usage"

    return InputLines
    
if __name__ == '__main__':
    totalList = processOptions()
    print totalList
        
    if len(totalList) != 0:
        srom_obj = SROMgenerator(totalList)
        srom_obj.initialize()
        srom_obj.computeSROM()
#        scatObj.checkGradient(scatObj.x_SROM,scatObj.p_SROM)
#
#        finalSize = scatObj.sampleSize
#        for size in range(5,finalSize+1,5):
#            scatObj.sampleSize = size
#            scatObj.initialize()
#            scatObj.computeSROM()