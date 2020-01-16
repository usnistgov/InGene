import os
from optparse import OptionParser
import sys

def processOptions():
    parser = OptionParser()
    #parser.add_option("-r", dest="response", help="Name of the file containing the inputs.", default="")
    parser.add_option("-i", dest="inputfile", help="Input file for landscapes.", default="")
    parser.add_option("-f", dest="flag", help="Flag for landscapce(l) or channel capacity(c).", default="b")

    [options, args] = parser.parse_args()

    # InputFile = open(options.inputfile, 'r')
    # InputLines = InputFile.readlines()

    return options.inputfile, options.flag

if __name__ == '__main__':
    input_file, flag = processOptions()

    data_fractions = [1,2,5,10]#,20]
    samples = list(range(1,6))

    duplicate_file = 'copy_'+input_file

    for df_i in range(0,len(data_fractions)):
        #sbound = min(len(samples),data_fractions[df_i])
        sbound = len(samples)
        for s_i in range(0,sbound):
            df = data_fractions[df_i]

            if df>=1:
                df_str = str(int(df))
            else:
                df_s = str(df).split('.')
                df_str = df_s[0]+'p'+df_s[1]

            s = samples[s_i]

            new_exp_tag = 'expressions'+df_str+'_'+str(s)

            new_trial_tag = 'trial'+df_str+'_'+str(s)

            f = open(input_file,'r')
            ilines = f.readlines()
            f.close()

            if df_i>0 or s_i>0:
                f = open(duplicate_file,'w')

                for l in ilines:
                    this_l = l.rstrip('\r\n')
                    if old_exp_tag in this_l:
                        this_l = this_l.replace(old_exp_tag,new_exp_tag)

                    if old_trial_tag in this_l:
                        this_l = this_l.replace(old_trial_tag,new_trial_tag)

                    print(this_l,file=f)

                f.close()

                os.system('mv '+duplicate_file+' '+input_file)

            os.system('python3 $INGENE/SEMIL/landmaker.py -i '+input_file+' -f '+flag)

            print(new_exp_tag,' completed.')

            old_exp_tag = new_exp_tag
            old_trial_tag = new_trial_tag

    #os.system('rm rinput.txt')
