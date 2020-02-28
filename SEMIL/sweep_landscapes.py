import os
from optparse import OptionParser
import sys
from SEMIL import SEMIL
import copy

def processOptions():
    parser = OptionParser()
    #parser.add_option("-r", dest="response", help="Name of the file containing the inputs.", default="")
    parser.add_option("-i", dest="inputfile", help="Input file for landscapes", default="")
    parser.add_option("-f", dest="flag", help="Flag for landscapce(l) or channel capacity(c)", default="l")
    parser.add_option("-s", dest="subsizes", help="string of subsample sizes (e.g. 1,2,5,10)", default="1")
    parser.add_option("-r", dest="replicates", help="number of replicates of each subsample size", default="1")

    [options, args] = parser.parse_args()

    InputFile = open(options.inputfile, 'r')
    InputLines = InputFile.readlines()

    sub_list = options.subsizes.rstrip('\r\n').split(',')

    return InputLines, options.flag, sub_list, int(options.replicates)

if __name__ == '__main__':
    input_lines, flag, data_fractions, reps = processOptions()

    current_dir = os.getcwd()
    print(current_dir)

    for l in input_lines:
        if 'Response file' in l:
            response_file = l.split(' = ')[1]
            response_file = response_file.rstrip('\r\n').split('.')[0]
        elif 'Casename' in l:
            casename = l.split(' = ')[1]
            casename = casename.rstrip('\r\n')

    old_response_tag = response_file
    old_case_tag = casename

    #data_fractions = [1,2,5,10]#,20]
    samples = list(range(1,reps+1))

    for df_i in range(0,len(data_fractions)):
        #sbound = min(len(samples),data_fractions[df_i])
        sbound = len(samples)
        for s_i in range(0,sbound):
            df_str = data_fractions[df_i]

            # if df>=1:
            #     df_str = str(int(df))
            # else:
            #     df_s = str(df).split('.')
            #     df_str = df_s[0]+'p'+df_s[1]

            s = samples[s_i]

            # new_exp_tag = 'expressions'+df_str+'_'+str(s)
            # new_trial_tag = 'trial'+df_str+'_'+str(s)

            new_response_tag = response_file+df_str+'_'+str(s)

            new_case_tag = casename+df_str+'_'+str(s)

            new_input_lines = []

            for l in input_lines:
                #this_l = l.rstrip('\r\n')
                if old_response_tag in l:
                    l = l.replace(old_response_tag,new_response_tag)

                if old_case_tag in l:
                    l = l.replace(old_case_tag,new_case_tag)

                new_input_lines.append(l)

            input_lines = copy.deepcopy(new_input_lines)

            # if df_i>0 or s_i>0:
            #     new_input_lines = []
            #
            #     for l in input_lines:
            #         #this_l = l.rstrip('\r\n')
            #         if old_exp_tag in l:
            #             l = l.replace(old_exp_tag,new_exp_tag)
            #
            #         if old_trial_tag in l:
            #             l = l.replace(old_trial_tag,new_trial_tag)
            #
            #         new_input_lines.append(l)
            #
            #     input_lines = copy.deepcopy(new_input_lines)

            SEMIL_instance = SEMIL(current_dir,input_lines,flag)

            if df_i>0 or s_i>0:
                SEMIL_instance.set_SROM_wts(wt_dict,wt_locs)
                #print(SEMIL_instance.wt_dict)
                #sys.stdout.flush()

            if flag=='l':
                SEMIL_instance.get_MI_landscape()

            if df_i==0 and s_i==0:
                wt_dict = copy.deepcopy(SEMIL_instance.wt_dict)
                wt_locs = copy.deepcopy(SEMIL_instance.SROM_obj.sample_locations)

            # f = open(input_file,'r')
            # ilines = f.readlines()
            # f.close()
            #
            # if df_i>0 or s_i>0:
            #     f = open(duplicate_file,'w')
            #
            #     for l in ilines:
            #         this_l = l.rstrip('\r\n')
            #         if old_exp_tag in this_l:
            #             this_l = this_l.replace(old_exp_tag,new_exp_tag)
            #
            #         if old_trial_tag in this_l:
            #             this_l = this_l.replace(old_trial_tag,new_trial_tag)
            #
            #         print(this_l,file=f)
            #
            #     f.close()
            #
            #     os.system('mv '+duplicate_file+' '+input_file)
            #
            # os.system('python3 $INGENE/SEMIL/landmaker.py -i '+input_file+' -f '+flag)

            print(new_response_tag,' completed.')

            old_response_tag = new_response_tag
            old_case_tag = new_case_tag

    #os.system('rm rinput.txt')
