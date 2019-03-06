# This script allows you to get the mean of a parameter,
# for example F_PL, from a previous run of quickmodel.py

# Usage:
#  python get_mean.py --input 'Src-Bosutinib-AB_mcmc_1.pickle'
#  or
#  python get_mean.py --input 'Src-Bosutinib-AB_mcmc_1.pickle' --parameter 'F_PL'

# Sonya M. Hanson

import json
import cPickle
import os

def get_mean(data_file, parameter='F_PL'):
    """
    Quick model for both spectra and single wavelength experiments
    ----------
    input_file : file_name
        File_name of pickle file to get mean from.
    parameter : string, default='F_PL'
        Parameter to get mean of. 
    """

    # load data from pickle file
    with open(r'%s'%data_file,'rb') as my_file:
        data = cPickle.load(my_file)

    # get t_equil from json file
    base_name = os.path.splitext(data_file)[0]
    json_base = base_name.replace('_mcmc','')
    json_file = json_base + '.json'

    with open(json_file) as f:
        json_data = json.load(f)
   
    t_equil = json_data['t_equil']

    # find mean of desired parameter removing data before equilibration
    my_mean = data[parameter][0][t_equil:].mean()

    # print and save this for later use
    print('The mean of parameter %s from data file %s is:' %(parameter,data_file))
    print(my_mean)

    mean_dict = {}
    mean_dict['data_file'] = data_file
    mean_dict[parameter] = my_mean

    with open('mean-%s.json'%parameter, 'w') as outfile:
                json.dump(mean_dict, outfile, sort_keys = True, indent = 4, ensure_ascii=False)

    print('Output file: mean-%s.json'%parameter)

def entry_point():

    import argparse

    parser = argparse.ArgumentParser(description="""Get mean from quickmodel run. Default parameter is F_PL.""")
    parser.add_argument("--input", help="The pickle file output by quickmodel for your experiment.",default=None)
    parser.add_argument("--parameter", help="The parameter you want to get the mean for.",default='F_PL')
    args = parser.parse_args()

    get_mean(data_file=args.input,parameter=args.parameter)

if __name__ == '__main__':
    entry_point()
