# Generate Numpy array of stated ligand concentration(Lstated) for logarithmic dilution along a row.
# 
# Example output:
# np.array([20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0001e-6], np.float64)
#
# Usage:	
#	python calculate_Lstated_array.py --n_wells 12 --h_conc 8e-06 --l_conc 2.53e-09 --dilution logarithmic
#       python calculate_Lstated_array.py --n_wells 6 --h_conc 100 --l_conc 10 --dilution linear
 

from __future__ import print_function, division
import numpy as np
import argparse
import sys

# Argument parser
parser = argparse.ArgumentParser(prog="calculate_Lstated_array")
parser.add_argument("--n_wells", dest="n_wells", help="Enter target number of different concentrations.", type=int)
parser.add_argument("--h_conc", dest= "highest_conc", help="The highest concentration in the series (Unit: M).", type=float)
parser.add_argument("--l_conc", dest="lowest_conc", help="The lowest concentration in the series (Unit: M).", type=float)
parser.add_argument("--dilution", dest="dilution_method", help="Defines logarithmic or linear dilution.Enter 'logarithmic' or 'linear.'.", type=str)
args = parser.parse_args()

n = args.n_wells
print("Number of wells: ", n)
highest_conc = args.highest_conc
print("Highest concentration (M): ", highest_conc)
lowest_conc = args.lowest_conc
print("Lowest concentration (M): ", lowest_conc)
dilution_method = args.dilution_method
print("Dilution method is {}.".format(dilution_method))


if dilution_method == "logarithmic":
   
     # Calculate dilution factor
    DF = (highest_conc / lowest_conc)**(1.0/(n-1))
    print("Dilution factor: {}".format(DF))

    # Calculate concentration of solute in each well
    Ltot = np.zeros(n)
    for i, conc  in enumerate(Ltot):
        Ltot[i] = highest_conc/(DF**i)

elif dilution_method == "linear":
   
    # Calculate interval, i.e. concentration difference between consecutive dilution points
    interval = (highest_conc - lowest_conc)/(n-1)
    print("Linear concentration interval(M): {} ".format(interval) )

    # Calculate linear dilution array
    Ltot = np.zeros(n)
    for i, conc  in enumerate(Ltot):
        Ltot[i] = highest_conc - i*interval

else:
    
    sys.exit("Error in --dilution argument! It must be specified as '--dilution linear' or '--dilution logarithmic' .")


# Create a string that codes for Lstated array:
array_string = "np.array(["
for i in range(n-1):
    array_string = array_string + str(Ltot[i]) + ","
array_string = array_string + str(Ltot[n-1]) + "], np.float64)"

# Print with matching format for input.py file
print("Copy this line to Lstated section of assaytools input.py file: ")
print(array_string)
