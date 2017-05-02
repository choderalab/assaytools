#!/usr/bin/env python

"""
Tools for assisting in parsing data to input into analysis procedure.
"""

#=============================================================================================
# Imports
#=============================================================================================

from assaytools import platereader
import string
import numpy as np

#=============================================================================================
# Functions and things to help with parsing
#=============================================================================================

def reorder2list(data,well):

    sorted_keys = sorted(well.keys(), key=lambda k:well[k])

    reorder_data = []

    for key in sorted_keys:
        try:
            reorder_data.append(data[key])
        except:
            pass

    reorder_data = [r.replace('OVER','70000') for r in reorder_data]

    reorder_data = np.asarray(reorder_data,np.float64)

    return reorder_data

ALPHABET = string.ascii_uppercase

well = dict()
for j in string.ascii_uppercase:
    for i in range(1,25):
        well['%s' %j + '%s' %i] = i

#=============================================================================================
# Parsing functions
#=============================================================================================

# This function requires an inputs dictionary with a single xml file input.

def get_data_using_inputs(inputs):

    """
    Parameters
    ----------
    inputs : dict
        Dictionary of input information
    """
    
    data = platereader.read_icontrol_xml(inputs['my_file'])
    
    complex_fluorescence = {}
    ligand_fluorescence = {}
    
    for i in range(0,15,2):
        protein_row = ALPHABET[i]
        buffer_row = ALPHABET[i+1]

        name = "%s-%s%s"%(inputs['ligand_order'][int(i/2)],protein_row,buffer_row)

        complex_fluorescence_data = platereader.select_data(data, inputs['section'], protein_row)
        ligand_fluorescence_data = platereader.select_data(data, inputs['section'], buffer_row)

        complex_fluorescence[name] = reorder2list(complex_fluorescence_data,well)
        ligand_fluorescence[name] = reorder2list(ligand_fluorescence_data,well)

    return [complex_fluorescence, ligand_fluorescence]
 
