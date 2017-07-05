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


def get_data_using_inputs(inputs):

    """
    Parses your data according to inputs dictionary. Requires you have a 'file_set' in your inputs.
    Parameters
    ----------
    inputs : dict
        Dictionary of input information
    """
    complex_fluorescence = {}
    ligand_fluorescence = {}
        
    for protein in inputs['file_set'].keys():

        #concatenate xmls into one dictionary with all ligand data

        my_file = []

        # Data just helps define our keys, which will be our data sections.
        # We are importing the first two files, because sometimes we have two files per
        # data collection, e.g. if we wanted a lot of types of reads 
        # (there is a limited number of type of read on the infinite per script).
        data = platereader.read_icontrol_xml(inputs['file_set']['%s'%protein][0])
        try:
            data.update(platereader.read_icontrol_xml(inputs['file_set']['%s'%protein][1]))
        except:
            pass
        
        for file in inputs['file_set']['%s'%protein]:
            my_file.append(file)
            new_dict = platereader.read_icontrol_xml(file)
            for key in data:
                try:
                    data[key].update(new_dict[key])
                except:
                    pass

        # Are there any experiments the user wants to skip analyzing?
        skipped_experiments=[]
        for i, ligand in enumerate(inputs['ligand_order']):
            if ligand == 'skip':
                skipped_experiments.append(i*2)
            else:
                continue

        skipped_rows=[]
        for i in skipped_experiments:
            skipped_rows.append(ALPHABET[i])
            skipped_rows.append(ALPHABET[i+1])
        print("Skipping analysis of rows: ", skipped_rows)

        for i in range(0,len(inputs['ligand_order']*2),2):

            if i in skipped_experiments:
                continue
            else:
                protein_row = ALPHABET[i]
                buffer_row = ALPHABET[i+1]

                name = "%s-%s-%s%s"%(protein,inputs['ligand_order'][int(i/2)],protein_row,buffer_row)
 
                # for spectra assays
                if 'wavelength' in inputs:
        
                    complex_fluorescence_data = platereader.select_data(data, inputs['section'], protein_row, wavelength = '%s' %inputs['wavelength'])
                    ligand_fluorescence_data = platereader.select_data(data, inputs['section'], buffer_row, wavelength = '%s' %inputs['wavelength'])

                # for single wavelength assays
                else:
                
                    complex_fluorescence_data = platereader.select_data(data, inputs['section'], protein_row)
                    ligand_fluorescence_data = platereader.select_data(data, inputs['section'], buffer_row)
                
                complex_fluorescence[name] = reorder2list(complex_fluorescence_data,well)
                ligand_fluorescence[name] = reorder2list(ligand_fluorescence_data,well)

    return [complex_fluorescence, ligand_fluorescence]

 
