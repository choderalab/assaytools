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
        sw_data = []

        # Data just helps define our keys, which will be our data sections.
        # We are importing the first two files, because sometimes we have two files per
        # data collection, e.g. if we wanted a lot of types of reads 
        # (there is a limited number of type of read on the infinite per script).
        data = platereader.read_icontrol_xml(inputs['file_set']['%s'%protein][0])
        try:
            data.update(platereader.read_icontrol_xml(inputs['file_set']['%s'%protein][1]))
        except:
            pass
        
        if 'single_well' in inputs:
            
            for file in inputs['file_set']['%s'%protein]:
                my_file.append(file)
                data_to_append = {key:None for key in data}
                new_dict = platereader.read_icontrol_xml(file)

                data_to_append.update(new_dict)

                sw_data.append(data_to_append) 
                #sw_data is a list of dataframes in the order of the list of files
            
            if inputs['protein_wells'][protein][0] not in range(1,25): #make sure not a two digit number
                if len(str(inputs['protein_wells'][protein][0]))==2 or 3: # is it e.g. 'A1' or 'A11'
                
                    if len(inputs['protein_wells'][protein]) == len(inputs['ligand_order']):  
                        pass

                    else:
                        print('***To define individual wells, you need to define the same number of wells as ligands.***')
                        break
            
            for i in range(0,len(inputs['ligand_order'])): # i defines ligand
                
                for j,protein_well_ID in enumerate(inputs['protein_wells'][protein]): # j defines protein
                
                    if str(protein_well_ID) in ALPHABET:
                        #print('Your ligands are in a row (1,2,3,4, etc. are different ligands)!')
                        protein_well = '%s%s'%(protein_well_ID,i)
                        buffer_well = '%s%s'%(inputs['buffer_wells'][protein][j],i)
                        
                        name = "%s-%s-%s%s"%(protein,inputs['ligand_order'][i],protein_well,buffer_well)
                        
                    elif len(str(protein_well_ID)) == 2 or 3:
                        if inputs['protein_wells'][protein][0] in range(1,25): # make sure not a two digit number
                            pass
                        else:
                            protein_well = protein_well_ID                         
                            buffer_well = inputs['buffer_wells'][protein][j]
                        
                            name = "%s-%s-%s%s"%(protein,inputs['ligand_order'][j],protein_well,buffer_well)
                        
                    else:
                        #print('Your ligands are in a column (A,B,C,D, etc. are different ligands)!')
                        protein_well = '%s%s'%(ALPHABET[i],protein_well_ID)
                        buffer_well = '%s%s'%(ALPHABET[i],inputs['buffer_wells'][protein][j])
                           
                        name = "%s-%s-%s%s"%(protein,inputs['ligand_order'][i],protein_well,buffer_well)
                        #print(name)

                    protein_data = []
                    buffer_data = []
                
                    for k in range(len(sw_data)): # k should define ligand concentration if your list of files was ordered correctly
                
                        # for spectra assays
                        if 'wavelength' in inputs:
                    
                            protein_data.append(float(platereader.select_data(sw_data[k], inputs['section'], [protein_well], 
                                                                    wavelength = '%s' %inputs['wavelength'])[protein_well]))
                            buffer_data.append(float(platereader.select_data(sw_data[k], inputs['section'], [buffer_well], 
                                                                    wavelength = '%s' %inputs['wavelength'])[buffer_well]))
        
                        # for single wavelength assays
                        else:
                
                            protein_data.append(float(platereader.select_data(sw_data[k], inputs['section'], [protein_well])[protein_well]))
                            buffer_data.append(float(platereader.select_data(sw_data[k], inputs['section'], [buffer_well])[buffer_well]))
    
                    complex_fluorescence[name] = np.asarray(protein_data)
                    ligand_fluorescence[name] = np.asarray(buffer_data)                 

        else:
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
                if ligand == None:
                    skipped_experiments.append(i*2)
                else:
                    continue

            skipped_rows=[]
            for i in skipped_experiments:
                skipped_rows.append(ALPHABET[i])
                skipped_rows.append(ALPHABET[i+1])
            if len(skipped_rows) != 0:
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
