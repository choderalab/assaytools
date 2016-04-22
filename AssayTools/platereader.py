#!/usr/bin/env python

"""
Tools for assisting in reading and manipulating data from plate readers.

"""

#=============================================================================================
# Imports
#=============================================================================================

import numpy as np
import re
import string
from lxml import etree

#=============================================================================================
# Tecan Infinite plate reader helper functions
#=============================================================================================

def read_icontrol_xml(filename):
    """
    Read a Tecan iControl XML-formatted file and return all section data.

    Parameters
    ----------
    filename : str
       The name of the XML file to be read.
    section_name : str
       The 'Name' attribute of the section to read.

    Returns
    -------
    data : dict
       data[well_label] is the well reading for well 'well_label' (e.g. data['A4'] = 324)

    Examples
    --------

    """
    # Parse XML file into nodes.
    root_node = etree.parse(filename)

    # Build a dict of section nodes.
    section_nodes = { section_node.get('Name') : section_node for section_node in root_node.xpath("/*/Section") }

    def extract_well_data(section_node):
        """Return a dict of all wells in a section.

        Parameters
        ----------
        section_node : xml node
           The XML node for a given section of the iControl file.

        Returns
        -------
        well_data : dict
           Returns all values for each section as a dictionary of dictionaries
           (e.g. well_data['280_TopRead'] is {'A1': 62213.0,'A10': 10505.0,...})
           Above is an example return for singlet data xml file, if input is a spectra data 
           xml file well_data will also include a dictionary for each well where the keys are the wavelength
           (e.g. well_data['em280'] is {'A1': {'280': '3216344', '285': '2587710'...}...})

        """

        # Get a list of all well nodes.
        well_nodes = section_node.xpath("*/Well")

        # Process all wells into data.
        well_data = dict()
        for well_node in well_nodes:
            if well_node.get('Type') == 'Single':
                well_name = well_node.get('Pos')
                for r in well_node:
                    read = r.text
                well_data[well_name] = read
            else:
                for well_node in well_nodes:
                    well_name = well_node.get('Pos')
                    wavelength_data = dict()
                    for wave in well_node:
                        wavelength = wave.attrib['WL']
                        wavelength_data[wavelength] = wave.text
                    well_data[well_name] = dict()
                    well_data[well_name] = wavelength_data

        return well_data

    # Process all sections.
    sections = dict()
    for (section_name, section_node) in section_nodes.iteritems():
        well_data = extract_well_data(section_node)
        sections[section_name] = dict()
        sections[section_name] = well_data

    return sections

def select_data(file_data,section_name,selection,*args,**kwargs):
    """
    Read a Tecan iControl XML-formatted file and extract a particular part (row, column, 
    well, or well selection) for a particular section.
    Parameters
    ----------
    input : str or dict
       EITHER the name of the XML file to be read OR the dict of an XML file made by read_icontrol_xml
    section_name : str
       The 'Name' attribute of the section to read.
    selection : str
       The selection to extract, for example 'A', '1', or ['A1'].
    wavelength (optional arg): str
       The wavelength to extract, for example '480'.
    Returns
    -------
    data : dictionary
       (e.g. {'A1': 471.0, 'A2': 418.0})
       or dictionary of dictionaries if spectra with no wavelength selected
       (e.g. {'A1': {'280': '3216344','285': '2587710'...}}
    Examples
    --------
    gefitinib_abl_singlet_A1 = select_data(singlet_file, '350_TopRead', ['A1'])
    gefitinib_abl_singlet_A = select_data(singlet_file, '350_TopRead', 'A')
    gefitinib_abl_singlet_1 = select_data(singlet_file, '350_TopRead', '1')
    bosutinib_abl_spectra_A1 = select_data(spectra_file, 'em280', ['A1'])
    bosutinib_abl_spectra_A_480 = select_data(spectra_file, 'em280', 'A', wavelength='480')
    """
    wavelength = kwargs.get('wavelength', None)
    
    #construct all possible wells and columns for 96 or 384 well plate
    rows =[]
    cols =[]
    for i in string.ascii_uppercase:
        rows.append('%s' % i)
    for i in range(1,25):
        cols.append('%s' % i)
    
    # read data from xml or dict
    
    if type(file_data) == dict:
        well_data = file_data
    else:
        well_data = read_icontrol_xml(file_data)

    section_data = well_data[section_name]
    
    # extract selection
    data = dict()
    for select in selection:
        if select in section_data:                          # if individual wells
            data[select] = section_data[select]
        elif any([selection == r for r in rows]):           # if row
            for col in cols:
                try:
                    data[selection + col] = section_data[selection + col]
                except KeyError:
                    continue
        elif any([selection == c for c in cols]):            # if column
            for row in rows:
                try:
                    data[row + selection] = section_data[row + selection]
                except KeyError:
                    continue
        else:
            print 'bad selection'

    # extract wavelength (only relevant for spectra)
    # if you don't include wavelength for spectra, data will be a dict of all wavelengths
    if type(data.itervalues().next()) == dict and wavelength != None:
        new_data = dict()
        for key in data.keys():
            new_data[key] = data[key][wavelength]
        data = new_data
            
    return data

def read_emission_spectra_text(filename):
    """
    Read text-formatted emission spectra.

    Parameters
    ----------
    filename : str
       The Tecan Infinite output filen to be read.

    Returns
    -------
    SRC_280 : numpy.array
    SRC_280_x : numpy.array
    SRC_280_x_num : numpy.array

    Examples
    --------

    """

    SRC_280 = np.genfromtxt(filename, dtype='str')
    SRC_280_x = SRC_280[0,:]
    SRC_280_x_num = re.findall(r'\d+', str(SRC_280_x )[1:-1])

    return [SRC_280, SRC_280_x, SRC_280_x_num]
