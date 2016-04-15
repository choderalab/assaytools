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
           well_data[well_name] is the value associated with well_name (e.g. well_data['D7'] = 423)

        """

        # Get a list of all well nodes.
        well_nodes = section_node.xpath("*/Well")

        # Process all wells into data.
        well_data = dict()
        for well_node in well_nodes:
            if well_node.get('Type') == 'Single':
                well_name = well_node.get('Pos')
                value = well_node.xpath("string()")
                well_data[well_name] = float(value)
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

def extract_data(filename,section_name,selection,*args,**kwargs):
    """
    Read a Tecan iControl XML-formatted file and extract a particular part (row, column, 
    well, or well selection) for a particular section.
    Parameters
    ----------
    filename : str
       The name of the XML file to be read.
    section_name : str
       The 'Name' attribute of the section to read.
    selection : str
       The selection to extract, for example 'A', '1', or ['A1'].
    wavelength (optional arg): str
       The wavelength to extract, for example '480'.
    Returns
    -------
    data : list of lists 
       or list of dictionaries if spectra with no wavelength selected
       data[0] is data and data[1] is the selection (e.g. [[471.0, 418.0], ['A1', 'A2']])
    Examples
    --------
    >>> gefitinib_abl_singlet_A1 = extract_data(singlet_file, '350_TopRead', ['A1'])
    >>> gefitinib_abl_singlet_A = extract_data(singlet_file, '350_TopRead', 'A')
    >>> gefitinib_abl_singlet_1 = extract_data(singlet_file, '350_TopRead', '1')
    >>> bosutinib_abl_spectra_A1 = extract_data(spectra_file, 'em280', ['A1'])
    >>> bosutinib_abl_spectra_A_480 = extract_data(spectra_file, 'em280', 'A', wavelength='480')
    """
    wavelength = kwargs.get('wavelength', None)
    
    #construct all possible wells and columns for 96 or 384 well plate
    rows =[]
    cols =[]
    for i in string.ascii_uppercase:
        rows.append('%s' % i)
    for i in range(1,25):
        cols.append('%s' % i)
    
    # import data from section
    well_data = read_icontrol_xml(filename)
    section_data = well_data[section_name]
    
    # extract selection
    data = []
    for select in selection:
        if select in section_data:                          # if individual wells
            data = [section_data[sele] for sele in selection]
        elif any([selection == r for r in rows]):           # if row
            for col in cols:
                try:
                    data.append(section_data[selection + col])
                except KeyError:
                    continue
        elif any([selection == c for c in cols]):            # if column
            for row in rows:
                try:
                    data.append(section_data[row + selection])
                except KeyError:
                    continue
        else:
            print 'bad selection'

    # extract wavelength (only relevant for spectra)
    # if you don't include wavelength for spectra, data[0] will be a dict of all wavelengths
    if type(data[0]) == dict and wavelength != None:
        new_data = []
        for i in range( len(data) ):
            new_data.append(data[i][wavelength])
        data = new_data
            
    return [data,selection]

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
