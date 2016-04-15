#!/usr/bin/env python

"""
Tools for assisting in reading and manipulating data from plate readers.

"""

#=============================================================================================
# Imports
#=============================================================================================

import numpy as np
import re
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
        sections[section_name]['well_data'] = well_data

    return sections

def identify_rows_and_cols(well_data):
    """
    Build a sorted list of unique rows and columns.

    Parameters
    ----------
    well_data : dict
       well_data[well_name] is the value associated with well_name (e.g. well_data['D7'] = 423)

    Returns
    -------
    rows : list of str
       Sorted list of unique rows.
    cols : list of str
       Sorted listof unique columns.

    Examples
    --------

    >>> well_data = { 'A1' : 1, 'A2' : 2, 'A3' : 3, 'B1' : 4, 'B2' : 5, 'B3', 6 }
    >>> [rows, cols] = identify_rows_and_cols(well_data)

    """

    rows = set()
    cols = set()

    for (well_name, value) in well_data.iteritems():
        # Parse well names like 'A7' and 'D8'
        m = re.match(r"([A-Z]+)([0-9]+)", well_name)

        row = m.group(1)
        col = int(m.group(2))

        rows.add(row)
        cols.add(col)

    rows = sorted(list(rows))
    cols = [ '%s' % entry for entry in sorted(list(cols)) ] # convert back to str after sorting numerically

    return [rows, cols]

def extract_rows(well_data, rows, cols):
    """
    Examples
    --------

    >>> well_data = { 'A1' : 1, 'A2' : 2, 'A3' : 3, 'B1' : 4, 'B2' : 5, 'B3', 6 }
    >>> [rows, cols] = identify_rows_and_cols(well_data)
    >>> row_data = extract_rows(well_data, rows, cols)

    """

    row_data = dict()
    for row in rows:
        row_data[row] = np.zeros([len(cols)], np.float64)
        for (index, col) in enumerate(cols):
            row_data[row][index] = well_data[row + col]
    return row_data

def extract_cols(well_data, rows, cols):
    """
    Examples
    --------

    >>> well_data = { 'A1' : 1, 'A2' : 2, 'A3' : 3, 'B1' : 4, 'B2' : 5, 'B3', 6 }
    >>> [rows, cols] = identify_rows_and_cols(well_data)
    >>> col_data = extract_cols(well_data, rows, cols)

    """

    col_data = dict()
    for col in cols:
        col_data[col] = np.zeros([len(rows)], np.float64)
        for (index, row) in enumerate(rows):
            col_data[col][index] = well_data[row + col]
    return col_data

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
