#!/usr/bin/env python

"""
Tools for assisting in reading and manipulating data from plate readers.

"""

#=============================================================================================
# Imports
#=============================================================================================

import numpy as np
import re
import pint

#=============================================================================================
# Tecan Infinite plate reader helper functions
#=============================================================================================

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
