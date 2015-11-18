# This script takes xml data file output from the Tecan Infinite m1000 Pro plate reader
# and makes quick and dirty images of the raw data.

#But with scans and not just singlet reads.

# The same procedure can be used to make matrices suitable for analysis using
# matrix = dataframe.values

# Made by Sonya Hanson, with some help from things that worked in xml2png.py
# Friday, June 20,2014

# Usage: python xml2png4scans.py *.xml

############ For future to combine with xml2png.py
#
#    for i, sect in enumerate(Sections):
#        reads = sect.xpath("*/Well")
#        parameters = root.xpath(path)[0]
#        if reads[0].attrib['Type'] == "Scan":
#
##############

import numpy as np
import matplotlib.pyplot as plt
from lxml import etree
import pandas as pd
import matplotlib.cm as cm
import sys
import os

# Define xml files.

xml_files = sys.argv[1:]

so_many = len(xml_files)
print "****This script is about to make png files for %s xml files. ****"  % so_many

# Define extract function that extracts parameters

def extract(taglist):
    result = []
    for p in taglist:
        print "Attempting to extract tag '%s'..." % p
        try:
            param = parameters.xpath("*[@Name='" + p + "']")[0]
            result.append( p + '=' + param.attrib['Value'])
        except:
            # tag not found
            result.append(None)

    return result

for file in xml_files:

    # Parse XML file.

    root = etree.parse(file)

    # Remove extension from xml filename.

    file_name = os.path.splitext(file)[0]

    # Extract plate type and barcode.

    plate = root.xpath("/*/Header/Parameters/Parameter[@Name='Plate']")[0]
    plate_type = plate.attrib['Value']

    try:
        bar = root.xpath("/*/Plate/BC")[0]
        barcode = bar.text
    except:
        barcode = 'no barcode'

    # Define Sections.

    Sections = root.xpath("/*/Section")
    much = len(Sections)
    print "****The xml file " + file + " has %s data sections:****" % much
    for sect in Sections:
        print sect.attrib['Name']

    for i, sect in enumerate(Sections):

        # Extract Parameters for this section.

        path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
        parameters = root.xpath(path)[0]

        # Parameters are extracted slightly differently depending on Absorbance or Fluorescence read.

        if  parameters[0].attrib['Value'] == "Absorbance":
            result = extract(["Mode", "Wavelength Start", "Wavelength End", "Wavelength Step Size"])
            title = '%s, %s, %s, %s' % tuple(result)

        else:
            result = extract(["Gain", "Excitation Wavelength", "Emission Wavelength", "Part of Plate", "Mode"])
            title = '%s, %s, %s, \n %s, %s' % tuple(result)

        print "****The %sth section has the parameters:****" %i
        print title

        # Extract Reads for this section.

        Sections = root.xpath("/*/Section")

        reads = root.xpath("/*/Section[@Name='" + sect.attrib['Name'] + "']/*/Well")

        wellIDs = [read.attrib['Pos'] for read in reads]

        data = [(float(s.text), float(s.attrib['WL']), r.attrib['Pos'])
                 for r in reads
                 for s in r]

        dataframe = pd.DataFrame(data, columns=['fluorescence','wavelength (nm)','Well'])

        dataframe_pivot = pd.pivot_table(dataframe, index = 'wavelength (nm)', columns= ['Well'])

        # Make plot, complete with separate png for each section.

        section_name = sect.attrib['Name']

        Alphabet = ['A','B','C','D','E','F','G','H']

        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
        for j,A in enumerate(Alphabet):
            for i in range(1,12):
                try:
                    dataframe_pivot.fluorescence.get(A + str(i)).plot(ax=axes[(j/3)%3,j%3], title=A, c=cm.hsv(i*15), ylim=[0,400])
                except:
                    print "****No row %s.****" %A

        fig.suptitle('%s \n %s \n Barcode = %s' % (title, plate_type, barcode), fontsize=14)
        fig.subplots_adjust(hspace=0.3)
        plt.savefig('%s_%s.png' % (file_name, section_name))
