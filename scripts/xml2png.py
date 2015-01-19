# This script takes xml data file output from the Tecan Infinite m1000 Pro plate reader 
# and makes quick and dirty images of the raw data.

# The same procedure can be used to make matrices suitable for analysis using
# matrix = dataframe.values

# Made by Sonya Hanson with significant help from JH Prinz.
# Wednesday, May 7, 2014

# Usage: python xml2png.py *.xml

###################
# To Do
#
# Set fixed Y axis.
# Define Wells for Legend.
# Define x and y axes.
# Figure out @Name='Mode' thing.
#
###################

import numpy as np
import matplotlib.pyplot as plt
from lxml import etree
import pandas as pd
import sys
import os

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

# Define get_wells_from_section function that extracts the data from each Section.
# It is written sort of strangely to ensure data is connected to the correct well.

def get_wells_from_section(path):
    reads = path.xpath("*/Well")
    wellIDs = [read.attrib['Pos'] for read in reads]

    data = [(float(s.text), r.attrib['Pos'])
         for r in reads
         for s in r]

    datalist = {
      well : value
      for (value, well) in data
    }

    welllist = [
                [
                 datalist[chr(64 + row) + str(col)]
                 if chr(64 + row) + str(col) in datalist else None
                 for row in range(1,17)
                ]
                for col in range(1,24)
                ]

    return welllist

def process_files(xml_files):
    """
    Main entry point.

    """
    # Define xml files.
    so_many = len(xml_files)
    print "****This script is about to make png files for %s xml files. ****"  % so_many

    for file in xml_files:

        # Parse XML file.
        root = etree.parse(file)

        # Remove extension from xml filename.
        file_name = os.path.splitext(file)[0]

        # Define Sections.
        Sections = root.xpath("/*/Section")
        much = len(Sections)
        print "****The xml file " + file + " has %s data sections:****" % much
        for sect in Sections:
            print sect.attrib['Name']

        data = []

        for i, sect in enumerate(Sections):
            # Extract Parameters for this section.
            path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
            parameters = root.xpath(path)[0]

            # Parameters are extracted slightly differently depending on Absorbance or Fluorescence read.
            if parameters[0].attrib['Value'] == "Absorbance":
                result = extract(["Mode", "Wavelength", "Part of Plate"])
                title = '%s, %s, %s' % tuple(result)

            else:
                result = extract(["Gain", "Excitation Wavelength", "Emission Wavelength", "Part of Plate", "Mode"])
                title = '%s, %s, %s, \n %s, %s' % tuple(result)

            print "****The %sth section has the parameters:****" %i
            print title

            # Extract Reads for this section.
            Sections = root.xpath("/*/Section")
            welllist = get_wells_from_section(sect)
            data.append({
                    'filename' : file_name,
                    'title' : title,
                    'dataframe' : pd.DataFrame(welllist)
                    })

        # Make plot, complete with subfigure for each section.
        fig, axes = plt.subplots(nrows=1,ncols=len(Sections), figsize=(20,4.5))

        for i, sect in enumerate(data):
            sect['dataframe'].plot(title = sect['title'] , ax = axes[i], legend=False )

        fig.tight_layout()
        fig.subplots_adjust(top=0.8)
        fig.suptitle("%s" % file_name, fontsize=18)

        plt.savefig('%s.png' % file_name)

    return

def entry_point():
    xml_files = sys.argv[1:]
    process_files(xml_files)

if __name__ == '__main__':
    xml_files = sys.argv[1:]
    process_files(xml_files)

