# This script takes xml data file output from the Tecan Infinite m1000 Pro plate reader
# and makes quick and dirty images of the raw data.

# But with scans and not just singlet reads.
# This script specifically combines four spectrum scripts (AB, CD, EF, GH) into a single dataframe and plot.

# The same procedure can be used to make matrices suitable for analysis using
# matrix = dataframe.values

# Made by Sonya Hanson, with some help from things that worked in xml2png.py and xml2png4scans.py
# Friday, November 18,2015

# Usage: python xml2png4scans-spectra.py *.xml

############ For future to combine with xml2png.py
#
#    for i, sect in enumerate(Sections):
#        reads = sect.xpath("*/Well")
#        parameters = root.xpath(path)[0]
#        if reads[0].attrib['Type'] == "Scan":
#
##############

import matplotlib.pyplot as plt
from lxml import etree
import pandas as pd
import matplotlib.cm as cm
import seaborn
import sys
import os

### Define extract function that extracts parameters

def extract(taglist):
    result = []
    for p in taglist:
        print "Attempting to extract tag '%s'..." % p
        try:
            param = parameters.xpath("*[@Name='" + p + "']")[0]
            result.append( p + '=' + param.attrib['Value'])
        except:
            ### tag not found
            result.append(None)

    return result

### Define an initial set of dataframes, one per each section

large_dataframe0 = pd.DataFrame()
large_dataframe1 = pd.DataFrame()
large_dataframe2 = pd.DataFrame()

def process_files(xml_files):
    """
    Main entry point.
    """

    ### Define xml files.

    xml_files = sys.argv[1:]

    so_many = len(xml_files)
    print "****This script is about to make png files for %s xml files. ****"  % so_many

    for file in xml_files:

        ### Parse XML file.

        root = etree.parse(file)

        ### Remove extension from xml filename.

        file_name = os.path.splitext(file)[0]

        ### Extract plate type and barcode.

        plate = root.xpath("/*/Header/Parameters/Parameter[@Name='Plate']")[0]
        plate_type = plate.attrib['Value']

        try:
            bar = root.xpath("/*/Plate/BC")[0]
            barcode = bar.text
        except:
            barcode = 'no barcode'

        ### Define Sections.

        Sections = root.xpath("/*/Section")
        much = len(Sections)
        print "****The xml file " + file + " has %s data sections:****" % much
        for sect in Sections:
            print sect.attrib['Name']

        for i, sect in enumerate(Sections):

            ### Extract Parameters for this section.

            path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
            parameters = root.xpath(path)[0]

            ### Parameters are extracted slightly differently depending on Absorbance or Fluorescence read.
            # Attach these to title1, title2, or title3, depending on section which will be the same for all 4 files.

            if  parameters[0].attrib['Value'] == "Absorbance":
                result = extract(["Mode", "Wavelength Start", "Wavelength End", "Wavelength Step Size"])
                globals()["title"+str(i)] = '%s, %s, %s, %s' % tuple(result)

            else:
                result = extract(["Gain", "Excitation Wavelength", "Emission Wavelength", "Part of Plate", "Mode"])
                globals()["title"+str(i)] = '%s, %s, %s, \n %s, %s' % tuple(result)

            print "****The %sth section has the parameters:****" %i
            print globals()["title"+str(i)]

            ### Extract Reads for this section.

            Sections = root.xpath("/*/Section")

            reads = root.xpath("/*/Section[@Name='" + sect.attrib['Name'] + "']/*/Well")

            wellIDs = [read.attrib['Pos'] for read in reads]

            data = [(s.text, float(s.attrib['WL']), r.attrib['Pos'])
                     for r in reads
                     for s in r]

            dataframe = pd.DataFrame(data, columns=['fluorescence','wavelength (nm)','Well'])
        
            ### dataframe_rep replaces 'OVER' (when fluorescence signal maxes out) with '3289277', an arbitrarily high number

            dataframe_rep = dataframe.replace({'OVER':'3289277'})

            dataframe_rep[['fluorescence']] = dataframe_rep[['fluorescence']].astype('float')

            ### Create large_dataframe1, large_dataframe2, and large_dataframe3 that collect data for each section
            ### as we run through cycle through sections and files.

            globals()["dataframe_pivot"+str(i)] = pd.pivot_table(dataframe_rep, index = 'wavelength (nm)', columns= ['Well'])
        
            print 'The max fluorescence value in this dataframe is %s'% globals()["dataframe_pivot"+str(i)].values.max()

            globals()["large_dataframe"+str(i)] = pd.concat([globals()["large_dataframe"+str(i)],globals()["dataframe_pivot"+str(i)]])

    ### Plot, making a separate png for each section.

    for i, sect in enumerate(Sections):

        section_name = sect.attrib['Name']
    
        path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
        parameters = root.xpath(path)[0]
    
        if  parameters[0].attrib['Value'] == "Absorbance":
            section_ylim = [0,0.2]
        else:
            section_ylim = [0,40000]

        Alphabet = ['A','B','C','D','E','F','G','H']

        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
        for j,A in enumerate(Alphabet):
            for k in range(1,12):
                try:
                    globals()["large_dataframe"+str(i)].fluorescence.get(A + str(k)).plot(ax=axes[(j/3)%3,j%3], title=A, c=cm.hsv(k*15), ylim=section_ylim, xlim=[240,800])
                except:
                    print "****No row %s.****" %A

        fig.suptitle('%s \n %s \n Barcode = %s' % (globals()["title"+str(i)], plate_type, barcode), fontsize=14)
        fig.subplots_adjust(hspace=0.3)
        plt.savefig('%s_%s.png' % (file_name, section_name))
        

    return

def entry_point():
    xml_files = sys.argv[1:]
    process_files(xml_files)

if __name__ == '__main__':
    xml_files = sys.argv[1:]
    process_files(xml_files)
