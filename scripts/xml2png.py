# A somewhat ugly, utilitarian script takes xml data file output from the Tecan Infinite m1000 Pro
# plate reader and allows for the quick visual inspection of raw data.
#
# Usage: python xml2png.py *.xml

# import math, xml, and dataframe libraries
import numpy as np
from lxml import etree
import pandas as pd
import string

# import libraries for making figures
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn

#import assaytools helper
from assaytools import platereader

# import libraries for interacting with arguments
import sys
import os
import argparse

# Define argparse stuff

parser = argparse.ArgumentParser(description="""Visualize your raw data by making a png from your xml:
> python xml2png.py --type spectra *.xml""")
parser.add_argument("files", nargs='*', help="xml file(s) to analyze")
parser.add_argument("--type", help="type of data file (spectra, singlet_96, singlet_384, scan)", choices=['spectra', 'singlet_96', 'singlet_384','scan'],default='singlet_96')
args = parser.parse_args()
print(args.files)
print("*** --type: analyzing %s file(s) ***" % args.type)

### Define extract function that extracts parameters

def extract(taglist, parameters): #["Mode", "Wavelength Start", "Wavelength End", "Wavelength Step Size"]
    result = []
    for p in taglist:
        print("Attempting to extract tag '%s'..." % p)
        try:
            param = parameters.xpath("*[@Name='" + p + "']")[0]
            result.append( p + '=' + param.attrib['Value'])
        except:
            ### tag not found
            result.append(None)

    return result

#############################################
 ###               singlet_96            ###
#############################################

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

def process_files_five(xml_files):
    """
    Main entry point.
    """

    so_many = len(xml_files)
    print("****This script is about to make png files for %s xml files. ****"  % so_many)

    for file in xml_files:

        # Parse XML file.

        root = etree.parse(file)

        # Remove extension from xml filename.

        file_name = os.path.splitext(file)[0]

        # Define Sections.

        Sections = root.xpath("/*/Section")
        much = len(Sections)
        print("****The xml file " + file + " has %s data sections:****" % much)
        for sect in Sections:
            print(sect.attrib['Name'])

        data = []

        for i, sect in enumerate(Sections):

           # Extract Parameters for this section.

            path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
            parameters = root.xpath(path)[0]

            # Parameters are extracted slightly differently depending on Absorbance or Fluorescence read.

            if  parameters[0].attrib['Value'] == "Absorbance":
                result = extract(["Mode", "Wavelength", "Part of Plate"], parameters)
                title = '%s, %s, %s' % tuple(result)

            else:
                result = extract(["Gain", "Excitation Wavelength", "Emission Wavelength", "Part of Plate", "Mode"], parameters)
                title = '%s, %s, %s, \n %s, %s' % tuple(result)

            print("****The %sth section has the parameters:****" % i)
            print(title)

            # Extract Reads for this section.

            Sections = root.xpath("/*/Section")

            welllist = get_wells_from_section(sect)

            data.append(
                {
                        'filename' : file_name,
                        'title' : title,
                        'dataframe' : pd.DataFrame(welllist, columns=list('ABCDEFGHIJKLMNOP'))
                }
            )

        # Make plot, complete with subfigure for each section.
        seaborn.set_palette("Paired", 10)
        seaborn.set_context("notebook", rc={"lines.linewidth": 2.5})

        fig, axes = plt.subplots(nrows=1,ncols=len(Sections), figsize=(60,4.5))

        for i, sect in enumerate(data):
            sect['dataframe'].plot(title = sect['title'] , ax = axes[i] )

        fig.tight_layout()
        fig.subplots_adjust(top=0.8)
        fig.suptitle("%s" % file_name, fontsize=18)

        plt.savefig('%s.png' % file_name)
        print('Look at how pretty your data is: %s.png' % file_name)

    return

#############################################
 ###               SPECTRA               ###
#############################################

### Define an initial set of dataframes, one per each section

large_dataframe0 = pd.DataFrame()
large_dataframe1 = pd.DataFrame()
large_dataframe2 = pd.DataFrame()
large_dataframe3 = pd.DataFrame()
large_dataframe4 = pd.DataFrame()
large_dataframe5 = pd.DataFrame()

def process_files_spectra(xml_files):
    """
    Main entry point.
    """

    so_many = len(xml_files)
    print("****This script is about to make png files for %s xml files. ****"  % so_many)

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
        print("****The xml file " + file + " has %s data sections:****" % much)
        for sect in Sections:
            print(sect.attrib['Name'])

        for i, sect in enumerate(Sections):

            ### Extract Parameters for this section.

            path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
            parameters = root.xpath(path)[0]

            ### Parameters are extracted slightly differently depending on Absorbance or Fluorescence read.
            # Attach these to title1, title2, or title3, depending on section which will be the same for all 4 files.

            if  parameters[0].attrib['Value'] == "Absorbance":
                result = extract(["Mode", "Wavelength Start", "Wavelength End", "Wavelength Step Size"], parameters)
                globals()["title"+str(i)] = '%s, %s, %s, %s' % tuple(result)

            else:
                result = extract(["Gain", "Excitation Wavelength", "Emission Wavelength", "Part of Plate", "Mode"], parameters)
                globals()["title"+str(i)] = '%s, %s, %s, \n %s, %s' % tuple(result)

            print("****The %sth section has the parameters:****" % i)
            print(globals()["title"+str(i)])

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

            print('The max fluorescence value in this dataframe is %s'% globals()["dataframe_pivot"+str(i)].values.max())

            globals()["large_dataframe"+str(i)] = pd.concat([globals()["large_dataframe"+str(i)],globals()["dataframe_pivot"+str(i)]])

            #print globals()["large_dataframe"+str(i)]

    ### Plot, making a separate png for each section.

    for i, sect in enumerate(Sections):

        section_name = sect.attrib['Name']

        path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
        parameters = root.xpath(path)[0]

        if  parameters[0].attrib['Value'] == "Absorbance":
            section_ylim = [0,0.2]
        else:
            section_ylim = [0,50000]

        Alphabet = ['A','B','C','D','E','F','G','H']

        fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
        for j,A in enumerate(Alphabet):
            for k in range(1,12):
                try:
                    globals()["large_dataframe"+str(i)].fluorescence.get(A + str(k)).plot(ax=axes[(j/3)%3,j%3], title=A, c=cm.hsv(k*15), ylim=section_ylim, xlim=[300,600])
                except:
                    print("****No row %s.****" %A)

        fig.suptitle('%s \n %s \n Barcode = %s' % (globals()["title"+str(i)], plate_type, barcode), fontsize=14)
        fig.subplots_adjust(hspace=0.3)
        plt.savefig('%s_%s.png' % (file_name, section_name))

    return

#############################################
 ###               SCAN                  ###
#############################################

def process_files_scan(xml_files):

    so_many = len(xml_files)
    print("****This script is about to make png files for %s xml files. ****"  % so_many)

    for file in xml_files:

        # Parse XML file.
        root = etree.parse(file)

        # Remove extension from xml filename.
        file_name = os.path.splitext(file)[0]

        # Extract plate type and barcode.
        plate = root.xpath("/*/Header/Parameters/Parameter[@Name='Plate']")[0]
        plate_type = plate.attrib['Value']

        bar = root.xpath("/*/Plate/BC")[0]
        barcode = bar.text

        # Define Sections.
        Sections = root.xpath("/*/Section")
        much = len(Sections)
        print("****The xml file " + file + " has %s data sections:****" % much)
        for sect in Sections:
            print(sect.attrib['Name'])

        data = []

        for i, sect in enumerate(Sections):

            # Extract Parameters for this section.
            path = "/*/Section[@Name='" + sect.attrib['Name'] + "']/Parameters"
            parameters = root.xpath(path)[0]

            # Parameters are extracted slightly differently depending on Absorbance or Fluorescence read.
            if  parameters[0].attrib['Value'] == "Absorbance":
                result = extract(["Mode", "Wavelength Start", "Wavelength End", "Wavelength Step Size"], parameters)
                title = '%s, %s, %s, %s' % tuple(result)

            else:
                result = extract(["Gain", "Excitation Wavelength", "Emission Wavelength", "Part of Plate", "Mode"], parameters)
                title = '%s, %s, %s, \n %s, %s' % tuple(result)

            print("****The %sth section has the parameters:****" %i)
            print(title)

            # Extract Reads for this section.
            Sections = root.xpath("/*/Section")

            reads = root.xpath("/*/*/*/Well")

            wellIDs = [read.attrib['Pos'] for read in reads]

            data = [(float(s.text), float(s.attrib['WL']), r.attrib['Pos'])
                    for r in reads
                    for s in r]

            dataframe = pd.DataFrame(data, columns=['fluorescence','wavelength (nm)','Well'])

            dataframe_pivot = pd.pivot_table(dataframe, index = 'wavelength (nm)', columns = ['Well'])

            # Make plot, complete with separate png for each section.
            section_name = sect.attrib['Name']

            fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('A' + str(i)).plot(ax=axes[0,0], title='A', c=cm.hsv(i*15))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('B' + str(i)).plot(ax=axes[0,1], title='B', c=cm.hsv(i*15))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('C' + str(i)).plot(ax=axes[0,2], title='C', c=cm.hsv(i*15))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('D' + str(i)).plot(ax=axes[1,0], title='D', c=cm.hsv(i*15))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('E' + str(i)).plot(ax=axes[1,1], title='E', c=cm.hsv(i*15))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('F' + str(i)).plot(ax=axes[1,2], title='F', c=cm.hsv(i*15))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('G' + str(i)).plot(ax=axes[2,0], title='G', c=cm.hsv(i*15))
            for i in range(1,12):
                dataframe_pivot.fluorescence.get('H' + str(i)).plot(ax=axes[2,1], title='H', c=cm.hsv(i*15))
            fig.suptitle('%s \n %s \n Barcode = %s' % (title, plate_type, barcode), fontsize=14)
            fig.subplots_adjust(hspace=0.3)
            plt.savefig('%s_%s.png' % (file_name, section_name))

    return

#############################################
 ###               singlet_384            ###
#############################################

def plot_singlet_one_section(data, section):

    #This assumes protein data is in row above buffer data

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 12))

    axes = axes.ravel()

    ALPHABET = string.ascii_uppercase

    well = dict()
    for j in string.ascii_uppercase:
        for i in range(1,25):
            well['%s' %j + '%s' %i] = i

    for i in range(0,15,2):
        protein_row = ALPHABET[i]
        buffer_row = ALPHABET[i+1]

        try:
            part1_data_protein = platereader.select_data(data, section, protein_row)
            part1_data_buffer = platereader.select_data(data, section, buffer_row)

            reorder_protein = reorder2list(part1_data_protein,well)
            reorder_buffer = reorder2list(part1_data_buffer,well)
            
        except:
            print('***no %s%s data***' %(protein_row,buffer_row) )
            continue    
        
        axes[i/2].set_color_cycle(['black','red'])

        if i/2 == 1:
            axes[i/2].plot(reorder_protein,marker='o',linestyle='None',label='protein+ligand')
            axes[i/2].plot(reorder_buffer,marker='o',linestyle='None',label='buffer+ligand')
            axes[i/2].legend(frameon=True)

        axes[i/2].plot(reorder_protein,marker='o',linestyle='None')
        axes[i/2].plot(reorder_buffer,marker='o',linestyle='None')
        axes[i/2].set_xticklabels(range(-4,25,5))
        axes[i/2].set_xlabel('Column Index', horizontalalignment='right',position=(1,1),fontsize=8)
        axes[i/2].set_ylabel('Fluorescence')
        axes[i/2].set_title('%s,%s' %(protein_row,buffer_row))

    fig.suptitle('%s' %section)

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

def process_files_singlet(xml_files):

    so_many = len(xml_files)
    print("****This script is about to make png files for %s xml files. ****"  % so_many)

    for file in xml_files:

        file_name = os.path.splitext(file)[0]

        data = platereader.read_icontrol_xml(file)

        print("****The xml file " + file + " has %s data sections:****" % len(data.keys()))
        print(data.keys())

        for key in data.keys():
            plot_singlet_one_section(data,key)

            #fig.tight_layout()
            #fig.suptitle("%s_%s" % (file_name,key), fontsize=18)

            plt.savefig('%s_%s.png' % (file_name,key))
            print('Look at how pretty your data is: %s_%s.png' % (file_name,key))
    return


def entry_point():
    xml_files = args.files
    if args.type == 'singlet_384':
        process_files_singlet(xml_files)
    if args.type == 'singlet_96':
        process_files_five(xml_files)
    if args.type == 'spectra':
        process_files_spectra(xml_files)
    if args.type == 'scan':
        process_files_scan(xml_files)

if __name__ == '__main__':
    xml_files = args.files
    if args.type == 'singlet_384':
        process_files_singlet(xml_files)
    if args.type == 'singlet_96':
        process_files_five(xml_files)
    if args.type == 'spectra':
        process_files_spectra(xml_files)
    if args.type == 'scan':
        process_files_scan(xml_files)
