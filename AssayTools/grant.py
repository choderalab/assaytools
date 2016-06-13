#!/usr/bin/env python

"""
Tools for assisting in reading, manipulating, and plotting data from plate readers for John's R01 due Feb 2016.

"""

#=============================================================================================
# Imports
#=============================================================================================

from lxml import etree
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

sns.set(style='white')
sns.set_context('talk')
sns.despine()

#=============================================================================================
# Tecan Infinite plate reader helper functions
#=============================================================================================

#THESE ARE ALL THE FUNCTIONS I MADE AND USED WHILE MAKING FIGURES FOR THE GRANT

#This function allows us to import xml format data files and convert them to a pandas dataframe
def xml2df(file):

    root = etree.parse(file)

    data = []

    reads = root.xpath("/*/Section[1]/*/Well")

    wellIDs = [read.attrib['Pos'] for read in reads]

    data = [(s.text, float(s.attrib['WL']), r.attrib['Pos'])
        for r in reads
        for s in r]

    dataframe = pd.DataFrame(data, columns=['fluorescence','wavelength (nm)','Well'])
            
    ### dataframe_rep replaces 'OVER' (when fluorescence signal maxes out) with '3289277', an arbitrarily high number

    dataframe_rep = dataframe.replace({'OVER':'3289277'})

    dataframe_rep[['fluorescence']] = dataframe_rep[['fluorescence']].astype('float')
            
    dataframe_pivot = pd.pivot_table(dataframe_rep, index = 'wavelength (nm)', columns = ['Well'])
    
    #Rearrange columns so they're in the right order
    cols =  dataframe_pivot['fluorescence'].columns.tolist()
    cols = [cols[0]] + cols[4:12] + cols[1:4] + [cols[12]] + cols[16:24] + cols[13:16]
    dataframe_reindex =  dataframe_pivot.reindex_axis(cols,level='Well',axis=1)
    
    return dataframe_reindex

#This function allows us to plot spectra
def plot_spectra_grid(file_set,ligands,*args,**kwargs):
    #Example Syntax:
    #   plot_spectra_grid(file_set,ligands,protein=['Src'],ligand=['Bosutinib'])
    #   plot_spectra_grid(file_set,ligands)
    #   plot_spectra_grid(file_set,ligands,output='Figure1')
    protein = kwargs.get('protein', None)
    ligand = kwargs.get('ligand', None)
    output = kwargs.get('output', None)
    
    if protein != None:
        these_proteins = protein
    else:
        these_proteins = file_set.keys()
    if ligand != None:
        these_ligands = ligand
    else:
        these_ligands = ligands
    print these_proteins
    print len(these_proteins)
    print these_ligands
    print len(these_ligands)
    
    grid = len(these_proteins) + len(these_ligands)
    if grid ==2:
        my_figsize = (12,8)
    else:
        my_figsize = (24,18)
    
    fig, axes = plt.subplots(nrows=len(these_ligands), ncols=len(these_proteins), figsize= my_figsize, sharey=True, sharex=True)
    
    for j,protein in enumerate(these_proteins):
        for k,ligand in enumerate(these_ligands):
            index = ligands.index(ligand)
            file = file_set[protein][index]
            
            if grid == 2:
                my_axes = None
            else:
                my_axes = axes[k,j]
    
            # pick a title
            title = "%s - %s" %(protein, ligand)
    
            # make a dataframe
            df = xml2df(file)
      
            # plot the spectra
            #fig = plt.figure();
            for i in range(11):
                df['fluorescence'].iloc[:,i].plot(ylim=(0,100000),xlim=(400,600),linewidth=3, ax = my_axes, c=cm.hsv(i*15), title=title);
                df['fluorescence'].iloc[:,11+i].plot(ylim=(0,100000),xlim=(400,600),legend=False, ax = my_axes, linewidth=4,c=cm.gray(i*15+50),fontsize =20, title=title);

            sns.despine()
            plt.yticks([])
            plt.tight_layout();
            
    if output != None:
        plt.savefig('%s.png'%output, dpi=1000)

#This function allows us to get wells from sections, parsing an xml file
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
                 for row in range(1,9)
                ]
                for col in range(1,13)
                ]
                
    return welllist

# This function makes a dataframe from our file once it's been parsed. 
def file2df(file,columns):
    root = etree.parse(file)
    #Just going to work with topread for now
    TopRead = root.xpath("/*/Section")[0]
    welllist = get_wells_from_section(TopRead)
    df = pd.DataFrame(welllist, columns = columns)
    return df

#This function allows us to plot the saturation curve at a single wavelength of the spectra
def plot_spectra2singlet(file_set,ligands,wavelength,output):
    
    fig, axes = plt.subplots(nrows=len(file_set), ncols=4, figsize=(22,22))
    
    proteins = file_set.keys()
    
    for j,protein in enumerate(file_set):
    
        files = file_set[protein]
        print file_set[protein]
    
        for i in range(len(files)):
        
            #Extract data from the xml file and make a dataframe
            df = xml2df(files[i])

            hardcode = wavelength #nm
     
            # This plots things.
            df.loc[hardcode][0:11].plot(ax = axes[j,i], xticks=[],linewidth=4)
            df.loc[hardcode][11:23].plot(ax = axes[j,i], xticks=[],linewidth=4,title ='%s - %s' %(proteins[j], ligands[i]))
        plt.text(4,15000,'wavelength %s nm'%hardcode,fontsize=20)
        
    plt.savefig(output,dpi=1000)
