# A somewhat ugly, utilitarian script takes xml data file output from the Tecan Infinite m1000 Pro 
# plate reader and allows for the quick visual inspection of raw data.
#
# Usage: python inputs.py
#        python quickmodel.py
#  or
#        python inputs.py
#        python quickmodel.py "/Users/hansons/Documents/github/fluorescence-assay-manuscript/data/singlet/DMSO-backfill/"

from assaytools import platereader
import matplotlib.pyplot as plt
import string
from glob import glob
import os
import string
import json
import numpy as np

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

def quick_model(inputs, path=None):

    if path!=None:
        xml_files = glob("%s/*.xml" % inputs['xml_file_path'])
    else:
        xml_files = glob("%s/*.xml" % path)

    #if args.path:
    #    xml_files = glob("%s/*.xml" % args.path)
    #else: 
    #    xml_files = glob("%s/*.xml" % inputs['xml_file_path'])
    
    for my_file in xml_files:
    
        file_name = os.path.splitext(my_file)[0]

        print file_name 
   
        data = platereader.read_icontrol_xml(my_file)

        well = dict()
        for j in string.ascii_uppercase:
            for i in range(1,25):
                well['%s' %j + '%s' %i] = i

        ALPHABET = string.ascii_uppercase
    
        for i in range(0,15,2):
            protein_row = ALPHABET[i]
            buffer_row = ALPHABET[i+1]
        
            name = "%s-%s%s"%(inputs['ligand_order'][i/2],protein_row,buffer_row)
            
            print name
                  
            metadata = {}
            metadata = dict(inputs)
            
            try:
                part1_data_protein = platereader.select_data(data, inputs['section'], protein_row)
                part1_data_buffer = platereader.select_data(data, inputs['section'], buffer_row)
            except:
                continue

            reorder_protein = reorder2list(part1_data_protein,well)
            reorder_buffer = reorder2list(part1_data_buffer,well)
    
            #these need to be changed so they are TAKEN FROM INPUTS!!!
        
            # Uncertainties in protein and ligand concentrations.
            dPstated = 0.35 * inputs['Pstated'] # protein concentration uncertainty
            dLstated = 0.08 * inputs['Lstated'] # ligand concentraiton uncertainty (due to gravimetric preparation and HP D300 dispensing)
            
            from assaytools import pymcmodels
            pymc_model = pymcmodels.make_model(inputs['Pstated'], dPstated, inputs['Lstated'], dLstated, 
               top_complex_fluorescence=reorder_protein,
               top_ligand_fluorescence=reorder_buffer,
               use_primary_inner_filter_correction=True, 
               use_secondary_inner_filter_correction=True, 
               assay_volume=inputs['assay_volume'], DG_prior='uniform')
                        
            import datetime
            my_datetime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
            
            mcmc = pymcmodels.run_mcmc(pymc_model, db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,my_datetime))
            
            map = pymcmodels.map_fit(pymc_model)
    
            pymcmodels.show_summary(pymc_model, map, mcmc)
    
            DeltaG = mcmc.DeltaG.trace().mean()
            dDeltaG = mcmc.DeltaG.trace().std()
    
            ## PLOT MODEL
            #from assaytools import plots
            #figure = plots.plot_measurements(Lstated, Pstated, pymc_model, mcmc=mcmc)
            #Code below inspired by import above, but didn't quite work to import it...
            plt.clf()
            plt.subplot(211)
            property_name = 'top_complex_fluorescence'
            complex = getattr(pymc_model, property_name)
            plt.semilogx(inputs['Lstated'], complex.value, 'ko',label='complex')
            property_name = 'top_ligand_fluorescence'
            ligand = getattr(pymc_model, property_name)
            plt.semilogx(inputs['Lstated'], ligand.value, 'ro',label='ligand')
            for top_complex_fluorescence_model in mcmc.top_complex_fluorescence_model.trace()[::10]:
                plt.semilogx(inputs['Lstated'], top_complex_fluorescence_model, 'k:')
            for top_ligand_fluorescence_model in mcmc.top_ligand_fluorescence_model.trace()[::10]:
                plt.semilogx(inputs['Lstated'], top_ligand_fluorescence_model, 'r:')
            plt.xlabel('$[L]_T$ (M)');
            plt.ylabel('fluorescence units');
            plt.legend(loc=0);
 
            ## PLOT TRACE
            plt.subplot(212)
            plt.hist(mcmc.DeltaG.trace(), 40, alpha=0.75, label="DeltaG = %.1f +- %.1f kT"%(DeltaG, dDeltaG))
            plt.axvline(x=DeltaG,color='blue')
            plt.legend(loc=0)
            plt.xlabel('$\Delta G$ ($k_B T$)');
            plt.ylabel('$P(\Delta G)$');

            plt.suptitle("%s: %s" % (name, my_datetime))
            plt.tight_layout()
            
            fig1 = plt.gcf()
            fig1.savefig('delG_%s-%s.png'%(name, my_datetime))
    
            np.save('DeltaG_%s-%s.npy'%(name, my_datetime),DeltaG)
            np.save('DeltaG_trace_%s-%s.npy'%(name, my_datetime),mcmc.DeltaG.trace())
            
            Kd = np.exp(mcmc.DeltaG.trace().mean())
            dKd = np.exp(mcmc.DeltaG.trace()).std()
            
            if (Kd < 1e-12):
                Kd_summary = "%.1f nM +- %.1f fM" % (Kd/1e-15, dKd/1e-15)
            elif (Kd < 1e-9):
                Kd_summary = "%.1f pM +- %.1f pM" % (Kd/1e-12, dKd/1e-12)
            elif (Kd < 1e-6):
                Kd_summary = "%.1f nM +- %.1f nM" % (Kd/1e-9, dKd/1e-9)
            elif (Kd < 1e-3):
                Kd_summary = "%.1f uM +- %.1f uM" % (Kd/1e-6, dKd/1e-6)
            elif (Kd < 1):
                Kd_summary = "%.1f mM +- %.1f mM" % (Kd/1e-3, dKd/1e-3)
            else:
                Kd_summary = "%.3e M +- %.3e M" % (Kd, dKd)
    
            outputs = {
                'raw_data_file'   : my_file,
                'name'            : name,
                'analysis'        : 'pymcmodels', #right now this is hardcoded, BOOO
                'outfiles'        : '%s_mcmc-%s.pickle, delG_%s-%s.png,DeltaG_%s-%s.npy,DeltaG_trace_%s-%s.npy'%(name,my_datetime,name,my_datetime,name,my_datetime,name,my_datetime),
                'DeltaG'          : "DeltaG = %.1f +- %.1f kT" % (DeltaG, dDeltaG),
                'Kd'              : Kd_summary,
                'datetime'        : my_datetime
            }
    
            metadata.update(outputs)
            
            metadata['Pstated'] = metadata['Pstated'].tolist()
            metadata['Lstated'] = metadata['Lstated'].tolist()

            with open('%s-%s.json'%(name,my_datetime), 'w') as outfile:
                json.dump(metadata, outfile, sort_keys = True, indent = 4, ensure_ascii=False)

def entry_point():

    import argparse

    # Define argparse stuff

    parser = argparse.ArgumentParser(description="""Analyze your fluorescence binding data by running assaytools on your xml files:
    > python inputs.py      
    > python quickmodel.py "/Users/hansons/Documents/github/fluorescence-assay-manuscript/data/singlet/DMSO-backfill/" """)
    parser.add_argument("path", nargs='*', help="path to xml file(s) to analyze",default=None)
    args = parser.parse_args()
    print args.path

    # Define inputs
    with open('inputs.json', 'r') as fp:
        inputs = json.load(fp)

    inputs['Lstated'] = np.asarray(inputs['Lstated'])
    inputs['Pstated'] = np.asarray(inputs['Pstated'])

    if args.path!=None:
        quick_model(inputs,path=args.path)
    else:
        quick_model(inputs)

if __name__ == '__main__':
    entry_point()



