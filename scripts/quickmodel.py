# A somewhat ugly, utilitarian script takes xml data file output from the Tecan Infinite m1000 Pro
# plate reader and allows for the quick visual inspection of raw data.
#
#    Usage:
#        python quickmodel.py --inputs 'inputs_example'
#        or
#        python quickmodel.py --inputs 'inputs_spectra_example' --type 'spectra'


from assaytools import parser
import matplotlib.pyplot as plt
import string
from glob import glob
import os
import json
import numpy as np
import traceback

import seaborn as sns
import pymbar

def quick_model(inputs, nsamples=1000, nthin=20):
    """
    Quick model for both spectra and single wavelength experiments

    Parameters
    ----------
    inputs : dict
        Dictionary of input information
    nsamples : int, optional, default=1000
        Number of MCMC samples to collect
    nthin : int, optional, default=20
        Thinning interval ; number of MCMC steps per sample collected
    """


    [complex_fluorescence, ligand_fluorescence] = parser.get_data_using_inputs(inputs)  

    for name in complex_fluorescence.keys():

            print(name)

            metadata = {}
            metadata = dict(inputs)

            #these need to be changed so they are TAKEN FROM INPUTS!!!

            # Uncertainties in protein and ligand concentrations.
            dPstated = 0.35 * inputs['Pstated'] # protein concentration uncertainty
            dLstated = 0.08 * inputs['Lstated'] # ligand concentraiton uncertainty (due to gravimetric preparation and HP D300 dispensing)

            from assaytools import pymcmodels
            pymc_model = pymcmodels.make_model(inputs['Pstated'], dPstated, inputs['Lstated'], dLstated,
               top_complex_fluorescence=complex_fluorescence[name],
               top_ligand_fluorescence=ligand_fluorescence[name],
               use_primary_inner_filter_correction=True,
               use_secondary_inner_filter_correction=True,
               assay_volume=inputs['assay_volume'], DG_prior='uniform')

            import datetime
            my_datetime = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
            my_datetime_filename = datetime.datetime.now().strftime("%Y-%m-%d %H%M")

            nburn = 0 # no longer need burn-in since we do automated equilibration detection
            niter = nthin*nsamples # former total simulation time
            # Note that nsamples and niter are not the same: Here nsamples is 
            # multiplied by nthin (default 20) so that your final mcmc samples will be the same 
            # as your nsamples, but the actual niter will be 20x that!
            mcmc = pymcmodels.run_mcmc(pymc_model,
                nthin=nthin, nburn=nburn, niter=niter,
                db = 'pickle', dbname = '%s_mcmc-%s.pickle'%(name,my_datetime))

            map = pymcmodels.map_fit(pymc_model)

            DeltaG_map = map.DeltaG.value
            DeltaG = mcmc.DeltaG.trace().mean()
            dDeltaG = mcmc.DeltaG.trace().std()

            ## DEFINE EQUILIBRATION
            #Calculate a mean and std from DeltaG trace after equil

            [t,g,Neff_max] = pymbar.timeseries.detectEquilibration(mcmc.DeltaG.trace())
            DeltaG_equil = mcmc.DeltaG.trace()[t:].mean()
            dDeltaG_equil = mcmc.DeltaG.trace()[t:].std()

            #This is so plotting works on the cluster
            plt.switch_backend('agg')

            ## PLOT MODEL
            #from assaytools import plots
            #figure = plots.plot_measurements(Lstated, Pstated, pymc_model, mcmc=mcmc)
            #Code below inspired by import above, but didn't quite work to import it...

            plt.clf()
            plt.figure(figsize=(8,8))

            plt.subplot(311)
            property_name = 'top_complex_fluorescence'
            complex = getattr(pymc_model, property_name)
            property_name = 'top_ligand_fluorescence'
            ligand = getattr(pymc_model, property_name)
            for top_complex_fluorescence_model in mcmc.top_complex_fluorescence_model.trace()[::10]:
                plt.semilogx(inputs['Lstated'], top_complex_fluorescence_model, marker='.',color='silver')
            for top_ligand_fluorescence_model in mcmc.top_ligand_fluorescence_model.trace()[::10]:
                plt.semilogx(inputs['Lstated'], top_ligand_fluorescence_model, marker='.',color='lightcoral', alpha=0.2)
            plt.semilogx(inputs['Lstated'], complex.value, 'ko',label='complex')
            plt.semilogx(inputs['Lstated'], ligand.value, marker='o',color='firebrick',linestyle='None',label='ligand')

            plt.xlabel('$[L]_T$ (M)');
            plt.ylabel('fluorescence units');
            plt.legend(loc=0);

            ## PLOT HISTOGRAM
            import matplotlib.patches as mpatches
            import matplotlib.lines as mlines

            interval = np.percentile(a=mcmc.DeltaG.trace()[t:], q=[2.5, 50.0, 97.5])
            [hist,bin_edges] = np.histogram(mcmc.DeltaG.trace()[t:],bins=40,normed=True)
            binwidth = np.abs(bin_edges[0]-bin_edges[1])

            #Print summary
            print( 'Delta G (95% credibility interval after equilibration):')
            print( '   %.3g [%.3g,%.3g] k_B T' %(interval[1],interval[0],interval[2]))
            print( 'Delta G (mean and std after equil):')
            print('   %.3g +- %.3g k_B T' %(DeltaG_equil,dDeltaG_equil) )

            #set colors for 95% interval
            clrs = [(0.7372549019607844, 0.5098039215686274, 0.7411764705882353) for xx in bin_edges]
            idxs = bin_edges.argsort()
            idxs = idxs[::-1]
            gray_before = idxs[bin_edges[idxs] < interval[0]]
            gray_after = idxs[bin_edges[idxs] > interval[2]]
            for idx in gray_before:
                clrs[idx] = (.5,.5,.5)
            for idx in gray_after:
                clrs[idx] = (.5,.5,.5)

            plt.subplot(312)

            plt.bar(bin_edges[:-1],hist,binwidth,color=clrs, edgecolor = "white");
            sns.kdeplot(mcmc.DeltaG.trace()[t:],bw=.4,color=(0.39215686274509803, 0.7098039215686275, 0.803921568627451),shade=False)
            plt.axvline(x=interval[0],color=(0.5,0.5,0.5),linestyle='--')
            plt.axvline(x=interval[1],color=(0.5,0.5,0.5),linestyle='--')
            plt.axvline(x=interval[2],color=(0.5,0.5,0.5),linestyle='--')
            plt.axvline(x=DeltaG_map,color='black')
            plt.xlabel('$\Delta G$ ($k_B T$)',fontsize=16);
            plt.ylabel('$P(\Delta G)$',fontsize=16);
            plt.xlim(-20,-8)
            hist_legend = mpatches.Patch(color=(0.7372549019607844, 0.5098039215686274, 0.7411764705882353),
                        label = '$\Delta G$ =  %.3g [%.3g,%.3g] $k_B T$'
                        %(interval[1],interval[0],interval[2]) )
            map_legend = mlines.Line2D([],[],color='black',label="MAP = %.1f $k_B T$"%DeltaG_map)
            plt.legend(handles=[hist_legend,map_legend],fontsize=14,loc=0,frameon=True)

            ## PLOT TRACE
            plt.subplot(313)
            plt.plot(range(0,t),mcmc.DeltaG.trace()[:t], 'go',label='equil. at %s'%t);
            plt.plot(range(t,len(mcmc.DeltaG.trace())),mcmc.DeltaG.trace()[t:], 'o');
            plt.xlabel('MCMC sample');
            plt.ylabel('$\Delta G$ ($k_B T$)');
            plt.legend(loc=2);

            plt.suptitle("%s: %s" % (name, my_datetime))
            plt.tight_layout()

            fig1 = plt.gcf()
            fig1.savefig('delG_%s-%s.png'%(name, my_datetime_filename))

            plt.close('all') # close all figures

            Kd = np.exp(DeltaG_equil)
            dKd = np.exp(mcmc.DeltaG.trace()[t:]).std()
            Kd_interval = np.exp(interval)

            if (Kd < 1e-12):
                Kd_summary_interval = '%.1f [%.1f,%.1f] fM' %(Kd_interval[1]/1e-15,Kd_interval[0]/1e-15,Kd_interval[2]/1e-15)
                Kd_summary = "%.1f fM +- %.1f fM" % (Kd/1e-15, dKd/1e-15)
            elif (Kd < 1e-9):
                Kd_summary_interval = '%.1f [%.1f,%.1f] pM' %(Kd_interval[1]/1e-12,Kd_interval[0]/1e-12,Kd_interval[2]/1e-12)
                Kd_summary = "%.1f pM +- %.1f pM" % (Kd/1e-12, dKd/1e-12)
            elif (Kd < 1e-6):
                Kd_summary_interval = '%.1f [%.1f,%.1f] nM' %(Kd_interval[1]/1e-9,Kd_interval[0]/1e-9,Kd_interval[2]/1e-9)
                Kd_summary = "%.1f nM +- %.1f nM" % (Kd/1e-9, dKd/1e-9)
            elif (Kd < 1e-3):
                Kd_summary_interval = '%.1f [%.1f,%.1f] uM' %(Kd_interval[1]/1e-6,Kd_interval[0]/1e-6,Kd_interval[2]/1e-6)
                Kd_summary = "%.1f uM +- %.1f uM" % (Kd/1e-6, dKd/1e-6)
            elif (Kd < 1):
                Kd_summary_interval = '%.1f [%.1f,%.1f] mM' %(Kd_interval[1]/1e-3,Kd_interval[0]/1e-3,Kd_interval[2]/1e-3)
                Kd_summary = "%.1f mM +- %.1f mM" % (Kd/1e-3, dKd/1e-3)
            else:
                Kd_summary_interval = '%.3e [%.3e,%.3e] M' %(Kd_interval[1],Kd_interval[0],Kd_interval[2])
                Kd_summary = "%.3e M +- %.3e M" % (Kd, dKd)

            print('Kd (95% credibility interval after equilibration):')
            print('   %s' %Kd_summary_interval)
            print('Kd (mean and std after equil):')
            print('   %s' %Kd_summary)

            outputs = {
                #'raw_data_file'   : my_file,
                'complex_fluorescence' : complex_fluorescence[name],
                'ligand_fluorescence'  : ligand_fluorescence[name],
                't_equil'              : t,
                'name'                 : name,
                'analysis'             : 'pymcmodels', #right now this is hardcoded, BOOO
                'outfiles'             : '%s_mcmc-%s.pickle, delG_%s-%s.pdf,DeltaG_%s-%s.npy,DeltaG_trace_%s-%s.npy'%(name,my_datetime,name,my_datetime,name,my_datetime,name,my_datetime),
                'DeltaG_cred_int'      : '$\Delta G$ =  %.3g [%.3g,%.3g] $k_B T$'  %(interval[1],interval[0],interval[2]),
                'DeltaG'               : "DeltaG = %.1f +- %.1f kT, MAP estimate = %.1f" % (DeltaG, dDeltaG, DeltaG_map),
                'Kd'                   : Kd_summary,
                'bin_edges'            : bin_edges,
                'hist'                 : hist,
                'binwidth'             : binwidth,
                'clrs'                 : clrs,
                'datetime'             : my_datetime
            }

            metadata.update(outputs)

            metadata['ligand_fluorescence'] = metadata['ligand_fluorescence'].tolist()
            metadata['complex_fluorescence'] = metadata['complex_fluorescence'].tolist()
            metadata['t_equil'] = metadata['t_equil'].tolist()
            metadata['bin_edges'] = metadata['bin_edges'].tolist()
            metadata['hist'] = metadata['hist'].tolist()
            metadata['Pstated'] = metadata['Pstated'].tolist()
            metadata['Lstated'] = metadata['Lstated'].tolist()

            import json
            with open('%s-%s.json'%(name,my_datetime), 'w') as outfile:
                json.dump(metadata, outfile, sort_keys = True, indent = 4, ensure_ascii=False)

def entry_point():

    import sys

    #This allows us to import local inputs.py
    sys.path.append('./')

    import argparse

    # Define argparse stuff

    parser = argparse.ArgumentParser(description="""Analyze your fluorescence binding data by running assaytools on your xml files:
    > python quickmodel.py --inputs 'inputs_example' """)
    parser.add_argument("--inputs", help="inputs file (python script, .py not needed)",default=None)
    parser.add_argument("--type", help="type of data (spectra, singlet)",choices=['spectra','singlet'],default='singlet')
    parser.add_argument("--nsamples", help="number of samples",default=1000, type=int)
    parser.add_argument("--nthin", help="thinning interval",default=20, type=int)
    args = parser.parse_args()

    # Define inputs
    # If an inputs###.py is defined, it is run and used. An inputs.json is also created.
    if args.inputs!=None:
        inputs_file = args.inputs
        name = 'inputs'
        #this is the equivalent of from inputs_file import name as inputs
        inputs = getattr(__import__(inputs_file, fromlist=[name]), name)
    if args.inputs==None:
        with open('inputs.json', 'r') as fp:
            inputs = json.load(fp)

    inputs['Lstated'] = np.asarray(inputs['Lstated'])
    inputs['Pstated'] = np.asarray(inputs['Pstated'])

    if args.type == 'singlet':
        quick_model(inputs, nsamples=args.nsamples, nthin=args.nthin)
    if args.type == 'spectra':
        quick_model(inputs, nsamples=args.nsamples, nthin=args.nthin)

if __name__ == '__main__':
    entry_point()
