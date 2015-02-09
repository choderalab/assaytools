#!/usr/bin/env python

"""
Tools for assisting in plotting data.

"""

#=============================================================================================
# Imports
#=============================================================================================

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

#=============================================================================================
# Plot fluorescence.
#=============================================================================================

def plot_measurements(Lstated, Pstated, pymc_model, map=None, mcmc=None, figsize=(10,5), subsample=10):
    """
    Plot pbserved fluorescence and absorbance.

    Parameters
    ----------
    Lstated : np.array
       Stated ligand concentrations (in M) corresponding to fluorescence measurements.
    Pstated : np.array
       Stated protein concentrations (in M) corresponding to fluorescence measurements.
    pymc_model : pymc model
       pymc model (created with pymcmodels.make_model) containing data
    map : pymc.MAP, optional, default=None
       If specified, will plot the MAP fit.
    mcmc : pymc.MCMC, optional, default=None
       If specified, will plot MCMC samples.
    figsize : (x,y) tuple, optional, default=(10,5)
       Figure size (in inches).
    subsample : int, optional, default=10
       For mcmc, degree to subsample trace data for plotting.

    Returns
    -------
    figure : matplotlib Figure
       The matplotlib-generated figure

    """

    # Determine how many plots we will generate
    nplots = 0
    if hasattr(pymc_model, 'top_complex_fluorescence') or hasattr(pymc_model, 'top_ligand_fluorescence'):
        nplots += 1
    if hasattr(pymc_model, 'bottom_complex_fluorescence') or hasattr(pymc_model, 'bottom_ligand_fluorescence'):
        nplots += 1
    if hasattr(pymc_model, 'ligand_ex_absorbance') or hasattr(pymc_model, 'ligand_em_absorbance'):
        nplots += 1

    figsize = (figsize[0], nplots*figsize[1])

    figure = plt.figure(figsize=figsize)
    plt.clf();
    plt.hold(True)

    # Plot sections.
    plot_index = 0
    if hasattr(pymc_model, 'top_complex_fluorescence') or hasattr(pymc_model, 'top_ligand_fluorescence'):
        subplot = plt.subplot(nplots, 1, plot_index+1)
        legend = list()

        property_name = 'top_complex_fluorescence'
        if hasattr(pymc_model, property_name):
            property = getattr(pymc_model, property_name)
            plt.semilogx(Lstated, property.value, 'ko');
            legend.append('complex');

        property_name = 'top_ligand_fluorescence'
        if hasattr(pymc_model, property_name):
            property = getattr(pymc_model, property_name)
            plt.semilogx(Lstated, property.value, 'ro');
            legend.append('ligand');

        if map and hasattr(pymc_model, 'top_complex_fluorescence'):
            plt.semilogx(Lstated, map.top_complex_fluorescence_model.value, 'k-')

        if map and hasattr(pymc_model, 'top_ligand_fluorescence'):
            plt.semilogx(Lstated, map.top_ligand_fluorescence_model.value, 'r-')

        if mcmc and hasattr(pymc_model, 'top_complex_fluorescence'):
            for top_complex_fluorescence_model in mcmc.top_complex_fluorescence_model.trace()[::subsample]:
                plt.semilogx(Lstated, top_complex_fluorescence_model, 'k:')

        if mcmc and hasattr(pymc_model, 'top_ligand_fluorescence'):
            for top_ligand_fluorescence_model in mcmc.top_ligand_fluorescence_model.trace()[::subsample]:
                plt.semilogx(Lstated, top_ligand_fluorescence_model, 'r:')

        plt.xlabel('$[L]_T$ (M)');
        plt.ylabel('fluorescence units');
        plt.title('top fluorescence');
        plt.legend(legend, loc='best');

        plot_index += 1

    if hasattr(pymc_model, 'bottom_complex_fluorescence') or hasattr(pymc_model, 'bottom_ligand_fluorescence'):
        subplot = plt.subplot(nplots, 1, plot_index+1)
        legend = list()

        property_name = 'bottom_complex_fluorescence'
        if hasattr(pymc_model, property_name):
            property = getattr(pymc_model, property_name)
            plt.semilogx(Lstated, property.value, 'ko');
            legend.append('complex');

        property_name = 'bottom_ligand_fluorescence'
        if hasattr(pymc_model, property_name):
            property = getattr(pymc_model, property_name)
            plt.semilogx(Lstated, property.value, 'ro');
            legend.append('ligand');

        if map and hasattr(pymc_model, 'bottom_complex_fluorescence'):
            plt.semilogx(Lstated, map.bottom_complex_fluorescence_model.value, 'k-')

        if map and hasattr(pymc_model, 'bottom_ligand_fluorescence'):
            plt.semilogx(Lstated, map.bottom_ligand_fluorescence_model.value, 'r-')

        if mcmc and hasattr(pymc_model, 'bottom_complex_fluorescence'):
            for bottom_complex_fluorescence_model in mcmc.bottom_complex_fluorescence_model.trace()[::subsample]:
                plt.semilogx(Lstated, bottom_complex_fluorescence_model, 'k:')

        if mcmc and hasattr(pymc_model, 'bottom_ligand_fluorescence'):
            for bottom_ligand_fluorescence_model in mcmc.bottom_ligand_fluorescence_model.trace()[::subsample]:
                plt.semilogx(Lstated, bottom_ligand_fluorescence_model, 'r:')

        plt.xlabel('$[L]_T$ (M)');
        plt.ylabel('fluorescence units');
        plt.title('bottom fluorescence');
        plt.legend(legend, loc='best');

        plot_index += 1

    if hasattr(pymc_model, 'ligand_ex_absorbance') or hasattr(pymc_model, 'ligand_em_absorbance'):
        subplot = plt.subplot(nplots, 1, plot_index+1)
        legend = list()

        property_name = 'ligand_ex_absorbance'
        if hasattr(pymc_model, property_name):
            property = getattr(pymc_model, property_name)
            plt.semilogx(Lstated, property.value, 'ko');
            legend.append('excitation');

        property_name = 'ligand_em_absorbance'
        if hasattr(pymc_model, property_name):
            property = getattr(pymc_model, property_name)
            plt.semilogx(Lstated, property.value, 'ro');
            legend.append('emission');

        if map and hasattr(pymc_model, 'ligand_ex_absorbance'):
            plt.semilogx(Lstated, map.ligand_ex_absorbance_model.value, 'k-')

        if map and hasattr(pymc_model, 'ligand_em_absorbance'):
            plt.semilogx(Lstated, map.ligand_em_absorbance_model.value, 'r-')

        if mcmc and hasattr(pymc_model, 'ligand_ex_absorbance'):
            for ligand_ex_absorbance_model in mcmc.ligand_ex_absorbance_model.trace()[::subsample]:
                plt.semilogx(Lstated, ligand_ex_absorbance_model, 'k:')

        if mcmc and hasattr(pymc_model, 'ligand_em_absorbance'):
            for ligand_em_absorbance_model in mcmc.ligand_em_absorbance_model.trace()[::subsample]:
                plt.semilogx(Lstated, ligand_em_absorbance_model, 'r:')

        plt.xlabel('$[L]_T$ (M)');
        plt.ylabel('absorbance');
        plt.title('ligand absorbance');
        plt.legend(legend, loc='best');
        oldaxis = plt.axis();
        plt.axis([oldaxis[0], oldaxis[1], 0, oldaxis[3]]);

        plot_index += 1

    return figure

def plot_mcmc_results(Lstated, Pstated, path_length, mcmc, subsample=10, figsize=(15,3)):
    """
    Plot the results of MCMC sampling.

    Parameters
    ----------
    Lstated : np.array
       Stated ligand concentrations (in M) corresponding to fluorescence measurements.
    Pstated : np.array
       Stated protein concentrations (in M) corresponding to fluorescence measurements.
    path_length : float
       The liquid height in the plate (in cm).
    mcmc : pymc.MCMC object
       The results of an MCMC simulation.
    figsize : (x,y) in inches, optional, default=(15,3)
       The figure size for each subplot.

    """

    from assaytools import pymcmodels

    # Determine how many samples we ended up with
    nsamples = mcmc.DeltaG.trace().size

    plt.rcParams['figure.figsize'] = figsize

    # Plot trace of DeltaG.
    figure = plt.figure(figsize=figsize);
    from assaytools.pymcmodels import DG_min, DG_max
    plt.plot(mcmc.DeltaG.trace(), 'o');
    plt.xlabel('MCMC sample');
    plt.ylabel('$\Delta G$ ($k_B T$)');
    plt.axis([0, nsamples, DG_min, DG_max]);

    # Plot histogram of DeltaG.
    figure = plt.figure(figsize=figsize);
    plt.hold(True)
    nbins = 40
    [N,x,patches] = plt.hist(mcmc.DeltaG.trace(), nbins);
    plt.hist(mcmc.DeltaG.trace(), nbins);
    plt.xlabel('$\Delta G$ ($k_B T$)');
    plt.ylabel('$P(\Delta G)$');
    plt.axis([DG_min, DG_max, 0, N.max()]);
    plt.title('histogram of estimates for binding free energy');
    DGmean = mcmc.DeltaG.trace().mean()
    plt.plot([DGmean, DGmean], [0, N.max()], 'r-');

    # Plot trace of true protein concentration.
    figure = plt.figure(figsize=figsize)
    plt.plot(mcmc.Ptrue.trace()*1e6, 'o');
    plt.xlabel('MCMC sample');
    plt.ylabel('$[P]$ ($\mu$M)');
    plt.title('Trace of true protein concentrations');

    # Plot trace of true protein concentration.
    figure = plt.figure(figsize=figsize)
    X = mcmc.Ptrue.trace() * 1e6;
    plt.hist(X.reshape(X.size,1), 40);
    plt.xlabel('$[P]$ ($\mu$M)');
    plt.ylabel('frequency');
    plt.title('Histogram of true protein concentrations');

    # Plot trace of true ligand concentration.
    figure = plt.figure(figsize=figsize)
    plt.semilogy(mcmc.Ltrue.trace()*1e6, 'o');
    plt.xlabel('MCMC sample');
    plt.ylabel('$[L]$ ($\mu$M)');
    plt.title('Trace of true ligand concentrations in complex experiment');

    # Plot histograms of true ligand concentration.
    figure = plt.figure(figsize=figsize)
    X = np.log10(mcmc.Ltrue.trace())
    plt.clf();
    plt.hold(True);
    for i in range(X.shape[1]):
        [N,x,patches] = plt.hist(X[:,i], histtype='stepfilled');
        plt.hist(X[:,i], histtype='stepfilled');
        plt.plot(np.log10(Lstated[i])*np.array([1,1]), [0, N.max()], 'k-', linewidth=3);
    plt.hold(False);
    plt.xlabel('$\log_{10} ([L]/M)$');
    plt.ylabel('frequency');
    plt.title('Ligand true concentrations for complex experiment');

    # Plot histograms of true ligand concentration.
    figure = plt.figure(figsize=figsize)
    X = np.log10(mcmc.Ltrue_control.trace())
    plt.clf()
    plt.hold(True)
    for i in range(X.shape[1]):
        [N,x,patches] = plt.hist(X[:,i], histtype='stepfilled');
        plt.hist(X[:,i], histtype='stepfilled');
        plt.plot(np.log10(Lstated[i])*np.array([1,1]), [0, N.max()], 'k-', linewidth=3);
    plt.hold(False)
    plt.xlabel('$\log_{10} ([L]/M)$');
    plt.ylabel('frequency');
    plt.title('Ligand true concentrations for ligand-only experiment');

    # Plot trace of intrinsic fluorescence parameters.
    figure = plt.figure(figsize=figsize)
    plt.semilogy(mcmc.F_PL.trace(), 'o', mcmc.F_L.trace(), 'o', mcmc.F_P.trace(), 'o', mcmc.F_buffer.trace(), 'o', mcmc.F_plate.trace(), 'o');
    plt.legend(['complex fluorescence', 'ligand fluorescence', 'protein fluorescence', 'buffer fluorescence', 'plate fluorescence']);
    plt.xlabel('MCMC sample');
    plt.ylabel('relative fluorescence intensity');
    plt.title('trace of relative fluorescence intensities of various species');

    # Plot trace of intrinsic fluorescence parameters.
    figure = plt.figure(figsize=figsize)
    plt.hold(True);
    nbins = 40;
    print mcmc.F_PL.trace().shape
    plt.hist(np.log10(mcmc.F_PL.trace()), nbins, histtype='stepfilled');
    plt.hist(np.log10(mcmc.F_L.trace()), nbins, histtype='stepfilled');
    plt.hist(np.log10(mcmc.F_P.trace()), nbins, histtype='stepfilled');
    plt.hist(np.log10(mcmc.F_buffer.trace()), nbins, histtype='stepfilled');
    plt.hist(np.log10(mcmc.F_plate.trace()), nbins, histtype='stepfilled');
    plt.hold(False);
    plt.legend(['complex fluorescence', 'ligand fluorescence', 'protein fluorescence', 'buffer fluorescence', 'plate fluorescence']);
    plt.xlabel('relative fluorescence intensity $log_{10} F$ ');
    plt.ylabel('frequency');
    plt.title('histogram of relative fluorescence intensities of various species');

    # Plot trace of measurement error.
    figure = plt.figure(figsize=figsize)
    plt.hold(True)
    plt.semilogy(np.exp(mcmc.sigma_top.trace()), 'ro');
    plt.semilogy(np.exp(mcmc.sigma_bottom.trace()), 'kv');
    plt.xlabel('MCMC sample');
    plt.ylabel('fluorescence measurement error');
    plt.legend(['top', 'bottom']);

    # Plot trace of absorbance measurement error.
    figure = plt.figure(figsize=figsize)
    plt.hold(True)
    plt.plot(np.exp(mcmc.sigma_abs.trace()), 'ro');
    plt.xlabel('MCMC sample');
    plt.ylabel('absorbance measurement error');

    # Plot extinction coefficient trace
    figure = plt.figure(figsize=figsize)
    plt.hold(True);
    plt.plot(mcmc.epsilon_ex.trace() / 1e3, 'ko');
    plt.plot(mcmc.epsilon_em.trace() / 1e3, 'ro');
    plt.legend(['$\epsilon_{ex}$', '$\epsilon_{em}$']);
    plt.xlabel('MCMC sample');
    plt.ylabel('ligand extinction coefficient $\epsilon / 10^3$ (1/M/cm)');
    plt.title('trace of ligand extinction coefficient $\epsilon$');

    # Plot extinction coefficient histogram
    figure = plt.figure(figsize=figsize)
    plt.hold(True);
    plt.hist(mcmc.epsilon_ex.trace() / 1e3, nbins, histtype='stepfilled');
    plt.hist(mcmc.epsilon_em.trace() / 1e3, nbins, histtype='stepfilled');
    plt.xlabel('ligand extinction coefficient $\epsilon / 10^3$ (1/M/cm)');
    plt.ylabel('frequency');
    plt.title('histogram of ligand extinction coefficient $\epsilon$');
    plt.legend(['$\epsilon_{ex}$', '$\epsilon_{em}$']);

    # Plot top inner filter effect
    figure = plt.figure(figsize=figsize)
    plt.clf();
    plt.hold(True)
    for (Ltrue, epsilon_ex, epsilon_em) in zip(mcmc.Ltrue.trace()[::subsample], mcmc.epsilon_ex.trace()[::subsample], mcmc.epsilon_em.trace()[::subsample]):
        primary_scaling = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, 0.0, path_length, Ltrue, geometry='top')
        secondary_scaling = pymcmodels.inner_filter_effect_attenuation(0.0, epsilon_em, path_length, Ltrue, geometry='top')
        scaling = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue, geometry='top')
        plt.semilogx(Lstated, primary_scaling, 'k-');
        plt.semilogx(Lstated, secondary_scaling, 'r-');
    plt.hold(False)
    plt.xlabel('stated $[L]_s$ (M)');
    plt.ylabel('inner filter effect');
    plt.title('top inner filter effect attenuation');
    plt.axis([Lstated.min(), Lstated.max(), 0, 1]);
    plt.legend(['primary', 'secondary'], loc='best');

    # Plot top inner filter effect
    figure = plt.figure(figsize=figsize)
    plt.clf();
    plt.hold(True)
    for (Ltrue, epsilon_ex, epsilon_em) in zip(mcmc.Ltrue.trace()[::subsample], mcmc.epsilon_ex.trace()[::subsample], mcmc.epsilon_em.trace()[::subsample]):
        primary_scaling = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, 0.0, path_length, Ltrue, geometry='bottom')
        secondary_scaling = pymcmodels.inner_filter_effect_attenuation(0.0, epsilon_em, path_length, Ltrue, geometry='bottom')
        scaling = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue, geometry='bottom')
        plt.semilogx(Lstated, primary_scaling, 'k-');
        plt.semilogx(Lstated, secondary_scaling, 'r-');
    plt.hold(False)
    plt.xlabel('stated $[L]_s$ (M)');
    plt.ylabel('inner filter effect');
    plt.title('bottom inner filter effect attenuation');
    plt.axis([Lstated.min(), Lstated.max(), 0, 1]);
    plt.legend(['primary', 'secondary'], loc='best');

