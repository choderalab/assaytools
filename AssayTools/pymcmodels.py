#!/usr/bin/env python

"""
A test of pymc for ITC.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy as np
import pymc
import pint

#=============================================================================================
# Physical constants
#=============================================================================================

Na = 6.02214179e23 # Avogadro's number (number/mol)
kB = Na * 1.3806504e-23 / 4184.0 # Boltzmann constant (kcal/mol/K)
C0 = 1.0 # standard concentration (M)

#=============================================================================================
# Parameters for MCMC sampling
#=============================================================================================

DG_min = np.log(1e-15) # kT, most favorable (negative) binding free energy possible; 1 fM
DG_max = +0 # kT, least favorable binding free energy possible
niter = 500000 # number of iterations
nburn = 50000 # number of burn-in iterations to discard
nthin = 500 # thinning interval

#=============================================================================================
# PyMC models
#=============================================================================================

# Create a pymc model
def make_model(Pstated, dPstated, Lstated, dLstated, Fobs_i, Fligand_i,
               DG_prior='uniform',
               concentration_priors='lognormal',
               use_primary_inner_filter_correction=True, assay_volume=100e-6, well_area=0.1586,
               epsilon=None, depsilon=None):
    """
    Build a PyMC model for an assay that consists of N wells of protein:ligand at various concentrations and an additional N wells of ligand in buffer, with the ligand at the same concentrations as the corresponding protein:ligand wells.

    Parameters
    ----------
    Pstated : numpy.array of N values
       Stated protein concentrations for all protein:ligand wells of assay. Units of molarity.
    dPstated : numpy.array of N values
       Absolute uncertainty in stated protein concentrations for all wells of assay. Units of molarity.
       Uncertainties currently cannot be zero.
    Lstated : numpy.array of N values
       Stated ligand concentrations for all protein:ligand and ligand wells of assay, which must be the same with and without protein. Units of molarity.
    dLstated : numpy.array of N values
       Absolute uncertainty in stated protein concentrations for all wells of assay. Units of molarity.
       Uncertainties currently cannot be zero
    DG_prior : str, optional, default='uniform'
       Prior to use for reduced free energy of binding (DG): 'uniform' (uniform over reasonable range), or 'chembl' (ChEMBL-inspired distribution); default: 'uniform'
    concentration_priors : str, optional, default='lognormal'
       Prior to use for protein and ligand concentrations. Available options are ['lognormal', 'normal'].
    use_primary_inner_filter_correction : bool, optional, default=True
       If true, will infer ligand extinction coefficient epsilon and apply primary inner filter correction to attenuate excitation light.
    assay_volume : float, optional, default=100e-6
       Assay volume. Units of L. Default 100 uL.
    well_area : float, optional, default=0.1586
       Well area. Units of cm^2. Default 0.1586 cm^2, for half-area plate.
    epsilon : float, optional, default=None
       Orthogonal measurement of ligand extinction coefficient at excitation wavelength. If None, will use a uniform prior.
    depsilon : float, optional, default=None
       Uncertainty (standard error) in measurement 'epsilon'.

    Returns
    -------
    pymc_model : dict
       A dict mapping variable names to onbjects that can be used as a PyMC model object.

    Examples
    --------
    Create a simple model

    >>> N = 12 # 12 wells per series of protein:ligand or ligand alone
    >>> Pstated = np.ones([N], np.float64) * 1e-6
    >>> Lstated = 20.0e-6 / np.array([10**(float(i)/2.0) for i in range(N)])
    >>> dPstated = 0.10 * Pstated
    >>> dLstated = 0.08 * Lstated
    >>> Fobs_i = array([ 689., 683., 664., 588., 207., 80., 28., 17., 10., 11., 10., 10.])
    >>> Fligand_i = array([ 174., 115., 57., 20., 7., 6., 6., 6., 6., 7., 6., 7.])
    >>> pymc_model = pymcmodels.make_model(Pstated, dPstated, Lstated, dLstated, Fobs_i, Fligand_i)

    """

    # Compute path length.
    path_length = assay_volume * 1000 / well_area # cm, needed for inner filter effect corrections

    # Compute number of samples.
    N = len(Lstated)

    # Check input.
    if (len(Pstated) != N):
        raise Exception('len(Pstated) [%d] must equal len(Lstated) [%d].' % (len(Pstated), len(Lstated)))
    if (len(dPstated) != N):
        raise Exception('len(dPstated) [%d] must equal len(Lstated) [%d].' % (len(dPstated), len(Lstated)))
    if (len(dLstated) != N):
        raise Exception('len(dLstated) [%d] must equal len(Lstated) [%d].' % (len(dLstated), len(Lstated)))

    # Create an empty dict to hold the model.
    pymc_model = dict()

    # Prior on binding free energies.
    if DG_prior == 'uniform':
        DeltaG = pymc.Uniform('DeltaG', lower=DG_min, upper=DG_max, value=-20) # binding free energy (kT), uniform over huge range
    elif DG_prior == 'chembl':
        DeltaG = pymc.Normal('DeltaG', mu=0, tau=1./(12.5**2)) # binding free energy (kT), using a Gaussian prior inspured by ChEMBL
    else:
        raise Exception("DG_prior = '%s' unknown. Must be one of 'DeltaG' or 'chembl'." % DG_prior)
    # Add to model.
    pymc_model['DeltaG'] = DeltaG

    # Create priors on true concentrations of protein and ligand.
    if concentration_priors == 'lognormal':
        Ptrue = pymc.Lognormal('Ptrue', mu=np.log(Pstated**2 / np.sqrt(dPstated**2 + Pstated**2)), tau=np.sqrt(np.log(1.0 + (dPstated/Pstated)**2))**(-2)) # protein concentration (M)
        Ltrue = pymc.Lognormal('Ltrue', mu=np.log(Lstated**2 / np.sqrt(dLstated**2 + Lstated**2)), tau=np.sqrt(np.log(1.0 + (dLstated/Lstated)**2))**(-2)) # ligand concentration (M)
        Ltrue_control = pymc.Lognormal('Ltrue_control', mu=np.log(Lstated**2 / np.sqrt(dLstated**2 + Lstated**2)), tau=np.sqrt(np.log(1.0 + (dLstated/Lstated)**2))**(-2)) # ligand concentration (M)
    elif concentration_priors == 'gaussian':
        # Warning: These priors could lead to negative concentrations.
        Ptrue = pymc.Normal('Ptrue', mu=Pstated, tau=dPstated**(-2)) # protein concentration (M)
        Ltrue = pymc.Normal('Ltrue', mu=Lstated, tau=dLstated**(-2)) # ligand concentration (M)
        Ltrue_control = pymc.Normal('Ltrue_control', mu=Lstated, tau=dLstated**(-2)) # ligand concentration (M)
    else:
        raise Exception("concentration_priors = '%s' unknown. Must be one of ['lognormal', 'normal']." % concentration_priors)
    # Add to model.
    pymc_model['Ptrue'] = Ptrue
    pymc_model['Ltrue'] = Ltrue
    pymc_model['Ltrue_control'] = Ltrue_control

    # extinction coefficient
    if use_primary_inner_filter_correction:
        if epsilon:
            epsilon = pymc.Lognormal('epsilon', mu=np.log(epsilon**2 / np.sqrt(depsilon**2 + epsilon**2)), tau=np.sqrt(np.log(1.0 + (depsilon/epsilon)**2))**(-2)) # prior is centered on measured extinction coefficient
            # TODO: Change this to lognormal, since epsilon cannot be negative
        else:
            epsilon = pymc.Uniform('epsilon', lower=0.0, upper=200e3) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm
            # TODO: Select a reasonable physical range for plausible extinction coefficients.
        # Add to model.
        pymc_model['epsilon'] = epsilon

    # Compute initial guesses.
    F_background_guess = min(Fobs_i.min(), Fligand_i.min())
    F_L_guess = ((Fligand_i.max() - F_background_guess) / Lstated.max())
    F_P_guess = 0.0
    F_PL_guess = ((Fobs_i.max() - F_background_guess) / Pstated)

    # Priors on fluorescence intensities of complexes (later divided by a factor of Pstated for scale).
    Fmax = max(Fobs_i.max(), Fligand_i.max())
    F_background = pymc.Uniform('F_background', lower=0.0, upper=Fmax, value=F_background_guess) # background fluorescence
    F_PL = pymc.Uniform('F_PL', lower=0.0, upper=2*Fmax/min(Pstated.max(),Lstated.max()), value=F_PL_guess) # complex fluorescence
    F_P = pymc.Uniform('F_P', lower=0.0, upper=2*(Fobs_i/Pstated).max(), value=F_P_guess) # protein fluorescence
    F_L = pymc.Uniform('F_L', lower=0.0, upper=2*(Fligand_i/Lstated).max(), value=F_L_guess) # ligand fluorescence
    # Add to model.
    pymc_model['F_background'] = F_background
    pymc_model['F_PL'] = F_PL
    pymc_model['F_P'] = F_P
    pymc_model['F_L'] = F_L

    # Unknown experimental measurement error.
    log_sigma = pymc.Uniform('log_sigma', lower=-10, upper=np.log(Fmax))
    @pymc.deterministic
    def precision(log_sigma=log_sigma): # measurement precision
        return np.exp(-2.0*log_sigma)
    # Add to model.
    pymc_model['log_sigma'] = log_sigma

    # Fluorescence model.
    from assaytools.bindingmodels import TwoComponentBindingModel
    @pymc.deterministic
    def Fmodel(F_background=F_background, F_PL=F_PL, F_P=F_P, F_L=F_L, Ptrue=Ptrue, Ltrue=Ltrue, DeltaG=DeltaG, epsilon=epsilon):
        if use_primary_inner_filter_correction:
            IF_i = np.exp(np.minimum(-epsilon*path_length*Ltrue[:], 0.0))
        else:
            IF_i = np.ones(N)
        [P_i, L_i, PL_i] = TwoComponentBindingModel.equilibrium_concentrations(DeltaG, Ptrue[:], Ltrue[:])
        Fmodel_i = IF_i[:]*(F_PL*PL_i + F_L*L_i + F_P*P_i + F_background)

        return Fmodel_i
    # Add to model.
    pymc_model['Fmodel'] = Fmodel

    # Fluorescence model, ligand only.
    @pymc.deterministic
    def Fligand(F_background=F_background, F_L=F_L, Ltrue_control=Ltrue_control, epsilon=epsilon):
        if use_primary_inner_filter_correction:
            IF_i = np.exp(np.minimum(-epsilon*path_length*Ltrue_control[:], 0.0))
        else:
            IF_i = np.ones(N)
        Fmodel_i = IF_i[:]*(F_L*Ltrue_control[:] + F_background)
        return Fmodel_i
    # Add to model.
    pymc_model['Fligand'] = Fligand

    # Experimental error on fluorescence observations.
    Fobs_model = pymc.Normal('Fobs_model', mu=Fmodel, tau=precision, size=[N], observed=True, value=Fobs_i) # observed data
    Fligand_model = pymc.Normal('Fligand_model', mu=Fligand, tau=precision, size=[N], observed=True, value=Fligand_i) # ligand only data
    # Add to model.
    pymc_model['Fobs_model'] = Fobs_model
    pymc_model['Fligand_model'] = Fligand_model

    # Promote this to a full-fledged PyMC model.
    pymc_model = pymc.Model(pymc_model)

    # Return the pymc model
    return pymc_model

