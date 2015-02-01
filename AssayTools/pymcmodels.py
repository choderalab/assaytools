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

def inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, concentration, geometry='top'):
    """
    Compute primary and secondar inner filter effect attenuation for top and bottom observation geometries.

    Parameters
    ----------
    epsilon_ex : float
       Exctinction coefficient at excitation wavelength. Units of 1/M/cm
    epsilon_em : float
       Extinction coefficient at emission wavelength. Units of 1/M/cm
    path_length : float
       Path length. Units of cm.
    concentration : float
       Concentration of species whose extinction coefficient is provided. Units of M.
    geometry : str, optional, default='top'
       Observation geometry, one of ['top', 'bottom'].

    Returns
    -------
    scaling : factor by which expected fluorescence is attenuated by primary and secondary inner filter effects

    """

    # Ensure concentration is a vector.
    if not hasattr(concentration, '__getitem__'):
        concentration = np.array([concentration])

    ELC_ex = epsilon_ex*path_length*concentration
    ELC_em = epsilon_em*path_length*concentration

    scaling = 1.0 # no attenuation

    if geometry == 'top':
        alpha = (ELC_ex + ELC_em)

        scaling = (1 - np.exp(-alpha)) / alpha
        # Handle alpha -> 0 case explicitly.
        indices = np.where(np.abs(alpha) < 0.01)
        scaling[indices] = 1. - alpha[indices]/2. + (alpha[indices]**2)/6. - (alpha[indices]**3)/24. + (alpha[indices]**4)/120.
    elif geometry == 'bottom':
        alpha = (ELC_ex - ELC_em)

        scaling = (1 - np.exp(-alpha)) / alpha
        # Handle alpha -> 0 case explicitly.
        indices = np.where(np.abs(alpha) < 0.01)
        scaling[indices] = 1. - alpha[indices]/2. + (alpha[indices]**2)/6. - (alpha[indices]**3)/24. + (alpha[indices]**4)/120.
        # Include additional term.
        scaling *= np.exp(-ELC_em)
    else:
        raise Exception("geometry '%s' unknown, must be one of ['top', 'bottom']" % geometry)

    return scaling

# Create a pymc model
def make_model(Pstated, dPstated, Lstated, dLstated,
               top_complex_fluorescence=None, top_ligand_fluorescence=None,
               bottom_complex_fluorescence=None, bottom_ligand_fluorescence=None,
               DG_prior='uniform',
               concentration_priors='lognormal',
               use_primary_inner_filter_correction=True,
               use_secondary_inner_filter_correction=True,
               assay_volume=100e-6, well_area=0.1586,
               epsilon_ex=None, depsilon_ex=None,
               epsilon_em=None, depsilon_em=None,
               ligand_ex_absorbance=None, ligand_em_absorbance=None):
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
    top_complex_fluorecence : numpy.array of N values, optional, default=None
       Fluorescence intensity (top) for protein:ligand mixture.
    top_ligand_fluorescence : numpy.array of N values, optional, default=None
       Fluorescence intensity (top) for ligand control.
    bottom_complex_fluorescence: numpy.array of N values, optional, default=None
       Fluorescence intensity (bottom) for protein:ligand mixture.
    bottom_ligand_fluorescence : numpy.array of N values, optional, default=None
       Fluorescence intensity (bottom) for ligand control.
    DG_prior : str, optional, default='uniform'
       Prior to use for reduced free energy of binding (DG): 'uniform' (uniform over reasonable range), or 'chembl' (ChEMBL-inspired distribution); default: 'uniform'
    concentration_priors : str, optional, default='lognormal'
       Prior to use for protein and ligand concentrations. Available options are ['lognormal', 'normal'].
    use_primary_inner_filter_correction : bool, optional, default=True
       If true, will infer ligand extinction coefficient epsilon and apply primary inner filter correction to attenuate excitation light.
    use_secondary_inner_filter_correction : bool, optional, default=True
       If true, will infer ligand extinction coefficient epsilon and apply secondary inner filter correction to attenuate excitation light.
    assay_volume : float, optional, default=100e-6
       Assay volume. Units of L. Default 100 uL.
    well_area : float, optional, default=0.1586
       Well area. Units of cm^2. Default 0.1586 cm^2, for half-area plate.
    epsilon_ex, depsilon_ex : float, optional, default=None
       Orthogonal measurement of ligand extinction coefficient at excitation wavelength (and uncertainty). If None, will use a uniform prior.
    epsilon_em, depsilon_em : float, optional, default=None
       Orthogonal measurement of ligand extinction coefficient at excitation wavelength (and uncertainty). If None, will use a uniform prior.
    ligand_ex_absorbance : np.array of N values, optional, default=None
       Ligand absorbance measurement for excitation wavelength.
    ligand_em_absorbance : np.array of N values, optional, default=None
       Ligand absorbance measurement for emission wavelength.

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
    >>> top_complex_fluoresnce = array([ 689., 683., 664., 588., 207., 80., 28., 17., 10., 11., 10., 10.])
    >>> top_ligand_fluorescence = array([ 174., 115., 57., 20., 7., 6., 6., 6., 6., 7., 6., 7.])
    >>> pymc_model = pymcmodels.make_model(Pstated, dPstated, Lstated, dLstated, top_complex_fluorescence=top_complex_fluorescence, bottom_complex_fluorescence=bottom_complex_fluorescence)

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
    model = dict()

    # Prior on binding free energies.
    if DG_prior == 'uniform':
        DeltaG = pymc.Uniform('DeltaG', lower=DG_min, upper=DG_max) # binding free energy (kT), uniform over huge range
    elif DG_prior == 'chembl':
        DeltaG = pymc.Normal('DeltaG', mu=0, tau=1./(12.5**2)) # binding free energy (kT), using a Gaussian prior inspured by ChEMBL
    else:
        raise Exception("DG_prior = '%s' unknown. Must be one of 'DeltaG' or 'chembl'." % DG_prior)
    # Add to model.
    model['DeltaG'] = DeltaG

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
    model['Ptrue'] = Ptrue
    model['Ltrue'] = Ltrue
    model['Ltrue_control'] = Ltrue_control

    # extinction coefficient
    model['epsilon_ex'] = 0.0
    model['epsilon_em'] = 0.0
    if use_primary_inner_filter_correction:
        if epsilon_ex:
            model['epsilon_ex'] = pymc.Lognormal('epsilon_ex', mu=np.log(epsilon_ex**2 / np.sqrt(depsilon_ex**2 + epsilon_ex**2)), tau=np.sqrt(np.log(1.0 + (depsilon_ex/epsilon_ex)**2))**(-2)) # prior is centered on measured extinction coefficient
        else:
            model['epsilon_ex'] = pymc.Uniform('epsilon_ex', lower=0.0, upper=1000e3, value=70000.0) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm

    if use_secondary_inner_filter_correction:
        if epsilon_em:
            model['epsilon_em'] = pymc.Lognormal('epsilon_em', mu=np.log(epsilon_em**2 / np.sqrt(depsilon_em**2 + epsilon_em**2)), tau=np.sqrt(np.log(1.0 + (depsilon_em/epsilon_em)**2))**(-2)) # prior is centered on measured extinction coefficient
        else:
            model['epsilon_em'] = pymc.Uniform('epsilon_em', lower=0.0, upper=1000e3, value=0.0) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm

    # Min and max observed fluorescence.
    Fmax = 0.0; Fmin = 1e6;
    if top_complex_fluorescence is not None:
        Fmax = max(Fmax, top_complex_fluorescence.max()); Fmin = min(Fmin, top_complex_fluorescence.min())
    if top_ligand_fluorescence is not None:
        Fmax = max(Fmax, top_ligand_fluorescence.max()); Fmin = min(Fmin, top_ligand_fluorescence.min())
    if bottom_complex_fluorescence is not None:
        Fmax = max(Fmax, bottom_complex_fluorescence.max()); Fmin = min(Fmin, bottom_complex_fluorescence.min())
    if bottom_ligand_fluorescence is not None:
        Fmax = max(Fmax, bottom_ligand_fluorescence.max()); Fmin = min(Fmin, bottom_ligand_fluorescence.min())

    # Compute initial guesses.
    # TODO: Generalize this for top and bottom.
    F_plate_guess = Fmin
    F_buffer_guess = Fmin / path_length
    F_L_guess = (top_ligand_fluorescence.max() - top_ligand_fluorescence.min()) / Lstated.max()
    F_P_guess = 0.0
    F_PL_guess = (Fmax - Fmin) / min(Pstated.max(), Lstated.max())
    log_gain_guess = - np.log((top_complex_fluorescence.max() - top_complex_fluorescence.min()) / (bottom_complex_fluorescence.max() - bottom_complex_fluorescence.min()))
    print "Fmin = %.1f, Fmax = %.1f" % (Fmin, Fmax)

    # Priors on fluorescence intensities of complexes (later divided by a factor of Pstated for scale).
    model['F_plate'] = pymc.Uniform('F_plate', lower=0.0, upper=Fmax, value=F_plate_guess) # plate fluorescence
    model['F_buffer'] = pymc.Uniform('F_buffer', lower=0.0, upper=Fmax/path_length, value=F_buffer_guess) # buffer fluorescence
    model['F_PL'] = pymc.Uniform('F_PL', lower=0.0, upper=2*Fmax/min(Pstated.max(),Lstated.max()), value=F_PL_guess) # complex fluorescence
    model['F_P'] = pymc.Uniform('F_P', lower=0.0, upper=2*(Fmax/Pstated).max(), value=F_P_guess) # protein fluorescence
    model['F_L'] = pymc.Uniform('F_L', lower=0.0, upper=2*(Fmax/Lstated).max(), value=F_L_guess) # ligand fluorescence

    # Unknown experimental measurement error.
    model['log_sigma_top'] = pymc.Uniform('log_sigma_top', lower=-10, upper=np.log(Fmax), value=np.log(5))
    model['sigma_top'] = pymc.Lambda('sigma_top', lambda log_sigma=model['log_sigma_top'] : np.exp(log_sigma) )
    model['precision_top'] = pymc.Lambda('precision_top', lambda log_sigma=model['log_sigma_top'] : np.exp(-2*log_sigma) )

    model['log_sigma_bottom'] = pymc.Uniform('log_sigma_bottom', lower=-10, upper=np.log(Fmax), value=np.log(5))
    model['sigma_bottom'] = pymc.Lambda('sigma_bottom', lambda log_sigma=model['log_sigma_bottom'] : np.exp(log_sigma) )
    model['precision_bottom'] = pymc.Lambda('precision_bottom', lambda log_sigma=model['log_sigma_bottom'] : np.exp(-2*log_sigma) )

    model['log_gain_bottom'] = pymc.Uniform('log_gain_bottom', lower=-6.0, upper=6.0, value=log_gain_guess) # plate material absorbance at emission wavelength
    model['gain_bottom'] = pymc.Lambda('gain_bottom', lambda log_gain_bottom=model['log_gain_bottom'] : np.exp(log_gain_bottom) )

    model['log_sigma_abs'] = pymc.Uniform('log_sigma_abs', lower=-10, upper=0, value=np.log(0.01))
    model['sigma_abs'] = pymc.Lambda('sigma_abs', lambda log_sigma=model['log_sigma_abs'] : np.exp(log_sigma) )
    model['precision_abs'] = pymc.Lambda('precision_abs', lambda log_sigma=model['log_sigma_abs'] : np.exp(-2*log_sigma) )

    # Fluorescence model.
    from assaytools.bindingmodels import TwoComponentBindingModel

    if top_complex_fluorescence is not None:
        @pymc.deterministic
        def top_complex_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                           F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                           Ptrue=Ptrue, Ltrue=Ltrue, DeltaG=DeltaG,
                                           epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            [P_i, L_i, PL_i] = TwoComponentBindingModel.equilibrium_concentrations(DeltaG, Ptrue[:], Ltrue[:])
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, L_i, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*L_i) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_PL*PL_i + F_L*L_i + F_P*P_i + F_buffer*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_complex_fluorescence_model'] = top_complex_fluorescence_model
        model['top_complex_fluorescence'] = pymc.Normal('top_complex_fluorescence',
                                                        mu=model['top_complex_fluorescence_model'], tau=model['precision_top'],
                                                        size=[N], observed=True, value=top_complex_fluorescence) # observed data

    if top_ligand_fluorescence is not None:
        @pymc.deterministic
        def top_ligand_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                          F_L=model['F_L'],
                                          Ltrue=Ltrue,
                                          epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*Ltrue) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_L*Ltrue + F_buffer*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_ligand_fluorescence_model'] = top_ligand_fluorescence_model
        model['top_ligand_fluorescence'] = pymc.Normal('top_ligand_fluorescence',
                                                       mu=model['top_ligand_fluorescence_model'], tau=model['precision_top'],
                                                       size=[N], observed=True, value=top_ligand_fluorescence) # observed data

    if bottom_complex_fluorescence is not None:
        @pymc.deterministic
        def bottom_complex_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                              F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                              Ptrue=Ptrue, Ltrue=Ltrue, DeltaG=DeltaG,
                                              epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em'],
                                              log_gain_bottom=model['log_gain_bottom']):
            [P_i, L_i, PL_i] = TwoComponentBindingModel.equilibrium_concentrations(DeltaG, Ptrue[:], Ltrue[:])
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, L_i, geometry='bottom')
            IF_i_plate = np.exp(-epsilon_ex*path_length*L_i) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_PL*PL_i + F_L*L_i + F_P*P_i + F_buffer*path_length)*np.exp(log_gain_bottom) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['bottom_complex_fluorescence_model'] = bottom_complex_fluorescence_model
        model['bottom_complex_fluorescence'] = pymc.Normal('bottom_complex_fluorescence',
                                                           mu=model['bottom_complex_fluorescence_model'], tau=model['precision_bottom'],
                                                           size=[N], observed=True, value=bottom_complex_fluorescence) # observed data

    if bottom_ligand_fluorescence is not None:
        @pymc.deterministic
        def bottom_ligand_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                             F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                             Ltrue=Ltrue,
                                             epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em'],
                                             log_gain_bottom=model['log_gain_bottom']):
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue, geometry='bottom')
            IF_i_plate = np.exp(-epsilon_ex*path_length*Ltrue) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_L*Ltrue + F_buffer*path_length)*np.exp(log_gain_bottom) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['bottom_ligand_fluorescence_model'] = bottom_ligand_fluorescence_model
        model['bottom_ligand_fluorescence'] = pymc.Normal('bottom_ligand_fluorescence',
                                                          mu=model['bottom_ligand_fluorescence_model'], tau=model['precision_bottom'],
                                                          size=[N], observed=True, value=bottom_ligand_fluorescence) # observed data

    if ligand_ex_absorbance is not None:
        model['plate_abs_ex'] = pymc.Uniform('plate_abs_ex', lower=0.0, upper=1.0, value=ligand_ex_absorbance.min())
        @pymc.deterministic
        def ligand_ex_absorbance_model(Ltrue=Ltrue,
                                       epsilon_ex=model['epsilon_ex'],
                                       plate_abs_ex=model['plate_abs_ex']):
            Fmodel_i = (1.0 - np.exp(-epsilon_ex*path_length*Ltrue)) + plate_abs_ex
            return Fmodel_i
        # Add to model.
        model['ligand_ex_absorbance_model'] = ligand_ex_absorbance_model
        model['ligand_ex_absorbance'] = pymc.Normal('ligand_ex_absorbance',
                                                    mu=model['ligand_ex_absorbance_model'], tau=model['precision_abs'],
                                                    size=[N], observed=True, value=ligand_ex_absorbance) # observed data

    if ligand_em_absorbance is not None:
        model['plate_abs_em'] = pymc.Uniform('plate_abs_em', lower=0.0, upper=1.0, value=ligand_em_absorbance.min())
        @pymc.deterministic
        def ligand_em_absorbance_model(Ltrue=Ltrue,
                                       epsilon_em=model['epsilon_em'],
                                       plate_abs_em=model['plate_abs_em']):
            Fmodel_i = (1.0 - np.exp(-epsilon_em*path_length*Ltrue)) + plate_abs_em
            return Fmodel_i
        # Add to model.
        model['ligand_em_absorbance_model'] = ligand_em_absorbance_model
        model['ligand_em_absorbance'] = pymc.Normal('ligand_em_absorbance',
                                                    mu=model['ligand_em_absorbance_model'], tau=model['precision_abs'],
                                                    size=[N], observed=True, value=ligand_em_absorbance) # observed data

    # Promote this to a full-fledged PyMC model.
    pymc_model = pymc.Model(model)

    # Return the pymc model
    return pymc_model

