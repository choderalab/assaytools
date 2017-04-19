"""
pymc models for the analysis of fluorescence assay data.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy as np
import pymc

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
# SUBROUTINES
#=============================================================================================

def LogNormalWrapper(name, mean, stddev, log_prefix='log_', size=1, observed=False, value=None):
    """
    Create a PyMC Normal stochastic, automatically converting parameters to be appropriate for a `LogNormal` distribution.
    Note that the resulting distribution is Normal, not LogNormal.
    This is appropriate if we want to describe the log of a volume or a fluorescence reading.

    Notes
    -----
    Everything is coerced into a 1D array.

    Parameters
    ----------
    mean : float or array-like (but not pymc variable)
       Mean of exp(X), where X is lognormal variate.
    stddev : float or array-like (but not pymc variable)
       Standard deviation of exp(X), where X is lognormal variate.
    log_prefix : str, optional, default='log_'
       Prefix appended to stochastic
    size : list of int, optional, default=1
       Size vector
    observed : bool, optional, default=False
       If True, the stochastic is fixed to this observed value.
    value : float, optional, default=None
       If observed=True, the observed value of the real quantity

    Returns
    -------
    stochastic : pymc.Normal
       Normal stochastic for random variable X.
       Name is prefixed with log_prefix
    deterministic : pymc.deterministic
       Deterministic encoding the exponentiated lognormal (real quantity)

    """
    if observed and (value is not None):
        # Compute parameters of lognormal distribution
        @pymc.deterministic(name=name + '_mu')
        def mu(mean=mean, stddev=stddev):
            mu = np.log(mean**2 / np.sqrt(stddev**2 + mean**2))
            return mu
        @pymc.deterministic(name=name + '_tau')
        def tau(mean=mean, stddev=stddev):
            tau = np.sqrt(np.log(1.0 + (stddev/mean)**2))**(-2)
            return tau
        stochastic = pymc.Normal(log_prefix + name, mu=mu, tau=tau, size=size, observed=observed, value=np.log(value))

        # Create deterministic from stochastic
        @pymc.deterministic(name=name)
        def deterministic(log_value=stochastic):
            return value

        return stochastic, deterministic

    #
    # Handle non-observed case
    #

    # Promote to array
    if size == 1:
        mean, stddev = np.array([mean]), np.array([stddev])
        size = [1]

    # Handle only strictly positive elements---all others are set to zero as constants
    try:
        nonzero_indices = np.where(mean > 0)[0]
        zero_indices = np.where(mean <= 0)[0]
    except:
        nonzero_indices = range(size[0])
        zero_indices = []
    nnonzero = len(nonzero_indices)
    nzeros = len(zero_indices)

    # Compute parameters of lognormal distribution
    @pymc.deterministic(name=name + '_mu')
    def mu(mean=mean, stddev=stddev):
        if not np.isscalar(mean):
            mean = mean[nonzero_indices]
        if not np.isscalar(stddev):
            stddev = stddev[nonzero_indices]
        mu = np.log(mean**2 / np.sqrt(stddev**2 + mean**2))
        return mu
    @pymc.deterministic(name=name + '_tau')
    def tau(mean=mean, stddev=stddev):
        if not np.isscalar(mean):
            mean = mean[nonzero_indices]
        if not np.isscalar(stddev):
            stddev = stddev[nonzero_indices]
        tau = np.sqrt(np.log(1.0 + (stddev/mean)**2))**(-2)
        return tau
    stochastic = pymc.Normal(log_prefix + name, mu=mu, tau=tau, size=nnonzero)

    # Create deterministic from stochastic
    @pymc.deterministic(name=name)
    def deterministic(log_value=stochastic):
        value = np.zeros(size, np.float64)
        value[nonzero_indices] = np.exp(log_value)
        return value

    return stochastic, deterministic

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
               quantum_yield_priors='lognormal',
               use_primary_inner_filter_correction=True,
               use_secondary_inner_filter_correction=True,
               assay_volume=100e-6, well_area=0.1586,
               epsilon_ex=None, depsilon_ex=None,
               epsilon_em=None, depsilon_em=None,
               ligand_ex_absorbance=None, ligand_em_absorbance=None,
               link_top_and_bottom_sigma=True,
               F_PL=None, dF_PL=None):
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
       Zero concentrations will automatically be handled as special cases.
    dLstated : numpy.array of N values
       Absolute uncertainty in stated protein concentrations for all wells of assay. Units of molarity.
       Uncertainties cannot be zero if the corresponding concentration is nonzero.
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
       Prior to use for protein and ligand concentrations. Available options are ['lognormal'].
    quantum_yield_priors : str, optional, default='lognormal'
       Prior to use for quantum yields. Available options are ['lognormal', 'uniform'].
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
    link_top_and_bottom_sigma : bool, optional, default=True
       If True, will link top and bottom fluorescence uncertainty sigma.
    F_PL, dF_PL : float, optional, default=None
       Prior on fluorescence quantum yield of protein-ligand complex, in units of RFU/M.
       Both F_PL and dF_PL must be specified for a lognormal prior to be used;
       otherwise a uniform prior is used.

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
    >>> top_complex_fluorescence = np.array([ 689., 683., 664., 588., 207., 80., 28., 17., 10., 11., 10., 10.], np.float32)
    >>> top_ligand_fluorescence = np.array([ 174., 115., 57., 20., 7., 6., 6., 6., 6., 7., 6., 7.], np.float32)
    >>> from pymcmodels import make_model
    >>> pymc_model = make_model(Pstated, dPstated, Lstated, dLstated, top_complex_fluorescence=top_complex_fluorescence, top_ligand_fluorescence=top_ligand_fluorescence)

    """

    # Compute path length.
    path_length = assay_volume * 1000 / well_area # cm, needed for inner filter effect corrections

    # Compute number of samples.
    N = len(Lstated)

    # Check input.
    # TODO: Check fluorescence and absorbance measurements for correct dimensions.
    if (len(Pstated) != N):
        raise Exception('len(Pstated) [%d] must equal len(Lstated) [%d].' % (len(Pstated), len(Lstated)))
    if (len(dPstated) != N):
        raise Exception('len(dPstated) [%d] must equal len(Lstated) [%d].' % (len(dPstated), len(Lstated)))
    if (len(dLstated) != N):
        raise Exception('len(dLstated) [%d] must equal len(Lstated) [%d].' % (len(dLstated), len(Lstated)))

    # Note whether we have top or bottom fluorescence measurements.
    top_fluorescence = (top_complex_fluorescence is not None) or (top_ligand_fluorescence is not None) # True if any top fluorescence measurements provided
    bottom_fluorescence = (bottom_complex_fluorescence is not None) or (bottom_ligand_fluorescence is not None) # True if any bottom fluorescence measurements provided

    # Create an empty dict to hold the model.
    model = dict()

    # Prior on binding free energies.
    if DG_prior == 'uniform':
        DeltaG = pymc.Uniform('DeltaG', lower=DG_min, upper=DG_max) # binding free energy (kT), uniform over huge range
    elif DG_prior == 'chembl':
        DeltaG = pymc.Normal('DeltaG', mu=0, tau=1./(12.5**2)) # binding free energy (kT), using a Gaussian prior inspured by ChEMBL
    else:
        raise Exception("DG_prior = '%s' unknown. Must be one of 'uniform' or 'chembl'." % DG_prior)
    # Add to model.
    model['DeltaG'] = DeltaG

    # Create priors on true concentrations of protein and ligand.
    if concentration_priors != 'lognormal':
        raise Exception("concentration_priors = '%s' unknown. Must be one of ['lognormal']." % concentration_priors)
    model['log_Ptrue'], model['Ptrue'] = LogNormalWrapper('Ptrue', mean=Pstated, stddev=dPstated, size=Pstated.shape) # protein concentration (M)
    model['log_Ltrue'], model['Ltrue'] = LogNormalWrapper('Ltrue', mean=Lstated, stddev=dLstated, size=Lstated.shape) # ligand concentration (M)
    model['log_Ltrue_control'], model['Ltrue_control'] = LogNormalWrapper('Ltrue_control', mean=Lstated, stddev=dLstated, size=Lstated.shape) # ligand concentration (M)

    # extinction coefficient
    model['epsilon_ex'] = 0.0
    if use_primary_inner_filter_correction:
        if epsilon_ex:
            model['log_epsilon_ex'], model['epsilon_ex'] = LogNormalWrapper('epsilon_ex', mean=epsilon_ex, stddev=depsilon_ex) # prior is centered on measured extinction coefficient
        else:
            model['epsilon_ex'] = pymc.Uniform('epsilon_ex', lower=0.0, upper=1000e3, value=70000.0) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm

    model['epsilon_em'] = 0.0
    if use_secondary_inner_filter_correction:
        if epsilon_em:
            model['log_epsilon_em'], model['epsilon_em'] = LogNormalWrapper('epsilon_em', mean=epsilon_em, stddev=depsilon_em) # prior is centered on measured extinction coefficient
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

    # Compute initial guesses for fluorescence quantum yield quantities.
    F_plate_guess = Fmin
    F_buffer_guess = Fmin / path_length
    F_L_guess = (Fmax - Fmin) / Lstated.max()
    F_P_guess = Fmin
    F_P_guess = Fmin / Pstated.min()
    F_PL_guess = (Fmax - Fmin) / min(Pstated.max(), Lstated.max())

    # Priors on fluorescence intensities of complexes (later divided by a factor of Pstated for scale).
    if quantum_yield_priors == 'uniform':
        model['F_plate'] = pymc.Uniform('F_plate', lower=0.0, upper=Fmax, value=F_plate_guess) # plate fluorescence
        model['F_buffer'] = pymc.Uniform('F_buffer', lower=0.0, upper=Fmax/path_length, value=F_buffer_guess) # buffer fluorescence
        if (F_PL is not None) and (dF_PL is not None):
            model['log_F_PL'], model['F_PL'] = LogNormalWrapper('F_PL', mean=F_PL, stddev=dF_PL)
        else:
            model['F_PL'] = pymc.Uniform('F_PL', lower=0.0, upper=2*Fmax/min(Pstated.max(),Lstated.max()), value=F_PL_guess) # complex fluorescence
        model['F_P'] = pymc.Uniform('F_P', lower=0.0, upper=2*(Fmax/Pstated).max(), value=F_P_guess) # protein fluorescence
        model['F_L'] = pymc.Uniform('F_L', lower=0.0, upper=2*(Fmax/Lstated).max(), value=F_L_guess) # ligand fluorescence
    elif quantum_yield_priors == 'lognormal':
        stddev = 1.0 # relative factor for stddev guess
        model['log_F_plate'], model['F_plate'] = LogNormalWrapper('F_plate', mean=F_plate_guess, stddev=stddev*F_plate_guess) # plate fluorescence
        model['log_F_buffer'], model['F_buffer'] = LogNormalWrapper('F_buffer', mean=F_buffer_guess, stddev=stddev*F_buffer_guess) # buffer fluorescence
        if (F_PL is not None) and (dF_PL is not None):
            model['log_F_PL'], model['F_PL'] = LogNormalWrapper('F_PL', mean=F_PL, stddev=dF_PL)
        else:
            model['log_F_PL'], model['F_PL'] = LogNormalWrapper('F_PL', mean=F_PL_guess, stddev=stddev*F_PL_guess) # complex fluorescence
        model['log_F_P'], model['F_P'] = LogNormalWrapper('F_P', mean=F_P_guess, stddev=stddev*F_P_guess) # protein fluorescence
        model['log_F_L'], model['F_L'] = LogNormalWrapper('F_L', mean=F_L_guess, stddev=stddev*F_L_guess) # ligand fluorescence
    else:
        raise Exception("quantum_yield_priors = '%s' unknown. Must be one of ['lognormal', 'uniform']." % quantum_yield_priors)

    # Unknown experimental measurement error.
    if top_fluorescence:
        model['log_sigma_top'] = pymc.Uniform('log_sigma_top', lower=-10, upper=np.log(Fmax), value=np.log(5))
        model['sigma_top'] = pymc.Lambda('sigma_top', lambda log_sigma=model['log_sigma_top'] : np.exp(log_sigma) )
        model['precision_top'] = pymc.Lambda('precision_top', lambda log_sigma=model['log_sigma_top'] : np.exp(-2*log_sigma) )

    if bottom_fluorescence:
        if top_fluorescence and bottom_fluorescence and link_top_and_bottom_sigma:
            # Use the same log_sigma for top and bottom fluorescence
            model['log_sigma_bottom'] = pymc.Lambda('log_sigma_bottom', lambda log_sigma_top=model['log_sigma_top'] : log_sigma_top )
        else:
            model['log_sigma_bottom'] = pymc.Uniform('log_sigma_bottom', lower=-10, upper=np.log(Fmax), value=np.log(5))
        model['sigma_bottom'] = pymc.Lambda('sigma_bottom', lambda log_sigma=model['log_sigma_bottom'] : np.exp(log_sigma) )
        model['precision_bottom'] = pymc.Lambda('precision_bottom', lambda log_sigma=model['log_sigma_bottom'] : np.exp(-2*log_sigma) )

    if top_fluorescence and bottom_fluorescence:
        # Gain that attenuates bottom fluorescence relative to top.
        # TODO: Replace this with plate absorbance?
        log_gain_guess = - np.log((top_complex_fluorescence.max() - top_complex_fluorescence.min()) / (bottom_complex_fluorescence.max() - bottom_complex_fluorescence.min()))
        model['log_gain_bottom'] = pymc.Uniform('log_gain_bottom', lower=-6.0, upper=6.0, value=log_gain_guess) # plate material absorbance at emission wavelength
        model['gain_bottom'] = pymc.Lambda('gain_bottom', lambda log_gain_bottom=model['log_gain_bottom'] : np.exp(log_gain_bottom) )
    elif (not top_fluorescence) and bottom_fluorescence:
        model['log_gain_bottom'] = 0.0 # no gain
        model['gain_bottom'] = pymc.Lambda('gain_bottom', lambda log_gain_bottom=model['log_gain_bottom'] : np.exp(log_gain_bottom) )

    if top_fluorescence:
        model['log_sigma_abs'] = pymc.Uniform('log_sigma_abs', lower=-10, upper=0, value=np.log(0.01))
        model['sigma_abs'] = pymc.Lambda('sigma_abs', lambda log_sigma=model['log_sigma_abs'] : np.exp(log_sigma) )
        model['precision_abs'] = pymc.Lambda('precision_abs', lambda log_sigma=model['log_sigma_abs'] : np.exp(-2*log_sigma) )

    # Fluorescence model.
    from assaytools.bindingmodels import TwoComponentBindingModel

    if top_complex_fluorescence is not None:
        @pymc.deterministic
        def top_complex_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                           F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                           Ptrue=model['Ptrue'], Ltrue=model['Ltrue'], DeltaG=model['DeltaG'],
                                           epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            [P_i, L_i, PL_i] = TwoComponentBindingModel.equilibrium_concentrations(DeltaG, Ptrue[:], Ltrue[:])
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, L_i, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*L_i) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_PL*PL_i + F_L*L_i + F_P*P_i + F_buffer*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_complex_fluorescence_model'] = top_complex_fluorescence_model
        model['log_top_complex_fluorescence'], model['top_complex_fluorescence'] = LogNormalWrapper('top_complex_fluorescence',
            mean=model['top_complex_fluorescence_model'], stddev=model['sigma_top'],
            size=[N], observed=True, value=top_complex_fluorescence) # observed data

    if top_ligand_fluorescence is not None:
        @pymc.deterministic
        def top_ligand_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                          F_L=model['F_L'],
                                          Ltrue=model['Ltrue'],
                                          epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em']):
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue, geometry='top')
            IF_i_plate = np.exp(-(epsilon_ex+epsilon_em)*path_length*Ltrue) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_L*Ltrue + F_buffer*path_length) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['top_ligand_fluorescence_model'] = top_ligand_fluorescence_model
        model['log_top_ligand_fluorescence'], model['top_ligand_fluorescence'] = LogNormalWrapper('top_ligand_fluorescence',
                                                       mean=model['top_ligand_fluorescence_model'], stddev=model['sigma_top'],
                                                       size=[N], observed=True, value=top_ligand_fluorescence) # observed data

    if bottom_complex_fluorescence is not None:
        @pymc.deterministic
        def bottom_complex_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                              F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                              Ptrue=model['Ptrue'], Ltrue=model['Ltrue'], DeltaG=model['DeltaG'],
                                              epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em'],
                                              log_gain_bottom=model['log_gain_bottom']):
            [P_i, L_i, PL_i] = TwoComponentBindingModel.equilibrium_concentrations(DeltaG, Ptrue[:], Ltrue[:])
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, L_i, geometry='bottom')
            IF_i_plate = np.exp(-epsilon_ex*path_length*L_i) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_PL*PL_i + F_L*L_i + F_P*P_i + F_buffer*path_length)*np.exp(log_gain_bottom) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['bottom_complex_fluorescence_model'] = bottom_complex_fluorescence_model
        model['log_bottom_complex_fluorescence'], model['bottom_complex_fluorescence'] = LogNormalWrapper('bottom_complex_fluorescence',
                                                           mean=model['bottom_complex_fluorescence_model'], stddev=model['sigma_bottom'],
                                                           size=[N], observed=True, value=bottom_complex_fluorescence) # observed data

    if bottom_ligand_fluorescence is not None:
        @pymc.deterministic
        def bottom_ligand_fluorescence_model(F_plate=model['F_plate'], F_buffer=model['F_buffer'],
                                             F_PL=model['F_PL'], F_P=model['F_P'], F_L=model['F_L'],
                                             Ltrue=model['Ltrue'],
                                             epsilon_ex=model['epsilon_ex'], epsilon_em=model['epsilon_em'],
                                             log_gain_bottom=model['log_gain_bottom']):
            IF_i = inner_filter_effect_attenuation(epsilon_ex, epsilon_em, path_length, Ltrue, geometry='bottom')
            IF_i_plate = np.exp(-epsilon_ex*path_length*Ltrue) # inner filter effect applied only to plate
            Fmodel_i = IF_i[:]*(F_L*Ltrue + F_buffer*path_length)*np.exp(log_gain_bottom) + IF_i_plate*F_plate
            return Fmodel_i
        # Add to model.
        model['bottom_ligand_fluorescence_model'] = bottom_ligand_fluorescence_model
        model['log_bottom_ligand_fluorescence'], model['bottom_ligand_fluorescence'] = LogNormalWrapper('bottom_ligand_fluorescence',
                                                          mean=model['bottom_ligand_fluorescence_model'], stddev=model['sigma_bottom'],
                                                          size=[N], observed=True, value=bottom_ligand_fluorescence) # observed data

    # Absorbance measurements

    if ligand_ex_absorbance is not None:
        model['plate_abs_ex'] = pymc.Uniform('plate_abs_ex', lower=0.0, upper=1.0, value=ligand_ex_absorbance.min())
        @pymc.deterministic
        def ligand_ex_absorbance_model(Ltrue=model['Ltrue'],
                                       epsilon_ex=model['epsilon_ex'],
                                       plate_abs_ex=model['epsilon_em']):
            Fmodel_i = (1.0 - np.exp(-epsilon_ex*path_length*Ltrue)) + plate_abs_ex
            return Fmodel_i
        # Add to model.
        model['ligand_ex_absorbance_model'] = ligand_ex_absorbance_model
        model['log_ligand_ex_absorbance'], model['ligand_ex_absorbance'] = LogNoramlWrapper('ligand_ex_absorbance',
                                                    mean=model['ligand_ex_absorbance_model'], stddev=model['sigma_abs'],
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
        model['log_ligand_em_absorbance'], model['ligand_em_absorbance'] = LogNormalWrapper('ligand_em_absorbance',
                                                    mean=model['ligand_em_absorbance_model'], stddev=model['sigma_abs'],
                                                    size=[N], observed=True, value=ligand_em_absorbance) # observed data

    # Promote this to a full-fledged PyMC model.
    pymc_model = pymc.Model(model)

    # Return the pymc model
    return pymc_model

def map_fit(pymc_model):
    """
    Find the maximum a posteriori (MAP) fit.

    Parameters
    ----------
    pymc_model : pymc model
       The pymc model to sample.

    Returns
    -------
    map : pymc.MAP
       The MAP fit.

    """
    map = pymc.MAP(pymc_model)
    ncycles = 50

    # DEBUG
    ncycles = 5

    for cycle in range(ncycles):
        if (cycle+1)%5==0: print('MAP fitting cycle %d/%d' % (cycle+1, ncycles))
        map.fit()

    return map

def run_mcmc_emcee(pymc_model, nwalkers=100, nburn=100, niter=1000, nthin=None):
    """
    Sample the model with pymc using sensible defaults and the emcee sampler.

    Parameters
    ----------
    pymc_model : pymc model
       The pymc model to sample.
    nwalkers: int
        The number ensemble walkers.
    nburn: int
        The number of initial steps that will be discarded.
    niter: int
        The number of MCMC iterations that are performed after the burn-in period.
    nthin: int
        The frequency with which to discard samples. If None, then nthin=nwalkers.

    Returns
    -------
    pymc_model : pymc.Model.Model
        The PyMC object used to initialize the class.
    mcmc_model: pymc.MCMC.MCMC
        The PyMC object that contains the MCMC traces.

    """
    if nthin is None:
        nthin = nwalkers

    import emcee

    def unpack_parameters(model):
        """
        Takes the parameters from the a pymc, which could be floats or arrays and unpacks them into a single array.

        Parameter
        ---------
        model: pymc.Model.Model or pymc.MCMC.MCMC
            The PyMC object from which the parameters will be extracted.

        Returns
        -------
        parameters: numpy.ndarray
            Vectorized form of the pymc model parameters
        """
        parameters = []
        for stoch in model.stochastics:
            if stoch.value.shape == ():
                parameters.append(float(stoch.value))
            else:
                for val in stoch.value:
                    parameters.append(val)
        return np.array(parameters)

    def log_post(parameters, model):
        """
        Packs a vectorized form of the parameters back into the pymc model to allow the evaluation of the log posterior.

        Parameters
        ----------
        parameters: numpy.ndarray
            Vectorized form of the pymc model parameters
        model: pymc.Model.Model
            The pymc that contains all the variables that need to be inferred.

        Returns
        -------
        logp: float
            The log of the posterior density.
        """
        # Basic error handling in case emcee proposes moves outside the domain of the variables.
        try:
            p_ind = 0
            for stoch in model.stochastics:
                if stoch.value.shape == ():
                    stoch.value = np.array(parameters[p_ind])
                    p_ind += 1
                else:
                    stoch.value = parameters[p_ind:(p_ind + len(stoch.value))]
                    p_ind += len(stoch.value)
            logp = model.logp
        except pymc.ZeroProbability:
            logp = -np.inf

        return logp

    # Find MAP. This will be used to initialize the emcee sampler.
    pymc.MAP(pymc_model).fit()

    # Perform a dummy run with pymc to initial trace data structures
    mcmc_model = pymc.MCMC(pymc_model)
    mcmc_model.sample(1)

    # Defining the emcee parameters by extracting the stochastic variables from pymc
    parameters = unpack_parameters(mcmc_model)
    ndim = len(parameters)

    # sample starting points for walkers around the MAP
    p0 = np.zeros((nwalkers, ndim))
    for walker in range(nwalkers):
        # Initializing walkers to be within about 20% of the (local) MAP estimate
        p0[walker, :] = parameters + np.random.normal(0, 0.2 * np.absolute(parameters))

    # Initiate emcee sampler by passing the likelihood function
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_post, args=[mcmc_model])

    # Burn-in
    pos, prob, state = sampler.run_mcmc(p0, nburn)
    sampler.reset()

    # Production
    sampler.run_mcmc(pos, niter)

    # Packing the trace of the parameters into the PyMC MCMC object
    p_ind = 0
    for stoch in mcmc_model.stochastics:
        if stoch.value.shape == ():
            trace = sampler.flatchain[:, p_ind]
            stoch.trace._trace[0] = trace[::nthin].copy()
            p_ind += 1
        else:
            trace = sampler.flatchain[:, p_ind:(p_ind + len(stoch.value))]
            stoch.trace._trace[0] = trace[::nthin].copy()
            p_ind += len(stoch.value)

    return mcmc_model, pymc_model

def run_mcmc(pymc_model, nthin=50, nburn=0, niter=20000, map=True, db='ram', dbname=None):
    """
    Sample the model with pymc. Initial values of the parameters can be chosen with a maximum a posteriori estimate.

    Parameters
    ----------
    pymc_model : pymc model
       The pymc model to sample.
    nthin: int
        The number of MCMC steps that constitute 1 iteration.
    nburn: int
        The number of MCMC iterations during the burn-in.
    niter: int
        The number of production iterations.
    map: bool
        Whether to initialize the parameters before MCMC with the maximum a posteriori estimate.
    db : str
        How to store model, default = 'ram' means not storing it. To store model use storage = 'pickle'. If not,
        supply the name of the database backend that will store the values of the stochastics and deterministics sampled
        during the MCMC loop.
   dbname : str
        name for storage object, default = None. To store model use e.g. dbname = 'my_mcmc.pickle'

    Returns
    -------
    mcmc : pymc.MCMC
       The MCMC samples.

    """
    # Find MAP:
    if map == True:
        pymc.MAP(pymc_model).fit()

    # Sample the model with pymc
    mcmc = pymc.MCMC(pymc_model, db=db, dbname=dbname, name='Sampler', verbose=True)

    mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'DeltaG'), proposal_sd=1.0, proposal_distribution='Normal')
    mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'log_F_PL'), proposal_sd=1.0, proposal_distribution='Normal')
    mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'log_F_P'), proposal_sd=1.0, proposal_distribution='Normal')
    mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'log_F_L'), proposal_sd=1.0, proposal_distribution='Normal')
    mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'log_F_plate'), proposal_sd=1.0, proposal_distribution='Normal')
    mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'log_F_buffer'), proposal_sd=1.0, proposal_distribution='Normal')
    if hasattr(pymc_model, 'epsilon_ex'):
        mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'epsilon_ex'), proposal_sd=10000.0, proposal_distribution='Normal')
    if hasattr(pymc_model, 'epsilon_em'):
        mcmc.use_step_method(pymc.Metropolis, getattr(pymc_model, 'epsilon_em'), proposal_sd=10000.0, proposal_distribution='Normal')

    mcmc.sample(iter=(nburn+niter), burn=nburn, thin=nthin, progress_bar=False, tune_throughout=False)

    return mcmc

def show_summary(pymc_model, mcmc, map):
    """
    Show summary statistics of MCMC and MAP estimates.

    Parameters
    ----------
    pymc_model : pymc model
       The pymc model to sample.
    map : pymc.MAP
       The MAP fit.
    mcmc : pymc.MCMC
       MCMC samples

    TODO
    ----
    * Automatically determine appropriate number of decimal places from statistical uncertainty.
    * Automatically adjust concentration units (e.g. pM, nM, uM) depending on estimated affinity.

    """

    # Compute summary statistics.
    DeltaG = map.DeltaG.value
    dDeltaG = mcmc.DeltaG.trace().std()
    Kd = np.exp(map.DeltaG.value)
    dKd = np.exp(mcmc.DeltaG.trace()).std()
    text = "DeltaG = %.1f +- %.1f kT\n" % (DeltaG, dDeltaG)
    if (Kd < 1e-12):
        text += "Kd = %.1f nM +- %.1f fM" % (Kd/1e-15, dKd/1e-15)
    elif (Kd < 1e-9):
        text += "Kd = %.1f pM +- %.1f pM" % (Kd/1e-12, dKd/1e-12)
    elif (Kd < 1e-6):
        text += "Kd = %.1f nM +- %.1f nM" % (Kd/1e-9, dKd/1e-9)
    elif (Kd < 1e-3):
        text += "Kd = %.1f uM +- %.1f uM" % (Kd/1e-6, dKd/1e-6)
    elif (Kd < 1):
        text += "Kd = %.1f mM +- %.1f mM" % (Kd/1e-3, dKd/1e-3)
    else:
        text += "Kd = %.3e M +- %.3e M" % (Kd, dKd)
    text += '\n'
    print(text)
