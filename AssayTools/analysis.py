"""
Classes for the analysis of fluorescence assay data.

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

#=============================================================================================
# PyMC models
#=============================================================================================

class CompetitiveBindingAnalysis(object):
    def __init__(self, solutions, wells, receptor_name, ligand_names,
        DeltaG_prior='uniform'):
        """
        Parameters
        ----------
        solutions : dict of str : Solution
            `solutions[name]`` is the Solution object corresponding to component `name`
        wells : autoprotocol.container.WellGroup
            Group of wells with `contents` and `measurements` properties defined
        receptor_name : str
            Name of receptor
        ligand_names : list of str
            Names of ligands
        DeltaG_prior : str, optional, default='uniform'
            Prior to use on DeltaG values for reaction ['uniform', 'chembl']

        """
        # Store data
        self.solutions = solutions
        self.wells = wells
        self.receptor_name = receptor_name
        self.ligand_names = ligand_names

        # Create the PyMC model
        model = dict()

        # Keep track of parameter name groups
        self.parameter_names = dict()

        # Construct priors on true solution concentrations
        self.parameter_names['concentrations'] = list()
        concentration_unit = 'moles/liter'
        for (name, solution) in solutions.iteritems():
            name += ' concentration'
            mu = solution.concentration.m_as(concentration_unit) # M, gaussian mean
            sigma = solution.uncertainty.m_as(concentration_unit) # M, gaussian mean
            if (mu == 0.0):
                true_concentration = 0.0
            else:
                true_concentration = pymc.Lognormal(name, mu=np.log(mu**2 / np.sqrt(sigma**2 + mu**2)), tau=np.sqrt(np.log(1.0 + (sigma/mu)**2))**(-2)) # protein concentration (M)
            model[name] = true_concentration
            self.parameter_names['concentrations'].append(name)

        # Construct binding reactions for competition assay and assign priors to binding affinities
        self.parameter_names['DeltaGs'] = list()
        from bindingmodels import GeneralBindingModel
        self.complex_names = list()
        self.reactions = list()
        for ligand_name in ligand_names:
            # Create complex name
            complex_name = receptor_name + ':' + ligand_name
            self.complex_names.append(complex_name)
            # Create the DeltaG prior
            name = 'DeltaG (%s + %s -> %s)' % (receptor_name, ligand_name, complex_name) # form the name of the pymc variable
            if DeltaG_prior == 'uniform':
                DeltaG = pymc.Uniform(name, lower=DG_min, upper=DG_max) # binding free energy (kT), uniform over huge range
            elif DeltaG_prior == 'chembl':
                DeltaG = pymc.Normal(name, mu=0, tau=1./(12.5**2)) # binding free energy (kT), using a Gaussian prior inspured by ChEMBL
            else:
                raise Exception("DeltaG_prior = '%s' unknown. Must be one of 'uniform' or 'chembl'." % DeltaG_prior)
            model[name] = DeltaG
            # Create the reaction
            reaction = (DeltaG, {complex_name : -1, receptor_name : +1, ligand_name : +1})
            # Append reactions.
            self.reactions.append(reaction)
            self.parameter_names['DeltaGs'].append(name)

        # Construct priors dispensed volumes in each well
        volume_unit = 'liters' # volume unit used throughout
        self.parameter_names['dispensed_volumes'] = list()
        self.parameter_names['well_volumes'] = list()
        for well in wells:
            component_variables= list()
            component_coefficients = list()
            for component in well.properties['contents']:
                name = 'volume of %s dispensed into well %s' % (component, well.humanize()) # TODO: Use plate name as well in case we may have well spanning multiple plates
                (volume, error) = well.properties['contents'][component]
                mu = volume.m_as(volume_unit)
                sigma = error.m_as(volume_unit)
                volume_dispensed = pymc.Lognormal(name, mu=np.log(mu**2 / np.sqrt(sigma**2 + mu**2)), tau=np.sqrt(np.log(1.0 + (sigma/mu)**2))**(-2)) # uL
                model[name] = volume_dispensed
                self.parameter_names['dispensed_volumes'].append(name)

                component_variables.append(volume_dispensed)
                component_coefficients.append(+1.0)

            name = 'volume of well %s' % well.humanize() # TODO: Use plate name too in case wells span multiple plates
            model[name] = pymc.LinearCombination(name, component_coefficients, component_variables)
            self.parameter_names['well_volumes'].append(name)

        # Create the PyMC Model object.
        self.model = pymc.Model(model)

        print('Model has %d stochastics and %d deterministics...' % (len(self.model.stochastics), len(self.model.deterministics)))

    def map_fit(self):
        """
        Find the maximum a posteriori (MAP) fit.

        Parameters
        ----------

        Returns
        -------
        map : pymc.MAP
           The MAP fit.

        """
        map = pymc.MAP(self.model)
        ncycles = 50
        nshow = 5

        # DEBUG
        ncycles = 5
        nshow = 1

        # Perform MAP fit
        print('Performing MAP fit...')
        for cycle in range(ncycles):
            if (cycle+1)%nshow==0: print('MAP fitting cycle %d/%d' % (cycle+1, ncycles))
            map.fit()

        return map

    def run_mcmc(self):
        """
        Sample the model with pymc using sensible defaults.

        Parameters
        ----------

        Returns
        -------
        mcmc : pymc.MCMC
           The MCMC samples.

        """

        # Sample the model with pymc
        mcmc = pymc.MCMC(self.model, db='ram', name='Sampler', verbose=True)
        nthin = 20
        nburn = nthin*10000
        niter = nthin*10000

        # DEBUG
        nburn = nthin*1000
        niter = nthin*1000

        for stochastic in self.model.stochastics:
            mcmc.use_step_method(pymc.Metropolis, stochastic, proposal_sd=1.0, proposal_distribution='Normal')

        print('Running MCMC...')
        mcmc.sample(iter=(nburn+niter), burn=nburn, thin=nthin, progress_bar=False, tune_throughout=False)

        return mcmc

    def show_summary(self, mcmc, map):
        """
        Show summary statistics of MCMC and MAP estimates.

        Parameters
        ----------
        map : pymc.MAP
           The MAP fit.
        mcmc : pymc.MCMC
           MCMC samples

        TODO
        ----
        * Automatically determine appropriate number of decimal places from statistical uncertainty.
        * Automatically adjust concentration units (e.g. pM, nM, uM) depending on estimated affinity.

        """

        # Compute summary statistics
        alpha = 0.95 # confidence interval width
        from scipy.stats import bayes_mvs
        for name in self.parameter_names['DeltaGs']:
            mle = getattr(map, name).value
            mean_cntr, var_cntr, std_cntr = bayes_mvs(getattr(mcmc, name).trace(), alpha=alpha)
            (center, (lower, upper)) = mean_cntr
            print("%64s : %5.1f [%5.1f, %5.1f] kT" % (mle, lower, upper))
