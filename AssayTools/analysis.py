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

        # Model solution concentrations
        for (name, solution) in solutions.iteritems():
            name += ' concentration'
            mu = solution.concentration.to('moles/liter').magnitude # M, gaussian mean
            sigma = solution.uncertainty.to('moles/liter').magnitude # M, gaussian mean
            if (mu == 0.0):
                true_concentration = 0.0
            else:
                true_concentration = pymc.Lognormal(name, mu=np.log(mu**2 / np.sqrt(sigma**2 + mu**2)), tau=np.sqrt(np.log(1.0 + (sigma/mu)**2))**(-2)) # protein concentration (M)
            model[name] = true_concentration

        # Construct reactions
        from bindingmodels import GeneralBindingModel
        reactions = list()
        for ligand_name in ligand_names:
            # Create the DeltaG prior.
            complex_name = receptor_name + ':' + ligand_name
            name = 'DeltaG ' + complex_name
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
            reactions.append(reaction)



        self.model = pymc.Model(model)

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

        # DEBUG
        ncycles = 5

        for cycle in range(ncycles):
            if (cycle+1)%5==0: print('MAP fitting cycle %d/%d' % (cycle+1, ncycles))
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

        # Compute summary statistics.
        DeltaG = map.DeltaG.value
        dDeltaG = mcmc.DeltaG.trace().std()
        Kd = np.exp(map.DeltaG.value)
        dKd = np.exp(mcmc.DeltaG.trace()).std()
        print "DeltaG = %.1f +- %.1f kT" % (DeltaG, dDeltaG)
        if (Kd < 1e-12):
            print "Kd = %.1f nM +- %.1f fM" % (Kd/1e-15, dKd/1e-15)
        elif (Kd < 1e-9):
            print "Kd = %.1f pM +- %.1f pM" % (Kd/1e-12, dKd/1e-12)
        elif (Kd < 1e-6):
            print "Kd = %.1f nM +- %.1f nM" % (Kd/1e-9, dKd/1e-9)
        elif (Kd < 1e-3):
            print "Kd = %.1f uM +- %.1f uM" % (Kd/1e-6, dKd/1e-6)
        elif (Kd < 1):
            print "Kd = %.1f mM +- %.1f mM" % (Kd/1e-3, dKd/1e-3)
        else:
            print "Kd = %.3e M +- %.3e M" % (Kd, dKd)
