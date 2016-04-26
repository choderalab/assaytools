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
        DeltaG_prior='uniform', inner_filter_effect=True):
        """
        Parameters
        ----------
        solutions : dict of str : Solution
            `solutions[name]` is the Solution object corresponding to component `name`
        wells : autoprotocol.container.WellGroup
            Group of wells with `contents` and `measurements` properties defined
        receptor_name : str
            Name of receptor
        ligand_names : list of str
            Names of ligands
        DeltaG_prior : str, optional, default='uniform'
            Prior to use on DeltaG values for reaction ['uniform', 'chembl']
        inner_filter_effect : bool, optional, default=True
            If True, use primary and secondary inner filter effects.
            If 'primary', only use primary inner filter effect.
            If False, no inner filter effect is used.

        TODO:
        * For well.humanize(), also include plate name in case we may have well spanning multiple plates

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
        for solution in solutions.values():
            name = 'true concentration of %s' % solution.name
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

        self.all_species = [self.receptor_name] + self.ligand_names + self.complex_names

        # Construct priors dispensed volumes in each well
        volume_unit = 'liters' # volume unit used throughout
        self.parameter_names['dispensed_volumes'] = list()
        self.parameter_names['well_volumes'] = list()
        for well in wells:
            component_variables = list()
            component_coefficients = list()
            for component in well.properties['contents']:
                name = 'volume of %s dispensed into well %s' % (component, well.humanize()) # TODO: Use plate name as well in case we may have well spanning multiple plates
                (volume, error) = well.properties['contents'][component]
                mu = volume.m_as(volume_unit)
                sigma = error.m_as(volume_unit)
                volume_dispensed = pymc.Lognormal(name, mu=np.log(mu**2 / np.sqrt(sigma**2 + mu**2)), tau=np.sqrt(np.log(1.0 + (sigma/mu)**2))**(-2)) # L
                model[name] = volume_dispensed
                self.parameter_names['dispensed_volumes'].append(name)

                component_variables.append(volume_dispensed)
                component_coefficients.append(+1.0)

            # Total volume in well
            name = 'volume of well %s' % well.humanize() # TODO: Use plate name too in case wells span multiple plates
            model[name] = pymc.LinearCombination(name, component_coefficients, component_variables)
            self.parameter_names['well_volumes'].append(name)

            # Total concentration of each species in well
            # TODO: Set total concentration of
            for component in well.properties['contents']:
                solution = solutions[component]
                if solution.concentration.m_as(concentration_unit) == 0.0: continue # skip stocks
                species = solution.species
                solution_concentration = model['true concentration of %s' % solution.name]
                solution_volume = model['volume of %s dispensed into well %s' % (component, well.humanize())]
                total_volume = model['volume of well %s' % well.humanize()]
                @pymc.deterministic
                def total_concentration(solution_concentration=solution_concentration, solution_volume=solution_volume, total_volume=total_volume):
                    return solution_concentration * solution_volume / total_volume
                name = 'total concentration of %s in well %s' % (species, well.humanize())
                model[name] = concentration

            # TODO: Equilibrium concentration of each component in well
            # NOTE: Don't include reactions with components that are not present
            # NOTE: Also, don't need to compute concentrations if there are zero or one species in the well.

        #
        # Spectroscopic measurements
        #

        # Determine all wavelengths and detection technologies in use
        all_wavelengths = set()
        fluorescence_wavelength_pairs = set()
        absorbance = False
        fluorescence = False
        fluorescence_top = False
        fluorescence_bottom = False
        for well in wells:
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                absorbance = True
                for wavelength in measurements['absorbance'].keys():
                    all_wavelengths.add(wavelength)
            if 'fluorescence' in measurements:
                fluorescence = True
                for (excitation_wavelength, emission_wavelength, geometry) in measurements['fluorescence'].keys():
                    if geometry == 'top': fluorescence_top = True
                    if geometry == 'bottom': fluorescence_bottom = True
                    all_wavelengths.add(excitation_wavelength)
                    all_wavelengths.add(emission_wavelength)
                    fluorescence_wavelength_pairs = (excitation_wavelength, emission_wavelength)
        print("all wavelengths in use:")
        print all_wavelengths
        self.all_wavelengths = all_wavelengths

        # Priors for spectroscopic parameters of each component
        if inner_filter_effect in ['primary', True]:
            self.parameter_names['extinction coefficients'] = list()
            for species in self.all_species:
                for wavelength in all_wavelengths:
                    name = 'extinction coefficient of %s at wavelength %s' % (species, wavelength)
                    extinction_coefficient = pymc.Uniform(name, lower=0.0, upper=1000e3, value=70000.0) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm
                    model[name] = extinction_coefficient
                    self.parameter_names['extinction coefficients'].append(name)

        #
        # Absorbance
        #

        # Prior on absorbance measurement error
        if absorbance:
            model['log absorbance error'] = pymc.Uniform('log absorbance error', lower=-10, upper=0, value=np.log(0.01))
            model['absorbance error'] = pymc.Lambda('absorbance error', lambda log_sigma=model['log absorbance error'] : np.exp(log_sigma) )
            model['absorbance precision'] = pymc.Lambda('absorbance precision', lambda log_sigma=model['log absorbance error'] : np.exp(-2*log_sigma) )
            self.parameter_names['absorbance'] = ['log absorbance error', 'absorbance error', 'absorbance precision']
            # Prior on plate absorbance at each wavelength
            for wavelength in all_wavelengths:
                name = 'plate absorbance at wavelength %s' % wavelength
                model[name] = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.0)
                self.parameter_names['absorbance'].append(name)

        # Absorbance measurements (for wells that have them)
        for well in wells:
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                for wavelength in measurements['absorbance'].keys():
                    extinction_coefficients = list()
                    concentrations = list()
                    for species in self.all_species:
                        extinction_coefficients.append( model['extinction coefficient of %s at wavelength %s' % (species, wavelength)] )
                        concentrations.append( model['concentration of %s in well %s' % (species, well.humanize())] )
                    plate_absorbance = model['plate absorbance at wavelength %s' % wavelength]

                    @pymc.deterministic
                    def absorbance_model(concentrations=concentrations, extinction_coefficients=extinction_coefficients, plate_absorbance=plate_absorbance):
                        ec = 0.0
                        for species in concentrations.keys():
                            ec += extinction_coefficients[species] * concentrations[species]
                        absorbance = (1.0 - np.exp(-ec * path_length)) + plate_absorbance
                        return absorbance
                    name = 'computed absorbance of well %s at wavelength %s' % (well.humanize(), wavelength)
                    model[name] = absorbance_model

                    measured_absorbance = measurements['absorbance'][wavelength]
                    name = 'measured absorbance of well %s at wavelength %s' % (well.humanize(), wavelength) # TODO: Include plate name
                    model[name] = pymc.Normal(name, mu=absorbance_model, tau=model['absorbance precision'], observed=True, value=measured_absorbance)

        #
        # Fluorescence
        #

        # TODO: We can have quantum yield instead?

        # Fluorescence quantum yields
        self.parameter_names['fluorescence'] = list()
        if fluorescence:
            for (excitation_wavelength, emission_wavelength) in fluorescence_wavelength_pairs:
                name = 'plate background fluorescence for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)
                quantum_yield = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.1)
                model[name] = quantum_yield
                self.parameter_names['fluorescence'].append(name)

                for species in self.all_species:
                    name = 'quantum yield of %s for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)
                    quantum_yield = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.1)
                    model[name] = quantum_yield
                    self.parameter_names['fluorescence'].append(name)

            # Fluorescence intensities * gains
            # TODO: If multiple gains are in use, slave them together through this intensity times a fixed gain factor.
            if fluorescence_top:
                max_top_fluorescence_intensity = 1.0e8 # TODO: Determine maximum possible fluorescence intensity

                name = 'top fluorescence illumination intensity'
                top_fluorescence_intensity = pymc.Uniform(name, lower=0.0, upper=max_top_fluorescence_intensity)
                model[name] = top_fluorescence_intensity
                self.parameter_names['fluorescence'].append(name)

            if fluorescence_bottom:
                max_bottom_fluorescence_intensity = 1.0e8 # TODO: Determine maximum possible fluorescence intensity

                name = 'bottom fluorescence illumination intensity'
                bottom_fluorescence_intensity = pymc.Uniform(name, lower=0.0, upper=max_bottom_fluorescence_intensity)
                model[name] = bottom_fluorescence_intensity
                self.parameter_names['fluorescence'].append(name)

        def all_initial_species_in_well(well):
            """List all initial species in the well"""
            all_species = list()
            for component in well.properties['contents']:
                solution = solutions[component]
                if solution.concentration.m_as(concentration_unit) == 0.0: continue # skip stocks
                species = solution.species
                all_species.append(species)
            return all_species

        def all_equilibrium_species_in_well(well):
            """List all species in the well after equilibrium has been reached"""
            all_species = list()
            for component in well.properties['contents']:
                solution = solutions[component]
                if solution.concentration.m_as(concentration_unit) == 0.0: continue # skip stocks
                species = solution.species
                all_species.append(species)
            return all_species

        # Fluorescence measurements (for wells that have them)
        for well in wells:
            measurements = well.properties['measurements']
            if 'fluorescence' in measurements:
                for (excitation_wavelength, emission_wavelength, geometry) in measurements['fluorescence'].keys():
                    # Extract extinction coefficients and concentrations for all species in well and pack them into lists
                    excitation_extinction_coefficients = list()
                    emission_extinction_coefficients = list()
                    concentrations = list()
                    quantum_yields = list()
                    for species in all_species_in_well(well): # TODO: Extract all equilibrium species in well
                        quantum_yields.append( model['quantum yield of %s for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)] )
                        extinction_coefficients.append( model['extinction coefficient of %s at wavelength %s' % (species, excitation_wavelength)] )
                        extinction_coefficients.append( model['extinction coefficient of %s at wavelength %s' % (species, emission_wavelength)] )
                        concentrations.append( model['concentration of %s in well %s' % (species, well.humanize())] )
                    plate_fluorescence = model['plate background fluorescence for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)]
                    top_illumination_intensity = model['top fluorescence illumination intensity']
                    bottom_illumination_intensity = model['bottom fluorescence illumination intensity']

                    @pymc.deterministic
                    def fluorescence_model(concentrations=concentrations, quantum_yields=quantum_yields,
                        excitation_extinction_coefficients=extinction_coefficients, emission_extinction_coefficients=emission_extinction_coefficients,
                        plate_fluorescence=plate_fluorescence,
                        top_illumination_intensity=top_illumination_intensity, bottom_illumination_intensity=bottom_illumination_intensity, geometry=geometry):

                        if geometry == 'top':
                            intensity = top_illumination_intensity
                        elif geometry == 'bottom':
                            intensity = bottom_illumination_intensity
                        fluorescence = intensity * plate_fluorescence # background
                        for (quantum_yield, concentration, excitation_extinction_coefficient, emission_extinction_coefficient) in zip(quantum_yields, concentrations, excitation_extinction_coefficients, emission_extinction_coefficients):
                            fluorescence += intensity * quantum_yield * concentration
                            # TODO: Add inner filter effects
                        return fluorescence
                    name = 'computed %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, well.humanize(), excitation_wavelength, emission_wavelength)
                    model[name] = fluorescence

                    measured_fluorescence = measurements['fluorescence'][(excitation_wavelength, emission_wavelength, geometry)]
                    name = 'measured %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, well.humanize(), excitation_wavelength, emission_wavelength)
                    model[name] = pymc.Normal(name, mu=fluorescence_model, tau=model['%s fluorescence precision' % geometry], observed=True, value=measured_fluorescence)

        # Add to model.



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
        ncycles = 1
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

    def show_summary(self, mcmc, map_fit):
        """
        Show summary statistics of MCMC and MAP estimates.

        Parameters
        ----------
        mcmc : pymc.MCMC
           MCMC samples
        map_fit : pymc.MAP
           The MAP fit.

        TODO
        ----
        * Automatically determine appropriate number of decimal places from statistical uncertainty.
        * Automatically adjust concentration units (e.g. pM, nM, uM) depending on estimated affinity.

        """

        # Compute summary statistics
        alpha = 0.95 # confidence interval width
        from scipy.stats import bayes_mvs
        for name in self.parameter_names['DeltaGs']:
            mle = getattr(map_fit, name).value
            mean_cntr, var_cntr, std_cntr = bayes_mvs(getattr(mcmc, name).trace(), alpha=alpha)
            (center, (lower, upper)) = mean_cntr
            print("%64s : %5.1f [%5.1f, %5.1f] kT" % (mle, lower, upper))
