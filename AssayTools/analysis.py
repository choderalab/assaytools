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
            `solutions[name]` is the Solution object corresponding to component `name`
        wells : autoprotocol.container.WellGroup
            Group of wells with `contents` and `measurements` properties defined
        receptor_name : str
            Name of receptor
        ligand_names : list of str
            Names of ligands
        DeltaG_prior : str, optional, default='uniform'
            Prior to use on DeltaG values for reaction ['uniform', 'chembl']

        TODO:
        * For well.humanize(), also include plate name in case we may have well spanning multiple plates

        """
        # Store data
        self.solutions = solutions
        self.wells = wells
        self.receptor_name = receptor_name
        self.ligand_names = list(ligand_names)

        # Set up internal data structures.
        self.model = dict() # the PyMC model; self.model[paramname] is the PyMC variable correspoding to 'paramname'
        self.parameter_names = dict() # dict to keep track of groups of related parameter names; self.parameter_names[groupname] is the list of PyMC variable names under 'groupname'

        # Construct different parts of the pymc model.
        self._create_solutions_model(solutions)
        self._create_competitive_binding_model(DeltaG_prior)
        self._create_dispensing_model()
        self._create_equilibrium_concentrations_model()
        self._create_extinction_coefficients_model()
        self._create_absorbance_model()
        self._create_fluorescence_model()

        # Create the PyMC Model object from the model dict
        self.model = pymc.Model(model)
        print('Model has %d stochastics and %d deterministics...' % (len(self.model.stochastics), len(self.model.deterministics)))

    def _create_solutions_model(self, solutions):
        """
        Create pymc model components for true concentrations of source receptor and ligand solutions.

        Populates the following fields:
        * parameter_names['concentrations'] : parameters associated with true concentrations of receptor and ligand solutions
        """
        self.parameter_names['solution concentrations'] = list()
        concentration_unit = 'moles/liter'
        for solution in solutions.values():
            name = 'true concentration of %s' % solution.name
            mu = solution.concentration.m_as(concentration_unit) # M, gaussian mean
            sigma = solution.uncertainty.m_as(concentration_unit) # M, gaussian mean
            if (mu == 0.0):
                true_concentration = 0.0
            else:
                true_concentration = pymc.Lognormal(name, mu=np.log(mu**2 / np.sqrt(sigma**2 + mu**2)), tau=np.sqrt(np.log(1.0 + (sigma/mu)**2))**(-2)) # protein concentration (M)
            self.model[name] = true_concentration
            self.parameter_names['solution concentrations'].append(name)

    def _create_competitive_binding_model(self, DeltaG_prior):
        """
        Create the binding free energy priors, binding reaction models, and list of all species whose concentrations will be tracked.

        Populates the following fields:
        * parameter_names['DeltaGs'] : parameters associated with DeltaG priors
        * reactions : reactions for GeneralBindingModel for competitive binding
        * complex_names : names of all complexes
        * all_species : names of all species (ligands, receptor, complexes)
        """
        self.parameter_names['DeltaGs'] = list()
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
            self.model[name] = DeltaG
            # Create the reaction
            reaction = (DeltaG, {complex_name : -1, receptor_name : +1, ligand_name : +1})
            # Append reactions.
            self.reactions.append(reaction)
            self.parameter_names['DeltaGs'].append(name)

        # Create a list of all species that may be present in assay plate
        self.all_species = [self.receptor_name] + self.ligand_names + self.complex_names

    def _create_dispensing_model(self):
        """
        Create nuisance parameters for dispensed volumes and actual concentrations of all species in each well.

        Populates the following fields:
        * parameter_names['dispensed_volumes'] : actual volumes dispensed into each well
        * parameter_names['well_volumes'] : actual well total volumes
        *
        """
        # Construct priors dispensed volumes in each well
        volume_unit = 'liters' # volume unit used throughout
        self.parameter_names['dispensed_volumes'] = list()
        self.parameter_names['well_volumes'] = list()
        for well in self.wells:
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
            self.model[name] = pymc.LinearCombination(name, component_coefficients, component_variables)
            self.parameter_names['well_volumes'].append(name)

            # Total concentration of each species in well
            self.parameters_names['well concentrations'] = dict()
            for component in well.properties['contents']:
                solution = self.solutions[component]
                if solution.concentration.m_as(concentration_unit) == 0.0: continue # skip stocks
                species = solution.species
                solution_concentration = self.model['true concentration of %s' % solution.name]
                solution_volume = model['volume of %s dispensed into well %s' % (component, well.humanize())]
                total_volume = model['volume of well %s' % well.humanize()]

                # TODO: Work with log volumes and log concentrations directly instead of lognormal distributions of normal volumes and concentrations?
                @pymc.deterministic
                def total_concentration(solution_concentration=solution_concentration, solution_volume=solution_volume, total_volume=total_volume):
                    return solution_concentration * solution_volume / total_volume
                name = 'total concentration of %s in well %s' % (species, well.humanize())
                self.model[name] = concentration
                self.parameter_names['well concentrations'].append(name)

                @pymc.deterministic
                def log_total_concentration(concentration=concentration):
                    return np.log(concentration)
                name = 'log total concentration of %s in well %s' % (species, well.humanize())
                self.model[name] = concentration
                self.parameter_names['well concentrations'].append(name)


    def all_initial_species_in_well(self, well):
        """
        List all initial species in the well


        """
        all_species = set()
        for component in well.properties['contents']:
            solution = self.solutions[component]
            if solution.concentration.m_as(concentration_unit) == 0.0: continue # skip stocks
            species = solution.species
            all_species.add(species)
        return all_species

    def all_equilibrium_species_in_well(self, well):
        """
        List all species present in the well after equilibrium has been reached

        """
        initial_species = self.all_initial_species_in_well(self, well)
        all_species = set(initial_species)

        previous_nspecies = len(all_species)
        while (current_nspecies != previous_nspecies):
            for reaction in self.reactions:
                # Decompose reaction into product and reactant sets
                (DeltaG, stoichiometry) = reaction
                reactants = set([ species for species in stoichiometry.keys() if stoichiometry[species] < 0 ])
                products = set([ species for species in stoichiometry.keys() if stoichiometry[species] > 0 ])
                # If reactants are present, add products; and vice-versa
                if reactants.issubset(all_species): all_species.update(products)
                if products.issubset(all_species): all_species.update(reactants)
                # Update count of number of species
                previous_nspecies = current_nspecies
                current_nspecies = len(all_species)

        return all_species

    def binding_reactions_for_well(self, well):
        """
        Determine subset of binding reactions relevant to a well.

        """
        # Determine all relevant species
        all_species = self.all_initial_species_in_well(well)
        reactions = list()
        for reaction in self.reactions:
            # Decompose reaction into product and reactant sets
            (DeltaG, stoichiometry) = reaction
            # Reactants and products should all be present in equilibrium concentration
            if set(stoichiometry.keys()).issubset(all_species)
                reactions.append(reaction)

        return reactions

    def log_total_concentrations_for_well(self, well):
        """
        Return list of pymc variables of log total concentrations of each species in the well.
        """
        log_total_concentations = dict()
        for component in well.properties['contents']:
            solution = self.solutions[component]
            if solution.concentration.m_as(concentration_unit) == 0.0: continue # skip stocks
            species = solution.species
            name = 'log total concentration of %s in well %s' % (species, well.humanize())
            log_total_concentration = self.model[name]
            log_total_concentrations[species] = log_total_concentration

        return log_total_concentrations

    def _create_equilibrium_concentrations_model(self):
        """
        Create model for equilibrium concentration of each species in each well.

        """
        from bindingmodels import GeneralBindingModel

        self.parameter_names['well concentrations'] = list()
        for well in self.wells():
            # Determine list of all species that can be in this well.
            all_species = self.all_equilibrium_species_in_well(well)

            # Determine relevant list of reactions for this well.
            reactions = self.binding_reactions_for_well(well)

            # Build list of total concentrations of each species in well.
            log_total_concentrations = self.log_total_concentrations_for_well(well)
            # TODO: Form conservation equations
            conservation_equations = [ (-6, {'RL' : +1, 'R' : +1}), (-6, {'RL' : +1, 'L' : +1}) ]

            # Compute equilibrium concentration of each component in well
            # TODO: Compute separately for each species.
            @deterministic
            def log_equilibrium_concentrations(log_total_concentrations=log_total_concentrations):
                return GeneralBindingModel.equilibrium_concentrations(reactions, conservation_equations)
            name = 'log equilibrium concentration of species %s in well %s' % (species, well.humanize())
            self.model[name] = log_equilibrium_concentrations


    def _create_extinction_coefficients_model(self, inner_filter_effect):
        """
        Determine all spectroscopic wavelengths in use and create model of extinction coefficients

        Parameters
        ----------
        inner_filter_effect : bool, optional, default=True
            If True, use primary and secondary inner filter effects.
            If 'primary', only use primary inner filter effect.
            If False, no inner filter effect is used.

        Populates the following fields:
        * parameter_names['extinction coefficients'] : all extinction coefficients
        * all_wavelengths : list of all wavelengths in use
        * absorbance : True if absorbance in use
        * fluorescence : True if fluorescence in use

        """
        #
        # Spectroscopic measurements
        #

        MAX_EXTINCTION_COEFFICIENT = 100e3

        # Determine all wavelengths and detection technologies in use
        self.all_wavelengths = set()
        for well in self.wells:
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                for wavelength in measurements['absorbance'].keys():
                    self.all_wavelengths.add(wavelength)
            if 'fluorescence' in measurements:
                for (excitation_wavelength, emission_wavelength, geometry) in measurements['fluorescence'].keys():
                    self.all_wavelengths.add(excitation_wavelength)
                    self.all_wavelengths.add(emission_wavelength)
        print("all wavelengths in use:")
        print self.all_wavelengths

        if (len(self.all_wavelengths) > 0):
            # Priors for spectroscopic parameters of each component
            if inner_filter_effect in ['primary', True]:
                self.parameter_names['extinction coefficients'] = list()
                for species in self.all_species:
                    for wavelength in self.all_wavelengths:
                        name = 'extinction coefficient of %s at wavelength %s' % (species, wavelength)
                        extinction_coefficient = pymc.Uniform(name, lower=0.0, upper=MAX_EXTINCTION_COEFFICIENT) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm
                        self.model[name] = extinction_coefficient
                        self.parameter_names['extinction coefficients'].append(name)

    def _create_absorbance_model(self):
        """
        Absorbance measurements.

        Populates the following fields

        """
        # Determine if absorbance is in use
        absorbance = False
        for well in self.wells:
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                absorbance = True

        # Prior on absorbance measurement error
        if absorbance:
            print('Absorbance measurements are available')
            self.model['log absorbance error'] = pymc.Uniform('log absorbance error', lower=-10, upper=0, value=np.log(0.01))
            self.model['absorbance error'] = pymc.Lambda('absorbance error', lambda log_sigma=self.model['log absorbance error'] : np.exp(log_sigma) )
            self.model['absorbance precision'] = pymc.Lambda('absorbance precision', lambda log_sigma=self.model['log absorbance error'] : np.exp(-2*log_sigma) )
            self.parameter_names['absorbance'] = ['log absorbance error', 'absorbance error', 'absorbance precision']
            # Prior on plate absorbance at each wavelength
            for wavelength in self.all_wavelengths:
                name = 'plate absorbance at wavelength %s' % wavelength
                self.model[name] = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.0)
                self.parameter_names['absorbance'].append(name)

        # Absorbance measurements (for wells that have them)
        for well in self.wells:
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                for wavelength in measurements['absorbance'].keys():
                    # Compile all parameters needed for absorbance model
                    extinction_coefficients = list()
                    concentrations = list()
                    for species in self.all_species:
                        extinction_coefficients.append( self.model['extinction coefficient of %s at wavelength %s' % (species, wavelength)] )
                        concentrations.append( self.model['concentration of %s in well %s' % (species, well.humanize())] )
                    plate_absorbance = self.model['plate absorbance at wavelength %s' % wavelength]

                    # Add computed absorbance model
                    @pymc.deterministic
                    def absorbance_model(concentrations=concentrations, extinction_coefficients=extinction_coefficients, plate_absorbance=plate_absorbance):
                        ec = 0.0
                        for species in concentrations.keys():
                            ec += extinction_coefficients[species] * concentrations[species]
                        absorbance = (1.0 - np.exp(-ec * path_length)) + plate_absorbance
                        return absorbance
                    name = 'computed absorbance of well %s at wavelength %s' % (well.humanize(), wavelength)
                    self.model[name] = absorbance_model
                    self.parameters['absorbance'].append(name)

                    # Add measured absorbance model
                    measured_absorbance = measurements['absorbance'][wavelength]
                    name = 'measured absorbance of well %s at wavelength %s' % (well.humanize(), wavelength) # TODO: Include plate name
                    self.model[name] = pymc.Normal(name, mu=absorbance_model, tau=self.model['absorbance precision'], observed=True, value=measured_absorbance)
                    self.parameters['absorbance'].append(name)

    def _create_fluorescence_model(self, inner_filter_effect):
        """
        Create model for fluorescence measurements.

        """

        # Determine if fluorescence is in use
        fluorescence = False
        fluorescence_top = False
        fluorescence_bottom = False
        fluorescence_wavelength_pairs = set()
        for well in self.wells:
            measurements = well.properties['measurements']
            if 'fluorescence' in measurements:
                fluorescence = True
                for (excitation_wavelength, emission_wavelength, geometry) in measurements['fluorescence'].keys():
                    if geometry == 'top': fluorescence_top = True
                    if geometry == 'bottom': fluorescence_bottom = True
                    fluorescence_wavelength_pairs = (excitation_wavelength, emission_wavelength)

        if fluorescence:
            self.parameter_names['fluorescence'] = list()
            for (excitation_wavelength, emission_wavelength) in fluorescence_wavelength_pairs:
                name = 'plate background fluorescence for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)
                quantum_yield = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.1)
                self.model[name] = quantum_yield
                self.parameter_names['fluorescence'].append(name)

                for species in self.all_species:
                    name = 'quantum yield of %s for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)
                    quantum_yield = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.1)
                    self.model[name] = quantum_yield
                    self.parameter_names['fluorescence'].append(name)

            # Fluorescence intensities * gains
            # TODO: If multiple gains are in use, slave them together through this intensity times a fixed gain factor.
            if fluorescence_top:
                max_top_fluorescence_intensity = 1.0e8 # TODO: Determine maximum possible fluorescence intensity

                name = 'top fluorescence illumination intensity'
                top_fluorescence_intensity = pymc.Uniform(name, lower=0.0, upper=max_top_fluorescence_intensity)
                self.model[name] = top_fluorescence_intensity
                self.parameter_names['fluorescence'].append(name)

            if fluorescence_bottom:
                max_bottom_fluorescence_intensity = 1.0e8 # TODO: Determine maximum possible fluorescence intensity

                name = 'bottom fluorescence illumination intensity'
                bottom_fluorescence_intensity = pymc.Uniform(name, lower=0.0, upper=max_bottom_fluorescence_intensity)
                self.model[name] = bottom_fluorescence_intensity
                self.parameter_names['fluorescence'].append(name)

        # Fluorescence measurements (for wells that have them)
        for well in self.wells:
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

                        intensity = (geometry == 'top') * top_illumination_intensity + (geometry == 'bottom') * bottom_illumination_intensity # select appropriate illumination intensity
                        fluorescence = intensity * plate_fluorescence # background
                        for (quantum_yield, concentration, excitation_extinction_coefficient, emission_extinction_coefficient) in zip(quantum_yields, concentrations, excitation_extinction_coefficients, emission_extinction_coefficients):
                            # TODO: Compute inner filter effects
                            inner_fileter_effect = 1.0
                            if self.inner_fileter_effect


                            fluorescence += intensity * quantum_yield * concentration * inner_filter_effect


                        return fluorescence
                    name = 'computed %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, well.humanize(), excitation_wavelength, emission_wavelength)
                    model[name] = fluorescence

                    measured_fluorescence = measurements['fluorescence'][(excitation_wavelength, emission_wavelength, geometry)]
                    name = 'measured %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, well.humanize(), excitation_wavelength, emission_wavelength)
                    model[name] = pymc.Normal(name, mu=fluorescence_model, tau=model['%s fluorescence precision' % geometry], observed=True, value=measured_fluorescence)


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

    def geneate_plots(self, mcmc, map=None, pdf_filename=None):
        """
        Generate interactive or PDF plots from MCMC trace.

        Parameters
        ----------
        mcmc : pymc.MCMC
           MCMC samples to plot
        map : pymc.MAP, optional, default=None
           Plot the maximum a posteriori (MAP) estimate if provided.
        pdf_filename : str, optional, default=None
           If specified, generate a PDF containing plots.

        """
        raise Exception('Not implemented yet.')
