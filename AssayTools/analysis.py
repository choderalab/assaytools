"""
Classes for the analysis of fluorescence assay data.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy as np
import pymc
from scipy.misc import logsumexp
from autoprotocol.unit import Unit

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

def LogNormalWrapper(name, mean, stddev):
    """
    Create a PyMC Normal stochastic, automatically converting parameters to be appropriate for a `LogNormal` distribution.
    Note that the resulting distribution is Normal, not LogNormal.
    This is appropriate if we want to describe the log of a volume or a fluorescence reading.

    Parameters
    ----------
    mean : float
       Mean of exp(X), where X is lognormal variate.
    stddev : float
       Standard deviation of exp(X), where X is lognormal variate.

    Returns
    -------
    stochastic : pymc.Normal
       Normal stochastic for random variable X

    """
    # Compute parameters of lognormal distribution
    mu = np.log(mean**2 / np.sqrt(stddev**2 + mean**2))
    tau = np.sqrt(np.log(1.0 + (stddev/mean)**2))**(-2)
    stochastic = pymc.Normal(name, mu=mu, tau=tau)
    return stochastic

def wellname(well):
    """
    Concatenate container and well names.

    Parameters
    ----------
    well : autoprotocol.container.Well
       The well

    Returns
    -------
    name : str
       [container name] [humanized well name]
    """
    return well.container.name + ' ' + well.humanize()

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

        """
        # Store data
        self.solutions = solutions
        self.wells = wells
        self.receptor_name = receptor_name
        self.ligand_names = list(ligand_names)

        # Set up internal data structures.
        self.model = dict() # the PyMC model; self.model[paramname] is the PyMC variable correspoding to 'paramname'
        self.parameter_names = dict() # dict to keep track of groups of related parameter names; self.parameter_names[groupname] is the list of PyMC variable names under 'groupname'

        # Construct components of the pymc model.
        self._create_solutions_model(solutions)
        self._create_dispensing_model()
        self._create_competitive_binding_model(receptor_name, DeltaG_prior)
        self._create_equilibrium_concentrations_model()
        self._create_extinction_coefficients_model()
        self._create_absorbance_model()
        self._create_fluorescence_model()

        # Create the PyMC Model object from the dictionary of pymc stochastics and deterministics.
        self.model = pymc.Model(self.model)
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
            if solution.species is None:
                continue # skip buffers or pure solvents
            name = 'log concentration of %s' % solution.shortname
            self.model[name] = LogNormalWrapper(name, mean=solution.concentration.m_as(concentration_unit), stddev=solution.uncertainty.m_as(concentration_unit))
            self.parameter_names['solution concentrations'].append(name)

    def _create_dispensing_model(self):
        """
        Create nuisance parameters for dispensed volumes and actual concentrations of all species in each well.

        Populates the following fields:
        * parameter_names['dispensed_volumes'] : actual volumes dispensed into each well
        * parameter_names['well_volumes'] : actual well total volumes
        *
        """
        volume_unit = 'liters' # volume unit used throughout
        self.parameter_names['dispensed volumes'] = list() # pymc variable names associated with dispensed volumes
        self.parameter_names['well volumes'] = list() # pymc variable names associated with well volumes
        self.parameter_names['well concentrations'] = list()

        for well in self.wells:
            # Volumes dispensed into each well
            log_volumes = list() # log volumes are in Liters
            for component in well.properties['contents']:
                name = 'log volume of %s dispensed into well %s' % (component, wellname(well))
                (volume, error) = well.properties['contents'][component]
                log_volume_dispensed = LogNormalWrapper(name, mean=volume.m_as(volume_unit), stddev=error.m_as(volume_unit))
                self.model[name] = log_volume_dispensed
                self.parameter_names['dispensed volumes'].append(name)
                log_volumes.append(log_volume_dispensed)

            # Total volume in well
            name = 'log volume of well %s' % wellname(well)
            @pymc.deterministic(name=name)
            def log_total_volume(log_volumes=log_volumes):
                return logsumexp(log_volumes)
            self.model[name] = log_total_volume
            self.parameter_names['well volumes'].append(name)

            # Total concentration of all species involving each component in well
            # TODO: In future, can we simply calculate initial concentrations of each species and use that in GeneralBindingModel instead?
            for component in well.properties['contents']:
                solution = self.solutions[component]
                if solution.species is None:
                    continue # skip buffers or pure solvents
                species = solution.species
                log_solution_concentration = self.model['log concentration of %s' % solution.shortname] # log concentrations are in molar
                log_solution_volume = self.model['log volume of %s dispensed into well %s' % (component, wellname(well))]
                log_total_volume = self.model['log volume of well %s' % wellname(well)]

                name = 'log total concentration of %s in well %s' % (species, wellname(well))
                @pymc.deterministic(name=name)
                def log_total_concentration(log_solution_concentration=log_solution_concentration, log_solution_volume=log_solution_volume, log_total_volume=log_total_volume):
                    return log_solution_concentration + log_solution_volume - log_total_volume
                self.model[name] = log_total_concentration
                self.parameter_names['well concentrations'].append(name)

    def _all_initial_species_in_well(self, well):
        """
        List all species initially added to the well, before association reactions occur.

        Parameters
        ----------
        well : autoprotocol.container.Well
           The well for which all species are to be listed.

        Returns
        -------
        all_species : list of str
           The names of all initial species added to the well

        """
        all_species = set()
        for component in well.properties['contents']:
            solution = self.solutions[component]
            if solution.species is None:
                continue # skip buffers and pure solvents
            species = solution.species
            all_species.add(species)
        return all_species

    def _binding_reactions_for_well(self, well):
        """
        Determine subset of binding reactions relevant to a well, as determined by
        which species were initially added to the well.

        Parameters
        ----------
        well : autoprotocol.container.Well
           The well for which all binding reactions are to be compiled.

        Returns
        -------
        reactions : list
           The subset of binding reactions from self.reactions that can occur in this well

        """
        # Determine all relevant species
        all_species = self._all_equilibrium_species_in_well(well)
        reactions = list()
        for reaction in self.reactions:
            # Decompose reaction into product and reactant sets
            (DeltaG, stoichiometry) = reaction
            # Reactants and products should all be present in equilibrium concentration
            if set(stoichiometry.keys()).issubset(all_species):
                reactions.append(reaction)

        return reactions

    def _log_total_concentrations_for_well(self, well):
        """
        Return list of pymc variables of log total concentrations of each species in the well.

        """
        log_total_concentations = dict()
        for component in well.properties['contents']:
            solution = self.solutions[component]
            if solution.species is None:
                continue # skip buffers and pure solvents
            species = solution.species
            name = 'log total concentration of %s in well %s' % (species, wellname(well))
            log_total_concentration = self.model[name]
            log_total_concentrations[species] = log_total_concentration

        return log_total_concentrations

    def _create_competitive_binding_model(self, receptor_name, DeltaG_prior):
        """
        Create the binding free energy priors, binding reaction models, and list of all species whose concentrations will be tracked.

        Populates the following fields:
        * parameter_names['DeltaGs'] : parameters associated with DeltaG priors
        * reactions : reactions for GeneralBindingModel with DeltaG parameters substituted for free energies
          formatted like [ ('DeltaG1', {'RL': -1, 'R' : +1, 'L' : +1}), ('DeltaG2', {'RP' : -1, 'R' : +1, 'P' : +1}) ]
        * conservation_equations : conservation relationships for GeneralBindingModel with species names substituted with log concentrations
          formatted like [ ('receptor', {'receptor' : +1, 'receptor:ligand' : +1}), ('ligand' : {'receptor:ligand' : +1, 'ligand' : +1}) ]
        * complex_names : names of all complexes
        * all_species : names of all species (ligands, receptor, complexes)
        """
        self.parameter_names['DeltaGs'] = list()
        self.complex_names = list()
        self.reactions = list() # reactions for GeneralBindingModel with pymc parameters substituted for free energies
        self.conservation_equations = list() # list of conservation relationships for GeneralBindingModel with pymc parameters substituted for log concentrations
        receptor_conservation_equation = { receptor_name : +1 } # Begin to populate receptor conservation equation
        for ligand_name in self.ligand_names:
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
            self.parameter_names['DeltaGs'].append(name)
            # Create the reaction for GeneralBindingModel
            self.reactions.append( (DeltaG, {complex_name : -1, receptor_name : +1, ligand_name : +1}) )
            # Create the conservation equation for GeneralBindingModel
            self.conservation_equations.append( (ligand_name, {ligand_name : +1, complex_name : +1}) )
            # Keep track of receptor conservation.
            receptor_conservation_equation[complex_name] = +1
        # Create the receptor conservation equation for GeneralBindingModel
        self.conservation_equations.append( (receptor_name, receptor_conservation_equation) )

        # Create a list of all species that may be present in assay plate
        self.all_species = [self.receptor_name] + self.ligand_names + self.complex_names

    def _conservation_equations_for_well(self, well):
        """
        Return list of conservation equations for each well.

        """
        # Get pymc variables corresponding to total concentrations of each initial species.
        conservation_equations = list()
        initial_species_in_well =  self._all_initial_species_in_well(well)
        for conservation_equation in self.conservation_equations:
            (species, stoichiometry) = conservation_equation
            if species in initial_species_in_well:
                name = 'log total concentration of %s in well %s' % (species, wellname(well))
                log_total_concentration = self.model[name]
                conservation_equations.append( (log_total_concentration, stoichiometry) )

        return conservation_equations

    def _all_equilibrium_species_in_well(self, well):
        """
        List all species present in the well after equilibrium has been reached.

        Parameters
        ----------
        well : autoprotocol.container.Well
           The well for which all species are to be listed.

        Returns
        -------
        all_species : list of str
           The names of all possible species that could be formed after equilibrium has been reached

        """
        initial_species = self._all_initial_species_in_well(well)
        all_species = set(initial_species)

        previous_nspecies = 0
        current_nspecies = len(all_species)
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

    def _create_equilibrium_concentrations_model(self):
        """
        Create model for equilibrium concentration of each species in each well.
        This function is agnostic to the exact binding model in use.

        """
        from bindingmodels import GeneralBindingModel

        self.parameter_names['well concentrations'] = list()
        for well in self.wells:
            # Determine list of all species that can be in this well
            all_initial_species     = self._all_initial_species_in_well(well)
            all_equilibrium_species = self._all_equilibrium_species_in_well(well)

            # Empty wells don't need models.
            if len(all_equilibrium_species) == 0:
                continue

            # Determine relevant list of binding reactions for this well.
            reactions = self._binding_reactions_for_well(well)
            if len(reactions)==0:
                # Initial concentrations will not change.
                for species in all_initial_species:
                    name = 'log concentration of %s in well %s' % (species, wellname(well))
                    parent_name = 'log total concentration of %s in well %s' % (species, wellname(well))
                    log_total_concentration = self.model[parent_name]
                    @pymc.deterministic(name=name)
                    def log_equilibrium_concentration(log_total_concentration=log_total_concentration):
                        return log_total_concentration
                    self.model[name] = log_equilibrium_concentration
                    self.parameter_names['well concentrations'].append(name)
                continue

            # Build list of conservation equations for this well
            conservation_equations = self._conservation_equations_for_well(well)
            if len(conservation_equations)==0:
                msg = 'No conservation equations for well %s\n' % wellname(well)
                msg += str(all_initial_species) + '\n'
                msg += str(all_equilibrium_species) + '\n'
                raise Exception(msg)

            # Compute equilibrium concentration of each component in well
            name = 'log equilibrium concentration of all species in well %s' % wellname(well)
            @pymc.deterministic(name=name)
            def log_equilibrium_concentrations(reactions=reactions, conservation_equations=conservation_equations):
                if len(reactions)==0 or len(conservation_equations)==0:
                    raise Exception(reactions, conservation_equations)
                return GeneralBindingModel.equilibrium_concentrations(reactions, conservation_equations)
            self.model[name] = log_equilibrium_concentrations

            # Separate out individual concentration components
            for species in all_equilibrium_species:
                name = 'log concentration of %s in well %s' % (species, wellname(well))
                @pymc.deterministic(name=name)
                def log_equilibrium_concentration(log_equilibrium_concentrations=log_equilibrium_concentrations):
                    return log_equilibrium_concentrations[species]
                self.model[name] = log_equilibrium_concentration
                self.parameter_names['well concentrations'].append(name)

    def _create_extinction_coefficients_model(self):
        """
        Determine all spectroscopic wavelengths in use and create model of extinction coefficients

        Populates the following fields:
        * parameter_names['extinction coefficients'] : all extinction coefficients
        * all_wavelengths : list of all wavelengths in use
        * absorbance : True if absorbance in use
        * fluorescence : True if fluorescence in use

        """
        #
        # Spectroscopic measurements
        #

        # TODO: Switch to log extinction coefficients and uniform prior in log extinction coefficients
        MAX_EXTINCTION_COEFFICIENT = Unit(100e3, '1/(moles/liter)/centimeter') # maximum physically reasonable extinction coefficient

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
        Amax = 0.0 # max absorbance
        for well in self.wells:
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                absorbance = True
                for wavelength in measurements['absorbance'].keys():
                    measured_absorbance = measurements['absorbance'][wavelength]
                    Amax = max(measured_absorbance, Amax)

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
                self.model[name] = pymc.Uniform(name, lower=0.0, upper=Amax, value=0.0)
                self.parameter_names['absorbance'].append(name)

        # Absorbance measurements (for wells that have them)
        for well in self.wells:
            # Compute path length
            log_well_volume = self.model['log volume of well %s' % wellname(well)]
            all_equilibrium_species_in_well = self._all_equilibrium_species_in_well(well)
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                for wavelength in measurements['absorbance'].keys():
                    # Compile all parameters needed for absorbance model
                    extinction_coefficients = list()
                    log_concentrations = list()
                    for species in all_equilibrium_species_in_well:
                        extinction_coefficients.append( self.model['extinction coefficient of %s at wavelength %s' % (species, wavelength)] )
                        log_concentrations.append( self.model['log concentration of %s in well %s' % (species, wellname(well))] )
                    plate_absorbance = self.model['plate absorbance at wavelength %s' % wavelength]

                    # Add computed absorbance model
                    name = 'computed absorbance of well %s at wavelength %s' % (wellname(well), wavelength)
                    @pymc.deterministic(name=name)
                    def absorbance_model(log_concentrations=log_concentrations, extinction_coefficients=extinction_coefficients, plate_absorbance=plate_absorbance, log_well_volume=log_well_volume):
                        well_volume = Unit(np.exp(log_well_volume), 'liters')
                        well_area = well.container.container_type.well_area
                        path_length = well_volume / well_area

                        absorbance = 0.0
                        for (extinction_coefficient, log_concentration) in zip(extinction_coefficients, log_concentrations):
                            extinction_coefficient = Unit(extinction_coefficient, '1/(moles/liter)/centimeter')
                            concentration = Unit(np.exp(log_concentration), 'moles/liter')
                            absorbance += extinction_coefficient * path_length * concentration
                        absorbance += plate_absorbance
                        return absorbance
                    self.model[name] = absorbance_model
                    self.parameter_names['absorbance'].append(name)

                    # Add measured absorbance model
                    measured_absorbance = measurements['absorbance'][wavelength]
                    name = 'measured absorbance of well %s at wavelength %s' % (wellname(well), wavelength) # TODO: Include plate name
                    self.model[name] = pymc.Normal(name, mu=absorbance_model, tau=self.model['absorbance precision'], observed=True, value=measured_absorbance)
                    self.parameter_names['absorbance'].append(name)

    def _create_fluorescence_model(self):
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
                    fluorescence_wavelength_pairs.add( (excitation_wavelength, emission_wavelength) )

        if fluorescence:
            self.parameter_names['fluorescence'] = list()
            for (excitation_wavelength, emission_wavelength) in fluorescence_wavelength_pairs:
                name = 'plate background fluorescence for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)
                quantum_yield = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.1)
                self.model[name] = quantum_yield
                self.parameter_names['fluorescence'].append(name)

                # TODO: If we had an estimate of quantum yield of some reference species, we could restraint these much better.
                for species in self.all_species:
                    name = 'quantum yield of %s for fluorescence excitation at %s and emission at %s' % (species, excitation_wavelength, emission_wavelength)
                    quantum_yield = pymc.Uniform(name, lower=0.0, upper=1.0, value=0.1)
                    self.model[name] = quantum_yield
                    self.parameter_names['fluorescence'].append(name)

            # Fluorescence intensities * gains
            # TODO: If multiple gains are in use, slave them together through this intensity times a fixed gain factor.
            if fluorescence_top:
                max_top_fluorescence_intensity = 1.0e8 # TODO: Determine maximum possible fluorescence intensity
                name = 'top fluorescence illumination intensity'
                self.model[name] = pymc.Uniform(name, lower=0.0, upper=max_top_fluorescence_intensity)
                self.parameter_names['fluorescence'].append(name)

                name = 'top fluorescence log uncertainty'
                self.model[name] = pymc.Uniform(name, lower=-10, upper=np.log(max_top_fluorescence_intensity), value=np.log(5))
                self.parameter_names['fluorescence'].append(name)

                name = 'top fluorescence precision'
                self.model[name] = pymc.Lambda(name, lambda log_sigma=self.model['top fluorescence log uncertainty'] : np.exp(-2*log_sigma))
                self.parameter_names['fluorescence'].append(name)

            if fluorescence_bottom:
                max_bottom_fluorescence_intensity = 1.0e8 # TODO: Determine maximum possible fluorescence intensity
                name = 'bottom fluorescence illumination intensity'
                self.model[name] = pymc.Uniform(name, lower=0.0, upper=max_bottom_fluorescence_intensity)
                self.parameter_names['fluorescence'].append(name)

                name = 'bottom fluorescence log uncertainty'
                self.model[name] = pymc.Uniform(name, lower=-10, upper=np.log(max_bottom_fluorescence_intensity), value=np.log(5))
                self.parameter_names['fluorescence'].append(name)

                name = 'bottom fluorescence precision'
                self.model[name] = pymc.Lambda(name, lambda log_sigma=self.model['bottom fluorescence log uncertainty'] : np.exp(-2*log_sigma))
                self.parameter_names['fluorescence'].append(name)

        # Fluorescence measurements (for wells that have them)
        for well in self.wells:
            log_well_volume = self.model['log volume of well %s' % wellname(well)]
            measurements = well.properties['measurements']
            if 'fluorescence' in measurements:
                for (excitation_wavelength, emission_wavelength, geometry) in measurements['fluorescence'].keys():
                    # Extract extinction coefficients and concentrations for all species in well and pack them into lists
                    quantum_yields = list()
                    excitation_extinction_coefficients = list()
                    emission_extinction_coefficients = list()
                    log_concentrations = list()
                    for species in self._all_equilibrium_species_in_well(well):
                        quantum_yields.append( self.model['quantum yield of %s for fluorescence excitation at %s and emission at %s' % (species, excitation_wavelength, emission_wavelength)] )
                        excitation_extinction_coefficients.append( self.model['extinction coefficient of %s at wavelength %s' % (species, excitation_wavelength)] )
                        emission_extinction_coefficients.append( self.model['extinction coefficient of %s at wavelength %s' % (species, emission_wavelength)] )
                        log_concentrations.append( self.model['log concentration of %s in well %s' % (species, wellname(well))] )
                    plate_fluorescence = self.model['plate background fluorescence for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)]
                    top_illumination_intensity = self.model['top fluorescence illumination intensity']
                    bottom_illumination_intensity = self.model['bottom fluorescence illumination intensity']

                    name = 'computed %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, wellname(well), excitation_wavelength, emission_wavelength)
                    @pymc.deterministic(name=name)
                    def fluorescence_model(log_well_volume=log_well_volume,
                        log_concentrations=log_concentrations, quantum_yields=quantum_yields,
                        excitation_extinction_coefficients=excitation_extinction_coefficients, emission_extinction_coefficients=emission_extinction_coefficients,
                        plate_fluorescence=plate_fluorescence,
                        top_illumination_intensity=top_illumination_intensity, bottom_illumination_intensity=bottom_illumination_intensity, geometry=geometry):

                        # Compute path length
                        well_volume = Unit(np.exp(log_well_volume), 'liters')
                        well_area = well.container.container_type.well_area
                        path_length = (well_volume / well_area).to('centimeter')

                        # Calculate attenuation due to inner filter effects
                        # TODO: We may need to fix some scaling factors in the extinction coefficient to go between log10 and ln-based absorbance/transmission.
                        inner_filter_effect_scaling = 1.0 # scaling factor applied from inner filter effect
                        ELC_excitation = 0.0 # sum of (extinction_coefficient * path_length * concentration) for all species at excitation wavelength
                        ELC_emission   = 0.0 # sum of (extinction_coefficient * path_length * concentration) for all species at emission wavelength
                        for (log_concentration, excitation_extinction_coefficient, emission_extinction_coefficient) in zip(log_concentrations, excitation_extinction_coefficients, emission_extinction_coefficients):
                            # Accumulate inner filter effects
                            excitation_extinction_coefficient = Unit(excitation_extinction_coefficient, '1/(moles/liter)/centimeter')
                            emission_extinction_coefficient   = Unit(emission_extinction_coefficient, '1/(moles/liter)/centimeter')
                            concentration = Unit(np.exp(log_concentration), 'moles/liter')

                            ELC_excitation += excitation_extinction_coefficient * path_length * concentration
                            ELC_emission   += emission_extinction_coefficient   * path_length * concentration
                        if geometry == 'top':
                            alpha = (ELC_excitation + ELC_emission) * np.log(10)
                        elif geometry == 'bottom':
                            alpha = (ELC_excitation - ELC_emission) * np.log(10)
                        inner_filter_effect_attenuation = (1 - np.exp(-alpha)) / alpha
                        # Handle alpha -> 0 case explicitly.
                        if np.abs(alpha) < 0.01:
                            inner_filter_effect_attenuation = 1. - alpha/2. + (alpha**2)/6. - (alpha**3)/24. + (alpha**4)/120.
                        if geometry == 'bottom':
                            # Include additional term for bottom detection geometry.
                            inner_filter_effect_attenuation *= np.exp(-ELC_emission)

                        # Compute incident intensity
                        intensity = (geometry == 'top') * top_illumination_intensity + (geometry == 'bottom') * bottom_illumination_intensity # select appropriate illumination intensity
                        intensity *= inner_filter_effect_attenuation # include inner filter effects

                        # Compute fluorescence intensity
                        fluorescence = intensity * plate_fluorescence # start by including plate background fluorescence
                        for (quantum_yield, log_concentration) in zip(quantum_yields, log_concentrations):
                            fluorescence += intensity * quantum_yield * np.exp(log_concentration)
                        return fluorescence
                    self.model[name] = fluorescence
                    self.parameter_names['fluorescence'].append(name)

                    measured_fluorescence = measurements['fluorescence'][(excitation_wavelength, emission_wavelength, geometry)]
                    name = 'measured %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, wellname(well), excitation_wavelength, emission_wavelength)
                    self.model[name] = pymc.Normal(name, mu=fluorescence_model, tau=self.model['%s fluorescence precision' % geometry], observed=True, value=measured_fluorescence)
                    self.parameter_names['fluorescence'].append(name)


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
