"""
Classes for the analysis of fluorescence assay data using pymc3

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy as np
import pymc3
from pymc3 import exp, log, cos, sin, tan, cosh, sinh, tanh, sqr, sqrt, erf, erfinv, dot
from pymc3 import maximum, minimum, sgn, ceil, floor
import theano
import scipy
import scipy.optimize
from theano.compile.ops import as_op
import theano.tensor
import sys

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

volume_unit = Unit(1.0, 'liter')
concentration_unit = Unit(1.0, 'moles/liter')

#=============================================================================================
# SUBROUTINES
#=============================================================================================

# NOTE: This appears in recent pymc3 dev versions, but was added after 3.0beta release
def logsumexp(x, axis=None):
     # Adapted from https://github.com/Theano/Theano/issues/1563
     x_max = theano.tensor.max(x, axis=axis, keepdims=True)
     return theano.tensor.log(theano.tensor.sum(theano.tensor.exp(x - x_max), axis=axis, keepdims=True)) + x_max

class LogNormalWrapper(pymc3.Lognormal):
    def __init__(self, mean=None, stddev=None, *args, **kwargs):
        """Specify Normal parameters, but create a lognormal distribution.

        Parameters
        ----------
        mean : float
            Mean of normal distribution.
        stddev : float
            Standard deviation of normal distribution

        """
        mu = log(mean**2 / np.sqrt(stddev**2 + mean**2))
        tau = sqrt(log(1.0 + (stddev/mean)**2))**(-2)
        super(LogNormalWrapper,self).__init__(mu=mu, tau=tau, *args, **kwargs)

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
# PyMC2 models
#=============================================================================================

class CompetitiveBindingAnalysis(object):
    def __init__(self, solutions, wells, receptor_name, DeltaG_prior='uniform'):
        """
        Parameters
        ----------
        solutions : dict of str : Solution
            `solutions[name]` is the Solution object corresponding to component `name`
        wells : autoprotocol.container.WellGroup
            Group of wells with `contents` and `measurements` properties defined
        receptor_name : str
            Name of receptor
        DeltaG_prior : str, optional, default='uniform'
            Prior to use on DeltaG values for reaction ['uniform', 'chembl']

        """
        # Store data
        self.solutions = solutions
        self.wells = wells
        self.receptor_name = receptor_name

        # Set up internal data structures.
        self.parameters = dict() # the PyMC model; self.parameters[paramname] is the PyMC variable correspoding to 'paramname'
        self.parameter_names = dict() # dict to keep track of groups of related parameter names; self.parameter_names[groupname] is the list of PyMC variable names under 'groupname'

        # DEBUG
        print('There are %d wells to analyze in the provided WellGroup' % len(wells))

        # Create the pymc3 model context
        with pymc3.Model() as model:
            self.model = model

            # Construct components of the pymc model.
            self._create_solutions_model()
            self._identify_ligand_names()
            self._create_dispensing_model()
            self._create_competitive_binding_model(receptor_name, DeltaG_prior)
            #self._create_equilibrium_concentrations_model()
            #self._create_extinction_coefficients_model()
            #self._create_absorbance_model()
            #self._create_fluorescence_model()

        # TODO: Print model parameter names
        print('Model variables:')
        for (index, var) in enumerate(self.model.vars):
            print('%8d %s' % (index, str(var)))

    def _create_solutions_model(self):
        """
        Create pymc model components for true concentrations of source receptor and ligand solutions.

        Populates the following fields:
        * parameter_names['concentrations'] : parameters associated with true concentrations of receptor and ligand solutions
        """

        # Determine solutions in use in plate
        solutions_in_use = set()
        for well in self.wells:
            for shortname in well.properties['contents']:
                solutions_in_use.add(shortname)
        print('Solutions in use: %s' % str(solutions_in_use))

        # Retain only solutions that appear in the plate
        self.solutions = { shortname : self.solutions[shortname] for shortname in solutions_in_use }

        self.parameter_names['solution concentrations'] = list()
        for solution in self.solutions.values():
            if solution.species is None:
                continue # skip buffers or pure solvents
            name = 'log concentration of %s' % solution.shortname
            self.parameters[name] = LogNormalWrapper(name, mean=solution.concentration.to_base_units().m, stddev=solution.uncertainty.to_base_units().m)
            self.parameter_names['solution concentrations'].append(name)

    def _identify_ligand_names(self):
        """
        Ligands are species that appear in wells that are not the receptor.

        """
        self.ligand_names = set()
        for solution in self.solutions.values():
            if solution.species is None:
                continue # skip buffers and pure solvents
            if solution.species == self.receptor_name:
                continue # sip receptor
            self.ligand_names.add(solution.species)
        self.ligand_names = list(self.ligand_names)

    def _create_dispensing_model(self):
        """
        Create nuisance parameters for dispensed volumes and actual concentrations of all species in each well.

        Populates the following fields:
        * parameter_names['dispensed_volumes'] : actual volumes dispensed into each well
        * parameter_names['well_volumes'] : actual well total volumes
        *
        """
        self.parameter_names['dispensed volumes'] = list() # pymc variable names associated with dispensed volumes
        self.parameter_names['well volumes'] = list() # pymc variable names associated with well volumes
        self.parameter_names['well concentrations'] = list()

        for well in self.wells:
            # Volumes dispensed into each well
            log_volumes = list() # log volumes are in Liters
            for component in well.properties['contents']:
                name = 'volume of %s dispensed into well %s' % (component, wellname(well))
                logname = 'log %s' % name
                (volume, error) = well.properties['contents'][component]
                log_volume_dispensed = LogNormalWrapper(logname, mean=volume.to_base_units().m, stddev=error.to_base_units().m)
                self.parameters[logname] = log_volume_dispensed
                #self.parameter_names['dispensed volumes'].append(logname)
                log_volumes.append(log_volume_dispensed)
                # Store real (non-log) value
                self.parameters[name] = pymc3.Deterministic(name, exp(self.parameters[logname]) / volume_unit.to_base_units().m)
                self.parameter_names['dispensed volumes'].append(name)

            # Total volume in well
            name = 'volume of well %s' % wellname(well)
            logname = 'log %s' % name
            self.parameters[logname] = pymc3.Deterministic(logname, logsumexp(log_volumes))
            #self.parameter_names['well volumes'].append(logname)
            # Store real (non-log) value
            self.parameters[name] = pymc3.Deterministic(name, exp(self.parameters[logname]) / volume_unit.to_base_units().m)
            self.parameter_names['well volumes'].append(name)

            # Total concentration of all species involving each component in well
            # TODO: In future, can we simply calculate initial concentrations of each species and use that in GeneralBindingModel instead?
            for component in well.properties['contents']:
                solution = self.solutions[component]
                if solution.species is None:
                    continue # skip buffers or pure solvents
                species = solution.species
                log_solution_concentration = self.parameters['log concentration of %s' % solution.shortname] # log concentrations are in molar
                log_solution_volume = self.parameters['log volume of %s dispensed into well %s' % (component, wellname(well))]
                log_total_volume = self.parameters['log volume of well %s' % wellname(well)]

                name = 'total concentration of %s in well %s' % (species, wellname(well))
                logname = 'log %s' % name
                self.parameters[logname] = pymc3.Deterministic(logname, log_solution_concentration + log_solution_volume - log_total_volume)
                # Store real (non-log) value
                self.parameters[name] = pymc3.Deterministic(name, exp(self.parameters[logname]) / concentration_unit.to_base_units().m)
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
            log_total_concentration = self.parameters[name]
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
        self.parameter_names['binding affinities'] = list()
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
                DeltaG = pymc3.Uniform(name, lower=DG_min, upper=DG_max) # binding free energy (kT), uniform over huge range
            elif DeltaG_prior == 'chembl':
                DeltaG = pymc3.Normal(name, mu=0, sd=12.5) # binding free energy (kT), using a Gaussian prior inspured by ChEMBL
            else:
                raise Exception("DeltaG_prior = '%s' unknown. Must be one of 'uniform' or 'chembl'." % DeltaG_prior)
            self.parameters[name] = DeltaG
            self.parameter_names['binding affinities'].append(name)
            # Create the reaction for GeneralBindingModel
            self.reactions.append( (DeltaG, {complex_name : -1, receptor_name : +1, ligand_name : +1}) )
            # Create the conservation equation for GeneralBindingModel
            self.conservation_equations.append( (ligand_name, {ligand_name : +1, complex_name : +1}) )
            # Keep track of receptor conservation.
            receptor_conservation_equation[complex_name] = +1
        # Create the receptor conservation equation for GeneralBindingModel
        self.conservation_equations.append( [receptor_name, receptor_conservation_equation] )

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
                log_total_concentration = self.parameters['log total concentration of %s in well %s' % (species, wellname(well))]
                conservation_equations.append( [log_total_concentration, stoichiometry] )

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
                    name = 'concentration of %s in well %s' % (species, wellname(well))
                    logname = 'log %s' % name
                    log_total_concentration = self.model['log total concentration of %s in well %s' % (species, wellname(well))]
                    print('%s slaved to parent' % logname)
                    self.model[logname] = pymc.Deterministic(logname, log_total_concentration)
                    #self.parameter_names['well concentrations'].append(logname)
                    # Store real (non-log) value
                    self.model[name] = pymc.Deterministic(name, exp(self.model[logname]))
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
            @as_op(itypes=[theano.tensor.lscalar, theano.tensor.dscalar, theano.tensor.dscalar], otypes=[theano.tensor.dvector])
            def log_equilibrium_concentrations(reactions=reactions, conservation_equations=conservation_equations):
                if len(reactions)==0 or len(conservation_equations)==0:
                    raise Exception(reactions, conservation_equations)
                solution = GeneralBindingModel.equilibrium_concentrations(reactions, conservation_equations)
                return solution

            self.parameters[name] = log_equilibrium_concentrations

            # Separate out individual concentration components
            for species in all_equilibrium_species:
                name = 'concentration of %s in well %s' % (species, wellname(well))
                logname = 'log %s' % name
                self.parameters[logname] = pymc3.Deterministic(logname, log_equilibrium_concentration[species])
                #self.parameter_names['well concentrations'].append(logname)
                # Store real (non-log) value
                self.parameters[name] = pymc3.Deterministic(name, exp(self.parameters[logname]))
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
        MIN_EXTINCTION_COEFFICIENT = Unit(0.1, '1/(moles/liter)/centimeter') # maximum physically reasonable extinction coefficient
        MAX_EXTINCTION_COEFFICIENT = Unit(100e3, '1/(moles/liter)/centimeter') # maximum physically reasonable extinction coefficient
        EXTINCTION_COEFFICIENT_GUESS = Unit(1.0, '1/(moles/liter)/centimeter') # maximum physically reasonable extinction coefficient

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
                    name = 'log extinction coefficient of %s at wavelength %s' % (species, wavelength)
                    log_extinction_coefficient = pymc3.Uniform(name, lower=np.log(MIN_EXTINCTION_COEFFICIENT.to_base_units().m), upper=np.log(MAX_EXTINCTION_COEFFICIENT.to_base_units().m), value=np.log(EXTINCTION_COEFFICIENT_GUESS.to_base_units().m)) # extinction coefficient or molar absorptivity for ligand, units of 1/M/cm
                    self.parameters[name] = log_extinction_coefficient
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
            self.parameters['log absorbance error'] = pymc3.Uniform('log absorbance error', lower=-10, upper=0, value=np.log(0.01))
            self.parameters['absorbance error'] = pymc3.Lambda('absorbance error', lambda log_sigma=self.parameters['log absorbance error'] : np.exp(log_sigma), trace=True)
            self.parameter_names['absorbance'] = ['log absorbance error', 'absorbance error']
            # Prior on plate absorbance at each wavelength
            for wavelength in self.all_wavelengths:
                name = 'plate absorbance at wavelength %s' % wavelength
                self.parameters[name] = pymc3.Uniform(name, lower=0.0, upper=Amax, value=0.0)
                self.parameter_names['absorbance'].append(name)

        # Absorbance measurements (for wells that have them)
        for well in self.wells:
            # Compute path length
            log_well_volume = self.parameters['log volume of well %s' % wellname(well)]
            all_equilibrium_species_in_well = self._all_equilibrium_species_in_well(well)
            measurements = well.properties['measurements']
            if 'absorbance' in measurements:
                for wavelength in measurements['absorbance'].keys():
                    # Compile all parameters needed for absorbance model
                    extinction_coefficients = [ self.parameters['log extinction coefficient of %s at wavelength %s' % (species, wavelength)] for species in all_equilibrium_species_in_well ]
                    log_concentrations = [ self.parameters['log extinction coefficient of %s at wavelength %s' % (species, wavelength)] for species in all_equilibrium_species_in_well ]
                    plate_absorbance = self.parameters['plate absorbance at wavelength %s' % wavelength]
                    log_well_area = np.log(well.container.container_type.well_area.to_base_units().m)

                    # Add computed absorbance model
                    name = 'computed absorbance of well %s at wavelength %s' % (wellname(well), wavelength)
                    @pymc3.deterministic(name=name)
                    def absorbance_model(log_concentrations=log_concentrations, log_extinction_coefficients=log_extinction_coefficients, plate_absorbance=plate_absorbance, log_well_volume=log_well_volume):
                        log_path_length = log_well_volume - log_well_area
                        absorbance = 0.0
                        for (log_extinction_coefficient, log_concentration) in zip(log_extinction_coefficients, log_concentrations):
                            absorbance += np.exp(log_extinction_coefficient + log_path_length + log_concentration)
                        absorbance += plate_absorbance
                        return absorbance
                    self.parameters[name] = absorbance_model
                    self.parameter_names['absorbance'].append(name)

                    # Add measured absorbance model
                    measured_absorbance = measurements['absorbance'][wavelength]
                    name = 'measured absorbance of well %s at wavelength %s' % (wellname(well), wavelength) # TODO: Include plate name
                    self.parameters[name] = pymc3.Normal(name, mu=absorbance_model, sd=self.parameters['absorbance error'], observed=True, value=measured_absorbance)
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
                quantum_yield = pymc3.Uniform(name, lower=0.0, upper=1.0, value=0.001)
                self.parameters[name] = quantum_yield
                self.parameter_names['fluorescence'].append(name)

                # TODO: If we had an estimate of quantum yield of some reference species, we could restraint these much better.
                for species in self.all_species:
                    name = 'quantum yield of %s for fluorescence excitation at %s and emission at %s' % (species, excitation_wavelength, emission_wavelength)
                    quantum_yield = pymc3.Uniform(name, lower=0.0, upper=1.0, value=0.1)
                    self.parameters[name] = quantum_yield
                    self.parameter_names['fluorescence'].append(name)

            MIN_LOG_FLUORESCENCE_INTENSITY = -10 # TODO: Determine minimum possible fluorescence intensity
            MAX_LOG_FLUORESCENCE_INTENSITY = +20 # TODO: Determine maximum possible fluorescence intensity
            MIN_LOG_FLUORESCENCE_UNCERTAINTY = -5
            MAX_LOG_FLUORESCENCE_UNCERTAINTY = +8 # TODO: MAKE THIS LARGER

            # Fluorescence intensities * gains
            # TODO: If multiple gains are in use, slave them together through this intensity times a fixed gain factor.
            if fluorescence_top:
                name = 'top fluorescence log illumination intensity'
                self.parameters[name] = pymc3.Uniform(name, lower=MIN_LOG_FLUORESCENCE_INTENSITY, upper=MAX_LOG_FLUORESCENCE_INTENSITY, value=0)
                self.parameter_names['fluorescence'].append(name)

                name = 'top fluorescence log uncertainty'
                self.parameters[name] = pymc3.Uniform(name, lower=MIN_LOG_FLUORESCENCE_UNCERTAINTY, upper=MAX_LOG_FLUORESCENCE_UNCERTAINTY, value=0)
                self.parameter_names['fluorescence'].append(name)

                name = 'top fluorescence error'
                self.parameters[name] = pymc3.Lambda(name, lambda log_sigma=self.parameters['top fluorescence error'] : np.exp(log_sigma))
                self.parameter_names['fluorescence'].append(name)

            if fluorescence_bottom:
                name = 'bottom fluorescence log illumination intensity'
                self.parameters[name] = pymc3.Uniform(name, lower=MIN_LOG_FLUORESCENCE_INTENSITY, upper=MAX_LOG_FLUORESCENCE_INTENSITY, value=0)
                self.parameter_names['fluorescence'].append(name)

                name = 'bottom fluorescence log uncertainty'
                self.parameters[name] = pymc3.Uniform(name, lower=MIN_LOG_FLUORESCENCE_UNCERTAINTY, upper=MAX_LOG_FLUORESCENCE_UNCERTAINTY, value=0)
                self.parameter_names['fluorescence'].append(name)

                name = 'bottom fluorescence error'
                self.parameters[name] = pymc3.Lambda(name, lambda log_sigma=self.parameters['bottom fluorescence error'] : np.exp(log_sigma))
                self.parameter_names['fluorescence'].append(name)

        # Fluorescence measurements (for wells that have them)
        for well in self.wells:
            log_well_volume = self.parameters['log volume of well %s' % wellname(well)]
            measurements = well.properties['measurements']
            if 'fluorescence' in measurements:
                for (excitation_wavelength, emission_wavelength, geometry) in measurements['fluorescence'].keys():
                    # Extract extinction coefficients and concentrations for all species in well and pack them into lists
                    all_species = self._all_equilibrium_species_in_well(well)
                    quantum_yields = [ self.parameters['quantum yield of %s for fluorescence excitation at %s and emission at %s' % (species, excitation_wavelength, emission_wavelength)] for species in all_species ]
                    log_excitation_extinction_coefficients = [ self.parameters['log extinction coefficient of %s at wavelength %s' % (species, excitation_wavelength)] for species in all_species ]
                    log_emission_extinction_coefficients = [ self.parameters['log extinction coefficient of %s at wavelength %s' % (species, emission_wavelength)] for species in all_species ]
                    log_concentrations = [ self.parameters['log concentration of %s in well %s' % (species, wellname(well))] for species in all_species ]
                    plate_fluorescence = self.parameters['plate background fluorescence for fluorescence excitation at %s and emission at %s' % (excitation_wavelength, emission_wavelength)]
                    top_log_illumination_intensity = self.parameters['top fluorescence log illumination intensity'] if fluorescence_top else 0
                    bottom_log_illumination_intensity = self.parameters['bottom fluorescence log illumination intensity'] if fluorescence_bottom else 0
                    name = 'computed %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, wellname(well), excitation_wavelength, emission_wavelength)
                    log_well_area = np.log(well.container.container_type.well_area.to_base_units().m)

                    @pymc3.deterministic(name=name, trace=True)
                    def fluorescence_model(log_well_volume=log_well_volume,
                        log_concentrations=log_concentrations, quantum_yields=quantum_yields,
                        log_excitation_extinction_coefficients=log_excitation_extinction_coefficients, log_emission_extinction_coefficients=log_emission_extinction_coefficients,
                        plate_fluorescence=plate_fluorescence,
                        top_log_illumination_intensity=top_log_illumination_intensity, bottom_log_illumination_intensity=bottom_log_illumination_intensity, geometry=geometry):

                        # TODO: Work entirely in log space?

                        # Compute path length
                        log_path_length = log_well_volume - log_well_area

                        # Calculate attenuation due to inner filter effects
                        # TODO: We may need to fix some scaling factors in the extinction coefficient to go between log10 and ln-based absorbance/transmission.
                        inner_filter_effect_scaling = 1.0 # scaling factor applied from inner filter effect
                        ELC_excitation = 0.0 # sum of (extinction_coefficient * path_length * concentration) for all species at excitation wavelength
                        ELC_emission   = 0.0 # sum of (extinction_coefficient * path_length * concentration) for all species at emission wavelength
                        for (log_concentration, log_excitation_extinction_coefficient, log_emission_extinction_coefficient) in zip(log_concentrations, log_excitation_extinction_coefficients, log_emission_extinction_coefficients):
                            # TODO: Accumulate using logsumexp?
                            ELC_excitation += np.exp(log_excitation_extinction_coefficient + log_path_length + log_concentration)
                            ELC_emission   += np.exp(log_emission_extinction_coefficient   + log_path_length + log_concentration)
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

                        # TODO: Fix this; ignoring for now
                        inner_filter_effect_attenuation = 1.0 # NO ATTENUATION

                        # Compute incident intensity
                        intensity = (geometry == 'top') * np.exp(top_log_illumination_intensity) + (geometry == 'bottom') * np.exp(bottom_log_illumination_intensity) # select appropriate illumination intensity
                        intensity *= inner_filter_effect_attenuation # include inner filter effects

                        # Compute fluorescence intensity
                        fluorescence = intensity * plate_fluorescence # start by including plate background fluorescence
                        for (quantum_yield, log_concentration) in zip(quantum_yields, log_concentrations):
                            fluorescence += intensity * quantum_yield * np.exp(log_concentration)

                        return fluorescence
                    self.parameters[name] = fluorescence_model
                    self.parameter_names['fluorescence'].append(name)

                    measured_fluorescence = measurements['fluorescence'][(excitation_wavelength, emission_wavelength, geometry)]
                    name = 'measured %s fluorescence of well %s at excitation wavelength %s and emission wavelength %s' % (geometry, wellname(well), excitation_wavelength, emission_wavelength)
                    self.parameters[name] = pymc3.Normal(name, mu=fluorescence_model, sd=self.parameters['%s fluorescence error' % geometry], observed=True, value=measured_fluorescence)
                    self.parameter_names['fluorescence'].append(name)


    def map_fit(self):
        """
        Find the maximum a posteriori (MAP) fit.

        Parameters
        ----------

        Returns
        -------
        map : pymc3.MAP
           The MAP fit.

        """
        print('Finding MAP fit...')
        print(self.model.unobserved_vars)
        #theano.compile.mode.Mode(linker='py', optimizer='fast_compile')
        with self.model:
            map_fit = pymc3.find_MAP(step=steps, fmin=scipy.optimize.fmin_l_bfgs_b, iprint=1) # Find starting value by optimization

        return map_fit

    def run_mcmc(self, map_fit=None, dbfilename='output'):
        """
        Sample the model with pymc using sensible defaults.

        Parameters
        ----------
        dbfilename : str, optional, default='output'
           Name of storage filename for database.

        Returns
        -------
        mcmc : pymc3.MCMC
           The MCMC samples.

        """

        ntune = 1000
        niter = 2000

        # Sample the model with pymc
        print('Sampling with MCMC...')

        #db = pymc3.backends.SQLite(name=dbfilename, model=self.model)
        db = None

        # Disable theano compilation
        #theano.compile.mode.Mode(linker='py', optimizer='fast_compile')

        with self.model:
            steps = list()
            for v in self.model.vars:
                #vars = [ v for v in self.model.vars ]
                #scaling = np.array([ 1.0 for v in self.model.vars ])
                #steps.append( pymc3.NUTS(vars=vars, scaling=scaling ) )
                #steps.append( pymc3.Metropolis(vars=vars, scaling=scaling ) )
                steps.append( pymc3.Metropolis(vars=[v], scaling=1.0 ) )

            trace = pymc3.sample(draws=niter, tune=ntune, step=steps, start=map_fit, trace=db, progressbar=True)

        return trace

    def show_summary(self, mcmc, map_fit=None):
        """
        Show summary statistics of MCMC and (optionally) MAP estimates.

        Parameters
        ----------
        mcmc : pymc3.MCMC
           MCMC samples
        map_fit : pymc3.MAP, optional, default=None
           The MAP fit.

        TODO
        ----
        * Automatically determine appropriate number of decimal places from statistical uncertainty.
        * Automatically adjust concentration units (e.g. pM, nM, uM) depending on estimated affinity.

        """

        # Compute summary statistics
        alpha = 0.95 # confidence interval width
        from scipy.stats import bayes_mvs
        for group in self.parameter_names:
            print(group)
            for name in self.parameter_names[group]:
                if map_fit:
                    mle = getattr(map_fit, name).value
                else:
                    mle = getattr(mcmc, name).trace().mean()
                mean_cntr, var_cntr, std_cntr = bayes_mvs(getattr(mcmc, name).trace(), alpha=alpha)
                (center, (lower, upper)) = mean_cntr
                print("%64s : %5.1f [%5.1f, %5.1f] kT" % (name, mle, lower, upper))
            print('')

    def generate_plots(self, trace, map_fit=None, pdf_filename=None):
        """
        Generate interactive or PDF plots from MCMC trace.

        Parameters
        ----------
        trace :
           MCMC sampling trace
        map_fit : pymc3.MAP, optional, default=None
           Plot the maximum a posteriori (MAP) estimate if provided.
        pdf_filename : str, optional, default=None
           If specified, generate a PDF containing plots.

        """
        pymc3.traceplot(trace)