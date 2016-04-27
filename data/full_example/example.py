"""
Initial crude example of new general analysis API that facilitates analysis of multiple experiments, competition assays, spectral assays, etc.

NOTE: This is just an initial crude long-form draft. All of the reusable components will be organized into tools to make general usage convenient.
"""

from autoprotocol.container import Well, WellGroup, Container
from autoprotocol.container_type import ContainerType
from autoprotocol.unit import Unit
import numpy as np

#
# Define assay plate container
#

# Define the container type for 4titude 4ti-0223.
# info: http://4ti.co.uk/microplates/black-clear-bottom/96-well/
# drawing: http://4ti.co.uk/files/1614/0542/7662/4ti-0223_Marketing_Drawing.pdf
# All arguments to ContainerType are required!
capabilities = ['pipette', 'spin', 'absorbance', 'fluorescence', 'luminescence', 'incubate', 'gel_separate', 'cover', 'seal', 'stamp', 'dispense']
well_diameter = Unit(6.30, "millimeters")
well_area = np.pi * (well_diameter/2)**2
container_type = ContainerType(name='4titude 4ti-0223', is_tube=False, well_count=96, well_depth_mm=Unit(11.15, 'millimeter'), well_volume_ul=Unit(300, 'milliliter'), well_coating='polystyrene', sterile=False, capabilities=capabilities, shortname='4ti-0223', col_count=12, dead_volume_ul=Unit(20,'milliliter'), safe_min_volume_ul=Unit(50, 'milliliter'))

# Generate the container a unique container ID
# This might come from a barcode in the future.
import uuid
id = str(uuid.uuid4())

# Define the container
container = Container(name="assay-plate", id=id, container_type=container_type)

# Initialize well properties for this container
for well in container.all_wells():
    well.set_volume(Unit(0.0, 'microliters')) # well starts empty
    well.set_properties({'contents' : dict()}) # contains the contents
    well.set_properties({'area' : well_area})

#
# Load information about DMSO stocks
# In future, this might come from a database query that returns JSON.
#

import csv
dmso_stocks_csv_filename = 'DMSOstocks-Sheet1.csv'
dmso_stocks = dict()
with open(dmso_stocks_csv_filename, 'rb') as csvfile:
     csvreader = csv.DictReader(csvfile, delimiter=',', quotechar='|')
     for row in csvreader:
         if row['id'] != '':
            for key in ['compound mass (mg)', 'purity', 'solvent mass (g)', 'molecular weight']:
                row[key] = float(row[key])
            dmso_stocks[row['id']] = row

#
# Define simple solutions
# Solutions are the fundamental liquids that have concentrations (and concentration uncertainties) for a single species.
# In future, we will associate densities, allow mixing of solutions, propagate uncertainty.
# Solutions have unknown *true* concentrations, except for buffer solutions which have exactly zero concentration.
#

class Buffer(object):
    def __init__(self, name):
        """
        """
        self.name = name
        self.component = None
        self.concentration = Unit(0.0, 'moles/liter')
        self.uncertainty = Unit(0.0, 'moles/liter')

DMSO = Buffer('DMSO')

class Solution(object):
    """A solution of a single component with a defined concentration and uncertainty.
    """
    def __init__(self, name, species, buffer, concentration, uncertainty):
        """
        Parameters
        ----------
        name : str
            A descrptive name for the solution.
        species : str
            The name of the species dissolved in buffer.
        buffer : Buffer
            The buffer the solution is prepared in.
        concentration : Unit
            The concentration of the species in the solution.
        uncertainty : Unit
            The concentration uncertainty of the species in the solution.
        """
        self.name = name
        self.component = component
        self.buffer = buffer
        self.concentration = concentration
        self.uncertainty = uncertainty

class ProteinSolution(Solution):
    """A protein solution in buffer prepared spectrophotometrically.
    """
    def __init__(self, name, species, buffer, absorbance, extinction_coefficient, molecular_weight, ul_protein_stock, ml_buffer):
        """
        Parameters
        ----------
        protein_name : str
            A descreptive name of the solution.
        species : str
            The protein name
        buffer : BufferSolution
            The corresponding buffer solution.
        absorbance : float
            Absorbance for 1 cm equivalent path length
        extinction_coefficient : float
            Extinction coefficient (1/M/cm)
        molecular_weight : float
            Molecular weight (g/mol or Daltons)
        ul_protein_stock : float
            Microliters of protein stock added to make solution
        ml_buffer : float
            Milliliters of buffer added to make solution

        """
        self.name = name + ' in ' + buffer.name
        self.species = species
        self.buffer = buffer
        concentration = (absorbance / extinction_coefficient) * (ul_protein_stock/1000.0) / ml_buffer # M
        spectrophotometer_CV = 0.10 # TODO: Specify in assumptions YAML file
        self.concentration = Unit(concentration, 'moles/liter')
        self.uncertainty = spectrophotometer_CV * self.concentration # TODO: Compute more precisely

class DMSOStockSolution(Solution):
    """A stock solution representing a compound dissolved in DMSO.
    """
    def __init__(self, dmso_stock):
        """
        Parameters
        ----------
        dmso_stock : dict
            The dictionary containing 'id', 'compound_name', 'compound mass (mg)', 'molecular weight', 'purity', 'solvent_mass'
        """
        self.name = '10 mM ' + dmso_stock['compound name'] + ' DMSO stock'
        self.species = dmso_stock['compound name']
        dmso_density = 1.1004 # g/cm3
        mass_uncertainty = 0.01 # TODO: Calculate from balance precision
        concentration = dmso_stock['compound mass (mg)'] / 1000 * dmso_stock['molecular weight'] * dmso_stock['purity'] / (dmso_stock['solvent mass (g)'] * dmso_density / 1000) # mol/liter
        self.concentration = Unit(concentration, 'moles/liter')
        self.uncertainty = mass_uncertainty * self.concentration
        self.buffer = DMSO

solutions = dict()
solutions['buffer'] = Buffer(name='20 mM Tris buffer')
solutions['Abl'] = ProteinSolution(name='1 uM Abl', species='Abl', buffer=solutions['buffer'], absorbance=4.24, extinction_coefficient=49850, molecular_weight=41293.2, ul_protein_stock=165.8, ml_buffer=14.0)
solutions['DMSO'] = DMSO
solutions['BOS'] = DMSOStockSolution(dmso_stocks['BOS001'])
solutions['BSI'] = DMSOStockSolution(dmso_stocks['BOI001'])
solutions['GEF'] = DMSOStockSolution(dmso_stocks['GEF001'])
solutions['ERL'] = DMSOStockSolution(dmso_stocks['ERL001'])
solutions['IMA'] = DMSOStockSolution(dmso_stocks['IMA001'])
receptor_name = 'Abl'
ligand_names = ['BOS', 'IMA']

#
# Dispense protein solution
#

Tecan_EVO_CV = 0.004 # for 100 uL volume; TODO: read these from YAML file
assay_volume = Unit(100.0, 'microliters')
for well in container.all_wells():
    well_name = well.humanize()
    contents = well.properties['contents']
    if well_name[0] in ['A', 'C', 'E', 'G']:
        contents['Abl'] = (assay_volume, Tecan_EVO_CV * assay_volume) # volume, error
    else:
        contents['buffer'] = (assay_volume, Tecan_EVO_CV * assay_volume) # volume, error
    well.set_volume(well.volume + assay_volume)

#
# Read ligand concentrations
#

d300_xml_filename = 'Src_Bos_Ima_96well_Mar2015 2015-03-07 1736.DATA.xml'

# Read HP D300 XML file
import xml.etree.ElementTree as ET
tree = ET.parse(d300_xml_filename)
root = tree.getroot()

# Read fluids
fluids = root.findall('./Fluids/Fluid')

# TODO: Rewrite fluid names to match stock names.
fluids[0].attrib['Name'] = 'GEF'
fluids[1].attrib['Name'] = 'IMA'
fluids[2].attrib['Name'] = 'DMSO'

def humanize_d300_well(row, column):
    """
    Return the humanized version of a D300 well index.
    """
    return chr(row + 97) + str(column+1)

# Read dispensed volumes
# TODO: Read only one plate
dispensed = root.findall('./Dispensed')[0]
volume_unit = dispensed.attrib['VolumeUnit'] # dispensed volume unit
wells = dispensed.findall('Plate/Well')
HP_D300_CV = 0.08 # CV for HP D300; TODO: Read this from a YAML file
for well in wells:
    row = int(well.attrib['R'])
    column = int(well.attrib['C'])
    dispensed_fluids = well.findall('Fluid')
    for dispensed_fluid in dispensed_fluids:
        dispensed_fluid_index = int(dispensed_fluid.attrib['Index'])
        fluid_name = fluids[dispensed_fluid_index].attrib['Name']
        total_volume = Unit(float(dispensed_fluid.attrib['TotalVolume']), volume_unit)
        detail = dispensed_fluid.findall('Detail')
        # Gather details
        time = detail[0].attrib['Time']
        cassette = int(detail[0].attrib['Cassette'])
        head = int(detail[0].attrib['Head'])
        dispensed_volume = Unit(float(detail[0].attrib['Volume']), volume_unit)

        # Retrieve well from container.
        well_name = humanize_d300_well(row, column)
        well = container.well(well_name)
        contents = well.properties['contents']
        contents[fluid_name] = (dispensed_volume, HP_D300_CV * dispensed_volume) # (volume, error)
        well.set_volume(well.volume + dispensed_volume)

#
# Load plate reader data.
#

filename = "Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml"
from assaytools import platereader
data = platereader.read_icontrol_xml(filename)

# Define wells for fluorescence assay (the verbose way; we'd provide convenience methods to format the plate)
for well in container.all_wells():
    well_name = well.humanize()

    # Attach plate reader data
    measurements = dict()
    for key in data:
        if key.startswith('Abs_'):
            # absorbance
            [prefix, wavelength] = key.split('_')
            wavelength = wavelength + ':nanometers'
            if 'absorbance' not in measurements:
                measurements['absorbance'] = dict()
            measurements['absorbance'][wavelength] = float(data[key][well_name])
        elif (key.endswith('_TopRead') or key.endswith('_BottomRead')):
            # top fluorescence read
            [wavelength, suffix] = key.split('_')
            excitation_wavelength = wavelength + ':nanometers'
            emission_wavelength = '450:nanometers'
            if key.endswith('_TopRead'):
                geometry = 'top'
            else:
                geometry = 'bottom'
            if 'fluorescence' not in measurements:
                measurements['fluorescence'] = dict()
            measurements['fluorescence'][(excitation_wavelength, emission_wavelength, geometry)] = float(data[key] [well_name])

    well.properties['measurements'] = measurements

#
# Actual API example
#

from assaytools.experiments import DMSOStockSolutions, BufferSolution, ProteinSolution, AssayPlate
from assaytools.analysis import CompetitiveBindingAnalysis

# Read DMSO stock solutions from inventory CSV file
dmso_stocks_csv_filename = 'DMSOstocks-Sheet1.csv'
solutions = DMSOStockSolutions(dmso_stocks_csv_filename)
ligand_solutions = (solutions['BOS001'], solutions['BSI001'], solutions['GEF001'], solutions['GEF001'])

# Define receptor and ligand species names.
receptor_species = 'Abl(D382N)'
ligand_species = [ ligand_solutions.species for solution in ligand_solutions ]

# Add buffer and protein stock solutions
solutions['buffer'] = Buffer(name='20 mM Tris buffer')
solutions[receptor_species] = ProteinSolution(name='1 uM Abl D382N', species=receptor_species, buffer=solutions['buffer'], absorbance=4.24, extinction_coefficient=49850, molecular_weight=41293.2, ul_protein_stock=165.8, ml_buffer=14.0)

# Populate the Container data structure with well contents and measurements
d300_xml_filename = 'LRL_Src_Bos_2rows_1_2 2015-09-11 1048.DATA.xml'
infinite_xml_filename = 'Abl_D382N_Bos_20160311_132205.xml'
plate = AssayPlate(protein_solution=solutions['Abl'], buffer_solution=solutions['buffer'], ligand_solutions=ligand_solutions,
    d300_xml_filename=d300_xml_filename, infinite_xml_filename=infinite_xml_filename)

# Create a model
experiment = CompetitiveBindingAnalysis(solutions=solutions, wells=plate.all_wells(), receptor_species=receptor_species, ligand_species=ligand_species)

# Fit the maximum a posteriori (MAP) estimate
map_fit = experiment.map_fit()

# Run some MCMC sampling and return the MCMC object
mcmc = experiment.run_mcmc()

# Show summary
experiment.show_summary(mcmc, map_fit)

# Generate plots
plots_filename = 'plots.pdf'
experiment.generate_plots(mcmc, pdf_filename=plots_filename)
