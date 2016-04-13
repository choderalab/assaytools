from autoprotocol.container import Well, WellGroup, Container
from autoprotocol.container_type import ContainerType
from autoprotocol.unit import Unit
import numpy as np

#
# Load DMSO stocks
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
# Define solutions
#

# TODO: Refine this manner of defining solutions.

class Solution(object):
    """A solution with a defined concentration and uncertainty.
    """
    def __init__(self, name, concentration, uncertainty):
        self.name = name
        self.concentration = concentration
        self.uncertainty = uncertainty

class ProteinSolution(Solution):
    """A protein solution prepared spectrophotometrically.
    """
    def __init__(self, name, absorbance, extinction_coefficient, molecular_weight, ul_protein_stock, ml_buffer):
        """
        Parameters
        ----------
        protein_name : str
            Name of the protein
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
        self.name = name
        concentration = (absorbance / extinction_coefficient) * (ul_protein_stock/1000.0) / ml_buffer # M
        spectrophotometer_CV = 0.10 # TODO: Specify in assumptions YAML file
        self.concentration = Unit(concentration, 'moles/liter')
        self.uncertainty = spectrophotometer_CV * self.concentration # TODO: Compute more precisely

class DMSOStockSolution(Solution):
    """A DMSO stock solution.
    """
    def __init__(self, dmso_stock):
        """
        Parameters
        ----------
        dmso_stock : dict
            The dictionary containing 'id', 'compound_name', 'compound mass (mg)', 'molecular weight', 'purity', 'solvent_mass'
        """
        self.name = '10 mM ' + dmso_stock['compound name'] + ' DMSO stock'
        dmso_density = 1.1004 # g/cm3
        mass_uncertainty = 0.01 # TODO: Calculate from balance precision
        concentration = dmso_stock['compound mass (mg)'] / 1000 * dmso_stock['molecular weight'] * dmso_stock['purity'] / (dmso_stock['solvent mass (g)'] * dmso_density / 1000) # mol/liter
        self.concentration = Unit(concentration, 'moles/liter')
        self.uncertainty = mass_uncertainty * concentration

solutions = dict()
solutions['Abl'] = ProteinSolution(name='1 uM Abl', absorbance=4.24, extinction_coefficient=49850, molecular_weight=41293.2, ul_protein_stock=165.8, ml_buffer=14.0)
solutions['BOS'] = DMSOStockSolution(dmso_stocks['BOS001'])
solutions['BSI'] = DMSOStockSolution(dmso_stocks['BOI001'])
solutions['GEF'] = DMSOStockSolution(dmso_stocks['GEF001'])
solutions['ERL'] = DMSOStockSolution(dmso_stocks['ERL001'])

#
# Read ligand concentrations
#



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

# Generate a unique container ID
import uuid
id = str(uuid.uuid4())

# Define the container
container = Container(name="assay-plate", id=id, container_type=container_type)

#
# Assay Parameters
#

Lstated = [ Unit(x, "moles/liter") for x in [20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0] ] # gefitinib concentration
Pstated = Unit(0.5e-6, "moles/liter")
well_volume = Unit(100, "milliliter")
receptor_name = 'Src'
ligand_name = 'gefitinib'

# Load data into well format.
filename = "Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml"
from assaytools import platereader
data = platereader.read_icontrol_xml(filename)

# Define wells for fluorescence assay (the verbose way; we'd provide convenience methods to format the plate)
ncolumns = 12
well_group = WellGroup([])
for column in range(ncolumns):
    for row in ['A', 'B']:
        well_name = row + str(column+1)
        well = container.well(well_name)
        well.set_volume(well_volume)
        well.set_properties({'area' : well_area})

        # Set concentrations of well components
        if row == 'A':
            well.set_properties({'concentrations' : {receptor_name : Pstated, ligand_name : Lstated[column]} })
            well.set_properties({'concentration_uncertainties' : {receptor_name  : 0.10 * Pstated, ligand_name : 0.08 * Lstated[column]} })
        elif row == 'B':
            well.set_properties({'concentrations' : {ligand_name : Lstated[column]} })
            well.set_properties({'concentration_uncertainties' : {ligand_name : 0.08 * Lstated[column]} })

        # Attach plate reader data
        for wavelength in ['280', '350', '480']:
            well.set_properties({' absorbance_' + wavelength + 'nm' : data['Abs_' + wavelength]['well_data'][well_name] })
        emission_wavelength = '480'
        for excitation_wavelength in ['280', '350']:
            well.set_properties({'fluorescence_top_ex' + excitation_wavelength + 'nm_em' + emission_wavelength + 'nm' : data[excitation_wavelength + '_TopRead']['well_data'][well_name]})
            well.set_properties({'fluorescence_bottom_ex' + excitation_wavelength + 'nm_em' + emission_wavelength + 'nm' : data[excitation_wavelength + '_BottomRead']['well_data'][well_name]})

        # Add to well group
        well_group.append(well)

# TODO: Analyze the well group.
# This is just a non-working stub for now.
# Create a model
from assaytools.analysis import CompetitiveBindingAnalysis
model = CompetitiveBindingAnalysis(wells=well_group, receptor=receptor_name, ligands=[ligand_name])
# fit the maximum a posteriori (MAP) estimate
map = model.map_fit()
# run some MCMC sampling and return the MCMC object
mcmc = model.run_mcmc()
