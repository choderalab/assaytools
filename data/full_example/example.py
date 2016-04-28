"""
Analyze Abl:Gefinitinib singlet assay.

"""

from autoprotocol.unit import Unit

#
# This information is different for each experiment
#

# TODO: Refactor into dict()

d300_xml_filename = 'Src_Bos_Ima_96well_Mar2015 2015-03-07 1736.DATA.xml' # HP D300 dispense simulated DATA file
infinite_xml_filename = 'Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml' # Tecan Infinite plate reader output data
dmso_stocks_csv_filename = 'DMSOstocks-Sheet1.csv' # CSV file of DMSO stock inventory
hpd300_fluids = ['GEF001', 'IMA001', 'DMSO'] # uuid of DMSO stocks from dmso_stocks_csv_filename (or 'DMSO' for pure DMSO) used to define HP D300 XML <Fluids> block
receptor_species = 'Abl(D382N)' # receptor name (just used for convenience)
protein_absorbance = 4.24 # absorbance reading of concentrated protein stock before dilution
protein_extinction_coefficient = 49850 # 1/M/cm extinction coefficient for protein
protein_molecular_weight = 41293.2 # g/mol protein molecular weight
ul_protein_stock = 165.8 # uL protein stock solution used to make 1 uM protein stock
ml_buffer = 14.0 # mL buffer used to make 1 uM protein stock
rows_to_analyze = ['A', 'B'] # rows to analyze
assay_volume = Unit(100.0, 'microliters') # quantity of protein or buffer dispensed into plate

#
# This is all boilerplate analysis for singlet assays, and can be wrapped into a helper function
#

# TODO: Refactor into helper class or function

# Read DMSO stock solutions from inventory CSV file
from assaytools.experiments import DMSOStockSolutions, DMSO, Buffer, ProteinSolution
solutions = DMSOStockSolutions(dmso_stocks_csv_filename) # all solutions from DMSO stocks inventory

# Enumerate all ligand species from DMSO stocks.
ligand_species = set( [ solution.species for solution in solutions.values() if (solution.species != None)] )

# Add buffer and protein stock solutions
solutions['buffer'] = Buffer(name='20 mM Tris buffer')
solutions[receptor_species] = ProteinSolution(name='1 uM %s' % receptor_species, species=receptor_species, buffer=solutions['buffer'],
absorbance=protein_absorbance, extinction_coefficient=protein_extinction_coefficient, molecular_weight=protein_molecular_weight, ul_protein_stock=ul_protein_stock, ml_buffer=ml_buffer)

# Populate the Container data structure with well contents and measurements
from assaytools.experiments import provision_assay_plate, dispense_evo, dispense_hpd300, read_infinite
plate = provision_assay_plate(name='assay-plate', plate_type='4titude 4ti-0223')
dispense_evo(plate, solution=solutions[receptor_species], volume=assay_volume, rows=['A', 'C', 'E', 'G'])
dispense_evo(plate, solution=solutions['buffer'], volume=assay_volume, rows=['B', 'D', 'F', 'H'])
dispense_hpd300(plate, solutions=[solutions[id] for id in hpd300_fluids], xml_filename=d300_xml_filename)
read_infinite(plate, xml_filename=infinite_xml_filename)

# Select specified rows for analysis.
from autoprotocol import WellGroup
well_group = WellGroup([well for well in plate.all_wells() if (well.humanize()[0] in rows_to_analyze)])

# Create a model
from assaytools.analysis import CompetitiveBindingAnalysis
experiment = CompetitiveBindingAnalysis(solutions=solutions, wells=well_group, receptor_name=receptor_species, ligand_names=ligand_species)

# Fit the maximum a posteriori (MAP) estimate
map_fit = experiment.map_fit()

# Run some MCMC sampling and return the MCMC object
mcmc = experiment.run_mcmc()

# Show summary
experiment.show_summary(mcmc, map_fit)

# Generate plots
plots_filename = 'plots.pdf'
experiment.generate_plots(mcmc, pdf_filename=plots_filename)
