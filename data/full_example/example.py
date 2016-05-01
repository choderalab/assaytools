"""
Analyze Abl:Gefinitinib singlet assay.

"""

from autoprotocol.unit import Unit
from assaytools.experiments import SingletAssay

#
# This information is different for each experiment.
# We use a 'dict' so that we can later store this information in a JSON database or something.
#

params = {
    'd300_xml_filename' : 'Src_Bos_Ima_96well_Mar2015 2015-03-07 1736.DATA.xml' # HP D300 dispense simulated DATA file
    'infinite_xml_filename' : 'Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml' # Tecan Infinite plate reader output data
    'dmso_stocks_csv_filename' : 'DMSOstocks-Sheet1.csv' # CSV file of DMSO stock inventory
    'hpd300_fluids' : ['GEF001', 'IMA001', 'DMSO'] # uuid of DMSO stocks from dmso_stocks_csv_filename (or 'DMSO' for pure DMSO) used to define HP D300 XML <Fluids> block
    'receptor_species' = 'Abl(D382N)' # receptor name (just used for convenience)
    'protein_absorbance' = 4.24 # absorbance reading of concentrated protein stock before dilution
    'protein_extinction_coefficient' = Unit(49850, '1/molar/centimeter') # 1/M/cm extinction coefficient for protein
    'protein_molecular_weight' = Unit(41293.2, 'daltons') # g/mol protein molecular weight
    'protein_stock_volume' = Unit(165.8, 'microliters') # uL protein stock solution used to make 1 uM protein stock
    'buffer_volume' = Unit(14.0, 'milliliters') # mL buffer used to make 1 uM protein stock
    'rows_to_analyze' = ['A', 'B'] # rows to analyze
    'assay_volume' = Unit(100.0, 'microliters') # quantity of protein or buffer dispensed into plate
}

# Create a single-point (singlet) assay.
assay = SingletAssay(**params)

# Fit the maximum a posteriori (MAP) estimate
map_fit = assay.experiment.map_fit()

# Run some MCMC sampling and return the MCMC object
mcmc = assay.experiment.run_mcmc()

# Show summary
assay.experiment.show_summary(mcmc, map_fit)

# Generate plots
plots_filename = 'plots.pdf'
assay.experiment.generate_plots(mcmc, pdf_filename=plots_filename)

#
# This is all boilerplate analysis for singlet assays, and can be wrapped into a helper function
#
