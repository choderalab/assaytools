"""
Analyze Abl:Gefinitinib singlet (single fluorescent inhibitor) assay.

"""

from autoprotocol.unit import Unit
from assaytools.experiments import SingleWavelengthAssay

#
# Single-wavelength competition assay, but only the data for a single probe will be used
#

params = {
    'd300_xml_filename' : 'data/Src_Bos_Ima_96well_Mar2015 2015-03-07 1736.DATA.xml', # HP D300 dispense simulated DATA file
    'infinite_xml_filename' : 'data/Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml', # Tecan Infinite plate reader output data
    'dmso_stocks_csv_filename' : 'data/DMSOstocks-Sheet1.csv', # CSV file of DMSO stock inventory
    'hpd300_fluids' : ['GEF001', 'IMA001', 'DMSO'], # uuid of DMSO stocks from dmso_stocks_csv_filename (or 'DMSO' for pure DMSO) used to define HP D300 XML <Fluids> block
    'hpd300_plate_index' : 1, # plate index for HP D300 dispense script
    'receptor_species' : 'Abl(D382N)', # receptor name (just used for convenience)
    'protein_absorbance' : 4.24, # absorbance reading of concentrated protein stock before dilution
    'protein_extinction_coefficient' : Unit(49850, '1/(moles/liter)/centimeter'), # 1/M/cm extinction coefficient for protein
    'protein_molecular_weight' : Unit(41293.2, 'daltons'), # g/mol protein molecular weight
    'protein_stock_volume' : Unit(165.8, 'microliters'), # uL protein stock solution used to make 1 uM protein stock
    'buffer_volume' : Unit(14.0, 'milliliters'), # mL buffer used to make 1 uM protein stock
    'rows_to_analyze' : ['A', 'B'], # rows to analyze
    'assay_volume' : Unit(100.0, 'microliters'), # quantity of protein or buffer dispensed into plate
    'measurements_to_analyze' : ['fluorescence top'], # which measurements to analyze (if specified -- this is optional)
    'wavelengths_to_analyze' : ['280:nanometers', '480:nanometers'], # which wavelengths to analyze (if specified -- this is optional)
}

# Create a single-wavelength assay.
assay = SingleWavelengthAssay(**params)

# Run some MCMC sampling and return the MCMC object
mcmc = assay.experiment.run_mcmc()

# Show summary
assay.experiment.show_summary(mcmc)

# Generate plots
plots_filename = 'plots.pdf'
assay.experiment.generate_plots(mcmc, pdf_filename=plots_filename)
