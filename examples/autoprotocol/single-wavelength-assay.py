"""
Analyze Abl:Gefinitinib singlet (single fluorescent inhibitor) assay.

"""

from autoprotocol.unit import Unit
from assaytools.experiments import SingletAssay

#
# Competition assay
#

params = {
    'd300_xml_filename' : 'data/Src_Bos_Ima_96well_Mar2015 2015-03-07 1736.DATA.xml', # HP D300 dispense simulated DATA file
    'infinite_xml_filename' : 'Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml', # Tecan Infinite plate reader output data
    'dmso_stocks_csv_filename' : 'DMSOstocks-Sheet1.csv', # CSV file of DMSO stock inventory
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

# Create a single-point (singlet) assay.
assay = SingletAssay(**params)

# Try emcee
import numpy as np
import emcee
import pymc
# This is the likelihood function for emcee
model = assay.experiment.model
from pymc.Node import ZeroProbability
from numpy import Inf
def lnprob(vals): # vals is a vector of parameter values to try
    # Set each random variable of the pymc model to the value
    # suggested by emcee
    for val, var in zip(vals, model.stochastics):
        var.value = val

    # Calculate logp
    try:
        logp = model.logp
    #except ZeroProbability as e:
    except Exception as e:
        print(e)
        logp = -Inf
    return logp

# emcee parameters
ndim = len(model.stochastics)
nwalkers = 2*ndim
# Find MAP
print('Finding MAP estimate...')
pymc.MAP(model).fit(iterlim=5)
start = np.empty(ndim)
for i, var in enumerate(model.stochastics):
    start[i] = var.value
# sample starting points for walkers around the MAP
p0 = np.random.randn(ndim * nwalkers).reshape((nwalkers, ndim)) + start
# instantiate sampler passing in the pymc likelihood function
print('Creating emcee sampler...')
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
# burn-in
print('Burn-in...')
pos, prob, state = sampler.run_mcmc(p0, 10)
sampler.reset()
# sample 10 * 500 = 5000
print('Production...')
sampler.run_mcmc(pos, 1000)
# Save samples back to pymc model
print('Copying back into PyMC model')
model = pymc.MCMC(model)
model.sample(1) # This call is to set up the chains
for i, var in enumerate(model.stochastics):
    var.trace._trace[0] = sampler.flatchain[:, i]
pymc.Matplot.plot(model)

# Fit the maximum a posteriori (MAP) estimate
#map_fit = assay.experiment.map_fit()
#map_fit = None

# Run some MCMC sampling and return the MCMC object
trace = assay.experiment.run_mcmc(map_fit=map_fit)

# Show summary
#assay.experiment.show_summary(mcmc, map_fit)

# Generate plots
#plots_filename = 'plots.pdf'
assay.experiment.generate_plots(trace, pdf_filename=plots_filename)
