# examples/direct-fluorescence-assay

Here is a series of ipython notebooks that walks you through simulation and analysis of experimental binding data using the direct fluorescence of compounds that increase fluorescence upon binding to their target protein.

 - `1 Simulating Experimental Fluorescence Binding Data.ipynb`
 - `2 MLE fit for two component binding - simulated and real data.ipynb`
 - `3a Bayesian fit txt file - SrcBosutinib.ipynb`
 - `3b Bayesian fit xml file - SrcGefitinib.ipynb`

## Command line tools

* `xml2png`

To plot raw data using `xml2png`, navigate into the data folder and run:

 `xml2png --type singlet_384 p38*.xml`

 or

 `xml2png --type singlet_96 Gef*.xml`

* `quickmodel`

To get a deltaG estimate from direct fluorescence data using our Bayesian modeling framework, via the `quickmodel` script run:

  `quickmodel --inputs 'inputs_p38_singlet'`

or
 
  `quickmodel --inputs 'inputs_Gef_WIP'`

or

  `quickmodel --inputs 'inputs_single_well'`

The three inputs.py files here are required to run these, but they can be changed if you want to run other types of analyses (different sections of the xml file, different wavelenghts, etc.).
 


