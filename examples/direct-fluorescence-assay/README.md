# examples/probe-assay

Here is a series of ipython notebooks that walks you through simulation and analysis of experimental binding data using the direct fluorescence of compounds that increase fluorescence upon binding to their target protein.

 - `1 Simulating Experimental Fluorescence Binding Data.ipynb`
 - `2 MLE fit for two component binding - simulated and real data.ipynb`
 - `3a Bayesian fit txt file - SrcBosutinib.ipynb`
 - `3b Bayesian fit xml file - SrcGefitinib.ipynb`

To plot raw data using `xml2png` navigate into the data folder and run:
 `xml2png --type singlet_384 p38*.xml`
 or
 `xml2png --type singlet_96 Gef*.xml`



