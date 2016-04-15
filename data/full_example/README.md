# full_example

----
## Files describing the Abl-Gefitinib Imatinib competition assay.

These are files describing two experiments:

* Abl-Gefitinib Imatinib competition assay performed on January 15, 2016.
* Abl spectra with Bosutinib, Bosutinib Isomer, Erlotinib, and Gefitinib performed on March 11, 2016.

Protein concentration data:

* `protein_data.txt`

Infinite data files:

* `Abl Gef Ima gain 120 bw1020 2016-01-19 16-22-45_plate_1.xml`- Abl Geftinib Imatinib Singlets (96-well)
* `Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml` - Abl Gefitnib Singlets (96-well)
* `Abl_D382N_Bos_20160311_132205.xml` - Rows A and B - Bosutinib Spectra Data for Abl
* `Abl_D382N_BosI_20160311_135952.xml` - Rows C and D - Bosutinib Isomer Spectra Data for Abl
* `Abl_D382N_Erl_20160311_143642.xml` - Rows E and F - Erlotinib Spectra Data for Abl
* `Abl_D382N_Gef_20160311_152340.xml` - Rows G and H - Gefitinib Spectra Data for Abl

D300 simulation reports:

* `Src_Bos_Ima_96well_Mar2015 2015-03-07 1736.DATA.xml` - D300 script used to dispense Gefitinib and Imatinib for Singlet Fluorescent and Singlet Imatinib Experiment
* `LRL_Src_Bos_2rows_1_2 2015-09-11 1048.DATA.xml` - D300 script used to dispense all ligands in Spectra assay

Compound stocks used:

* `DMSOstocks-Sheet1.csv`
   * For the competition assay GEF999 and IMA001 were used
   * For the spectra assay BOS001, BOI001, GEF001, and ERL001 were used.

## API NOTES

`example.py` contains a verbose example of using the new general `assaytoos.analysis.CompetitiveBindingAnalysis` API to analyze spectroscopic measurements using one or more components.

Using API looks like this:

```python
from assaytools.analysis import CompetitiveBindingAnalysis
# Define an experiment by specifying source solutions, WellGroup to analyze, and names of receptor and ligand components for competitive binding model.
experiment = CompetitiveBindingAnalysis(solutions=solutions, wells=container.all_wells(), receptor_name=['Abl'], ligand_names=['bosutinib', 'erlotinib', 'gefinitb', 'bosutinib isomer'])
# Determine the maximum a posteriori (MAP) estimate.
map_fit = experiment.map_fit()
# Run MCMC sampling and return the MCMC object.
mcmc = experiment.run_mcmc()
```

### Specifying solutions

`Solution` objects describe the main source solution stocks from which wells are filled.
`Solution` objects specify at most the concentration of a *single* component (such as protein or ligand) and its corresponding uncertainty.
Each `Solution` has a `name`, `concentration`, and concentration `uncertainty`.

There are three types of `Solution`:
* `BufferSolution`: A buffer that contains no compound species.
* `ProteinSolution`: A single protein species prepared spectrophotometrically from a high-concentration stock solution whose `absorbance` measurement was taken, with corresponding concentration calculated automatically from `extinction_coefficient` and `molecular_weight`. A quantity `ul_protein_stock` microliters of high-concentration stock solution is added to `ml_buffer` milliliters of buffer, with the source `BufferSolution` specified as the `buffer` attribute.
* `DMSOStockSolution`: A DMSO stock solution prepared graviemetrically. This is currently constructed from a `dict` object specifying a compound stock solution through a number of required parameters (`{'compound name', 'compound mass (mg)', 'molecular weight', 'purity', 'solvent mass (g)'}`)

The `solutions` argument of `CompetitiveBindingAnalysis` is a `dict` specifying a solution key associated with a `Solution` object, such as
```python
solutions = dict()
solutions['buffer'] = BufferSolution(name='20 mM Tris buffer')
solutions['Abl'] = ProteinSolution(name='1 uM Abl', buffer=solutions['buffer'], absorbance=4.24, extinction_coefficient=49850, molecular_weight=41293.2, ul_protein_stock=165.8, ml_buffer=14.0)
solutions['BOS'] = DMSOStockSolution(dmso_stocks['BOS001'])
```

### Specifying wells

The [`autoprotocol.container` API](http://pythonhosted.org/autoprotocol/_modules/autoprotocol/container.html) is used to create a `WellGroup` for analysis of any spectroscopic data available within the specified wells. If a whole `Container` is specified (as `plate`, for example), all wells can be selected with [`container.all_wells()`](http://pythonhosted.org/autoprotocol/autoprotocol.html#autoprotocol.container.Container.all_wells).

Each [`Well`](http://pythonhosted.org/autoprotocol/autoprotocol.html#autoprotocol.container.Well) has a `properties` attribute that is a `dict` allowing additional well properties to be specified. We use this to hold two important additional pieces of information:
* `contents`: A `dict` specifying the contents of each well, in terms of the volume of each `Solution` that was added
* `measurements`: A `dict` specifying any spectroscopic measurements made of the well, which need not be made for every well of the `WellGroup` analyzed

#### Specifying well contents

The `well.properties['contents']` is a `dict` specifying how much of each `Solution` was added to the well:
```python
from autoprotocol.unit import Unit
well.properties['contents'] = { 'Abl' : Unit(100, 'microliters'), 'BOS' : Unit(2, 'microliters')}
```
Solutions that were not added need to be specified.

#### Specifying well measurements

The `well.properties['measurements']` is a `dict` specifying which measurements, if any, were performed for the corresponding well.
The `analysis` framework supports arbitrary sets of measurements for individual wells---the same set of measurements need not be made across all wells in the specified `WellSet` used for analysis.

Supported measurements include:
* `absorbance` : the absorbance of a well at one or more wavelengths
* `fluorescence` : the fluorescence of a well at one or more excitation/emission wavelengths and geometries (`{top, bottom}`)

##### Absorbance measurements

Absorbance measurements are specified with the corresponding wavelength
```python
from autoprotocol.unit import Unit
well.properties['measurements']['absorbance'] = { '280:nanometers' : 0.437 }
```
***NOTE: `Unit` objects, such as `Unit(280, 'nanometers')`, cannot be used as not unique keys, so we are forced to use the string representations, such as '280:nanometers', as keys.***

##### Fluorescence measurements

Fluorescence measurements are specified with the corresponding excitation wavelength, emission wavelength, and detection geometry as a tuple:
```python
from autoprotocol.unit import Unit
well.properties['measurements']['fluorescence'] = {
    ('280:nanometers', '450:nanometers', 'top') : 12425,
    ('280:nanometers', '450:nanometers', 'bottom') : 1425
    }
```

#### Analysis

Well properties may therefore look like this:
```python
{'measurements': {'absorbance': {'480:nanometers': 0.0413, '350:nanometers': 0.0723, '280:nanometers': 0.639}, 'fluorescence': {('280:nanometers', '450:nanometers', 'bottom'): 17775.0, ('350:nanometers', '450:nanometers', 'top'): 89.0, ('280:nanometers', '450:nanometers', 'top'): 8671.0, ('350:nanometers', '450:nanometers', 'bottom'): 172.0}}, 'contents': {'buffer': Unit(100.0, 'microliter'), 'DMSO': Unit(300.0, 'nanoliter'), 'IMA001': Unit(100.0, 'nanoliter')}, 'area': Unit(31.1724531052, 'millimeter ** 2')}
```
