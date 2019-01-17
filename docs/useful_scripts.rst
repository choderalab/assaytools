**************
Useful Scripts
**************

.. highlight:: bash

xml2png
=======

Converts xml data file output from the Tecan Infinite M1000 Pro plate reader to png plot of fluorescence and absorbance values. It allows for the quick visual inspection of raw experimental results.

::

    $ xml2png *.xml --singlet 'singlet_96'


quickmodel
==========

Builds quick Bayesian model of both spectra and single wavelength two component binding experiments.

As input, it requires xml output files of the experiment form plate reader and a python script(`inputs.py`) that includes all experimental design details.

1. Run `calculate_L_stated_array` to generate ligand concentration array and copy it into `inputs.py`.
2. Construct `inputs.py` script based on experimental design.
3. Run`quickmodel`.

::

    $ quickmodel --inputs 'inputs' --type 'singlet' --nsamples 10000


inputs.py
---------

`inputs.py` should be manually constructed to record experimental design details, following the layout of of `inputs_example.py` file.
Ligand concentration array (`Lstated` section) can be constructed using `calculate_L_stated_array.py` script.

**Sections of `inputs.py`**

- `xml_file_path` : relative path to xml plate reader output files.
- `file_set` : option to group multiple experimental sets with a dictionary key.
- `ligand_order` : List of ligand names per each experiment set (one protein, one buffer row). If `None` Python object is specified as ligand name in this list, `quickmodel` analysis will skip the analysis of that experiment. For example:
::

    'ligand_order'  :  [None, None, 'ligand3', 'ligand4']


- `section` : Data section label of Tecan Infinite M1000 Pro plate reader as specified in its method.
- `wavelength` : Emission wavelength picked for analysis (nm).
- `Lstated` : Experimental value of ligand concentration (M), in NumPy array form. It can be constructed using `calculate_L_stated_array.py` script.
- `L_error` : Estimated % error in stated ligand concentrations.
- `Pstated` : Experimental value of protein concentration (M).
- `P_error` : Estimated % error in stated protein concentration.
- `assay_volume` : Volume of each assay sample (L).
- `well_area` : Area of microtiter plate well (cm^2)


calculate_L_stated_array.py
---------------------------

Generates Numpy array of stated ligand concentration (Lstated) for logarithmic or linear dilution along a row, adjusted for true experimental ligand concentration. This numpy array is necessary to construct `inputs.py` file for `quickmodel.py` analysis.
Provide information on how ligand titration is constructed: number of wells in each titration (`--n_wells`), highest and lowest ligand concentrations in molar units (--h_conc and --l_conc), target ligand stock concentration and true ligand stock concentration in molar units (--target_stock_conc and --true_stock_conc), and serial dilution mode (`--dilution`, linear or logarithmic) as inputs.

::

    $ calculate_Lstated_array --n_wells 12 --h_conc 8e-06 --l_conc 2.53e-09 --target_stock_conc 0.010 --true_stock_conc 0.0100344 --dilution logarithmic

The numpy array this script prints out must be directly copied to `Lstated` section of `inputs.py`.




