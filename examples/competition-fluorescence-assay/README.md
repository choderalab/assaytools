# examples/competition-fluorescence-assay

Here is a series of ipython notebooks that walks you through simulation and analysis of experimental binding data using non-fluorescent ligands that compete off a compound that increases fluorescence upon binding a target protein.

-`1a-modelling-CompetitiveBinding-FluorescenceBindingAssay.ipynb` - Predicting binding curves and fluorescence signal for competitive binding experiment of one fluorescent and one non-fluorescent ligand using standart experimental assumptions (no ligand depletion). HSA - dansyl amide - phenylbutazone is used as an example.

-`1b-modelling-CompetitiveBinding-AnalyticalSolution.ipynb` - Exact solution to 2 ligand competition equilibrium expressions is used to predict saturation curves and fluorescence signal. HSA - dansyl amide - phenylbutazone is used as an example.

-`2a-competition-assay-modeling.ipynb` - A notebook exploring how concentrations of a competitive ligand effect the fluorescence signal. ** This will be replaced. **

-`2b MLE fit for three component binding - simulated data.ipynb` - A notebook that generates simulated competition assay data and then analyzes it using the three component binding model and least squares.

-`3-Competition-Assay-Data-Plotting.ipynb` - A notebook that imports and plots competition assay data. ** This needs to be cleaned up and analysis added to it. **

## Command line tools 

To plot raw data using `xml2png` navigate into the data folder and run:

`xml2png --type singlet_96 *.xml`


