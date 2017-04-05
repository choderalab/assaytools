import numpy as np
from scipy import optimize
from assaytools import pymcmodels
from assaytools.bindingmodels import TwoComponentBindingModel


class AssaySimulator(object):
    """
    Class to predict fluorescence data using an AssayTools pymc model.

    Example
    -------
    Decide what concentrations of the ligand (in M) you would like to predict the fluorescence for. For instance:
    >>> l_total = 10 ** (np.linspace(-10, -5, num=12))

    Given PyMC simulation data from assaytools, initialize the class
    >>> sim_assay = AssaySimulator(pymc_data=data, l_total=l_total, inner_filter=False)

    By default, if no sample_index is specified, assaytools parameters will be drawn from the PyMC posterior. If you'd
    like the use the assaytools expectation values, then
    >>> sim_assay.set_mean_parameters()

    Predict the fluorescence given the assaytools parameters.
    >>> fluorescence = simulate_fluorescence()

    You can see how the fluorescence is expected to be at specific protein concentrations (in M) with
    >>>  fluorescence = simulate_fluorescence(p_total=10E-6)

    """
    def __init__(self, pymc_data, l_total, p_total=None, sample_index=None, inner_filter=True, geometry='top',
                 assay_volume=100E-6, well_area=0.3969):
        """
        Parameters
        ----------
        pymc_data: dict
            Dictionary of pymc variables that have been simulated.
        l_total: numpy.ndarray
            Array of ligand concentrations in M.
        sample_index: int
            The iteration number from which pymc model parameters will be selected. If None, a random draw will be made.
        p_total: float
            Concentration of the protein in M.
        inner_filter: bool
            Whether to account for the inner filter effect.
        geometry: str
            Where the fluorescence is measured. Either 'top' or 'bottom'.
        assay_volume: float
            The volume of each well.
        well_area: float
            The surface area of each well that is visible to the detector.

        """
        # Non pymc paramters
        self.path_length = assay_volume * 1000 / well_area
        self.inner_filter = inner_filter
        self.geometry = geometry

        # Unpack the pymc parameters for a given MCMC iteration index.
        if sample_index is None:
            sample_index = np.random.choice(len(pymc_data['DeltaG'][0]), size=1)

        self.DeltaG = pymc_data['DeltaG'][0][sample_index]
        self.l_total = l_total
        if p_total is None:
            self.p_total = pymc_data["Ptrue"][0][sample_index,:].mean()
        else:
            self.p_total = p_total

        # Pymc fluorescence parameters
        self.F_buffer = pymc_data['F_buffer'][0][sample_index]
        self.F_plate = pymc_data['F_plate'][0][sample_index]
        self.F_PL = pymc_data['F_PL'][0][sample_index]
        self.F_L = pymc_data['F_L'][0][sample_index]
        self.F_P = pymc_data['F_P'][0][sample_index]

        if self.inner_filter:
            epsilon_em = pymc_data['epsilon_em'][0][sample_index]
            epsilon_ex = pymc_data['epsilon_ex'][0][sample_index]
            self.IF_i = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, self.path_length, self.l_total, self.geometry)
            self.IF_i_plate = np.exp(-epsilon_ex * self.path_length * self.l_total)

        # Select the parameters for the noise based on where the fluorescence was taken
        if self.geometry == 'top':
            self.sigma = pymc_data['sigma_top'][0][sample_index]
        else:
            self.sigma = pymc_data['sigma_bottom'][0][sample_index]

    def set_mean_parameters(self, pymc_data, t_equil=None):
        """
        Set the assaytools parameters to their expectation values, with the exception of p_total and l_total.

        Parameters
        -----------
        pymc_data: dict
            Dictionary of assaytools pymc variables that have been simulated.
        t_equil: int
            The iteration number from which the pymc MCMC samples are considered "equilibrated". i.e. the end of the
            burn-in period.
        """
        if t_equil is None:
            t_equil = 0

        self.DeltaG = np.mean(pymc_data['DeltaG'][0][t_equil:])

        # Pymc fluorescence parameters
        self.F_buffer = np.mean(pymc_data['F_buffer'][0][t_equil:])
        self.F_plate = np.mean(pymc_data['F_plate'][0][t_equil:])
        self.F_PL = np.mean(pymc_data['F_PL'][0][t_equil:])
        self.F_L = np.mean(pymc_data['F_L'][0][t_equil:])
        self.F_P = np.mean(pymc_data['F_P'][0][t_equil:])

        if self.inner_filter:
            epsilon_em = np.mean(pymc_data['epsilon_em'][0][t_equil:])
            epsilon_ex = np.mean(pymc_data['epsilon_ex'][0][t_equil:])
            self.IF_i = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, self.path_length, self.l_total, self.geometry)
            self.IF_i_plate = np.exp(-epsilon_ex * self.path_length * self.l_total)

        # Select the parameters for the noise based on where the fluorescence was taken
        if self.geometry == 'top':
            self.sigma = np.mean(pymc_data['sigma_top'][0][t_equil:])
        else:
            self.sigma = np.mean(pymc_data['sigma_bottom'][0][t_equil:])

    def simulate_fluorescence(self, DeltaG=None, p_total=None, noisy=True):
        """
        Predict the fluorescence of the complex using AssayTools posterior for a given protein concentration
        and range of ligand concentrations.

        Parameters
        ----------
        DeltaG: float
            The binding free energy in thermal units.
        p_total: float
            Concentration of protein in M.
        noisy: bool
            Whether to add detector noise to the predicted fluorescence.

        Returns
        --------
        Fmodel: numpy.ndarray
            The expected fluorescence given the supplied model parameters
        """
        if p_total is None:
            p_total = self.p_total
        if DeltaG is None:
            DeltaG = self.DeltaG

        # Predict the concentrations of the complex, free protein and free ligand
        [P_free, L_free, PL] = TwoComponentBindingModel.equilibrium_concentrations(DeltaG, p_total, self.l_total)

        # Model the fluorescence
        if self.inner_filter:
            Fmodel = self.IF_i * (self.F_PL * PL + self.F_L * L_free + self.F_P * P_free + self.F_buffer * self.path_length) + self.IF_i_plate * self.F_plate
            if noisy:
                Fmodel += np.random.normal(loc=0.0, scale=self.sigma, size=len(self.l_total))
        else:
            Fmodel = self.F_PL * PL + self.F_L * L_free + self.F_P * P_free + self.F_buffer * self.path_length
            if noisy:
                Fmodel += np.random.normal(loc=0.0, scale=self.sigma, size=len(self.l_total))

        return Fmodel

    def fit_deltaG(self, p_total=None):
        """
        Estimate the binding free energy given the fluorescence model parameters using least-squares fitting.

        Parameters
        ----------
        p_total: float
            Concentration of protein in M.

        Returns
        -------
        fit: float
            The estimated binding free energy in thermal units
        """

        if p_total is None:
            p_total = self.p_total

        # The fluorescence data that will be fit to
        target = self.simulate_fluorescence(p_total)

        def sum_of_squares(DeltaG, target=target):
            """
            The sum of squares between model fluorescence and the target
            """
            model = self.simulate_fluorescence(DeltaG, p_total, noisy=False)
            return np.sum((model - target)**2)

        # Start the initial guess within about 10% of the "true" value
        guess = self.DeltaG + np.random.normal(loc=0, scale=0.1 * np.abs(self.DeltaG))

        fit = optimize.minimize(sum_of_squares, guess, method='BFGS')

        return fit.x[0]
