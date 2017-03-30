import numpy as np
from scipy import optimize
from assaytools import pymcmodels
from assaytools.bindingmodels import TwoComponentBindingModel


class AssaySimulator(object):
    """
    Class to simulate fluorescence assays using an AssayTools pymc model
    """
    def __init__(self, pymc_data, L_total, sample_index, P_total=None, inner_filter=True, geometry='top', assay_volume=100E-6, well_area=0.3969):
        """
        Parameters
        ----------
        pymc_data: dict
            Dictionary of pymc variables that have been simulated.
        L_total: numpy.ndarray
            Array of ligand concentrations in M.
        sample_index: int
            The iteration number from which pymc model parameters will be selected.
        P_total: float
            Concentration of protein in M.
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
        # Unpack the pymc parameters for a given MCMC iteration index.
        self.DeltaG = pymc_data['DeltaG'][0][sample_index]
        self.L_total = L_total
        if P_total is None:
            self.P_total = 1E-6
        else:
            self.P_total = P_total

        # Pymc fluorescence parameters
        self.F_buffer = pymc_data['F_buffer'][0][sample_index]
        self.F_plate = pymc_data['F_plate'][0][sample_index]
        self.F_PL = pymc_data['F_PL'][0][sample_index]
        self.F_L = pymc_data['F_L'][0][sample_index]
        self.F_P = pymc_data['F_P'][0][sample_index]

        if self.inner_filter:
            epsilon_em = pymc_data['epsilon_em'][0][sample_index]
            epsilon_ex = pymc_data['epsilon_ex'][0][sample_index]
            self.IF_i = pymcmodels.inner_filter_effect_attenuation(epsilon_ex, epsilon_em, self.path_length, self.L_total, geometry)
            self.IF_i_plate = np.exp(-epsilon_ex * self.path_length * self.L_total)

        # Select the parameters for the noise based on where the fluorescence was taken
        if geometry == 'top':
            self.sigma = pymc_data['sigma_top'][0][sample_index]
        else:
            self.sigma = pymc_data['sigma_bottom'][0][sample_index]

    def simulate_fluorescence(self, DeltaG=None, P_total=None, noisy=True):
        """
        Predict the fluorescence of the complex using AssayTools posterior for a given protein concentration
        and range of ligand concentrations.

        Parameters
        ----------
        DeltaG: float
            The binding free energy in thermal units.
        P_total: float
            Concentration of protein in M.
        noisy: bool
            Whether to add detector noise to the predicted fluorescence.

        Returns
        --------
        Fmodel: numpy.ndarray
            The expected fluorescence given the supplied model parameters
        """
        if P_total is None:
            P_total = self.P_total
        if DeltaG is None:
            DeltaG = self.DeltaG

        # Predict the concentrations of the complex, free protein and free ligand
        [P_free, L_free, PL] = TwoComponentBindingModel.equilibrium_concentrations(DeltaG, P_total, self.L_total)

        # Model the fluorescence
        if self.inner_filter:
            Fmodel = self.IF_i * (self.F_PL * PL + self.F_L * L_free + self.F_P * P_free + self.F_buffer * self.path_length) + self.IF_i_plate * self.F_plate
            if noisy:
                Fmodel += np.random.normal(loc=0.0, scale=self.sigma, size=len(self.L_total))
        else:
            Fmodel = self.F_PL * PL + self.F_L * L_free + self.F_P * P_free + self.F_buffer * self.path_length
            if noisy:
                Fmodel += np.random.normal(loc=0.0, scale=self.sigma, size=len(self.L_total))

        return Fmodel

    def fit_deltaG(self, P_total=None):
        """
        Estimate the binding free energy given the fluorescence model parameters using least-squares fitting.

        Parameters
        ----------
        P_total: float
            Concentration of protein in M.

        Returns
        -------
        fit: float
            The estimated binding free energy in thermal units
        """

        if P_total is None:
            P_total = self.P_total

        # The fluorescence data that will be fit to
        target = self.simulate_fluorescence(P_total)

        def sum_of_squares(DeltaG, target=target):
            """
            The sum of squares between model fluorescence and the target
            """
            model = self.simulate_fluorescence(DeltaG, P_total, noisy=False)
            return np.sum((model - target)**2)

        # Start the initial guess within about 10% of the "true" value
        guess = self.DeltaG + np.random.normal(loc=0, scale=0.1 * np.abs(self.DeltaG) )

        fit = optimize.minimize(sum_of_squares, guess, method='BFGS')

        return fit.x[0]