import numpy as np
from scipy import optimize
from assaytools import pymcmodels
from assaytools.bindingmodels import TwoComponentBindingModel


class AssaySimulator(object):
    """
    Class to predict fluorescence data using an AssayTools pymc model. Upon initialization of this class, fitted
    fluorescence parameters from assaytools are generated.

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
            Dictionary of assaytools pymc variables that have been sampled with MCMC.
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

    def set_mean_parameters(self, pymc_data, t_equil=0):
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
            Fmodel = self.F_PL * PL + self.F_L * L_free + self.F_P * P_free + self.F_buffer * self.path_length + self.F_plate
            if noisy:
                Fmodel += np.random.normal(loc=0.0, scale=self.sigma, size=len(self.l_total))

        return Fmodel

    def generate_deltag_estimates(self, nsamples=100):
        """
        Generate samples of fitter DeltaGs. This function draws multiple sets of fluorescence titration data, each with
        different amounts of random noise added. For draw,the binding free energy is estimated using least squares
        regression.

        Parameters
        ----------
        nsamples: int
            The number of fluorescence data draws

        Returns
        -------
        estimates: numpy.ndarray
            Binding free energies estimates, where each has been fitted using least squares regression
        """
        estimates = np.zeros(nsamples)
        for sample in range(nsamples):
            # The fluorescence data that will be fit to. Random noise is added the fluorescence data with noisy=True.
            target = self.simulate_fluorescence(noisy=True)

            def sum_of_squares(g):
                """
                The sum of squares between model fluorescence and the target
                """
                guess = self.simulate_fluorescence(DeltaG=g, noisy=False)
                return np.sum((guess - target) ** 2)

            # The initial guess within about 10% of the "true" value
            initial_guess = self.DeltaG + np.random.normal(loc=0, scale=0.1 * np.abs(self.DeltaG))
            fit = optimize.minimize(sum_of_squares, initial_guess, method='BFGS')
            estimates[sample] = fit.x[0]

        return estimates


def predict_assay_error(pymc_data, l_total, p_total, DeltaG=None, nsamples=10, nposterior_samples=100, t_equil=0, **kwargs):
    """
    Function to predict the expected coefficient of variation, variance, and squared bias of estimated binding free
    energies from assaytools PyMC data at different protein and ligand concentrations.

    Parameters
    ----------
    pymc_data: dict
        Dictionary of assaytools pymc variables that have been sampled with MCMC.
    DeltaG: float
        The binding free energy in thermal units.
    l_total: float or numpy.ndarray
        Array of ligand concentrations in M.
    p_total: float or numpy.ndarray
        Array of protein concentrations in M.
    nsamples: int
        The number of times fluorescence data will be drawn at all ligand concentrations for a given PyMC posterior
        sample.
    nposterior_samples:
        The number of assaytools posterior samples to draw. For each posterior sample, nsamples of fluorescence
        titration profiles will be drawn.
    t_equil: int
        The index from which it the assaytools PyMC MCMC samples are equilibrated. i.e. the end of the burnin-period.

    Returns
    -------
    (cv, var, bias): (np.ndarray, np.ndarray, np.ndarray)
        the mean coefficien of variation (cv), variance (var) and bias of the fitted free energy for each protein
        concentration.

    Example
    -------
    Run the assaytools PyMC pipeline and save the MCMC samples, here called 'pymc_data'. It's worth discarding the
    initial samples where the MCMC appears unequilibrated. The end of the burn-in time can be estimated with
    pymbar.timeseries by looking at trace of the free energy:
    >>> from pymbar import timeseries
    >>> (t_equil, g, N_eff) = timeseries.detectEquilibration(pymc_data['DeltaG'][0], fast=True, nskip=1)
    We'll be using t_equil below.

    We can now estimate how the protein concentration used in the assay will affect the coefficient of variation,
    variance, and bias of the fitted free energy.

    First, choose the ligand concentrations that can be used in the assay:
    >>> L_total = 10 ** (np.linspace(-10, -5, num=12))

    Next, pick the range of protein concentrations for which the error of the binding free energy will be estimated:
    >>> P_totals =10 ** (np.linspace(-10, -1, num=20))

    Estimating the variance error meaures for the above concentrations. Note how we're inputing the burn-in time, as
    well parsing the assay parameters.
    >>> (CV, Var, Bias2) = predict_assay_error(pymc_data, L_total, P_totals, t_equil=t_equil, geometry='top',
    >>>                                        assay_volume=100E-6, well_area=0.3969)

    Finally, it's a question of picking the protein concentration that reduces whichever metric you prefer.
    """
    # Pre-assigning the coefficient of variation, the bias, and variance respectively.
    CV = []
    Bias2 = []
    Var = []
    
    delta_gs = pymc_data['DeltaG'][0]           # All of the binding free energy samples
    pymc_indices = range(t_equil, len(delta_gs))     # The indices of the PyMC samples that are after the burn-in period:
    
    if DeltaG is None:
        DeltaG = np.mean(delta_gs)            # Used to determine the coefficient of variation.
        
    mean_delta_g = DeltaG

    # Looping over protein concentrations, drawing from the posterior and fitting the affinity:
    for p in range(len(p_total)):
        estimates = []
        bias_squared = []
        # Draw parameters from the posterior
        for i in range(nposterior_samples):
            ind = np.random.choice(pymc_indices,1)[0]
            simulator = AssaySimulator(pymc_data=pymc_data, l_total=l_total, sample_index=ind,  p_total=p_total[p], **kwargs)
            simulator.DeltaG = DeltaG
            # Generate fitted DeltaG estimates with stochastic noise added before each fit.
            estimates_per_posterior_sample = simulator.generate_deltag_estimates(nsamples=nsamples)
            bias_squared.append(np.mean((estimates_per_posterior_sample - simulator.DeltaG)**2))
            estimates.append(estimates_per_posterior_sample)
        # Collect the bias, variance, and coefficient of variantion
        bias_squared = np.array(bias_squared)
        Bias2.append(np.mean(bias_squared))
        estimates = np.array(estimates)
        Var.append(np.var(estimates))
        CV.append(np.std(estimates)/np.abs(mean_delta_g) * 100)

    return np.array(CV), np.array(Var), np.array(Bias2)
