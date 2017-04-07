from assaytools.assay_simulator import AssaySimulator

def TestAssaySimulator(object):
    """
    A collection of functions to simulate the AssayTools simulator.
    """
    def test_initialization(model_data):
        """
        Test the initialization of the class

        pymc_data: dict
            Dictionary of pymc variables that have been simulated.

        Returns
        -------
        simulator: assaytools.AssaySimulator
            Class to sample from posterior in Assay Tools.
        """

        L_total = 10 ** (np.linspace(-10, -5, num=12))
        return AssaySimulator(pymc_data=model_data, L_total=L_total, sample_index=6000)
