"""

numpyro models for analysis of fluorescence assay data

"""

# =============================================================================================
# IMPORTS
# =============================================================================================

import abc
import numpy as np
import numpyro
import numpyro.distributions as dist
import arviz as az
import jax
import jax.numpy as jnp
from jax import random
from numpyro.infer import MCMC, NUTS, Predictive
from .numpyrobindingmodels import TwoComponentBindingModel
import matplotlib.pyplot as plt

# =============================================================================================
# Physical constants
# =============================================================================================

Na = 6.02214179e23  # Avogadro's number (number/mol)
kB = Na * 1.3806504e-23 / 4184.0  # Boltzmann constant (kcal/mol/K)
C0 = 1.0  # standard concentration (M)

# =============================================================================================
# Parameters for MCMC sampling
# =============================================================================================

DG_min = np.log(
    1e-15
)  # kT, most favorable (negative) binding free energy possible; 1 fM
DG_max = +0  # kT, least favorable binding free energy possible
niter = 500000  # number of iterations
nburn = 50000  # number of burn-in iterations to discard
nthin = 500  # thinning interval

# =============================================================================================
# numpyro submodels
# i'm assuming that we have LogNormal priors for now
# =============================================================================================

# construct a lognorm dist with proper loc and scale
def construct_lognorm(loc, scale):
    u = jnp.log(loc ** 2 / jnp.sqrt(loc ** 2 + scale ** 2))
    sig = jnp.log(1 + scale ** 2 / loc ** 2)
    return dist.LogNormal(loc=u, scale=sig)


# dispense an amount of liquid with expected concentration mu and variance var
def dispense(mu, var, name=""):
    return numpyro.sample(f"dispense_{name}", construct_lognorm(loc=mu, scale=var))


# sample a unknown value e.g. extinction coefficient or quantum yield
def hidden(min=0, max=10e6, name=""):
    return numpyro.sample(name, dist.Uniform(low=min, high=max))


# =============================================================================================
# numpyro base modules
# =============================================================================================


class MCMCModel:
    def __init__(
        self,
        model,
        mcmc_args={"num_warmup": nburn, "num_samples": niter, "thinning": nthin},
    ):
        self.model = model
        self.mcmc = None
        self.mcmc_args = mcmc_args
        self.params = None
        rng_key = random.PRNGKey(0)  # TODO: make option to choose random seed
        self.rng_key_infer, self.rng_key_predict = random.split(rng_key)

    def run_mcmc(self, *args, **kwargs):
        nuts_kernel = NUTS(self.model)
        self.mcmc = MCMC(nuts_kernel, **self.mcmc_args)
        self.mcmc.run(self.rng_key_infer, *args, **kwargs)
        self.sample_params()

    def sample_params(self):
        self.params = self.mcmc.get_samples()

    def predict(self, *args, **kwargs):
        predictor = Predictive(self.model, self.params)
        return predictor(self.rng_key_predict, *args, **kwargs)

    def plot_results(self):
        self.mcmc.print_summary()
        data = az.from_numpyro(self.mcmc)
        az.plot_trace(data, compact=True)


# =============================================================================================
# numpyro models
# =============================================================================================

# toy model for fluorescence curve
# assume we have precise solutions of each concentration of protein/ligand available
# and only the complex fluoresces
def toy_model(Pstated, dPstated, Lstated, dLstated, fluorescence=None):
    # binding free energy (kT), using a Uniform prior
    dG = numpyro.sample("dG", dist.Uniform(DG_min, DG_max))

    # we use our stated concentrations as priors
    # we use plate as each dispense should be conditionally independent
    with numpyro.plate("dispense_plate", len(Pstated)):
        Ptrue = dispense(Pstated, dPstated, name="Ptrue")
        Ltrue = dispense(Lstated, dLstated, name="Ltrue")

    # compute equilibrium concentrations using model
    [P_i, L_i, PL_i] = TwoComponentBindingModel.equilibrium_concentrations(
        dG, Ptrue, Ltrue
    )

    # scale/noise factors
    # all the bounds are sorta arbitrary
    f_var = hidden(
        max=100, name="f_measure_var"
    )  # variance for fluorescence measurement
    f_gain = hidden(max=1e13, name="f_gain")  # gain / quantum yield
    f_background = hidden(max=100, name="f_background")  # baseline signal

    # assume that only the complex fluoresces
    # e.g. FRET experiment at a particular wavelength
    F_PL_i = f_background + f_gain * PL_i

    # assume each measurement also has conditionally independent error
    with numpyro.plate("measure_plate", len(F_PL_i)):
        measurement = numpyro.sample(
            "measure_fluorescence",
            dist.Normal(loc=F_PL_i, scale=f_var),
            obs=fluorescence,
        )

    return measurement
