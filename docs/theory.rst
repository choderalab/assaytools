.. _theory:

******
Theory
******

This section describes the theory behind the various Bayesian model fitting schemes in `AssayTools`.

Bayesian analysis
=================

`AssayTools` uses `Bayesian inference <https://en.wikipedia.org/wiki/Bayesian_inference>`_ to infer unknown parameters (such as ligand binding free energies) from experimental spectroscopic data.
It does this to allow the complete uncertainty---in the form of the joint distribution of all unknown parameters---to be rigorously characterized.
The most common way to summarize these results is generally to extract confidence intervals of individual parameters, but much more sophisticated analyses---such as examining the correlation between parameters---are also possible.

The Bayesian analysis scheme is intended to be modular, and the user can select whether certain effects (such as :ref:`inner filter effects <inner-filter-effects>`) are incorporated into the model.
Below, we describe the components of the model.
If an effect that carries unknown nuisance parameters (such as extinction coefficients for :ref:`inner filter effects <inner-filter-effects>`), these nuisance parameters carry additional prior terms along with them and are inferred as part of the inference procedure.

Unknown parameters
------------------
.. _parameters:

For convenience, we define the unknown parameters in the model for reference:

* the *true total ligand concentration* :math:`L_{true}` in the well (including all species involving the ligand)
* the *true receptor concentration* :math:`R_{true}` (including all species involving the receptor)

Data
----
.. _data:

For each experiment, the data is given by a set of observations for each well.
Each well is associated with some properties:

* the *volume* :math:`V` of sample in the well (mostly buffer)
* a total concentration :math:`[R]_T` of *receptor* added to the well
* a total concentration :math:`[L]_T` of *ligand* added to the well (or potentially multiple ligands)
* the *well area* :math:`A` with the assumption that the well is cylindrical (allowing the path length :math:`l` to be computed)

and one or more experimental measurements:

* a *top fluorescence measurement* (returning toward the incident excitation light) :math:`F_{top}`
* a *bottom fluorescence measurement* (proceeding in the same direction as the incident excitation light) :math:`F_{bottom}`
* an *absorbance measurement* :math:`A`

Priors
------
.. _priors:

Each unknown parameter in the model is assigned a *prior* distribution that reflects the state of knowledge we have of its value before including the effects of the observed data.

Concentrations
^^^^^^^^^^^^^^

While we design the experiment to dispense the *intended* amount of protein and ligand into each well, the true amount dispensed into the well will vary due to random pipetting error.
The *true* concentrations of protein :math:`R_{true}` and ligand :math:`L_{true}` in each well are therefore unknown.
Because we propagate the pipetting error along with the intended concentrations, we have the intended ("stated") protein concentration :math:`P_{stated}` and its standard error :math:`\delta P_{stated}` as input.
Similarly, the stated ligand concentration :math:`L_{stated}` and its error :math:`\delta L_{stated}` are also known.

If we assume the dispensing process is free of bias, the simplest distribution that fits the stated concentration and its standard deviation without making additional assumptions is a Gaussian.

We assign these true concentrations for the receptor :math:`R_{true}` and ligand :math:`L_{true}` a prior distribution.
If ``concentration_priors`` is set to ``gaussian``, this is precisely what is used

.. math::

   R_{true} &\sim N(R_{stated}, \delta R_{stated}) \\
   L_{true} &\sim N(L_{stated}, \delta L_{stated}) \\

This is expressed in the ``pymc`` model as ::

  Ptrue = pymc.Normal('Ptrue', mu=Pstated, tau=dPstated**(-2)) # protein concentration (M)
  Ltrue = pymc.Normal('Ltrue', mu=Lstated, tau=dLstated**(-2)) # ligand concentration (M)

.. note:: ``pymc`` uses the *precision* :math:`\tau \equiv \sigma^{-2}` instead of the variance :math:`\sigma^2` as a parameter of the normal distribution.

Gaussian priors have the unfortunate drawback that there is a small but nonzero probability that these concentrations would be negative, leading to nonsensical (unphysical) concentrations.
To avoid this, we generally use a lognormal distribution (selected by ``concentration_priors='lognormal'``).
The lognormal priors are expressed in the ``pymc`` model as ::

  Ptrue = pymc.Lognormal('Ptrue', mu=np.log(Pstated**2 / np.sqrt(dPstated**2 + Pstated**2)), tau=np.sqrt(np.log(1.0 + (dPstated/Pstated)**2))**(-2)) # protein concentration (M)
  Ltrue = pymc.Lognormal('Ltrue', mu=np.log(Lstated**2 / np.sqrt(dLstated**2 + Lstated**2)), tau=np.sqrt(np.log(1.0 + (dLstated/Lstated)**2))**(-2)) # ligand concentration (M)

.. note:: The parameters of a *lognormal distribution* differ from those of a normal distribution by the relationship `described here <https://en.wikipedia.org/wiki/Log-normal_distribution>`_. The parameters above ensure that the mean concentration is the stated concentration and the standard deviation is its experimental uncertainty.  The relationship between the mean and variance of the normal distribution :math:`\mu_N, \sigma_N^2` and the parameters for the lognormal distribution is given by:
.. math::

   \mu_{LN} &= \ln \frac{\mu_N^2}{\sqrt{\mu_N^2 + \sigma_N^2}} \\
   \sigma^2_{LN} &= \ln \left[ 1 + \left( \frac{\sigma_N}{\mu_N}\right)^2 \right] \\
   \tau_{LN} &= \ln \left[ 1 + \left( \frac{\sigma_N}{\mu_N}\right)^2 \right]^{-1}


Binding free energy
^^^^^^^^^^^^^^^^^^^

The ligand binding free energy :math:`\Delta G` is unknown, and presumed to either be unknown over a large uniform range with the ``uniform`` prior

.. math::

   \Delta G \sim U(-\Delta G_{min}, +\Delta G_{max})

This is expressed in the ``pymc`` model as ::

  DeltaG = pymc.Uniform('DeltaG', lower=DG_min, upper=DG_max) # binding free energy (kT), uniform over huge range

This uniform prior has the drawback that affinities near the extreme measurable ranges are simply unknown with equal likelihood out to absurd extreme values.

We can attenuate the posterior probabilities at extreme affinities by using a prior inspired by the range of data recorded in `ChEMBL <https://www.ebi.ac.uk/chembl/>`_ via the ``chembl`` prior, with a Gaussian form

.. math::

   \Delta G &\sim N(0, \sigma^2) \\
   \sigma &= 12.5 \: \mathrm{kcal/mol}

This is expressed in the ``pymc`` model as ::

  DeltaG = pymc.Normal('DeltaG', mu=0, tau=1./(12.5**2)) # binding free energy (kT), using a Gaussian prior inspired by ChEMBL

Components
----------

We now discuss the various modular components of the Bayesian inference scheme.

Fluorescence
^^^^^^^^^^^^
.. _fluorescence:

Fluorescence can be measured from either the top (from which the plate is illuminated in the Tecan Infinite M1000PRO), bottom, or both.
Observed fluorescence depends on the concentration of each species :math:`X_i`, the quantum efficiencies of each species at the excitation/emission wavelengths :math:`q_i(ex,em)`, and the intrinsic background fluorescence of the plate :math:`F_{plate}`

.. math::

   F_{top} = I_0 \left[ \sum_{i \in components} q_i(ex,em) [X_i] + F_{plate,top} \right]

   F_{bottom} = I_0 \left[ \sum_{i \in components} q_i(ex,em) [X_i] + F_{plate,bottom} \right]

Here, :math:`I_0` is the incident excitation intensity.

Absorbance
^^^^^^^^^^
.. _absorbance:

The absorbance is determined by the the extinction coefficient of each component (`R`, `L`, `RL` for simple two-component binding) at the excitation wavelength, as well as any intrinsic absorbance of the plate at that wavelength.

.. math::

   A = \sum_{i \in components} \epsilon_{ex,i} l [X_i] + A_{plate}

Inner filter effects
^^^^^^^^^^^^^^^^^^^^
.. _inner-filter-effects:

At high ligand concentrations, if the ligand has significant absorbance at the excitation wavelength, the amount of light reaching the bottom of the well is less than the light reaching the top of the well.
This is called the *primary inner filter effect*, and has the net effect of attenuating the observed fluorescence.
Similarly, the *secondary inner filter effect* is caused by significant absorbance at the emission wavelength.
When both effects are combined, the net attenuation effect depends on the geometry of excitation and detection:

.. note:: Add figure illustrating inner filter effects.

.. note:: Derive attenuation factors.

Measured extinction coefficients
""""""""""""""""""""""""""""""""

If the extinction coefficients have been measured, we have a measurement :math:`\epsilon` and corresponding standard error :math:`\delta \epsilon` available.
Because extinction coefficients must be positive, we use a lognormal distirbution to model the true extinction coefficients about the measured value

.. math::

  \epsilon \sim \mathrm{LN}(\mu, \tau) \\
  \mu = \ln \frac{\epsilon^2}{\sqrt{\epsilon^2 + \delta \epsilon^2}} \\
  \tau = \ln \left[ 1 + \left( \frac{\delta \epsilon}{\epsilon}\right)^2 \right]^{-1}

This is modeled in the code as ::

  model['epsilon_ex'] = pymc.Lognormal('epsilon_ex', mu=np.log(epsilon_ex**2 / np.sqrt(depsilon_ex**2 + epsilon_ex**2)), tau=np.sqrt(np.log(1.0 + (depsilon_ex/epsilon_ex)**2))**(-2)) # prior is centered on measured extinction coefficient
  model['epsilon_em'] = pymc.Lognormal('epsilon_em', mu=np.log(epsilon_em**2 / np.sqrt(depsilon_em**2 + epsilon_em**2)), tau=np.sqrt(np.log(1.0 + (depsilon_em/epsilon_em)**2))**(-2)) # prior is centered on measured extinction coefficient

Inferred extinction coefficients
""""""""""""""""""""""""""""""""

If the extinction coefficients

Binding models
==============

`AssayTools` has a variety of binding models implemented.
Though the user must currently specify the model to be fit to the data, we plan to include the ability to automatically select the most appropriate binding model automatically using `reversible-jump Monte Carlo (RJMC) <https://en.wikipedia.org/wiki/Reversible-jump_Markov_chain_Monte_Carlo>`_, which also permits `Bayesian hypothesis testing <https://en.wikipedia.org/wiki/Bayes_factor>`_.
All binding models are subclasses of the :class:`BindingModel` abstract base class, and users can implement their own binding models as subclasses.

Two-component binding model
---------------------------

A two-component binding model is implemented in :class:`assaytools.bindingmodels.TwoComponentBinding`.
When it is known that receptor `R` associates with ligand `L` in a 1:1 fashion, we can write the dissociation constant :math:`K_d` in terms of the equilibrium concentrations of each species:

.. math::

   K_d = \frac{[R][L]}{[RL]}

Incorporating conservation of mass constraints

.. math::

   [R]_T &= [R] + [RL] \\
   [L]_T &= [L] + [RL]

we can eliminate the unknown concentrations of free receptor :math:`[R]` and free ligand :math:`[L]` to obtain an expression for the complex concentration :math:`[RL]` in terms of fixed quantities (dissociation constant :math:`K_d` and total concentrations :math:`[R]_T` and :math:`[L]_T`):

.. math::

   K_d = \frac{([R]_T - [RL]) ([L]_T - [RL])}{[RL]}

   [RL] K_d = ([R]_T - [RL]) ([L]_T - [RL])

   0 = [RL]^2 - ([R]_T + [L]_T + K_d) [RL] + [R]_T [L]_T

This quadratic equation has closed-form solution, with only one branch of the solution giving :math:`0 < [RL] < \min([R]_T, [L]_t)`:

.. math::

   K_d = \frac{1}{2} \left[ ([R]_T + [L]_T + K_d) - \sqrt{([R]_T + [L]_T + K_d)^2 - 4 [R]_T [L]_T} \right]

Note that this form is not always numerically stable since :math:`[R]_T`, :math:`[L]_T`, and :math:`K_d` may differ by orders of magnitude, leading to slightly negative numbers inside the square-root.
`AssayTools` uses the logarithms of these quantities instead, and guards against negative values inside the square root.

Competitive binding model
-------------------------

When working with N ligands :math:`L_n` that bind a single receptor :math:`R`, we utilize a competitive binding model implemented in :class:`assaytools.bindingmodels.CompetitiveBindingModel`.
Here, the dissociation constants :math:`K_n` are defined as

.. math::

   K_n = \frac{[R][L_n]}{[RL_n]}

with corresponding conservation of mass constraints

.. math::

   [R]_T &= [R] + \sum_{n=1}^N [RL_n] \\
   [L_n]_T &= [L_n] + [RL_n], n = 1,\ldots, N

The solution must also satisfy some constraints:

.. math::

   0 \le [RL_n] \le \min([L_n], [R]_T) \:\:,\:\: n = 1,\ldots,N

   \sum_{n=1}^N [RL_n] \le [R]_T

We can rearrange these expressions to give

.. math::

   [R][L_n] - [RL_n] K_n = 0  \:\:,\:\: n = 1, \ldots, N

and eliminate :math:`[RL_n]` and :math:`[R]` to give

.. math::

   \left( [R]_T - \sum_{n=1}^N [RL_n] \right) * ([L_n]_T - [RL_n]) - [RL_n] K_n = 0  \:\:,\:\: n = 1, \ldots, N

This leads to a coupled series of equations that cannot easily be solved in closed form, but are straightforward to solve numerically using the solver :func:`scipy.optimize.fsolve`, starting from an initial guess that ensures the constraints remain satisfied.
