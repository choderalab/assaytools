.. _theory:

******
Theory
******

This section describes the theory behind the various Bayesian model fitting schemes in `AssayTools`.

Bayesian analysis
=================

`AssayTools` uses `Bayesian inference `https://en.wikipedia.org/wiki/Bayesian_inference>`_ to infer unknown parameters (such as ligand binding free energies) from experimental spectroscopic data.
It does this to allow the complete uncertainty---in the form of the joint distribution of all unknown parameters---to be rigorously characterized.
The most common way to summarize these results is generally to extract confidence intervals of individual parameters, but much more sophisticated analyses---such as examining the correlation between parameters---are also possible.

Unknown parameters
------------------
.. _parameters:

For convenience, we define the unknown parameters in the model for reference:
* `Ltrue` - true ligand concentration in the well
* `Ptrue` - true protein concentration in the well

Data
----
.. _data:

The data are given as...

.. math::

   (a + b)^2 = a^2 + 2ab + b^2

   (a - b)^2 = a^2 - 2ab + b^2

Priors
------
.. _priors:

Each unknown parameter in the model is assigned a *prior* distribution that reflects the state of knowledge we have of its value before including the effects of the observed data.

Concentrations
^^^^^^^^^^^^^^

While we design the experiment to dispense the *intended* amount of protein and ligand into each well, the true amount dispensed into the well will vary due to random pipetting error.
The *true* concentrations of protein (`Ptrue`) and ligand (`Ltrue`) in each well are therefore unknown.
Because we propagate the pipetting error along with the intended concentrations, we have the intended ("stated") protein concentration (`Pstated`) and its standard error (`dPstated`) as input.
Similarly, the stated ligand concentration (`Lstated`) and its error (`dLstated`) are also known.

If we assume the dispensing process is free of bias, the simplest distribution that fits the stated concentration and its standard deviation without making additional assumptions is a Gaussian.

We assign these true concentrations `Ptrue` (protein) and `Ltrue` (ligand) a prior distribution.
If `concentration_priors` is set to `gaussian`, this is precisely what is used ::

  Ptrue = pymc.Normal('Ptrue', mu=Pstated, tau=dPstated**(-2)) # protein concentration (M)
  Ltrue = pymc.Normal('Ltrue', mu=Lstated, tau=dLstated**(-2)) # ligand concentration (M)

.. note:: `pymc` uses the precision `tau = sigma**(-2)` instead of the standard deviation `sigma` as a parameter of the distribution.

Gaussian priors have the unfortunate drawback that there is a small but nonzero probability that these concentrations would be negative, leading to nonsensical (unphysical) concentrations.
To avoid this, we generally use a lognormal distribution (selected by `concentration_priors='lognormal'`.
The lognormal priors are defined as ::

  Ptrue = pymc.Lognormal('Ptrue', mu=np.log(Pstated**2 / np.sqrt(dPstated**2 + Pstated**2)), tau=np.sqrt(np.log(1.0 + (dPstated/Pstated)**2))**(-2)) # protein concentration (M)
  Ltrue = pymc.Lognormal('Ltrue', mu=np.log(Lstated**2 / np.sqrt(dLstated**2 + Lstated**2)), tau=np.sqrt(np.log(1.0 + (dLstated/Lstated)**2))**(-2)) # ligand concentration (M)

.. note:: The parameters of a lognormal distribution differ from those of a normal distribution by the relationship `described here <https://en.wikipedia.org/wiki/Log-normal_distribution>`_. The parameters above ensure that the mean concentration is the stated concentration and the standard deviation is its experimental uncertainty.

Binding free energy
^^^^^^^^^^^^^^^^^^^

The ligand binding free energy `DeltaG` is unknown.

Inner filter effects
^^^^^^^^^^^^^^^^^^^^

At high ligand concentrations, if the ligand has significant absorbance at the excitation wavelength, the amount of light reaching the bottom of the well is less than the light reaching the top of the well.
This is called the *primary inner filter effect*, and has the net effect of attenuating the observed fluorescence.
Similarly, the *secondary inner filter effect* is caused by significant absorbance at the emission wavelength.
When both effects are combined, the net attenuation effect depends on the geometry of excitation and detection:

If we are allowing for primary inner filter effects, in which incident excitation light is absorbed by the ligand, we use a lognormal distribution for the ligand extinction coefficient at the excitation wavelength `epsilon_ex` ::

  model['epsilon_ex'] = pymc.Lognormal('epsilon_ex', mu=np.log(epsilon_ex**2 / np.sqrt(depsilon_ex**2 + epsilon_ex**2)), tau=np.sqrt(np.log(1.0 + (depsilon_ex/epsilon_ex)**2))**(-2)) # prior is centered on measured extinction coefficient

If we are allowing for secondary inner filter effects, in which emission light is absorbed by the ligand, we use a lognormal distribution for the ligand extinction coefficient at the emission wavelength `epsilon_ex` ::

  model['epsilon_em'] = pymc.Lognormal('epsilon_em', mu=np.log(epsilon_em**2 / np.sqrt(depsilon_em**2 + epsilon_em**2)), tau=np.sqrt(np.log(1.0 + (depsilon_em/epsilon_em)**2))**(-2)) # prior is centered on measured extinction coefficient

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

This leads to a coupled series of equations that cannot easily be solved in closed form, but are straightforward to solve numerically using the `scipy` solver `fsolve`, starting from an initial guess that ensures the constraints remain satisfied.
