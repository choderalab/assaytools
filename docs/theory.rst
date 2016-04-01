.. _theory:

******
Theory
******

Bayesian sampling
=================

`AssayTools` uses a Bayesian model to infer unknown parameters (such as the free energy of ligand binding) from experimental data.
It does this to allow the complete uncertainty---in the form of the joint distribution of all unknown parameters---to be characterized.
The most common way to summarize these results is generally confidence intervals of individual parameters, but much more sophisticated analyses are also possible.

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
