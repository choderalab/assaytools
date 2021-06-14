#!/usr/bin/env python

"""
Various ligand binding models for use in assays.

"""

# =============================================================================================
# IMPORTS
# =============================================================================================

import numpy as np
import jax.numpy as jnp
import copy

from math import sqrt, exp, log

# =============================================================================================
# Physical constants
# =============================================================================================

Na = 6.02214179e23  # Avogadro's number (number/mol)
kB = Na * 1.3806504e-23 / 4184.0  # Boltzmann constant (kcal/mol/K)
C0 = 1.0  # standard concentration (M)

# =============================================================================================
# Binding models
# =============================================================================================


class BindingModel(object):
    """
    Abstract base class for reaction models.

    """

    def __init__(self):
        pass


# =============================================================================================
# Two-component binding model
# =============================================================================================


class TwoComponentBindingModel(BindingModel):
    """
    Simple two-component association.

    """

    @classmethod
    def equilibrium_concentrations(cls, DeltaG, Ptot, Ltot):
        """
        Compute equilibrium concentrations for simple two-component association.

        Parameters
        ----------
        DeltaG : float
           Reduced free energy of binding (in units of kT)
        Ptot : float or numpy array
           Total protein concentration summed over bound and unbound species, molarity.
        Ltot : float or numpy array
           Total ligand concentration summed over bound and unbound speciesl, molarity.

        Returns
        -------
        P : float or numpy array with same dimensions as Ptot
           Free protein concentration, molarity.
        L : float or numpy array with same dimensions as Ptot
           Free ligand concentration, molarity.
        PL : float or numpy array with same dimensions as Ptot
           Bound complex concentration, molarity.

        """
        # Handle only strictly positive elements---all others are set to zero as constants
        try:
            nonzero_indices = jnp.where(Ltot > 0)[0]
            zero_indices = jnp.where(Ltot <= 0)[0]
        except:
            nonzero_indices = jnp.array(range(Ltot.shape[0]))
            zero_indices = jnp.array([])
        nnonzero = len(nonzero_indices)
        nzeros = len(zero_indices)

        # Numerically stable variant
        dtype = jnp.float32
        Ptot = Ptot.astype(dtype)  # promote to dtype
        Ltot = Ltot.astype(dtype)  # promote to dtype
        PL = jnp.zeros(Ptot.shape, dtype)
        logP = jnp.log(jnp.take(Ptot, nonzero_indices))
        logL = jnp.log(jnp.take(Ltot, nonzero_indices))
        logPLK = jnp.logaddexp(jnp.logaddexp(logP, logL), DeltaG)
        PLK = jnp.exp(logPLK)
        sqrt_arg = 1.0 - jnp.exp(jnp.log(4.0) + logP + logL - 2.0 * logPLK)
        sqrt_arg = jnp.where(sqrt_arg >= 0.0, sqrt_arg, 0)  # ensure always positive
        PL = PL.at[nonzero_indices].set(
            0.5 * PLK * (1.0 - jnp.sqrt(sqrt_arg))
        )  # complex concentration (M)

        # Compute remaining concentrations.
        P = Ptot - PL
        # free protein concentration in sample cell after n injections (M)
        L = Ltot - PL
        # free ligand concentration in sample cell after n injections (M)

        # Ensure all concentrations are within limits, correcting cases where numerical issues cause problems.
        PL = jnp.where(PL >= 0.0, PL, 0.0)  # complex cannot have negative concentration
        P = jnp.where(P >= 0.0, P, 0.0)
        L = jnp.where(L >= 0.0, L, 0.0)

        """
        # Check all concentrations are nonnegative
        # this check doesn't work with jax as it requires concrete values (no tracer)
        # but P, L, PL all have tracers

        assert jnp.all(P >= 0)
        assert jnp.all(L >= 0)
        assert jnp.all(PL >= 0)
        """

        return [P, L, PL]
