#!/usr/bin/env python

"""
Various ligand binding models for use in assays.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy as np
import copy
import time

import scipy.optimize
import scipy.integrate

from math import sqrt, exp, log

#=============================================================================================
# Physical constants
#=============================================================================================

Na = 6.02214179e23 # Avogadro's number (number/mol)
kB = Na * 1.3806504e-23 / 4184.0 # Boltzmann constant (kcal/mol/K)
C0 = 1.0 # standard concentration (M)

#=============================================================================================
# Binding models
#=============================================================================================

class BindingModel(object):
    """
    Abstract base class for reaction models.

    """

    def __init__(self):
        pass


#=============================================================================================
# LRU cache that supports mutable args
#=============================================================================================

import typing
from collections import OrderedDict
from functools import lru_cache, wraps

def lru_cache2(size=128):

    cache = OrderedDict()

    def make_key(x):
        if isinstance(x, typing.Mapping):
            return (id(type(x)), *((k, make_key(v)) for k, v in sorted(x.items())))
        elif isinstance(x, typing.Set):
            return (id(type(x)), *map(make_key, sorted(x)))
        elif isinstance(x, typing.Sequence):
            return (id(type(x)), *map(make_key, x))
        else:
            try:
                hash(x)
            except TypeError:
                return id(type(x)), str(x)
            else:
                return id(type(x)), x

    def decorator(fn):
        @wraps(fn)
        def wrapped(*args, **kwargs):
            key = make_key((id(fn), args, kwargs))
            try:
                return cache[key]
            except KeyError:
                res = fn(*args, **kwargs)
                cache[key] = res
                while len(cache) >= size:
                    cache.popitem(last=False)
                return res
        return wrapped
    return decorator

#=============================================================================================
# Two-component binding model
#=============================================================================================

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
      initial_time = time.time()

      # Handle only strictly positive elements---all others are set to zero as constants
      try:
          nonzero_indices = np.where(Ltot > 0)[0]
          zero_indices = np.where(Ltot <= 0)[0]
      except:
          nonzero_indices = range(size[0])
          zero_indices = []
      nnonzero = len(nonzero_indices)
      nzeros = len(zero_indices)

      # Original form:
      #Kd = np.exp(DeltaG)
      #sqrt_arg = (Ptot + Ltot + Kd)**2 - 4*Ptot*Ltot
      #sqrt_arg[sqrt_arg < 0.0] = 0.0
      #PL = 0.5 * ((Ptot + Ltot + Kd) - np.sqrt(sqrt_arg));  # complex concentration (M)


      # Numerically stable variant?
      PL = np.zeros(Ptot.shape)
      logP = np.log(Ptot[nonzero_indices])
      logL = np.log(Ltot[nonzero_indices])
      logPLK = np.logaddexp(np.logaddexp(logP, logL), DeltaG)
      PLK = np.exp(logPLK);
      sqrt_arg = 1.0 - np.exp(np.log(4.0) + logP + logL - 2*logPLK);
      sqrt_arg[sqrt_arg < 0.0] = 0.0 # ensure always positive
      PL[nonzero_indices] = 0.5 * PLK * (1.0 - np.sqrt(sqrt_arg));  # complex concentration (M)

      # Another variant
      #PL = 2*Ptot*Ltot / ((Ptot+Ltot+Kd) + np.sqrt((Ptot + Ltot + Kd)**2 - 4*Ptot*Ltot));  # complex concentration (M)
      # Yet another numerically stable variant?
      #logPLK = np.logaddexp(np.log(Ptot + Ltot),  DeltaG);
      #PLK = np.exp(logPLK);
      #xy = np.exp(np.log(Ptot) + np.log(Ltot) - 2.0*logPLK);
      #chi = 1.0 - 4.0 * xy;
      #chi[chi < 0.0] = 0.0 # prevent square roots of negative numbers
      #PL = 0.5 * PLK * (1 - np.sqrt(chi))

      # Ensure all concentrations are within limits, correcting cases where numerical issues cause problems.
      PL[PL < 0.0] = 0.0 # complex cannot have negative concentration
      #PL_max = np.minimum(Ptot, Ltot)
      #indices = np.where(PL > PL_max)
      #PL[indices] = PL_max[indices]

      # Compute remaining concentrations.
      P = Ptot - PL; # free protein concentration in sample cell after n injections (M)
      L = Ltot - PL; # free ligand concentration in sample cell after n injections (M)

      return [P, L, PL]

#=============================================================================================
# General robust competitive binding model
#=============================================================================================

class GeneralBindingModel(BindingModel):
   """
   General robust binding model for one protein and arbitrary number of competitive ligands.

   """

   _time = 0
   _ncalls = 0

   @classmethod
   #@lru_cache2(size=128)
   def equilibrium_concentrations(cls, reactions, conservation_equations, tol=1.0e-8):
      """
      Compute the equilibrium concentrations of each complex species for a general set of binding reactions.

      Parameters
      ----------
      reactions : list
          List of binding reactions.
          Each binding reaction is encoded as a tuple of (log equilibrium constant, dict of stoichiometry)
          Example: K_d = [RL] / ([R] [L]) becomes [ (-10, {'RL': -1, 'R' : +1, 'L' : +1}) ]
      conservation_equations : list
          List of mass conservation laws.
          Each mass conservation law is encoded as a tuple of (log total concentration, dict of stoichiometry of all species)
          Example: [R]tot = 10^-6 M = [RL] + [R] and [L]tot = 10^-6 M = [RL] + [L] becomes [ (-6, {'RL' : +1, 'R' : +1}), (-6, {'RL' : +1, 'L' : +1}) ]
      tol : float, optional, default=1.0e-8
          Solution tolerance.

      Returns
      -------
      log_concentrations : dict of str : float
          log_concentrations[species] is the log concentration of specified species

      Examples
      --------

      Simple 1:1 association

      >>> reactions = [ (-10, {'RL': -1, 'R' : +1, 'L' : +1}) ]
      >>> conservation_equations = [ (-6, {'RL' : +1, 'R' : +1}), (-6, {'RL' : +1, 'L' : +1}) ]
      >>> log_concentrations = GeneralBindingModel.equilibrium_concentrations(reactions, conservation_equations)

     Competitive 1:1 association

      >>> reactions = [ (-10, {'RL': -1, 'R' : +1, 'L' : +1}), (-5, {'RP' : -1, 'R' : +1, 'P' : +1}) ]
      >>> conservation_equations = [ (-6, {'RL' : +1, 'RP' : +1, 'R' : +1}), (-6, {'RL' : +1, 'L' : +1}), (-5, {'RP' : +1, 'P' : +1}) ]
      >>> log_concentrations = GeneralBindingModel.equilibrium_concentrations(reactions, conservation_equations)

      TODO
      ----
      * Can we allow the caller to specify initial conditions instead of conservation laws?

      """

      nreactions = len(reactions)
      nconservation = len(conservation_equations)
      nequations = nreactions + nconservation

      # Determine names of all species.
      all_species = set()
      for (log_equilibrium_constant, reaction) in reactions:
          for species in reaction.keys():
              all_species.add(species)
      all_species = list(all_species) # order is now fixed
      nspecies = len(all_species)

      # Construct function with appropriate roots.
      # TODO: Meta-program efficient form of target function for speed?
      # TODO: Form the logK, logQ, S, and C  matrices in numpy outside of this function so that operations can be numpy-accelerated
      # http://assaytools.readthedocs.io/en/latest/theory.html#general-binding-model
      # Rewrite theory docs in these matrices as well.

      # Form equations as numpy
      logK = np.zeros([nreactions], np.float64)
      S = np.zeros([nreactions, nspecies], np.float64)
      C = np.zeros([nconservation, nspecies], np.float64)
      logQ = np.zeros([nconservation], np.float64)
      equation_index = 0
      # Reactions
      for (reaction_index, (log_equilibrium_constant, reaction)) in enumerate(reactions):
          logK[reaction_index] = log_equilibrium_constant
          for (species_index, species) in enumerate(all_species):
              if species in reaction:
                  S[reaction_index, species_index] = reaction[species]
      # Conservation of mass
      for (equation_index, (log_total_concentration, conservation_equation)) in enumerate(conservation_equations):
          logQ[equation_index] = log_total_concentration
          for (species_index, species) in enumerate(all_species):
              if species in conservation_equation:
                  C[equation_index, species_index] = conservation_equation[species]

      # Define target function
      from scipy.misc import logsumexp
      def ftarget_numpy(logX):
          target = np.zeros([nequations], np.float64)
          jacobian = np.zeros([nequations, nspecies], np.float64)

          # Reactions
          target[0:nreactions] = - logK[:] + np.dot(S[:,:], logX[:])
          jacobian[0:nreactions,:] = S[:,:]
          # Conservation
          target[nreactions:] = - logQ[:] + logsumexp(C + logX, axis=1)
          for equation_index in range(nconservation):
              nonzero_indices = np.where(C[equation_index,:] != 0)[0]
              logsum = logsumexp(np.log(C[equation_index, nonzero_indices]) + logX[nonzero_indices])
              jacobian[nreactions+equation_index,:] = C[equation_index,:] * np.exp(logX[:] - logsum)

          return (target, jacobian)

      # Construct function with appropriate roots.
      def ftarget(X):
          target = np.zeros([nequations], np.float64)
          jacobian = np.zeros([nequations, nspecies], np.float64)
          equation_index = 0
          # Reactions
          for (log_equilibrium_constant, reaction) in reactions:
              target[equation_index] = - log_equilibrium_constant
              for (species_index, species) in enumerate(all_species):
                  if species in reaction:
                      stoichiometry = reaction[species]
                      target[equation_index] += stoichiometry * X[species_index]
                      jacobian[equation_index][species_index] = stoichiometry
              equation_index += 1
          # Conservation of mass
          from scipy.misc import logsumexp
          for (log_total_concentration, conservation_equation) in conservation_equations:
              target[equation_index] = - log_total_concentration
              log_concentrations = list()
              for (species_index, species) in enumerate(all_species):
                  if species in conservation_equation:
                      stoichiometry = conservation_equation[species]
                      log_concentrations.append(X[species_index] + np.log(stoichiometry))
              log_concentrations = np.array(log_concentrations)
              logsum = logsumexp(log_concentrations)
              target[equation_index] += logsum
              for (species_index, species) in enumerate(all_species):
                  log_concentrations = list()
                  if species in conservation_equation:
                      stoichiometry = conservation_equation[species]
                      jacobian[equation_index, species_index] = stoichiometry * np.exp(X[species_index] - logsum)
              equation_index += 1

          return (target, jacobian)

      # Construct initial guess
      # We assume that all matter is equally spread out among all species via the conservation equations
      from scipy.misc import logsumexp
      LOG_ZERO = -100
      X = LOG_ZERO * np.ones([nspecies], np.float64)
      for (species_index, species) in enumerate(all_species):
          for (log_total_concentration, conservation_equation) in conservation_equations:
              log_total_stoichiometry = np.log(np.sum([stoichiometry for stoichiometry in conservation_equation.values()]))
              if species in conservation_equation:
                  stoichiometry = conservation_equation[species]
                  X[species_index] = logsumexp([X[species_index], log_total_concentration + np.log(stoichiometry) - log_total_stoichiometry])

      # Solve
      from scipy.optimize import root
      options = {'xtol' : tol}
      initial_time = time.time()
      sol = root(ftarget, X, method='lm', jac=True, tol=tol, options=options)
      final_time = time.time()
      cls._time += (final_time - initial_time)
      cls._ncalls += 1
      if (cls._ncalls % 1000 == 0):
          print('%8d calls : %8.3f s total (%8.3f ms/call)' % (cls._ncalls, cls._time, cls._time/cls._ncalls*1000))

      if (sol.success == False):
          msg  = "root-finder failed to converge:\n"
          msg += str(sol)
          raise Exception(msg)

      log_concentrations = { all_species[index] : sol.x[index] for index in range(nspecies) }
      return log_concentrations
