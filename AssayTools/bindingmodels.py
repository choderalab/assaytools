#!/usr/bin/env python

"""
Various ligand binding models for use in assays.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy as np
import copy

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

      # Original form:
      #Kd = np.exp(DeltaG)
      #sqrt_arg = (Ptot + Ltot + Kd)**2 - 4*Ptot*Ltot
      #sqrt_arg[sqrt_arg < 0.0] = 0.0
      #PL = 0.5 * ((Ptot + Ltot + Kd) - np.sqrt(sqrt_arg));  # complex concentration (M)

      # Numerically stable variant?
      logP = np.log(Ptot)
      logL = np.log(Ltot)
      logPLK = np.logaddexp(np.logaddexp(logP, logL), DeltaG)
      PLK = np.exp(logPLK);
      sqrt_arg = 1.0 - np.exp(np.log(4.0) + logP + logL - 2*logPLK);
      sqrt_arg[sqrt_arg < 0.0] = 0.0 # ensure always positive
      PL = 0.5 * PLK * (1.0 - np.sqrt(sqrt_arg));  # complex concentration (M)

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

   @classmethod
   def equilibrium_concentrations(cls, reactions, conservation_equations, tol=1.0e-10):
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

      # Solve
      from scipy.optimize import root
      X = np.zeros([nspecies], np.float64)
      sol = root(ftarget, X, jac=True, tol=tol)
      if (sol.success == False):
          msg  = "root-finder failed to converge:\n"
          msg += str(sol)
          raise Exception(msg)

      log_concentrations = { all_species[index] : sol.x[index] for index in range(nspecies) }
      return log_concentrations
