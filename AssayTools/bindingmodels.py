#!/usr/bin/env python

"""
Various ligand binding models for use in assays.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy as np
import copy
import pint

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
# Competitive binding model
#=============================================================================================

class CompetitiveBindingModel(BindingModel):
   """
   Competitive binding model for one protein and arbitrary number of competitive ligands.

   """

   @classmethod
   def equilibrium_concentrations(cls, Ka_n, C0_R, C0_Ln, V, c0=None):
      """
      Compute the equilibrium concentrations of each complex species for N ligands competitively binding to a receptor.

      ARGUMENTS

      Ka_n (numpy N-array of float) - Ka_n[n] is the association constant for receptor and ligand species n (1/M)
      x_R (float) - the total number of moles of receptor in the sample volume
      x_n (numpy N-array of float) - x_n[n] is the total number of moles of ligand species n in the sample volume
      V (float) - the total sample volume (L)

      RETURNS

      C_n (numpy N-array of float) - C_n[n] is the concentration of complex of receptor with ligand species n

      EXAMPLES

      >>> V = 1.4303e-3 # volume (L)
      >>> x_R = V * 510.e-3 # receptor
      >>> x_Ln = np.array([V * 8.6e-6, 200.e-6 * 55.e-6]) # ligands
      >>> Ka_n = np.array([1./(400.e-9), 1./(2.e-11)]) # association constants
      >>> from bindingmodels import CompetitiveBindingModel
      >>> C_PLn = CompetitiveBindingModel.equilibrium_concentrations(Ka_n, x_R, x_Ln, V)

      NOTES

      Each complex concentration C_n must obey the relation

      Ka_n[n] = C_RLn[n] / (C_R * C_Ln[n])           for n = 1..N

      with conservation of mass constraints

      V * (C_Ln[n] + C_RLn[n]) = x_Ln[n]             for n = 1..N



      and

      V * (C_R + C_RLn[:].sum()) = x_R

      along with the constraints

      0 <= V * C_RLn[n] <= min(x_Ln[n], x_R)         for n = 1..N
      V * C_RLn[:].sum() <= x_R

      We can rearrange these expressions to give

      V * C_R * C_Ln[n] * Ka_n[n] - V * C_RLn[n] = 0

      and eliminate C_Ln[n] and C_R to give

      V * (x_R/V - C_RLn[:].sum()) * (x_Ln[n]/V - C_RLn[n]) * Ka_n[n] - V * C_RLn[n] = 0    for n = 1..N

      """

      x_R = C0_R * V
      x_Ln = C0_Ln * V

      nspecies = Ka_n.size
      #print "x_R = ", x_R
      #print "x_Ln = ", x_Ln
      #print "x_Ln / V = ", x_Ln / V
      #print "Ka_n = ", Ka_n

      # Define optimization functions
      def func(C_RLn):
         f_n = V * (x_R/V - C_RLn[:].sum()) * (x_Ln[:]/V - C_RLn[:]) * Ka_n[:] - V * C_RLn[:]
         #print "f_n = ", f_n
         return f_n

      def fprime(C_RLn):
         nspecies = C_RLn.size
         G_nm = np.zeros([nspecies,nspecies], np.float64) # G_nm[n,m] is the derivative of func[n] with respect to C_RLn[m]
         for n in range(nspecies):
            G_nm[n,:] = - V * (x_Ln[:]/V - C_RLn[:]) * Ka_n[:]
            G_nm[n,n] -= V * (Ka_n[n] * (x_R/V - C_RLn[:].sum()) + 1.0)
         return G_nm

      def sfunc(s):
         #print "s = ", s
         f_n = V * (x_R/V - (s[:]**2).sum()) * (x_Ln[:]/V - s[:]**2) * Ka_n[:] - V * s[:]**2
         #print "f_n = ", f_n
         return f_n

      def sfprime(s):
         nspecies = s.size
         G_nm = np.zeros([nspecies,nspecies], np.float64) # G_nm[n,m] is the derivative of func[n] with respect to C_RLn[m]
         for n in range(nspecies):
            G_nm[n,:] = - V * (x_Ln[:]/V - s[:]**2) * Ka_n[:]
            G_nm[n,n] -= V * (Ka_n[n] * (x_R/V - (s[:]**2).sum()) + 1.0)
            G_nm[n,:] *= 2. * s[n]
         return G_nm

      # Allocate storage for complexes
      # Compute equilibrium concentrations.
      #x0 = np.zeros([nspecies], np.float64)
      #x0 = (x_Ln / V).copy()
      #x = scipy.optimize.fsolve(func, x0, fprime=fprime)
      #C_RLn = x

      #x0 = np.sqrt(x_Ln / V).copy()
      #x = scipy.optimize.fsolve(sfunc, x0, fprime=sfprime)
      #C_RLn = x**2

      def objective(x):
         f_n = func(x)
         G_nm = fprime(x)

         obj = (f_n**2).sum()
         grad = 0.0 * f_n
         for n in range(f_n.size):
            grad += 2 * f_n[n] * G_nm[n,:]

         return (obj, grad)

      #x0 = np.zeros([nspecies], np.float64)
      #bounds = list()
      #for n in range(nspecies):
      #   m = min(C0_R, C0_Ln[n])
      #   bounds.append( (0., m) )
      #[x, a, b] = scipy.optimize.fmin_l_bfgs_b(objective, x0, bounds=bounds)
      #C_RLn = x

      def ode(c_n, t, Ka_n, x_Ln, x_R):
         dc_n = - c_n[:] + Ka_n[:] * (x_Ln[:]/V - c_n[:]) * (x_R/V - c_n[:].sum())
         return dc_n

      def odegrad(c_n, t, Ka_n, x_Ln, x_R):
         N = c_n.size
         d2c = np.zeros([N,N], np.float64)
         for n in range(N):
            d2c[n,:] = -Ka_n[n] * (x_Ln[n]/V - c_n[n])
            d2c[n,n] += -(Ka_n[n] * (x_R/V - c_n[:].sum()) + 1.0)
         return d2c

      #if c0 is None: c0 = np.zeros([nspecies], np.float64)
      #maxtime = 100.0 * (x_R/V) / Ka_n.max()
      #time = [0, maxtime / 2.0, maxtime]
      #c = scipy.integrate.odeint(ode, c0, time, Dfun=odegrad, args=(Ka_n, x_Ln, x_R))
      #C_RLn = c[-1,:]

      #c = np.zeros([nspecies], np.float64)
      #maxtime = 1.0 / Ka_n.min()
      #maxtime = 1.0 / ((x_R/V) * Ka_n.min())
      #maxtime = 1.0
      #time = [0, maxtime]
      #c = scipy.optimize.fsolve(ode, c, fprime=odegrad, args=(0.0, Ka_n, x_Ln, x_R), xtol=1.0e-6)
      #c = scipy.integrate.odeint(ode, c, time, Dfun=odegrad, args=(Ka_n, x_Ln, x_R), mxstep=50000)
      #c = c[-1,:]
      #C_RLn = c

      #print "C_RLn = ", C_RLn
      #print ""

      c = np.zeros([nspecies], np.float64)
      sorted_indices = np.argsort(-x_Ln)
      for n in range(nspecies):
         indices = sorted_indices[0:n+1]
         #c[indices] = scipy.optimize.fsolve(ode, c[indices], fprime=odegrad, args=(0.0, Ka_n[indices], x_Ln[indices], x_R), xtol=1.0e-6, warning=False)
         c[indices] = scipy.optimize.fsolve(ode, c[indices], fprime=odegrad, args=(0.0, Ka_n[indices], x_Ln[indices], x_R), xtol=1.0e-6)
      C_RLn = c

      return C_RLn


#=============================================================================================
# General robust competitive binding model
#=============================================================================================

class GeneralBindingModel(BindingModel):
   """
   General robust binding model for one protein and arbitrary number of competitive ligands.

   """

   @classmethod
   def equilibrium_concentrations(cls, reactions, conservation_equations):
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
      sol = root(ftarget, X, jac=True)
      if (sol.success == False):
          msg  = "root-finder failed to converge:\n"
          msg += str(sol)
          raise Exception(msg)

      log_concentrations = { all_species[index] : sol.x[index] for index in range(nspecies) }
      return log_concentrations
