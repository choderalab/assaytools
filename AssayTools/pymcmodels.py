#!/usr/bin/env python

"""
A test of pymc for ITC.

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import numpy
import pymc

from math import sqrt, exp, log

#=============================================================================================
# Physical constants
#=============================================================================================

Na = 6.02214179e23 # Avogadro's number (number/mol)
kB = Na * 1.3806504e-23 / 4184.0 # Boltzmann constant (kcal/mol/K)
C0 = 1.0 # standard concentration (M)

#=============================================================================================
# Experimental parameters and data
#=============================================================================================

# ABRF-MIRG'02 dataset 10
V0 = 1.4301e-3 # volume of calorimeter sample cell (L)
V0 = V0 - 0.044e-3 # Tellinghuisen volume correction for VP-ITC (L)
DeltaV = 8.e-6 # injection volume (L)
P0_stated = 32.e-6 # protein stated concentration (M)
Ls_stated = 384.e-6 # ligand syringe stated concentration (M)
temperature = 298.15 # temperature (K)
dP0 = 0.10 * P0_stated # uncertainty in protein stated concentration (M) - 10% error
dLs = 0.10  * Ls_stated # uncertainty in ligand stated concentration (M) - 10% error
q_n = numpy.array([
    -13.343, -13.417, -13.279, -13.199, -13.118, -12.781, -12.600, -12.124, -11.633, -10.921, -10.009, -8.810, 
    -7.661, -6.272, -5.163, -4.228, -3.519, -3.055, -2.599, -2.512, -2.197, -2.096, -2.087, -1.959, -1.776, -1.879,
    -1.894, -1.813, -1.740, -1.810]) # integrated heats of injection (kcal/mol injectant)
q_n = q_n * DeltaV * Ls_stated * 1000.0 # convert injection heats to cal/injection
beta = 1.0 / (kB * temperature) # inverse temperature 1/(kcal/mol)

#=============================================================================================
# Model
#=============================================================================================

# Determine number of observations.
N = q_n.size

# Determine guesses for initial values
log_sigma_guess = log(q_n[N-4:N].std()) # cal/injection
DeltaG_guess = -8.3 # kcal/mol
DeltaH_guess = -12.0 # kcal/mol
DeltaH_0_guess = q_n[N-1] # cal/injection

# Determine min and max range for log_sigma
log_sigma_min = log_sigma_guess - 10
log_sigma_max = log_sigma_guess + 5

# Determine range for priors for thermodynamic parameters.
DeltaG_min = -40. # (kcal/mol)
DeltaG_max = +40. # (kcal/mol)
DeltaH_min = -100. # (kcal/mol)
DeltaH_max = +100. # (kcal/mol)
heat_interval = q_n.max() - q_n.min()
DeltaH_0_min = q_n.min() - heat_interval # (cal/mol)
DeltaH_0_max = q_n.max() + heat_interval # (cal/mol)

# Define priors.
#P0 = pymc.Normal('P0', mu=P0_stated, tau=1.0/dP0**2, value=P0_stated)
#Ls = pymc.Normal('Ls', mu=Ls_stated, tau=1.0/dLs**2, value=Ls_stated)

P0 = pymc.Lognormal('P0', mu=log(P0_stated), tau=1.0/log(1.0+(dP0/P0_stated)**2), value=P0_stated)
Ls = pymc.Lognormal('Ls', mu=log(Ls_stated), tau=1.0/log(1.0+(dLs/Ls_stated)**2), value=Ls_stated)

log_sigma = pymc.Uniform('log_sigma', lower=log_sigma_min, upper=log_sigma_max, value=log_sigma_guess)
DeltaG = pymc.Uniform('DeltaG', lower=DeltaG_min, upper=DeltaG_max, value=DeltaG_guess)
DeltaH = pymc.Uniform('DeltaH', lower=DeltaH_min, upper=DeltaH_max, value=DeltaH_guess)
DeltaH_0 = pymc.Uniform('DeltaH_0', lower=DeltaH_0_min, upper=DeltaH_0_max, value=DeltaH_0_guess)
        
#=============================================================================================
# Initial guess for thermodynaic parameters
#=============================================================================================

@pymc.deterministic
def expected_injection_heats(P0=P0, Ls=Ls, DeltaG=DeltaG, DeltaH=DeltaH, DeltaH_0=DeltaH_0, q_n_obs=q_n):
    """
    Expected heats of injection for two-component binding model.
    
    ARGUMENTS
    
    DeltaG - free energy of binding (kcal/mol)
    DeltaH - enthalpy of binding (kcal/mol)
    DeltaH_0 - heat of injection (cal/mol)
    
    """

    debug = False

    Kd = exp(beta * DeltaG) * C0 # dissociation constant (M)

    # Compute complex concentrations.
    Pn = numpy.zeros([N], numpy.float64) # Pn[n] is the protein concentration in sample cell after n injections (M)
    Ln = numpy.zeros([N], numpy.float64) # Ln[n] is the ligand concentration in sample cell after n injections (M)
    PLn = numpy.zeros([N], numpy.float64) # PLn[n] is the complex concentration in sample cell after n injections (M)
    for n in range(N):
        # Instantaneous injection model (perfusion)
        d = 1.0 - (DeltaV / V0) # dilution factor (dimensionless)
        P = V0 * P0 * d**(n+1) # total quantity of protein in sample cell after n injections (mol)
        L = V0 * Ls * (1. - d**(n+1)) # total quantity of ligand in sample cell after n injections (mol)
        PLn[n] = 0.5/V0 * ((P + L + Kd*V0) - sqrt((P + L + Kd*V0)**2 - 4*P*L));  # complex concentration (M)
        Pn[n] = P/V0 - PLn[n]; # free protein concentration in sample cell after n injections (M)
        Ln[n] = L/V0 - PLn[n]; # free ligand concentration in sample cell after n injections (M)
        
    # Compute expected injection heats.
    q_n = numpy.zeros([N], numpy.float64) # q_n_model[n] is the expected heat from injection n
    # Instantaneous injection model (perfusion)
    d = 1.0 - (DeltaV / V0) # dilution factor (dimensionless)
    q_n[0] = (1000.0*DeltaH) * V0 * PLn[0] + DeltaH_0 # first injection
    for n in range(1,N):
        q_n[n] = (1000.0*DeltaH) * V0 * (PLn[n] - d*PLn[n-1]) + DeltaH_0 # subsequent injections

    # Debug output
    if debug:
        print "DeltaG = %6.1f kcal/mol ; DeltaH = %6.1f kcal/mol ; DeltaH_0 = %6.1f ucal/injection" % (DeltaG, DeltaH, DeltaH_0*1e6)
        for n in range(N):
            print "%6.1f" % (PLn[n]*1e6),
        print ""
        for n in range(N):
            print "%6.1f" % (q_n[n]*1e6),
        print ""
        for n in range(N):
            print "%6.1f" % (q_n_obs[n]*1e6),
        print ""
        print ""
    
    return q_n

@pymc.deterministic
def tau(log_sigma=log_sigma):
    """
    Injection heat measurement precision.

    """
    return exp(-2.0*log_sigma)

# Define observed data.
observed_injection_heats = pymc.Normal('q_n', size=N, mu=expected_injection_heats, tau=tau, observed=True, value=q_n)
