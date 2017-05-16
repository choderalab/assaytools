import numpy as np
import tensorflow as tf
from scipy.misc import logsumexp
import copy

from assaytools.bindingmodels import GeneralBindingModel

decimal = 7 # number of decimal places to expect agreement
nsamples = 100


def equilibrium_concentrations(logK, S, logQ, C, tol=1.0e-8):
   """
   Compute the equilibrium concentrations of each complex species for a general set of binding reactions.

   Parameters
   ----------
   logK : numpy.array with shape [nreactions]
       Corresponding log equilibrium constants for each reaction
   S : numpy.array with shape [nreactions,nspecies]
       reaction_stoichiometry[reaction,species] the stoichiometry of each species
   logQ : numpy.array with shape [nconservation]
       log_concentrations[conservation] is the log concentration of overall species 'conservation'
   C : numpy.array with shape [nconservation,nspecies]
       conservation_equations[conservation,species] is the number of 'species' involved in overall species 'conservation'
   tol : float, optional, default=1.0e-8
       Solution tolerance.

   Returns
   -------
   log_species_concentrations : numpy.array with shape [nspecies]
       log_concentrations[species] is the log concentration of specified species

   Examples
   --------

   Simple 1:1 association.

   Species are RL, R, L
   Conservation equations are for R, L

   >>> log_equilibrium_constants = np.array([-10])
   >>> reaction_stoichiometry = np.array([[-1, +1, +1]])
   >>> log_concentrations = np.array([-6, -6])
   >>> conservation_equations = np.array([[+1, +1, 0], [+1, 0, +1]])
   >>> log_species_concentrations = equilibrium_concentrations(log_equilibrium_constants, reaction_stoichiometry, log_concentrations, conservation_equations)

   Competitive 1:1 association

   Species are RL, RP, R, P, L
   Conservation equations are for R, P, L

   >>> log_equilibrium_constants = np.array([-10, -5])
   >>> reaction_stoichiometry = np.array([[-1, 0, +1, 0, +1], [0, -1, +1, +1, 0]])
   >>> log_concentrations = np.array([-6, -5, -6])
   >>> conservation_equations = np.array([[+1, +1, +1, 0, 0], [0, +1, 0, +1, 0], [+1, 0, 0, 0, +1]])
   >>> log_species_concentrations = equilibrium_concentrations(log_equilibrium_constants, reaction_stoichiometry, log_concentrations, conservation_equations)

   TODO
   ----
   * Can we allow the caller to specify initial conditions instead of conservation laws?

   """

   # Get dimensions
   [nreactions, nspecies] = S.shape
   [nconservation, nspecies] = C.shape
   nequations = nreactions + nconservation

   # Construct function with appropriate roots.
   from scipy.misc import logsumexp
   def ftarget(logX):
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

       return [target, jacobian]

   # Construct initial guess
   # We assume that all matter is equally spread out among all species via the conservation equations
   from scipy.misc import logsumexp
   LOG_ZERO = -100
   X = LOG_ZERO * np.ones([nspecies], np.float64)
   for species_index in range(nspecies):
       for conservation_equation in range(nconservation):
           log_total_concentration = logQ[conservation_equation]
           log_total_stoichiometry = np.log(C[conservation_equation,:].sum())
           stoichiometry = C[conservation_equation,species_index]
           if (stoichiometry > 0):
               X[species_index] = logsumexp([X[species_index], log_total_concentration + np.log(stoichiometry) - log_total_stoichiometry])

   # Solve
   from scipy.optimize import root
   options = {'xtol' : tol}
   sol = root(ftarget, X, method='lm', jac=True, tol=tol, options=options)
   if (sol.success == False):
       msg  = "root-finder failed to converge:\n"
       msg += str(sol)
       raise Exception(msg)

   return X

#
# Simple 1:1 binding
#


# RL R L
logK = np.array([-10], np.float64)
S = np.array([[-1, +1, +1]], np.float64)
# R L
logQ = np.array([-6, -6], np.float64)
C = np.array([[+1, +1, 0], [+1, 0, +1]], np.float64)

# Get dimensions
[nreactions, nspecies] = S.shape
[nconservation, nspecies] = C.shape
nequations = nreactions + nconservation

#
# TensorFlow model
#

# Input concentrations
#tf_logK = tf.Variable(logK, tf.float64)
#tf_logQ = tf.Variable(logQ, tf.float64)
tf_S = tf.constant(S, tf.float64)
tf_C = tf.constant(C, tf.float64)

# Input and output
tf_logK = tf.placeholder(tf.float64, shape=nreactions, name='tf_logK')
tf_logQ = tf.placeholder(tf.float64, shape=nconservation, name='tf_logQ')

# Model
tf_eq = tf.py_func(equilibrium_concentrations, [tf_logK, tf_S, tf_logQ, tf_C], tf.float64)

#
# Initialize
#

# training loop
init = tf.global_variables_initializer()
sess = tf.Session()
sess.run(init) # reset values

decimal = 7 # number of decimal places to expect agreement
nsamples = 10

#
# Simple 1:1 binding
#

for sample in range(nsamples):
    log_affinity_true = np.random.uniform(-15, 0)
    log_Rtot_true = np.random.uniform(-15, 0)
    log_Ltot_true = np.random.uniform(-15, 0)

    logK = np.array([log_affinity_true])
    logQ = np.array([log_Rtot_true, log_Ltot_true])

    log_concentrations = sess.run(tf_eq, {tf_logK : logK, tf_logQ : logQ})

    print(log_concentrations)

    log_affinity = (log_concentrations[1] + log_concentrations[2] - log_concentrations[0])
    np.testing.assert_almost_equal(log_affinity, log_affinity_true, decimal=decimal, err_msg="log_affinity was %f, expected %f." % (log_affinity, log_affinity_true))
    log_Rtot = logsumexp([log_concentrations[1], log_concentrations[0]])
    np.testing.assert_almost_equal(log_Rtot, log_Rtot_true, decimal=decimal, err_msg="log [R]tot was %f, expected %f" % (log_Rtot, log_Rtot_true))
    log_Ltot = logsumexp([log_concentrations[2], log_concentrations[0]])
    np.testing.assert_almost_equal(log_Ltot, log_Ltot_true, decimal=decimal, err_msg="log [L]tot was %f, expected %f" % (log_Ltot, log_Ltot_true))
