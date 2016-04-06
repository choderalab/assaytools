"""
Test binding models
"""


import numpy as np
from scipy.misc import logsumexp
import copy

def test_general_binding_model():
    from assaytools.bindingmodels import GeneralBindingModel

    decimal = 7 # number of decimal places to expect agreement

    # Simple 1:1 binding
    log_affinity_true = -10 # ln (Kd / M)
    log_Rtot_true = -6 # ln ([R]tot / M)
    log_Ltot_true = -6 # ln ([L]tot / M)

    reactions = [ (log_affinity_true, {'RL': -1, 'R' : +1, 'L' : +1}) ]
    conservation_equations = [ (log_Rtot_true, {'RL' : +1, 'R' : +1}), (log_Ltot_true, {'RL' : +1, 'L' : +1}) ]
    log_concentrations = GeneralBindingModel.equilibrium_concentrations(reactions, conservation_equations)

    log_affinity = (log_concentrations['R'] + log_concentrations['L'] - log_concentrations['RL'])
    np.testing.assert_almost_equal(log_affinity, -10, decimal=decimal, err_msg="log_affinity was %f, expected %f." % (log_affinity, log_affinity_true))
    log_Rtot = logsumexp([log_concentrations['R'], log_concentrations['RL']])
    np.testing.assert_almost_equal(log_Rtot, log_Rtot_true, decimal=decimal, err_msg="log [R]tot was %f, expected %f" % (log_Rtot, log_Rtot_true))
    log_Ltot = logsumexp([log_concentrations['L'], log_concentrations['RL']])
    np.testing.assert_almost_equal(log_Ltot, log_Ltot_true, decimal=decimal, err_msg="log [L]tot was %f, expected %f" % (log_Ltot, log_Ltot_true))

if __name__ == "__main__":
   test_general_binding_model()
