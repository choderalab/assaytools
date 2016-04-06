# Simple 1:1 binding
reactions = [ (-10, {'RL': -1, 'R' : +1, 'L' : +1}) ]
conservation_equations = [ (-6, {'RL' : +1, 'R' : +1}), (-6, {'RL' : +1, 'L' : +1}) ]
from assaytools.bindingmodels import GeneralBindingModel
log_concentrations = GeneralBindingModel.equilibrium_concentrations(reactions, conservation_equations)
print(log_concentrations)

# check solution
import numpy as np
print("log(K_d) = -10")
print(log_concentrations['R'] + log_concentrations['L'] - log_concentrations['RL'])
print("[R]tot = %f" % np.exp(-6))
print(np.exp(log_concentrations['R']) + np.exp(log_concentrations['RL']))
print("[L]tot = %f" % np.exp(-6))
print(np.exp(log_concentrations['L']) + np.exp(log_concentrations['RL']))
