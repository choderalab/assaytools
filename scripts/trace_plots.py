import matplotlib
matplotlib.use('Agg')

# NOTE right now this does not work for py3 because pickle files 

import cPickle
import matplotlib.pyplot as plt
import numpy as np

data_file = 'p38-Erlotinib-EF_mcmc-2017-05-05 13:16.pickle'
###ABOVE INPUT PICKLE FILE NAME NEEDS MANUAL EDITTING

with open(r'%s'%data_file,'rb') as my_file:
    data = cPickle.load(my_file)

import traceback
import matplotlib
colors = matplotlib.cm.viridis(np.linspace(0, 0.85, len(data.keys())))

fig, axs = plt.subplots(len(data.keys())/2,2, figsize=(20, 80), facecolor='w', edgecolor='k')
axs = axs.ravel()

for i, key in enumerate(data.keys()):
    try:
        axs[i].plot(data['DeltaG'][0],data[key][0],color=colors[i],marker='.',linestyle = 'None');
        axs[i].set_title('%s'%key,fontsize=16);
        plt.tight_layout()
    except Exception as e:
        # There is often an error here for _state_, which it can't plot.
        # It is fine, just uncomment the line below if you want to see it.
        #print(traceback.print_exc())
        print('%s has no [0]'%key)

    plt.savefig('trace_p38-Erlotinib-EF_May5.png')
    ###ABOVE OUTPUT PNG NAME NEEDS MANUAL EDITTING

plt.clf()
fig, axs = plt.subplots(len(data.keys())/2,2, figsize=(20, 80), facecolor='w', edgecolor='k')
axs = axs.ravel()
for i, key in enumerate(data.keys()):
    try:
        axs[i].plot(data[key][0],color=colors[i],marker='.',linestyle = 'None');
        axs[i].set_title('%s'%key,fontsize=16);
        plt.tight_layout()
    except Exception as e:
        # There is often an error here for _state_, which it can't plot.
        # It is fine, just uncomment the line below if you want to see it.
        #print(traceback.print_exc())
        print('%s has no [0]'%key)

    plt.savefig('DeltaGvall_p38-Erlotinib-EF_May5.png')
    ###ABOVE OUTPUT PNG NAME NEEDS MANUAL EDITTING

