# This is a script version of DMSO-backfill.ipynb

from assaytools import platereader
import matplotlib.pyplot as plt
import string
import numpy as np

part1 = "p38_singlet1_20160420_153238.xml"
part2 = "p38_singlet2_20160420_154750.xml"

part1_data = platereader.read_icontrol_xml(part1)
part2_data = platereader.read_icontrol_xml(part2)

part1_data_C = platereader.select_data(part1_data, '280_480_TOP_120', 'C')
part1_data_D = platereader.select_data(part1_data, '280_480_TOP_120', 'D')
part1_data_K = platereader.select_data(part1_data, '280_480_TOP_120', 'K')
part1_data_L = platereader.select_data(part1_data, '280_480_TOP_120', 'L')

well = dict()
for j in string.ascii_uppercase:
    for i in range(1,25):
        well['%s' %j + '%s' %i] = i

Lstated = np.array([20.0e-6,14.0e-6,9.82e-6,6.88e-6,4.82e-6,3.38e-6,2.37e-6,1.66e-6,1.16e-6,0.815e-6,0.571e-6,0.4e-6,0.28e-6,0.196e-6,0.138e-6,0.0964e-6,0.0676e-6,0.0474e-6,0.0320e-6,0.0240e-6,0.0160e-6,0.0120e-6,0.008e-6,0.00001e-6], np.float64) # ligand concentration, M

# Stated concentrations of protein and ligand.
Pstated = 0.5e-6 * np.ones([24],np.float64) # protein concentration, M

# Assay configuration details
import math
assay_volume = 50e-6 # assay volume, L
well_area = 0.1369 # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
path_length = assay_volume * 1000 / well_area # cm, needed for inner filter effect corrections     

# Uncertainties in protein and ligand concentrations.
dPstated = 0.35 * Pstated # protein concentration uncertainty
dLstated = 0.08 * Lstated # ligand concentraiton uncertainty (due to gravimetric preparation and HP D300 dispensing)

def reorder2list(data,wells):
    
    sorted_keys = sorted(well.keys(), key=lambda k:well[k])
    
    reorder_data = []
    
    for key in sorted_keys:
        try:
            reorder_data.append(data[key])
        except:
            pass
        
    reorder_data = np.asarray(reorder_data,np.float64)
    
    return reorder_data

def quick_model(protein_data, buffer_data,name):
    reorder_protein = reorder2list(protein_data,well)
    reorder_buffer = reorder2list(buffer_data,well)
    
    print name   
 
    from assaytools import pymcmodels
    pymc_model = pymcmodels.make_model(Pstated, dPstated, Lstated, dLstated, 
               top_complex_fluorescence=reorder_protein,
               top_ligand_fluorescence=reorder_buffer,
               use_primary_inner_filter_correction=True, 
               use_secondary_inner_filter_correction=True, 
               assay_volume=assay_volume, DG_prior='uniform')
    
    mcmc = pymcmodels.run_mcmc(pymc_model)
    
    from assaytools import plots
    figure = plots.plot_measurements(Lstated, Pstated, pymc_model, mcmc=mcmc)
    
    map = pymcmodels.map_fit(pymc_model)
    
    pymcmodels.show_summary(pymc_model, map, mcmc)
    
    DeltaG = map.DeltaG.value
    
    np.save('DeltaG_%s.npy'%name,DeltaG)
    np.save('DeltaG_trace_%s.npy'%name,mcmc.DeltaG.trace())
    

quick_model(part1_data_C,part1_data_D,'BSI_w_backfill')

quick_model(part1_data_K,part1_data_L,'BSI')

DeltaG_trace_BSI_w_backfill = np.load('DeltaG_trace_BSI_w_backfill.npy')
DeltaG_trace_BSI = np.load('DeltaG_trace_BSI.npy')

plt.clf()

plt.hist(DeltaG_trace_BSI, 40, alpha=0.75, label='Bsi no backfill')
plt.hist(DeltaG_trace_BSI_w_backfill, 40, alpha=0.75, label='Bsi with backfill')
plt.plot([DeltaG_trace_BSI_w_backfill.mean(),DeltaG_trace_BSI_w_backfill.mean()],[0, 350],'g-')
plt.plot([DeltaG_trace_BSI.mean(),DeltaG_trace_BSI.mean()],[0, 350],'b-')
#plt.plot([literature_Bos,literature_Bos],[0, 350],'r-',label='IUPHARM data')
plt.xlabel('$\Delta G$ ($k_B T$)');
plt.ylabel('$P(\Delta G)$');
plt.title('histogram of estimates for binding free energy');
plt.legend(loc=0);

plt.savefig('BSI_w_backfill.png')
