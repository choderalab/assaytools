import json
import numpy as np

inputs = {
    'xml_file_path' :  "./data/single_wavelength_copy",
    'section'       :  '280_480_TOP_120',
    'ligand_order'  :  ['Bosutinib','Bosutinib Isomer','Erlotinib','Gefitinib','Ponatinib','Lapatinib','Saracatinib','Vandetanib'],
    'Lstated'       :  np.array([20.0e-6,14.0e-6,9.82e-6,6.88e-6,4.82e-6,3.38e-6,2.37e-6,1.66e-6,1.16e-6,0.815e-6,0.571e-6,0.4e-6,0.28e-6,0.196e-6,0.138e-6,0.0964e-6,0.0676e-6,0.0474e-6,0.0320e-6,0.0240e-6,0.0160e-6,0.0120e-6,0.008e-6,0.0], np.float64), # ligand concentration, M
    'protein'       :  'p38',
    'Pstated'       :  0.5e-6 * np.ones([24],np.float64), # protein concentration, M
    'assay_volume'  :  50e-6, # assay volume, L
    'well_area'     :  0.1369, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)
