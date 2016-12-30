import json
import numpy as np

inputs = {
    'xml_file_path' :  "./data/Abl_Gef_Ima_copy/",
    'section'       :  '280_BottomRead',
    'ligand_order'  :  ['Gef-Ima','Gef-Ima','Gef-Ima','Gef-Ima'],
    'Lstated'       :  np.array([20.0e-6,9.15e-6,4.18e-6,1.91e-6,0.875e-6,0.4e-6,0.183e-6,0.0837e-6,0.0383e-6,0.0175e-6,0.008e-6,0.0001e-6], np.float64), # ligand concentration
    'protein'       :  'Abl',
    'Pstated'       :  0.5e-6 * np.ones([12],np.float64), # protein concentration, M
    'assay_volume'  :  100e-6, # assay volume, L
    'well_area'     :  0.3969, # well area, cm^2 for 4ti-0203 [http://4ti.co.uk/files/3113/4217/2464/4ti-0201.pdf]
    }

inputs['Lstated'] = inputs['Lstated'].tolist()
inputs['Pstated'] = inputs['Pstated'].tolist()

with open('inputs.json', 'w') as fp:
    json.dump(inputs, fp)

