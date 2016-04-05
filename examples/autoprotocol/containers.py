from autoprotocol.container import Well, Container
from autoprotocol.container_type import ContainerType
from autoprotocol.unit import Unit
import numpy as np

# Define the container type for 4titude 4ti-0110.
# info: http://4ti.co.uk/microplates/storage-plates/96-well-microplates/
# drawing: http://4ti.co.uk/files/6813/7510/4696/4ti_0110.pdf
# All arguments to ContainerType are required!
capabilities = ['pipette', 'spin', 'incubate', 'gel_separate', 'cover', 'stamp', 'seal', 'dispense'] # capabilities of container
well_diameter = Unit(6.50, "millimeters")
container_type = ContainerType(name='4titude 4ti-0110', is_tube=False, well_count=96, well_depth_mm=Unit(11.60, 'millimeter'), well_volume_ul=Unit(300, 'milliliter'), well_coating='polystyrene', sterile=False, capabilities=capabilities, shortname='4ti-0110', col_count=12, dead_volume_ul=Unit(5,'milliliter'), safe_min_volume_ul=Unit(10, 'milliliter'))

# Define the container type for 4titude 4ti-0223.
# info: http://4ti.co.uk/microplates/black-clear-bottom/96-well/
# drawing: http://4ti.co.uk/files/1614/0542/7662/4ti-0223_Marketing_Drawing.pdf
# All arguments to ContainerType are required!
capabilities = ['pipette', 'spin', 'absorbance', 'fluorescence', 'luminescence', 'incubate', 'gel_separate', 'cover', 'seal', 'stamp', 'dispense']
well_diameter = Unit(6.30, "millimeters")
well_area = np.pi * (well_diameter/2)**2
container_type = ContainerType(name='4titude 4ti-0223', is_tube=False, well_count=96, well_depth_mm=Unit(11.15, 'millimeter'), well_volume_ul=Unit(300, 'milliliter'), well_coating='polystyrene', sterile=False, capabilities=capabilities, shortname='4ti-0223', col_count=12, dead_volume_ul=Unit(20,'milliliter'), safe_min_volume_ul=Unit(50, 'milliliter'))
