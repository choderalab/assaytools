from autoprotocol.container import Well, WellGroup, Container
from autoprotocol.container_type import ContainerType
from autoprotocol.unit import Unit
import numpy as np
import random, string

#
# Define simple solutions
# Solutions are the fundamental liquids that have concentrations (and concentration uncertainties) for a single species.
# In future, we will associate densities, allow mixing of solutions, propagate uncertainty.
# Solutions have unknown *true* concentrations, except for buffer solutions which have exactly zero concentration.
#

class Buffer(object):
    def __init__(self, name):
        """
        """
        self.name = name
        self.species = None
        self.concentration = Unit(0.0, 'moles/liter')
        self.uncertainty = Unit(0.0, 'moles/liter')

DMSO = Buffer('DMSO')

class Solution(object):
    """A solution of a single component with a defined concentration and uncertainty.
    """
    def __init__(self, name, species, buffer, concentration, uncertainty):
        """
        Parameters
        ----------
        name : str
            A descrptive name for the solution.
        species : str
            The name of the species dissolved in buffer.
        buffer : Buffer
            The buffer the solution is prepared in.
        concentration : Unit
            The concentration of the species in the solution.
        uncertainty : Unit
            The concentration uncertainty of the species in the solution.
        """
        self.name = name
        self.species = species
        self.buffer = buffer
        self.concentration = concentration
        self.uncertainty = uncertainty

class ProteinSolution(Solution):
    """A protein solution in buffer prepared spectrophotometrically.
    """
    def __init__(self, name, species, buffer, absorbance, extinction_coefficient, molecular_weight, ul_protein_stock, ml_buffer):
        """
        Parameters
        ----------
        protein_name : str
            A descreptive name of the solution.
        species : str
            The protein name
        buffer : BufferSolution
            The corresponding buffer solution.
        absorbance : float
            Absorbance for 1 cm equivalent path length
        extinction_coefficient : float
            Extinction coefficient (1/M/cm)
        molecular_weight : float
            Molecular weight (g/mol or Daltons)
        ul_protein_stock : float
            Microliters of protein stock added to make solution
        ml_buffer : float
            Milliliters of buffer added to make solution

        """
        self.name = name + ' in ' + buffer.name
        self.species = species
        self.buffer = buffer
        concentration = (absorbance / extinction_coefficient) * (ul_protein_stock/1000.0) / ml_buffer # M
        spectrophotometer_CV = 0.10 # TODO: Specify in assumptions YAML file
        self.concentration = Unit(concentration, 'moles/liter')
        self.uncertainty = spectrophotometer_CV * self.concentration # TODO: Compute more precisely

class DMSOStockSolution(Solution):
    """A stock solution representing a compound dissolved in DMSO.
    """
    def __init__(self, dmso_stock):
        """
        Parameters
        ----------
        dmso_stock : dict
            The dictionary containing 'id', 'compound_name', 'compound mass (mg)', 'molecular weight', 'purity', 'solvent_mass'
        """
        self.name = '10 mM ' + dmso_stock['compound name'] + ' DMSO stock'
        self.species = dmso_stock['compound name']
        dmso_density = 1.1004 # g/cm3
        mass_uncertainty = 0.01 # TODO: Calculate from balance precision
        concentration = dmso_stock['compound mass (mg)'] / 1000 * dmso_stock['molecular weight'] * dmso_stock['purity'] / (dmso_stock['solvent mass (g)'] * dmso_density / 1000) # mol/liter
        self.concentration = Unit(concentration, 'moles/liter')
        self.uncertainty = mass_uncertainty * self.concentration
        self.buffer = DMSO

def DMSOStockSolutions(dmso_stocks_csv_filename):
    """
    Create DMSO stock solutions.
    """

    #
    # Load information about DMSO stocks
    # In future, this might come from a database query that returns JSON.
    #

    import csv
    dmso_stocks_csv_filename = 'DMSOstocks-Sheet1.csv'
    dmso_stocks = dict()
    with open(dmso_stocks_csv_filename, 'rb') as csvfile:
         csvreader = csv.DictReader(csvfile, delimiter=',', quotechar='|')
         for row in csvreader:
             if row['id'] != '':
                for key in ['compound mass (mg)', 'purity', 'solvent mass (g)', 'molecular weight']:
                    row[key] = float(row[key])
                dmso_stocks[row['id']] = row

    solutions = dict()
    solutions['DMSO'] = DMSO # Add DMSO
    for id in dmso_stocks:
        solutions[id] = DMSOStockSolution(dmso_stocks[id])

    return solutions

#
# Define container types
#

def define_container_types():
    #
    # Define assay plate container
    #

    container_types = dict()

    # Define the container type for 4titude 4ti-0223.
    # info: http://4ti.co.uk/microplates/black-clear-bottom/96-well/
    # drawing: http://4ti.co.uk/files/1614/0542/7662/4ti-0223_Marketing_Drawing.pdf
    # All arguments to ContainerType are required!
    capabilities = ['pipette', 'spin', 'absorbance', 'fluorescence', 'luminescence', 'incubate', 'gel_separate', 'cover', 'seal', 'stamp', 'dispense']
    container_type = ContainerType(name='4titude 4ti-0223', is_tube=False, well_count=96, well_depth_mm=Unit(11.15, 'millimeter'), well_volume_ul=Unit(300, 'milliliter'), well_coating='polystyrene', sterile=False, capabilities=capabilities, shortname='4ti-0223', col_count=12, dead_volume_ul=Unit(20,'milliliter'), safe_min_volume_ul=Unit(50, 'milliliter'))

    # Attach well area.
    well_diameter = Unit(6.30, "millimeters")
    well_area = np.pi * (well_diameter/2)**2
    setattr(container_type, 'well_area', well_area)

    container_types[container_type.name] = container_type

    return container_types

container_types = define_container_types()

#
# Helper functions for our common assays
#

def generate_uuid(size=6, chars=(string.ascii_uppercase + string.digits)):
    """
    Generate convenient universally unique id (UUID)

    Parameters
    ----------
    size : int, optional, default=6
       Number of alphanumeric characters to generate.
    chars : list of chars, optional, default is all uppercase characters and digits
       Characters to use for generating UUIDs

    NOTE
    ----
    This is not really universally unique, but it is convenient.

    """
    return ''.join(random.choice(chars) for _ in range(size))

def provision_assay_plate(name, plate_type='4titude 4ti-0223', id=None):
    """
    Provision a new assay plate.

    Parameters
    ----------
    name : str
       The name of the container
    plate_type : str, optional, default='4titude 4ti-0223'
       The name of the plate type used to retrieve the `container_type` from library
    id : str, optional, default=None
       Unless `id` is specified, a unique container ID will be autogenerated.
    """
    if id == None:
        id = generate_uuid

    # Define the container
    container_type = container_types[plate_type]
    container = Container(name="assay-plate", id=id, container_type=container_type)

    # Initialize well properties for this container
    for well in container.all_wells():
        well.set_volume(Unit(0.0, 'microliters')) # well starts empty

    return container

def dispense_evo(container, solution, volume, rows):
    """
    Dispense a given volume into the specified rows with a Tecan EVO pipetting robot.

    Parameters
    ----------
    container : autoprotocol.containers.Container
       The container to dispense the specified solution into
    solution : Solution
       The solution to dispense
    volume : Unit compatible with 'microliters'
       The volume to dispense
    rows : list of char
       The rows to disepnse into (e.g. ['A', 'C', 'E', 'G'])

    TODO
    ----
    * Are `contents` stored by solution name?
    * Generalize beyond only specifying 'rows'
    * Handle case where the same solution is dispensed multiple times.
    """
    CV = 0.004 # for 100 uL volume; TODO: compute this automatically based on dispense volume
    for well in container.all_wells():
        if 'contents' not in well.properties:
            well.properties['contents'] = dict()
        contents = well.properties['contents']

        well_name = well.humanize()
        if well_name[0] in rows:
            contents[solution.name] = (volume, CV * volume) # volume, error
            well.set_volume(well.volume + volume)

def dispense_hpd300(container, solutions, xml_filename, plate_index=0):
    """
    Dispense one or more solutions into a container using an HP D300.

    Parameters
    ----------
    container :

    solutions : list of Solutions
       List of solutils corresponding to the <Fluids/> block in the HP D300 XML filename
    xml_filename : str
       HP D300 simulation DATA XML filename
    plate_index : int, optional, default=0
       If the XML file contains multiple <Plate> tags, use the specified index.
    """
    CV = 0.08 # CV for HP D300

    # Read HP D300 XML file
    import xml.etree.ElementTree as ET
    tree = ET.parse(xml_filename)
    root = tree.getroot()

    # Read fluids
    fluids = root.findall('./Fluids/Fluid')

    # Rewrite fluid names to match stock names.
    print fluids
    print solutions
    for (index, solution) in enumerate(solutions):
        fluids[index].attrib['Name'] = solution.name

    def humanize_d300_well(row, column):
        """
        Return the humanized version of a D300 well index.
        """
        return chr(row + 97) + str(column+1)

    # Read dispensed volumes from the specified Plate entry
    dispensed = root.findall('./Dispensed')[0]
    volume_unit = dispensed.attrib['VolumeUnit'] # dispensed volume unit
    hpd300_wells = dispensed.findall("Plate[@Index='%d']/Well" % plate_index)
    for hpd300_well in hpd300_wells:
        # Get corresponding container well.
        row_index = int(hpd300_well.attrib['R'])
        column_index = int(hpd300_well.attrib['C'])
        well_name = humanize_d300_well(row_index, column_index)

        # Retrieve well from container.
        well = container.well(well_name)
        if 'contents' not in well.properties:
            well.properties['contents'] = dict()
        contents = well.properties['contents']

        # Handle dispensed fluids.
        # TODO: Deal with case where we might dispense the same fluid twice.
        dispensed_fluids = hpd300_well.findall('Fluid')
        for dispensed_fluid in dispensed_fluids:
            dispensed_fluid_index = int(dispensed_fluid.attrib['Index'])
            fluid_name = fluids[dispensed_fluid_index].attrib['Name']
            total_volume = Unit(float(dispensed_fluid.attrib['TotalVolume']), volume_unit)
            detail = dispensed_fluid.findall('Detail')

            # Gather details
            time = detail[0].attrib['Time']
            cassette = int(detail[0].attrib['Cassette'])
            head = int(detail[0].attrib['Head'])
            dispensed_volume = Unit(float(detail[0].attrib['Volume']), volume_unit)

            # Fill well
            contents[fluid_name] = (dispensed_volume, CV * dispensed_volume) # (volume, error)
            well.set_volume(well.volume + dispensed_volume)


def read_infinite(container, xml_filename):
    """
    Read measurements from a Tecan Infinite reader XML file.

    Measurements are stored in `well.properties['measurements']` under
    * `absorbance` : e.g. { '280:nanometers' : 0.425 }
    * `fluorescence` : e.g. { ('280:nanometers', '350:nanometers', 'top') : 15363 }

    TODO
    ----
    * Eventually, we can read all components of the Infinite file directly here.
    """
    # TODO: Replace read_icontrol_xml with direct processing of XML file
    from assaytools import platereader
    data = platereader.read_icontrol_xml(xml_filename)

    for well in container.all_wells():
        well_name = well.humanize()

        # Attach plate reader data
        # TODO: Only process wells for which measurements are available
        if 'measurements' not in well.properties:
            well.properties['measurements'] = dict()
        measurements = well.properties['measurements']

        for key in data:
            if key.startswith('Abs_'):
                # absorbance
                [prefix, wavelength] = key.split('_')
                wavelength = wavelength + ':nanometers'
                if 'absorbance' not in measurements:
                    measurements['absorbance'] = dict()
                measurements['absorbance'][wavelength] = float(data[key][well_name])
            elif (key.endswith('_TopRead') or key.endswith('_BottomRead')):
                # top fluorescence read
                [wavelength, suffix] = key.split('_')
                excitation_wavelength = wavelength + ':nanometers'
                emission_wavelength = '450:nanometers'
                if key.endswith('_TopRead'):
                    geometry = 'top'
                else:
                    geometry = 'bottom'
                if 'fluorescence' not in measurements:
                    measurements['fluorescence'] = dict()
                measurements['fluorescence'][(excitation_wavelength, emission_wavelength, geometry)] = float(data[key] [well_name])


class Assay(object):
    """
    Assay base class
    """
    pass

params = {
    'd300_xml_filename' : 'Src_Bos_Ima_96well_Mar2015 2015-03-07 1736.DATA.xml' # HP D300 dispense simulated DATA file
    'infinite_xml_filename' : 'Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml' # Tecan Infinite plate reader output data
    'dmso_stocks_csv_filename' : 'DMSOstocks-Sheet1.csv' # CSV file of DMSO stock inventory
    'hpd300_fluids' : ['GEF001', 'IMA001', 'DMSO'] # uuid of DMSO stocks from dmso_stocks_csv_filename (or 'DMSO' for pure DMSO) used to define HP D300 XML <Fluids> block
    'receptor_species' = 'Abl(D382N)' # receptor name (just used for convenience)
    'protein_absorbance' = 4.24 # absorbance reading of concentrated protein stock before dilution
    'protein_extinction_coefficient' = Unit(49850, '1/molar/centimeter') # 1/M/cm extinction coefficient for protein
    'protein_molecular_weight' = Unit(41293.2, 'daltons') # g/mol protein molecular weight
    'protein_stock_volume' = Unit(165.8, 'microliters') # uL protein stock solution used to make 1 uM protein stock
    'buffer_volume' = Unit(14.0, 'milliliters') # mL buffer used to make 1 uM protein stock
    'rows_to_analyze' = ['A', 'B'] # rows to analyze
    'assay_volume' = Unit(100.0, 'microliters') # quantity of protein or buffer dispensed into plate
}

class SingletAssay(Assay):
    def __init__(self,
        d300_xml_filename,
        infinite_xml_filename,
        dmso_stocks_csv_filename,
        hpd300_fluids,
        receptor_species,
        protein_absorbance,
        protein_extinction_coefficient,
        protein_molecular_weight,
        protein_stock_volume,
        buffer_volume,
        rows_to_analyze,
        assay_volume):
        """
        Set up a single-point assay.

        Parameters
        ----------
        d300_xml_filename : str
            HP D300 dispense simulated DATA file
        infinite_xml_filename : str
            Tecan Infinite plate reader output data, e.g. 'Abl Gef gain 120 bw1020 2016-01-19 15-59-53_plate_1.xml'
        dmso_stocks_csv_filename : str
            CSV file of DMSO stock inventory, e.g. 'DMSOstocks-Sheet1.csv'
        hpd300_fluids : list of str
            uuids of DMSO stocks from dmso_stocks_csv_filename (or 'DMSO' for pure DMSO) used to define HP D300 XML <Fluids> block, e.g. ['GEF001', 'IMA001', 'DMSO']
        receptor_species
            Name of receptor species, e.g. 'Abl(D382N)'
        protein_absorbance
            Absorbance reading of concentrated protein stock before dilution
        protein_extinction_coefficient : Unit compatible with 1/molar/centimeters
            Extinction coefficient for protein, e.g. Unit(49850, '1/molar/centimeter')
        protein_molecular_weight : Unit compatible with daltons
            Protein molecular weight
        protein_stock_volume : Unit compatible with microliters
            Volume of high-concentration protein stock solution used to make 1 uM protein stock
        buffer_volume : Unit compatible with milliliters
            Volume of buffer used to make ~1 uM protein stock used to fill wells
        rows_to_analyze : list
            Rows to analyze, e.g. ['A', 'B']
        assay_volume : Unit compatible with microliters
            Quantity of protein or buffer dispensed into plate

        """
        # Read DMSO stock solutions from inventory CSV file
        from assaytools.experiments import DMSOStockSolutions, DMSO, Buffer, ProteinSolution
        solutions = DMSOStockSolutions(dmso_stocks_csv_filename) # all solutions from DMSO stocks inventory

        # Enumerate all ligand species from DMSO stocks.
        ligand_species = set( [ solution.species for solution in solutions.values() if (solution.species != None)] )

        # Add buffer and protein stock solutions
        solutions['buffer'] = Buffer(name='20 mM Tris buffer')
        solutions[receptor_species] = ProteinSolution(name='1 uM %s' % receptor_species, species=receptor_species, buffer=solutions['buffer'],
        absorbance=protein_absorbance, extinction_coefficient=protein_extinction_coefficient, molecular_weight=protein_molecular_weight, ul_protein_stock=ul_protein_stock, ml_buffer=ml_buffer)

        # Populate the Container data structure with well contents and measurements
        from assaytools.experiments import provision_assay_plate, dispense_evo, dispense_hpd300, read_infinite
        plate = provision_assay_plate(name='assay-plate', plate_type='4titude 4ti-0223')
        dispense_evo(plate, solution=solutions[receptor_species], volume=assay_volume, rows=['A', 'C', 'E', 'G'])
        dispense_evo(plate, solution=solutions['buffer'], volume=assay_volume, rows=['B', 'D', 'F', 'H'])
        dispense_hpd300(plate, solutions=[solutions[id] for id in hpd300_fluids], xml_filename=d300_xml_filename)
        read_infinite(plate, xml_filename=infinite_xml_filename)

        # Select specified rows for analysis.
        from autoprotocol import WellGroup
        well_group = WellGroup([well for well in plate.all_wells() if (well.humanize()[0] in rows_to_analyze)])

        # Create a model
        from assaytools.analysis import CompetitiveBindingAnalysis
        self.experiment = CompetitiveBindingAnalysis(solutions=solutions, wells=well_group, receptor_name=receptor_species, ligand_names=ligand_species)
