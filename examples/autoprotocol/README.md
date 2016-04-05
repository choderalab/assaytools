# Using AssayTools with autoprotocol

The examples here explore the use of `AssayTools` with [autoprotocol](http://autoprotocol.org/) as implemented in [`autoprotocol-python`](https://github.com/autoprotocol/autoprotocol-python).

[`autoprotcol-python`](https://github.com/autoprotocol/autoprotocol-python) can be installed via `conda`:
```bash
$ conda install -c omnia autoprotocol
```

## Manifest

* `example.py` - a (verbose) example of manually defining a plate's contents and data for the `autoprotocol` API
* `containers.py` - some common containers used in the laboratory

## Proposal

As discussed in our meeting, we'd like to generalize the model to permit each well to have arbitrary concentrations of multiple components, and to attach a competitive binding model for each ligand species.  Later, we could even provide a set of models to enable automated model selection.

For now, I'd like to propose an API for a general experiment, as well as a few convenience subclasses that make it easy to set up experiments of the kind we are currently performing.

### Specifying wells

We first need some way to describe, for each well:
* total well volume
* liquid level height (which might be computed from the well geometry and volume)
* concentrations (and uncertainties) of all components in the well

I wonder if we might want to piggypack on [`autoprotocol-python`](https://github.com/autoprotocol/autoprotocol-python) to specify these in a form that will be compatible with [autoprotocol](http://autoprotocol.org).

We can use the [Well](http://autoprotocol-python.readthedocs.org/en/latest/autoprotocol.html#id2) object to describe all of this information by standardizing the names of specific additional well properties (which can be added via [`well.add_properties()`](http://autoprotocol-python.readthedocs.org/en/latest/autoprotocol.html#well-add-properties) or set through [`well.set_properties()`](http://autoprotocol-python.readthedocs.org/en/latest/autoprotocol.html#well-set-properties)).  We can use a
[WellGroup](http://autoprotocol-python.readthedocs.org/en/latest/autoprotocol.html#container-wellgroup) to specify a collection of wells to pass to the Bayesian model for analysis.

An example of specifying a well using this approach would be
```python
# autoprotocol imports
from autoprotocol.container import Well, Container
from autoprotocol.container_type import ContainerType
from autoprotocol.unit import Unit

# Define the container type for 4titude 4ti-0110.
# info: http://4ti.co.uk/microplates/storage-plates/96-well-microplates/
# drawing: http://4ti.co.uk/files/6813/7510/4696/4ti_0110.pdf
# All arguments to ContainerType are required!
capabilities = ['spin', 'incubate'] # capabilities of container
well_diameter_mm = 6.50
#container_type = ContainerType(name='4titude 4ti-0110', is_tube=False, well_count=96, well_depth_mm=11.60, well_volume_ul=300, well_coating='polystyrene', sterile=False, capabilities=capabilities, shortname='4ti-0110', col_count=12, dead_volume_ul=5, safe_min_volume_ul=10)
container_type = ContainerType(name='4titude 4ti-0110', is_tube=False, well_count=96, well_depth_mm=Unit(11.60, 'millimeter'), well_volume_ul=Unit(300, 'milliliter'), well_coating='polystyrene', sterile=False, capabilities=capabilities, shortname='4ti-0110', col_count=12, dead_volume_ul=Unit(5,'milliliter'), safe_min_volume_ul=Unit(10, 'milliliter'))

# Generate a unique container ID
import uuid
id = str(uuid.uuid4())

# Define the container
container = Container(id=id, container_type=container_type)

# Define wells (the verbose way; we'd provide convenience methods to format the plate)
well = container.well("A1") # retrieve a specific well from the plate
well.set_volume(Unit(100, "milliliter")) # set the well volume

# Set concentrations of well components
concentrations = {
    'Src' : Unit(0.5, "micromoles/liter"),
    'bosutinib' : Unit(20, "micromoles/liter"),
    'imatinib' : Unit(10, "micromoles/liter")
}
well.set_properties({'concentrations' : concentrations})

# We would also need the optical path length through liquid, or well area.
well.set_properties({'area' : Unit(35, "millimeters**2")})
```
Since [`Unit`](http://autoprotocol-python.readthedocs.org/en/latest/autoprotocol.html#autoprotocol-unit) does not provide automatic unit conversion capabilities, we would have to work out a way to automatically convert between `simtk.unit` and `autoprotocol.unit.Unit`.

### Specifying plate reader data

[autoprotocol](http://autoprotocol.org) doesn't specify any particular format for the data output from plate reads, so we're stuck with coming up with our own.  We might want to create a class that ingests an Infinite XML file and provides a pythonic way to access data, e.g.

```python
>>> from assaytools.platereaders import InfiniteDataset
>>> data = InfiniteDataset(filename='data/2015-01-20 18-14-25_plate_1.xml')
>>> print data.welldata('A6')
{'TopFluorescenceRead' : 12451, 'BottomFluorescenceRead' : 74252, 'Absorbance' : 0.5632 }
```
where `InfiniteDataset` would be a subclass of a more generic `PlateReaderDataset`.

Because we haven't been super consistent about naming schemes (though we should be!) we may want to provide a mapping dict either here or in the analysis stage that maps the names in the XML file to names that we want in the dict returned by `data.welldata(id)`.

### General API for Bayesian analysis
We could specify the components of a well as a `dict` of concentrations.
```python
# Create a model
from assaytools.analysis import CompetitiveBindingAnalysis
receptor = 'Src' # receptor name (only one allowed)
ligands = ['bosutinib', 'bosutinib-isomer', 'imatinib'] # competitive ligand names
model = CompetitiveBindingAnalysis(wells=wellgroup, data=welldata, receptor=receptor, ligands=ligands)
# fit the maximum a posteriori (MAP) estimate
map = model.map_fit()
# run some MCMC sampling and return the MCMC object
mcmc = model.run_mcmc()
```


## Other notes

* [Transcriptic container types](https://developers.transcriptic.com/v1.0/docs/containers)
