AssayTools
==========

*Assay modeling and Bayesian analysis made easy.*

AssayTools is a python library that allows users to model automated assays and analyze assay data using powerful Bayesian techniques that allow complete characterization of the uncertainty in fit models.
With AssayTools, you can

 - Create models of experimental assays (e.g. fluorescence or absorbance assays of ligand binding from titration curves prepared by automated liquid handlers)
 - Analyze data from these assays using powerful Bayesian techniques
 - Derive parameters for these experimental assays from real data to use in modeling
 - Model new assay configurations to determine expected accuracy and bias across a range of ligand affinities

The library also ships with a command-line application for visualizing [Tecan Infinite](http://lifesciences.tecan.com/products/reader_and_washer/microplate_readers/infinite_m1000_pro) plate reader XML output.
When you install AssayTools, the script will be installed under the name ``xml2png``.

* [Download the code](https://github.com/choderalab/assaytools/releases/latest)
* [See it in action](https://github.com/choderalab/assaytools/tree/master/examples)
* [Get involved](https://github.com/assaytools/assaytools/issues)

.. raw:: html

  <!-- robot video -->
  <br/>
  <iframe id="player" type="text/html" width="500" height="300" style="display:block; margin:auto"
  src="http://www.youtube.com/embed/H6zIcMst4lY"/>
  frameborder="0"></iframe>
  <br/>

.. raw:: html

   <div style="display:none">


--------------------------------------------------------------------------------

Documentation
-------------

.. toctree::
   :maxdepth: 1

   getting_started
   theory

API Reference
-------------

.. toctree::
   :maxdepth: 1

   analysis

--------------------------------------------------------------------------------

   </div>

License
-------
AssayTools is licensed under the Lesser GNU General Public License (LGPL v2.1+).
