AssayTools
==========

*Assay modeling and Bayesian analysis made easy.*

AssayTools is a python library that allows users to model automated assays and analyze assay data using powerful Bayesian techniques that allow complete characterization of the uncertainty in fit models.
With AssayTools, you can

 - Create models of experimental assays (e.g. fluorescence assays of ligand binding from titration curves prepared by automated liquid handlers)
 - Analyze data from these assays usijng powerful Bayesian techniques
 - Derive parameters for these experimental assays from real data to use in modeling
 - Model new assay configurations to determine expected accuracy and bias across a range of ligand affinities

The library also ships with a command-line application for visualizing Tecan Infinite plate XML output.
When you install AssayTools, the script will be installed under the name ``xml2png``.

.. raw:: html

  <div>
      <h2 style="display: inline; float:left; margin-left:5em">
          <a href="https://github.com/choderalab/assaytools/releases/latest">
          Download the Code</a>
      </h2>
      <h2 style="display: inline; float:right; margin-right:5em">
         <a href="https://github.com/choderalab/assaytools/tree/master/examples/fluorescence-binding-assay">
          See it in Action</a>
      </h2>
      <div style="clear:both"></div>
      <h2 style="display: inline; float:right; margin-right:7em">
          <a href="https://github.com/assaytools/assaytools/issues">
          Get Involved</a>
      </h2>
  </div>
  <br/>
  <iframe id="player" type="text/html" width="500" height="300" style="display:block; margin:auto"
  src="http://www.youtube.com/embed/H6zIcMst4lY"/>
  frameborder="0"></iframe>



.. raw:: html

   <div style="display:none">


--------------------------------------------------------------------------------

Documentation
-------------

.. toctree::
   :maxdepth: 1

   getting_started
   theory
..   examples/index
   whatsnew
   faq
..   Discussion Forums <http://discourse.mdtraj.org>

API Reference
-------------

.. toctree::
   :maxdepth: 1

..   load_functions
..   atom_selection
   analysis
..   api/trajectory_files
..   api/reporters
..   api/utils
   xml2png


Developing
----------

.. toctree::
   :maxdepth: 1

   building_docs


--------------------------------------------------------------------------------

.. raw:: html

   </div>

License
-------
AssayTools is licensed under the Lesser GNU General Public License (LGPL v2.1+).
