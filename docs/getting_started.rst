.. _getting-started:

************
Installation
************

.. highlight:: bash

Install with conda
------------------
.. _install-with-conda:

To install `AssayTools` and its dependencies with `conda <conda.pydata.org>`_`, use the following commands:

::

   $ conda install -c omnia assaytools

Install from source
-------------------
.. _install-from-source

You can install the latest development version of ``AssayTools`` from github via ``pip``:

::

   $ pip install git+https://github.com/choderalab/assaytools.git

Testing Your Installation
=========================
.. _testing-your-installation

Running the tests is a great way to verify that everything is working.
The test suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can install via ``conda``:

::

   $ conda install nose

Then, to run the tests:

::

   $ nosetests -vv assaytools
