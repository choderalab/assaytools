.. _getting-started:

************
Installation
************

Install with conda
------------------
.. _install-with-conda:

To install `AssayTools` with `conda <conda.pydata.org>`_`, use the following commands ::

  $ conda install -c omnia assaytools

.. note::
   ``conda`` will automatically all of the tricky dependencies from binary packages automatically!
   This includes :mod:`pytables`, :mod:`numpy`, :mod:`scipy`, :mod:`pandas`, and :mod:`pymc`!
   The easiest way to get ``conda`` is with the `Anaconda python distribution <https://store.continuum.io/cshop/anaconda/>`_, or its smaller version `Miniconda <http://conda.pydata.org/miniconda.html>`_.

Install from source
-------------------

You can install the latest development version of ``AssayTools`` from github via ``pip``:

  $ pip install git+https://github.com/choderalab/assaytools.git

Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working.
The test suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick up via ``conda`` if you don't already have it. ::

  $ conda install nose

Then, to run the tests:

  $ nosetests -vv assaytools
