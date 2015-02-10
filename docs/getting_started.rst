.. _getting-started:

************
Installation
************

Platforms
=========

AssayTools runs on CPython 2.7, 3.3, and 3.4.

Install with Conda
------------------
.. _install-with-conda:

`conda <http://www.continuum.io/blog/conda>`_ is a python package manager built for scientific python. Unlike ``easy_install`` or ``pip``, it handles binaries and binary dependencies, which are critical for most scientific workflows. If you're a ``conda`` user, you can install MDTraj by adding the omnia channel. If you're not a conda user, you should look into it. ::

To install MDTraj with conda, use the following commands ::

  $ conda config --add channels http://conda.binstar.org/choderalab
  $ conda install assaytools

.. note:: ``conda`` will automatically all of the tricky dependencies from binary packages automatically! This includes pytables / numpy / scipy / pandas / pymc! The easiest way to get conda is with the `Anaconda python distribution <https://store.continuum.io/cshop/anaconda/>`_.

Install from Source
-------------------
Clone the source code repository from github::

  $ git clone git://github.com/choderalab/assaytools.git

If you don't have ``git``, you can download the source code as a zip file from
https://github.com/rmcgibbo/mdtraj/archive/master.zip. Then, in the directory containing the source code, you can install it with. ::

  $ python setup.py install

Dependencies
============

To use mdtraj, the following libraries and software will need to be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Linux and Mac machines. Windows is not
        well supported.

    `Python <http://python.org>`_ >= 2.6
        The development package (``python-dev`` or ``python-devel``
        on most Linux distributions) is recommended.

    `NumPy <http://numpy.scipy.org/>`_ >= 1.6.0
        Numpy is the base package for numerical computing in python.

    `PyMC <http://pymc-devs.github.io/pymc/>`_ >= 2.3.3
        PyMC is used as the underlying Markov chain Monte Carlo (MCMC) based Bayesian inference engine for AssayTools.

Optional packages:

    `Pandas <http://pandas.pydata.org>`_ >= 0.9.0
        Some functionality requires pandas.

Avoid Hassles with Anaconda or Canopy
-------------------------------------

The easiest way to get all of the dependencies is to install one of the
pre-packaged scientific python distributes like `Enthought's Canopy
<https://www.enthought.com/products/canopy/>`_ or `Continuum's Anaconda
<https://store.continuum.io/>`_. These distributions already contain all of
the dependences, and are distributed via 1-click installers for Windows, Mac
and Linux.

.. note:: We recommend Continuum's Anaconda.

Manually Installing the Dependencies
------------------------------------

Oh, please don't!

Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``pip`` if you don't already have it. ::

  pip install nose

Then, to run the tests, open a python shell and do ::

  >>> import assaytools
  >>> assaytools.test()

From the source directory, you can also run the tests with ``nosetests`` on
the command line.
