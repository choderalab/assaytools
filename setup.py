"""AssayTools: Assay modeling and Bayesian analysis made easy.

AssayTools is a python library that allows users to model automated assays and analyze assay data using powerful Bayesian techniques that allow complete characterization of the uncertainty in fit models.

"""

from __future__ import print_function, absolute_import
DOCLINES = __doc__.split("\n")

import sys
from setuptools import setup, Extension
sys.path.insert(0, '.')
from basesetup import (find_packages, write_version_py, build_ext,
                       StaticLibrary, CompilerDetection)

try:
    import numpy
except ImportError:
    print('Building and running assaytools requires numpy', file=sys.stderr)
    sys.exit(1)

try:
    # add an optional command line flag --no-install-deps to setup.py
    # to turn off setuptools automatic downloading of dependencies
    sys.argv.remove('--no-install-deps')
    no_install_deps = True
except ValueError:
    no_install_deps = False


setup_kwargs = {}
if 'setuptools' in sys.modules:
    setup_kwargs['zip_safe'] = False
    setup_kwargs['entry_points'] = {'console_scripts':
              ['xml2png = assaytools.scripts.xml2png:entry_point',
               'quickmodel = assaytools.scripts.quickmodel:entry_point',
               'ipnbdoctest = assaytools.scripts.ipnbdoctest:entry_point']}

    if sys.version_info[0] == 2:
        # required to fix cythoninze() for old versions of setuptools on
        # python 2
        m = sys.modules['setuptools.extension']
        m.Extension.__dict__ = m._Extension.__dict__
else:
    setup_kwargs['scripts'] = ['scripts/xml2png.py','scripts/quickmodel.py']


##########################
VERSION = "0.1.0"
ISRELEASED = False
__version__ = VERSION
##########################


CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)
Programming Language :: C
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

write_version_py(VERSION, ISRELEASED)
setup(name='assaytools',
      author='John D. Chodera',
      author_email='john.chodera@choderalab.org',
      description=DOCLINES[0],
      long_description="\n".join(DOCLINES[2:]),
      version=__version__,
      license='LGPLv2.1+',
      url='http://assaytools.choderalab.org',
      download_url = "https://github.com/choderalab/assaytools/releases/latest",
      platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
      classifiers=CLASSIFIERS.splitlines(),
      packages=find_packages(),
      package_dir={'assaytools': 'AssayTools', 'assaytools.scripts': 'scripts'},
      cmdclass={'build_ext': build_ext},
      package_data={'assaytools.html': ['static/*']},
      install_requires=[
#        'numpy',
#        'pandas',
#        'pymc',
#        'matplotlib',
#        'seaborn',
#        'lxml'
        ],
      **setup_kwargs)
