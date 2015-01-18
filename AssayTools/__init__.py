"""AssayTools is a python library that allows users to model automated assays and analyze assay data using powerful Bayesian techniques that allow complete characterization of the uncertainty in fit models.
"""

from assaytools import models

def test(label='full', verbose=2):
    """Run tests for mdtraj using nose.

    Parameters
    ----------
    label : {'fast', 'full'}
        Identifies the tests to run. The fast tests take about 10 seconds,
        and the full test suite takes about two minutes (as of this writing).
    verbose : int, optional
        Verbosity value for test outputs, in the range 1-10. Default is 2.
    """
    import assaytools
    from assaytools.testing.nosetester import AssaytoolsTester
    tester = AssaytoolsTester(assaytools)
    return tester.test(label=label, verbose=verbose, extra_argv=('--exe',))
# prevent nose from discovering this function, or otherwise when its run
# the test suite in an infinite loop
test.__test__ = False


def capi():
    import os
    import sys
    module_path = sys.modules['assaytools'].__path__[0]
    return {
        'lib_dir':  os.path.join(module_path, 'core', 'lib'),
        'include_dir': os.path.join(module_path, 'core', 'lib'),
    }
