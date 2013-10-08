
# FIXME How do we get line numbers of the failing tests from nosetests!?

# To run just these tests:
#
#     nosetests nipype.interfaces.minc

import os

from nipype.testing import (assert_equal, assert_true, assert_raises,
                            assert_not_equal, skipif)

import nipype.interfaces.minc as minc
from nipype.interfaces.minc import check_minc, no_minc

@skipif(no_minc)
def test_mincversion():
    ver = minc.Info.version()

    if ver:
        # If ver is None then MINC is not installed.

        # FIXME these version restrictions could be relaxed; they
        # refer to the latest version available that was installed by hand
        # instead of those available in the Debian ecosystem.

        yield assert_true, ver['minc']                  >= '2.2.00'
        yield assert_true, ver['libminc']               >= '2.2.00'
        yield assert_true, ver['netcdf'].split(' ')[0]  >= '4.1.3'
        yield assert_true, ver['hdf5']                  >= '1.8.8'
