
# FIXME How do we get line numbers of the failing tests from nosetests!?

import os

from nipype.testing import (assert_equal, assert_true, assert_raises,
                            assert_not_equal, skipif)

import nipype.interfaces.minc as minc
from nipype.interfaces.minc import check_minc, no_minc

# FIXME write some tests!
