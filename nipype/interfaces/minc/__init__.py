# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The minc module provides classes for interfacing with the `MINC
<http://www.bic.mni.mcgill.ca/ServicesSoftware/MINC>`_ command line tools.

Top-level namespace for minc.
"""

from .base import (Info, check_minc, no_minc)
from .utils import (ToRaw, Convert, Copy,
                    ToEcat, Dump, Average,
                    Calc, Extract, Blob,
                    Blur, Math)
