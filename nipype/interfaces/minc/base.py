# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The minc module provides classes for interfacing with the `MINC
<http://www.bic.mni.mcgill.ca/ServicesSoftware/MINC>`_ command line tools.  This
was written to work with MINC version 2.2.00.

TODO Documentation, docstrings, etc.

Author: Carlo Hamalainen <carlo@carlo-hamalainen.net>
        http://carlo-hamalainen.net
"""


# TODO Check exit-code behaviour of minc commands.

from nipype.interfaces.base import (
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    StdOutCommandLineInputSpec,
    StdOutCommandLine,
    File,
    InputMultiPath,
    traits,
)

import os

import warnings
warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)

def check_minc():
    """Returns True if and only if MINC is installed.'
    """

    return Info.version() is not None

def no_minc():
    """Returns True if and only if MINC is *not* installed.
    """
    return not check_minc()

class Info(object):
    """Handle MINC version information.

    version refers to the version of MINC on the system
    """

    @staticmethod
    def version():
        """Check for minc version on the system

        Parameters
        ----------
        None

        Returns
        -------
        version : dict
           Version number as dict or None if MINC not found

        """
        try:
            clout = CommandLine(command='mincinfo',
                                args='-version',
                                terminal_output='allatonce').run()
        except IOError:
            return None

        out = clout.runtime.stdout

        def read_program_version(s):
            if 'program' in s:
                return s.split(':')[1].strip()
            return None

        def read_libminc_version(s):
            if 'libminc' in s:
                return s.split(':')[1].strip()
            return None

        def read_netcdf_version(s):
            if 'netcdf' in s:
                return ' '.join(s.split(':')[1:]).strip()
            return None

        def read_hdf5_version(s):
            if 'HDF5' in s:
                return s.split(':')[1].strip()
            return None

        versions = {'minc':    None,
                    'libminc': None,
                    'netcdf':  None,
                    'hdf5':    None,
                   }

        for l in out.split('\n'):
            for (name, f) in [('minc',      read_program_version),
                              ('libminc',   read_libminc_version),
                              ('netcdf',    read_netcdf_version),
                              ('hdf5',      read_hdf5_version),
                             ]:
                if f(l) is not None: versions[name] = f(l)

        return versions
