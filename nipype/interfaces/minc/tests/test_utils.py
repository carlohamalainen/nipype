
# FIXME How do we get line numbers of the failing tests from nosetests!?

import os
import tempfile

from nipype.testing import (assert_equal, assert_true, assert_raises,
                            assert_not_equal, skipif)

import nipype.interfaces.minc as minc
from nipype.interfaces.minc import check_minc, no_minc
import numpy as np

def touch(suffix='.mnc'):
    f = tempfile.NamedTemporaryFile(suffix='.mnc', delete=False)
    f.close()
    return f.name

def remove(filename):
    try:
        os.remove(filename)
    except OSError:
        pass
 
# FIXME skip import if pyminc not installed?
from pyminc.volumes.factory import *

def write_minc(filename, data):
    # FIXME doc

    nr_dims = len(data.shape)
    assert 1 <= nr_dims <= 3

    dim_names = ['xspace', 'yspace', 'zspace'][:nr_dims][::-1]
    starts = [-5, -5, -5][:nr_dims][::-1]
    steps = [0.1, 0.2, 0.3][:nr_dims][::-1]

    vol = volumeFromDescription(filename, dim_names, data.shape, starts, steps, volumeType='int')
    vol.data = data
    vol.writeFile()
    vol.closeVolume()

def read_minc(filename):
    # FIXME doc
    return volumeFromFile(filename).data

@skipif(no_minc)
def test_mincaverage():
    # FIXME doc
    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = touch()
    file2 = touch()
    output_file = touch()

    shape = (40, 30)

    data1 = np.random.random(shape)
    data2 = np.random.random(shape)
    data_avg_numpy = np.average(np.array([data1, data2]), axis=0)

    write_minc(file1, data1)
    write_minc(file2, data2)

    a = minc.AverageTask(input_files=[file1, file2], output_file=output_file, clobber=True)
    a.run()

    data_avg_minc = read_minc(output_file)

    remove(file1)
    remove(file2)
    remove(output_file)

    yield assert_true, np.max(np.abs(data_avg_numpy - data_avg_minc)) < 1e-07
