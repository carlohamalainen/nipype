import os
import tempfile

from nipype.testing import (assert_equal, assert_true, assert_raises,
                            assert_not_equal, skipif, example_data)

import nipype.interfaces.minc as minc
from nipype.interfaces.minc import check_minc, no_minc
import numpy as np
import hashlib

def create_empty_temp_file(suffix='.mnc'):
    """
    Create an empty temporary file with the given suffix. It is the caller's
    responsibility to delete the file when done.
    """

    f = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
    f.close()
    return f.name

def remove(filename):
    """
    Delete a file; fail silently.
    """

    try:
        os.remove(filename)
    except OSError:
        pass

# FIXME skip import if pyminc not installed?
from pyminc.volumes.factory import *

def write_minc(filename, data):
    """
    Write a 2D or 3D NumPy array to the given file.
    """

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
    """
    Read the main data variable of a MINC file.
    """

    return volumeFromFile(filename).data

@skipif(no_minc)
def test_mincaverage_basic():
    """
    Test mincaverage, basic usage.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = create_empty_temp_file()
    file2 = create_empty_temp_file()
    output_file = create_empty_temp_file()

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


@skipif(no_minc)
def test_mincaverage_filelist():
    """
    Test mincaverage, pass list of files using -filelist. Also set output format to 'double'.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = create_empty_temp_file()
    file2 = create_empty_temp_file()
    output_file = create_empty_temp_file()

    filelist = create_empty_temp_file(suffix='.txt')
    open(filelist, 'w').write(file1 + '\n' + file2 + '\n')

    shape = (40, 30)

    data1 = np.random.random(shape)
    data2 = np.random.random(shape)
    data_avg_numpy = np.average(np.array([data1, data2]), axis=0)

    write_minc(file1, data1)
    write_minc(file2, data2)

    a = minc.AverageTask(filelist=filelist, output_file=output_file, clobber=True, format_double=True)
    a.run()

    data_avg_minc = read_minc(output_file)

    remove(file1)
    remove(file2)
    remove(output_file)

    yield assert_true, np.max(np.abs(data_avg_numpy - data_avg_minc)) < 1e-07

@skipif(no_minc)
def test_mincextract_ascii():
    """
    Test extract.
    """

    input_file  = example_data('minc/minc_test_2D_00.mnc')
    output_file = create_empty_temp_file(suffix='.raw')

    e = minc.ExtractTask(input_file=input_file, output_file=output_file)
    e.run()

    # remove(output_file)
    output = open(output_file).read()

    yield assert_true, len(output.split('\n')) == 201

    # FIXME will this have issues with roundoff on different platforms?
    yield assert_true, hashlib.md5(output).hexdigest() == '89941a49f2536cc114a5cfd4b5df0dd2'

def test_mincextract_bytes():
    """
    Test extract.
    """

    input_file  = example_data('minc/minc_test_2D_00.mnc')
    output_file = create_empty_temp_file(suffix='.raw')

    e = minc.ExtractTask(input_file=input_file, output_file=output_file, write_byte=True)
    e.run()

    # remove(output_file)
    output = open(output_file).read()

    yield assert_true, len(output) == 200

    # FIXME will this have issues with roundoff on different platforms?
    print hashlib.md5(output).hexdigest()
    yield assert_true, hashlib.md5(output).hexdigest() == '19afc43ef27c969afaa23e41bcefb78e'
