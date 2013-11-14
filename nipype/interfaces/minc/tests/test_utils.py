import os
import tempfile

from nipype.testing import (assert_equal, assert_true, assert_raises,
                            assert_not_equal, skipif, example_data)

import nipype.interfaces.minc as minc
from nipype.interfaces.minc import check_minc, no_minc
import numpy as np
import hashlib

from nipype.interfaces.base import (isdefined)

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

    a = minc.Average(input_files=[file1, file2], output_file=output_file)
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

    a = minc.Average(filelist=filelist, output_file=output_file, format_double=True)
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

    e = minc.Extract(input_file=input_file, output_file=output_file)
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

    e = minc.Extract(input_file=input_file, output_file=output_file, write_byte=True)
    e.run()

    # remove(output_file)
    output = open(output_file).read()

    yield assert_true, len(output) == 200

    # FIXME will this have issues with roundoff on different platforms?
    print hashlib.md5(output).hexdigest()
    yield assert_true, hashlib.md5(output).hexdigest() == '19afc43ef27c969afaa23e41bcefb78e'

@skipif(no_minc)
def test_minccalc_add():
    """
    Test mincaverage, add two files together.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = create_empty_temp_file()
    file2 = create_empty_temp_file()
    output_file = create_empty_temp_file()

    shape = (40, 30)

    data1 = np.random.random(shape)
    data2 = np.random.random(shape)

    write_minc(file1, data1)
    write_minc(file2, data2)

    a = minc.Calc(input_files=[file1, file2], output_file=output_file, expression='A[0] + A[1]')
    a.run()

    data_minccalc_output = read_minc(output_file)

    remove(file1)
    remove(file2)
    remove(output_file)

    yield assert_true, np.max(np.abs((data1 + data2) - data_minccalc_output)) < 1e-07

@skipif(no_minc)
def test_minccalc_sub():
    """
    Test mincaverage, subtract two files.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = create_empty_temp_file()
    file2 = create_empty_temp_file()
    output_file = create_empty_temp_file()

    shape = (40, 30)

    data1 = np.random.random(shape)
    data2 = np.random.random(shape)

    write_minc(file1, data1)
    write_minc(file2, data2)

    a = minc.Calc(input_files=[file1, file2], output_file=output_file, expression='A[0] - A[1]')
    a.run()

    data_minccalc_output = read_minc(output_file)

    remove(file1)
    remove(file2)
    remove(output_file)

    yield assert_true, np.max(np.abs((data1 - data2) - data_minccalc_output)) < 1e-07

@skipif(no_minc)
def test_minccalc_sumsquares():
    """
    Test mincaverage, sum of squares.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    shape = (40, 30)

    nr_files = 5

    filenames = [create_empty_temp_file() for _ in range(nr_files)]
    data      = [np.random.random(shape)  for _ in range(nr_files)]

    for i in range(nr_files): write_minc(filenames[i], data[i])

    output_file = create_empty_temp_file()

    a = minc.Calc(input_files=filenames, output_file=output_file,
                      expression='total = 0; for {i in [0:(len(A)-1)]} { total = total + A[i]^2 }; total')
    a.run()

    data_minccalc_output = read_minc(output_file)

    for f in filenames: remove(f)
    remove(output_file)

    sum_squares = np.sum([x**2 for x in data], axis=0)
    yield assert_true, np.max(np.abs(sum_squares - data_minccalc_output)) < 1e-06

@skipif(no_minc)
def test_minccalc_std():
    expression_string = r"""s0 = s1 = s2 = 0;

    for { i in [0:len(A)) } {
        v=A[i];
        if (!isnan(v)) {
            s0 = s0 + 1;
            s1 = s1 + v;
            s2 = s2 + v*v;
        }
    };

    if (s0 > 1) {
        sqrt((s2 - s1*s1/s0) / (s0-1));
    }
    else {
        NaN;
    };
    """

    expression_file = create_empty_temp_file(suffix='.txt')
    open(expression_file, 'w').write(expression_string)

    np.random.seed(0) # FIXME fix the seed for tests?

    shape = (1, 1)

    nr_files = 1000

    filenames = [create_empty_temp_file() for _ in range(nr_files)]
    data      = [np.random.random(shape)  for _ in range(nr_files)]

    for i in range(nr_files): write_minc(filenames[i], data[i])

    output_file = create_empty_temp_file()

    a = minc.Calc(input_files=filenames, output_file=output_file,
                      expfile=expression_file)
    a.run()

    data_minccalc_output = read_minc(output_file)

    data_std = np.std(np.array(data), axis=0)

    for f in filenames: remove(f)
    remove(output_file)
    remove(expression_file)

    yield assert_true, np.abs(data_std - data_minccalc_output) < 0.001

@skipif(no_minc)
def test_mincblur_basic():
    """
    Blur a file.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = create_empty_temp_file()

    shape = (40, 30, 50)

    data1 = np.random.random(shape)

    write_minc(file1, data1)

    b = minc.Blur(input_file=file1, fwhm=1)

    r = b.run()

    outputs = r.outputs.get()

    output_blur = read_minc(outputs['output_file'])

    for f in outputs.values():
        if isdefined(f): remove(f)

    # FIXME Fairly crude test - check the average of the output.
    # Will this cause problems on different architectures/compilers/etc?
    yield assert_true, np.average(output_blur) - 0.19848914388   < 1e-10
    yield assert_true, np.std(output_blur)     - 0.173123419395  < 1e-10

@skipif(no_minc)
def test_mincblur_partial():
    """
    Blur a file and calculate partial derivatives.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = create_empty_temp_file()

    shape = (40, 30, 50)

    data1 = np.random.random(shape)

    write_minc(file1, data1)

    b = minc.Blur(input_file=file1, fwhm=1, partial=True)

    r = b.run()

    outputs = r.outputs.get()

    output_blur = read_minc(outputs['output_file'])
    output_dxyz = read_minc(outputs['partial_dxyz'])
    output_dx   = read_minc(outputs['partial_dx'])
    output_dy   = read_minc(outputs['partial_dy'])
    output_dz   = read_minc(outputs['partial_dz'])

    for f in outputs.values():
        if isdefined(f): remove(f)

    # FIXME Fairly crude test - check the average of the output.
    # Will this cause problems on different architectures/compilers/etc?
    yield assert_true, np.average(output_blur) -  0.19848914388     < 1e-10
    yield assert_true, np.std(output_blur)     -  0.173123419395    < 1e-10
    yield assert_true, np.average(output_dxyz) -  0.211823779297    < 1e-10
    yield assert_true, np.average(output_dx)   - -1.11356467009e-05 < 1e-10
    yield assert_true, np.average(output_dy)   - -1.33981078863e-05 < 1e-10
    yield assert_true, np.average(output_dz)   - -8.20965766907e-06 < 1e-10

@skipif(no_minc)
def test_mincmath_add():
    """
    Test mincmath, add two files together.
    """

    np.random.seed(0) # FIXME fix the seed for tests?

    file1 = create_empty_temp_file()
    file2 = create_empty_temp_file()
    output_file = create_empty_temp_file()

    shape = (40, 30)

    data1 = np.random.random(shape)
    data2 = np.random.random(shape)

    write_minc(file1, data1)
    write_minc(file2, data2)

    a = minc.Math(input_files=[file1, file2], output_file=output_file, calc_add=True)

    a.run()

    data_mincmath_output = read_minc(output_file)

    remove(file1)
    remove(file2)
    remove(output_file)

    yield assert_true, np.max(np.abs((data1 + data2) - data_mincmath_output)) < 1e-07
