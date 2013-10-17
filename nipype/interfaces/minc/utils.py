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

# FIXME double-check behaviour of usedefault=True on ranges with values specified, and also int traits. Don't want
# command line options to be specified if the user doesn't set the value.

from nipype.interfaces.base import (
    TraitedSpec,
    CommandLineInputSpec,
    CommandLine,
    StdOutCommandLineInputSpec,
    StdOutCommandLine,
    File,
    InputMultiPath,
    traits,
    isdefined,
)

import os

import warnings
warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)

class ExtractInputSpec(StdOutCommandLineInputSpec):
    """
    For the MINC command mincextract.
    """

    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    position=-1)

    _xor_write = ('write_ascii', 'write_ascii', 'write_byte',
                  'write_short', 'write_int', 'write_long',
                  'write_float', 'write_double', 'write_signed',
                  'write_unsigned',)

    write_ascii = traits.Bool(
                desc='Write out data as ascii strings (default)',
                argstr='-ascii',
                xor=_xor_write)

    write_byte = traits.Bool(
                desc='Write out data as bytes',
                argstr='-byte',
                xor=_xor_write)

    write_short = traits.Bool(
                desc='Write out data as short integers',
                argstr='-short',
                xor=_xor_write)

    write_int = traits.Bool(
                desc='Write out data as 32-bit integers',
                argstr='-int',
                xor=_xor_write)

    write_long = traits.Bool(
                desc='Superseded by -int',
                argstr='-long',
                xor=_xor_write)

    write_float = traits.Bool(
                desc='Write out data as single precision floating-point values',
                argstr='-float',
                xor=_xor_write)

    write_double = traits.Bool(
                desc='Write out data as double precision floating-point values',
                argstr='-double',
                xor=_xor_write)

    _xor_signed = ('write_signed', 'write_unsigned')

    write_signed = traits.Bool(
                desc='Write out signed data',
                argstr='-signed',
                xor=_xor_signed)

    write_unsigned = traits.Bool(
                desc='Write out unsigned data',
                argstr='-unsigned',
                xor=_xor_signed)

    write_range = traits.Tuple(
                traits.Float, traits.Float, argstr='-range %s %s',
                desc='Specify the range of output values\nDefault value: 1.79769e+308 1.79769e+308',)

    _xor_normalize = ('normalize', 'nonormalize',)

    normalize = traits.Bool(
                    desc='Normalize integer pixel values to file max and min',
                    argstr='-normalize',
                    xor=_xor_normalize)

    nonormalize = traits.Bool(
                    desc='Turn off pixel normalization',
                    argstr='-nonormalize',
                    xor=_xor_normalize)

    image_range = traits.Tuple(
                    traits.Float, traits.Float,
                    desc='Specify the range of real image values for normalization',
                    argstr='-image_range %s %s')

    image_minimum = traits.Float(
                    desc='Specify the minimum real image value for normalization. Default value: 1.79769e+308',
                    argstr='-image_minimum %s')

    image_maximum = traits.Float(
                    desc='Specify the maximum real image value for normalization. Default value: 1.79769e+308',
                    argstr='-image_maximum %s')

    start = InputMultiPath(
                    traits.Int,
                    desc='Specifies corner of hyperslab (C conventions for indices)',
                    sep=',',
                    argstr='-start %s',)

    count = InputMultiPath(
                    traits.Int,
                    desc='Specifies edge lengths of hyperslab to read',
                    sep=',',
                    argstr='-count %s',)

    # FIXME make sure that len(start) == len(count)?

    _xor_flip = ('flip_positive_direction', 'flip_negative_direction', 'flip_any_direction')

    flip_positive_direction = traits.Bool(desc='Flip images to always have positive direction.', argstr='-positive_direction', xor=_xor_flip)
    flip_negative_direction = traits.Bool(desc='Flip images to always have negative direction.', argstr='-negative_direction', xor=_xor_flip)
    flip_any_direction      = traits.Bool(desc='Do not flip images (Default).',                  argstr='-any_direction',      xor=_xor_flip)

    _xor_x_flip = ('flip_x_positive', 'flip_x_negative', 'flip_x_any')

    flip_x_positive = traits.Bool(desc='Flip images to give positive xspace:step value (left-to-right).', argstr='+xdirection',    xor=_xor_x_flip)
    flip_x_negative = traits.Bool(desc='Flip images to give negative xspace:step value (right-to-left).', argstr='-xdirection',    xor=_xor_x_flip)
    flip_x_any      = traits.Bool(desc='Don\'t flip images along x-axis (default).',                      argstr='-xanydirection', xor=_xor_x_flip)

    _xor_y_flip = ('flip_y_positive', 'flip_y_negative', 'flip_y_any')

    flip_y_positive = traits.Bool(desc='Flip images to give positive yspace:step value (post-to-ant).',   argstr='+ydirection',    xor=_xor_y_flip)
    flip_y_negative = traits.Bool(desc='Flip images to give negative yspace:step value (ant-to-post).',   argstr='-ydirection',    xor=_xor_y_flip)
    flip_y_any      = traits.Bool(desc='Don\'t flip images along y-axis (default).',                      argstr='-yanydirection', xor=_xor_y_flip)

    _xor_z_flip = ('flip_z_positive', 'flip_z_negative', 'flip_z_any')

    flip_z_positive = traits.Bool(desc='Flip images to give positive zspace:step value (inf-to-sup).',  argstr='+zdirection',       xor=_xor_z_flip)
    flip_z_negative = traits.Bool(desc='Flip images to give negative zspace:step value (sup-to-inf).',  argstr='-zdirection',       xor=_xor_z_flip)
    flip_z_any      = traits.Bool(desc='Don\'t flip images along z-axis (default).',                    argstr='-zanydirection',    xor=_xor_z_flip)

class ExtractOutputSpec(TraitedSpec):
    # FIXME Not sure if I'm defining the outout specs correctly.

    output_file = File(
                    desc='output file',
                    exists=True,
                    genfile=True,)

class ExtractTask(StdOutCommandLine):
    input_spec  = ExtractInputSpec
    output_spec = ExtractOutputSpec
    cmd = 'mincextract'

    def _gen_outfilename(self):
        """
        Extract foo.mnc to foo.raw.
        """

        output_file = self.inputs.output_file

        if isdefined(output_file):
            return output_file
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '.raw'

class ToRawInputSpec(StdOutCommandLineInputSpec):
    """
    For the MINC command minctoraw.
    """

    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    # FIXME xor on the write_ options

    write_byte = traits.Bool(
                desc='Write out data as bytes',
                argstr='-byte',)

    write_short = traits.Bool(
                desc='Write out data as short integers',
                argstr='-short',)

    write_int = traits.Bool(
                desc='Write out data as 32-bit integers',
                argstr='-int',)

    write_long = traits.Bool(
                desc='Superseded by -int',
                argstr='-long',)

    write_float = traits.Bool(
                desc='Write out data as single precision floating-point values',
                argstr='-float',)

    write_double = traits.Bool(
                desc='Write out data as double precision floating-point values',
                argstr='-double',)

    # FIXME xor option on signed/unsigned?

    write_signed = traits.Bool(
                desc='Write out signed data',
                argstr='-signed',)

    write_unsigned = traits.Bool(
                desc='Write out unsigned data',
                argstr='-unsigned',)

    write_range = traits.Tuple(
                traits.Float, traits.Float, argstr='-range %s %s',
                desc='Specify the range of output values\nDefault value: 1.79769e+308 1.79769e+308',) # FIXME minctoraw output is missing a negative?

    _xor_normalize = ('normalize', 'nonormalize',)

    normalize = traits.Bool(
                    desc='Normalize integer pixel values to file max and min',
                    argstr='-normalize',
                    xor=_xor_normalize)

    nonormalize = traits.Bool(
                    desc='Turn off pixel normalization',
                    argstr='-nonormalize',
                    xor=_xor_normalize)

class ToRawOutputSpec(TraitedSpec):
    # FIXME Not sure if I'm defining the outout specs correctly.

    output_file = File(
                    desc='output file',
                    exists=True,
                    genfile=True,)

class ToRawTask(StdOutCommandLine):
    input_spec  = ToRawInputSpec
    output_spec = ToRawOutputSpec
    cmd = 'minctoraw'

    def _gen_outfilename(self):
        # FIXME check isdefined as in the Extract code.
        """
        Convert foo.mnc to foo.raw.
        """
        return os.path.splitext(self.inputs.input_file)[0] + '.raw'

class ConvertInputSpec(CommandLineInputSpec):
    input_file = File(
                    desc='input file for converting',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    mandatory=True,
                    genfile=False,
                    argstr='%s',
                    position=-1,)

    clobber = traits.Bool(
                desc='Overwrite existing file.',
                argstr='-clobber',)

    two = traits.Bool(
                desc='Create a MINC 2 output file.',
                argstr='-2',)

    template = traits.Bool(
                desc='Create a template file.',
                argstr='-template',)

    compression = traits.Enum(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                            argstr='-compress %s',
                            desc='Set the compression level, from 0 (disabled) to 9 (maximum).',)

    chunk = traits.Range(low=0,
                        desc='Set the target block size for chunking (0 default, >1 block size).',
                        value=0,
                        usedefault=False,
                        argstr='-chunk %d',)

class ConvertOutputSpec(TraitedSpec):
    # FIXME Am I defining the output spec correctly?
    output_file = File(
                    desc='output file',
                    exists=True,)

class ConvertTask(CommandLine):
    input_spec  = ConvertInputSpec
    output_spec = ConvertOutputSpec
    cmd = 'mincconvert'

    def _list_outputs(self):
        # FIXME seems generic, is this necessary?
        outputs = self.output_spec().get()
        outputs['output_file'] = self.inputs.output_file
        return outputs

class CopyInputSpec(CommandLineInputSpec):
    """
    Implement minccopy? Its man page says:

        NOTE: This program is intended primarily for use with scripts such
        as mincedit.  It does not follow the typical design rules of most MINC
        command-line tools and therefore should be  used only with caution.

    It doesn't run on a standard MNC file, e.g.

        $ minccopy /home/carlo/tmp/foo.mnc /tmp/foo_copy.mnc
        (from miopen): Can't write compressed file
        ncvarid: ncid -1: NetCDF: Not a valid ID

    """

    input_file = File(
                    desc='input file to copy',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    mandatory=True,
                    genfile=False,
                    argstr='%s',
                    position=-1,)

    pixel_values = traits.Bool(
                desc='Copy pixel values as is.',
                argstr='-pixel_values',)

    real_values = traits.Bool(
                desc='Copy real pixel intensities (default).',
                argstr='-real_values',
                usedefault=False,)

class CopyOutputSpec(TraitedSpec):
    # FIXME Am I defining the output spec correctly?
    output_file = File(
                    desc='output file',
                    exists=True,)

class CopyTask(CommandLine):
    input_spec  = CopyInputSpec
    output_spec = CopyOutputSpec
    cmd = 'minccopy'

    def _list_outputs(self):
        # FIXME seems generic, is this necessary?
        outputs = self.output_spec().get()
        outputs['output_file'] = self.inputs.output_file
        return outputs

class ToEcatInputSpec(CommandLineInputSpec):
    input_file = File(
                    desc='input file to convert',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    mandatory=False,
                    genfile=True,
                    argstr='%s',
                    position=-1,)

    ignore_patient_variable = traits.Bool(
                    desc='Ignore informations from the minc patient variable.',
                    argstr='-ignore_patient_variable',)

    ignore_study_variable = traits.Bool(
                    desc='Ignore informations from the minc study variable.',
                    argstr='-ignore_study_variable',)

    ignore_acquisition_variable = traits.Bool(
                    desc='Ignore informations from the minc acquisition variable.',
                    argstr='-ignore_acquisition_variable',)

    ignore_ecat_acquisition_variable = traits.Bool(
                    desc='Ignore informations from the minc ecat_acquisition variable.',
                    argstr='-ignore_ecat_acquisition_variable',)

    ignore_ecat_main = traits.Bool(
                    desc='Ignore informations from the minc ecat-main variable.',
                    argstr='-ignore_ecat_main',)

    ignore_ecat_subheader_variable = traits.Bool(
                    desc='Ignore informations from the minc ecat-subhdr variable.',
                    argstr='-ignore_ecat_subheader_variable',)

    no_decay_corr_fctr = traits.Bool(
                    desc='Do not compute the decay correction factors',
                    argstr='-no_decay_corr_fctr',)

    voxels_as_integers = traits.Bool(
                    desc='Voxel values are treated as integers, scale and calibration factors are set to unity',
                    argstr='-label',)

class ToEcatOutputSpec(TraitedSpec):
    # FIXME Am I defining the output spec correctly?
    output_file = File(
                    desc='output file',
                    exists=True,)

class ToEcatTask(CommandLine):
    input_spec  = ToEcatInputSpec
    output_spec = ToEcatOutputSpec
    cmd = 'minctoecat'

    def _list_outputs(self):
        # FIXME seems generic, is this necessary?
        outputs = self.output_spec().get()
        outputs['output_file'] = self.inputs.output_file
        return outputs

    def _gen_filename(self, name):
        if name == 'output_file':
            return os.path.splitext(self.inputs.input_file)[0] + '.v'
        return None

class DumpInputSpec(StdOutCommandLineInputSpec):
    """
    For the MINC command mincdump.
    """

    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    _xor_coords_or_header = ('coordinate_data', 'header_data',)

    coordinate_data = traits.Bool(
                    desc='Coordinate variable data and header information',
                    argstr='-c',
                    xor=_xor_coords_or_header,)

    header_data = traits.Bool(
                    desc='Header information only, no data',
                    argstr='-h',
                    xor=_xor_coords_or_header,)

    _xor_annotations = ('annotations_brief', 'annotations_full',)

    # FIXME Instead of an enum, make a separate Bool trait called fortran_indices and another
    # called c_indices?

    annotations_brief = traits.Enum('c', 'f',
                            argstr='-b %s',
                            desc='Brief annotations for C or Fortran indices in data',
                            xor=_xor_annotations)

    annotations_full = traits.Enum('c', 'f',
                            argstr='-f %s',
                            desc='Full annotations for C or Fortran indices in data',
                            xor=_xor_annotations)

    variables = InputMultiPath(
                            traits.Str,
                            desc='Output data for specified variables only',
                            sep=',',
                            argstr='-v %s',)

    line_length = traits.Range(low=0,
                        desc='Line length maximum in data section (default 80)',
                        value=0,
                        usedefault=False,
                        argstr='-l %d',)

    netcdf_name = traits.Str(
                        desc='Name for netCDF (default derived from file name)',
                        argstr='-n %s',)

    precision = traits.Either(
                        traits.Int(),
                        traits.Tuple(traits.Int, traits.Int),
                        desc='Display floating-point values with less precision',
                        argstr='%s',)

class DumpOutputSpec(TraitedSpec):
    # FIXME Not sure if I'm defining the outout specs correctly.

    output_file = File(
                    desc='output file',
                    exists=True,
                    genfile=True,)

class DumpTask(StdOutCommandLine):
    """
    >>> from os.path import splitext
    >>> from nipype.testing import example_data
    >>> input_file = example_data('minc/minc_test_2D_01.mnc')
    >>> output_file = splitext(input_file)[0] + '.txt'
    >>> toraw = DumpTask(input_file=input_file)
    >>> toraw.cmdline == 'mincdump %s > %s' % (input_file, output_file,)
    True
    """

    input_spec  = DumpInputSpec
    output_spec = DumpOutputSpec
    cmd = 'mincdump'

    def _format_arg(self, name, spec, value):
        if name == 'precision':
            if isinstance(value, int):
                return '-p %d' % value
            elif isinstance(value, tuple):
                return '-p %d,%d' % (value[0], value[1],)
            else:
                raise NotImplemented # FIXME some other exception?
        return super(DumpTask, self)._format_arg(name, spec, value)

    # FIXME Are we forced to send outout to a file? Can we pipe it
    # to another minc command directly?
    def _gen_outfilename(self):
        # FIXME allow user to specify an output file as in the Extract task.
        """
        Dump foo.mnc to foo.txt.
        """
        return os.path.splitext(self.inputs.input_file)[0] + '.txt'

class AverageInputSpec(CommandLineInputSpec):
    _xor_input_files = ('input_files', 'filelist')

    input_files = InputMultiPath(
                    traits.File,
                    desc='input file(s) for averaging',
                    exists=True,
                    mandatory=True,
                    sep=' ', # FIXME test with files that contain spaces - does InputMultiPath do the right thing?
                    argstr='%s',
                    position=-2, # FIXME test with multiple files, is order ok?
                    xor=_xor_input_files) # FIXME test this xor

    output_file = File(
                    desc='output file',
                    mandatory=True,
                    genfile=False,
                    argstr='%s',
                    position=-1,)

    two = traits.Bool(desc='Produce a MINC 2.0 format output file', argstr='-2')

    _xor_clobber = ('clobber', 'no_clobber')

    clobber     = traits.Bool(desc='Overwrite existing file.',                  argstr='-clobber',      xor=_xor_clobber)
    no_clobber  = traits.Bool(desc='Don\'t overwrite existing file (default).', argstr='-noclobber',    xor=_xor_clobber)

    _xor_verbose = ('verbose', 'quiet',)

    verbose = traits.Bool(desc='Print out log messages (default).', argstr='-verbose',  xor=_xor_verbose)
    quiet   = traits.Bool(desc='Do not print out log messages.',    argstr='-quiet',    xor=_xor_verbose)

    debug   = traits.Bool(desc='Print out debugging messages.', argstr='-debug')

    # FIXME How to handle stdin option ('-') here? Not relevant?
    filelist = traits.File(desc='Specify the name of a file containing input file names.', argstr='-filelist %s', mandatory=True, xor=_xor_input_files)

    _xor_check_dimensions = ('check_dimensions', 'no_check_dimensions',)

    check_dimensions    = traits.Bool(desc='Check that dimension info matches across files (default).', argstr='-check_dimensions',     xor=_xor_check_dimensions)
    no_check_dimensions = traits.Bool(desc='Do not check dimension info.',                              argstr='-nocheck_dimensions',   xor=_xor_check_dimensions)


    # FIXME mincaverage seems to accept more than one of these options; I assume
    # that it takes the last one, and it makes more sense for these to be
    # put into an xor case.

    _xor_format = ('format_filetype', 'format_byte', 'format_short',
                   'format_int', 'format_long', 'format_float', 'format_double',
                   'format_signed', 'format_unsigned',)

    format_filetype     = traits.Bool(desc='Use data type of first file (default).',                    argstr='-filetype', xor=_xor_format)
    format_byte         = traits.Bool(desc='Write out byte data.',                                      argstr='-byte',     xor=_xor_format)
    format_short        = traits.Bool(desc='Write out short integer data.',                             argstr='-short',    xor=_xor_format)
    format_int          = traits.Bool(desc='Write out 32-bit integer data.',                            argstr='-int',      xor=_xor_format)
    format_long         = traits.Bool(desc='Superseded by -int.',                                       argstr='-long',     xor=_xor_format)
    format_float        = traits.Bool(desc='Write out single-precision floating-point data.',           argstr='-float',    xor=_xor_format)
    format_double       = traits.Bool(desc='Write out double-precision floating-point data.',           argstr='-double',   xor=_xor_format)
    format_signed       = traits.Bool(desc='Write signed integer data.',                                argstr='-signed',   xor=_xor_format)
    format_unsigned     = traits.Bool(desc='Write unsigned integer data (default if type specified).',  argstr='-unsigned', xor=_xor_format) # FIXME mark with default=?


    max_buffer_size_in_kb = traits.Range(
                                low=0,
                                desc='Specify the maximum size of the internal buffers (in kbytes).',
                                value=4096,
                                usedefault=False,
                                argstr='-max_buffer_size_in_kb %d',)

    _xor_normalize = ('normalize', 'nonormalize',)
    normalize   = traits.Bool(desc='Normalize data sets for mean intensity.', argstr='-normalize', xor=_xor_normalize)
    nonormalize = traits.Bool(desc='Do not normalize data sets (default).',   argstr='-nonormalize', xor=_xor_normalize, default=True) # FIXME check default=? behaviour

    voxel_range = traits.Tuple(
                traits.Int, traits.Int, argstr='-range %d %d',
                desc='Valid range for output data.',) # FIXME ensure min <= max??? Ditto for other range parameters.

    sdfile = traits.File(
                desc='Specify an output sd file (default=none).',
                argstr='-sdfile %s',)

    _xor_copy_header = ('copy_header, no_copy_header')

    copy_header     = traits.Bool(desc='Copy all of the header from the first file (default for one file).',            argstr='-copy_header',   xor=_xor_copy_header)
    no_copy_header  = traits.Bool(desc='Do not copy all of the header from the first file (default for many files)).',  argstr='-nocopy_header', xor=_xor_copy_header)

    avgdim = traits.Str(desc='Specify a dimension along which we wish to average.', argstr='-avgdim %s')

    binarize = traits.Bool(desc='Binarize the volume by looking for values in a given range.', argstr='-binarize')

    binrange = traits.Tuple(
                traits.Float, traits.Float, argstr='-binrange %s %s',
                desc='Specify a range for binarization. Default value: 1.79769e+308 -1.79769e+308.') # FIXME shouldn't that be -1.79769e+308 1.79769e+308? Min then max?

    binvalue = traits.Float(desc='Specify a target value (+/- 0.5) for binarization. Default value: -1.79769e+308', argstr='-binvalue %s')
		
    weights = InputMultiPath(
                            traits.Str,
                            desc='Specify weights for averaging ("<w1>,<w2>,...").',
                            sep=',',
                            argstr='-weights %s',)

    width_weighted = traits.Bool(desc='Weight by dimension widths when -avgdim is used.', argstr='-width_weighted', requires=('avgdim',))

class AverageOutputSpec(TraitedSpec):
    # FIXME Am I defining the output spec correctly?
    output_file = File(
                    desc='output file',
                    exists=True,)

class AverageTask(CommandLine):
    input_spec  = AverageInputSpec
    output_spec = AverageOutputSpec
    cmd = 'mincaverage'

    def _list_outputs(self):
        # FIXME seems generic, is this necessary?
        outputs = self.output_spec().get()
        outputs['output_file'] = self.inputs.output_file
        return outputs

class BlobInputSpec(CommandLineInputSpec):
    input_file = File(
                    desc='input file to blob',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    mandatory=True,
                    genfile=False,
                    argstr='%s',
                    position=-1,)

    trace       = traits.Bool(desc='compute the trace (approximate growth and shrinkage) -- FAST',  argstr='-trace')
    determinant = traits.Bool(desc='compute the determinant (exact growth and shrinkage) -- SLOW',  argstr='-determinant')
    translation = traits.Bool(desc='compute translation (structure displacement)',                  argstr='-translation')
    magnitude   = traits.Bool(desc='compute the magnitude of the displacement vector',              argstr='-magnitude')

class BlobOutputSpec(TraitedSpec):
    # FIXME Am I defining the output spec correctly?
    output_file = File(
                    desc='output file',
                    exists=True,)

class BlobTask(CommandLine):
    input_spec  = BlobInputSpec
    output_spec = BlobOutputSpec
    cmd = 'mincblob'

    def _list_outputs(self):
        # FIXME seems generic, is this necessary?
        outputs = self.output_spec().get()
        outputs['output_file'] = self.inputs.output_file
        return outputs

class CalcInputSpec(CommandLineInputSpec):
    _xor_input_files = ('input_files', 'filelist')

    input_files = InputMultiPath(
                    traits.File,
                    desc='input file(s) for calculation',
                    exists=True,
                    mandatory=True,
                    sep=' ', # FIXME test with files that contain spaces - does InputMultiPath do the right thing?
                    argstr='%s',
                    position=-2,) # FIXME test with multiple files, is order ok?


    two = traits.Bool(desc='Produce a MINC 2.0 format output file', argstr='-2')

    _xor_clobber = ('clobber', 'no_clobber')

    clobber     = traits.Bool(desc='Overwrite existing file.',                  argstr='-clobber',      xor=_xor_clobber)
    no_clobber  = traits.Bool(desc='Don\'t overwrite existing file (default).', argstr='-noclobber',    xor=_xor_clobber)

    _xor_verbose = ('verbose', 'quiet',)

    verbose = traits.Bool(desc='Print out log messages (default).', argstr='-verbose',  xor=_xor_verbose)
    quiet   = traits.Bool(desc='Do not print out log messages.',    argstr='-quiet',    xor=_xor_verbose)

    debug   = traits.Bool(desc='Print out debugging messages.', argstr='-debug')

    # FIXME How to handle stdin option here? Not relevant?
    filelist = traits.File(desc='Specify the name of a file containing input file names (- for stdin).', argstr='-filelist %s', mandatory=True, xor=_xor_input_files)

    _xor_copy_header = ('copy_header, no_copy_header')

    copy_header     = traits.Bool(desc='Copy all of the header from the first file.',            argstr='-copy_header',   xor=_xor_copy_header)
    no_copy_header  = traits.Bool(desc='Do not copy all of the header from the first file.',  argstr='-nocopy_header', xor=_xor_copy_header)

    # FIXME mincaverage seems to accept more than one of these options; I assume
    # that it takes the last one, and it makes more sense for these to be
    # put into an xor case.

    _xor_format = ('format_filetype', 'format_byte', 'format_short',
                   'format_int', 'format_long', 'format_float', 'format_double',
                   'format_signed', 'format_unsigned',)

    format_filetype     = traits.Bool(desc='Use data type of first file (default).',                    argstr='-filetype', xor=_xor_format)
    format_byte         = traits.Bool(desc='Write out byte data.',                                      argstr='-byte',     xor=_xor_format)
    format_short        = traits.Bool(desc='Write out short integer data.',                             argstr='-short',    xor=_xor_format)
    format_int          = traits.Bool(desc='Write out 32-bit integer data.',                            argstr='-int',      xor=_xor_format)
    format_long         = traits.Bool(desc='Superseded by -int.',                                       argstr='-long',     xor=_xor_format)
    format_float        = traits.Bool(desc='Write out single-precision floating-point data.',           argstr='-float',    xor=_xor_format)
    format_double       = traits.Bool(desc='Write out double-precision floating-point data.',           argstr='-double',   xor=_xor_format)
    format_signed       = traits.Bool(desc='Write signed integer data.',                                argstr='-signed',   xor=_xor_format)
    format_unsigned     = traits.Bool(desc='Write unsigned integer data (default if type specified).',  argstr='-unsigned', xor=_xor_format) # FIXME mark with default=?

    voxel_range = traits.Tuple(
                traits.Int, traits.Int, argstr='-range %d %d',
                desc='Valid range for output data.',)

    max_buffer_size_in_kb = traits.Range(
                                low = 0,
                                desc='Specify the maximum size of the internal buffers (in kbytes).',
                                value=0,
                                usedefault=False,
                                argstr='-max_buffer_size_in_kb %d',)

    _xor_check_dimensions = ('check_dimensions', 'no_check_dimensions',)

    check_dimensions    = traits.Bool(desc='Check that files have matching dimensions (default).', argstr='-check_dimensions',     xor=_xor_check_dimensions)
    no_check_dimensions = traits.Bool(desc='Do not check that files have matching dimensions.',                              argstr='-nocheck_dimensions',   xor=_xor_check_dimensions)

    # FIXME are ignore_nan and propagate_nan mutually exclusive?
    ignore_nan = traits.Bool(desc='Ignore invalid data (NaN) for accumulations.', argstr='-ignore_nan')

    propagate_nan = traits.Bool(desc='Invalid data in any file at a voxel produces a NaN (default).', argstr='-propagate_nan')

    # FIXME Double-check that these are mutually exclusive?
    _xor_nan_zero_illegal = ('output_nan', 'output_zero', 'output_illegal_value')

    output_nan      = traits.Bool(desc='Output NaN when an illegal operation is done (default).',                           argstr='-nan',              xor=_xor_nan_zero_illegal)
    output_zero     = traits.Bool(desc='Output zero when an illegal operation is done.',                                    argstr='-zero',             xor=_xor_nan_zero_illegal)
    output_illegal  = traits.Bool(desc='Value to write out when an illegal operation is done. Default value: 1.79769e+308', argstr='-illegal_value',    xor=_xor_nan_zero_illegal)

    _xor_expression = ('expression', 'expfile')

    expression = traits.Str(desc='Expression to use in calculations.',    argstr='-expression %s', xor=_xor_expression, mandatory=True)
    expfile    = traits.File(desc='Name of file containing expression.',  argstr='-expfile %s',    xor=_xor_expression, mandatory=True)

    # FIXME test this one, the argstr will probably need tweaking, see _format_arg.
    outfiles = traits.List(
                traits.Tuple(traits.Str, traits.File, argstr='-outfile %s %s',
                desc='List of (symbol, file) tuples indicating that output should be written to the specified file, taking values from the symbol which should be created in the expression (see the EXAMPLES section). If this option is given, then all non-option arguments are taken as input files. This option can be used multiple times for multiple output files.'))

    eval_width = traits.Int(200, desc='Number of voxels to evaluate simultaneously.', argstr='-eval_width %s', usedefault=False)


class CalcOutputSpec(TraitedSpec):
    # FIXME Am I defining the output spec correctly?
    output_file = File(
                    desc='output file',
                    exists=True,)

class CalcTask(CommandLine):
    input_spec  = CalcInputSpec
    output_spec = CalcOutputSpec
    cmd = 'minccalc'

    def _list_outputs(self):
        # FIXME seems generic, is this necessary?
        outputs = self.output_spec().get()
        outputs['output_file'] = self.inputs.output_file
        return outputs

