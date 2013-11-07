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

# FIXME -filelist options accept '-' for stdin. Can we support this?

# FIXME Can we check the arguments to "-range min max" to avoid min > max?

# FIXME output_file(s) should be optional, in line with the Nipype convention.

# FIXME check all interface .help() outputs for out_file vs output_file.

# FIXME check that the _list_outputs functions do a return at the end.

# FIXME factor out common definitions, e.g. max_buffer_size_in_kb?

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
import os.path

import warnings
warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)

def aggregate_filename(files, new_suffix):
    """
    Try to work out a sensible name given a set of files that have
    been combined in some way (e.g. averaged). If we can't work out a
    sensible prefix, we use the first filename in the list.

    Examples
    --------

    >>> import nipype.interfaces.minc as minc
    >>> minc.utils.aggregate_filename(['/tmp/foo1.mnc', '/tmp/foo2.mnc', '/tmp/foo3.mnc'], 'averaged')
    '/tmp/foo_averaged.mnc'

    >>> minc.utils.aggregate_filename(['/tmp/foo1.mnc', '/tmp/blah1.mnc'], 'averaged')
    '/tmp/foo1_averaged.mnc'

    """

    path            = os.path.split(files[0])[0]
    names           = [os.path.splitext(os.path.split(x)[1])[0] for x in files]
    common_prefix   = os.path.commonprefix(names)

    if common_prefix == '':
        return os.path.abspath(os.path.join(path, os.path.splitext(files[0])[0] + '_' + new_suffix + '.mnc'))
    else:
        return os.path.abspath(os.path.join(path, common_prefix + '_' + new_suffix + '.mnc'))


class ExtractInputSpec(StdOutCommandLineInputSpec):
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
                desc='Write out data as ascii strings (default).',
                argstr='-ascii',
                xor=_xor_write)

    write_byte = traits.Bool(
                desc='Write out data as bytes.',
                argstr='-byte',
                xor=_xor_write)

    write_short = traits.Bool(
                desc='Write out data as short integers.',
                argstr='-short',
                xor=_xor_write)

    write_int = traits.Bool(
                desc='Write out data as 32-bit integers.',
                argstr='-int',
                xor=_xor_write)

    write_long = traits.Bool(
                desc='Superseded by write_int.',
                argstr='-long',
                xor=_xor_write)

    write_float = traits.Bool(
                desc='Write out data as single precision floating-point values.',
                argstr='-float',
                xor=_xor_write)

    write_double = traits.Bool(
                desc='Write out data as double precision floating-point values.',
                argstr='-double',
                xor=_xor_write)

    _xor_signed = ('write_signed', 'write_unsigned')

    write_signed = traits.Bool(
                desc='Write out signed data.',
                argstr='-signed',
                xor=_xor_signed)

    write_unsigned = traits.Bool(
                desc='Write out unsigned data.',
                argstr='-unsigned',
                xor=_xor_signed)

    write_range = traits.Tuple(
                traits.Float, traits.Float, argstr='-range %s %s',
                desc='Specify the range of output values\nDefault value: 1.79769e+308 1.79769e+308.',)

    _xor_normalize = ('normalize', 'nonormalize',)

    normalize = traits.Bool(
                    desc='Normalize integer pixel values to file max and min.',
                    argstr='-normalize',
                    xor=_xor_normalize)

    nonormalize = traits.Bool(
                    desc='Turn off pixel normalization.',
                    argstr='-nonormalize',
                    xor=_xor_normalize)

    image_range = traits.Tuple(
                    traits.Float, traits.Float,
                    desc='Specify the range of real image values for normalization.',
                    argstr='-image_range %s %s')

    image_minimum = traits.Float(
                    desc='Specify the minimum real image value for normalization. Default value: 1.79769e+308.',
                    argstr='-image_minimum %s')

    image_maximum = traits.Float(
                    desc='Specify the maximum real image value for normalization. Default value: 1.79769e+308.',
                    argstr='-image_maximum %s')

    start = InputMultiPath(
                    traits.Int,
                    desc='Specifies corner of hyperslab (C conventions for indices).',
                    sep=',',
                    argstr='-start %s',)

    count = InputMultiPath(
                    traits.Int,
                    desc='Specifies edge lengths of hyperslab to read.',
                    sep=',',
                    argstr='-count %s',)

    # FIXME Can we make sure that len(start) == len(count)?

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
    output_file = File(desc='output file in raw/text format', exists=True)

class Extract(StdOutCommandLine):
    """Dump a hyperslab of MINC file data.

    Examples
    --------

    >>> from nipype.interfaces.minc import Extract
    >>> from nipype.testing import minc2Dfile

    >>> extract = Extract(input_file=minc2Dfile)
    >>> extract.run() # doctest: +SKIP

    >>> extract = Extract(input_file=minc2Dfile, start=[3, 10, 5], count=[4, 4, 4]) # extract a 4x4x4 slab at offset [3, 10, 5]
    >>> extract.run() # doctest: +SKIP
    """

    input_spec  = ExtractInputSpec
    output_spec = ExtractOutputSpec
    _cmd = 'mincextract'

    # FIXME Does this play nicely with a workflow?
    def _gen_outfilename(self):
        """
        If the user specified output_file then return that, otherwise
        return the full path to the input file with the extension
        changed to '.raw'.
        """

        output_file = self.inputs.output_file

        if isdefined(output_file):
            return output_file
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '.raw'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

class ToRawInputSpec(StdOutCommandLineInputSpec):
    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    position=-1)

    _xor_write = ('write_byte', 'write_short', 'write_int',
                  'write_long', 'write_float', 'write_double')

    write_byte = traits.Bool(
                desc='Write out data as bytes.',
                argstr='-byte',
                xor=_xor_write)

    write_short = traits.Bool(
                desc='Write out data as short integers.',
                argstr='-short',
                xor=_xor_write)

    write_int = traits.Bool(
                desc='Write out data as 32-bit integers.',
                argstr='-int',
                xor=_xor_write)

    write_long = traits.Bool(
                desc='Superseded by write_int.',
                argstr='-long',
                xor=_xor_write)

    write_float = traits.Bool(
                desc='Write out data as single precision floating-point values.',
                argstr='-float',
                xor=_xor_write)

    write_double = traits.Bool(
                desc='Write out data as double precision floating-point values.',
                argstr='-double',
                xor=_xor_write)

    _xor_signed = ('write_signed', 'write_unsigned')

    write_signed = traits.Bool(
                desc='Write out signed data.',
                argstr='-signed',
                xor=_xor_signed)

    write_unsigned = traits.Bool(
                desc='Write out unsigned data.',
                argstr='-unsigned',
                xor=_xor_signed)

    write_range = traits.Tuple(
                traits.Float, traits.Float, argstr='-range %s %s',
                desc='Specify the range of output values\nDefault value: 1.79769e+308 1.79769e+308.',)

    _xor_normalize = ('normalize', 'nonormalize',)

    normalize = traits.Bool(
                    desc='Normalize integer pixel values to file max and min.',
                    argstr='-normalize',
                    xor=_xor_normalize)

    nonormalize = traits.Bool(
                    desc='Turn off pixel normalization.',
                    argstr='-nonormalize',
                    xor=_xor_normalize)

class ToRawOutputSpec(TraitedSpec):
    output_file = File(desc='output file in raw format', exists=True)

class ToRaw(StdOutCommandLine):
    """Dump a chunk of MINC file data. This program is largely
    superceded by mincextract (see Extract).

    Examples
    --------

    >>> from nipype.interfaces.minc import ToRaw
    >>> from nipype.testing import minc2Dfile

    >>> toraw = ToRaw(input_file=minc2Dfile)
    >>> toraw.run() # doctest: +SKIP

    >>> toraw = ToRaw(input_file=minc2Dfile, write_range=(0, 100))
    >>> toraw.run() # doctest: +SKIP
    """

    input_spec  = ToRawInputSpec
    output_spec = ToRawOutputSpec
    _cmd = 'minctoraw'

    def _gen_outfilename(self):
        """
        If the user specified output_file then return that, otherwise
        return the full path to the input file with the extension
        changed to '.raw'.
        """

        output_file = self.inputs.output_file

        if isdefined(output_file):
            return output_file
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '.raw'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

class ConvertInputSpec(CommandLineInputSpec):
    input_file = File(
                    desc='input file for converting',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    genfile=True,
                    argstr='%s',
                    position=-1,)

    clobber = traits.Bool(
                desc='Overwrite existing file.',
                argstr='-clobber', usedefault=True, default_value=True)

    two = traits.Bool(
                desc='Create a MINC 2 output file.',
                argstr='-2',)

    template = traits.Bool(
                desc='Create a template file. The dimensions, variables, and attributes of the input file are preserved but all data it set to zero.',
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
    output_file = File(desc='output file', exists=True)

class Convert(CommandLine):
    """convert between MINC 1 to MINC 2 format.

    Examples
    --------

    >>> from nipype.interfaces.minc import ToRaw
    >>> from nipype.testing import minc2Dfile
    >>> c = Convert(input_file=minc2Dfile, output_file='/tmp/out.mnc', two=True) # Convert to MINC2 format.
    >>> c.run() # doctest: +SKIP
    """


    input_spec  = ConvertInputSpec
    output_spec = ConvertOutputSpec
    _cmd = 'mincconvert'

    def _gen_filename(self, name):
        if name == 'output_file':
            output_file = self.inputs.output_file

            if isdefined(output_file):
                return os.path.abspath(output_file)
            else:
                return aggregate_filename([self.inputs.input_file], 'convert_output')
        else:
            raise NotImplemented

    def _gen_outfilename(self):
        return self._gen_filename('output_file')

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

class CopyInputSpec(CommandLineInputSpec):
    input_file = File(
                    desc='input file to copy',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    genfile=True,
                    argstr='%s',
                    position=-1,)

    _xor_pixel = ('pixel_values', 'real_values')

    pixel_values = traits.Bool(
                desc='Copy pixel values as is.',
                argstr='-pixel_values',
                xor=_xor_pixel)

    real_values = traits.Bool(
                desc='Copy real pixel intensities (default).',
                argstr='-real_values',
                xor=_xor_pixel)

class CopyOutputSpec(TraitedSpec):
    output_file = File(desc='output file', exists=True)

class Copy(CommandLine):
    """
    Copy image values from one MINC file to another. Both the input
    and output files must exist, and the images in both files must
    have an equal number dimensions and equal dimension lengths.

    NOTE: This program is intended primarily for use with scripts
    such as mincedit. It does not follow the typical design rules of
    most MINC command-line tools and therefore should be used only
    with caution.
    """

    input_spec  = CopyInputSpec
    output_spec = CopyOutputSpec
    _cmd = 'minccopy'

    def _gen_filename(self, name):
        if name == 'output_file':
            output_file = self.inputs.output_file

            if isdefined(output_file):
                return os.path.abspath(output_file)
            else:
                return aggregate_filename([self.inputs.input_file], 'copy_output')
        else:
            raise NotImplemented

    def _gen_outfilename(self):
        return self._gen_filename('output_file')

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
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
    output_file = File(desc='output file', exists=True)

class ToEcat(CommandLine):
    """Convert a 2D image, a 3D volumes or a 4D dynamic volumes
    written in MINC file format to a 2D, 3D or 4D Ecat7 file.

    Examples
    --------

    >>> from nipype.interfaces.minc import ToEcat
    >>> from nipype.testing import minc2Dfile

    >>> c = ToEcat(input_file=minc2Dfile)
    >>> c.run() # doctest: +SKIP

    >>> c = ToEcat(input_file=minc2Dfile, voxels_as_integers=True)
    >>> c.run() # doctest: +SKIP

    """

    input_spec  = ToEcatInputSpec
    output_spec = ToEcatOutputSpec
    _cmd = 'minctoecat'

    def _gen_outfilename(self):
        """
        If the user specified output_file then return that, otherwise
        return the full path to the input file with the extension
        changed to '.v'.
        """

        output_file = self.inputs.output_file

        if isdefined(output_file):
            return output_file
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '.v'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

class DumpInputSpec(StdOutCommandLineInputSpec):
    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    position=-1)

    _xor_coords_or_header = ('coordinate_data', 'header_data',)

    coordinate_data = traits.Bool(
                    desc='Coordinate variable data and header information.',
                    argstr='-c',
                    xor=_xor_coords_or_header)

    header_data = traits.Bool(
                    desc='Header information only, no data.',
                    argstr='-h',
                    xor=_xor_coords_or_header)

    _xor_annotations = ('annotations_brief', 'annotations_full',)

    annotations_brief = traits.Enum('c', 'f',
                            argstr='-b %s',
                            desc='Brief annotations for C or Fortran indices in data.',
                            xor=_xor_annotations)

    annotations_full = traits.Enum('c', 'f',
                            argstr='-f %s',
                            desc='Full annotations for C or Fortran indices in data.',
                            xor=_xor_annotations)

    variables = InputMultiPath(
                            traits.Str,
                            desc='Output data for specified variables only.',
                            sep=',',
                            argstr='-v %s')

    line_length = traits.Range(low=0,
                        desc='Line length maximum in data section (default 80).',
                        value=80,
                        usedefault=False,
                        argstr='-l %d')

    netcdf_name = traits.Str(
                        desc='Name for netCDF (default derived from file name).',
                        argstr='-n %s')

    precision = traits.Either(
                        traits.Int(),
                        traits.Tuple(traits.Int, traits.Int),
                        desc='Display floating-point values with less precision',
                        argstr='%s',) # See _format_arg in Dump for actual formatting.

class DumpOutputSpec(TraitedSpec):
    output_file = File(desc='output file', exists=True)

class Dump(StdOutCommandLine):
    """Dump a MINC file. Typically used in conjunction with mincgen (see Gen).

    Examples
    --------

    >>> from nipype.interfaces.minc import Dump
    >>> from nipype.testing import minc2Dfile

    >>> dump = Dump(input_file=minc2Dfile)
    >>> dump.run() # doctest: +SKIP

    >>> dump = Dump(input_file=minc2Dfile, output_file='/tmp/out.txt', precision=(3, 4))
    >>> dump.run() # doctest: +SKIP

    """

    input_spec  = DumpInputSpec
    output_spec = DumpOutputSpec
    _cmd = 'mincdump'

    def _format_arg(self, name, spec, value):
        if name == 'precision':
            if isinstance(value, int):
                return '-p %d' % value
            elif isinstance(value, tuple) and isinstance(value[0], int) and isinstance(value[1], int):
                return '-p %d,%d' % (value[0], value[1],)
            else:
                raise ValueError, 'Invalid precision argument: ' + str(value)
        return super(Dump, self)._format_arg(name, spec, value)

    def _gen_outfilename(self):
        """
        If the user specified output_file then return that, otherwise
        return the full path to the input file with the extension
        changed to '.txt'.
        """

        output_file = self.inputs.output_file

        if isdefined(output_file):
            return output_file
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '.txt'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

class AverageInputSpec(CommandLineInputSpec):
    _xor_input_files = ('input_files', 'filelist')

    input_files = InputMultiPath(
                    traits.File,
                    desc='input file(s) for averaging',
                    exists=True,
                    mandatory=True,
                    sep=' ',
                    argstr='%s',
                    position=-2,
                    xor=_xor_input_files)

    output_file = File(
                    desc='output file',
                    genfile=True,
                    argstr='%s',
                    position=-1,)

    two = traits.Bool(desc='Produce a MINC 2.0 format output file.', argstr='-2')

    clobber = traits.Bool(desc='Overwrite existing file.', argstr='-clobber', usedefault=True, default_value=True)

    _xor_verbose = ('verbose', 'quiet',)

    verbose = traits.Bool(desc='Print out log messages (default).', argstr='-verbose',  xor=_xor_verbose)
    quiet   = traits.Bool(desc='Do not print out log messages.',    argstr='-quiet',    xor=_xor_verbose)

    debug   = traits.Bool(desc='Print out debugging messages.', argstr='-debug')

    filelist = traits.File(desc='Specify the name of a file containing input file names.', argstr='-filelist %s', exists=True, mandatory=True, xor=_xor_input_files)

    _xor_check_dimensions = ('check_dimensions', 'no_check_dimensions',)

    check_dimensions    = traits.Bool(desc='Check that dimension info matches across files (default).', argstr='-check_dimensions',     xor=_xor_check_dimensions)
    no_check_dimensions = traits.Bool(desc='Do not check dimension info.',                              argstr='-nocheck_dimensions',   xor=_xor_check_dimensions)

    _xor_format = ('format_filetype', 'format_byte', 'format_short',
                   'format_int', 'format_long', 'format_float', 'format_double',
                   'format_signed', 'format_unsigned',)

    format_filetype     = traits.Bool(desc='Use data type of first file (default).',            argstr='-filetype', xor=_xor_format)
    format_byte         = traits.Bool(desc='Write out byte data.',                              argstr='-byte',     xor=_xor_format)
    format_short        = traits.Bool(desc='Write out short integer data.',                     argstr='-short',    xor=_xor_format)
    format_int          = traits.Bool(desc='Write out 32-bit integer data.',                    argstr='-int',      xor=_xor_format)
    format_long         = traits.Bool(desc='Superseded by -int.',                               argstr='-long',     xor=_xor_format)
    format_float        = traits.Bool(desc='Write out single-precision floating-point data.',   argstr='-float',    xor=_xor_format)
    format_double       = traits.Bool(desc='Write out double-precision floating-point data.',   argstr='-double',   xor=_xor_format)
    format_signed       = traits.Bool(desc='Write signed integer data.',                        argstr='-signed',   xor=_xor_format)
    format_unsigned     = traits.Bool(desc='Write unsigned integer data (default).',            argstr='-unsigned', xor=_xor_format)

    max_buffer_size_in_kb = traits.Range(
                                low=0,
                                desc='Specify the maximum size of the internal buffers (in kbytes).',
                                value=4096,
                                argstr='-max_buffer_size_in_kb %d',)

    _xor_normalize = ('normalize', 'nonormalize',)

    normalize   = traits.Bool(desc='Normalize data sets for mean intensity.', argstr='-normalize',   xor=_xor_normalize)
    nonormalize = traits.Bool(desc='Do not normalize data sets (default).',   argstr='-nonormalize', xor=_xor_normalize)

    voxel_range = traits.Tuple(
                traits.Int, traits.Int, argstr='-range %d %d',
                desc='Valid range for output data.')

    sdfile = traits.File(
                desc='Specify an output sd file (default=none).',
                argstr='-sdfile %s')

    _xor_copy_header = ('copy_header', 'no_copy_header')

    copy_header     = traits.Bool(desc='Copy all of the header from the first file (default for one file).',            argstr='-copy_header',   xor=_xor_copy_header)
    no_copy_header  = traits.Bool(desc='Do not copy all of the header from the first file (default for many files)).',  argstr='-nocopy_header', xor=_xor_copy_header)

    avgdim = traits.Str(desc='Specify a dimension along which we wish to average.', argstr='-avgdim %s')

    binarize = traits.Bool(desc='Binarize the volume by looking for values in a given range.', argstr='-binarize')

    binrange = traits.Tuple(
                traits.Float, traits.Float, argstr='-binrange %s %s',
                desc='Specify a range for binarization. Default value: 1.79769e+308 -1.79769e+308.')

    binvalue = traits.Float(desc='Specify a target value (+/- 0.5) for binarization. Default value: -1.79769e+308', argstr='-binvalue %s')
		
    weights = InputMultiPath(
                            traits.Str,
                            desc='Specify weights for averaging ("<w1>,<w2>,...").',
                            sep=',',
                            argstr='-weights %s',)

    width_weighted = traits.Bool(desc='Weight by dimension widths when -avgdim is used.', argstr='-width_weighted', requires=('avgdim',))

class AverageOutputSpec(TraitedSpec):
    output_file = File(desc='output file', exists=True)

class Average(CommandLine):
    """Average a number of MINC files.

    Examples
    --------

    >>> from nipype.interfaces.minc import Dump
    >>> from nipype.testing import minc2Dfile, nonempty_minc_data

    >>> files = [nonempty_minc_data(i) for i in range(3)]
    >>> average = Average(input_files=files, output_file='/tmp/tmp.mnc')
    >>> dump.run() # doctest: +SKIP

    """

    input_spec  = AverageInputSpec
    output_spec = AverageOutputSpec
    _cmd = 'mincaverage'

    def _gen_filename(self, name):
        if name == 'output_file':
            output_file = self.inputs.output_file

            if isdefined(output_file):
                return os.path.abspath(output_file)
            else:
                return aggregate_filename(self.inputs.input_files, 'averaged')
        else:
            raise NotImplemented

    def _gen_outfilename(self):
        return self._gen_filename('output_file')

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
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
                    genfile=True,
                    argstr='%s',
                    position=-1,)

    trace       = traits.Bool(desc='compute the trace (approximate growth and shrinkage) -- FAST',  argstr='-trace')
    determinant = traits.Bool(desc='compute the determinant (exact growth and shrinkage) -- SLOW',  argstr='-determinant')
    translation = traits.Bool(desc='compute translation (structure displacement)',                  argstr='-translation')
    magnitude   = traits.Bool(desc='compute the magnitude of the displacement vector',              argstr='-magnitude')

class BlobOutputSpec(TraitedSpec):
    output_file = File(desc='output file', exists=True)

class Blob(CommandLine):
    """Calculate blobs from minc deformation grids.

    Examples
    --------

    >>> from nipype.interfaces.minc import Blob
    >>> from nipype.testing import minc2Dfile

    >>> blob = Blob(input_file=minc2Dfile, output_file='/tmp/tmp.mnc', trace=True)
    >>> blob.run() # doctest: +SKIP
    """

    input_spec  = BlobInputSpec
    output_spec = BlobOutputSpec
    _cmd = 'mincblob'

    def _gen_filename(self, name):
        if name == 'output_file':
            output_file = self.inputs.output_file

            if isdefined(output_file):
                return os.path.abspath(output_file)
            else:
                return aggregate_filename([self.inputs.input_file], 'blob')
        else:
            raise NotImplemented

    def _gen_outfilename(self):
        return self._gen_filename('output_file')

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

class CalcInputSpec(CommandLineInputSpec):
    _xor_input_files = ('input_files', 'filelist')

    input_files = InputMultiPath(
                    traits.File,
                    desc='input file(s) for calculation',
                    exists=True,
                    mandatory=True,
                    sep=' ',
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file',
                    genfile=True,
                    argstr='%s',
                    position=-1,)

    two = traits.Bool(desc='Produce a MINC 2.0 format output file.', argstr='-2')

    clobber = traits.Bool(desc='Overwrite existing file.', argstr='-clobber', usedefault=True, default_value=True)

    _xor_verbose = ('verbose', 'quiet',)

    verbose = traits.Bool(desc='Print out log messages (default).', argstr='-verbose',  xor=_xor_verbose)
    quiet   = traits.Bool(desc='Do not print out log messages.',    argstr='-quiet',    xor=_xor_verbose)

    debug   = traits.Bool(desc='Print out debugging messages.', argstr='-debug')

    filelist = traits.File(desc='Specify the name of a file containing input file names.', argstr='-filelist %s', mandatory=True, xor=_xor_input_files)

    _xor_copy_header = ('copy_header', 'no_copy_header')

    copy_header     = traits.Bool(desc='Copy all of the header from the first file.',         argstr='-copy_header',   xor=_xor_copy_header)
    no_copy_header  = traits.Bool(desc='Do not copy all of the header from the first file.',  argstr='-nocopy_header', xor=_xor_copy_header)

    _xor_format = ('format_filetype', 'format_byte', 'format_short',
                   'format_int', 'format_long', 'format_float', 'format_double',
                   'format_signed', 'format_unsigned',)

    format_filetype     = traits.Bool(desc='Use data type of first file (default).',            argstr='-filetype', xor=_xor_format)
    format_byte         = traits.Bool(desc='Write out byte data.',                              argstr='-byte',     xor=_xor_format)
    format_short        = traits.Bool(desc='Write out short integer data.',                     argstr='-short',    xor=_xor_format)
    format_int          = traits.Bool(desc='Write out 32-bit integer data.',                    argstr='-int',      xor=_xor_format)
    format_long         = traits.Bool(desc='Superseded by -int.',                               argstr='-long',     xor=_xor_format)
    format_float        = traits.Bool(desc='Write out single-precision floating-point data.',   argstr='-float',    xor=_xor_format)
    format_double       = traits.Bool(desc='Write out double-precision floating-point data.',   argstr='-double',   xor=_xor_format)
    format_signed       = traits.Bool(desc='Write signed integer data.',                        argstr='-signed',   xor=_xor_format)
    format_unsigned     = traits.Bool(desc='Write unsigned integer data (default).',            argstr='-unsigned', xor=_xor_format)

    voxel_range = traits.Tuple(
                traits.Int, traits.Int, argstr='-range %d %d',
                desc='Valid range for output data.',)

    max_buffer_size_in_kb = traits.Range(
                                low = 0,
                                desc='Specify the maximum size of the internal buffers (in kbytes).',
                                value=0,
                                usedefault=False,
                                argstr='-max_buffer_size_in_kb %d')

    _xor_check_dimensions = ('check_dimensions', 'no_check_dimensions',)

    check_dimensions    = traits.Bool(desc='Check that files have matching dimensions (default).',  argstr='-check_dimensions',     xor=_xor_check_dimensions)
    no_check_dimensions = traits.Bool(desc='Do not check that files have matching dimensions.',     argstr='-nocheck_dimensions',   xor=_xor_check_dimensions)

    # FIXME Is it sensible to use ignore_nan and propagate_nan at the same time? Document this.
    ignore_nan = traits.Bool(desc='Ignore invalid data (NaN) for accumulations.', argstr='-ignore_nan')
    propagate_nan = traits.Bool(desc='Invalid data in any file at a voxel produces a NaN (default).', argstr='-propagate_nan')

    # FIXME Double-check that these are mutually exclusive?
    _xor_nan_zero_illegal = ('output_nan', 'output_zero', 'output_illegal_value')

    output_nan      = traits.Bool(desc='Output NaN when an illegal operation is done (default).',                           argstr='-nan',              xor=_xor_nan_zero_illegal)
    output_zero     = traits.Bool(desc='Output zero when an illegal operation is done.',                                    argstr='-zero',             xor=_xor_nan_zero_illegal)
    output_illegal  = traits.Bool(desc='Value to write out when an illegal operation is done. Default value: 1.79769e+308', argstr='-illegal_value',    xor=_xor_nan_zero_illegal)

    _xor_expression = ('expression', 'expfile')

    expression = traits.Str(desc='Expression to use in calculations.',    argstr='-expression \'%s\'', xor=_xor_expression, mandatory=True)
    expfile    = traits.File(desc='Name of file containing expression.',  argstr='-expfile %s',    xor=_xor_expression, mandatory=True)

    # FIXME test this one, the argstr will probably need tweaking, see _format_arg.
    outfiles = traits.List(
                traits.Tuple(traits.Str, traits.File, argstr='-outfile %s %s',
                desc='List of (symbol, file) tuples indicating that output should be written to the specified file, taking values from the symbol which should be created in the expression (see the EXAMPLES section). If this option is given, then all non-option arguments are taken as input files. This option can be used multiple times for multiple output files.'))

    eval_width = traits.Int(200, desc='Number of voxels to evaluate simultaneously.', argstr='-eval_width %s', usedefault=False)

class CalcOutputSpec(TraitedSpec):
    output_file = File(desc='output file', exists=True)

class Calc(CommandLine):
    """Compute an expression using MINC files as input.

    Examples
    --------

    >>> from nipype.interfaces.minc import Calc
    >>> from nipype.testing import minc2Dfile, nonempty_minc_data

    >>> file0 = nonempty_minc_data(0)
    >>> file1 = nonempty_minc_data(1)
    >>> calc = Calc(input_files=[file0, file1], output_file='/tmp/calc.mnc', expression='A[0] + A[1]') # add files together
    >>> calc.run() # doctest: +SKIP
    """

    input_spec  = CalcInputSpec
    output_spec = CalcOutputSpec
    _cmd = 'minccalc'

    def _gen_filename(self, name):
        if name == 'output_file':
            output_file = self.inputs.output_file

            if isdefined(output_file):
                return os.path.abspath(output_file)
            else:
                return aggregate_filename([self.inputs.input_file], 'calc_output')
        else:
            raise NotImplemented

    def _gen_outfilename(self):
        return self._gen_filename('output_file')

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs


# FIXME mincbbox produces output like
#
#   -5.000000 -5.000000 -5.000000    4.800000 2.800000 8.800000
#
# so perhaps this would be better returned as a pair of Python
# lists instead of sending to an output file?

class BBoxInputSpec(StdOutCommandLineInputSpec):
    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file = File(
                    desc='output file containing bounding box corners',
                    position=-1)

    threshold = traits.Int(0, desc='VIO_Real value threshold for bounding box. Default value: 0.', argstr='-threshold')

    _xor_one_two = ('one_line', 'two_lines')

    one_line  = traits.Bool(desc='Output on one line (default): start_x y z width_x y z', argstr='-one_line',  xor=_xor_one_two)
    two_lines = traits.Bool(desc='Output on two lines: start_x y z \n width_x y z',       argstr='-two_lines', xor=_xor_one_two)

    format_mincresample = traits.Bool(desc='Output format for mincresample: (-step x y z -start x y z -nelements x y z',    argstr='-mincresample')
    format_mincreshape  = traits.Bool(desc='Output format for mincreshape: (-start x,y,z -count dx,dy,dz',                  argstr='-mincreshape')
    format_minccrop     = traits.Bool(desc='Output format for minccrop: (-xlim x1 x2 -ylim y1 y2 -zlim z1 z2',              argstr='-minccrop')

    # FIXME Not implemented, will clash with our parsing of the output?
    # Command-specific options:
    # Options for logging progress. Default = -verbose.
    #  -verbose:      Write messages indicating progress
    #  -quiet:        Do not write log messages
    #  -debug:        Print out debug info.

class BBoxOutputSpec(TraitedSpec):
    output_file = File(desc='output file containing bounding box corners', exists=True)

class BBox(StdOutCommandLine):
    """Determine a bounding box.

    FIXME doctests
    """

    input_spec  = BBoxInputSpec
    output_spec = BBoxOutputSpec
    _cmd = 'mincbbox'

    # FIXME Does this play nicely with a workflow?
    def _gen_outfilename(self):
        output_file = self.inputs.output_file

        if isdefined(output_file):
            return output_file
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '_bbox.txt'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs



class BeastInputSpec(CommandLineInputSpec):
    """

    TODO:

    Command-specific options:
     -verbose:          Enable verbose output.
     -positive:         Specify mask of positive segmentation (inside mask) instead of the default mask.
     -output_selection: Specify file to output selected files.
     -count:            Specify file to output the patch count.
     -mask:             Specify a segmentation mask instead of the the default mask.
     -no_mask:          Do not apply a segmentation mask. Perform the segmentation over the entire image.
     -no_positive:      Do not apply a positive mask.
    Generic options for all commands:
     -help:             Print summary of command-line options and abort
     -version:          Print version number of program and exit
    Copyright (C) 2011	Simon Fristed Eskildsen, Vladimir Fonov,
                Pierrick Coupe, Jose V. Manjon

    This program comes with ABSOLUTELY NO WARRANTY; for details type 'cat COPYING'.
    This is free software, and you are welcome to redistribute it under certain
    conditions; type 'cat COPYING' for details.

    Usage: mincbeast [options] <library dir> <input> <output>
           mincbeast -help

    Get this example to work?

    https://github.com/BIC-MNI/BEaST/blob/master/README.library


        2.3 Source the minc-toolkit (if installed):
        $ source /opt/minc/minc-toolkit-config.sh

        2.4 Generate library by running:
        $ beast_prepareADNIlib -flip <ADNI download directory> <BEaST library directory>
        Example:
        $ sudo beast_prepareADNIlib -flip Downloads/ADNI /opt/minc/share/beast-library-1.1

        3. Test the setup
        3.1 Normalize your data
        $ beast_normalize -modeldir /opt/minc/share/icbm152_model_09c input.mnc normal.mnc normal.xfm
        3.2 Run BEaST
        $ mincbeast /opt/minc/share/beast-library-1.1 normal.mnc brainmask.mnc -conf /opt/minc/share/beast-library-1.1/default.2mm.conf -same_res
    """

    probability_map = traits.Bool(desc='Output the probability map instead of crisp mask.', argstr='-probability')
    flip_images = traits.Bool(desc='Flip images around the mid-sagittal plane to increase patch count.', argstr='-flip')
    load_moments = traits.Bool(desc='Do not calculate moments instead use precalculated library moments. (for optimization purposes)', argstr='-load_moments')
    fill_holes = traits.Bool(desc='Fill holes in the binary output.', argstr='-fill')
    median_filter = traits.Bool(desc='Apply a median filter on the probability map.', argstr='-median')
    nlm_filter = traits.Bool(desc='Apply an NLM filter on the probability map (experimental).', argstr='-nlm_filter')

    clobber = traits.Bool(desc='Overwrite existing file.', argstr='-clobber', usedefault=True, default_value=True)

    configuration_file = traits.File(desc='Specify configuration file.', argstr='-configuration %s')

    voxel_size = traits.Int(4, desc='Specify voxel size for calculations (4, 2, or 1). Default value: 4. Assumes no multiscale. Use configuration file for multiscale.', argstr='-voxel_size %s')

    abspath = traits.Bool(desc='File paths in the library are absolute (default is relative to library root).',
                          argstr='-abspath', usedefault=True, default_value=True)

    patch_size = traits.Int(1, desc='Specify patch size for single scale approach. Default value: 1.', argstr='-patch_size %s')

    search_area = traits.Int(2, desc='Specify size of search area for single scale approach. Default value: 2.', argstr='-search_area %s')

    confidence_level_alpha      = traits.Float(0.5, desc='Specify confidence level Alpha. Default value: 0.5',          argstr='-alpha %s')
    smoothness_factor_beta      = traits.Float(0.5, desc='Specify smoothness factor Beta. Default value: 0.25',         argstr='-beta %s')
    threshold_patch_selection   = traits.Float(0.95, desc='Specify threshold for patch selection. Default value: 0.95', argstr='-threshold %s')
    number_selected_images      = traits.Int(20, desc='Specify number of selected images. Default value: 20',           argstr='-selection_num %s')

    same_resolution = traits.Bool(desc='Output final mask with the same resolution as input file.', argstr='-same_resolution')

    library_dir = traits.Directory(desc='library directory', position=-3, argstr='%s', mandatory=True)
    input_file  = traits.File(desc='input file',             position=-2, argstr='%s', mandatory=True)
    output_file = traits.File(desc='output file',            position=-1, argstr='%s', genfile=True)

class BeastOutputSpec(TraitedSpec):
    output_file = File(desc='output file in raw/text format', exists=True)

class Beast(CommandLine):
    """FIXME
    """

    input_spec  = BeastInputSpec
    output_spec = BeastOutputSpec
    _cmd = 'mincbeast'

    # FIXME Does this play nicely with a workflow?
    def _gen_outfilename(self):
        output_file = self.inputs.output_file

        if isdefined(output_file):
            return output_file
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '_mask.mnc'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

class PikInputSpec(CommandLineInputSpec):
    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    _xor_image_type = ('jpg', 'png')

    jpg = traits.Bool(desc='Output a jpg file.',            xor=_xor_image_type)
    png = traits.Bool(desc='Output a png file (default).',  xor=_xor_image_type)

    output_file = File(
                    desc='output file',
                    argstr='%s',
                    position=-1)

    clobber = traits.Bool(
                desc='Overwrite existing file.',
                argstr='-clobber', usedefault=True, default_value=True)

    # FIXME not implemented: --verbose
    #                        --fake
    #                        --lookup    ==> arguments to pass to minclookup

    scale = traits.Int(2, desc='Scaling factor for resulting image, by default images are output at twice their original resolution.', argstr='--scale %s')

    width = traits.Int(desc='Autoscale the resulting image to have a fixed image width (in pixels).', argstr='--width %s')

    depth = traits.Enum(8, 16, desc='Bitdepth for resulting image 8 or 16 (MSB machines only!)', argstr='--depth %s')

    _xor_title = ('title_string', 'title_with_filename')

    title = traits.Either(
                    traits.Bool(desc='Use input filename as title in resulting image.'),
                    traits.Str(desc='Add a title to the resulting image.'),
                    argstr='%s') # see _format_arg for actual arg string

    title_size = traits.Int(desc='Font point size for the title.', argstr='--title_size %s', requires=['title'])

    annotated_bar = traits.Bool(desc='create an annotated bar to match the image (use height of the output image)', argstr='--anot_bar')

    # FIXME tuple of floats? Not voxel values? Man page doesn't specify.
    minc_range = traits.Tuple(
                    traits.Float, traits.Float,
                    desc='Valid range of values for MINC file.',
                    argstr='--range %s %s')

    _xor_image_range = ('image_range', 'auto_range')

    image_range = traits.Tuple(
                    traits.Float, traits.Float,
                    desc='Range of image values to use for pixel intensity.',
                    argstr='--image_range %s %s',
                    xor=_xor_image_range)

    auto_range = traits.Bool(desc='Automatically determine image range using a 5 and 95% PcT. (histogram)', argstr='--auto_range', xor=_xor_image_range)

    start = traits.Int(desc='Slice number to get. (note this is in voxel co-ordinates).', argstr='--slice %s') # FIXME Int is correct?

    _xor_slice = ('slice_z', 'slice_y', 'slice_x')

    slice_z = traits.Bool(desc='Get an axial/transverse (z) slice.',    argstr='-z', xor=_xor_slice)
    slice_y = traits.Bool(desc='Get a coronal (y) slice.',              argstr='-y', xor=_xor_slice)
    slice_x = traits.Bool(desc='Get a sagittal (x) slice.',             argstr='-x', xor=_xor_slice) # FIXME typo in man page? sagital?

    triplanar = traits.Bool(desc='Create a triplanar view of the input file.', argstr='--triplanar')
    tile_size = traits.Int(desc='Pixel size for each image in a triplanar.', argstr='--tilesize')

    _xor_sagittal_offset = ('sagittal_offset', 'sagittal_offset_perc')

    sagittal_offset = traits.Int(desc='Offset the sagittal slice from the centre.', argstr='--sagittal_offset')
    sagittal_offset_perc = traits.Range(low=0, high=100,
                        desc='Offset the sagittal slice by a percentage from the centre.',
                        argstr='--sagittal_offset_perc %d',)

    _xor_vertical_horizontal = ('vertical_triplanar_view', 'horizontal_triplanar_view')

    vertical_triplanar_view   = traits.Bool(desc='Create a vertical triplanar view (Default).', argstr='--vertical',   xor=_xor_vertical_horizontal)
    horizontal_triplanar_view = traits.Bool(desc='Create a horizontal triplanar view.',         argstr='--horizontal', xor=_xor_vertical_horizontal)

class PikOutputSpec(TraitedSpec):
    output_file = File(desc='output image', exists=True)

class Pik(CommandLine):
    """FIXME
    """

    input_spec  = PikInputSpec
    output_spec = PikOutputSpec
    _cmd = 'mincpik'

    # FIXME Does this play nicely with a workflow?
    def _gen_outfilename(self):
        output_file = self.inputs.output_file

        if isdefined(output_file):
            assert not isdefined(self.inputs.png) # FIXME make a warning instead?
            assert not isdefined(self.inputs.jpg) # FIXME make a warning instead?
            return output_file
        else:
            b = os.path.splitext(self.inputs.input_file)[0]

            if isdefined(self.inputs.png) and self.inputs.png:
                return b + '.png'
            elif isdefined(self.inputs.jpg) and self.inputs.jpg:
                return b + '.jpg'
            else:
                # By default we'll write a png file.
                return b + '.png'

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_file'] = os.path.abspath(self._gen_outfilename())
        return outputs

    def _format_arg(self, name, spec, value):
        if name == 'title':
            if isinstance(value, bool) and value:
                return '--title'
            elif isinstance(value, str):
                return '--title --title_text %s' % (value,)
            else:
                raise ValueError, 'Unknown value for "title" argument: ' + str(value)
        return super(Pik, self)._format_arg(name, spec, value)


    @property
    def cmdline(self):
        output_file = self.inputs.output_file

        if isdefined(output_file):
            return super(Pik, self).cmdline
        else:
            # FIXME this seems like a bit of a hack. Can we force output_file
            # to show up in cmdline by default, even if it isn't specified in
            # the instantiation of Pik?
            return '%s %s' % (super(Pik, self).cmdline, self._gen_outfilename())

class BlurInputSpec(CommandLineInputSpec):
    """ FIXME not implemented
    -verbose:     Write messages indicating progress
    -quiet:       Do not write log messages
    -debug:       Print out debug info.
    -version:     Print out version info and exit.
    """

    input_file = File(
                    desc='input file',
                    exists=True,
                    mandatory=True,
                    argstr='%s',
                    position=-2,)

    output_file_base = File(
                    desc='output file base',
                    argstr='%s',
                    position=-1)

    clobber = traits.Bool(
                desc='Overwrite existing file.',
                argstr='-clobber', usedefault=True, default_value=True)

    _xor_kernel = ('gaussian', 'rect')

    gaussian = traits.Bool(desc='Use a gaussian smoothing kernel (default).', argstr='-gaussian', xor=_xor_kernel)
    rect     = traits.Bool(desc='Use a rect (box) smoothing kernel.',         argstr='-rect',     xor=_xor_kernel)

    gradient   = traits.Bool(desc='Create the gradient magnitude volume as well.',                          argstr='-gradient')
    partial    = traits.Bool(desc='Create the partial derivative and gradient magnitude volumes as well.',  argstr='-partial')

    no_apodize = traits.Bool(desc='Do not apodize the data before blurring.',                               argstr='-no_apodize')

    _xor_main_options = ('fwhm', 'fwhm3d', 'standard_dev')

    fwhm = traits.Float(
                0,
                desc='Full-width-half-maximum of gaussian kernel. Default value: 0.',
                argstr='-fwhm %s',
                xor=_xor_main_options,
                mandatory=True)

    standard_dev = traits.Float(
                    0,
                    desc='Standard deviation of gaussian kernel. Default value: 0.',
                    argstr='-standarddev %s',
                    xor=_xor_main_options,
                    mandatory=True)

    fwhm3d = traits.Tuple(
                    traits.Float, traits.Float, traits.Float,
                    argstr='-3dfwhm %s %s %s',
                    desc='Full-width-half-maximum of gaussian kernel. Default value: -1.79769e+308 -1.79769e+308 -1.79769e+308.',
                    xor=_xor_main_options,
                    mandatory=True)

    dimensions = traits.Enum(1, 2, 3, desc='Number of dimensions to blur (either 1,2 or 3). Default value: 3.', argstr='-dimensions %s', default=3)

class BlurOutputSpec(TraitedSpec):
    output_file = File(desc='Blurred output file.', exists=True)

    gradient_dxyz   = File(desc='Gradient dxyz.')
    partial_dx      = File(desc='Partial gradient dx.')
    partial_dy      = File(desc='Partial gradient dy.')
    partial_dz      = File(desc='Partial gradient dz.')
    partial_dxyz    = File(desc='Partial gradient dxyz.')

class Blur(StdOutCommandLine):
    """
    Convolve an input volume with a Gaussian blurring kernel of
    user-defined width.  Optionally, the first partial derivatives
    and the gradient magnitude volume can be calculated.

    Examples
    --------

    >>> from nipype.interfaces.minc import Blur
    >>> from nipype.testing import minc3Dfile

    (1) Blur  an  input  volume with a 6mm fwhm isotropic Gaussian
    blurring kernel:

    >>> blur = Blur(input_file=minc3Dfile, fwhm=6, output_file_base='/tmp/out_6')
    >>> extract.run() # doctest: +SKIP

    mincblur will create /tmp/out_6_blur.mnc.

    (2) Calculate the blurred and gradient magnitude data:

    >>> blur = Blur(input_file=minc3Dfile, fwhm=6, gradient=True, output_file_base='/tmp/out_6')
    >>> extract.run() # doctest: +SKIP

    will create /tmp/out_6_blur.mnc and /tmp/out_6_dxyz.mnc.

    (3) Calculate the blurred data, the partial derivative volumes
    and  the gradient magnitude for the same data:

    >>> blur = Blur(input_file=minc3Dfile, fwhm=6, partial=True, output_file_base='/tmp/out_6')
    >>> extract.run() # doctest: +SKIP

    will create /tmp/out_6_blur.mnc, /tmp/out_6_dx.mnc,
    /tmp/out_6_dy.mnc, /tmp/out_6_dz.mnc and /tmp/out_6_dxyz.mnc.
    """

    input_spec  = BlurInputSpec
    output_spec = BlurOutputSpec
    _cmd = 'mincblur'

    def _gen_output_base(self):
        output_file_base = self.inputs.output_file_base

        if isdefined(output_file_base):
            return output_file_base
        else:
            return os.path.splitext(self.inputs.input_file)[0] + '_bluroutput'

    def _list_outputs(self):
        outputs = self.output_spec().get()

        output_file_base = self._gen_output_base()

        outputs['output_file'] = output_file_base + '_blur.mnc'

        if isdefined(self.inputs.gradient):
            outputs['gradient_dxyz'] = output_file_base + '_dxyz.mnc'

        if isdefined(self.inputs.partial):
            outputs['partial_dx']   = output_file_base + '_dx.mnc'
            outputs['partial_dy']   = output_file_base + '_dy.mnc'
            outputs['partial_dz']   = output_file_base + '_dz.mnc'
            outputs['partial_dxyz'] = output_file_base + '_dxyz.mnc'

        return outputs

    @property
    def cmdline(self):
        output_file_base = self.inputs.output_file_base

        if isdefined(output_file_base):
            return super(Blur, self).cmdline
        else:
            # FIXME this seems like a bit of a hack. Can we force output_file
            # to show up in cmdline by default, even if it isn't specified in
            # the instantiation of Pik?
            return '%s %s' % (super(Blur, self).cmdline, self._gen_output_base())

"""
Command-specific options:
General options:
 -dimension:             Specify a dimension along which we wish to perform a calculation.
 -ignore_nan:            Ignore invalid data (NaN) for accumulations.
 -propagate_nan:         Invalid data in any file at a voxel produces a NaN (default).
 -nan:                   Output NaN when an illegal operation is done (default).
 -zero:                  Output zero when an illegal operation is done.
 -illegal_value:         Value to write out when an illegal operation is done.
		Default value: 1.79769e+308
Options for specifying constants:
 -constant:              Specify a constant argument.
		Default value: 1.79769e+308
 -const:                 Synonym for -constant.
		Default value: 1.79769e+308
 -const2:                Specify two constant arguments.
		Default value: 1.79769e+308 1.79769e+308
Operations:
 -invert:                Calculate 1/x at each voxel (use -constant for c/x).
 -sqrt:                  Take square root of a volume.
 -square:                Take square of a volume.
 -abs:                   Take absolute value of a volume.
 -max:                   Synonym for -maximum.
 -maximum:               Find maximum of N volumes.
 -minimum:               Find minimum of N volumes.
 -exp:                   Calculate c2*exp(c1*x). The constants c1 and c2 default to 1.
 -log:                   Calculate log(x/c2)/c1. The constants c1 and c2 default to 1.
 -scale:                 Scale a volume: volume * c1 + c2.
 -clamp:                 Clamp a volume to lie between two values.
 -segment:               Segment a volume using range of -const2: within range = 1, outside range = 0.
 -nsegment:              Opposite of -segment: within range = 0, outside range = 1.
 -percentdiff:           Percent difference between 2 volumes, thresholded (const def=0.0).
 -pd:                    Synonym for -percentdiff.
 -isnan:                 Test for NaN values in vol1.
 -nisnan:                Negation of -isnan.
 -count_valid:           Count the number of valid values in N volumes.
Generic options for all commands:
 -help:                  Print summary of command-line options and abort
 -version:               Print version number of program and exit

Usage: mincmath [options] [<in1.mnc> ...] <out.mnc>
       mincmath -help

"""

class MathInputSpec(CommandLineInputSpec):
    """
    FIXME Not implemented
    -verbose:               Print out log messages (default).
    -quiet:                 Do not print out log messages.
    -debug:                 Print out debugging messages.

    """

    _xor_input_files = ('input_files', 'filelist')

    input_files = InputMultiPath(
                    traits.File,
                    desc='input file(s) for calculation',
                    exists=True,
                    mandatory=True,
                    sep=' ',
                    argstr='%s',
                    position=-2,
                    xor=_xor_input_files)

    filelist = traits.File(desc='Specify the name of a file containing input file names.', argstr='-filelist %s', exists=True, mandatory=True, xor=_xor_input_files)

    clobber = traits.Bool(
                desc='Overwrite existing file.',
                argstr='-clobber', usedefault=True, default_value=True)

    two = traits.Bool(
                desc='Create a MINC 2 output file.',
                argstr='-2',)

    _xor_copy_header = ('copy_header', 'no_copy_header')

    copy_header     = traits.Bool(desc='Copy all of the header from the first file (default for one file).',            argstr='-copy_header',   xor=_xor_copy_header)
    no_copy_header  = traits.Bool(desc='Do not copy all of the header from the first file (default for many files)).',  argstr='-nocopy_header', xor=_xor_copy_header)

    _xor_format = ('format_filetype', 'format_byte', 'format_short',
                   'format_int', 'format_long', 'format_float', 'format_double',
                   'format_signed', 'format_unsigned',)

    format_filetype     = traits.Bool(desc='Use data type of first file (default).',            argstr='-filetype', xor=_xor_format)
    format_byte         = traits.Bool(desc='Write out byte data.',                              argstr='-byte',     xor=_xor_format)
    format_short        = traits.Bool(desc='Write out short integer data.',                     argstr='-short',    xor=_xor_format)
    format_int          = traits.Bool(desc='Write out 32-bit integer data.',                    argstr='-int',      xor=_xor_format)
    format_long         = traits.Bool(desc='Superseded by -int.',                               argstr='-long',     xor=_xor_format)
    format_float        = traits.Bool(desc='Write out single-precision floating-point data.',   argstr='-float',    xor=_xor_format)
    format_double       = traits.Bool(desc='Write out double-precision floating-point data.',   argstr='-double',   xor=_xor_format)
    format_signed       = traits.Bool(desc='Write signed integer data.',                        argstr='-signed',   xor=_xor_format)
    format_unsigned     = traits.Bool(desc='Write unsigned integer data (default).',            argstr='-unsigned', xor=_xor_format)

    voxel_range = traits.Tuple(
                traits.Int, traits.Int, argstr='-range %d %d',
                desc='Valid range for output data.')

    max_buffer_size_in_kb = traits.Range(
                                low=0,
                                desc='Specify the maximum size of the internal buffers (in kbytes).',
                                value=4096,
                                argstr='-max_buffer_size_in_kb %d',)

    _xor_check_dimensions = ('check_dimensions', 'no_check_dimensions',)

    check_dimensions    = traits.Bool(desc='Check that dimension info matches across files (default).', argstr='-check_dimensions',     xor=_xor_check_dimensions)
    no_check_dimensions = traits.Bool(desc='Do not check dimension info.',                              argstr='-nocheck_dimensions',   xor=_xor_check_dimensions)

    test_gt = traits.Either(traits.Bool(), traits.Float(), desc='Test for vol1 > vol2 or vol1 > constant.',             argstr='%s')
    test_lt = traits.Either(traits.Bool(), traits.Float(), desc='Test for vol1 < vol2 or vol1 < constant.',             argstr='%s')
    test_eq = traits.Either(traits.Bool(), traits.Float(), desc='Test for integer vol1 == vol2 or vol1 == constant.',   argstr='%s')
    test_ne = traits.Either(traits.Bool(), traits.Float(), desc='Test for integer vol1 != vol2 or vol1 != const.',      argstr='%s')
    test_ge = traits.Either(traits.Bool(), traits.Float(), desc='Test for vol1 >= vol2 or vol1 >= const.',              argstr='%s')
    test_le = traits.Either(traits.Bool(), traits.Float(), desc='Test for vol1 <= vol2 or vol1 <= const.',              argstr='%s')

    calc_and = traits.Bool(desc='Calculate vol1 && vol2 (&& ...).', argstr='-and')
    calc_or  = traits.Bool(desc='Calculate vol1 || vol2 (|| ...).', argstr='-or')
    calc_not = traits.Bool(desc='Calculate !vol1.',                 argstr='-not')

    calc_add = traits.Either(traits.Bool(), traits.Float(), desc='Add N volumes or volume + constant.',         argstr='%s')
    calc_sub = traits.Either(traits.Bool(), traits.Float(), desc='Subtract 2 volumes or volume - constant.',    argstr='%s')
    calc_mul = traits.Either(traits.Bool(), traits.Float(), desc='Multiply N volumes or volume * constant.',    argstr='%s')
    calc_div = traits.Either(traits.Bool(), traits.Float(), desc='Divide 2 volumes or volume / constant.',    argstr='%s')

class MathOutputSpec(TraitedSpec):
    output_file = File(desc='output file', exists=True)

class Math(StdOutCommandLine):
    input_spec  = MathInputSpec
    output_spec = MathOutputSpec
    _cmd = 'mincmath'

    def _format_arg(self, name, spec, value):
        if name == 'test_gt':
            if isinstance(value, bool):
                return '-gt'
            elif isinstance(value, float):
                return '-gt -const %s' % value
            else:
                raise ValueError, 'Invalid gt argument: ' + str(value)
        if name == 'test_lt':
            if isinstance(value, bool):
                return '-lt'
            elif isinstance(value, float):
                return '-lt -const %s' % value
            else:
                raise ValueError, 'Invalid lt argument: ' + str(value)
        if name == 'test_eq':
            if isinstance(value, bool):
                return '-eq'
            elif isinstance(value, float):
                return '-eq -const %s' % value
            else:
                raise ValueError, 'Invalid eq argument: ' + str(value)
        if name == 'test_ne':
            if isinstance(value, bool):
                return '-ne'
            elif isinstance(value, float):
                return '-ne -const %s' % value
            else:
                raise ValueError, 'Invalid ne argument: ' + str(value)
        if name == 'test_ge':
            if isinstance(value, bool):
                return '-ge'
            elif isinstance(value, float):
                return '-ge -const %s' % value
            else:
                raise ValueError, 'Invalid ge argument: ' + str(value)
        if name == 'test_le':
            if isinstance(value, bool):
                return '-le'
            elif isinstance(value, float):
                return '-le -const %s' % value
            else:
                raise ValueError, 'Invalid le argument: ' + str(value)
        if name == 'calc_add':
            if isinstance(value, bool):
                return '-add'
            elif isinstance(value, float):
                return '-add -const %s' % value
            else:
                raise ValueError, 'Invalid add argument: ' + str(value)
        if name == 'calc_sub':
            if isinstance(value, bool):
                return '-sub'
            elif isinstance(value, float):
                return '-sub -const %s' % value
            else:
                raise ValueError, 'Invalid sub argument: ' + str(value)
        if name == 'calc_mul':
            if isinstance(value, bool):
                return '-mult'
            elif isinstance(value, float):
                return '-mult -const %s' % value
            else:
                raise ValueError, 'Invalid mul argument: ' + str(value)
        if name == 'calc_div':
            if isinstance(value, bool):
                return '-div'
            elif isinstance(value, float):
                return '-div -const %s' % value
            else:
                raise ValueError, 'Invalid div argument: ' + str(value)

        return super(Math, self)._format_arg(name, spec, value)

# TODO from volgenmodel:
# mincnorm  ??? Not in my installation of MINC.
# mincpik
# mincmath
# mincresample
