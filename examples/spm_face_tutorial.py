"""
Using SPM for analysis
=======================

The spm_auditory_tutorial.py recreates the classical workflow described in the SPM8 manual (http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf)
using auditory dataset that can be downloaded from http://www.fil.ion.ucl.ac.uk/spm/data/face_rep/face_rep_SPM5.html:

    python spm_tutorial.py

"""
from copy import deepcopy


"""Import necessary modules from nipype."""

import nipype.interfaces.io as nio           # Data i/o 
import nipype.interfaces.spm as spm          # spm
import nipype.interfaces.matlab as mlab      # how to run matlab
import nipype.interfaces.fsl as fsl          # fsl
import nipype.interfaces.utility as util     # utility 
import nipype.pipeline.engine as pe          # pypeline engine
import nipype.algorithms.modelgen as model   # model specification
import nipype.externals.pynifti as ni
import os                                    # system functions

"""

Preliminaries
-------------

Confirm package dependencies are installed.  (This is only for the
tutorial, rarely would you put this in your own code.)
"""

from nipype.utils.misc import package_check

package_check('numpy', '1.3', 'tutorial1')
package_check('scipy', '0.7', 'tutorial1')
package_check('networkx', '1.0', 'tutorial1')
package_check('IPython', '0.10', 'tutorial1')

"""Set any package specific configuration. The output file format
for FSL routines is being set to uncompressed NIFTI and a specific
version of matlab is being used. The uncompressed format is required
because SPM does not handle compressed NIFTI.
"""

# Tell fsl to generate all output in uncompressed nifti format
print fsl.Info.version()
fsl.FSLCommand.set_default_output_type('NIFTI')

# Set the way matlab should be called
mlab.MatlabCommand.set_default_matlab_cmd("matlab -nodesktop -nosplash")


"""
Setting up workflows
--------------------

In this tutorial we will be setting up a hierarchical workflow for spm
analysis. This will demonstrate how pre-defined workflows can be setup
and shared across users, projects and labs.


Setup preprocessing workflow
----------------------------

This is a generic preprocessing workflow that can be used by different analyses

"""

preproc = pe.Workflow(name='preproc')

"""Use :class:`nipype.interfaces.spm.Realign` for motion correction
and register all images to the mean image.
"""

realign = pe.Node(interface=spm.Realign(), name="realign")

slice_timing = pe.Node(interface=spm.SliceTiming(), name="slice_timing")


"""Use :class:`nipype.interfaces.spm.Coregister` to perform a rigid
body registration of the functional data to the structural data.
"""

coregister = pe.Node(interface=spm.Coregister(), name="coregister")
coregister.inputs.jobtype = 'estimate'



segment = pe.Node(interface=spm.Segment(), name="segment")
# Uncomment the following line for faster execution
#segment.inputs.gaussians_per_class = [1, 1, 1, 4]

"""Warp functional and structural data to SPM's T1 template using
:class:`nipype.interfaces.spm.Normalize`.  The tutorial data set
includes the template image, T1.nii.
"""

normalize_func = pe.Node(interface=spm.Normalize(), name = "normalize_func")
normalize_func.inputs.jobtype = "write"

normalize_struc = pe.Node(interface=spm.Normalize(), name = "normalize_struc")
normalize_struc.inputs.jobtype = "write"


"""Smooth the functional data using
:class:`nipype.interfaces.spm.Smooth`.
"""

smooth = pe.Node(interface=spm.Smooth(), name = "smooth")


def get_vox_dims(volume):
    if isinstance(volume, list):
        volume = volume[0]
    nii = ni.load(volume)
    hdr = nii.get_header()
    voxdims = hdr.get_zooms()
    return [float(voxdims[0]), float(voxdims[1]), float(voxdims[2])]

preproc.connect([(realign,coregister,[('mean_image', 'target')]),
                 (segment, normalize_func, [('transformation_mat','parameter_file')]),
                 (segment, normalize_struc, [('transformation_mat','parameter_file'),
                                             ('modulated_input_image', 'apply_to_files'),
                                             (('modulated_input_image', get_vox_dims), 'write_voxel_sizes')]),
                 (realign, slice_timing, [('realigned_files', 'in_files')]),
                 (slice_timing, normalize_func, [('timecorrected_files', 'apply_to_files'),
                                            (('timecorrected_files', get_vox_dims), 'write_voxel_sizes')]),
                 (normalize_func, smooth, [('normalized_files', 'in_files')]),
                 ])


"""
Set up analysis workflow
------------------------

"""

l1analysis = pe.Workflow(name='analysis')

"""Generate SPM-specific design information using
:class:`nipype.interfaces.spm.SpecifyModel`.
"""

modelspec = pe.Node(interface=model.SpecifyModel(), name= "modelspec")

"""Generate a first level SPM.mat file for analysis
:class:`nipype.interfaces.spm.Level1Design`.
"""

level1design = pe.Node(interface=spm.Level1Design(), name= "level1design")

"""Use :class:`nipype.interfaces.spm.EstimateModel` to determine the
parameters of the model.
"""

level1estimate = pe.Node(interface=spm.EstimateModel(), name="level1estimate")
level1estimate.inputs.estimation_method = {'Classical' : 1}

threshold = pe.Node(interface=spm.Threshold(), name="threshold")


"""Use :class:`nipype.interfaces.spm.EstimateContrast` to estimate the
first level contrasts specified in a few steps above.
"""

contrastestimate = pe.Node(interface = spm.EstimateContrast(), name="contrastestimate")

l1analysis.connect([(modelspec,level1design,[('session_info','session_info')]),
                  (level1design,level1estimate,[('spm_mat_file','spm_mat_file')]),
                  (level1estimate,contrastestimate,[('spm_mat_file','spm_mat_file'),
                                                  ('beta_images','beta_images'),
                                                  ('residual_image','residual_image')]),
                  (contrastestimate, threshold,[('spm_mat_file','spm_mat_file'),
                                                    ('spmT_images', 'spmT_images')])
                  ])

"""
Preproc + Analysis pipeline
---------------------------

"""

def makelist(item):
    return [item]

l1pipeline = pe.Workflow(name='firstlevel')
l1pipeline.connect([(preproc, l1analysis, [('realign.realignment_parameters',
                                            'modelspec.realignment_parameters'),
                                           (('smooth.smoothed_files', makelist),
                                            'modelspec.functional_runs')])
                  ])


"""
Data specific components
------------------------

The nipype tutorial contains data for two subjects.  Subject data
is in two subdirectories, ``s1`` and ``s2``.  Each subject directory
contains four functional volumes: f3.nii, f5.nii, f7.nii, f10.nii. And
one anatomical volume named struct.nii.

Below we set some variables to inform the ``datasource`` about the
layout of our data.  We specify the location of the data, the subject
sub-directories and a dictionary that maps each run to a mnemonic (or
field) for the run type (``struct`` or ``func``).  These fields become
the output fields of the ``datasource`` node in the pipeline.

In the example below, run 'f3' is of type 'func' and gets mapped to a
nifti filename through a template '%s.nii'. So 'f3' would become
'f3.nii'.

"""

# Specify the location of the data downloaded from http://www.fil.ion.ucl.ac.uk/spm/data/face_rep/face_rep_SPM5.html
data_dir = os.path.abspath('spm_face_data')
# Specify the subject directories
subject_list = ['M03953']
# Map field names to individual subject runs.
info = dict(func=[['RawEPI', 'subject_id', 5, ["_%04d"%i for i in range(6,357)]]],
            struct=[['Structural', 'subject_id', 7, '']])

infosource = pe.Node(interface=util.IdentityInterface(fields=['subject_id']), name="infosource")

"""Here we set up iteration over all the subjects. The following line
is a particular example of the flexibility of the system.  The
``datasource`` attribute ``iterables`` tells the pipeline engine that
it should repeat the analysis on each of the items in the
``subject_list``. In the current example, the entire first level
preprocessing and estimation will be repeated for each subject
contained in subject_list.
"""

infosource.iterables = ('subject_id', subject_list)

"""
Now we create a :class:`nipype.interfaces.io.DataGrabber` object and
fill in the information from above about the layout of our data.  The
:class:`nipype.pipeline.NodeWrapper` module wraps the interface object
and provides additional housekeeping and pipeline specific
functionality.
"""

datasource = pe.Node(interface=nio.DataGrabber(infields=['subject_id'],
                                               outfields=['func', 'struct']),
                     name = 'datasource')
datasource.inputs.base_directory = data_dir
datasource.inputs.template = '%s/s%s_%04d%s.img'
datasource.inputs.template_args = info



"""
Experimental paradigm specific components
-----------------------------------------

Here we create a structure that provides information
about the experimental paradigm. This is used by the
:class:`nipype.interfaces.spm.SpecifyModel` to create the information
necessary to generate an SPM design matrix.
"""

from nipype.interfaces.base import Bunch
from scipy.io.matlab import loadmat

mat = loadmat(os.path.join(data_dir, "sots.mat"))
sot = mat['sot'][0]
itemlag = mat['itemlag'][0]

subjectinfo = [Bunch(conditions=['N1', 'N2', 'F1', 'F2'],
                            onsets=[sot[0], sot[1], sot[2], sot[3]],
                            durations=[[0], [0], [0], [0]],
                            amplitudes=None,
                            tmod=None,
                            pmod=None,
                            regressor_names=None,
                            regressors=None)]

#TODO fix the parametric one
subjectinfo_param = [Bunch(conditions=['N1', 'N2', 'F1', 'F2'],
                            onsets=[sot[0], sot[1], sot[2], sot[3]],
                            durations=[[0], [0], [0], [0]],
                            amplitudes=None,
                            tmod=None,
                            pmod=[None,
                                  Bunch(name=['Lag'],
                                        param=itemlag[1].tolist(),
                                        poly=[2]),
                                  None,
                                  Bunch(name=['Lag'],
                                        param=itemlag[3].tolist(),
                                        poly=[2])],
                            regressor_names=None,
                            regressors=None)]

"""Setup the contrast structure that needs to be evaluated. This is a
list of lists. The inner list specifies the contrasts and has the
following format - [Name,Stat,[list of condition names],[weights on
those conditions]. The condition names must match the `names` listed
in the `subjectinfo` function described above.
"""

cont1 = ('positive effect of condition_1','T', ['N1','N2','F1','F2'], [1,1,1,1])
cont2 = ('positive effect of Fame_1','T', ['N1','N2','F1','F2'],[1,1,-1,-1])
cont3 = ('positive effect of Rep_1','T', ['N1','N2','F1','F2'],[1,-1,1,-1])
cont4 = ('positive interaction: Fame xRep_1','T', ['N1','N2','F1','F2'],[1,-1,-1,1])

fcont1 = ('Average effect of condition', 'F', [cont1])
fcont2 = ('Main effect of Fame', 'F', [cont2])
fcont3 = ('Main effect of Rep', 'F', [cont3])
fcont4 = ('Interaction: Fame x Rep', 'F', [cont4])

contrasts = [cont1, cont2, cont3, cont4, fcont1, fcont2, fcont3, fcont4]

"""
parametric f-contrast
"""
cont1 = ('Famous_lag1','T', ['F2xLag^1'],[1])
cont2 = ('Famous_lag2','T', ['F2xLag^2'],[1])
fcont1 = ('Famous Lag', 'F', [cont1, cont2])
paramcontrasts = [cont1, cont2, fcont1]

num_slices = 24
TR = 2.

slice_timingref = l1pipeline.inputs.preproc.slice_timing
slice_timingref.num_slices = num_slices
slice_timingref.time_repetition = TR
slice_timingref.time_acquisition = TR - TR/float(num_slices)
slice_timingref.slice_order = range(num_slices,0,-1)
slice_timingref.ref_slice = num_slices/2

l1pipeline.inputs.preproc.smooth.fwhm = [8, 8, 8]

# set up node specific inputs
modelspecref = l1pipeline.inputs.analysis.modelspec
modelspecref.input_units             = 'scans'
modelspecref.output_units            = 'scans'
modelspecref.time_repetition         = TR

l1designref = l1pipeline.inputs.analysis.level1design
l1designref.timing_units       = modelspecref.output_units
l1designref.interscan_interval = modelspecref.time_repetition
l1designref.microtime_resolution = slice_timingref.num_slices
l1designref.microtime_onset = slice_timingref.ref_slice
l1designref.bases = {'hrf':{'derivs': [1,1]}}

#l1designref.factor_info = [dict(name='Fame', levels = 2),
#                           dict(name = 'Rep', levels = 2)]

l1pipeline.inputs.analysis.modelspec.subject_info = subjectinfo
l1pipeline.inputs.analysis.contrastestimate.contrasts = contrasts
l1pipeline.inputs.analysis.threshold.contrast_index = 1

paramanalysis = l1analysis.clone(name='paramanalysis')
l1pipeline.connect([(preproc, paramanalysis, [('realign.realignment_parameters',
                                               'modelspec.realignment_parameters'),
                                              (('smooth.smoothed_files', makelist),
                                               'modelspec.functional_runs')])
                  ])

paramanalysis.inputs.level1design.bases = {'hrf':{'derivs': [0,0]}}
paramanalysis.inputs.modelspec.subject_info = subjectinfo_param
paramanalysis.inputs.contrastestimate.contrasts = paramcontrasts

                 
"""
Setup the pipeline
------------------

The nodes created above do not describe the flow of data. They merely
describe the parameters used for each function. In this section we
setup the connections between the nodes such that appropriate outputs
from nodes are piped into appropriate inputs of other nodes.

Use the :class:`nipype.pipeline.engine.Pipeline` to create a
graph-based execution pipeline for first level analysis. The config
options tells the pipeline engine to use `workdir` as the disk
location to use when running the processes and keeping their
outputs. The `use_parameterized_dirs` tells the engine to create
sub-directories under `workdir` corresponding to the iterables in the
pipeline. Thus for this pipeline there will be subject specific
sub-directories.

The ``nipype.pipeline.engine.Pipeline.connect`` function creates the
links between the processes, i.e., how data should flow in and out of
the processing nodes.
"""

level1 = pe.Workflow(name="level1")
level1.base_dir = os.path.abspath('spm_face_tutorial/workingdir')

level1.connect([(infosource, datasource, [('subject_id', 'subject_id')]),
                (datasource,l1pipeline,[('func','preproc.realign.in_files'),
                                        ('struct', 'preproc.coregister.source'),
                                        ('struct', 'preproc.segment.data')]),
                (infosource,l1pipeline,[('subject_id','analysis.modelspec.subject_id')]),
                ])


"""

Setup storage results
---------------------

Use :class:`nipype.interfaces.io.DataSink` to store selected outputs
from the pipeline in a specific location. This allows the user to
selectively choose important output bits from the analysis and keep
them.

The first step is to create a datasink node and then to connect
outputs from the modules above to storage locations. These take the
following form directory_name[.[@]subdir] where parts between [] are
optional. For example 'realign.@mean' below creates a directory called
realign in 'l1output/subject_id/' and stores the mean image output
from the Realign process in the realign directory. If the @ is left
out, then a sub-directory with the name 'mean' would be created and
the mean image would be copied to that directory.
"""

datasink = pe.Node(interface=nio.DataSink(), name="datasink")
datasink.inputs.base_directory = os.path.abspath('spm_auditory_tutorial/l1output')

def getstripdir(subject_id):
    return os.path.join(os.path.abspath('spm_auditory_tutorial/workingdir'),'_subject_id_%s' % subject_id)

# store relevant outputs from various stages of the 1st level analysis
level1.connect([(infosource, datasink,[('subject_id','container'),
                                       (('subject_id', getstripdir),'strip_dir')]),
                (l1pipeline, datasink,[('analysis.contrastestimate.con_images','contrasts.@con'),
                                       ('analysis.contrastestimate.spmT_images','contrasts.@T'),
                                       ('paramanalysis.contrastestimate.con_images','paramcontrasts.@con'),
                                       ('paramanalysis.contrastestimate.spmT_images','paramcontrasts.@T')]),
                ])


"""
Execute the pipeline
--------------------

The code discussed above sets up all the necessary data structures
with appropriate parameters and the connectivity between the
processes, but does not generate any output. To actually run the
analysis on the data the ``nipype.pipeline.engine.Pipeline.Run``
function needs to be called.
"""

if __name__ == '__main__':
    level1.run()
    level1.write_graph()
