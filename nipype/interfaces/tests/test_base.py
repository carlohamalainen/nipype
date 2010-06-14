# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import nipype.interfaces.base as nii
from nipype.testing import assert_equal, assert_not_equal, assert_raises, assert_true
import os
import tempfile
import shutil
from nipype.testing import (assert_equal, assert_not_equal, assert_raises,
                            assert_true, assert_false, with_setup, package_check, skipif)
import nipype.interfaces.base as nib
from nipype.interfaces.base import Undefined
from nipype.interfaces.base import InterfaceResult

#test Bunch
def test_bunch():
    b = nii.Bunch()
    yield assert_equal, b.__dict__,{}
    b = nii.Bunch(a=1,b=[2,3])
    yield assert_equal, b.__dict__,{'a': 1, 'b': [2,3]}

def test_bunch_attribute():
    b = nii.Bunch(a=1,b=[2,3],c=None)
    yield assert_equal, b.a ,1
    yield assert_equal, b.b, [2,3]
    yield assert_equal, b.c, None

def test_bunch_repr():
    b = nii.Bunch(b=2,c=3,a=dict(n=1,m=2))
    yield assert_equal, repr(b), "Bunch(a={'m': 2, 'n': 1}, b=2, c=3)"

def test_bunch_methods():
    b = nii.Bunch(a=2)
    b.update(a=3)
    newb = b.dictcopy()
    yield assert_equal, b.a, 3
    yield assert_equal, b.get('a'), 3
    yield assert_equal, b.get('badkey', 'otherthing'), 'otherthing'
    yield assert_not_equal, b, newb
    yield assert_equal, type(dict()), type(newb)
    yield assert_equal, newb['a'], 3

def test_bunch_hash():
    # NOTE: Since the path to the json file is included in the Bunch,
    # the hash will be unique to each machine.
    pth = os.path.split(os.path.abspath(__file__))[0]
    json_pth = os.path.join(pth, 'realign_json.json')
    b = nii.Bunch(infile = json_pth, 
                  otherthing = 'blue',
                  yat = True)
    newbdict, bhash = b._get_bunch_hash()
    yield assert_equal, bhash, 'ddcc7b4ec5675df8cf317a48bd1857fa'
    # Make sure the hash stored in the json file for `infile` is correct.
    jshash = nii.md5()
    fp = file(json_pth)
    jshash.update(fp.read())
    fp.close()
    yield assert_equal, newbdict['infile'][0][1], jshash.hexdigest()
    yield assert_equal, newbdict['yat'], True


# create a temp file
#global tmp_infile, tmp_dir
#tmp_infile = None
#tmp_dir = None
def setup_file():
    #global tmp_infile, tmp_dir
    tmp_dir = tempfile.mkdtemp()
    tmp_infile = os.path.join(tmp_dir, 'foo.txt')
    open(tmp_infile, 'w').writelines('123456789')
    return tmp_infile

def teardown_file(tmp_dir):
    shutil.rmtree(tmp_dir)


def test_TraitedSpec():
    yield assert_true, nib.TraitedSpec().hashval
    yield assert_equal, nib.TraitedSpec().__repr__(), '\n\n'
    
    class spec(nib.TraitedSpec):
        foo = nib.traits.Int
        goo = nib.traits.Float(usedefault=True)


    yield assert_equal, spec().foo, Undefined
    yield assert_equal, spec().goo, 0.0
    specfunc = lambda x : spec(hoo=x)
    yield assert_raises, nib.traits.TraitError, specfunc, 1
    infields = spec(foo=1)
    hashval = ({'foo': 1, 'goo': 0.0}, '278d33656cb1a3e4f49b0e27dfed1b7d')
    yield assert_equal, infields.hashval[0], hashval[0]
    yield assert_equal, infields.hashval[1], hashval[1]
    yield assert_equal, infields.__repr__(), '\nfoo = 1\ngoo = 0.0\n'


def checknose():
    """check version of nose for known incompatability"""
    mod = __import__('nose')
    if mod.__versioninfo__[1] <= 11:
        return 0
    else:
        return 1

@skipif(checknose)
def test_TraitedSpec_withFile():
    tmp_infile = setup_file()
    tmpd, nme = os.path.split(tmp_infile)
    yield assert_true, os.path.exists(tmp_infile)
    class spec2(nib.TraitedSpec):
        moo = nib.File(exists=True)
        doo = nib.traits.List(nib.File(exists=True))
    infields = spec2(moo=tmp_infile,doo=[tmp_infile])
    yield assert_equal, infields.hashval[1], '8c227fb727c32e00cd816c31d8fea9b9'
    teardown_file(tmpd)
    
def test_Interface():
    yield assert_equal, nib.Interface.input_spec, None
    yield assert_equal, nib.Interface.output_spec, None
    yield assert_raises, NotImplementedError, nib.Interface
    yield assert_raises, NotImplementedError, nib.Interface.help
    yield assert_raises, NotImplementedError, nib.Interface._inputs_help
    yield assert_raises, NotImplementedError, nib.Interface._outputs_help
    yield assert_raises, NotImplementedError, nib.Interface._outputs

    class DerivedInterface(nib.Interface):
        def __init__(self):
            pass
        
    nif = DerivedInterface()
    yield assert_raises, NotImplementedError, nif.run
    yield assert_raises, NotImplementedError, nif.aggregate_outputs
    yield assert_raises, NotImplementedError, nif._list_outputs
    yield assert_raises, NotImplementedError, nif._get_filecopy_info

def test_BaseInterface():
    yield assert_equal, nib.BaseInterface.help(), None
    yield assert_equal, nib.BaseInterface._get_filecopy_info(), []
    

    class InputSpec(nib.TraitedSpec):
        foo = nib.traits.Int(desc='a random int')
        goo = nib.traits.Int(desc='a random int', mandatory=True)
        hoo = nib.traits.Int(desc='a random int', usedefault=True)
        zoo = nib.File(desc='a file', copyfile=False)
        woo = nib.File(desc='a file', copyfile=True)
    class OutputSpec(nib.TraitedSpec):
        foo = nib.traits.Int(desc='a random int')
    class DerivedInterface(nib.BaseInterface):
        input_spec = InputSpec
        
    yield assert_equal, DerivedInterface.help(), None
    yield assert_equal, DerivedInterface()._outputs(), None
    yield assert_equal, DerivedInterface._get_filecopy_info()[0]['key'], 'woo'
    yield assert_true, DerivedInterface._get_filecopy_info()[0]['copy']
    yield assert_equal, DerivedInterface._get_filecopy_info()[1]['key'], 'zoo'
    yield assert_false, DerivedInterface._get_filecopy_info()[1]['copy']
    yield assert_equal, DerivedInterface().inputs.foo, Undefined
    yield assert_raises, ValueError, DerivedInterface()._check_mandatory_inputs
    yield assert_equal, DerivedInterface(goo=1)._check_mandatory_inputs(), None
    yield assert_raises, ValueError, DerivedInterface().run
    yield assert_raises, NotImplementedError, DerivedInterface(goo=1).run
    
    class DerivedInterface2(DerivedInterface):
        output_spec = OutputSpec
        def _run_interface(self, runtime):
            return runtime

    yield assert_equal, DerivedInterface2.help(), None
    yield assert_equal, DerivedInterface2()._outputs().foo, Undefined
    yield assert_raises, Exception, DerivedInterface2(goo=1).run

    class DerivedInterface3(DerivedInterface2):
        def _run_interface(self, runtime):
            runtime.returncode = 1
            return runtime

    yield assert_equal, DerivedInterface3(goo=1).run().outputs, None
    
    class DerivedInterface4(DerivedInterface):
        def _run_interface(self, runtime):
            runtime.returncode = 0
            return runtime
    yield assert_equal, DerivedInterface4(goo=1).run().outputs, None

    class DerivedInterface5(DerivedInterface2):
        def _run_interface(self, runtime):
            runtime.returncode = 0
            return runtime
    yield assert_raises, NotImplementedError, DerivedInterface5(goo=1).run

    nib.BaseInterface.input_spec = None
    yield assert_raises, Exception, nib.BaseInterface

def test_Commandline():
    yield assert_raises, Exception, nib.CommandLine
    ci = nib.CommandLine(command='which')
    yield assert_equal, ci.cmd, 'which'
    yield assert_equal, ci.inputs.args, Undefined
    ci2 = nib.CommandLine(command='which', args='ls')
    yield assert_equal, ci2.cmdline, 'which ls'
    ci3 = nib.CommandLine(command='echo')
    ci3.inputs.environ = {'MYENV' : 'foo'}
    res = ci3.run()
    yield assert_equal, res.runtime.environ['MYENV'], 'foo'
    yield assert_equal, res.outputs, None

    class CommandLineInputSpec1(nib.TraitedSpec):
        foo = nib.traits.Str(argstr='%s', desc='a str')
        goo = nib.traits.Bool(argstr='-g', desc='a bool', position=0)
        hoo = nib.traits.List(argstr='-l %s', desc='a list')
        moo = nib.traits.List(argstr='-i %d...', desc='a repeated list',
                              position=-1)
        noo = nib.traits.Int(argstr='-x %d', desc='an int')
        roo = nib.traits.Str(desc='not on command line')
    nib.CommandLine.input_spec = CommandLineInputSpec1
    ci4 = nib.CommandLine(command='cmd')
    ci4.inputs.foo = 'foo'
    ci4.inputs.goo = True
    ci4.inputs.hoo = ['a', 'b']
    ci4.inputs.moo = [1, 2, 3]
    ci4.inputs.noo = 0
    ci4.inputs.roo = 'hello'
    cmd = ci4._parse_inputs()
    yield assert_equal, cmd[0], '-g'
    yield assert_equal, cmd[-1], '-i 1 -i 2 -i 3'
    yield assert_true, 'hello' not in ' '.join(cmd)
    
    class CommandLineInputSpec2(nib.TraitedSpec):
        foo = nib.File(argstr='%s', desc='a str', genfile=True)
    nib.CommandLine.input_spec = CommandLineInputSpec2
    ci5 = nib.CommandLine(command='cmd')
    yield assert_raises, NotImplementedError, ci5._parse_inputs

    class DerivedClass(nib.CommandLine):
        input_spec = CommandLineInputSpec2
        def _gen_filename(self, name):
            return 'filename'
        
    ci6 = DerivedClass(command='cmd')
    yield assert_equal, ci6._parse_inputs()[0], 'filename'
    nib.CommandLine.input_spec = nib.CommandLineInputSpec
