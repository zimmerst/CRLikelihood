'''
Created on Dec 15, 2015

@author: zimmer
'''
import sys, os
import numpy as np
from dsphs.pointlike.utilities import keyword_options, path
from dsphs.pointlike.Models import *
from skymaps import SkyDir,SkyIntegrator,Background,CompositeSkySpectrum,DiffuseFunction,EffectiveArea, \
                    Exposure,IsotropicSpectrum,IsotropicConstant,IsotropicPowerLaw 
import copy
import collections

class PointlikeException(Exception): pass

class PointSource(object):
    """ combine name, skydir, model """
    def __init__(self,skydir,name,model=None,free_parameters=True,leave_parameters=False):
        self.name    = name
        self.skydir = skydir
        self.model  = PowerLaw() if model is None else model
        #if not free_parameters:
        if not leave_parameters:
            for i in xrange(len(self.model.free)): self.model.free[i] = free_parameters
        self.duplicate = False
    def __str__(self):
        return '\n'.join(['\n',
                                '='*60,
                                'Name:\t\t%s'%(self.name),
                                'R.A. (J2000):\t%.5f'%(self.skydir.ra()),
                                'Dec. (J2000):\t%.5f'%(self.skydir.dec()),
                                'Model:\t\t%s'%(self.model.full_name()),
                                '\t'+self.model.__str__(indent='\t'), 
                                ])

    def copy(self):
        """ Create a deep copy of the point source. """
        return PointSource(SkyDir(self.skydir.ra(),self.skydir.dec()),
                           self.name,self.model.copy(),leave_parameters=True)

class DiffuseSource(object):
    """ Associate a spatial model with a spectral scaling model."""
    __counter = 0

    def __init__(self,diffuse_model,scaling_model,name=None):

        self.dmodel = diffuse_model
        self.smodel = scaling_model

        self.smodel.background = True

        if name is None:
            self.name = 'Diffuse Source %d'%(DiffuseSource.__counter)
            DiffuseSource.__counter += 1
        else: self.name = name

        if not isinstance(self.dmodel,collections.Iterable):
            self.dmodel = [self.dmodel]
   
    def __str__(self): return '\n'.join((self.name,'\t'+self.dmodel.__str__(),
            '\t'+self.smodel.__str__()))

    def __getstate__(self):
        """ this is a poor man's pickeling. Convert all the IsotropicSpectrum and DiffuseFunction
            objects into pickelable filenames. This implementation should be more robust. 
            
            You should be able to pickle Diffuse sources:

                >>> import pickle
                >>> import os
                >>> from uw.utilities import path
                >>> from uw.like.Models import Constant, FileFunction
                >>> from skymaps import IsotropicPowerLaw

            Pickle isotropic:

                >>> iso = path.expand('$GLAST_EXT/diffuseModels/v2r0p1/isotrop_2year_P76_source_v1.txt')
                >>> ds = DiffuseSource(IsotropicSpectrum(iso),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))

            After pickling the object, the original object should be uncahgned:

                >>> print type(ds.dmodel[0])
                <class 'skymaps.IsotropicSpectrum'>
            
            Pickle galactic diffuse:

                >>> gal = path.expand('$GLAST_EXT/diffuseModels/v2r0p1/ring_2year_P76_v0.fits')
                >>> ds = DiffuseSource(DiffuseFunction(gal),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))
            
            Pickle Isotropic PowerLaw:

                >>> ds = DiffuseSource(IsotropicPowerLaw(),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))

            Pickle Isotropic Constant:

                >>> ds = DiffuseSource(IsotropicConstant(),Constant())
                >>> pickle.dump(ds, open(os.devnull,'w'))
            
        """
        d=copy.copy(self.__dict__)

        def convert_spectrum(spectrum):
            if isinstance(spectrum,IsotropicSpectrum) or isinstance(spectrum,DiffuseFunction):
                return (type(spectrum),spectrum.name())
            elif isinstance(spectrum,IsotropicPowerLaw):
                return (type(spectrum),spectrum.flux(),spectrum.index())
            elif isinstance(spectrum,IsotropicConstant):
                return (type(spectrum),spectrum.constant())
            else:
                # unrecognized type
                return spectrum

        d['dmodel'] = map(convert_spectrum,d['dmodel'])
        return d

    def __setstate__(self,state):
        """ recreate the dmodel. """
        self.__dict__ = state
        def unconvert_spectrum(spectrum):
            if type(spectrum) == tuple:
                return spectrum[0](*spectrum[1:])
            else:
                return spectrum

        self.dmodel = map(unconvert_spectrum,self.dmodel)

    def copy(self):
        """ Make a copy of a diffuse source. 
        
            First, create the DS

                >>> from uw.like.Models import Constant
                >>> from skymaps import IsotropicPowerLaw
                >>> ds = DiffuseSource(IsotropicConstant(),Constant())

            Nowe, we can copy it:

                >>> ds_copy = ds.copy()
                >>> print type(ds.dmodel[0])
                <class 'skymaps.IsotropicConstant'>
                >>> print type(ds_copy.dmodel[0])
                <class 'skymaps.IsotropicConstant'>
        """
        return DiffuseSource(
            name = self.name,
            diffuse_model = self.dmodel,
            scaling_model = self.smodel.copy())

class ExtendedSource(DiffuseSource):
    """ Class inherting from DiffuseSource but implementing a spatial source. 
        The main difference is the requirement of a spatial model to accomany 
        a spectral model. """

    defaults = (
        ('name',None,'The name of the extended source.'),
        ('model',None,'a Model object.'),
        ('spatial_model',None,"""The spatial model to use. 
                                       This is a SpatialModel object."""),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,**kwargs):
        """ Make the naming consistent with the PointSource object so that
            extended sources 'feel' like point sources.

            """
        keyword_options.process(self, kwargs)

        if self.model == None: self.model = PowerLaw()
        if self.spatial_model == None: self.spatial_model = Disk()

        if not isinstance(self.spatial_model,SpatialModel):
            raise Exception("The diffuse_model passed to an Extended Source must inherit from SpatialModel.")

        super(ExtendedSource,self).__init__(
            diffuse_model = self.spatial_model,
            scaling_model = self.model,
            name          = self.name)

        self.model.background = False

    @property
    def skydir(self): return self.spatial_model.center

    @property
    def smodel(self): 
        """ No reason to keep a model & smodel. """
        return self.model

    @smodel.setter
    def smodel(self, value): self.model = value


    def __str__(self,indent=''):
        return indent+('\n'+indent).join(['\n',
                          '='*60,
                          'Name:\t\t%s'%(self.name),
                          'R.A. (J2000):\t\t%.5f'%(self.spatial_model.center.ra()),
                          'Dec. (J2000):\t\t%.5f'%(self.spatial_model.center.dec()),
                          'Model:\t\t%s'%(self.model.full_name()),
                          '\t'+self.model.__str__(indent='\t'), 
                          'SpatialModel:\t%s'%(self.spatial_model.full_name()),
                          '\t'+self.spatial_model.__str__(indent='\t')
                         ])

    def copy(self):
        """ Create a deep copy of an extended source. """
        return ExtendedSource(name=self.name,
                              spatial_model=self.spatial_model.copy(),
                              model=self.model.copy())


###=========================================================================###

class Singleton(object):
    """Implement a singleton class."""

    __instances = {}

    def __init__(self,constructor,key=None,*args,**kwargs):
        """Constructor -- a class name; note -- not a string!
            Key            -- an optional key; constructor name used otherwise.
                                **keys are unique, not constructors!**
            Additional args and kwargs are passed to constructor."""

        self.add(constructor,key,*args,**kwargs)

    def add(self,constructor,key=None,*args,**kwargs):
        inst = Singleton._Singleton__instances
        key  = str(constructor) if key is None else key
        if key not in inst.keys(): inst[key] = constructor(*args,**kwargs)
        return inst[key]

    def __call__(self,key):
        inst = Singleton._Singleton__instances
        key  = str(key)
        if key in inst: return inst[key]

class Singleton2(Singleton):

    def __init__(self):
        pass

def get_diffuse_source(spatialModel='ConstantValue',
                              spatialModelFile=None,
                              spectralModel='PowerLaw',
                              spectralModelFile=None,
                              name=None,
                              diffdir = None):

    """ Return a DiffuseSource instance suitable for
         instantiating a child of ROIDiffuseModel.

         NB -- don't support front/back distinction atm.

         The list of supported models is currently very short, but covers
         the usual cases for modeling diffuse backgrounds.  Additional
         use cases can be developed on an ad hoc basis.
         
         Arguments:
         
         spatialModel -- an XML-style keyword.  Valid options are
                               1) ConstantValue (isotropic)
                               2) MapCubeFunction (from a FITS file)
                               
         spatialModelFile -- if a mapcube is specified, its location
         
         spectralModel -- This can be either an XML-style keyword or an
                                instance of Model.
                                If an XML-style keyword, valid options are
                                1) FileFunction
                                2) PowerLaw
                                3) Constant
                               
         spectralModelFile -- if a tabular function is specified,
                                     its location

         name -- a name for the ol' model

         
         diffdir -- if the XML files specify paths relative to some
                        directory, set this variable appropriately
    """

    if (diffdir is not None):
        if spatialModelFile is not None:
            spatialModelFile = os.path.join(diffdir,spatialModelFile)
        if spectralModelFile is not None:
            spectralModelFile = os.path.join(diffdir,spectralModelFile)

    # check input sanity
    if not isinstance(spectralModel,Model):
        if (spectralModelFile is not None):
            if not os.path.exists(path.expand(spectralModelFile)):
                raise Exception('Could not find the ASCII file specified for FileFunction')
        elif not (spectralModel == 'PowerLaw' or spectralModel == 'Constant'):
            raise NotImplementedError,'Must provide one of the understood spectral models.'
        else:
            pass

    if spatialModel=='MapCubeFunction':
        if (spatialModelFile is None) or (not os.path.exists(path.expand(spatialModelFile))):
            raise Exception('Could not find the FITS file specified for MapCubeFunction (file = %s).' % spatialModelFile)
    elif spatialModel != 'ConstantValue':
        raise NotImplementedError,'Must provide one of the understood spatial models.'
    else:
        pass                  

    ston = Singleton2()
    dmodel = None; smodel = None


    # deal with isotropic models
    if spatialModel=='ConstantValue':
        if isinstance(spectralModel,Model):
            smodel=spectralModel
            dmodel=IsotropicConstant()
        elif spectralModelFile is not None:
            smodel = FileFunction(normalization=1, file=spectralModelFile)
            dmodel = IsotropicConstant()
        elif spectralModel == 'PowerLaw':
            # use Sreekumar-like defaults
            smodel = PowerLaw(index=2.1)
            smodel.set_flux(1.5e-5, emin=100, emax=N.inf)

            dmodel = IsotropicConstant()
        else:
            raise Exception("Unable to parse input.")

    # deal with mapcubes
    else:
        if spectralModel == 'FileFunction':
            dmodel1 = IsotropicSpectrum(spectralModelFile)
            dmodel2 = ston.add(DiffuseFunction,spatialModelFile,spatialModelFile)
            dmodel  = CompositeSkySpectrum(dmodel1,dmodel2)
            dmodel.saveme1 = dmodel1; dmodel.saveme2 = dmodel2
            smodel  = Constant()
        else:
            dmodel = ston.add(DiffuseFunction,path.expand(spatialModelFile),path.expand(spatialModelFile))
            dmodel.filename=spatialModelFile
            if spectralModel == 'PowerLaw':
                smodel = ScalingPowerLaw()
            elif spectralModel == 'Constant':
                smodel = Constant()
            else:
                smodel = spectralModel

    if (dmodel is None) or (smodel is None):
         raise Exception('Was unable to parse input.')

    return DiffuseSource(dmodel,smodel,name)


