'''
Created on Dec 16, 2015

@author: zimmer
@brief:  base class for CR model definitions. 
'''
import subprocess, os, pyfits
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.stats.morestats import pdf_fromgamma

from scipy.interpolate import interp1d
import numpy as np

try:
    from skymaps import SkyDir
except ImportError:
    print 'Could not find skymaps, trying FSSC init'
    from pyLikelihood import SkyDir

from dsphs.pointlike.Models import FileFunction

from uw.like.SpatialModels import InterpProfile #FIXME: remove dependency!

from dsphs.utils.par2cmd import par2cmd
from dsphs.utils.tools import yaml_load

def array2DtoProfile(array2D,pixelsize=1.,interpolate=False):
    # get dimensions
    n, m = np.shape(array2D)
    x = np.arange(0,n/2.)
    y = array2D[m/2,n/2:]
    x*= np.abs(pixelsize)
    return (x,y if not interpolate else interp1d(x,y))

class TabulatedProfile(InterpProfile):
    ''' 
    to be compatible with Fabio's files
    '''
    def __init__(self,datfile,center=None,kind='cubic'):
        r_in_degrees, pdf = np.loadtxt(datfile,unpack=True)
        kwargs = {}
        kwargs['r_in_degrees']=r_in_degrees 
        kwargs['pdf']=pdf
        kwargs['kind']=kind
        kwargs['center']=center
        super(TabulatedProfile,self).__init__(**kwargs)
    
class CRModel(object):
    '''
    Note a CR model is defined by a spectrum, a spatial model and a unique name.
    '''

    def __init__(self, params,datadir=None,**kwargs):
        '''
        Constructor
        '''
        self.datadir            = datadir
        self.spectrum           = FileFunction(normalization=1, file=os.path.join(datadir,params['spectrum']))
        self.center             = SkyDir(194.9531, 27.9807,SkyDir.EQUATORIAL)
        self.spatialModel       = TabulatedProfile(os.path.join(datadir,params['profile']),center=self.center)
        self.convolvedTemplate  = None
        if 'name' in params:
            self.name               = params['name']
        self.__dict__.update(kwargs)

    def __repr__(self):
        ''' 
        the string representation of the object
        '''
        return self.name
        
    def quickConvolution(self,parfile,verbose=True,cleanup=True):
        '''
        performs a quick gtsrcmap call on the model
        '''
        if not self.convolvedTemplate is None:
            print 'Warning, convolved template attribute already set, will overwrite model definition.'
        templateName = "%s_spatial.fits"%self.name
        self.spatialModel.save_template(templateName,npix=500)
        os.environ['SRC_TEMPLATE']=os.path.abspath(templateName)
        cmd,parDict = par2cmd(parfile,returnDict=True)
        # replace the datadir
        while "$(DATADIR)" in cmd:
            cmd = cmd.replace("$(DATADIR)",os.path.abspath(os.getenv("DATADIR")))
        print cmd
        pts = subprocess.Popen(cmd,shell=True)
        pts.communicate()
        if cleanup:
            os.remove(templateName)
            # add more cleanup files.
        # need to cleanup stuff...
        ofile = parDict['outfile']
        while '"' in ofile:
            ofile = ofile.replace('"',"")
        self.convolvedTemplate = os.path.abspath(ofile)
        return self.convolvedTemplate
    
    def findRadius(self,scaled_srcmap,algorithm='delta_average'):
        '''
        depending on the algorithm, return the disk radius that matches best
        '''
        print 'INFO: Algorithm: %s'%algorithm
        def val2str(v): return ("%1.1f"%v).replace(".","d")
        def str2val(s): return float(s.replace("d","."))
        
        def average(a,ebounds=None): 
            return np.average(a,axis=0)
            
        def weighted_energy(a,ebounds,index=-2.0):
            def scaled_slice(slice,e,e0):
                return slice*np.power(e/e0,index)
            e_min, e_max = ebounds['E_MIN'],ebounds['E_MAX']
            energies = np.power(10,np.log10(e_min)+(np.log10(e_max)-np.log10(e_min))/2.)
            e0 = e_min[0]
            Z_a = np.sum(np.array([scaled_slice(a[i],e,e0) for i,e in enumerate(e_min)]),axis=0)
            return Z_a
        # implement other algorithms as we need.
        def weighted_energy_spectra(a,ebounds):
            def scaled_slice(a,e,e0):
                return a*self.spectrum.i_flux(e0,e,e_weight=0,cgs=False)#add integral here.
            e_min, e_max = ebounds['E_MIN'],ebounds['E_MAX']
            energies = np.power(10,np.log10(e_min)+(np.log10(e_max)-np.log10(e_min))/2.)
            e0 = e_min[0]
            Z_a = np.sum(np.array([scaled_slice(a[i],e,e0) for i,e in enumerate(e_min)]),axis=0)
            return Z_a
          
        if self.convolvedTemplate is None:
            raise Exception("can't find convolved template. Have you run the convolution yet?")
        print 'loading convolved template %s'%self.convolvedTemplate
        target, tHeader = pyfits.getdata(self.convolvedTemplate,"target",header=True)
        ebounds = pyfits.getdata(self.convolvedTemplate,"ebounds")
        target2D = eval(algorithm)(np.array(target,dtype=float),ebounds=ebounds)
        target_profile = array2DtoProfile(target2D, pixelsize=np.abs(tHeader['CDELT1']), interpolate=True)
        target_profile_interp = target_profile[1]
        matching = None
        disk_radii = np.arange(0.1,1.1,0.1)
        foms = np.nan*np.ones( len(disk_radii) )
        for i in range(0,10):
            val = (.1+float(i)/10.)
            reference, rHeader = pyfits.getdata(scaled_srcmap,val2str(val),header=True)
            ref2D = eval(algorithm)(np.array(reference,dtype=float),ebounds=ebounds)
            ref_profile = array2DtoProfile(ref2D, pixelsize=np.abs(rHeader['CDELT1']), interpolate=True)
            ref_profile_interp = ref_profile[1]
            x = target_profile[0]
            fractional_difference = np.abs(ref_profile_interp(x) - target_profile_interp(x))/ref_profile_interp(x)
            fom = np.sum(fractional_difference)
            foms[i] = fom
            print '*DEBUG* val, rel difference (sum): ',val, fom
            row = np.array([val,fom])
            if matching is None:
                matching = row
            else:
                matching = np.vstack((matching,row))
        # best key
        # minimize the fractional difference 
        return val2str(matching.T[0][np.argmin(matching.T[1])]),(disk_radii,foms)
    
def init_models(configfile,datadir=None):
    '''
    the work horse - here all models are instantiated and returned as list
    '''
    modelfile = yaml_load(configfile['models'])
    models = [CRModel(modelfile[m],name=m,datadir=datadir) for m in modelfile]
    return models