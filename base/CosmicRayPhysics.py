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

try:
    from skymaps import SkyDir
except ImportError:
    print 'Could not find skymaps, trying FSSC init'
    from pyLikelihood import SkyDir

from dsphs.pointlike.Models import FileFunction

from uw.like.SpatialModels import InterpProfile #FIXME: remove dependency!

from dsphs.utils.par2cmd import par2cmd
from dsphs.utils.tools import yaml_load

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
        def val2str(val): return ("%1.1f"%val).replace(".","d")
            
        def delta_average(a,b,ebounds=None):
            avg_a = np.average(a,axis=0)
            avg_b = np.average(b,axis=0)
            return np.abs(avg_a - avg_b)
        def delta_weighted_energy(a,b,ebounds,index=-2.0):
            def scaled_slice(slice,e,e0):
                return slice*np.power(e/e0,index)
            e_min, e_max = ebounds['E_MIN'],ebounds['E_MAX']
            energies = np.power(10,np.log10(e_min)+(np.log10(e_max)-np.log10(e_min))/2.)
            e0 = e_min[0]
            Z_a = np.sum(np.array([scaled_slice(a[i],e,e0) for i,e in enumerate(e_min)]))
            Z_b = np.sum(np.array([scaled_slice(b[i],e,e0) for i,e in enumerate(e_min)]))
            return np.abs(Z_a - Z_b)
        # implement other algorithms as we need.
        
        if self.convolvedTemplate is None:
            raise Exception("can't find convolved template. Have you run the convolution yet?")
        print 'loading convolved template %s'%self.convolvedTemplate
        target = np.array(pyfits.getdata(self.convolvedTemplate,"target"),dtype=float)
        ebounds = pyfits.getdata(self.convolvedTemplate,"ebounds")
        matching = None
        for i in range(0,10):
            val = (.1+float(i)/10.)
            reference = np.array(pyfits.getdata(scaled_srcmap,val2str(val)),dtype=float)
            fom = np.sum(eval(algorithm)(target,reference,ebounds=ebounds)) # to convert a map to one single value
            print '*DEBUG* val, fom: ',val, fom
            row = np.array([val,fom])
            if matching is None:
                matching = row
            else:
                matching = np.vstack((matching,row))
        # best key
        return val2str(matching.T[0][np.argmin(matching.T[1])])
    
def init_models(configfile,datadir=None):
    '''
    the work horse - here all models are instantiated and returned as list
    '''
    modelfile = yaml_load(configfile['models'])
    models = [CRModel(modelfile[m],name=m,datadir=datadir) for m in modelfile]
    return models