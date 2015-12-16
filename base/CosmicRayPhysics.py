'''
Created on Dec 16, 2015

@author: zimmer
@brief:  base class for CR model definitions. 
'''
import subprocess, os, pyfits

from skymaps import SkyDir
from uw.like.Models import FileFunction
from uw.like.SpatialModels import InterpProfile

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

    def __init__(self, params):
        '''
        Constructor
        '''
        self.spectrum           = FileFunction(normalization=1, file=params['spectrum'])
        self.center             = SkyDir(194.9531, 27.9807,SkyDir.EQUATORIAL)
        self.spatialModel       = TabulatedProfile(params['profile'],center=self.center)
        self.convolvedTemplate  = None
        self.name               = params['name']
        
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
        print cmd
        pts = subprocess.Popen(cmd,shell=True)
        pts.communicate()
        if cleanup:
            os.remove(templateName)
            # add more cleanup files.
        self.convolvedTemplate = parDict['outfile']
        return self.convolvedTemplate
    
    def findRadius(self,scaled_srcmap,algorithm='delta_average'):
        '''
        depending on the algorithm, return the disk radius that matches best
        '''
        def val2str(val): return ("%1.1f"%val).replace(".","d")
            
        def delta_average(a,b):
            avg_a = np.average(a,axis=0)
            avg_b = np.average(b,axis=0)
            return np.abs(avg_a - avg_b)
        # implement other algorithms as we need.
        
        if self.convolvedTemplate is None:
            raise Exception("can't find convolved template. Have you run the convolution yet?")
        target = np.array(pyfits.getdata(self.convolvedTemplate,"target"),dtype=float)
        matching = None
        for i in range(0,10):
            val = (.1+float(i)/10.)
            reference = np.array(pyfits.getdata(scaled_srcmap,val2str(val)),dtype=float)
            fom = algorithm(target,reference)
            row = np.array([val,fom])
            if matching is None:
                matching = row
            else:
                matching = np.vstack((matching,row))
        # best key
        return val2str(matching.T[0][np.argmin(matching.T[1])])
    
def init_models(configfile):
    '''
    the work horse - here all models are instantiated and returned as list
    '''
    modelfile = yaml_load(configfile['models'])
    models = [CRModel(modelfile[m],name=m) for m in modelfile]
    return models