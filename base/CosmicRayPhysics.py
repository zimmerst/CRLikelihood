'''
Created on Dec 16, 2015

@author: zimmer
@brief:  base class for CR model definitions. 
'''
from __builtin__ import None
from uw.like.Models import *
from uw.like.SpatialModels import *

class CRModel(object):
    '''
    Note a CR model is defined by a spectrum, a spatial model and a unique name.
    '''

    def __init__(self, params):
        '''
        Constructor
        '''
        self.spectrum = None
        self.spatialModel = None
        self.name = None
        
    def __repr__(self):
        ''' 
        the string representation of the object
        '''
        return self.name
    
    def quickConvolution(self):
        '''
        performs a quick gtsrcmap call on the model
        '''
        pass
    
    def findRadius(self,scaled_srcmap,algorithm='delta'):
        '''
        depending on the algorithm, return the disk radius that matches best
        '''
        pass

def init_models(configfile):
    '''
    the work horse - here all models are instantiated and returned as list
    '''
    pass