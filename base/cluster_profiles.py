'''
Created on Apr 7, 2015

@author: zimmer
'''

import numpy as np
# astropy -> deals with units :)
from astropy import units as u
# some science tools stuff.
from dsphs.pointlike.SpatialModels import InterpProfile
SMALL_ANALYTIC_EXTENSION=1e-10
SMALL_NUMERIC_EXTENSION=1e-3

def Mpc2Deg(r_mpc,dlum=100):
    # Coma is 100 Mpc away!
    return np.arctan(r_mpc/dlum)*(180./np.pi)*u.deg

class CoredProfile(InterpProfile):
    # these are bogus parameters!
    #default_p = [1.0,0.9,0.1]
    #param_names = ['sigma','index','core']
    #default_limits = [[SMALL_ANALYTIC_EXTENSION,10.],[SMALL_NUMERIC_EXTENSION,10.],[SMALL_ANALYTIC_EXTENSION,10.]]
    #log = [False,False,False]
    #steps = [0.04,0.01,0.04]

    def __init__(self,rvir=None,core=None,index=None,center=None,kind='cubic'):
        def I_r(r,p=.9,theta_c=0.1):
            a = 1+(r/theta_c)**2
            return np.power(a,-p)
        radii = np.power(10,np.linspace(-2, np.log10(1.5),1000)) # log sampling
        allowed_radii = radii[np.where(radii<=rvir)]
        kwargs = {}
        kwargs['r_in_degrees']=allowed_radii
        kwargs['pdf']=I_r(allowed_radii,p=index,theta_c=core)
        kwargs['kind']=kind
        kwargs['center']=center
        super(CoredProfile,self).__init__(**kwargs)
                
    def setp(self,key,value):
        if key in self.__dict__:
            self.__dict__[key]=value

def getRadialProfileFromSquaredMap(data,radial_bins=None,errors=False):
    ''' 
        computes a radial profile from a squared map 
        data : array n x n containing the data to plot
        returns: array with radial bins, profile and optional errors (1sigma)
        radial_bins : provide custom binning
        errors : if set to True, returns errors on top of values alone.
    '''
    # FIXME: handle even pixels #
    dimensions = np.shape(data)
    if dimensions[0]!=dimensions[1]:
        raise Exception("could not broadcast shapes, expect squared data")
    if not (dimensions[0] % 2 != 0):
        raise NotImplementedError("can only handle odd data so far!")
    center = tuple([val/2 for val in dimensions])
    x0, y0 = center
    r_px = lambda x,y : np.sqrt((x-x0)**2+(y-y0)**2)
    r_arr = np.zeros(np.shape(data))
    for i,x in enumerate(data[0][:]):
        for j,y in enumerate(data[:][0]):
            r_arr[i][j]=r_px(i,j)
            
    if radial_bins is None:
        radial_bins = np.linspace(0,np.max(r_arr))
    
    r_prof = np.zeros(np.shape(radial_bins))
    r_prof_err = np.zeros(np.shape(radial_bins))
    for i,radius in enumerate(radial_bins):
        mask = None
        if i<len(radial_bins)-1:
            try:
                mask = np.where(r_arr>radius) and np.where(r_arr<=radial_bins[i+1])
            except IndexError,IE:
                print 'shape of radius',np.shape(radial_bins),' index failure: ',i
                raise IndexError(IE)
        else:
            mask = np.where(r_arr>radius)
        r_prof[i]=np.average(data[mask])
        r_prof_err[i]=np.std(data[mask])#np.sqrt(np.sum(data[mask]**2))
    ret = [radial_bins,r_prof]
    if errors:
        ret.append(r_prof_err)
    return ret
