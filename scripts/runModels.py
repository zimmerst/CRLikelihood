'''
Created on Dec 15, 2015

@brief: will run likelihood calculations on a set of models.
@author: zimmer
'''

#!/usr/bin/env python
from optparse import OptionParser
import yaml, copy, os 
from os.path import join
import numpy as np

from dsphs.base.CosmicRayPhysics    import CRModel, init_models
from dsphs.base.sed                 import SED
from dsphs.base.lnlfn               import ProfileLnL, ProfileLimit, LnLFn
from dsphs.utils.tools              import yaml_dump, update_dict

usage = "Usage: %prog  [options] config.yaml"
description = "run all model comparisons"
parser = OptionParser(usage=usage,description=description)
parser.add_option("--text_output",action='store_true',default=False)
parser.add_option("--seddir",default=None)
(opts, args) = parser.parse_args()

config = {'jsigma': None}
update_dict(config,yaml.load(open(args[0])),add_new_keys=True)

datadir = os.getenv("DATADIR","data")
print 'initializing CR models...'
models_to_test = init_models(config,datadir) # this should load and init all models
print 'found %i models'%len(models_to_test)

# r holds all results.
r = dict(
        models = copy.copy([str(m) for m in models_to_test]), # note: we are only interested in the names, list contains full obj. reference
        ulimits68 = np.nan*np.ones( len(models_to_test) ),  # upper limits w/o nuisance J-factor
        fluxes68 = np.nan*np.ones( len(models_to_test) ),   # flux upper limits w/o nuisance J-factor    
        ulimits95 = np.nan*np.ones( len(models_to_test) ),  # upper limits w/o nuisance J-factor
        fluxes95 = np.nan*np.ones( len(models_to_test) ),   # flux upper limits w/o nuisance J-factor    
        ulimits99 = np.nan*np.ones( len(models_to_test) ),  # 99% upper limits w/o nuisance J-factor
        fluxes99 = np.nan*np.ones( len(models_to_test) ),   # 99% flux upper limits w/o nuisance J-factor
        mle = np.nan*np.ones( len(models_to_test) ),      # MLE w/o nuisance J-factor
        ts = np.nan*np.ones( len(models_to_test) ),       # TS value, only for curiosity!
        )
data = {}

for i,model in enumerate(models_to_test):
    if not isinstance(model,CRModel): 
        raise Exception("Model is not a CRModel instance, giving up!")
    print 'working on model %s'%str(model)
    # find associated disk
    print 'convolve template with PSF'
    model.quickConvolution(config['parfile'],verbose=True,cleanup=False)
    print 'find matching disk'
    matching_disk = model.findRadius(join(datadir,config['scaled_disk_srcmap']),algorithm='delta_average')
    if matching_disk is None:
        print "Error: could not find associated radius with profile %s"%str(model)
    print 'best matching radius : %s'%matching_disk
    # load matching SED file
    sedfile = join(datadir,config['sed'][matching_disk])
    print 'Loading SED file:', sedfile
    sed = SED.read_yaml(sedfile)
    # global logLike
    print 'get global logLike given the spectral form'
    pd = sed.get_global_logLike(model.spectrum, config['jsigma'])
    (norm, flux, lnl, p1lnl, flnl) = pd[0]
    p2lnl = LnLFn(pd[1][0],pd[1][3])(norm)
    lnlx = norm; lnly = lnl
    print 'Calculating limits.'
    try:
        lnlfn = LnLFn(lnlx,lnly)
        p1lnlfn = LnLFn(lnlx,p1lnl)
        p2lnlfn = LnLFn(lnlx,p2lnl)
        r['models'][i]=str(model)
        r['mle'][i] = lnlfn.mle()
        r['ulimits68'][i]  = ProfileLimit( lnlx, lnly).getLimit( 0.32 )        
        r['fluxes86'][i]  = ProfileLimit( flux, flnl).getLimit( 0.32 )
        r['ulimits95'][i]  = ProfileLimit( lnlx, lnly).getLimit( 0.05 )        
        r['fluxes95'][i]  = ProfileLimit( flux, flnl).getLimit( 0.05 )
        r['ulimits99'][i]  = ProfileLimit( lnlx, lnly).getLimit( 0.01 )
        r['fluxes99'][i]  = ProfileLimit( flux, flnl).getLimit( 0.01 )
        r['ts'][i] = float(2*(p1lnlfn(p1lnlfn.mle()) - p1lnlfn(0)))
    except Exception, message:
        print message
        continue
    print "Writing output",config['outfile']
    yaml_dump(r,config['outfile'],update=True) # keep track of things!
print 'Done.'
