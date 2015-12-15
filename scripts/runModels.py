'''
Created on Dec 15, 2015

@brief: will run likelihood calculations on a set of models.
@author: zimmer
'''

#!/usr/bin/env python
import os, tempfile, glob, shutil, time, copy, sys
from os.path import join, splitext, basename, dirname
import numpy as np
import yaml
import pyfits
from scipy.optimize import fmin, brentq
from tempfile import NamedTemporaryFile

from dsphs.base.target import load_targets_from_config
import dsphs.utils.tools as tools
from dsphs.utils.tools import mplconfig, pwd, mkscratch, mkdir
from uw.like.Models import FileFunction, PowerLaw

mplconfig()
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')
from dsphs.base.sed import SED
from dsphs.base.lnlfn import ProfileLnL, ProfileLimit, LnLFn
from LikelihoodState import LikelihoodState

datadir = os.getenv("DATADIR","data")
diskFile = "srcmap_scaled_disk.fits"

workdir=join(basedir,"analysis/%s/%s/"%(target,opts.spectrum))

jsigma = target.get_safe_jsigma(profile) # this is a friggin' remnant...

matching_disk = None
# PSEUDOCODE!
models_to_test = init_models() # this should load and init all models
for model in models_to_test:
    # find associated disk
    model.quickConvolution(verbose=True)
    matching_disk = model.findRadius(os.path.join(datadir,diskFile))
    if matching_disk is None:
        print "Error: could not find associated radius with profile %s"%str(model)
    # load matching SED file
    sedfile = os.path.join(datadir,"sed_%s"%matching_disk)
    print 'Loading ', sedfile
    sed = SED.read_yaml(sedfile)
    # calculate flux Limit
    pd = sed.get_global_logLike(model.spectrum, jsigma)
    (norm, flux, lnl, p1lnl, flnl) = pd[0]
    p2lnl = LnLFn(pd[1][0],pd[1][3])(norm)
    m1lnl = np.zeros(np.shape(p2lnl))
    profile_data = dict(norm=norm,flux=flux,lnl=lnl,p1lnl=p1lnl,
                        flnl=flnl,p2lnl=p2lnl,m1lnl=m1lnl)
    data = np.vstack( (norm, flux, lnl, p1lnl, flnl, p2lnl, m1lnl) ).T
    lnlfile = join(workdir,'%s_%s.yaml'%(target,str(model)))
    print "Writing %s"%lnlfile  
    tools.yaml_dump(dict(lnldata = profile_data,param=spparam,
                         jvalue=jvalue,jsigma=jsigma,
                         jnominal=jnominal,rjvalues=rjvalues),lnlfile)