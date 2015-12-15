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
from dsphs.pointlike.Models import FileFunction, PowerLaw

mplconfig()
import matplotlib
try:             os.environ['DISPLAY']
except KeyError: matplotlib.use('Agg')

from dsphs.base.sed import SED

from dsphs.base.lnlfn import ProfileLnL, ProfileLimit, LnLFn

from LikelihoodState import LikelihoodState

def FileSpectrum(infile,jvalue,z=0,sleep=None):
    # Input of the data file format is:
    # Column 1: Energy [GeV]
    # Column 2: ~Differential Flux...
    #           (sigmav/(8*pi*mass^2)) * Sum( sigmav_i/sigmav * dN/dE ) * 1e29
    #           [cm^3 s^-1 / GeV^2 ] * [ photons/GeV ]
    data = np.loadtxt(infile,unpack=True)
    # Deal with files containing NaN...
    if np.isnan(np.sum(data[1])):
        print "WARNING: NaN found... %s"%basename(infile)
        return None
    data = data.compress( data[1] > 0, axis=1) # Clip zeros
    # Convert energy (MeV)
    energy = data[0] * 1e3 
    energy*=np.power(1+z,-1)

    # Convert to dN/dE (ph/MeV/s/cm^2)
    dnde = data[1] * 1e-29 * 1e-3 * jvalue 
    temp = NamedTemporaryFile(dir=mkscratch())
    filename = temp.name
    np.savetxt(filename, np.array((energy,dnde)).T, fmt="%-12g")
    spec = FileFunction(file=filename)
    if sleep is not None: tools.sleep(sleep)
    return spec

from optparse import OptionParser
usage = "Usage: %prog  [options] config.yaml"
description = "Full instance of dwarf analysis"
parser = OptionParser(usage=usage,description=description)
parser.add_option("-d","--dwarf",default="draco")
parser.add_option("-t","--target",default="draco")
parser.add_option("-s","--spectrum",default='pmssm')
parser.add_option("-p","--param", default=0, type='float')
parser.add_option("--mc",action='store_true')
parser.add_option("--plot",action='store_true')
parser.add_option("--text_output",action='store_true',default=False)
parser.add_option("--seddir",default=None)
(opts, args) = parser.parse_args()

config = { 'profile' : None, 'seed' : 0, 'randomize' : {'jvalues' : False} ,'decay':False}
tools.update_dict(config,yaml.load(open(args[0])),add_new_keys=True)
basedir = pwd()

profile = config['profile'] #if 'profile' in config else None

# Parse the SED file
sed_index = config.get('sed_index',2.0)
sedfile = join(basedir,"analysis/%s/pwl/%s/%s_sed.yaml"%(opts.target,sed_index,opts.target))

# logLike2 -- standard likelihood
# logLike3 -- likelihood with user-defined SED input file

sedfiles = [sedfile]
suffixes = ['logLike2']
decay = False
if config['decay']:
    print '*DECAY mode requested*'
    decay = True
if opts.seddir is not None:
    sedfiles += [join(opts.seddir,"analysis/%s/pwl/%s/%s_sed.yaml"%
                      (opts.target,sed_index,opts.target))]
    suffixes += ['logLike3']
    
for sedfile, suffix in zip(sedfiles,suffixes):

    TARGETS = load_targets_from_config(config)
    jnominal = TARGETS[opts.target].get_safe_jvalue(profile)
    target = TARGETS[opts.target]

    workdir=join(basedir,"analysis/%s/%s/"%(target,opts.spectrum))
    
    # Form of the likelihood function
    fntypes = ['lgauss','lgauss_like','lgauss']

    # Methods
    methods = ['profile','profile','marginalize']

    # Apply correction to J factor
    rescale = [True, False, False]
    
    # Create an array of target sets with randomized J values
    rtargets = []

    # Generate target sets with randomized J values
    for t, r in zip(fntypes,rescale):
        if config['randomize']['jvalue']:
            rtargets.append(TARGETS.randomize_jvalues(profile=profile,
                                                      seed=config['seed'],
                                                      fntype=t,rescale=r))
        else:
            rtargets.append(TARGETS)

    target = TARGETS[opts.target]
    z = target.z # should be defined!
    z = 0 ## override z=0 for all objects
    jvalue = target.get_safe_jvalue(profile)
    jsigma = target.get_safe_jsigma(profile)
    print "J-Value: %.3f +/- %.3f"%(np.log10(jvalue),jsigma)

    rjvalues = {}
    for ft, rtgts in zip(fntypes,rtargets):
        rjvalues[ft] = rtgts[opts.target].get_safe_jvalue(profile)
    
    jfile = join(workdir,'jvalues.yaml')
    out = open(jfile,'w')
    out.write(yaml.dump(dict(jvalue=jvalue,jsigma=jsigma,jnominal=jnominal,rjvalues=rjvalues)))
    out.close()


    spparam = None
    
    # Setup the spectrum
    if opts.spectrum == 'pwl':
        if decay:
            print 'WARNING Decay not defined for powerlaw. exiting.'
            sys.exit()
        indices = config['spectra'][opts.spectrum]
        spectra = [ PowerLaw(index=idx) for idx in indices]
        spparam = indices
        outfiles = [ "%.1f/%s_%s.dat"%(idx,target,suffix) for idx in indices ]
    elif opts.spectrum.startswith("IC"):
        masses = config['spectra'][opts.spectrum]#filesIC.keys()
        spectra = [ FileSpectrum(m=...) for m in masses]
        outfiles = [ "%.1f/%s_logLike2.dat"%(mass,target) mass in masses ]
        spparam = masses
    else:
        specdir = join(basedir,"spectra/%05d"%opts.param)
        infiles = sorted(glob.glob(join(specdir,'out-yield-ex*.dat')))
        spectra = [ FileSpectrum(f,jvalue,sleep=0.5) for f in infiles] 
        outfiles = [ "%.1f/"%opts.param + splitext(basename(f))[0]+"_logLike.dat" for f in infiles ]

    profile_data = {}
    spparam = ['%.1f'%t for t in spparam]
    if opts.spectrum in DMFitFunction.channels() or opts.spectrum.startswith("IC"):
        print 'Output of spectra we want to create'
        for i,m in enumerate(masses):
            print m,str(spectra[i])
    
    print 'Loading ', sedfile
    sed = SED.read_yaml(sedfile)
    # Parse the spectral file
    for i, (spec, outfile) in enumerate(zip(spectra, outfiles)):
        # Deal with files containing NaN
        if spec is None: continue
        pd = []        
        for ft, m, rtgts in zip(fntypes,methods,rtargets):

            target = rtgts[opts.target]
            
            if isinstance(spec,DMFitFunction) or isinstance(spec,PowerLaw):
                if not decay:
                    spec.setp('norm',target.get_safe_jvalue(profile))
            pd.append(sed.get_global_logLike(spec, jsigma, m, ft))
            
#            norm,flux,lnl,plnl,flnl,plnl2 = \
 #               sed.get_global_logLike(spec, jvalue, jsigma, fn)

        (norm, flux, lnl, p1lnl, flnl) = pd[0]
        p2lnl = LnLFn(pd[1][0],pd[1][3])(norm)
        m1lnl = np.zeros(np.shape(p2lnl))
        ### FIXME! TRY WITH THIS ONE FIRST! --  still broken!
        #m1lnl = LnLFn(pd[2][0],pd[2][3])(norm)

        profile_data[spparam[i]] = dict(norm=norm,flux=flux,lnl=lnl,p1lnl=p1lnl,
                                        flnl=flnl,p2lnl=p2lnl,m1lnl=m1lnl)
        
        # Write logLike to text file
        if not opts.text_output: continue
        
        #outfile=join(workdir,splitext(basename(infile))[0]+"_logLike.dat")
        outfile=join(workdir,outfile)
        mkdir(dirname(outfile))
        print "Writing %s"%outfile
        data = np.vstack( (norm, flux, lnl, p1lnl, flnl, p2lnl, m1lnl) ).T
        out = open(outfile,'w')
        out.write('%-15s %-15s %-15s %-15s %-15s %-15s %-15s\n'%('# normPar','flux[ph/cm2/s]','logLike',
                                                                 'p1logLike','flogLike','p2logLike','m1logLike'))
        np.savetxt(out,data,fmt='%-15.6g %-15.6g %-15.6g %-15.6g %-15.6g %-15.6g %-15.6g')
        out.close()

    lnlfile = join(workdir,'%s_%s.yaml'%(target,suffix))
    print "Writing %s"%lnlfile
    tools.yaml_dump(dict(lnldata = profile_data,param=spparam,
                         jvalue=jvalue,jsigma=jsigma,
                         jnominal=jnominal,rjvalues=rjvalues),lnlfile)
    
    """
    if opts.plot:
        import pylab as plt

        sed_emin, sed_emax = sed.get_emins(),sed.get_emaxs()
        sed_energy = np.sqrt(sed_emin*sed_emax)
        
        sed_eflux = sed.get_efluxes()
        sed_lnl=sed.get_logLikes()
        sed_ul = sed.get_ulimits(alpha=0.05)
        
        x = np.logspace(-3,5,100)
        xmin = x.min(); xmax = x.max()
        
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        color='k'; index=2.0
        kwargs = dict(linestyle='none', color=color)
        ax1.errorbar(sed_energy, sed_ul, 
        yerr=[0.5*sed_ul, plt.zeros(len(sed_ul))],
        lolims=True, lw=1, 
        **kwargs)
        ax1.errorbar(sed_energy, sed_ul, 
        xerr=[sed_energy-sed_emin, sed_emax-sed_energy],
        capsize=0, lw=1.25, 
        **kwargs)
        ax1.plot(plt.nan,plt.nan,lw=2,color=color,
        label=r'$\mathrm{dN/dE \propto E^{-%g}}$'%index)
        
        ecut = energy>min(sed_emin)
        ax1.plot(energy[ecut],energy[ecut]**2*dnde[ecut],'-')
        
        ax1.set_xlabel('Energy (MeV)')
        ax1.set_ylabel(r'$\mathrm{\int\ E\,\frac{dN}{dE} dE\ (MeV\,cm^{-2}\,s^{-1})}$')
        plt.legend(loc='upper left')
    """
