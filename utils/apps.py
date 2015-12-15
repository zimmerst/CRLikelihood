#!/usr/bin/env python
import subprocess as sub
import os
import shutil
from os import getenv
from os.path import join, isfile
import numpy as np
import pyfits
from GtApp import GtApp
from tempfile import NamedTemporaryFile, mkdtemp

from dsphs.utils.lsf import bsub
from dsphs.utils.tools import mkdir,pwd

class XGtApp(GtApp):
    def command(self, do_timing=True):
        use_xrootd = bool(int(os.getenv("USE_XROOTD","0")))
        if use_xrootd:
            xrdPath = os.getenv("XRDPRE","None")
            if not isfile(xrdPath):
                raise RuntimeWarning("cannot find xrootd preload path, fall back to posix")
                use_xrootd = False
        chatter = 2
        try:
            chatter = self.pars['chatter']
        except KeyError:
            pass
        if use_xrootd:
            self.app = "%s -g %s"%(xrdPath,self.app)
        if do_timing and os.name == 'posix' and chatter != 0:
            return 'time -p ' + self.app + self.pars()
        else:
            return self.app + self.pars()

def gttsmap(**kwargs):
    app=GtApp("gttsmap")
    return app.run(**kwargs)

def gtexpmap(**kwargs):
    app=GtApp("gtexpmap")
    return app.run(**kwargs)

def gtselect(**kwargs):
    app=XGtApp('gtselect','dataSubselector')
    return app.run(**kwargs)

def gtmktime(**kwargs):
    app = GtApp('gtmktime')
    return app.run(**kwargs)

def gtltcube(**kwargs):
    app = GtApp('gtltcube','Likelihood')
    return app.run(**kwargs)

def gtltsum(**kwargs):
    app=GtApp('gtltsum')
    return app.run(**kwargs)

def gtdiffrsp(**kwargs):
    app=GtApp('gtdiffrsp')
    return app.run(**kwargs)

def gtbin(**kwargs):
    app=GtApp('gtbin')
    return app.run(**kwargs)

def gtmodel(**kwargs):
    app=GtApp('gtmodel')
    return app.run(**kwargs)

def gtincept(**kwargs):
    evfile  = kwargs['evfile']
    tmpfile = NamedTemporaryFile(dir='/scratch/').name
    gtifile = kwargs['gtifile']
    outfile = kwargs['outfile']
    scfile  = kwargs['scfile']
    chatter = kwargs['chatter'] if 'chatter' in kwargs else 2

    # Note, add on gtis to 'evfile'. This is a bit distructive,
    # but should cause no real harm.
    if chatter > 2: print "Copying %s to %s"%(evfile,tmpfile)
    shutil.copy(evfile,tmpfile)
    e = pyfits.open(tmpfile, mode='update')

    # Temporarily address issue https://jira.slac.stanford.edu/browse/OBS-20
    if len(e)==4 and e[2].name=='GTI' and e[3].name=='GTI' and e[3].data==None:
        del(e[3])
        
    if chatter > 2: print "Updating GTIs..."
    g = pyfits.open(gtifile)
    e['GTI'] = g['GTI']

    # And re-write header info (inaccuracy due to pil conversion float to str)
    for ehdu, ghdu in zip(e,g):
        ehdu.header['TSTART'] = ghdu.header['TSTART']
        ehdu.header['TSTOP']  = ghdu.header['TSTOP']

    e.flush()

    gtmktime(evfile=tmpfile,
             outfile=outfile,
             scfile=scfile,
             # Here, apply no additional gtmktime cuts
             filter='T', 
             roicut='no',
             # This will cut on previously added GTI cuts
             apply_filter='yes',
             chatter=chatter
             )
    
    os.remove(tmpfile)

def gtltcuber(queue='long',**kwargs):
    """ Submit a set of ltcubes to the batch and collect them at the end 
    evfile  : List of LAT ft1 fits files
    scfile  : LAT spacecraft file
    outfile : Output summed ltcube
    tmin    : Minimum time
    tmax    : Maximum time"""

    workdir = mkdtemp(dir='./')
    tmp = mkdir('/tmp/%s/'%os.environ['USER'])
    os.environ['PFILES']=tmp+';'+os.environ['PFILES'].split(';')[-1]

    tmin = kwargs.get('tmin',0)
    tmax = kwargs.get('tmax',0)
    zmax = kwargs.get('zmax',0)
    dcostheta = kwargs.get('dcostheta',180)
    binsz = kwargs.get('binsz',1)
    chatter = kwargs.get('chatter',2)
    sleep = kwargs.get('sleep','5m')

    if kwargs['evfile'].startswith('@'):
        evfiles = np.loadtxt(kwargs['evfile'].strip('@'),dtype='str')
    else:
        evfiles = [ kwargs['evfile'] ]

    #### Calculated the ltcubes ####
    jobname = 'gtltcube.%s'%os.path.basename(workdir)
    lst, cmnds,logs = [],[],[]
    i = 0
    for evfile in evfiles:
        #### gtltcube doesn't treat tmin/tmax sensibly ####
        header = pyfits.open(evfile)[0].header
        tstart = int(float(header['TSTART']))
        tstop  = int(float(header['TSTOP']))
        if tstart >= tmax or tstop <= tmin:
            print "Skipping %s"%evfile
            continue
        tbegin = tmin if tmin > tstart else tstart
        tend   = tmax if tmax < tstop else tstop

        i += 1
        ltfile = join(workdir,'ltcube_%03i.fits'%i)
        select = join(workdir,'select_%03i.fits'%i)
        params = dict( evfile=evfile,
                       select=select,
                       ltfile=ltfile,
                       scfile=kwargs['scfile'],
                       tmin=tbegin, tmax=tend,
                       zmax=zmax
                       )
        print params
        select = """gtselect \
 infile=%(evfile)s \
 outfile=%(select)s \
 tmin=%(tmin)s tmax=%(tmax)s zmax=%(zmax)s \
 ra=0 dec=0 rad=180 emin=1 emax=1e6 \
 chatter=4 clobber=yes;
"""%params

        #### tmin and tmax ignored by ltcube ####
        ltcube = """gtltcube\
 evfile="%(select)s" \
 outfile="%(ltfile)s" \
 scfile="%(scfile)s" \
 dcostheta=0.025 binsz=1 \
 zmax=%(zmax)d tmin=%(tmin)d tmax=%(tmax)d \
 chatter=4 clobber=yes;
"""%params
        cmnd = select + '\n\n' + ltcube
        cmnds.append(cmnd)
        lst.append(ltfile)
    bsub(jobname,cmnds,logfile=None,q=queue,sleep=sleep)
    ltcube_lst = join(workdir,'ltcubes.lst')
    np.savetxt(ltcube_lst,lst,fmt="%s")

    #### Combine ltcubes ####
    depend = "done(%s)"%jobname
    jobname='merge.%s'%os.path.basename(workdir)
    script = join(workdir,'merge.py')
    params = dict( script = script,
                   workdir=workdir)
    #cmnd ="""python %(script)s; rm -rf %(workdir)s"""%params
    cmnd ="""python %(script)s"""%params

    open(script,'w').write("""\
from dsphs.utils.apps import gtltsumr;
import os
gtltsumr("@%s","%s")
""" % (ltcube_lst,kwargs['outfile']) )
    bsub(jobname,[cmnd],q=queue,w=depend)

def gtltsumr(infiles, outfile):
    """ Run gtltsum recursively to sum all livetime cubes. """
    if isinstance(infiles,str) and infiles.startswith('@'):
        infiles = np.loadtxt(infiles.strip('@'),dtype='str').tolist()
        # For single entry in lst
        if isinstance(infiles,str): infiles = [infiles]

    if len(infiles) == 1:
        shutil.copy(infiles[0], outfile)
        return

    infiles = list(infiles) # in case generator

    if len(infiles) < 2: raise Exception('Must sum >= 2 livetime cubes')

    temp = NamedTemporaryFile(suffix='.fits',delete=False)
    tempfile = temp.name

    gtltsumr.i=1

    def sum_ltcube(infile1,infile2):
        print 'Merging file %s and %s (%s/%s)' % (infile1, infile2,gtltsumr.i,len(infiles))
        gtltsumr.i+=1
        if infile1 == outfile:
            gtltsum(infile1=infile1, infile2=infile2, outfile=tempfile)
            return tempfile
        else:
            gtltsum(infile1=infile1, infile2=infile2, outfile=outfile)
            return outfile
    try:
        accum=reduce(sum_ltcube, infiles)
    except RuntimeError, msg:
        if os.path.exists(outfile): os.remove(outfile)
        raise RuntimeError(msg)

    if accum == tempfile:
        shutil.move(tempfile, outfile)
    else:
        os.remove(tempfile)

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()
