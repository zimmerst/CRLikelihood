#!/usr/bin/env python
import subprocess as sub
import numpy as np
from numpy import cos,sin,arccos,pi
import os, re, time, copy
import numbers
import shutil, shlex
import pyfits
import yaml
try:
    from skymaps import SkyDir
except ImportError:
    print 'Could not find skymaps, trying FSSC init'
    from pyLikelihood import SkyDir

def stversion():
    if 'VERSION' in os.environ:
        return os.environ['VERSION'].split('-')
    else:
        raise Exception('VERSION environment variable is not defined.')

def relocalizeSpatialMap(iname,oname,center,proj='AIT',galactic=True):
    safe_copy(iname, oname, sleep=10) # make a safe copy first, rest is pyfits 'magic'
    pf = pyfits.open(oname,mode='update')
    # clean up header of existing stuff
    h = pf[0].header
    ctypes = ["GLON","GLAT"] if galactic else ["RA","DEC"]
    ctypes_new = []
    for ct in ctypes:
        while (len(ct)+len(proj))<=8: ct+="-" # need to fill w/ -
        ct+=proj.upper()
        ctypes_new.append(ct)
    crvals = [center.l(),center.b()] if galactic else [center.ra(),center.dec()]
    my_map = {"CTYPE":ctypes_new,"CRVAL":crvals}
    for key in ['CTYPE','CRVAL']:
        for i in [0,1]:
            hkey = "%s%i"%(key,i+1)
            h.update(hkey,my_map[key][i])
    pf.flush() # write changes to file.
    
def makeSafeName(srcname):
    rep = {".":"d","+":"p","-":"n"}
    for key in rep:
        srcname = srcname.replace(key,rep[key])
    return srcname

def pwd():
    # Careful, won't work after a call to os.chdir...
    return os.environ['PWD']

def mkdir(dir):
    if not os.path.exists(dir):  os.makedirs(dir)
    return dir

def config():
    npconfig()
    mplconfig()

def mkscratch():
    if os.path.exists('/scratch/'):    
        return(mkdir('/scratch/%s/'%os.environ['USER']))
    elif os.path.exists('/tmp/'):
        return(mkdir('/tmp/%s/'%os.environ['USER']))
    else:
        raise Exception('...')

def mplconfig():
    #try:             mplconfigdir = os.environ['MPLCONFIGDIR']
    #except KeyError: 
    #    if os.path.exists('/scratch'):
    #        mplconfigdir = '/scratch/%s/'%os.environ['USER']
    #    else:
    #        mplconfigdir = '/tmp/%s/'%os.environ['USER']
    mplconfigdir = '/tmp/%s/'%os.environ['USER']
    os.environ['MPLCONFIGDIR'] = mplconfigdir
    os.environ['MATPLOTLIBRC'] = '%s/.matplotlib/'%os.environ['HOME']
    return mkdir(mplconfigdir)

def npconfig():
    np.seterr(all='ignore')

def round_figures(x, n=5):
    """Returns x rounded to n significant figures."""
    if np.asarray(x).ndim == 0: 
        return round(x, int(n - np.ceil(np.nan_to_num(np.log10(abs(x))))))
    else:
        return np.array( [round(_x, int(n - np.ceil(np.nan_to_num(np.log10(abs(_x)))))) for _x in x] )

def bootstrap(infile,outfile,seed=0,scale=1.):
    if not os.path.splitext(infile)[-1] == '.fits':
        raise Exception("Cannot bootstrap non-fits file.")
    
    # If seed set to None, just copy the file
    if seed is None:
        print "Copying: %s"%infile
        shutil.copyfile(infile,outfile)
        return

    print "Bootstrap: %s"%infile
    hdu = pyfits.open(infile)

    if scale is None: scale = 1.
    ninput = len(hdu[1].data)
    noutput = int(ninput/scale)

    print "\t Original: %i"%ninput

    # Do the bootstrapping
    np.random.seed(int(seed))
    index = np.random.randint(0,ninput,noutput)
    hdu[1].data = hdu[1].data[index]

    print "\t Bootstrap: %i"%len(hdu[1].data)
    hdu.writeto(outfile)
    return

class astroserver(object):
    """ Wrapper around the glast astroserver.
    Pass in command-line args as kwargs changing
    '_' to '-'. Checks kwargs f"""

    def __init__(self):
        self.exe = "/u/gl/glast/astroserver/prod/astro"
        p = sub.Popen((self.exe + " --help").split(), stdout=sub.PIPE,stderr=sub.PIPE)
        stdout,stderr=p.communicate()
        self.opts = stdout

    def __call__(self,arg,**kwargs):
        command = self.exe
        for key,val in kwargs.items():
            kwarg = key.replace('_','-')
            if kwarg not in self.opts: raise Exception("%s\n%s"%(self.exe,self.opts))
            command += " -%s %s"%(kwarg,val) 
        command += " %s"%arg
        return command

def isRoman(string):
    pattern = '^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$'
    return not (re.search(pattern, string) is None)

def isArabic(string):
    pattern = '^[0-9]'
    return not (re.search(pattern, string) is None)

def isNumeral(string):
    return isRoman(string) or isArabic(string)

isnum = lambda x: isinstance(x, numbers.Real)

def nancut(a, axis=None):
    if len(a.shape) > 1:
        cut = ~np.isnan(np.sum(a,axis=axis))
    else:
        cut = ~np.isnan(a)
    return a[cut]

def isDMFit(string):
    from uw.darkmatter.spectral import DMFitFunction
    l = DMFitFunction.channel_mapping.values()
    return string in [item for sublist in l for item in sublist]

def event_mask(all_events,allowed_events):
    """ Takes unique event identifier
    and compares it to a set of allowed"""
    return np.array( [ x in allowed_events for x in all_events ])

def hms2decimal(hms):
    DEGREE = 360.
    HOUR = 24.
    MINUTE = 60.
    SECOND = 3600.

    if not isinstance(hms,str):
        raise Exception("...")
    
    if 'h' in hms:
        hour,minute,second = np.array(re.split('[hms]',hms))[:3].astype(float)
        #sign = -1 if hour < 0 else 1
        decimal = (hour + minute * 1./MINUTE + second * 1./SECOND)*(DEGREE/HOUR)
    else:
        sign = -1 if hms.startswith('-') else 1
        degree,minute,second = np.array(re.split('[dms]',hms))[:3].astype(float)
        decimal = abs(degree) + minute * 1./MINUTE + second * 1./SECOND
        decimal *= sign

    return decimal

#adapted from iLat by J. Chiang 
def Angdist(x):
    """ Angular distance in radians corresponding to a cosinus  """ 
    x = np.asarray(x)
    if x.ndim == 0: x = x.reshape((1))
    angdist=np.zeros(len(x))
    #try: 
    #    angdist=np.zeros(len(x))
    #except TypeError:
    #    angdist = 0
    if any(x > 1.00001):
        print "Error: x must be smaller than 1"
    else:
        angdist += arccos(x)*(x<1)
        angdist += np.pi/2 * (1-x.astype(int))*(x>=1)
	return angdist 

#adapted from iLat by J. Chiang 
def angdist( ra1, dec1, ra2, dec2 ):
	"""Angular separation in degrees between sky coordinates"""
	ra1 = ra1*np.pi/180.;	dec1 = dec1*np.pi/180.
	ra2 = ra2*np.pi/180.;	dec2 = dec2*np.pi/180.
	mu = (cos(dec1)*cos(ra1)*cos(dec2)*cos(ra2)
		  + cos(dec1)*sin(ra1)*cos(dec2)*sin(ra2) + sin(dec1)*sin(dec2))
	return Angdist(mu)*180./np.pi

def random_glon(catalog, glats, dist=1., roisize=1., seed=0, samplesize=None):
    """ Choose random glons with fixed glats
    """
    print "Selecting random GLONs"
    from uw.like.roi_catalogs import BaseCatalog2FGL
    from uw.like.roi_extended import ExtendedSource
    from skymaps import SkyDir

    np.random.seed(seed)
    if not size is None:
        # pick exactly n from glats
        if size>glats:
            print "WARNING: more fields requested than the positions that are available"
        glats = np.random.choice(glats,size)
    chunk = 10*len(glats)

    if isinstance(catalog, BaseCatalog2FGL):
        print "Reading pointlike catalog..."
        skydir = SkyDir(0,0,SkyDir.GALACTIC)
        GLON, GLAT = [],[]
        for source in catalog.get_sources(skydir,180.):
            GLON.append(source.skydir.l())
            GLAT.append(source.skydir.b())
        GLON,GLAT = np.asarray((GLON,GLAT))
    else: 
        print "Reading FITS catalog..."
        hdulist = pyfits.open(catalog)
        # GLON [0,360]; GLAT [-90,90]; GC [0,0]
        GLON = hdulist[1].data.field('GLON')
        GLAT = hdulist[1].data.field('GLAT')

    points = []
    for i,b in enumerate(glats):
        v = np.random.uniform(size=chunk)    
        glon=np.degrees(2*pi*v)

        # Check points don't overlap with catalog or themselves
        for l in glon:
            if len(points) > 0:
                exist = np.array(points)
                if sum(angdist(l,b,exist[:,0],exist[:,1]) < 2*roisize ) != 0:
                    continue
            if sum(angdist(l,b,GLON,GLAT) < dist) == 0:
                points.append([l,b])
                break
    points = np.array(points)
    if len(points) != len(glats):
        raise Exception("Length Mismatch")
    print 'Random GLONS : ', len(points)
    return points

def random_skydir(catalog, size=100, glatcut=None, dist=1., roisize=1., seed=0, samplesize=None):
    """
    catalog - can be either point to a fits file or a pointlike catalog object
    """
    print "Selecting random points on the sky"
    from uw.like.roi_catalogs import BaseCatalog2FGL
    from uw.like.roi_extended import ExtendedSource
    from skymaps import SkyDir
    if not samplesize is None:
        size = int(samplesize)
        print '*custom hack: required a set size of %i*'%samplesize
    np.random.seed(seed)
    chunk = 10*size

    if isinstance(catalog, BaseCatalog2FGL):
        print "Reading pointlike catalog..."
        skydir = SkyDir(0,0,SkyDir.GALACTIC)
        GLON, GLAT, RADIUS = [],[],[]
        for source in catalog.get_sources(skydir,180.):
            GLON.append(source.skydir.l())
            GLAT.append(source.skydir.b())
            # Extra padding around extended sources
            if isinstance(source, ExtendedSource): RADIUS.append(4.0)
            else:                                  RADIUS.append(0.0)
        GLON,GLAT,RADIUS = np.asarray((GLON,GLAT,RADIUS))
    else: 
        print "Reading FITS catalog..."
        hdulist = pyfits.open(catalog)
        # GLON [0,360]; GLAT [-90,90]; GC [0,0]
        GLON = hdulist[1].data.field('GLON')
        GLAT = hdulist[1].data.field('GLAT')
        RADIUS = np.zeros(len(GLON))

    print len(GLON),len(GLAT)
     
    points = []
    i = 0
    while len(points) < size:
        # Create random numbers
        u = np.random.uniform(size=chunk); v = np.random.uniform(size=chunk)
        # In radians
        #glon=2*pi*v - pi
        #glat=arccos(2*u - 1) - pi/2
         
        #  glon [0,360];  glat [-90,90]
        glon=np.degrees(2*pi*v)
        glat=np.degrees(np.arccos(2*u - 1)) - 90
     
        # Cut on glat
        if glatcut:
            cut = (glat>=glatcut)|(glat<=-glatcut)
            glat = glat[cut]
            glon = glon[cut]
     
        # Check points don't overlap with catalog or themselves
        for l,b in zip(glon,glat):
            if len(points) > 0:
                exist = np.array(points)
                if sum(angdist(l,b,exist[:,0],exist[:,1]) < 2*roisize ) == 0:
                    if sum( angdist(l,b,GLON,GLAT) < (dist + RADIUS) ) == 0:
                        points.append([l,b])
            else:
                #if sum(angdist(l,b,glon,glat) < dist ) == 1:
                #    if sum(angdist(l,b,GLON,GLAT) < dist) == 0:
                #        points.append([l,b])
                if sum(angdist(l,b,GLON,GLAT) < dist) == 0:
                    points.append([l,b])

        print 'Chunk %i : '%i,len(points)
        i += 1
    points = np.array(points)[:size]
    print 'Random Points : ', len(points)
    return points

def coverage_range(spectrum,param):
    from uw.darkmatter.spectral import DMFitFunction
    if spectrum in DMFitFunction.channels():
        logpar = np.log10(param)
        # mu+mu-
        if DMFitFunction.channel2int(spectrum) == 2:
            slope = 1.25
            ymin,ymax = -27.3,-25.3
        # tau+tau-
        elif DMFitFunction.channel2int(spectrum) == 3:
            slope = 1.25
            ymin,ymax = -28.,-26.
        # b-bbar
        elif DMFitFunction.channel2int(spectrum) == 4:
            slope = 1.
            ymin,ymax = -28.,-26.
        # W+W-
        elif DMFitFunction.channel2int(spectrum) == 7:
            slope = 1.1
            ymin,ymax = -28.,-26.
        return 10**(slope*logpar+ymin),10**(slope*logpar+ymax)
    else: 
        return np.nan,np.nan

def random_sigmav(spectrum,param,sigmav=None):
    if sigmav is None:
        low,high = np.log10(coverage_range(spectrum,param))
    else:
        low,high = np.log10(sigmav[0]), np.log10(sigmav[1])
    sigmav = 10**( np.random.uniform(low,high) )
    return sigmav


def parse_sleep(sleep):
    MINUTE=60
    HOUR=60*MINUTE
    DAY=24*HOUR
    WEEK=7*DAY
    if isinstance(sleep,float) or isinstance(sleep,int):
        return sleep
    elif isinstance(sleep,str):
        try: return float(sleep)
        except ValueError: pass
        
        if sleep.endswith('s'): return float(sleep.strip('s'))
        elif sleep.endswith('m'): return float(sleep.strip('m'))*MINUTE
        elif sleep.endswith('h'): return float(sleep.strip('h'))*HOUR
        elif sleep.endswith('d'): return float(sleep.strip('d'))*DAY
        elif sleep.endswith('w'): return float(sleep.strip('w'))*WEEK
        else: raise ValueError
    else:
        raise ValueError

def sleep(sleep):
    return time.sleep(parse_sleep(sleep))

def get_resources():
    import resource
    usage=resource.getrusage(resource.RUSAGE_SELF)
    return '''usertime=%s systime=%s mem=%s mb
           '''%(usage[0],usage[1],
                (usage[2]*resource.getpagesize())/1000000.0 )

def safe_copy(infile, outfile, sleep=10, attempts=10):
    infile = infile.replace("@","") if infile.startswith("@") else infile
    # Try not to step on any toes....
    sleep = parse_sleep(sleep)
    if infile.startswith("root:"):
        print 'file is on xrootd - switching to XRD library'
        cmnd = "xrdcp %s %s"%(infile,outfile)
    else:
        cmnd = "cp %s %s"%(infile,outfile)
    i = 0
    print "Attempting to copy file..."
    while i < attempts: 
        status = sub.call(shlex.split(cmnd))
        if status == 0: 
            return status
        else:
            print "%i - Copy failed; sleep %ss"%(i,sleep)
            time.sleep(sleep)
        i += 1
    raise Exception("File *really* doesn't exist: %s"%infile)

def yaml_load(filename):
    """ Load a yaml file (use libyaml when available) """
    if not os.path.exists(filename):
        raise Exception('File does not exist: %s'%(filename))
    if hasattr(yaml,'CLoader'):
        Loader = yaml.CLoader
    else: 
        Loader = yaml.Loader

    try:
        ret = yaml.load(open(filename),Loader=Loader)
    except:
        ret = yaml.load(open(filename),Loader=yaml.Loader)
    return ret
    
def yaml_dump(x, filename):    
    """ Dump object to a yaml file (use libyaml when available) 
        x        : output to dump to the file
        filename : output file (can be file-type or path string)
    """
    
    if hasattr(yaml,'CDumper'):
        Dumper = yaml.CDumper
    else: 
        Dumper = yaml.Dumper

    if isinstance(filename, basestring):
        out = open(filename,'w')
    elif isinstance(filename, file):
        out = filename
    else:
        raise Exception("Unrecognized file: ",filename)
    out.write( yaml.dump(x) )
    out.close()

def expand_list_arg(args,delimiter=',',argtype=str):

    if args is None: return args    
    elif not isinstance(args,list):
        o = args.split(delimiter)
    else:
        o = []
        for arg in args: o += arg.split(delimiter)

    o = [argtype(t) for t in o]
        
    return o
    
def update_config(infile, outfile, **kwargs):
    config = yaml_load(infile)
    update_dict(config,kwargs,add_new_keys=True)
    yaml_dump(config, outfile)

    #import yaml
    #config = yaml.load(open(infile))
    #config.update(kwargs)
    #out = open(outfile,'w')
    #out.write(yaml.dump(config))
    #out.close()

def update_dict(d0,d1,add_new_keys=False,append=False):
    """Recursively update the contents of python dictionary d0 with
    the contents of python dictionary d1."""

    if d0 is None or d1 is None: return
    
    for k, v in d0.iteritems():

        if not k in d1: continue

        if isinstance(v,dict) and isinstance(d1[k],dict):
            update_dict(d0[k],d1[k],add_new_keys,append)
        elif isinstance(v,list) and isinstance(d1[k],str):
            d0[k] = d1[k].split(',')            
        elif isinstance(v,np.ndarray) and append:
            d0[k] = np.concatenate((v,d1[k]))
        else: d0[k] = d1[k]

    if add_new_keys:
        for k, v in d1.iteritems(): 
            if not k in d0: d0[k] = d1[k]
    
def load_config(defaults,config=None,**kwargs):
    """Create a configuration dictionary.  The defaults tuple is used
    to create a dictionary in which valid key values are predefined.
    The config dictionary and kwargs are then used to update the
    values in the default configuration dictionary."""

    o = {}
    for item in defaults:
        
        item_list = [None,None,'',None,str]
        item_list[:len(item)] = item        
        key, value, comment, groupname, item_type = item_list

        if len(item) == 1:
            raise Exception('Option tuple must have at least one element.')
                    
        if value is None and (item_type == list or item_type == dict):
            value = item_type()
            
        keypath = key.split('.')

        if len(keypath) > 1:
            groupname = keypath[0]
            key = keypath[1]
                    
        if groupname:
            group = o.setdefault(groupname,{})
            group[key] = value
        else:
            o[key] = value
            
    update_dict(o,config)
    update_dict(o,kwargs)

    return o

def addTargetSetToCatalog(catalog,config):
    # returns a catalog which includes the target set
    from uw.like.roi_catalogs import SourceCatalog
    from skymaps import SkyDir
    from dsphs.base.target import load_targets_from_config

    if not isinstance(catalog,SourceCatalog):
        raise Exception("before merging, catalog must be source catalog")
    if config['randomize']['profile'] is None: 
        config['randomize']['profile']='point'
        print RuntimeWarning("no randomize.profile chosen, falling back to point")
    if config['randomize']['vetoyaml'] is None:
        # do nothing. 
        return catalog
    cfg = copy.deepcopy(config)
    del config
    # work entirely on 2nd config.
    cfg['targetfile'] = cfg['randomize']['vetoyaml'] 
    target_set = load_targets_from_config(cfg) # this one load what's in vetoyaml
    target_set.set_profile(cfg['randomize']['profile'])
    print 'loaded %i sources with profile %s'%(len(target_set),cfg['randomize']['profile'])
    tlist = [t.get_source_model(profile=cfg['randomize']['profile']) for t in target_set]
    ps, ds = [], []
    if cfg['randomize']['profile']=="point": ps = tlist
    else: ds = tlist
    #print ps, ds
    print 'merging %i sources with catalog containing %i sources'%(len(target_set),
                                                                 len(catalog.get_sources(SkyDir(0,0,SkyDir.GALACTIC),radius=180.)))
    cnew = catalog.merge_lists(SkyDir(0,0,SkyDir.GALACTIC),radius=180.,user_point_list=ps,user_diffuse_list=ds)
    catalog.sources = cnew[0]+cnew[1] # joint catalogs
    del cfg # cleanup
    print 'size after merger %i'%len(catalog.sources)
    return catalog

def tolist(x):
    """ convenience function that takes in a 
        nested structure of lists and dictionaries
        and converts everything to its base objects.
        This is useful for dupming a file to yaml.

        (a) numpy arrays into python lists

            >>> type(tolist(np.asarray(123))) == int
            True
            >>> tolist(np.asarray([1,2,3])) == [1,2,3]
            True

        (b) numpy strings into python strings.

            >>> tolist([np.asarray('cat')])==['cat']
            True

        (c) an ordered dict to a dict

            >>> ordered=OrderedDict(a=1, b=2)
            >>> type(tolist(ordered)) == dict
            True

        (d) converts unicode to regular strings

            >>> type(u'a') == str
            False
            >>> type(tolist(u'a')) == str
            True

        (e) converts numbers & bools in strings to real represntation,
            (i.e. '123' -> 123)

            >>> type(tolist(np.asarray('123'))) == int
            True
            >>> type(tolist('123')) == int
            True
            >>> tolist('False') == False
            True
    """
    if isinstance(x,list):
        return map(tolist,x)
    elif isinstance(x,dict):
        return dict((tolist(k),tolist(v)) for k,v in x.items())
    elif isinstance(x,np.ndarray) or \
            isinstance(x,np.number):
        # note, call tolist again to convert strings of numbers to numbers
        return tolist(x.tolist())
#    elif isinstance(x,PhaseRange):
#        return x.tolist(dense=True)
#    elif isinstance(x,OrderedDict):
#        return dict(x)
    elif isinstance(x,basestring) or isinstance(x,np.str):
        x=str(x) # convert unicode & numpy strings 
        try:
            return int(x)
        except:
            try:
                return float(x)
            except:
                if x == 'True': return True
                elif x == 'False': return False
                else: return x
    else:
        return x

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()


