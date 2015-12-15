#!/usr/bin/env python
import re, os
from os.path import join
import collections
import yaml, copy
import numpy as np

import xml.etree.cElementTree as ElementTree
try:
    from skymaps import SkyDir
except ImportError:
    print 'Could not find skymaps, trying FSSC init'
    from pyLikelihood import SkyDir

from dsphs.base.cluster_profiles import CoredProfile

from dsphs.pointlike.Sources import PointSource, ExtendedSource # good
from dsphs.pointlike.SpatialModels import SpatialModel, Gaussian, Disk, SpatialMap
from dsphs.pointlike.Models import PowerLaw, BrokenPowerLaw, SmoothBrokenPowerLaw

from dsphs.utils.set import *
from dsphs.utils.tools import isRoman, isArabic, isNumeral, isnum, hms2decimal
from dsphs.utils.tools import random_skydir, random_glon, yaml_dump
from dsphs.utils.dmfit import channel2int

from dsphs.pointlike.xml_parsers import write_sources


class Target(SetItem):
    def __init__(self, name, **kwargs):
        # Setup all the naming
        super(Target,self).__init__(name)
        self.force_extvalue=None # force the extension!
        self.__dict__.update(kwargs)

        for x in ['ra','dec','l','b']:
            try: self.__dict__.pop(x)
            except KeyError: pass

        if 'codename' in kwargs: self.codename = kwargs['codename']
        else: self._create_codename()

        if 'abbreviation' in kwargs: self.abbreviation = kwargs['abbreviation']
        else: self._create_abbreviation()

        if 'nickname' in kwargs: self.nickname = kwargs['nickname']
        else: self._create_nickname()

        if 'color' in kwargs: self.color = kwargs['color']
        else: self.color='black'

        if 'skydir' in kwargs:
            if not isinstance(kwargs['skydir'],SkyDir):
                print "skydir must be of type skymaps.SkyDir"
            self.skydir = kwargs['skydir']
        elif 'ra' in kwargs and 'dec' in kwargs:
            ra, dec = kwargs['ra'],kwargs['dec']
            if isinstance(ra,str) and 'h' in ra:
                self.skydir = SkyDir(hms2decimal(ra),hms2decimal(dec),SkyDir.EQUATORIAL)
            else:
                self.skydir = SkyDir(ra,dec,SkyDir.EQUATORIAL)
        elif 'l' in kwargs and 'b' in kwargs:
            self.skydir = SkyDir(kwargs['l'],kwargs['b'],SkyDir.GALACTIC)
        else:
            self.skydir = SkyDir(0 ,0, SkyDir.GALACTIC)
            #raise Exception("Unable to create skydir.")

        if 'distance' in kwargs: self.distance = kwargs['distance']
        else: self.distance = None

        if 'dsigma' in kwargs: self.dsigma = kwargs['dsigma']
        else: self.dsigma = None

        if 'extension' in kwargs: self.extension = kwargs['extension']
        else: self.extension = {'point': None}

        if 'extsigma' in kwargs: self.extension = kwargs['extsigma']
        else: self.extsigma = {'point': None}

        if 'refs' in kwargs: self.refs = kwargs['refs']
        else: self.refs = None

    def __str__(self):
        return self.codename
 
    def set_name(self,name,codename=None):
        self.name = name
        self._create_codename(codename)
        self._create_abbreviation()
        self._create_nickname()
 
    def _create_codename(self, codename = None):
        """ For creating directors and filenames """
        if codename is not None: 
            self.codename = codename
        else:
            split_name = self.name.split()
            # Try (and probably fail) to do something clever ...
            self.codename = ''.join("%s_"%(s.lower() if not isNumeral(s) else s) 
                                    for s in split_name).strip('_')

    def _create_abbreviation(self, abbreviation = None):
        """ Fill the abbreviation """
        if abbreviation is not None:
            # Hopefully this is given
            self.abbreviation = abbreviation
        else:
            # Try (and probably fail) to do something clever ...
            split_name = self.name.split()
            if isArabic(split_name[-1]):
                number = " %s"%split_name.pop(-1)
            else:
                number = ''

            if len(split_name) > 1:
                self.abbreviation = split_name[0][0] + split_name[1][:2]
            else:
                self.abbreviation = split_name[0][:3]
            self.abbreviation += number
        self.abbr = self.abbreviation


    def _create_nickname(self, nickname = None):
        """ Create a shortened name for plot legends etc."""
        
        if nickname is not None:
            # Hopefully this is given
            self.nickname = nickname
        else:
            # Try (and probably fail) to do something clever ...
            split_name = self.name.split()
            self.nickname=self.name.replace(' ','')
            if len(self.nickname) > 9:
                if len(split_name) == 1:
                    #Remove the vowels
                    self.nickname = re.sub('[aeiou]','',self.nickname)
                elif len(split_name) == 2 and not isNumeral(split_name[1]) :
                    #Take the first name
                    self.nickname = split_name[0]
                else:
                    #Push everything together
                    self.nickname = split_name[0][0].join(s for s in split_name[1:])

    def ra(self): return self.skydir.ra()
    def dec(self): return self.skydir.dec()

    def l(self): return self.skydir.l()
    def b(self): return self.skydir.b()

    def glon(self): return self.skydir.l()
    def glat(self): return self.skydir.b()

    def get_extension(self, profile='point'):
        return self.extension[profile]

    def get_safe_extension(self, profile='point'):
        return np.nan_to_num(self.extension[profile])
    def force_extension(self,ext):
        self.force_extvalue=ext

    def toDict(self):

        import inspect
        
        odict = {}
        for k,v in self.__dict__.iteritems():

            if isinstance(v,SkyDir):
                odict[k] = {'ra' : v.ra(), 'dec' : v.dec(), 'l' : v.l(), 'b' : v.b() }
            elif not k.startswith('__') and not inspect.ismethod(v):
                odict[k] = v
        return odict
    
    def plot_kwargs(self):
        kwargs = dict(color = self.color,
                      label = self.name,
                      lw = 1.5,
                      linestyle = '-',
                      )
        return kwargs


class DMTarget(Target):
    """ Subclass of Target for sources with dark matter parameters 
    (e.g., DM spectra, sigmav, jfactors, etc.)
    """

    def __init__(self, name, **kwargs):
        super(DMTarget,self).__init__(name, **kwargs)

        self.default_jvalue = 1e18 
        self.default_jsigma = 0.20
        self.z = 0.
        # Dark matter characteristics
        if 'profile' in kwargs: self.profile = kwargs['profile']
        else: self.profile = 'point'

        self._setup_profiles(kwargs['profiles'] if 'profiles' in kwargs else None)

    def _setup_profiles(self, profiles):
        if profiles is None:
            point = dict( jvalue = self.default_jvalue,
                          jsigma = self.default_jsigma,
                          extension = 0,
                          extsigma  = 0,
                          )
            self.profiles = dict( point = point )
        else:
            self.profiles = profiles


    def get_profile_value(self, key, profile=None):
        if profile is None: profile = self.profile
        if key is None: return self.profiles[profile]
        try:       return self.profiles[profile][key]
        except KeyError: return np.nan

    def set_profile_value(self, key, value, profile=None):
        if profile is None: profile = self.profile
        if key not in self.profiles[profile]: raise KeyError
        self.profiles[profile][key] = value

    def get_profile_ra(self,profile=None):
        return self.get_profile_value(key='ra',profile=profile)

    def get_profile_dec(self,profile=None):
        return self.get_profile_value(key='dec',profile=profile)
    
    def get_profile_glat(self,profile=None):
        return self.get_profile_value(key='b',profile=profile)
    
    def get_profile_glon(self,profile=None):
        return self.get_profile_value(key='l',profile=profile)
    
    def get_jvalue(self, profile=None):
        jvalue = self.get_profile_value(key='jvalue',profile=profile)
        if jvalue < 100: jvalue = 10**jvalue
        return jvalue

    def get_safe_jvalue(self, profile=None):
        jvalue = self.get_jvalue(profile)
        return self.default_jvalue if np.isnan(jvalue) else jvalue 

    def set_jvalue(self, value, profile=None):
        self.set_profile_value(key='jvalue',value=value,profile=profile)

    # Eventually we want to convert over to "jfactor" instead of "jvalue"...
    get_jfactor = get_jvalue
    get_safe_jfactor = get_safe_jvalue
    set_jfactor = set_jvalue

    

    def get_jsigma(self, profile=None):
        return self.get_profile_value('jsigma',profile)

    def get_safe_jsigma(self, profile=None):
        jsigma = self.get_jsigma(profile)
        return self.default_jsigma if np.isnan(jsigma) else jsigma

    def get_extension(self, profile=None):
        extension = self.get_profile_value('extension',profile)
        if self.force_extvalue:
            extension = self.force_extvalue
        return 10**extension

    def get_safe_extension(self, profile='point'):
        spatial_model = self.get_spatial_model(profile=profile)
        # spatial_model is None if sigma is not defined
        if spatial_model is None: return None
        elif isinstance(spatial_model,SkyDir): return 0.0
        else:                     return spatial_model.sigma

    def get_extsigma(self, profile=None):
        return self.get_profile_value('extsigma',profile)

    def get_safe_extsigma(self, profile=None):
        jsigma = self.get_extsigma(profile)
        return 0 if np.isnan(extsigma) else extsigma

    def get_profile_dir(self, profile=None):
        ra_dir = self.get_profile_ra(profile)
        dec_dir = self.get_profile_dec(profile)
        if np.isnan(ra_dir) or np.isnan(dec_dir):
            pass
        else:
            print "Using RA,DEC = %4.2f %4.2f"%(ra_dir,dec_dir)
            return SkyDir(ra_dir,dec_dir)
        glon_dir = self.get_profile_glon(profile)
        glat_dir = self.get_profile_glat(profile)
        if np.isnan(glon_dir) or np.isnan(glat_dir):
            pass
        else:
            print "Using L,B = %4.2f %4.2f"%(glon_dir,glat_dir)
            return SkyDir(glon_dir,glat_dir,SkyDir.GALACTIC)
        print "Using default direction"
        return self.skydir


    def get_spectral_model(self, spectrum='bbsrc', param=100, **kwargs):
        """ Allow the user to be too clever..."""
        if spectrum == 'pwl':
            index = param
            if 'sigmav' in kwargs: print "Popping 'sigmav' from kwargs..."
            if 'sigmav' in kwargs: kwargs.pop('sigmav')
            if 'norm' not in kwargs: kwargs['norm'] = 1e-11
            spectral_model = PowerLaw(index=index, **kwargs)
        elif spectrum == 'BrokenPowerLaw':
            if 'norm' not in kwargs: kwargs['norm'] = 1e-11
            pars = {'Index_1':None,'Index_2':None,'E_break':None}
            if isinstance(param,list):
                pars = dict(zip(["Index_1","Index_2","E_break"]),param)
            elif isinstance(param, dict):
                pars.update(param)
            else:
                raise Exception("could not parse BrokenPowerLaw target model")
            pars.setdefault("norm",kwargs['norm'])
            spectral_model = BrokenPowerLaw(**pars)
            return spectral_model
        elif spectrum == 'SmoothBrokenPowerLaw':
            if 'norm' not in kwargs: kwargs['norm'] = 1e-11
            pars = {'Index_1':None,'Index_2':None,'E_break':None,"beta":None}
            if isinstance(param,list):
                pars = dict(zip(["Index_1","Index_2","E_break","beta"]),param)
            elif isinstance(param, dict):
                pars.update(param)
            else:
                raise Exception("could not parse BrokenPowerLaw target model")
            pars.setdefault("norm",kwargs['norm'])
            pars.setdefault("beta",0.1)
            spectral_model = SmoothBrokenPowerLaw(**pars)
            return spectral_model
        else:            
            raise Exception("Cannot understand parsed Spectrum.")
        #             channel0=DMFitFunction.channel2int(spectrum)
        #             mass = param
        #             if 'norm' in kwargs: 
        #                 print "Moving 'norm' to sigmav..."
        #                 kwargs['sigmav'] = kwargs['norm']
        #             kwargs['norm'] = self.get_safe_jvalue()             
        #             spectral_model = DMFitFunction(channel0=channel0, mass=mass, **kwargs)
        spectral_model.set_default_limits(oomp_limits=True)
        scale = spectral_model.get_scale(0)
        spectral_model.set_limits(0,1e-10*scale,1e10*scale,scale=scale)
        return spectral_model

    # ADW: Temporary chnage
    #def get_spatial_model(self, profile='point', radius=0.5, npix=500):
    def get_spatial_model(self, profile='point', truncate=0.5, npix=500,galactic=True):

        profile_name = profile.lower().split('_')[0]        
        center = self.get_profile_dir(profile)

        if profile_name == 'point':
            return center
        elif profile_name == 'disk':
            spatial_model = Disk(center=center)
            truncate = 2.*truncate
        elif profile_name == 'gauss':
            spatial_model = Gaussian(center=center)
            # don't need this here.
            #         elif profile_name == 'nfw':
            #             spatial_model = NFW(center=center)    
            #         elif profile_name == 'einasto':
            #             spatial_model = Einasto(center=center)
            #         elif profile_name == 'burkert':
            #             spatial_model = Burkert(center=center)
        elif profile_name == 'cored':
            core = self.get_profile_value("core",profile)
            rvir = self.get_profile_value("rvir",profile)
            idx  = self.get_profile_value("index",profile)
            kwargs = {"center":center,"rvir":rvir,
                      "index":idx,"core":core}
            print 'FIXME: cored profile pars: {}'.format(kwargs)
            spatial_model = CoredProfile(**kwargs)
            return spatial_model
        elif profile_name.startswith("spatialmap"):
            self.extension = np.nan
            spatial_model = SpatialMap(file=self.get_profile_value("filename", profile))
            return spatial_model 
        else:
            raise Exception("Unrecognized spatial profile: %s"%profile)

        # Make sure the extension is within lower limit
        extension = self.get_extension(profile)
        if np.isnan(extension):   return center
        else:                     spatial_model.setp('sigma',extension)

        # Monkey patch some parameters for writing FITS file
        def effective_edge(energy=None): return truncate
        spatial_model.effective_edge = effective_edge
        def template_diameter(): return 4*truncate
        spatial_model.template_diameter = template_diameter
        def template_npix(): return npix
        spatial_model.template_npix = template_npix
        return spatial_model

    def get_source_model(self, profile='point',npix=500,truncate=0.5,galactic=True,**kwargs):
        spectral_model = self.get_spectral_model(**kwargs)
        spatial_model = self.get_spatial_model(profile,truncate,npix)
        if isinstance(spatial_model,SpatialModel):
            if not galactic:
                spatial_model.change_coordsystem(SkyDir.EQUATORIAL)
            else:   
                spatial_model.change_coordsystem(SkyDir.GALACTIC)
            s = ExtendedSource(
                name=self.codename,
                model=spectral_model,
                spatial_model=spatial_model)
        else:
            s = PointSource(
                name=self.codename,
                model=spectral_model,
                skydir=spatial_model)
        s.model.free[1] = False
        return s


    def toTAB(self, **kwargs):
        s = "%(target)s & %(glon)s & %(glat)s & %(dist)s & %(jfactor)s & %(alpha)s & %(notes)s \\\\"
        refs = ""
        for k in sorted(self.refs,reverse=True):
            refs += "[%s|%s], "%(k,self.refs[k])
            refs.strip(", ")
        jfactor = "$%.1f \pm %.2f$"%(np.log10(self.get_jvalue()),self.get_jsigma())
        if np.isnan(self.get_jvalue()): jfactor = "--"
        params = dict(target=self.name,
                      distance="%s \+/- %s"%(self.distance,self.dsigma),
                      glon="%.1f"%self.l(),
                      glat="%.1f"%self.b(),
                      dist="%.0f"%self.distance,
                      jfactor=jfactor,
                      alpha="--",
                      notes="--",
                      )
        return s%params
        

class TargetSet(Set):
    """ Subclass of Set that allows access by name, codename,
    abbreviation, or nickname."""

    def __init__(self, *items):
        super(TargetSet,self).__init__(*items)

    def add(self, *items):
        """ Add one or more items to the collection. 
        """
        super(TargetSet,self).add(*items)
        for item in self.list:
            self.dict[item.codename] = item
        self._update_names()

    def pop(self, key):
        """ Pop an item from the collection. 
        """
        item = self.dict[key]
        super(TargetSet,self).pop(item.name)
        self.dict.pop(item.codename)
        self._update_names()
        return item

    def _update_names(self):
        codenames,nicknames,abbreviations = [],[],[]
        for item in self.list:
            codenames.append(item.codename)
            nicknames.append(item.nickname)
            abbreviations.append(item.abbreviation)
        self.codenames=codenames
        self.nicknames=nicknames
        self.abbreviations=abbreviations

    def skydirs(self):
        return [target.skydir for target in self.list]

    def glons(self):
        return [target.glon() for target in self.list]

    def glats(self):
        return [target.glat() for target in self.list]

    def set_profile(self, profile):
        for target in self.list: target.profile = profile
    def force_extension(self, extension):
        for target in self.list: target.force_extension(extension)
        print 'INFO: setting custom extension to targets %s'%str(extension)

    def randomize(self, catalog, **kwargs):
        """ Choose random positions for glat and glon 
        avoiding sources in the input catalog.        """
        randloc = random_skydir(catalog,
                                size=len(self),
                                **kwargs)
        newset = copy.deepcopy(self)
        for i,loc in enumerate(randloc):
            skydir = SkyDir(loc[0], loc[1], SkyDir.GALACTIC)
            newset[i].skydir = skydir
        return newset

    def slide(self, catalog, **kwargs):
        """ Hold glat constant and slide to random glon
        avoiding sources in the input catalog.        """
        glats = self.glats()
        newloc = random_glon(catalog, glats, **kwargs)
        newset = copy.deepcopy(self)
        for i,loc in enumerate(newloc):
            skydir = SkyDir(loc[0], loc[1], SkyDir.GALACTIC)
            newset[i].skydir = skydir
        return newset

    def randomize_jvalues(self, profile=None, seed=0, fntype='lgauss', rescale=False):
        """ Randomize the nominal J-factor in accord with the 
        unconstrained ensemble...
        http://indico.cern.ch/getFile.py/access?contribId=0&resId=0&materialId=2&confId=126652

        fntype : Form of the likleihood for J.  This will be used to
        determine the form of the distribution from which randomized J
        values will be drawn.
        
        """
        np.random.seed(seed)
        newset = copy.deepcopy(self)
        for i,target in enumerate(newset.list):
            jvalue = target.get_jvalue(profile)
            jsigma = target.get_jsigma(profile)

            if fntype=='lgauss':            
                scale = 10**( np.log(10)*jsigma**2 )
                mean = jvalue
                if rescale: mean *= scale                
                rjvalue = np.random.lognormal(mean=np.log(mean),sigma=jsigma*np.log(10))
            elif fntype=='lgauss_like':
                rjvalue = np.random.lognormal(mean=np.log(jvalue),sigma=jsigma*np.log(10))
            else:
                raise Exception('Unknown fntype: ' + fntype)


            #rjvalue = 10**np.random.normal(np.log10(jvalue),jsigma)

            newset[i].set_jvalue(rjvalue,profile)

        return newset

    def loadYAML(self, filename, ttype=Target):
        """ Load a TargetSet from a yaml file.
        """
        yamldict = yaml.load(open(filename))
        self.loadDICT(yamldict,ttype)

    def loadDICT(self, dict, ttype=Target):
        """ Load a TargetSet from a python dictionary.
        """
        for dwarf in sorted(dict.keys()):
            self.add(ttype(**dict[dwarf]))

    def writeREG(self,filename,coordsys='galactic',ds9args=''):
        #for the items in the set, create a regionstring
        f = open(fname,'w')
        if len(ds9args)!=0:
            f.write(ds9args+'\n')
        for target in self.list:
            f.write(item.make_regionstring(coordsys=coordsys)+'\n')
        f.close()
        print '*INFO* created regionfile %s'%fname

    def writeXML(self, filename=None):
        if filename is None: filename="%s.xml"%self.name
        point_sources, extended_sources = [],[]
        for target in self.list:
            model = target.get_source_model()
            if isinstance(model,PointSource):
                point_sources.append(model)
            elif isinstance(model,ExtendedSource):
                extended_sources.append(model)
            else:
                raise Exception("Unrecognized target model: %s"%(type(model)))
        write_sources(point_sources, extended_sources, 
                      filename, convert_extended=True, expand_env_vars=True)

    def writeTAB(self, filename):
        """ For making a confluence table"""
        f = open(filename,'w')
        s = "|%(target)s|%(distance)s|%(glon)s|%(glat)s|%(ra)s|%(dec)s|%(refs)s|"
        params = dict(target="| TARGET ",
                      distance="| DISTANCE (kpc) ",
                      glon="| GLON (deg.) ",
                      glat="| GLAT (deg.) ",
                      ra="| RA (deg.) ",
                      dec="| DEC (deg.) ",
                      refs="| REFERENCE |"
                      )
        f.write(s%params+'\n')

        for target in self.list:
            refs = ""
            for k in sorted(target.refs,reverse=True):
                refs += "[%s|%s], "%(k,target.refs[k])
            refs.strip(", ")
            
            params = dict(target=target.name,
                          distance="%s \+/- %s"%(target.distance,target.dsigma),
                          glon="%.3f"%target.l(),
                          glat="%.3f"%target.b(),
                          ra="%.3f"%target.ra(),
                          dec="%.3f"%target.dec(),
                          refs=refs
                          )
            f.write(s%params+'\n')

    def writeYAML(self, filename):
        odict = {}
        for target in self.list:
            odict[target.codename] = target.toDict() 
        yaml_dump(odict,filename)



targetdir=join(os.environ['DSPHSROOT'],'targets')

def load_srcs_from_xml(xmlfile):

    et = ElementTree.ElementTree(file=xmlfile)
    tree = et.getroot()
    srcs = {}

    for c in tree.findall('source'):

        s = {}

        s['name'] = c.attrib['name']
        s['spectrum_type'] = c.find('spectrum').attrib['type']
        s['spatialModel_type'] = c.find('spatialModel').attrib['type']
        if s['spatialModel_type'] == 'SkyDirFunction': continue

        for p in c.findall('spatialModel/parameter'):
            s[p.attrib['name']] = {}            
            for k, v in p.attrib.iteritems():
                s[p.attrib['name']][k] = v

        for p in c.findall('spectrum/parameter'):
            s[p.attrib['name']] = {}            
            for k, v in p.attrib.iteritems():
                s[p.attrib['name']][k] = v
        
        srcs[s['name']] = s
    
    return srcs

def load_targets_from_config(config):
    if isinstance(config,basestring):
        config = yaml.load(open(config))
    dirname  = config.get('targetdir', None)
    filename = config.get('targetfile', None)
    default  = config.get('targetdefault',None)
    return load_targets(filename=filename,dirname=dirname,default=default)

def load_targets(filename=None, dirname=None, default=None):
    """
    Load a TargetSet from a filename and directory. Merge J-factors
    with default file if provided.

    filename: Input yaml file containing target set
    dirname : Directory of file
    default : Input yaml file containing default target set
    """
    if filename is None: filename = 'dwarfs.yaml'
    
    if dirname is None:  
        if not os.path.isabs(filename):
            dirname = targetdir
            filename = join(dirname,filename)
    targetdict = yaml.load(open(filename))

    if default: 
        if not os.path.isabs(default):
            default = join(dirname,default)
        defaultdict = yaml.load(open(default))
        # ADW: Only join J-factors for now...
        print "Merging with target list: %s"%default
        #targetdict = merge_target_dicts(defaultdict,targetdict)
        targetdict = merge_target_dicts(defaultdict,targetdict,keys=['profiles'])
    
    targets = TargetSet()
    targets.loadDICT(targetdict, DMTarget)
    return targets

def merge_target_dicts(d,u,keys=None):
    """ Merge two target dicts.
    """
    # ADW: Explicitly update nested dicts
    # http://stackoverflow.com/questions/3232943/
    def update(d, u):
        for k, v in u.iteritems():
            if isinstance(v, collections.Mapping):
                r = update(d.get(k, {}), v)
                d[k] = r
            else:
                d[k] = u[k]
        return d

    ret = dict(d)
    for name,target in u.items():
        print "Updating target: %s"%name
        if keys is None:
            if name not in ret:
                ret[name] = target
            else:
                update(ret[name],target)
        else:
            # Should throw exception if key doesn't exist
            for key in keys:
                update(ret[name][key],target[key])

    return ret

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()

