#!/usr/bin/env python
import numpy as np
from os.path import splitext

import inspect, collections
from scipy.optimize import fmin

from dsphs.utils.set import SetItem, Set
import dsphs.utils.tools as tools

from dsphs.base.lnlfn import LnLFn, ProfileLimit, ProfileLnL, MarginalLnL, BayesianLimit

class SEDBin(SetItem):
    def __init__(self,name,emin,emax,**kwargs):
        super(SEDBin,self).__init__(name)
        self.nobs = np.nan
        self.binTs = np.nan
        self.emin,self.emax = emin,emax
        self.pulimit,self.bulimit = np.nan, np.nan
        self.pllimit,self.bllimit = np.nan, np.nan
        self.__dict__.update(**kwargs)
        self.update_lnlfn()
        self.update_ulimit()
        self.update_llimit()
        self.update_mle()
    def get_value(self, value):
        return self.__dict__[value]

    def update_ulimit(self,alpha=0.05,ulimit='pulimit'):
        if ulimit=='pulimit':
            self.update_pulimit(alpha)
        elif ulimit=='bulimit':
            self.update_bulimit(alpha)
        else:
            raise ValueError('Limit not implemented: %s'%ulimit)

    def update_pulimit(self,alpha=0.05):
        try:
            self.pulimit = ProfileLimit(self.eflux,self.logLike).getLimit(alpha)
        except ValueError as e:
            print "WARNING: Unable to compute upper limit: ", e
            self.pulimit = np.nan

    def update_bulimit(self,alpha=0.05):
        try:
            self.bulimit = BayesianLimit(self.eflux,self.logLike).getLimit(alpha)
        except ValueError as e:
            print "WARNING: Unable to compute upper limit: ", e
            self.bulimit = np.nan

    def update_llimit(self,alpha=0.05,llimit='pllimit'):
        if llimit=='pllimit':
            self.update_pllimit(alpha)
        else:
            raise ValueError('Limit not implemented: %s'%llimit)

    def update_pllimit(self,alpha=0.05):
        try:
            self.pllimit = ProfileLimit(self.eflux,self.logLike).getLowerLimit(alpha)
        except ValueError as e:
            print "WARNING: Unable to compute lower limit: ", e
            self.pllimit = np.nan

    def update_lnlfn(self):
        self.lnlfn = LnLFn(self.eflux,self.logLike)

    def update_mle(self):
        self.mle = self.lnlfn.mle()

    def todict(self):
        return dict( emin=self.emin, emax=self.emax,
                     flux=self.flux, eflux=self.eflux,
                     logLike=self.logLike,
                     eflux2npred=self.eflux2npred,
                     nobs = self.nobs, binTs = self.binTs
                     )

    def __str__(self):
        np.set_printoptions(precision=3)
        s = "SED bin (%s):\n"%self.name
        s+= "  energy: [%g, %g]\n"%(self.emin,self.emax)
        s+= "  eflux:   %s\n"%str(self.eflux)
        s+= "  logLike: %s"%str(self.logLike)
        np.set_printoptions()
        return s

class SED(Set):
    def __init__(self, sed):
        bins = []
        if isinstance(sed, basestring):
            # Doesn't work yet...
            ext = splitext(sed)[-1]
            if ext == '.yaml':  self.__init__(SED.read_yaml(sed))
            elif ext == '.txt': self.__init__(SED.read_text(sed))
            elif ext == '.fits': self.__init__(SED.read_fits(sed))
            else: raise TypeError("Unrecognized file type...%s"%sed)
        elif isinstance(sed, collections.Iterable):
            for i,s in enumerate(sed):
                if isinstance(s,SEDBin): bins.append(s)
                elif isinstance(s,dict): bins.append(SEDBin('bin%i'%i,**s))
                else: raise TypeError("Unrecognized bin type...%s"%type(bin))
            super(SED,self).__init__(*bins)

    def get_values(self,value):
        return np.asarray([ item.get_value(value) for item in self.vals()])

    def get_params(self):
        return self.get_values('param')

    def get_emins(self):
        return self.get_values('emin')

    def get_emaxs(self):
        return self.get_values('emax')

    def get_energies(self):
        return np.sqrt(self.get_emins()*self.get_emaxs())

    def get_efluxes(self):
        return self.get_values('eflux')
        
    def get_logLikes(self):
        return self.get_values('logLike')

    def get_eflux2npred(self):
        return self.get_values('eflux2npred')
    
    def get_ulimits(self,alpha=0.05,ulimit='pulimit'):
        for b in self.vals(): b.update_ulimit(alpha,ulimit)
        return self.get_values(ulimit)

    def get_llimits(self,alpha=0.05,llimit='pllimit'):
        for b in self.vals(): b.update_llimit(alpha,llimit)
        return self.get_values(llimit)

    def get_mles(self):
        return self.get_values('mle')

    def get_global_bin_logLike(self, x, spectrum, jsigma = None):

        emin,emax = self.get_emins(),self.get_emaxs()
        energy = np.sqrt(emin*emax)
        eflux,logLike = self.get_efluxes(), self.get_logLikes()
        limits = self.get_ulimits(alpha=0.05)

        x = np.logspace(-3,5,250)
        xmin = x.min(); xmax = x.max()

        start_flux = min(limits/energy)
        spectrum.set_flux(start_flux,emin=emin[0],emax=emax[-1])
     
        flux = start_flux * x
        norm = spectrum[0] * x
        
        # Predicted number of counts
        pred = np.array([spectrum.i_flux(_emin,_emax,e_weight=1) for _emin,_emax in zip(emin,emax)])
        # Create a list of likelihoods
        likes = [ LnLFn(f,l-max(l)) for f,l in zip(eflux,logLike) ]
        # Because of the way the spline works add epsilon
        eflux_min = eflux.min() * ( 1 + 1e-4) 
        like = lambda x: np.array([lnlfn(x*p) if x*p > eflux_min else lnlfn(eflux_min) for lnlfn,p in zip(likes,pred)])
        return like

    def get_global_logLike(self, spectrum, jsigma = None, method='profile',fntype='lgauss'):
        """ Get the global log-likelihood derived from the input spectrum
        spectrum : Input pointlike Model-class object
        jsigma   : Uncertainty on j-factor
        returns  :
           (norm, flux, lnl, plnl, flnl)
        """
        
        emin,emax = self.get_emins(),self.get_emaxs()
        energy = np.sqrt(emin*emax)
        eflux,logLike = self.get_efluxes(), self.get_logLikes()
        limits = self.get_ulimits(alpha=0.05)

        x = np.logspace(-3,5,250)
        xmin = x.min(); xmax = x.max()

        start_flux = min(limits/energy)
        spectrum.set_flux(start_flux,emin=emin[0],emax=emax[-1])
     
        flux = start_flux * x
        norm = spectrum[0] * x
        
        # Predicted energy flux
        pred = np.array([spectrum.i_flux(_emin,_emax,e_weight=1) for _emin,_emax in zip(emin,emax)])
        # Create a list of likelihoods
        likes = [ LnLFn(f,l-max(l)) for f,l in zip(eflux,logLike) ]
        # Because of the way the spline works add epsilon
        eflux_min = eflux.min() * ( 1 + 1e-4) 
        like = lambda x: sum([lnlfn(x*p) if x*p > eflux_min else lnlfn(eflux_min) for lnlfn,p in zip(likes,pred)])
        # Global log-likelihood using nominal value of nuisance parameter
        lnl = np.array([like(_x) for _x in x])
        
        # Log-likelihood after profiling over nuisance parameter
        if jsigma:

            if method == 'profile':            
                plnl = ProfileLnL(x,lnl,sigma=jsigma,fntype=fntype)(x)
            elif method == 'marginalize':
                plnl = MarginalLnL(x,lnl,sigma=jsigma,fntype=fntype)(x)                
        else:
            plnl = lnl
            
        # Log-likelihood as a function of flux rather than normParam
        flnl = lnl

        return norm, flux, lnl, plnl, flnl

    def get_global_limit(self, spectrum, jsigma = None, alpha=0.05):
        """ Return limit from global likelihood fit """
        norm, flux, lnl, plnl, flnl = self.get_global_logLike(spectrum,jsigma)
        if jsigma is None: return ProfileLimit(norm,lnl).getLimit(alpha)
        else:              return ProfileLimit(norm,plnl).getLimit(alpha)

    def get_global_mle(self, spectrum, jsigma = None):
        """ return factor to multiply spectrum by """
        norm, flux, lnl, plnl, flnl = self.get_global_logLike(spectrum,jsigma)
        if jsigma is None: return LnLFn(norm,lnl).mle()
        else:              return LnLFn(norm,plnl).mle()

    def get_bin_limit(self,spectrum, jsigma = None, alpha=0.05):
        """ Return limit from bumping into first bin. """
        emin,emax = self.get_emins(),self.get_emaxs()
        energy = np.sqrt(emin*emax)
        eflux,logLike = self.get_efluxes(), self.get_logLikes()
        limits = self.get_ulimits(alpha=alpha)

        start_flux = min(limits/energy)
        spectrum.set_flux(start_flux,emin=emin[0],emax=emax[-1])

        # Predicted number of counts
        pred = np.array([spectrum.i_flux(_emin,_emax,e_weight=1) for _emin,_emax in zip(emin,emax)])

        fn = lambda x: sum(np.where(x*pred<limits, limits - x*pred, np.inf))
        x = fmin(fn, 1, disp=False)

        ulimit = spectrum[0] * x
        return float(ulimit)

    def tolist(self):
        l = []
        for i in self.vals(): l.append(i.todict())
        return l

    @classmethod
    def read_yaml(cls,filename):
        """ Read an SED object from a yaml file """
        sed = tools.yaml_load(filename)
        return SED(sed)
    
    def write_yaml(self, filename):
        """ Write an SED object to a yaml file """
        #if not isinstance(cls, SED): raise TypeError("Not SED object")
        tools.yaml_dump(self.tolist(),filename)

    @classmethod
    def read_text(cls,filename):
        """ Read an SED object from a plain text file. """
        data = np.loadtxt(filename,unpack=True)
        emins, emaxs = np.unique(data[0]),np.unique(data[1])
        eflux = data[2].reshape(len(emins),-1)
        logLike = data[3].reshape(len(emins),-1)
        sed = [ dict(emin=emin,emax=emax,logLike=logLike[i],eflux=eflux[i])
                for i,(emin,emax) in enumerate(zip(emins,emaxs)) ]
        return SED(sed)

    def write_text(self, filename):
        """ Write an SED object to a plain text file. """
        f = open(filename,'w')
        header = '### %s\n'%filename
        #header += ''.join('%-15s'%s for s in 
        #                  ['#[0]emin(MeV)','[1]emax(MeV)','[2]eflux(MeV/cm^2/s)','[3]logLikelihood\n'])
        header += '#[0]emin(MeV) [1]emax(MeV) [2]eflux(MeV/cm^2/s) [3]deltaLogLike\n'
        f.write(header)

        logLikes = self.get_logLikes()
        efluxes = self.get_efluxes()
        emins = self.get_emins().repeat(logLikes.shape[1])
        emaxs = self.get_emaxs().repeat(logLikes.shape[1])
        data = np.vstack([emins,emaxs,efluxes.ravel(),logLikes.ravel()]).T
        np.savetxt(f,data,fmt='%-15.9g')
        f.close()

    @classmethod
    def read_fits(cls, filename):
        """ Read an SED object from a FITS file. """
        import pyfits
        f = pyfits.open(filename)
        logLike = f[0].data
        emins = f['ENERGIES'].data['EMIN']
        emaxs = f['ENERGIES'].data['EMAX']
        eflux = f['EFLUX'].data['EFLUX']
        sed = [ dict(emin=emin,emax=emax,logLike=logLike[i],eflux=eflux) 
                for i,(emin,emax) in enumerate(zip(emins,emaxs)) ]
        f.close()
        return SED(sed)

    def write_fits(self, filename):
        """ Write an SED object to a FITS file. """
        import pyfits
        like_hdu = pyfits.PrimaryHDU(self.get_logLikes())

        columns = [
            pyfits.Column(name='EMIN',format='D',
                          unit='MeV',array=self.get_emins()),
            pyfits.Column(name='EMAX',format='D',
                          unit='MeV',array=self.get_emaxs())
        ]
        energy_hdu = pyfits.new_table(columns)
        energy_hdu.name = 'ENERGIES'

        columns = [
            pyfits.Column(name='EFLUX',format='D',
                          unit='MeV/cm**2/s**2',array=self.get_efluxes()[0])
        ]
        eflux_hdu = pyfits.new_table(columns)
        eflux_hdu.name = 'EFLUX'
        
        hdulist = pyfits.HDUList([like_hdu,energy_hdu,eflux_hdu])
        hdulist.writeto(filename,clobber=True)
        
if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()

