#!/usr/bin/env python
"""
Covienent function for IRF interface.
Adapted from plot_irfs.py by Jim Chiang.
"""

import pyIrfLoader
from dsphs.pointlike.utilities.quantile import Quantile
from scipy.integrate import quad
import numpy as np

pyIrfLoader.Loader_go()

class Irf(object):
    _factory = pyIrfLoader.IrfsFactory_instance()
    def __init__(self, irfsName, inc=0, phi=0):
        self._irfs = self._factory.create(irfsName)
        self._inc = inc
        self._phi = phi

class Psf(Irf):
    def __init__(self, irfsName, energy=100., inc=0, phi=0):
        Irf.__init__(self, irfsName, inc, phi)
        self._psf = self._irfs.psf()
        self._energy = energy
    def __call__(self, separation):
        psf, energy, inc, phi = self._psf, self._energy, self._inc, self._phi
        try:
            y = []
            for x in separation:
                y.append(psf.value(x, energy, inc, phi))
            return np.array(y)
        except TypeError:
            return psf.value(separation, energy, inc, phi)

    def quantile(self, q):
        integrand = lambda r: self(r)*2*np.pi*np.radians(r)
        quantile = Quantile(integrand, 0, np.inf)
        return quantile(q)

class Aeff(Irf):        
    def __init__(self, irfsName, inc=0, phi=0):
        Irf.__init__(self, irfsName, inc, phi)
        self._aeff = self._irfs.aeff()
    def __call__(self, energy):
        aeff, inc, phi = self._aeff, self._inc, self._phi
        try:
            y = []
            for x in energy:
                y.append(aeff.value(x, inc, phi))
            return np.array(y)
        except TypeError:
            return aeff.value(energy, inc, phi)

if __name__ == "__main__":
    from optparse import OptionParser
    usage = "Usage: %prog  [options] input"
    description = "python script"
    parser = OptionParser(usage=usage,description=description)
    (opts, args) = parser.parse_args()
    
    psf = Psf("P7REP_CLEAN_V10",energy=1e3)
    psf.quantile(.68)
    psf.quantile(.95)
