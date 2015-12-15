#!/usr/bin/env python
# $Header: /nfs/slac/g/glast/ground/cvs/users/kadrlica/dsphs/like/lnlfn.py,v 1.17 2014/11/29 05:01:10 mdwood Exp $

"""
@author Matthew Wood <mdwood@slac.stanford.edu>
@author Alex Drlica-Wagner <kadrlica@slac.stanford.edu>
"""

__author__   = "Matthew Wood <mdwood@slac.stanford.edu>"
__date__     = "$Date: 2014/11/29 05:01:10 $"

import sys
import scipy as sp
import numpy as np
import math
import scipy.optimize as opt
from scipy.interpolate import UnivariateSpline
import scipy.special as spf
from scipy.integrate import quad
import scipy.stats as stats
import copy

uniform = stats.uniform
norm = stats.norm
lognorm = stats.lognorm

"""
A note on the highly confusing scipy.stats.lognorm function...
The three inputs to this function are:
s           : This is the variance of the underlying gaussian distribution
scale = 1.0 : This is the mean of the linear-space lognormal distribution.
              The mean of the underlying normal distribution occurs at ln(scale)
loc = 0     : This linearly shifts the distribution in x (DO NOT USE)

The convention is different for numpy.random.lognormal
mean        : This is the mean of the underlying normal distribution (so mean = log(scale))
sigma       : This is the standard deviation of the underlying normal distribution (so sigma = s)

For random sampling:
numpy.random.lognormal(mean, sigma, size)
mean        : This is the mean of the underlying normal distribution (so mean = exp(scale))
sigma       : This is the standard deviation of the underlying normal distribution (so sigma = s)

scipy.stats.lognorm.rvs(s, scale, loc, size)
s           : This is the standard deviation of the underlying normal distribution
scale       : This is the mean of the generated random sample scale = exp(mean)

Remember, pdf in log space is

plot( log(x), stats.lognorm(sigma,scale=exp(mean)).pdf(x)*x )
"""

def norm(x,mu,sigma=1.0):
    """ Scipy norm function """
    return stats.norm(loc=mu,scale=sigma).pdf(x)

def ln_norm(x,mu,sigma=1.0):
    """ Natural log of scipy norm function truncated at zero """
    return np.log(stats.norm(loc=mu,scale=sigma).pdf(x))

def lognorm(x,mu,sigma=1.0):
    """ Log-normal function from scipy """
    return stats.lognorm(sigma,scale=mu).pdf(x)

def log10norm(x,mu,sigma=1.0):
    """ Scale scipy lognorm from natural log to base 10 
    x     : input parameter
    mu    : mean of the underlying log10 gaussian
    sigma : variance of underlying log10 gaussian
    """
    return stats.lognorm(sigma*np.log(10),scale=mu).pdf(x)

def ln_log10norm(x,mu,sigma=1.0):
    """ Natural log of base 10 lognormal """
    return np.log(stats.lognorm(sigma*np.log(10),scale=mu).pdf(x))

def gauss(x,mu,sigma=1.0):
    s2 = sigma*sigma
    return 1./np.sqrt(2*s2*np.pi)*np.exp(-(x-mu)*(x-mu)/(2*s2))

def lngauss(x,mu,sigma=1.0):
    s2 = sigma*sigma
    return -0.5*np.log(2*s2*np.pi) - np.power(x-mu,2)/(2*s2)

def lgauss(x,mu,sigma=1.0,logpdf=False):
    """ Log10 normal distribution...

    x     : Parameter of interest for scanning the pdf
    mu    : Peak of the lognormal distribution (mean of the underlying normal distribution is log10(mu)
    sigma : Standard deviation of the underlying normal distribution
    """
    x = np.array(x,ndmin=1)

    lmu = np.log10(mu)
    s2 = sigma*sigma
 
    lx = np.zeros(x.shape)
    v = np.zeros(x.shape)

    lx[x>0] = np.log10(x[x>0])

    v = 1./np.sqrt(2*s2*np.pi)*np.exp(-(lx-lmu)**2/(2*s2))

    if not logpdf: v /= (x*np.log(10.))

    v[x<=0] = -np.inf

    return v

def lnlgauss(x,mu,sigma=1.0,logpdf=False):

    x = np.array(x,ndmin=1)

    lmu = np.log10(mu)
    s2 = sigma*sigma

    lx = np.zeros(x.shape)
    v = np.zeros(x.shape)

    lx[x>0] = np.log10(x[x>0])

    v = -0.5*np.log(2*s2*np.pi) - np.power(lx-lmu,2)/(2*s2) 
    if not logpdf: v -= 2.302585*lx + np.log(np.log(10.))

    v[x<=0] = -np.inf

    return v

# Possible correction factor:
# ( np.exp( 0.5*(np.log(10)*sigma)**2)**2 )

class LnLFn(object):
    """
    Helper class for interpolating a 1-D log-likelihood function from a
    set of tabulated values.  
    """
    def __init__(self,x,y):

        msk = np.isfinite(y)
        x = np.array(x[msk],ndmin=1)
        y = np.array(y[msk],ndmin=1)

        self._xmin = x[0]
        self._xmax = x[-1]
        self._ymin = y[0]
        self._ymax = y[-1]
        self._dydx_lo = (y[1]-y[0])/(x[1]-x[0])
        self._dydx_hi = (y[-1]-y[-2])/(x[-1]-x[-2])

        self._fn = UnivariateSpline(x,y,s=0,k=2)

    @property
    def xmin(self):
        return self._xmin

    @property
    def xmax(self):
        return self._xmax

    def __call__(self,x):

        xnew = np.asarray(x)
        nx = xnew.ndim

        below_bounds = xnew < self._xmin
        above_bounds = xnew > self._xmax

        # Handle case where x is a scalar
        if xnew.ndim == 0:

            if xnew > self._xmax:
                return self._ymax + (x-self._xmax)*self._dydx_hi
            elif xnew < self._xmin:
                return self._ymin + (x-self._xmin)*self._dydx_lo
            else:
                return self._fn(xnew)
        else:

            dxhi = np.asarray(xnew-self._xmax)
            dxlo = np.asarray(xnew-self._xmin)

            # UnivariateSpline will only accept 1-D arrays so this
            # passes a flattened version of the array.
            y = self._fn(xnew.ravel())
            y.resize(xnew.shape)
            
            # If y is a rank-0 array convert it to rank-1
            if y.ndim == 0:
                y = y[np.newaxis,...]

            y[above_bounds] = (self._ymax + dxhi[above_bounds]*self._dydx_hi)
            y[below_bounds] = (self._ymin + dxlo[below_bounds]*self._dydx_lo)
            return y

    def mle(self):
        """ Maximum likelihood estimator """
        x0 = self._xmax/4.
        ret = opt.fmin(lambda x: np.where(self._xmax>x>0, -self(x), np.inf)[0], 
                       x0, disp=False)

        mle = float(ret[0])
        return mle

    def fn_mle(self):
        return self(self.mle())

class LnLFnBase(object):
    """
    Base class for the evaluation of likelihood of the form:

    L(x,y|z') = L_z(x*y|z')*L_y(y)

    where x is the parameter of interest, y is a nuisance parameter,
    and L_z is a likelihood constraining z = x*y.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.

    fntype : string 
  
       Function with which to model the nuisance parameter.  Valid options are:

       gauss         : normal distribution with pdf truncated at zero.
       lgauss        : log-normal posterior
       lgauss_like   : log-normal likelihood
       lgauss_logpdf : log-normal posterior in the basis of the log of the nuisance

    sigma : float
        Width parameter of the nuisance parameter distribution.
    
    """
    def __init__(self,lnlx,lnly,sigma=0.1,fntype='lgauss'):

        lnlx = np.array(lnlx,ndmin=1)
        lnly = np.array(lnly,ndmin=1)

        self._xmin = min(lnlx) #lnlx[0]
        self._xmax = max(lnlx) #lnlx[-1]
        self._sigma=sigma
        self._fntype=fntype
        self._pdfnorm=1
        
        if len(lnlx) < 2 or len(lnly) < 2:
            raise ValueError("")
        if isinstance(fntype, basestring):

            # Posterior
            if fntype == 'lgauss':
                self._ypdf = lgauss
                self._ylnpdf = lnlgauss
            elif fntype == 'lgauss_like':
                self._ypdf = lambda x, y, s: lgauss(y,x,s)
                self._ylnpdf = lambda x, y, s: lnlgauss(y,x,s)
            elif fntype == 'lgauss_logpdf':
                self._ypdf = lambda x, y, s: lgauss(x,y,s,logpdf=True)
                self._ylnpdf = lambda x, y, s: lnlgauss(x,y,s,logpdf=True)
            elif fntype == 'gauss':
                self._ypdf = gauss
                self._ylnpdf = lngauss
                self._pdfnorm = 1-0.5*(1+spf.erf(-1.0/(math.sqrt(2.0)*self._sigma)))
            elif fntype == 'norm':
                self._ypdf = norm
                self._ylnpdf = ln_norm
                self._pdfnorm = quad(norm, 0, np.inf, args=(1.0,self._sigma) )[0]
            else:
                raise ValueError("fntype = %s is not supported." %(fntype))
        else:
            self._ypdf = fntype
            self._ylnpdf = lambda x, *args, **kwargs: np.log(fntype(x,*args,**kwargs))
            self._pdfnorm = quad(self._ypdf, 0, np.inf, args=(1.0, self._sigma) )[0]

        self._fn_lnl = LnLFn(lnlx,lnly)

    def like(self,x,y):
        """Evaluate the 2-D likelihood in the x/y parameter space.
        The dimension of the two input arrays should be the same.

        Parameters
        ----------
        x : array_like
        Array of coordinates in the `x` parameter.
        
        y : array_like       
        Array of coordinates in the `y` nuisance parameter.
        """
        
        z = self._fn_lnl(x*y)
        return np.exp(z)*self._ypdf(y,1.0,self._sigma)/self._pdfnorm

    def lnlike(self,x,y):

        return self._fn_lnl(x*y) + \
            self._ylnpdf(y,1.0,self._sigma) - np.log(self._pdfnorm)
        

class MarginalLnL(LnLFnBase):
    """
    MarginalLnL(lnlx,lnly,sigma=0.1,fntype='lgauss')

    Generate the marginal likelihood for the parameter x given a
    likelihood L(z) where z = x*y and y is a nuisance parameter.
    The marginalization is calculated by evaluating the 1-D integral:

    L(x) = \int L(x*y)P(y)dy

    where P(y) is the p.d.f. for y which is assumed to have the
    physical domain [0,inf].  The 1-D log-likelihood for z evaluated
    at y=y_{0} is passed as a set of tabulated values using the arrays
    `lnlx` and `lnly`.  

    This class returns a function that can be used to evaluate the
    marginal likelihood.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.

    fntype : string 
  
       Function with which to model the nuisance parameter.  Options
       are normal ('gauss') and log-normal ('lgauss').  In the case of
       the normal distribution the p.d.f. is truncated at zero.

    sigma : float
        Width parameter of the nuisance parameter distribution.

    Examples
    --------
    >>> import DMLimits as dmlim

    >>> lnlx = np.linspace(0,10,100)
    >>> lnly = lnlx*lnlx
    >>> mlnl = dmlim.MarginalLnL(lnlx,lnly,sigma=0.2,fntype='lgauss')
    >>> mlnly = mlnl(lnlx)


    """

    def __init__(self,lnlx,lnly,sigma=0.1,fntype='lgauss'):
        super(MarginalLnL,self).__init__(lnlx,lnly,sigma,fntype)

        logx = np.linspace(np.log(self._xmin),np.log(self._xmax),1000)
        x = np.exp(logx)
        logz = self.marginal_lnlike(x)
        self._mfn = UnivariateSpline(x,logz,s=0,k=2)

    def nuisancePDF(self,x):
        return self._ypdf(x,1.0,self._sigma)/self._pdfnorm

    def mle(self):
        """ Maximum likelihood estimator """

        xmax = self._fn_lnl._xmax
        x0 = self._fn_lnl.mle()
        ret = opt.fmin(lambda x: np.where(xmax>x>0, -self(x), np.inf)[0], 
                       x0, disp=False)
                       
        mle = float(ret[0])
        return mle
        
    def marginal_lnlike(self,x):

        yedge = np.logspace(-10*self._sigma,10.*self._sigma,1001)
        yw = yedge[1:]-yedge[:-1]
        yc = 0.5*(yedge[1:]+yedge[:-1])

        s = self.like(x[:,np.newaxis],yc[np.newaxis,:])
        z = np.sum(s*yw,axis=1)
        msk = z>0
        dlogzdx = (np.log(z[msk][-1]) - np.log(z[msk][-2]))/(x[msk][-1]-x[msk][-2])

        logz = np.zeros(z.shape)
        logz[msk] = np.log(z[msk])
        logz[~msk] = logz[msk][-1] + (x[~msk] - x[msk][-1])*dlogzdx

        return logz

    def __call__(self,x):
        """Evaluate the log-likelihood for the x parameter
        marginalizing over the nuisance parameter y."""

        x = np.array(x,ndmin=1)
        return self._mfn(x)

class ProfileLnL(LnLFnBase):
    """
    ProfileLnL(lnlx,lnly,sigma=0.1,fntype='lgauss')

    Profile likelihood for the parameter x given a
    likelihood L(z) where z = x*y and y is a nuisance parameter.
    The profile is calculated by evaluating the supremum:

    L(x) = \sup L(x*y)P(y)dy

    where P(y) is the p.d.f. for y which is assumed to have the
    physical domain [0,inf].  The 1-D log-likelihood for z evaluated
    at y=y_{0} is passed as a set of tabulated values using the arrays
    `lnlx` and `lnly`.  

    This class returns a function that can be used to evaluate the
    profile likelihood.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.

    fntype : string 
  
       Function with which to model the nuisance parameter.  Options
       are normal ('gauss') and log-normal ('lgauss').  In the case of
       the normal distribution the p.d.f. is truncated at zero.

    sigma : float
        Width parameter of the nuisance parameter distribution.

    Examples
    --------
    >>> import lnlfn

    >>> lnlx = np.linspace(0,10,100)
    >>> lnly = lnlx*lnlx
    >>> plnl = lnlfn.ProfileLnL(lnlx,lnly,sigma=0.2,fntype='lgauss')
    >>> plnly = plnl(lnlx)

    """

    def __init__(self,lnlx,lnly,sigma=0.1,fntype='lgauss'):
        super(ProfileLnL,self).__init__(lnlx,lnly,sigma,fntype)
    
    def mle(self):
        """ Maximum likelihood estimator """

        xmax = self._fn_lnl._xmax
        x0 = self._fn_lnl.mle()
        ret = opt.fmin(lambda x: np.where(xmax>x>0, -self(x)[1], np.inf)[0], x0, disp=False)
                       
        mle = float(ret[0])
        return mle

    def mle2(self):
        """ Maximum likelihood estimator """

        xmax = self._fn_lnl._xmax
        x0 = self._fn_lnl.mle()
        mle = opt.fmin(lambda x: -self.lnlike(x[0],x[1]), np.array([x0,1.0]), 
                       disp=False)
        return mle

    def profile_lnlike(self,x):
        """Evaluate the profiled log-likelihood.
           x : array-like
        """

        x = np.array(x,ndmin=1)

        z = []
        y = []
        
        for xtmp in x:
            
            fn = lambda t: -self.lnlike(xtmp,t)
            ytmp = opt.fmin(fn,1.0,disp=False)[0]
#            ytmp = opt.minimize_scalar(fn,bracket=(0.1,10.0)).x
            ztmp = self.lnlike(xtmp,ytmp)
            z.append(ztmp)
            y.append(ytmp)

        return np.array(y), np.array(z)
    
    def __call__(self,x):
        """Evaluate the profiled log-likelihood.
           x : array-like
        """
        return np.squeeze(self.profile_lnlike(x)[1])

class BayesianLimit(object):
    """
    BayesianLimit(lnlx,lnly,prior='uniform')

    Evaluate the upper limit on a parameter using Bayesian
    methodology.  The assumed prior can be configured with the `prior`
    option.  A 1-D log-likelihood curve is defined using the input
    arrays `lnlx` and `lnly`.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.
       
    prior : string

       Set the form of the prior.
      
       uniform : Flat prior w/ P(mu) = 1 for mu >=0
       poisson : Non-informative prior for poisson mean w/ P(mu) =
       1/sqrt(mu)

    Examples
    --------
    >>> import DMLimits as dmlim

    >>> lnlx = np.linspace(0,10,100)
    >>> lnly = lnlx*lnlx
    >>> blim = dmlim.BayesianLimit(lnlx,lnly)
    >>> print blim.getLimit(0.05)  
    """

    def __init__(self,lnlx,lnly,prior='uniform'):
        self._fn_lpdf = LnLFn(lnlx,lnly)
        self._fn_mle = self._fn_lpdf.mle()

        if prior == 'uniform':
            self._fn_pdf = lambda t: np.exp(self._fn_lpdf(t))
        elif prior == 'poisson':
            self._fn_pdf = lambda t: np.exp(self._fn_lpdf(t))*t**-0.5
        else:
            raise Exception('Unrecognized prior: ' + prior)

        dx = lnlx[1:]-lnlx[:-1]
        cdf = np.cumsum(np.concatenate(([0],self._fn_pdf(lnlx[:-1]+dx*0.5)*dx)))
        cdf /= cdf[-1]

        msk = cdf < 1-1E-10
        self._fn_cdf = UnivariateSpline(lnlx,cdf,s=0)        
        self._icdf = UnivariateSpline(cdf[msk],lnlx[msk],s=0,k=1)        


    def getLimit(self,alpha=0.05):
        """Evaluate the upper limit corresponding to a C.L. of (1-alpha)%.

        Parameters
        ----------
        alpha : Upper limit confidence level.
        """

        return float(self._icdf(1-alpha))

class ProfileLimit(object):
    """
    ProfileLimit(lnlx,lnly)

    Evaluate the upper limit on a parameter using the MINOS
    methodology.  A 1-D log-likelihood curve for the parameter is
    defined using the input arrays `lnlx` and `lnly`.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.

    Examples
    --------
    >>> import DMLimits as dmlim

    >>> lnlx = np.linspace(0,10,100)
    >>> lnly = lnlx*lnlx
    >>> plim = dmlim.ProfileLimit(lnlx,lnly)
    >>> print plim.getLimit(0.05)  
    """

    def __init__(self,lnlx,lnly):
        """ Initialize the profile likelihood function.
        Some complications having to do with tolerance
        and zeros """
        
        self._fn = LnLFn(lnlx,lnly)
        self._xmin = self._fn.xmin
        self._xmax = self._fn.xmax

        # Find xmin bounded by zero
        #ftmp = lambda x: np.where(x>0, -self._fn(x), np.inf)[0]
        ftmp = lambda x: np.where(x>self._xmin, -self._fn(x), np.inf)[0]
        # starting guess
        #x0 = np.sqrt(min(lnlx[lnlx > 0])*self._xmax)
        #self._mle = opt.fmin(ftmp, x0, disp=False, xtol=1e-10*self._xmin)[0]
        self._mle = self._fn.mle()
        #self._xmin = opt.fminbound(lambda x: -self._fn(x),lnlx[0],lnlx[-1])
        self._fmax = self._fn(self._mle)
        #print self._xmin, self._fmax

    def getLimit(self,alpha):
        """Evaluate the upper limit corresponding to a C.L. of (1-alpha)%.

        Parameters
        ----------
        alpha : Upper limit confidence level.
        """
        dlnl = pow(math.sqrt(2.)*spf.erfinv(1-2*alpha),2.)/2.  
        rf = lambda x: self._fn(x)+dlnl-self._fmax
        return opt.brentq(rf,self._mle,self._xmax,xtol=1e-10*self._mle)
    
    def getUpperLimit(self,alpha):
        return self.getLimit(alpha)

    def getLowerLimit(self,alpha):
        """Evaluate the lower limit corresponding to a C.L. of (1-alpha)%.

        Parameters
        ----------
        alpha : Lower limit confidence level.
        """
        dlnl = pow(math.sqrt(2.)*spf.erfinv(1-2*alpha),2.)/2.  
        rf = lambda x: self._fn(x)+dlnl-self._fmax
        return opt.brentq(rf,self._xmin,self._mle,xtol=1e-10*self._mle)

class StackedLnL(object):
    """
    StackLnL()
    
    Stack multiple profile log-likelihoods in order to calculate
    a combined log-likelihood curve.
    """
    def __init__(self):
        self.lnl_fn = []
        self.xlim = -1

    def addlnl(self,lnlx,lnly):
        """Add a likelihood profile to the stacked likelihood.

        Parameters
        ----------
        lnlx : array_like
           Array of points at which the log-likelihood is evaluated.
            
        lnly : array_like       
           Array of log-likelihood values evalated at the points in
           `lnlx`.
            """
        fn = LnLFn(lnlx,lnly)
        self.lnl_fn.append(fn)
        if self.xlim < 0:
            self.xlim = max(lnlx)
        elif max(lnlx) < self.xlim:
            self.xlim = max(lnlx)

    def eval(self,x):
        s = 0
        for fn in self.lnl_fn:
            s += fn(x)
        return s

    def __call__(self,x):
        return self.eval(x)

    def mle(self):
        """ Maximum likelihood estimator """
        x0 = self.xlim/4.
        ret = opt.fmin(lambda x: np.where(self.xlim>x>0, -self(x), np.inf)[0], x0, disp=False)
        mle = float(ret[0])
        return mle

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import scipy.stats
    import matplotlib as mpl

    def poisson_lnl(counts, model, nuisance=1.0):
        epsilon = 1e-3
        model = np.where(model <= 0, epsilon, model)/nuisance
        loglike = counts*np.log(model) - model
        return loglike - np.max(loglike)

    # Gaussian likelihood
    x0 = 1
    sigma = 0.5

    # Log-normal nuisance
    fntype = 'lgauss'
    jsigma = 0.32

    lnlx = np.linspace(0.0,20.0,200)
    lnly = -(lnlx-x0)**2

    #lnlx = np.linspace(0.0,20.0,200)
    #lnly = poisson_lnl(5,lnlx)

    mlnl = MarginalLnL(lnlx,lnly,sigma=jsigma,fntype=fntype)
    plnl = ProfileLnL(lnlx,lnly,sigma=jsigma,fntype=fntype)

    mlnly = np.array(mlnl(lnlx))
    plnly = np.array(plnl(lnlx))

    mlnl_fn = LnLFn(lnlx,mlnly)
    plnl_fn = LnLFn(lnlx,plnly)

    mlnl_xmax = opt.fminbound(lambda x: -mlnl_fn(x),lnlx[0],lnlx[-1])
    plnl_xmax = opt.fminbound(lambda x: -plnl_fn(x),lnlx[0],lnlx[-1])

    mlnl_ymax = mlnl_fn(mlnl_xmax)
    plnl_ymax = plnl_fn(plnl_xmax)

    print 'mlnl xmax ', mlnl_xmax
    print 'plnl xmax ', plnl_xmax

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel(r'$-\log\mathcal{L}$')
    ax.grid(True)

    ax.set_ylim(0,10)
    ax.set_xlim(0,10)

    plt.plot(lnlx, -lnly, 'r-', linewidth=2)
    plt.plot(lnlx, -mlnly+mlnl_ymax, 'b--', linewidth=2,label='Marginal LnL')
    plt.plot(lnlx, -plnly+plnl_ymax, 'g-.', linewidth=2,label='Profile LnL')

    ax.legend()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(True)

    xlo = 0.1
    xhi = 15.
    ylo = 0.01
    yhi = 4

    t1,t2 = sp.mgrid[xlo:xhi:.01, ylo:yhi:.01]

    #z = mlnl.like(t1,t2)
    #lnz = mlnl.lnlike(t1,t2)

    z = plnl.like(t1,t2)
    lnz = plnl.lnlike(t1,t2)

    x = np.linspace(xlo,xhi,100)
    y = [plnl.get_xmax(t) for t in x]

    plt.plot(x,y,linewidth=2,color='r',linestyle='--')

    cs = plt.contour(-lnz.transpose(), 
                     extent=[xlo,xhi,ylo,yhi],levels=[0.35,1.0,2.3,4.6,6.],
                     cmap=mpl.cm.jet, origin='lower',vmax=3)
    plt.clabel(cs, fontsize=9, inline=1)

    fig = plt.figure()

    y = np.linspace(0.1,3.0)

    z1 = [plnl.lnlike(2.0,t) for t in y]
    z2 = [plnl.lnlike(1.0,t) for t in y]
    z3 = [plnl.lnlike(0.5,t) for t in y]
    z4 = [plnl.lnlike(0.1,t) for t in y]

    plt.plot(y,z1)
    plt.plot(y,z2)
    plt.plot(y,z3)
    plt.plot(y,z4)

    print '95% C.L. Profile Limit:                 ', 
    print ProfileLimit(lnlx,lnly).getLimit(0.05)
    print '95% C.L. Profile Limit (Marginalized):  ', 
    print ProfileLimit(lnlx,mlnly).getLimit(0.05)
    print '95% C.L. Bayesian Limit:                ', 
    print BayesianLimit(lnlx,lnly).getLimit(0.05)
    print '95% C.L. Bayesian Limit (Marginalized): ', 
    print BayesianLimit(lnlx,mlnly).getLimit(0.05)


    #import numpy as np
    #import scipy.stats as stats
    # 
    #mean=1.0;    # Mean of underlying normal distribution (mu or ln(x_0))
    # 
    ## Random samples
    #gauss_sample = np.exp(np.random.normal(loc=mean,scale=jsigma,size=1e4))
    #scipy_sample = stats.lognorm.rvs(jsigma,loc=0,scale=np.exp(mean),size=1e4); 
    #numpy_sample = np.random.lognormal(mean=mean,sigma=jsigma,size=1e4)
    # 
    #print 'gauss:',stats.lognorm.fit(gauss_sample,floc=0)
    #print 'scipy:',stats.lognorm.fit(scipy_sample,floc=0)
    #print 'numpy:',stats.lognorm.fit(numpy_sample,floc=0)
    # 
    #plt.figure()
    #plt.hist(gauss_sample,bins=100,alpha=0.5, fc='k', label='gauss',normed=True);
    #plt.hist(scipy_sample,bins=100,alpha=0.5, fc='b', label='scipy',normed=True); 
    #plt.hist(numpy_sample,bins=100,alpha=0.5, fc='r', label='numpy',normed=True);
    #x = np.linspace(0,20,1e4)
    #plt.plot(x,stats.lognorm(jsigma, scale=np.exp(mean), loc=0).pdf(x),c='b',lw=2,label='scipy')
    #plt.legend()
    #
    # In log-space remember to multiply the pdf by x
    #plt.hist(np.log(scipy_sample),bins=100,alpha=0.5, fc='k', label='gauss',normed=True);
    #plt.plot(np.log(x),stats.lognorm(jsigma, scale=np.exp(mean)*x, loc=0).pdf(x),c='b',lw=2,label='scipy')   
    # 
    ## Do the same thing in log10 space
    #gauss_sample = 10**(np.random.normal(loc=mean,scale=jsigma,size=1e4))
    #scipy_sample = stats.lognorm.rvs(jsigma*np.log(10),loc=0,scale=10**mean,size=1e4); 
    #numpy_sample = np.random.lognormal(mean=mean*np.log(10),sigma=jsigma*np.log(10),size=1e4)
    # 
    #print 'gauss:',stats.lognorm.fit(gauss_sample,floc=0)
    #print 'scipy:',stats.lognorm.fit(scipy_sample,floc=0)
    #print 'numpy:',stats.lognorm.fit(numpy_sample,floc=0)
    # 
    #plt.figure()
    #plt.hist(gauss_sample,bins=100,alpha=0.5, fc='k', label='gauss',normed=True);
    #plt.hist(scipy_sample,bins=100,alpha=0.5, fc='b', label='scipy',normed=True); 
    #plt.hist(numpy_sample,bins=100,alpha=0.5, fc='r', label='numpy',normed=True);
    #x = np.linspace(0,20,1e4)
    #plt.plot(x,stats.lognorm(jsigma*np.log(10), scale=10**mean, loc=0).pdf(x),c='b',lw=2,label='scipy')
    # 
    #plt.plot(x[1:],lgauss(x[1:], mu=10**mean, sigma=jsigma),c='k',lw=2,label='lgauss')
    #plt.legend()

    plt.show()

