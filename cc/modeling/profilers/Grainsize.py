# -*- coding: utf-8 -*-

"""
Module for calculating the grain size distribution. 

Author: R. Lombaert

"""

import sys
import numpy as np
from astropy import constants as cst
from astropy import units as u
from scipy.integrate import trapz

from cc.modeling.profilers import Profiler
from cc.modeling.profilers import Velocity, Mdot


def MRN(r,a,w,v,mdot_dust,sd,a_min=0.005e-4,a_max=0.25e-4,returnA=0):

    '''
    Calculate the MRN grain size distribution for an AGB wind. 
    
    This assumes an average drift velocity independent of grain size, and a 
    constant dust-to-gas ratio throughout the wind.
    
    The MRN distribution is taken from Mathis, Rumpl and Nordsieck (1977), and 
    reads: 
    n_d(a,r) = A(r) n_Htot(r) a^-3.5
    
    A(r) is a scaling factor that is derived from the dust mass-loss rate. This 
    causes the H-number density to drop out of the equation, so this information
    is not needed:
    A(r) = md / [nHtot 16/3 pi^2 sd r^2 Int((vg + w) a^-0.5 da)]
    
    The MRN distribution assumes minimum grain size of .005 micron and maximum
    grain size of 0.25 micron. While the method allows going outside this range
    (e.g. such as in the case of calculating a local value just outside the 
    range when solving a differential equation), it is recommended changing them
    if the grain size distribution extends beyond it. These are the values used
    to derive A(r).
    
    @param r: The radial grid (cm)
    @type r: array
    @param a: The grain size grid (cm)
    @type a: array
    @param w: The drift velocity profile (cm/s) as a function of r and a.
    @type w: Drift()
    @param v: The gas velocity (cm/s). Either as a profile, or as a cst
    @type v: float/Velocity()
    @param mdot_dust: The dust mass-loss rate (msun/yr). Either as a profile, or
                      as a cst
    @type mdot_dust: float/Mdot()
    @param sd: The average specific density of the dust grains in g/cm3.
    @type sd: float
    
    @keyword a_min: The minimum grain size in cm for MRN
    
                    (default: 0.005e-4)
    @type a_min: float
    @keyword a_max: The maximum grain size in cm for MRN
    
                    (default: 0.25e-4)
    @type a_max: float
    @keyword returnA: Return the scale factor separately as well. This is the 
                      second element of the tuple. Note that this does not 
                      contain Htot, ie A = A_local/nhtot
                      
                      (default: 0)
    @type returnA: bool
    
    @return: The grain size distribution as a function of a
    @rtype: array
    
    '''
    
    #-- Check if the mass-loss rate is given as a constant or a profile
    mddi = mdot_dust.eval(r) if isinstance(mdot_dust,Mdot.Mdot) else mdot_dust
    mddi = (mddi*u.Msun/u.yr).to(u.g/u.s).value
    vi = v.eval(r) if isinstance(v,Velocity.Velocity) else v
    
    #-- Calculate the integration over a for A. Use the drift's default a grid,
    #   because A(r) does not depend on the specific a requested here.
    term1 = 2*(np.sqrt(a_max)-np.sqrt(a_min))*vi
    term2 = trapz(x=w.a,y=w.eval(x=r)*w.a**-0.5)

    #-- Calculate the factors for A
    denominator = 16*np.pi**2*sd*r**2*(term1+term2)
    numerator = mddi*3.
    A = numerator/denominator
    
    #-- Calculate the distribution
    #   Note that we do not need the TOTAL hydrogen density, since it drops out 
    #   through A where nHtot ends up in the denominator
    nd = np.outer(A,a**-3.5)
    
    if returnA: return (nd,A)
    
    return nd
    
    
    
class Distribution(Profiler.Profiler2D): 
    
    '''
    An interface for a grain size distribution as a function of grain size and
    radius, and for given minimum and maximum grain sizes.
    
    '''
    
    def __init__(self,r,a,func,*args,**kwargs):
    
        '''
        Create an instance of the Distribution() class. Requires a radial
        and a grain size grid.
        
        The function can also be given as an interpolation object.
        
        The optional args and kwargs give the additional arguments for the 
        temperature function, which are ignored in case func is an interpolation 
        object.
        
        @param r: The radial points (cm)
        @type r: array
        @param a: The grain size points (cm). Must include minimum and maximum
                  grain sizes, unless manually declared in kwargs.
        @type a: array
        @param func: The function that describes the grain size distribution. 
                     Can be given as an interpolation object.
        @type func: function/interpolation object
                
        @keyword args: Additional parameters passed to the functions when eval
                       or diff are called. 
                       
                       (default: [])
        @type args: tuple
        @keyword kwargs: Additional keywords passed to the functions when eval
                         or diff are called. 
                       
                         (default: {})
        @type kwargs: dict
        
        '''

        #-- If the function is given as a string, retrieve it from the local 
        #   module, or from the Profiler module.
        if isinstance(func,str):
            try:
                func = getattr(sys.modules[__name__],func)
            except AttributeError:
                func = getattr(Profiler,func)
                
        if not kwargs.has_key('a_min'):
            kwargs['a_min'] = a[0]
        if not kwargs.has_key('a_max'):
            kwargs['a_max'] = a[-1]
            
        #-- Do not name func, dfunc, etc in function call, or *args won't work            
        super(Distribution, self).__init__(r,a,func,*args,**kwargs)
        
        self.r = self.x
        self.a = self.y
        self.nd = self.z
        