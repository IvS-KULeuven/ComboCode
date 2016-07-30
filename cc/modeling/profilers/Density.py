# -*- coding: utf-8 -*-

"""
Module for calculating density structures and molecular abundances as a 
function of radius. 

Author: R. Lombaert

"""

import sys
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d
from astropy import constants as cst
from astropy import units as u

from cc.modeling.profilers import Velocity, Profiler, Mdot


def mean_molecular_weight(fH,fHe):

    '''
    Calculate the mean molecular weight, given the fractional abundances
    fH = n(H)/n(H2)
    and
    fHe = n(He)/n(Htot)
    with
    n(Htot) = n(H) + n(H2)
    
    Derived from: rho_gas / n(gas) = mu * m_H
    
    @param fH: The fractional abundance of atomic hydrogen with respect to
               molecular hydrogen. fH = n(H)/n(H2). Can be a constant or 
               a radially dependent value.
    @type fH: float/array
    @param fHe: The fractional abundance of helium with respect to total 
                hydrogen number density. fHe = n(He)/(n(H)+2*n(H2)). Can be
                a constant or a radially dependent value.
    @type fHe: float/array
    
    @return: The mean molecular weight in units of m_H
    @rtype: float
    
    '''
    
    mu = (fH+2.+4.*fHe*(fH+2.))/(fH+1.+fHe*(fH+2.))
    
    return mu
    


def dens_mdot(r,mdot,v):

    '''
    Calculate the mass density profile based on mass-loss rate (Msun/yr) and 
    expansion velocity (cm/s). 
    
    Based on conservation of mass.
    
    @param r: The radial points (cm)
    @type r: array
    @param mdot: The mass-loss-rate profile or as a single value (constant mass 
                 loss) (Msun/yr)
    @type mdot: Mdot()/float
    @param v: The velocity profile or as a single value (constant velocity) 
              (cm/s)
    @type v: Velocity()/float
    
    @return: The density profile (g/cm3)
    @rtype: array/float
    
    '''
    
    #-- Check if v or mdot are given as floats, in which case create a constant
    #   profile
    vi = v.eval(r) if isinstance(v,Velocity.Velocity) else v
    mdoti = mdot.eval(r) if isinstance(mdot,Mdot.Mdot) else mdot
    
    #-- Calculate the density profile.
    denominator = 4.*np.pi*r**2.*vi
    nominator = (mdoti*u.Msun/u.yr).to(u.g/u.s).value
    return nominator/denominator



class Density(Profiler.Profiler):

    '''
    An interface for the density profile as a function of radius.
    
    '''
    
    def __init__(self,r,func=dens_mdot,dfunc=None,order=3,dust=0,*args,\
                 **kwargs):
    
        '''
        Create an instance of the Density() class. Requires a radial grid, a 
        function for calculating the density profile, and its arguments.
        
        The default is the mass-loss rate density profile. Requires mdot and v
        to be given either as constants or Profile()s. 
       
        @param r: The radial points (cm)
        @type r: array
        
        @keyword func: The function for calculating the density profile. 
        
                       (default: dens_mdot)
        @type func: func
        @keyword dfunc: The function that describes the derivative of the  
                        profile with respect to r. Can be given as an interp1d 
                        object. If None, a generic central difference is taken  
                        and interpolated.
        
                        (default: None)
        @type dfunc: function/interp1d object
        @keyword order: Order of the spline interpolation. Only relevant if 
                        spline is requested instead of linear. Default is cubic
                        interpolation.
                        
                        (default: 3)
        @type order: int
        @keyword dust: This is a dust density profile rather than a gas density
                       
                       (default: 0)
        @type dust: bool
        
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
                
        #-- Do not name func, dfunc, etc in function call, or *args won't work           
        super(Density, self).__init__(r,func,dfunc,order,*args,**kwargs)
        
        self.r = self.x
        self.rho = self.y
        self.drhodr = self.dydx
        
        self.dust = dust
        
        #-- Some constants
        self.mh = cst.m_p.cgs.value      #in g, mass hydrogen atom
        
        #-- Define the number density dictionary and its interpolators
        self.n = dict()
        self.fac = dict()
    
    
    
    def getMeanMolecularWeight(self):
    
        '''
        Return the mean molecular weight for this Density profile. 
        
        Requires calcNumberDensity to have been called before.
        
        @return: The mean molecular weight in units of m_H
        @rtype: float
        
        '''
        
        return self.mu

        
    
    def calcNumberDensity(self,fH=0.,fHe=0.,sd=None,a=None,gsd=None,order=3): 
    
        '''
        In case of gas, calculate the H_2, H, He and total number density 
        profiles. In case of dust, calculate the dust number density. 
        
        Takes into account fractional abundances of atomic hydrogen and helium,
        and the specific density/average grain size, respectively.
        
        Note the specific definitions of fH and fHe. These are implemented this 
        way to work well with the energy balance calculation and the definition
        of dT/dR as well as H_dg.
        
        @keyword fH: The fractional abundance of atomic hydrogen with respect to
                     molecular hydrogen. fH = n(H)/n(H2). Can be a constant or 
                     a radially dependent value.
                     
                     (default: 0)
        @type fH: float/array
        @keyword fHe: The fractional abundance of helium with respect to total 
                      hydrogen number density. fHe = n(He)/(n(H)+2*n(H2)). Can be
                      a constant or a radially dependent value.
                      
                      (default: 0)
        @type fHe: float/array
        @keyword sd: The average specific density of the dust grains in g/cm3.
                     Only relevant if this is a dust density profile. 
                     
                     (default: None)
        @type sd: float     
        @keyword a: The average grain size (cm)
        
                    (default: None)
        @type a: float
           
        '''
        
        if self.dust: 
            self.sd = sd
            self.a = a
            denominator = sd*4./3.*np.pi*a**3
            self.n['Dust'] = self.rho/denominator
            self.fac['Dust'] = 1./denominator
            return
            
        #-- Remember fractional abundances
        self.fH = fH
        self.fHe = fHe
        
        #-- Set the mean molecular weight
        self.mu = mean_molecular_weight(fH,fHe)
        
        #-- Total gas number density is 
        #   n_H + 4*n_He = n(H2) (n(H)/n(H2) + 2)(1 + 4*n(He)/n_H)
        self.n['Gas'] = self.rho/self.mh
        self.fac['Gas'] = 1./self.mh
        
        #-- n(H2) taking into account fractional abundances of H and He
        denominator = self.mh*(fH+2.)*(1.+4.*fHe)
        self.n['H2'] = self.rho/denominator
        self.fac['H2'] = 1./denominator
        
        #-- Atomic hydrogen fractional abundance given with respect to H2
        self.n['H'] = fH*self.n['H2']
        self.fac['H'] = fH/denominator
        
        #-- Total hydrogen abundance n_H = n(H) + 2*n(H2)
        self.n['Htot'] = self.n['H'] + 2.*self.n['H2']
        self.fac['Htot'] = fH/denominator + 2./denominator

        #-- He abundance given by fractional abundance with respect to Htot
        self.n['He'] = fHe*self.n['Htot']
        self.fac['He'] = fHe*(fH/denominator + 2./denominator)
        
        
    
    def eval(self,r=None,dtype='dens',warn=1):
    
        '''
        Evaluate the density function at a coordinate point.
        
        r can be any value or array. If func is an interp1d object, it is
        in principle limited by the r-range of the interpolator. However, 
        linear extrapolation is enabled, but it is advised not to extend much 
        beyond the given r-range.
                
        @keyword r: The coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type r: array/float
        @keyword dtype: The type of density requested. Default is the mass 
                        density. Alternatives: nh2, nh, nhtot, nhe, ngas. 
        
                        (default: dens)
        @type dtype: str
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The profile evaluated at x
        @rtype: array/float
        
        '''
        
        #-- Call the Profiler method in the normal dens case
        dtype = dtype.lower()
        rho = super(Density,self).eval(x=r,warn=warn)
        if not dtype in ['nh2', 'nh', 'nhtot', 'nhe', 'ngas', 'ndust']:
            return rho
        
        #-- Check if number densities have already been calculated
        if not self.n: 
            raise ValueError('Number densities not yet calculated.')
            
        #-- Alternatively, return the default requested number density ...
        dtype = dtype[1:].capitalize()
        if r is None:
            return self.n[dtype]
        
        #-- Otherwise calc the requested number density for the new grid
        return self.fac[dtype]*rho



class DustToGas(Profiler.Profiler):

    '''
    An interface for the dust-to-gas ratio profile as a function of radius.
    
    '''
    
    def __init__(self,r,func=Profiler.constant,dfunc=None,order=3,\
                 *args,**kwargs):
    
        '''
        Create an instance of the DustToGas() class. Requires a radial grid, a 
        function, and its arguments. 
        
        Default is a constant dust-to-gas ratio, in which case d2g must be given
       
        @param r: The radial points (cm)
        @type r: array
        @keyword func: The function for calculating the density profile. 
        
                       (default: dens_mdot)
        @type func: func
        @keyword dfunc: The function that describes the derivative of the  
                        profile with respect to r. Can be given as an interp1d 
                        object. If None, a generic central difference is taken  
                        and interpolated.
        
                        (default: None)
        @type dfunc: function/interp1d object
        @keyword order: Order of the spline interpolation. Default is cubic.
                        
                        (default: 3)
        @type order: int
        
        @keyword args: Additional parameters passed to the functions when eval
                       or diff are called. 
                       
                       (default: [])
        @type args: tuple
        @keyword kwargs: Additional keywords passed to the functions when eval
                         or diff are called. 
                       
                         (default: {})
        @type kwargs: dict
                
        '''

        #-- Do not name func, dfunc, etc in function call, or *args won't work        
        super(DustToGas, self).__init__(r,func,dfunc,order,*args,**kwargs)
        
        self.r = self.x
        self.d2g = self.y
        self.dd2gdr = self.dydx
        