# -*- coding: utf-8 -*-

"""
Module for calculating the velocity profile as a function of radius. 

Author: R. Lombaert

"""

import numpy as np
import sys
from astropy import constants as cst
from astropy import units as u
from scipy.integrate import trapz

from cc.modeling.profilers import Profiler
from cc.modeling.profilers import Mdot



def driftRPDF(r,a,l,v,mdot,opac,sd,radiance,T=None,P=0,alpha=0.,mu=2.,\
              w_thermal='none'):

    '''
    Calculate the drift velocity as a function of radius and grain size by 
    equating the radiation pressure force and drag force. 
    
    If the grain size is given as a constant value, the profile is calculated
    only for an average grain size, and is essentially 1d, in which case the
    drift does not depend on grain size.
    
    For now in the approximation that the thermal velocities are small compared
    to the drift. This is only valid in the outer, cool regions of the wind. 
    
    The function integrates the luminosity and the opacity over the wavelength 
    grid to estimate the component for the radiation pressure.
    
    Inclusion of the thermal velocity requires T profile. Mean molecular weight
    assumed to be one in that case. 
    
    @param r: The radial points (cm)
    @type r: array/float
    @param a: The grain size grid (cm)
    @type a: array/float
    @param l: The wavelength grid (cm)
    @type l: array/float
    @param v: The velocity profile (cm/s)
    @type v: float/Velocity()
    @param mdot: The mass-loss profile (Msun/yr)
    @type mdot: float/Mdot()
    @param opac: The opacity profile (cm2/g), must include l-dependence. Can be
                 extinction, absorption, or scattering, but note that the
                 calculation assumes this is an absorption opacity.
    @type opac: Opacity()
    @param sd: The specific density of the dust species
    @type sd: float
    @param radiance: The luminosity profile (ergs/s), must include l-dependence
    @type radiance: Radiance()

    @keyword T: The temperature profile. Only relevant if w_thermal is not 
                none.
                
                (default: None)
    @type T: Temperature()
    @keyword P: The porosity of the grain (or equivalent, to represent 
                effective grain surface). Default is the spherical case. Is
                given as the ratio of vacuum per unit volume of a grain.
                
                (default: 0)
    @type P: float
    @keyword alpha: Sticking coefficient (Kwok 1975, Van Marle 2011). Default is
                    current value in MCP/GASTRoNOoM. Should be 0.25 according to
                    Kwok and van Marle. 1-alpha gives the fraction of elastic
                    collisions.
    
                    (default: 0.)
    @type alpha: float
    @keyword mu: The mean molecular weight. Default is a wind with purely H2. 
                 To get the GASTRoNOoM value of 1.4, set fH to 1.5.
    
                 (default: 2.0)
    @type mu: float
    @keyword w_thermal: Type of thermal velocity. 
                           - sqrt(2kT/m): most probably speed df(v)/dr = 0
                           - sqrt(8/pi kT/m): mean speed (f(v) * v dv)/f(v) dv
                           - sqrt(3kT/m): rms speed sqrt((f(v) * v2 dv)/f(v) dv)
                           - Kwok version: 3/4 sqrt(3kT/m)
                           - Epstein: 4/3 sqrt(8kT/mpi)
                        Choose from: epstein, kwok, prob, mean, rms, none 
                         
                        (default: none)
    @type w_thermal: str
    
    @return: The drift velocity (cm/s) as a function of r and a
    @rtype: array/float
    
    '''
    
    #-- Define some constants 
    c = cst.c.cgs.value
    m_p = cst.m_p.cgs.value
    k_b = cst.k_B.cgs.value
    w_thermal = w_thermal.lower()
    
    #-- Check whether v and mdot are constants
    vi = v.eval(r) if isinstance(v,Velocity) else v
    mdoti = mdot.eval(r) if isinstance(mdot,Mdot.Mdot) else mdot
    mdoti = (mdoti*u.Msun/u.yr).to(u.g/u.s).value
    
    #-- Integrate the opacity and luminosity profiles over wavelength
    #   For this: calculate the emitting surface of the central source
    L = radiance.getLuminosity(l=l,ftype='flambda')
    OpacL = trapz(x=l,y=opac.eval(l)*L)
    
    #-- Calculate the drift for each grain size
    #   1) create the 2d array with r-dependent and a-dependent 1d arrays
    #   2) Then add in anything that's constant
    #   Note that Q(a) = kappa*4/3*a*sd*(1-P)^(2/3)
    arr = np.outer(vi/mdoti,a)
    vK = np.sqrt(arr/(c*(1-alpha))*OpacL*4./3.*sd*(1.-P)**(2./3.))
    
    #-- Calculate the thermal term (see Decin 2006) or return vK (default)
    if not w_thermal in ['kwok','mean','rms','prob','epstein']: 
        return vK

    #-- RMS, ie sqrt(Int(f(v) v^2 dv)/Int(f(v) dv))
    if w_thermal == 'rms':
        vT = (3.*k_b*T.eval(r)/(mu*m_p))**0.5
        
    #-- Most probable speed: d(f(v)/dv = 0)
    elif w_thermal == 'prob':
        vT = (2.*k_b*T.eval(r)/(mu*m_p))**0.5
        
    #-- Mean, ie Int(f(v) v dv)/Int(f(v) dv) ==> likely the most appropriate
    elif w_thermal == 'mean':
        vT = (8.*k_b*T.eval(r)/(mu*m_p*np.pi))**0.5
    
    #-- Implementation by Epstein 1924, with the correct factor 4./3.
    elif w_thermal == 'epstein': 
        vT = 4/3.*(8.*k_b*T.eval(r)/(mu*m_p*np.pi))**0.5
        
    #-- kwok: Unknown factor 3/4 and rms velocity
    else: 
        vT = 0.75*(3.*k_b*T.eval(r)/(mu*m_p))**0.5
    
    #--  Make sure the right axis of vT is multiplied with vK^-1
    vT = np.transpose([vT])
    xT = 0.5*(vT/vK)**2.
    factor = ((1.+xT**2)**0.5-xT)**0.5
    
    return factor*vK



def vbeta(r,r0,v0,vinf,beta):

    '''
    The beta-power law description of the velocity profile. 
    
    @param r: The radial points (cm)
    @type r: array
    @param r0: The inner radius (typically at dust condensation radius, cm)
    @type r0: float    
    @param v0: The initial velocity at r0 (cm/s)
    @type v0: float    
    @param vinf: The terminal wind velocity (cm/s)
    @type vinf: float    
    @param beta: The exponent of the beta law
    @type beta: float    
    
    @return: The beta law velocity profile (cm/s)
    @rtype: array
    
    '''
    
    return v0 + (vinf-v0)*(1-r0/r)**beta



def dvbetadr(r,r0,v0,vinf,beta):

    '''
    The derivative of the vbeta function: analytic. 
    
    @param r: The radial points (cm)
    @type r: array  
    @param r0: The inner radius (typically at dust condensation radius, cm)
    @type r0: float    
    @param v0: The initial velocity at r0 (cm/s)
    @type v0: float    
    @param vinf: The terminal wind velocity (cm/s)
    @type vinf: float    
    @param beta: The exponent of the beta law
    @type beta: float    
    
    @return: The analytic derivative of the beta law velocity profile (cm^2/s)
    @rtype: array
    
    '''
    
    return beta*(vinf-v0)*(1.-r0/r)**(beta-1)*(r0*r**(-2.))



def vdust(r,v,w,*args,**kwargs):
    
    '''
    Calculate the dust velocity profile from a gas velocity and a drift velocity
    profile.
    
    Requires averaging the drift velocity over grain size. Several normalisation
    possibilities available. See Drift().avgDrift.
    
    @param r: The radial points (cm)
    @type r: array
    @param v: The velocity profile (cm/s)
    @type v: float/Velocity()
    @param w: The drift velocity profile (cm/s) as a function of r and a.
    @type w: float/Drift() 
                
    @keyword args: Additional parameters passed to the avgDrift method
                   
                   (default: [])
    @type args: tuple
    @keyword kwargs: Additional keywords passed to the avgDrift method
                   
                     (default: {})
    @type kwargs: dict    
    
    @return: The dust velocity profile (cm/s)
    @rtype: array
    
    '''
    
    vi = v.eval(r) if isinstance(v,Velocity) else v
    wi = w.avgDrift(r,*args,**kwargs) if isinstance(w,Drift) else w
    
    return vi + wi
    

    
class Velocity(Profiler.Profiler): 
    
    '''
    An interface for a velocity profile and its derivative
    
    '''
    
    def __init__(self,r,func,dfunc=None,order=3,*args,**kwargs):
    
        '''
        Create an instance of the Velocity() class. Requires a radial grid and 
        a velocity function.
        
        The function can also be given as an interp1d object.
        
        The optional args and kwargs give the additional arguments for the 
        velocity function, which are ignored in case func is an interp1d object.
        
        @param r: The radial points (cm)
        @type r: array
        @param func: The function that describes the velocity profile. Can be 
                     given as an interp1d object.
        @type func: function
        
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
        
        #-- If the function is given as a string, retrieve it from the local 
        #   module, or from the Profiler module.
        if isinstance(func,str):
            try:
                func = getattr(sys.modules[__name__],func)
            except AttributeError:
                func = getattr(Profiler,func)
        
        #-- In case of vbeta, use the analytic function for the derivative
        if func == vbeta:
            dfunc = dvbetadr        

        #-- Do not name func, dfunc, etc in function call, or *args won't work
        super(Velocity, self).__init__(r,func,dfunc,order,*args,**kwargs)
        
        self.r = self.x
        self.v = self.y
        self.dvdr = self.dydx
        
        
        
class Drift(Profiler.Profiler2D):

    '''
    An interface for the drift velocity as a function of radius and grain size.
    
    '''
    
    def __init__(self,r,a,func,*args,**kwargs):
    
        '''
        Create an instance of the Drift() class. Requires a radial grid, a grain
        size grid and a velocity function.
        
        The function can also be given as an interpolator object.
        
        The optional args and kwargs give the additional arguments for the 
        drift function, which are ignored in case func is an interpolator object
        
        @param r: The radial points (cm)
        @type r: array
        @param a: The grain size grid (cm)
        @type a: array
        @param func: The function that describes the drift profile. Can be 
                     given as an interpolator object.
        @type func: function
        
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
        
        super(Drift, self).__init__(r,a,func,*args,**kwargs)
        self.r = self.x
        self.a = self.y
        self.w = self.z

    
    
    def avgDrift(self,r=None,norm_type='standard',nd=None,warn=1):
        
        '''
        Calculate the average drift velocity as a function of radius, 
        normalising with respect to the average grain size (requires the 
        grain size distribution). 
        
        @keyword r: The radial points (cm). If None, the default r grid is used.
        
                    (default: None)
        @type r: float/array
        @keyword norm_type: The weighting of normalisation for the average drift 
                            velocity: 
                                - collisional: Int(w*nd*a^2*da)/Int(nd*a^2*da)
                                - dens: Int(w*nd*a^3*da)/Int(nd*a^3*da)
                                - nd: Int(w*nd*da)/Int(nd*da)
                                - a: Int(w*a*da)/Int(a*da)
                                - standard: Int(w*da)/Int(da)
                        
                            (default: 'standard')
        @type norm_type: str
        @keyword nd: The grain size distribution. Only needed when norm_type is
                     dens or nd.
                     
                     (default: None)
        @type nd: Distribution()
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                   
        @return: The average drift velocity evaluated at r (cm/s). Array is
                 reduced by one dimension since a drops out as variable.
        @rtype: array/float
        
        ''' 
        
        if self.a.size == 1: 
            return np.squeeze(self.eval(r,warn=warn))
            
        norm_type = norm_type.lower()
        #-- Normalise over grain size
        if norm_type == 'a':
            wsum = trapz(y=self.eval(x=r,warn=warn)*self.a,x=self.a,axis=1)
            norm = trapz(y=self.a*self.a,x=self.a)
                        
        #-- Normalize over grain surface * nd (collisional drift)
        elif norm_type == 'collisional':
            #-- Constant pi drops out due to normalisation
            cs_tot = nd.eval(r,self.a,warn=warn)*self.a**2
            wsum = trapz(y=self.eval(x=r,warn=warn)*cs_tot,x=self.a,axis=1)
            norm = trapz(y=cs_tot,x=self.a,axis=1)
            
        #-- Normalise over number density
        elif norm_type == 'nd':
            ndens = nd.eval(r,self.a,warn=warn)
            wsum = trapz(y=self.eval(x=r,warn=warn)*ndens,x=self.a,axis=1)
            norm = trapz(y=ndens,x=self.a,axis=1)
        
        #-- Normalise over dust density
        elif norm_type == 'dens':
            #-- Constant 4/3*pi*spec_dens drops out due to normalisation.
            rho = nd.eval(r,self.a,warn=warn)*self.a**3
            wsum = trapz(y=self.eval(x=r,warn=warn)*rho,x=self.a,axis=1)
            norm = trapz(y=rho,x=self.a,axis=1)
            
        #-- Default is 'standard', so if neither a or nd or dens, do standard.
        else:
            wsum = trapz(y=self.eval(x=r,warn=warn),x=self.a,axis=1)
            norm = trapz(y=np.ones_like(self.a),x=self.a)

        return wsum/norm

        