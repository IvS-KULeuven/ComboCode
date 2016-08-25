# -*- coding: utf-8 -*-

"""
Module for calculating the spectral radiance as a function of wavelength. 

Provides the alternative expressions such as flux density, luminosity, etc as 
well. 

Note that I avoid the term specific or total intensity, as there is confusion 
between the precise definition of the term. Spectral intensity refers to 
radiance integrated over a surface, ie in units of erg/s/Hz/sr. Hence, I 
consistently refer to radiance when specific intensity is concerned.

Also note that radiance is also known as the brightness. That term however 
refers to the subjective brightness of an object, and is not very appropriate in
this context.

For an overview of terminology:
https://en.wikipedia.org/wiki/Radiance

For an overview on relation between brightness (or specific intensity), flux
density, and luminosity:
http://www.cv.nrao.edu/course/astr534/Brightness.html

Author: R. Lombaert

"""

import sys
import numpy as np
from astropy import constants as cst

from cc.modeling.profilers import Profiler



def blackbody(f,T):

    '''
    Calculate a blackbody spectrum (ie brightness) as a function of 
    frequency and for a given temperature.
    
    Gives the energy per unit time per unit area of emitting surface in the 
    normal direction per unit solid angle per unit frequency.
    
    Note that GASTRoNOoM works in brightness! Not flux density!
    
    @param f: The frequency grid (cm)
    @type f: array/float
    @param T: The temperature value (not array, K) 
    @type T: float
    
    @return: The blackbody spectrum in erg/s/cm2/Hz
    @rtype: array/float
    
    '''
    
    #-- Define some constant parameters
    h = cst.h.cgs.value         # in erg*s, Planck constant
    c = cst.c.cgs.value         # in cm/s
    k = cst.k_B.cgs.value       # in erg/K, Boltzmann constant

    return 2*h*f**3/c**2/(np.exp((h*f)/(k*T))-1.)



class Radiance(Profiler.Profiler): 
    
    '''
    An interface for a spectral radiance and its derivative.
    
    '''
    
    def __init__(self,func,f=None,l=None,dfunc=None,order=3,*args,**kwargs):
    
        '''
        Create an instance of the Radiance() class. Requires a frequency grid  
        and a spectral radiance function.
        
        The function can also be given as an interpolation object.
        
        The optional args and kwargs give the additional arguments for the 
        temperature function, which are ignored in case func is an interp1d 
        object.
        
        @param func: The function that describes the intensity profile. Can be 
                     given as an interp1d object.
        @type func: function
        
        @keyword f: The frequency points (Hz). This or wavelength array must
                    be given!  
                    
                    (default: None)
        @type f: array
        @keyword l: The wavelength points (cm). This or frequency array must
                    be given!  
                    
                    (default: None)
        @type l: array        
        @keyword dfunc: The function that describes the derivative of the  
                        profile with respect to r. Can be given as an interp 
                        object. If None, a generic central difference is taken  
                        and interpolated.
        
                        (default: None)
        @type dfunc: function/interpolation object
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
        
        if not (l is None or f is None): 
            raise ValueError('Both wavelength and frequency given. ' + \
                             'Define only one.')
        
        #-- If the function is given as a string, retrieve it from the local 
        #   module, or from the Profiler module.
        if isinstance(func,str):
            try:
                func = getattr(sys.modules[__name__],func)
            except AttributeError:
                func = getattr(Profiler,func)
        
        #-- If wavelength is given, calculate frequency and reverse order
        self.c = cst.c.cgs.value
        if f is None: 
            f = self.c/l[::-1]
        
        #-- Do not name func, dfunc, etc in function call, or *args won't work
        super(Radiance, self).__init__(f,func,dfunc,order,*args,**kwargs)

        #-- Remember frequency. Do not remember wavelength, since the order
        #   would be wrong.
        self.f = self.x
        self.R = self.y
        self.dRdl = self.dydx
        
        #-- Remember Lsun and Rsun for convenience
        self.Lsun = cst.L_sun.cgs.value
        self.Rsun = cst.R_sun.cgs.value
        
    
    
    def setSurface(self,rstar):
    
        '''
        Define the surface area of the star. 
        
        @param rstar: The stellar radius (cm)
        @type rstar: float
        
        '''
        
        self.surface = 4*np.pi*rstar**2
    
    
    
    def eval(self,f=None,l=None,ftype='Rnu',warn=1):
        
        '''
        Evaluate the radiance function at a coordinate point.
        
        l/f can be any value or array. If func is an interpolation object, it is
        in principle limited by the l/f-range of the interpolator. However, 
        extrapolation is enabled, but it is advised not to extend much 
        beyond the given l/f-range.
        
        Can be evaluated both in frequency and wavelength (l or f), and with 
        respect to both frequency and wavelength (ftype). 
        
        @keyword f: The frequency grid (Hz). Passed to self.eval(). Can be None
                    in which case the default grid is used. Can also be given as
                    wavelength 
                    
                    (default: None)
        @type f: array/float
        @keyword l: The wavelength grid (cm). Converted to frequency before 
                    evaluation.
                    
                    (default: None)
        @type l: array/float  
        @keyword ftype: Return the spectral radiance with respect to wavelength
                        or frequency. Given as '*nu' or '*lambda', where * is 
                        R, I, F, ... Default is always '*nu', also
                        when ftype is not recognized. Rlambda is calculated from
                        Rlambda = Rnu*nu/lambda. 
                        
                        (default: 'Rnu')
        @type ftype: str 
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The profile evaluated at l/f (erg/s/cm2/Hz/sr or erg/s/cm3/sr)
        @rtype: array/float
        
        '''
        
        #-- Make sure only one independent variable is given.
        if not (l is None or f is None): 
            msg = 'Both wavelength and frequency given. Define only one.'
            raise ValueError(msg)
        
        #-- Wavelength is given, so calculate frequency from it. Check for l 
        #   instead of f. This way f can still be None, and it will be passed
        #   as such to eval()
        if not l is None: 
            f = self.c/l
        #-- Define wavelength as well in case radiance wrt wavelength is needed
        else:
            l = self.c/self.f if f is None else self.c/f
            
        #-- Note that the frequency grid is not reversed, since the original 
        #   order of the input wavelength should be maintained.
        radiance = super(Radiance,self).eval(x=f,warn=warn)
        
        #-- Check if radiance wrt wavelength is requested
        if ftype.lower()[1:] == 'lambda':
            radiance = radiance*f/l
            
        return radiance
    
    
    
    def diff(self,f=None,l=None,warn=1):
        
        '''
        Evaluate the derivative of the radiance function at a coordinate point.
        
        For now only possible with respect to frequency.
        
        x can be any value or array. If func is an interpolation object, it is
        in principle limited by the x-range of the interpolator. However, 
        linear extrapolation is enabled, but it is advised not to extend much 
        beyond the given x-range.
                
        @keyword f: The frequency grid (Hz). Passed to self.eval(). Can be None
                    in which case the default grid is used. Can also be given as
                    wavelength 
                    
                    (default: None)
        @type f: array/float
        @keyword l: The wavelength grid (cm). Converted to frequency before 
                    evaluation.
                    
                    (default: None)
        @type l: array/float  
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        
        @return: The derivative evaluated at l or f
        @rtype: array/float
        
        '''
        
        #-- Make sure only one independent variable is given.
        if not (l is None or f is None): 
            msg = 'Both wavelength and frequency given. Define only one.'
            raise ValueError(msg)
        
        #-- Wavelength is given, so calculate frequency from it. Check for l 
        #   instead of f. This way f can still be None, and it will be passed
        #   as such to eval()
        if f is None: 
            f = self.c/l
        
        #-- Note that the frequency grid is not reversed, since the original 
        #   order of the input wavelength should be maintained.
        return super(Radiance,self).diff(x=f,warn=warn)
        
        
        
    def getLuminosity(self,surface=None,f=None,l=None,ftype='Lnu',warn=1):
    
        '''
        Return the spectral luminosity, i.e. integrated over the surface
        and solid angle, as a function of frequency or wavelength.
        
        Integrated over frequency, this leads to the (total) luminosity, or 
        also the "flux". Do not integrate over wavelength! Units not right for 
        that.
        
        Note that the integration over solid angle leads to the factor of pi. 
        Fnu = Int(cos(theta) dOmega) = Int Int (cos(theta)sin(theta)dtheta dpsi)
        which has primitive function sin^2(t)/, psi in [0,2pi], theta in 
        [0,pi/2]. Note that this is the outgoing flux, along the direction of 
        the ray (opposite the direction would be [pi/2,pi], since theta is wrt
        the normal of the surface).

        Can be evaluated both in frequency and wavelength (l or f), and with 
        respect to both frequency and wavelength (ftype). 
        
        @keyword surface: The surface over which the energy is being emitted
                          (cm^2). If None, self.setSurface must have been called
                          before. self.surface gives the stellar surface area.
                          
                          (default: None)
        @type surface: float
        @keyword f: The frequency grid (Hz). Passed to self.eval(). Can be None
                    in which case the default grid is used. Can also be given as
                    wavelength 
                    
                    (default: None)
        @type f: array/float
        @keyword l: The wavelength grid (cm). Converted to frequency before 
                    evaluation.
                    
                    (default: None)
        @type l: array/float  
        @keyword ftype: Return spectral luminosity with respect to wavelength
                        or frequency. Given as 'Lnu' or 'Llambda'. Default is 
                        always 'Lnu', also when ftype is not recognized. Llambda
                        is calculated from Llambda = Lnu*nu/lambda.
                        
                        (default: 'Lnu')
        @type ftype: str 
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                        
        @return: The spectral luminosity (erg/s/Hz or erg/s/cm)
        @rtype: array
        
        '''
        
        #-- Retrieve radiance, and convert to luminosity.
        radiance = self.eval(l=l,f=f,ftype=ftype,warn=warn)
        if surface is None: surface = self.surface
        return np.pi*radiance*surface
    
    
    
    def getFluxDensity(self,f=None,l=None,ftype='Fnu',warn=1):
    
        '''
        Return the spectral flux density, i.e. integrated over solid 
        angle, as a function of wavelength.
        
        Integrated over wavelength, this leads to the (total) flux density. 
        For integration over wavelength, convert to Flambda = Fnu*nu/lambda 
        first.
        
        Note that the integration over solid angle leads to the factor of pi. 
        Fnu = Int(cos(theta) dOmega) = Int Int (cos(theta)sin(theta)dtheta dpsi)
        which has primitive function sin^2(t)/2, psi in [0,2pi], theta in 
        [0,pi/2]. Note that this is the outgoing flux, along the direction of 
        the ray (opposite the direction would be [pi/2,pi], since theta is wrt
        the normal of the surface).
        
        Can be evaluated both in frequency and wavelength (l or f), and with 
        respect to both frequency and wavelength (ftype). 
        
        @keyword f: The frequency grid (Hz). Passed to self.eval(). Can be None
                    in which case the default grid is used. Can also be given as
                    wavelength 
                    
                    (default: None)
        @type f: array/float
        @keyword l: The wavelength grid (cm). Converted to frequency before 
                    evaluation.
                    
                    (default: None)
        @type l: array/float  
        @keyword ftype: Return the spectral flux dens with respect to wavelength
                        or frequency. Given as 'Fnu' or 'Flambda'. Default is 
                        always 'Fnu', also when ftype is not recognized. Flambda
                        is calculated from Flambda = Fnu*nu/lambda.
                        
                        (default: 'Fnu')
        @type ftype: str 
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The spectral flux density (erg/s/cm2/Hz or erg/s/cm3)
        @rtype: array
        
        '''
        
        #-- Retrieve radiance, and convert to flux density.
        radiance = self.eval(l=l,f=f,ftype=ftype,warn=warn)
        return np.pi*radiance
    
    
    
    def getIntensity(self,surface=None,f=None,l=None,ftype='Inu',warn=1):
    
        '''
        Return the spectral intensity, i.e. integrated over the surface, as a 
        function of wavelength or frequency.
        
        Integrated over frequency, this leads to the (total) intensity. Do not
        integrate over wavelength! Units not right for that.
        
        Not to be confused with the specific intensity, which is the radiance or
        brightness of the source! 

        Can be evaluated both in frequency and wavelength (l or f), and with 
        respect to both frequency and wavelength (ftype). 
        
        @keyword surface: The surface over which the energy is being emitted
                          (cm^2). If None, self.setSurface must have been called
                          before. self.surface gives the stellar surface area.
                          
                          (default: None)
        @type surface: float
        @keyword f: The frequency grid (Hz). Passed to self.eval(). Can be None
                    in which case the default grid is used. Can also be given as
                    wavelength 
                    
                    (default: None)
        @type f: array/float
        @keyword l: The wavelength grid (cm). Converted to frequency before 
                    evaluation.
                    
                    (default: None)
        @type l: array/float  
        @keyword ftype: Return the spectral intensity with respect to wavelength
                        or frequency. Given as 'Inu' or 'Ilambda'. Default is 
                        always 'Inu', also when ftype is not recognized. Ilambda
                        is calculated from Ilambda = Inu*nu/lambda.
                        
                        (default: 'Inu')
        @type ftype: str 
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The spectral intensity (erg/s/Hz/sr or erg/s/cm/sr)
        @rtype: array
        
        '''
        
        radiance = self.eval(l=l,f=f,ftype=ftype,warn=warn)
        if surface is None: surface = self.surface
        return radiance*surface
        