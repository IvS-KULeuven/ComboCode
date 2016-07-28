# -*- coding: utf-8 -*-

"""
Module for calculating the temperature profile as a function of radius. 

Author: R. Lombaert

"""

import sys
import numpy as np

from cc.modeling.profilers import Profiler



def Teps(r,T0,r0,epsilon):

    '''
    Define a temperature power law.
    
    The functional form is: 
    Tg(r) = T0 * (R0/r)^epsilon
        
    @param r: The radial grid in cm
    @type r: array
    @param T0: The initial temperature at Ri in K
    @type T0: float
    @param r0: The inner radius in cm
    @type r0: float
    @param epsilon: The exponent of the temp law
    @type epsilon: float
    
    @return: The temperature grid (K)
    @rtype: array
    
    '''
    
    return T0 * (r/r0)**-epsilon
    
    

def Tdust(r,T0,r0,s=1):
    
    '''
    Return a dust temperature power law of the form as suggested by 
    observational evidence.
    
    See Thesis p32, where power is -2/(4+s) in 
    T(r) = T_eff*(2*r/R_STAR)**(-2/(4+s))
    and s is typically 1 in the Rayleigh limit and optically thin case, with
    T_eff = T0, R_STAR = r0.
    
    @param r: The radial grid for the power law in Rstar
    @type r: array
    @param T0: The stellar effective temperature
    @type T0: float
    @param r0: The initial value for the radius, typically the stellar radius.
    @type r0
    
    @keyword s: The s parameter in the power law T(r) given above.
    
                (default: 1)
    @type s: int

    @return: The temperature grid (K)
    @rtype: array
    
    '''
    
    return T0*(2.*r/r0)**(-2./(4+s))



class Temperature(Profiler.Profiler): 
    
    '''
    An interface for a temperature profile and its derivative
    
    '''
    
    def __init__(self,r,func,dfunc=None,order=3,*args,**kwargs):
    
        '''
        Create an instance of the Temperature() class. Requires a radial grid  
        and a temperature function.
        
        The function can also be given as an interp1d object.
        
        The optional args and kwargs give the additional arguments for the 
        temperature function, which are ignored in case func is an interp1d 
        object.
        
        @param r: The radial points (cm)
        @type r: array
        @param func: The function that describes the temperature profile. Can be 
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
                
        #-- Do not name func, dfunc, etc in function call, or *args won't work        
        super(Temperature, self).__init__(r,func,dfunc,order,*args,**kwargs)
        
        self.r = self.x
        self.T = self.y
        self.dTdr = self.dydx
        