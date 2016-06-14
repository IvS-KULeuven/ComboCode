# -*- coding: utf-8 -*-

"""
Module for calculating the opacity. 

Author: R. Lombaert

"""

import sys
import numpy as np
from astropy import constants as cst

from cc.modeling.profilers import Profiler
from cc.modeling.profilers import Density    
from cc.tools.io import KappaReader as KR


def read_opacity(x,species,index=1,order=3): 

    '''
    Read the opacities with the KappaReader and interpolate. 
    
    Returns a tuple with the original grid, and the interpolation object. Used 
    by the Opacity class to set an interpolation object for a .opac/.particle
    file through species given in usr/Dust.dat.
    
    @param x: The x grid requested for the interpolation. Note that this is a 
              dummy variable to allow Profiler to work with this function. 
    @type x: array
    
    @param species: The dust species (from Dust.dat)
    @type species: string
            
    @keyword index: The index of the kappas in the .opacity/.particle file. 
                    0: extinction, 1: absorption, 2: scattering
                    
                    (default: 1)
    @type index: int
    @keyword order: Order of the spline interpolation. Default is cubic.
                    
                    (default: 3)
    @type order: int
            
    @return: The interpolator for the mass extinction/absorption/scattering
             coefficients. (in cgs!)
    @rtype: spline1d 
    
    '''
    
    kr = KR.KappaReader()
    return (kr.getWavelength(species),kr.getKappas(species,index),
            kr.interpolate(species,index,order))
    


class Opacity(Profiler.Profiler): 
    
    '''
    An interface for an opacity profile and any type of extinction quantity.
    
    '''
    
    def __init__(self,l,func,dfunc=None,order=3,*args,**kwargs):
    
        '''
        Create an instance of the Opacity() class. Requires a wavelength grid.
        
        The function can also be given as an interpolation object.
        
        The optional args and kwargs give the additional arguments for the 
        opacity function, which are ignored in case func is an interpolation 
        object.
        
        Eval and diff work with the mass extinction coefficient in cm2/g. 
        
        In case func refers to an interpolation object in KappaReader, the
        args/kwargs should always contain an index indicating whether extinction
        (default), absorption or scattering is requested:
            - 0: extinction, 
            - 1: absorption, 
            - 2: scattering
        
        @param l: The wavelength points (cm)
        @type l: array
        @param func: The function that describes the opacity profile. Can be 
                     given as an interpolation object.
        @type func: function
        
        @keyword dfunc: Function that describes the derivative of the profile 
                        with respect to x. Can be given as an interpolation 
                        object. If None, a generic central difference is taken & 
                        interpolated with a spline of which the order can be 
                        chosen.
        
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

        #-- If the function is given as a string, retrieve it from the local 
        #   module, or from the Profiler module.
        if isinstance(func,str):
            try:
                func = getattr(sys.modules[__name__],func)
            except AttributeError:
                func = getattr(Profiler,func)
                
        #-- Do not name func, dfunc, etc in function call, or *args won't work
        super(Opacity, self).__init__(l,func,dfunc,order,*args,**kwargs)
        
        self.l = self.x
        self.kappa = self.y
        if kwargs.has_key('index'):
            self.index = kwargs['index']
        else:
            self.index = 0
        
        
    
    def getEfficiency(self,l,a,sd,P=0.):
    
        '''
        Calculate the extinction/absorption/scattering (depending on index) 
        efficiency for given grain size and porosity on a wavelength grid.
        
        Follows: Q = (extinction cross section) / (geometric cross section)
        where geometric cross section is the effective surface, hence 
        pi*a^2/(1-P)^(2/3), taking into account the porosity of the grain.
        
        @param l: The wavelength points (cm)
        @type l: array
        @param a: The grain size (cm)
        @type a: float
        @param sd: The specific density of the dust species (g/cm3)
        @type sd: float
        @keyword P: The porosity of the grain (or equivalent, to represent 
                    effective grain surface). Default is the spherical case. Is
                    given as the ratio of vacuum per unit volume of a grain.
                    
                    (default: 0)
        @type P: float
        
        '''
        
        kappas = self.eval(l)
        return kappas * 4./3. * sd * a * (1-P)**(2./3.)
    
    
    