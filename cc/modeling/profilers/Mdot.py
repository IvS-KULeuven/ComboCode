# -*- coding: utf-8 -*-

"""
Module for calculating mass loss structures as a function of radius. 

Author: R. Lombaert

"""

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d
from astropy import constants as cst
from astropy import units as u

from cc.modeling.profilers import Profiler



class Mdot(Profiler.Profiler):

    '''
    An interface for the mass-loss rate profile as a function of radius.
    
    '''
    
    def __init__(self,r,func=Profiler.constant,dfunc=None,order=3,\
                 *args,**kwargs):
    
        '''
        Create an instance of the Mdot() class. Requires a radial grid, a 
        function, and its arguments. 
        
        Default is a constant mass-loss rate, in which case mdot must be given
        in Msun/yr.
       
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
        
        super(Mdot, self).__init__(r,func=func,dfunc=dfunc,order=order,\
                                   *args,**kwargs)
        
        self.r = self.x
        self.mdot = self.y
        self.dmdotdr = self.dydx

        