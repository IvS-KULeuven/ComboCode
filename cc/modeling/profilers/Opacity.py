# -*- coding: utf-8 -*-

"""
Module for calculating the opacity. 

Author: R. Lombaert

"""

import collections
import numpy as np
from astropy import constants as cst
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d

from cc.modeling.profilers import Profiler
from cc.modeling.profilers import Density    
from cc.tools.readers import KappaReader as KR
from cc.tools.numerical import Operators as op
from cc.data import Data


def read_opacity(x,species,index=1,unit='cm',*args,**kwargs): 

    '''
    Read the opacities with the KappaReader and interpolate. 
    
    Returns a tuple with the original grid, and the interpolation object. Used 
    by the Opacity class to set an interpolation object for a .opac/.particle
    file through species given in usr/Dust.dat.
    
    Additional args and kwargs can be passed to the interpolate method.
    
    @param x: The x grid requested for the interpolation. Note that this is a 
              dummy variable to allow Profiler to work with this function. 
    @type x: array
    @param species: The dust species (from Dust.dat)
    @type species: string
            
    @keyword index: The index of the kappas in the .opacity/.particle file. 
                    0: extinction, 1: absorption, 2: scattering
                    
                    (default: 1)
    @type index: int
    @keyword unit: The unit of the wavelength. Can be given as u.Unit() 
                   object or as a string representation of those objects.
                   Can range from length, to frequency, and energy
    
                   (default: cm)
    @type unit: str/u.Unit()

    @return: The interpolator for the mass extinction/absorption/scattering
             coefficients. (in cgs!)
    @rtype: spline1d 
    
    '''
    
    kr = KR.KappaReader()
    return (kr.getWavelength(species,unit=unit),kr.getKappas(species,index),
            kr.interpolate(species,index,unit=unit,*args,**kwargs))
    


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

        Note that if func is an interpolator object, the original input x and y
        grids can be passed as additional keywords xin/lin and yin, which would 
        then be arrays. Otherwise, the x and the interpolator(x) are set as 
        xin/lin and yin. xin/lin and yin are ignored if func is a function, even
        if it returns an interpolator (in which case the original grids are 
        known)
        
        Extrapolation is done outside of the "original" (depending on if xin/yin
        are given) wavelength ranges according to a power law ~l^-alpha, with 
        alpha given upon eval() calls.
        
        @param l: The wavelength points (cm)
        @type l: array
        @param func: The function that describes the opacity profile. Can be 
                     given as an interpolation object, in which case lin and yin
                     keywords can be passed as arrays for the original grids.
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
        
        #-- Make sure lin is recognized
        if kwargs.has_key('lin'): 
            kwargs['xin'] = kwargs['lin']
            del kwargs['lin']
        
        #-- Do not name func, dfunc, etc in function call, or *args won't work
        super(Opacity, self).__init__(l,func,dfunc,order,*args,**kwargs)
        
        
        #-- Remember the index of the opacities file. If relevant.
        if kwargs.has_key('index'):
            self.index = kwargs['index']
        else:
            self.index = 0
        
        #-- Set the default alpha. Not used if func is not an interpolator.
        self.alpha = 2.
        
        #-- If the requested function is not an interpolator of a file, do not
        #   change anything.
        if not self.interp_func:
            self.l = self.x
            self.kappa = self.y
            return
        
        #-- For Opacity(), in case of an interpolator object as function, 
        #   the extrapolation must be re-done for it to be safe.
        #   The eval() method is rewritten here to function with a
        #   power law extrapolation, for which the settings can be chosen. 
        #   The diff method must not be rewritten since it is based on the y
        #   profile.
        self.l = None
        self.kappa = None
        self.y = None
        
        #-- Run the class' own eval() method, with l=self.x, which will not be 
        #   equal to the class' self.l attribute, and thus it will evaluated
        #   properly. Evaluated with alpha==2. Reevaluated if alpha is different
        self.kappa = self.eval(self.x,alpha=2.)
        self.y = self.kappa
        
        #-- Now also redo the derivative method, in case it is an interpolator
        if self.interp_dfunc: 
            self.dfunc = spline1d(x=self.x,y=op.diff_central(self.y,self.x),\
                                  k=self.order)
        
            #-- Evaluate the derivative with the default grid
            self.dydx = self.dfunc(self.x,*self._dargs,**self._dkwargs)
            
        #-- Now set l, dydl.
        self.l = self.x
        self.dydl = self.dydx
        


    def __call__(self,l=None,*args,**kwargs):
    
        '''
        Evaluate the profile function at a coordinate point.
        
        l can be any value or array. No warning is printed in case of 
        extrapolation, which is done with a power law at small or large l values
        by the eval() function.
        
        The default y-grid is returned if l is None.
        
        Passes the call to eval. Additional keywords for the power law 
        extrapolation can be passed this way as well.
                
        @keyword l: The primary coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type l: array/float
                
        @return: The profile evaluated at l
        @rtype: array/float
        
        '''
        
        return self.eval(l,*args,**kwargs)
        
        
    
    def eval(self,l=None,alpha=2.):
        
        '''
        Evaluate the profile function at a coordinate point.
        
        l can be any value or array. No warning is printed in case of 
        extrapolation, which is done with a power law at small or large l values
        
        The default y-grid is returned if l is None and alpha is 2 (the default)
      
        @keyword l: The coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type l: array/float
        @keyword alpha: The exponent of the wavelength-dependent power law
                        extrapolation, such that kappa ~ lambda^-alpha
        
                        (default: 2.)
        @type alpha: float
                
        @return: The profile evaluated at l
        @rtype: array/float
        
        '''
        
        #-- Return self.y since l was given as None
        if l is None and alpha == self.alpha:
            return self.y
        
        #-- l can still be None. Apparently different alpha requested. So calc.
        if l is None: 
            l = self.l
        
        #-- First retrieve the original profile. Don't warn since we re-do the
        #   extrapolation anyway. If this is the first call from the __init__
        #   method of Opacity(), this will be the standard grid.
        y = super(Opacity,self).eval(l,warn=0)
        
        #-- If self.func is not an interpolator, no warnings needed, and no 
        #   extrapolation needed either.
        if not self.interp_func: return y
        
        #-- Determine the regions where extrapolation is done, i.e. outside the
        #   original l-grid's range. 
        lmin, lmax = self.xin[0], self.xin[-1]
        ymin, ymax = self.yin[0], self.yin[-1]
        
        #-- Replace the extrapolated values with the new power law. Make sure l
        #   and y are an array for this.
        larr, y = Data.arrayify(l), Data.arrayify(y)
        y[larr<lmin] = ymin*(larr[larr<lmin]/lmin)**(-alpha)
        y[larr>lmax] = ymax*(larr[larr>lmax]/lmax)**(-alpha)
        
        return y if isinstance(l,collections.Iterable) else y[0]



    def diff(self,l=None,alpha=2.):
        
        '''
        Evaluate the derivative of the profile function at a coordinate point.
        
        l can be any value or array. No warning is printed in case of 
        extrapolation, which is done with a power law at small or large l values
        
        The default y-grid is returned if l is None and alpha is 2 (the default)
                
        @keyword l: The coordinate point(s). If None, the default 
                    coordinate grid is used.
                    
                    (default: None)
        @type l: array/float
        @keyword alpha: The exponent of the wavelength-dependent power law
                        extrapolation, such that kappa ~ lambda^-alpha
        
                        (default: 2.)
        @type alpha: float
                
        @return: The derivative evaluated at l
        @rtype: array/float
        
        '''
        
        #-- Return self.y since l was given as None
        if l is None and alpha == self.alpha:
            return self.dydx
        
        #-- l can still be None. Apparently different alpha requested. So calc.
        if l is None: 
            l = self.l
        
        #-- First retrieve the original profile. Don't warn since we re-do the
        #   extrapolation anyway. This will have the wrong alpha if alpha is not
        #   2, but we re-do the extrapolation anyway. The original profile will
        #   be untouched.
        dydl = super(Opacity,self).diff(l,warn=0)
        
        #-- If self.func or self.dfunc is not an interpolator, no warnings 
        #   needed, and no extrapolation needed either.
        if not (self.interp_dfunc and self.interp_func): return dydl
        
        #-- Determine the regions where extrapolation is done, i.e. outside the
        #   original l-grid's range. 
        lmin, lmax = self.xin[0], self.xin[-1]
        ymin, ymax = self.yin[0], self.yin[-1]
        
        #-- Replace the extrapolated values with the new power law. Make sure l
        #   and y are an array for this.
        larr, dydl = Data.arrayify(l), Data.arrayify(dydl)
        dydl[larr<lmin] = ymin*(larr[larr<lmin]/lmin)**(-alpha-1.)*-alpha/lmin
        dydl[larr>lmax] = ymax*(larr[larr>lmax]/lmax)**(-alpha-1.)*-alpha/lmax
        
        return dydl if isinstance(l,collections.Iterable) else dydl[0]
        
        
        
    def getEfficiency(self,l,a,sd,P=0.,alpha=2.):
    
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
        @keyword alpha: The exponent of the wavelength-dependent power law
                        extrapolation, such that kappa ~ lambda^-alpha
        
                        (default: 2.)
        @type alpha: float
        
        @return: The dimensionless extinction efficiencies
        @rtype: array
                
        '''
        
        kappas = self.eval(l,alpha=alpha)
        return kappas * 4./3. * sd * a * (1-P)**(2./3.)
    
    
    