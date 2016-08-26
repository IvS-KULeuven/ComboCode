# -*- coding: utf-8 -*-

"""
Module for calculating the temperature profile as a function of radius. 

Author: R. Lombaert

"""

import numpy as np
import collections 
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d

from cc.data import Data
from cc.modeling.profilers import Profiler
from cc.tools.numerical import Operators as op



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
    
    def __init__(self,r,func,dfunc=None,order=3,inner=0,inner_eps=0.5,\
                 *args,**kwargs):
    
        '''
        Create an instance of the Temperature() class. Requires a radial grid  
        and a temperature function.
        
        The function can also be given as an interp1d object.
        
        The optional args and kwargs give the additional arguments for the 
        temperature function, which are ignored in case func is an interp1d 
        object.
        
        An additional option concerns the extrapolation to smaller radii, e.g.
        in the inner wind between the stellar surface and the dust condensation
        radius. In this case, the eval/diff methods will differ between inner 
        wind (r<r0) and the rest of the wind (r>r0) in terms of evaluation. 
        
        At r>r0: The given func (or interpolator) is used.
        At r<r0: A Teps power law is used, for which 1 out of Tstar or 
        epsilon can be defined. The default sets epsilon == 0.5 (the optically
        thin case), but inner_epsilon can be passed to the initialisation if 
        needed. r0, T0 must be defined in kwargs either for the
        main wind's Teps law, or as an additional argument if a file is read.
        
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
        @keyword inner: Applies a power law to the inner region, as a means of
                        extrapolation. Off by default, but can be turned on. 
                        In this case, r0 is taken from kwargs, ie the r0 from 
                        the power law of the wind (if Teps is chosen), or as an 
                        additional keyword if a file is read.
                        
                        (default: 0)
        @type inner: bool
        @keyword inner_eps: The exponent for the power law in the inner wind. 
        
                            (default: 0.5)
        @type inner_eps: float
        
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
        super(Temperature, self).__init__(r,func,dfunc,order,*args,**kwargs)
        
        self.inner = inner
        self.inner_eps = inner_eps

        #-- No extrapolation to the inner wind requested
        if not self.inner: 
            self.r = self.x
            self.T = self.y
            self.dTdr = self.dydx
            return
            
        #-- Alternatively, a power law extrapolation is requested. Extract r0,T0
        self.r0 = kwargs.get('r0')
        self.T0 = kwargs.get('T0')
    
        self.r = None
        self.T = None
        self.y = None
        
        #-- Run the class' own eval() method, with l=self.x, which will not be 
        #   equal to the class' self.r attribute, and thus it will evaluated
        #   properly. Evaluated with inner_eps==0.5. Reevaluated if alpha is
        #   different upon eval call.
        self.T = self.eval(self.x,inner_eps=self.inner_eps)
        self.y = self.T
        
        #-- Now also redo the derivative method. Then call diff to set the 
        #   standard derivative of the instance
        if self.interp_dfunc: 
            self.dfunc = spline1d(x=self.x,y=op.diff_central(self.y,self.x),\
                                  k=self.order)
        
        #-- Evaluate the derivative with the default grid
        self.dydx = self.diff(self.x,inner_eps=0.5)
        
        #-- Now set l, dydl.
        self.r = self.x
        self.dTdr = self.dydx
        


    def __call__(self,r=None,warn=1,*args,**kwargs):
    
        '''
        Evaluate the profile function at a coordinate point.
        
        r can be any value or array. 
        
        The default y-grid is returned if r is None.
        
        Passes the call to eval. Additional keywords for the power law 
        extrapolation can be passed this way as well.
                
        @keyword r: The primary coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type r: array/float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                        
        @return: The profile evaluated at r
        @rtype: array/float
        
        '''
        
        return self.eval(r,warn=warn,*args,**kwargs)
        
        
    
    def eval(self,r=None,warn=1,inner_eps=None):
        
        '''
        Evaluate the profile function at a coordinate point.
        
        r can be any value or array. 
        
        The default y-grid is returned if r is None and inner_eps is 0.5 
        (the default).
      
        @keyword r: The coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type r: array/float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        @keyword inner_eps: The exponent of the Teps power law inner wind
                            extrapolation. If default, the inner_eps defined
                            upon initialisation is used.
        
                            (default: None)
        @type inner_eps: float
        
        @return: The profile evaluated at r
        @rtype: array/float
        
        '''

        #-- First retrieve the original profile. If this is the first call from 
        #   the __init__ method of Opacity(), this will be the standard grid.        
        y = super(Temperature,self).eval(r,warn=warn)
        
        #-- No inner power law requested, just pass on the original 
        if not self.inner: 
            return y
            
        #-- Return self.y since r was given as None, if the eps is correct
        if inner_eps is None: inner_eps = self.inner_eps
        if r is None and inner_eps == self.inner_eps:
            return self.y
        
        #-- r can still be None. Apparently different inner_eps requested.
        #   So calc the profile anew with the inner wind law. Need r defined.
        if r is None: 
            r = self.r
        
        #-- Replace the extrapolated values in the inner wind with the new power
        #   law. Make sure r is an array for this.
        rarr = Data.arrayify(r)
        y[rarr<self.r0] = Teps(rarr[rarr<self.r0],T0=self.T0,r0=self.r0,\
                               epsilon=inner_eps)
        
        return y if isinstance(r,collections.Iterable) else y[0]



    def diff(self,r=None,warn=1,inner_eps=None):
        
        '''
        Evaluate the derivative of the profile function at a coordinate point.
        
        r can be any value or array. 
        
        The default y-grid is returned if r is None and inner_eps is 0.5 
        (the default).
        
        @keyword r: The coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type r: array/float        
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        @keyword inner_eps: The exponent of the Teps power law inner wind
                            extrapolation. If default, the inner_eps defined
                            upon initialisation is used.
        
                            (default: None)
        @type inner_eps: float
                
        @return: The derivative evaluated at r
        @rtype: array/float
        
        '''
        
        #-- First retrieve the original profile. If this is the first call from 
        #   the __init__ method of Opacity(), this will be the standard grid.        
        dydx = super(Temperature,self).diff(r,warn=warn)
        
        #-- No inner power law requested, just pass on the original 
        if not self.inner: 
            return dydx
            
        #-- Return self.dydr since r was given as None, if the eps is correct
        if inner_eps is None: inner_eps = self.inner_eps
        if r is None and inner_eps == self.inner_eps:
            return self.y
        
        #-- r can still be None. Apparently different inner_eps requested.
        #   So calc the profile anew with the inner wind law. Need r defined.
        if r is None: 
            r = self.r
            
        #-- Replace the extrapolated values in the inner wind with the new power
        #   law. Make sure r is an array for this.
        rarr = Data.arrayify(r)
        dy_fac = Teps(rarr[rarr<self.r0],T0=self.T0,r0=self.r0,\
                      epsilon=inner_eps-1)
        dydx[rarr<self.r0] = -dy_fac*inner_eps/self.r0
        
        return dydx if isinstance(r,collections.Iterable) else dydx[0]