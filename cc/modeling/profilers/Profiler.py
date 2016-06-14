# -*- coding: utf-8 -*-

"""
Tools for making profiles of any kind.

Author: R. Lombaert

"""

import numpy as np
from scipy.interpolate import interp1d, interp2d, UnivariateSpline, BivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d
from scipy.interpolate import RectBivariateSpline as spline2d
import os

from cc.tools.numerical import Operators as op
from cc.tools.io import DataIO



def waterFraction1StepProfiler(model_id,path_gastronoom,fraction,rfrac):

    '''
    Create a 1-step fractional profile for water.
    
    The original water abundance profile is taken from the output of the 
    original model without fractional abundances. 
    
    These fraction profiles can be used for CHANGE_ABUNDANCE_FRACTION in mline
    
    @param model_id: The model id of the original cooling model
    @type model_id: string
    @param path_gastronoom: The model subfolder in ~/GASTRoNOoM/
    @type path_gastronoom: string
    @param fraction: the fraction used
    @type fraction: float
    @param rfrac: the radius at the step to the fractional abundance [cm]
    @type rfrac: float
    
    '''
    
    rfrac = float(rfrac)
    fraction = float(fraction)
    filename = os.path.join(cc.path.gastronoom,path_gastronoom,'models',\
                            model_id,'coolfgr_all%s.dat'%model_id)
    rad = Gastronoom.getGastronoomOutput(filename=filename,keyword='RADIUS',\
                                         return_array=1)
    fraction_profile = np.ones(len(rad))
    step_index = np.argmin(abs(rad-rfrac))
    fraction_profile[step_index:] = fraction
    output_filename = os.path.join(cc.path.gastronoom,path_gastronoom,\
                                   'profiles',\
                                   'water_fractions_%s_%.2f_r%.3e.dat'\
                                   %(model_id,fraction,rfrac))
    DataIO.writeCols(output_filename,[rad,fraction_profile])



def interp_file(x=None,read_func=DataIO.readCols,xcol=0,ycol=1,itype='spline',\
                ikwargs={'ext': 3,'k': 3},*args,**kwargs):

    '''
    Read an arbitrary file and return an interpolation object for that file. 
    
    Extra arguments are passed to the chosen read_func. This must include a 
    filename. Order of args/kwargs follows that of the requested function.

    @keyword x: The x grid requested for the interpolation. Note that this is a 
                dummy variable to allow Profiler to work with this function. x 
                can however be used for the interpolation if xcol is None, but 
                that requires x to have the same dimensions as the y variable.
                
                (default: None)
    @type x: array
    @keyword read_func: The function used to read the file. Default is 
                        straightforward function for reading columns. 
                        Alternatives include getInputData, getMCMaxOutput, etc.
                        Can be given as a string, in which case the function 
                        must be defined in DataIO.
                        
                        (default: DataIO.readCols)
    @type read_func: function/str
    @keyword xcol: The column index of the X independent variable
    
                   (default: 0)
    @type xcol: int
    @keyword ycol: The column index of the Y dependent variable
    
                   (default: 1)
    @type ycol: int
    @keyword itype: The type of interpolator used. Either 'linear' or 'spline'.
                
                    (default: 'spline')
    @type itype: str
    @keyword ikwargs: Extra keyword arguments for the interpolation. Default 
                      extrapolates by returning boundary values, and 
                      interpolates with spline of order 3.
                    
                      (default: {'ext':3,'k':3})
    @type ikwargs: dict
    
    @return: The interpolation object for the file.
    @rtype: interpolation object
    
    '''
    
    #-- Select the interpolation type
    itype = itype.lower()
    if itype == 'linear':
        interp = interp1d
    else:
        interp = spline1d
    
    #-- Check how the read_func is defined and set it
    if isinstance(read_func,str):
        read_func = getattr(DataIO,read_func)
    
    #-- Read the file and select columns
    data = read_func(*args,**kwargs)
    if not xcol is None: x = data[xcol]
    y = data[ycol]
    
    #-- Create the interpolation object and return
    return (x,y,interp(x=x,y=y,**ikwargs))



def step(x,ylow,yhigh,xstep):

    '''
    Step function. At x <= xstep: ylow, at x > xstep: yhigh.
    
    @param x: x grid (can be array or float)
    @type x: array/float
    @param ylow: y value for x < xstep
    @type ylow: float
    @param yhigh: y value for x > xstep
    @type yhigh: float 
    @param xstep: The boundary value between high and low end
    @type xstep: float
    
    '''
    
    return np.piecewise(x,(x>xstep,x<=xstep),(yhigh,ylow))
    


def constant2D(x,y,cfunc,axis=1,*args,**kwargs):

    '''
    Define a 2D array where the variable is constant with respect to one of the 
    two axes. 
    
    The function can be passed any arbitrary args and kwargs. 
    
    @param x: The primary coordinate on axis 0
    @type x: array
    
    @param y: The secondary coordinate on axis 1. 
    @type y: array
    @param cfunc: The function for the variable axis. Can be an interpolator. 
    @type cfunc: function/interpolator
    
    @keyword axis: The axis that is constant. x is 0, y is 1.
                
                   (default: 1)
    @type axis: int
    
    @keyword args: args to be passed to function. Should be empty in case of
                   interpolator.
                    
                   (default: ())
    @type args: dict
    @keyword kwargs: kwargs to be passed to function. Should be empty in case of
                     interpolator.
                    
                     (default: {})
    @type kwargs: dict
    
    @return: The 2d profile (x.size,y.size)
    @rtype: array    
    
    '''
    
    #-- Default is axis 1, so check if it is 0, otherwise assume 1
    #   In case of 0, evaluate y as variable, otherwise x.
    fvar = cfunc(y if axis == 0 else x,*args,**kwargs)
    
    #-- Set the cst axis to 1. Then multiply the two arrays
    fcst = np.ones_like(x if axis == 0 else y)
    z = np.outer(fcst if axis == 0 else fvar, fvar if axis == 0 else fcst)
    
    return z



def constant(x,c=None,*args,**kwargs): 

    '''
    Define a function for a constant profile. 
    
    @param x: The radial points
    @type x: array
    
    @keyword c: The constant value. If None, it is taken from args/kwargs. This
                allows arbitrary naming of the constant. Except of course either
                x or y. A constant must always be given!
                
                (default: None)
    @type c: float
    
    @keyword args: In case the constant is defined under a different name
                    
                   (default: ())
    @type args: dict
    @keyword kwargs: In case the constant is defined under a different name
                    
                     (default: {})
    @type kwargs: dict
    
    @return: The constant profile, either 1d (x.size) or 2d (x.size,y.size)
    @rtype: array
    
    '''
    
    if c is None: 
        if len(args) + len(kwargs) > 1: 
            raise TypeError("constant() got unexpected arguments")
        elif len(args) + len(kwargs) == 0: 
            raise TypeError("constant() got no arguments")
        c = (args+tuple(kwargs.values()))[0]
    
    return np.zeros_like(x)+c



def zero(x,*args,**kwargs): 

    '''
    Define a function for a zero profile. 
    
    @param x: The coordinate points
    @type x: array
    
    @keyword args: To catch any extra keywords. They are not used.
                    
                   (default: ())
    @type args: dict    
    @keyword kwargs: To catch any extra keywords. They are not used.
                    
                     (default: {})
    @type kwargs: dict
    
    @return: The constant profile, either 1d (x.size) or 2d (x.size,y.size)
    @rtype: array
    
    '''
    
    return np.zeros_like(x)



def zero2D(x,y=None,*args,**kwargs): 

    '''
    Define a function for a zero profile. 
    
    @param x: The coordinate points
    @type x: array
    @keyword y: The secondary points. If None, a 1d array is returned. If not
                None a 2d array is returned.
              
                (default: None)
    @type y: array
    
    @keyword args: To catch any extra keywords. They are not used.
                    
                   (default: ())
    @type args: dict    
    @keyword kwargs: To catch any extra keywords. They are not used.
                    
                     (default: {})
    @type kwargs: dict
    
    @return: The constant profile, either 1d (x.size) or 2d (x.size,y.size)
    @rtype: array
    
    '''
    
    if not y is None: x = np.outer(x,y)
    print x.shape
    
    return np.zeros_like(x)

    
    
class Profiler(object): 
    
    '''
    An interface for creating profiles, and allowing to evaluate them and 
    calculate the central difference with respect to a coordinate grid.
    
    '''
    
    #-- Note that inheriting classes depend on this order of arguments
    def __init__(self,x,func=interp_file,dfunc=None,order=3,*args,\
                 **kwargs):
    
        '''
        Create an instance of the Profiler() class. Requires a coordinate grid 
        and a function object for the profile. A function for the derivative is
        optional. The functions can also be given as an interpolation object.
        
        The optional args and kwargs give the additional arguments for the 
        two function, which are ignored in case func is an interpolation object.
        
        The default coordinate grid is evaluated for both the function and the
        derivative. They are saved in self.y and self.dydx. Alternatively, new
        evaluations can be attained through eval and diff.
        
        @param x: The default coordinate points, minimum three points. In the 
                  case of an interpolation function, this is the default grid
                  returned by the instance. The original x/y of the 
                  interpolated profile are saved as xori/yori in the object.
        @type x: array
        
        @keyword func: The function that describes the profile with respect to 
                       x. Can be given as an interp1d object. Default is a read
                       function that interpolates data and returns the 
                       interpolator object. If interpolation object, x and 
                       eval(x) are assumed to be the original grids
                       
                       (default: interp_file)
        @type func: function/interp1d object
        @keyword dfunc: Function that describes the derivative of the profile 
                        with respect to x. Can be given as an interpolation 
                        object. If None, a generic central difference is taken & 
                        interpolated with a spline of which the order can be 
                        chosen.
        
                        (default: None)
        @type dfunc: function/interpolation object
        @keyword order: Order of the spline interpolation of the derivative. 
                        Default is cubic. Not used for the interpolation if func
                        returns an interpolation object. Use read_order in that
                        case.
                        
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
        
        #-- Check len of coordinate input. Cannot be less than three for the 
        #   derivative.
        if len(x) < 3:
            raise ValueError('Coordinate grid must have more than 2 elements.')
                
        #-- set functional args, remember spline order. args/kwargs for 
        #   func and dfunc are saved separately. They are removed if either are
        #   interpolation objects. They are always accessible, the _** variables
        #   are passed to the evaluation.
        self.args = args
        self.kwargs = kwargs
        self._args = self.args
        self._dargs = self.args
        self._dkwargs = self.kwargs
        self._kwargs = self.kwargs
        self.order = order
        
        #-- By default no interpolation, so leave this off.
        self.interp_func = 0
        self.interp_dfunc = 0
        
        #-- Defaults for xori/yori, not used in case of normal functions
        self.xin = np.empty(0)
        self.yin = np.empty(0)
        
        #-- Evaluate the default grid with function. Set x as None first so it
        #   can actually evaluate. Set x once derivative has been evaluated. 
        self.x = None
        if not (isinstance(func,interp1d) or isinstance(func,UnivariateSpline)):
            #-- Evaluate the function, and check what is returned: array or 
            #   interpolation object
            y = func(x,*args,**kwargs)
            
            #-- Interpolation object: so set that as this instance's func. In 
            #   this case, the variable x passed to the class is the default x
            #   grid of this instance. The original x/y-grid is saved as xi, yi
            if isinstance(y,tuple):
                self.func = y[2]
                self.xin = y[0]
                self.yin = y[1]
                self._args = []
                self._kwargs = {}
                self.interp_func = 1
                self.y = self.func(x)
            else: 
                self.func = func
                self.y = y
                
        #-- func is an interpolation object, so just run the normal evaluation.
        #   Set _args/_kwargs to empty, so none are ever passed to the interpol
        else:
            self.func = func
            self._args = []
            self._kwargs = {}
            self.interp_func = 1
            self.y = self.func(x)
            self.yin = self.y
            self.xin = x
                                                    
        #-- Set the derivative function, resorting to default if needed
        if not dfunc is None:
            self.dfunc = dfunc
        elif self.func == constant:
            self.dfunc = zero
        else: 
            #-- Extend array slightly to allow odeint to succeed.
            #   Need better fix for this.
            #x0 = x[0]-(x[1]-x[0])#*0.5
            #xn = x[-1]+(x[-1]-x[-2])#*0.5
            #x_ext = np.hstack([[x0],x,[xn]])
            
            #-- Evaluate the function, and set up an interpolator for 
            #   central difference. The interpolator will extrapolate 
            #   beyond the given x range. This is necessary for odeint to work.
            #   Usually x-range is not exceeded much. 
            self.dfunc = spline1d(x=x,y=op.diff_central(self.y,x),\
                                  k=self.order)
        
        if (isinstance(self.dfunc,interp1d) \
                or isinstance(self.dfunc,UnivariateSpline)):
            self._dargs = []
            self._dkwargs = {}
            self.interp_dfunc = 1
        
        #-- Evaluate the derivative with the default grid
        self.dydx = self.dfunc(x,*self._dargs,**self._dkwargs)
        
        #-- Now set x.
        self.x = x
    
    
    
    def __call__(self,x=None,warn=1):
    
        '''
        Evaluate the profile function at a coordinate point.
        
        x can be any value or array. If func is an interpolation object, it is
        in principle limited by the x-range of the interpolator. It is advised
        not to extend much beyond the given x-range. No warning is printed here.
        Use eval() if that functionality is preferred.
        
        The default y-grid is returned if x is None.
        
        Passes the call to eval.
                
        @keyword x: The primary coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type x: array/float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The profile evaluated at x
        @rtype: array/float
        
        '''
        
        return self.eval(x,warn)
        
        
    
    def eval(self,x=None,warn=1):
        
        '''
        Evaluate the profile function at a coordinate point.
        
        x can be any value or array. If func is an interpolation object, it is
        in principle limited by the x-range of the interpolator. It is advised 
        not to extend much beyond the given x-range.
        
        The default y-grid is returned if x is None.a
      
        @keyword x: The coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type x: array/float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The profile evaluated at x
        @rtype: array/float
        
        '''
        
        if self.interp_func and warn:
            #-- Are all requested values in range of the original grid?
            if np.any((x>self.xin[-1])|(x<self.xin[0])):
                m = 'Warning! There were values outside of interpolation range.'
                print(m)
        
        if x is None:
            return self.y
        
        return self.func(x,*self._args,**self._kwargs)
            #-- Remove 2nd dimension if it only contains one element.
            #return np.squeeze(profile)
    
    
    
    def diff(self,x=None,warn=1):
        
        '''
        Evaluate the derivative of the profile function at a coordinate point.
        
        x can be any value or array. If func is an interpolation object, it is
        in principle limited by the x-range of the interpolator. It is advised 
        not to extend much beyond the given x-range.
                
        @keyword x: The coordinate point(s). If None, the default 
                    coordinate grid is used.
                    
                    (default: None)
        @type x: array/float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The derivative evaluated at x
        @rtype: array/float
        
        '''
        
        
        if self.interp_dfunc and warn:
            #-- Are all requested values in range of the original grid?
            if np.any((x>self.x[-1])|(x<self.x[0])):
                m = 'Warning! There were values outside of interpolation range.'
                print(m)
        
        if x is None:
            return self.dydx
        
        return self.dfunc(x,*self._dargs,**self._dkwargs)
            
            
            
class Profiler2D(object): 
    
    '''
    An interface for creating 2D profiles, and allowing to evaluate them and 
    calculate the central difference with respect to a 2D coordinate grid.
    
    '''
    
    #-- Note that inheriting classes depend on this order of arguments
    def __init__(self,x,y,func,*args,**kwargs):
    
        '''
        Create an instance of the Profiler2D() class. Requires 2 coordinate 
        arrays and a function object for profile.
       
        The function can also be given as an interpolation object.
        
        The optional args and kwargs give the additional arguments for the 
        function, which are ignored in case func is an interpolation object.
        
        The default coordinate grids are both evaluated for the function. 
        They are saved in self.z, as an array of dimensions (x.size,y.size). 
        Alternatively, new evaluations can be attained through eval and diff.
        
        @param x: The default coordinates of the primary independent variable. 
                  Minimum three points.
        @type x: array
        @param y: The default coordinates of the secondary independent variable. 
                  Minimum three points.
        @type y: array
        @param func: The function that describes the profile with respect to x 
                     and y. Can be given as a 2D interpolation object.
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

        #-- set functional args, remember spline order. args/kwargs for 
        #   func are saved separately. They are removed if either are
        #   interpolation objects. They are always accessible, the _** variables
        #   are passed to the evaluation.
        self.args = args
        self.kwargs = kwargs
        self._args = args
        self._kwargs = kwargs
        self.func = func
        
        if isinstance(self.func,interp2d) \
                or isinstance(self.func,BivariateSpline):
            self.interp_func = 1
            self._args = []
            self._kwargs = {}
        else:
            self.interp_func = 0
        
        #-- Evaluate the default grid with function. Set x as None first so it
        #   can actually evaluate. Set x once derivative has been evaluated.
        self.x = None
        self.y = None
        self.z = self.func(x,y,*self._args,**self._kwargs)
        
        #-- Now set x and y.
        self.x = x
        self.y = y 
        
    
    def __call__(self,x=None,y=None,warn=1):
    
        '''
        Evaluate the profile function at a coordinate point.
        
        x/y can be any value or array. If func is an interpolation object, it is
        in principle limited by the x/y-range of the interpolator. It is advised
        not to extend much beyond the given x/y-range.
        
        If one of the two variables is None, it is replaced by the default grid.
                
        Passes on the call to self.eval().
        
        @keyword x: The primary coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type x: array/float
        @keyword y: The secondary coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type y: array/float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        
        @return: The profile evaluated at x and y
        @rtype: array/float
        
        '''
    
        return self.eval(x=x,y=y,warn=warn)
    
    
    
    def eval(self,x=None,y=None,warn=1):
        
        '''
        Evaluate the profile function at a coordinate point.
        
        x/y can be any value or array. If func is an interpolation object, it is
        in principle limited by the x/y-range of the interpolator. It is advised
        not to extend much beyond the given x/y-range.
        
        If one of the two variables is None, it is replaced by the default grid.
                
        @keyword x: The primary coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type x: array/float
        @keyword y: The secondary coordinate point(s). If None, the default 
                    coordinate grid is used.
        
                    (default: None)
        @type y: array/float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
                
        @return: The profile evaluated at x and y
        @rtype: array/float
        
        '''
        
        if x is None and y is None: 
            return self.z
        
        if x is None:
            x = self.x
            
        if y is None:
            y = self.y
# Ignoring equal arrays for now. Likely faster this way.
#         if np.array_equal(x,self.x) and np.array_equal(y,self.y):
#             return self.z
            
        if self.interp_func and warn:
            #-- Check whether this is either the first time running eval or 
            #   whether all requested values lie in range of the original grid
            if np.any((x>self.x[-1])|(x<self.x[0])) \
                    or np.any((y>self.y[-1])|(y<self.y[0])):
                m = 'Warning! There were values outside of 2D interpolation '+\
                    'range.'
                print(m)
        
        return self.func(x,y,*self._args,**self._kwargs)
    
    
    