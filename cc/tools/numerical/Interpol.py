# -*- coding: utf-8 -*-

"""
Tools for inter- and extrapolation.

Author: R. Lombaert

"""

from scipy import array, hstack
from scipy import exp
from scipy.optimize import leastsq
from scipy import isnan

from cc.plotting import Plotting2


def getResiduals(p,x,y,func='power'):
    
    """
    Calculate residuals between function evaluation and given y-values.
    
    @param p: parameters for function evaluation
    @type p: list
    @param x: the x-grid evaluated by function
    @type x: array
    @param y: The given y-values
    @type y: array
    
    @keyword func: type of function f(x), can be ['power','exp','power10',\
                   'square','linear'] for now.
                      
                   (default: 'power')
    @type func: string
    @return: residuals as a function of the x-grid
    @rtype: array
    
    """    
    
    x,y,p = array(x),array(y),array(p)
    err = y - pEval(x,p,func)                        
    return err
    


def pEval(x,p,func='power',x_unit=None):
    
    """
    Evaluate a function prescription.
    
    @param x: The x-grid for function evaluation
    @type x: array
    @param p: parameters for function evaluation
    @type p: list

    @keyword func: type of function f(x), can be ['power','exp','power10',\
                   'square','linear'] for now.
                       
                   (default: 'power')
    @type func: string
    @keyword x_unit: the unit of the x array, only relevant if func=='planck'. 
                     Units other than micron not yet available.
    
                     (default: None)
    @type x_unit: string
    
    @return: f(x)
    @rtype: array
    
    """
    
    x,p = array(x),array(p)
    if func.lower() == 'power':
        return p[0] + p[2]*(x+p[1])**(p[3])
    if func.lower() == 'exp':
        return p[0] + p[2]*(exp(p[3]*(x+p[1])))
    elif func.lower() == 'power10':
        return p[0] + p[2]*(10**(p[3]*(x+p[1])))
    elif func.lower() == 'square':
        return p[0] + p[2]*(x+p[1]) + p[3]*(x+p[1])**2
    elif func.lower() == 'linear':
        return p[0] + p[2]*(x + p[1])
  
        
        
def fitFunction(x_in,y_in,x_out,show=0,func='linear',initial=[1,0,0.5,-1.5],\
                maxfev=20000):
     
    '''
    Fit a template function to input data.
    
    @param x_in: The input x values
    @type x_in: list/array
    @param y_in: The input y values
    @type y_in: list/array
    @param x_out: The output x grid
    @type x_out: list/array

    @keyword show: show the fit result
    
                   (default: 0)
    @type show: bool
    @keyword func: The function template (see pEval)
    
                   (default: 'linear')
    @type func: string
    @keyword initial: initialisation parameters for your function, number of
                      parameters depends on type of function (see pEval)
    
                      (default: [1,0,0.5,-1.5])
    @type initial: list(float)
    @keyword maxfev: Maximum number of calls to the function for the leastsq 
                     method. Zero means (1+number of elements in x_in)*100
                     
                     (default: 20000)
    @return: The output y values
    @rtype: list/array
    
    '''
    
    p_lsq = leastsq(getResiduals,initial,args=(x_in,y_in,func),\
                    maxfev=20000)[0]
    y_out = pEval(x_out,p_lsq,func=func)
    if show:
        print p_lsq
        Plotting2.plotCols(x=[x_in,hstack([x_out,x_in])],\
                           y=[y_in,hstack([y_out,pEval(x_in,p_lsq,func=func)])],\
                           xlogscale=1,ylogscale=1)
    return y_out
          
          
              
def linInterpol(x_in,y_in,x):
    
    """
    Linear interpolation between two points.
    
    @param x_in: 2 x-points for interpolation
    @type x_in: list
    @param y_in: 2 y-points for interpolation corresponding to x_in
    @type y_in: list
    @param x: x value for interpolation
    @type x: float/array
    
    @return: f(x), where f is linear interpolation from x_in and y_in
    @rtype: float/array
    
    """
    
    x = array(x)
    try:
        return y_in[0] + (y_in[0]-y_in[1])/(x_in[0]-x_in[1]) * (x-x_in[0])
    except ZeroDivisionError:
        print 'Identical x-coordinates were submitted: Division by zero. ' + \
              'Aborting.'
        return
