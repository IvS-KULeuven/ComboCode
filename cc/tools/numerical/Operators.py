# -*- coding: utf-8 -*-

"""
Module with several methods that function as operators for numerical grids.

"""

import numpy as np

def diff_central(y,x):

    '''
    Calculate the numerical central difference of a function with respect to a
    single variable. 
    
    This is essentially what np.gradient does, but allows a non-equidistant 
    independent variable.
    
    @param x: The independent variable
    @type x: array
    @param y: The variable dependent on x
    @type y: array
    
    @return: diff_central(y)/diff_central(x)
    @rtype: array

    '''
    
    z1 = np.hstack((y[0], y[:-1]))
    z2 = np.hstack((y[1:], y[-1]))
    dx1 = np.hstack((0, np.diff(x)))
    dx2 = np.hstack((np.diff(x), 0))
    return (z2-z1) / (dx2+dx1)