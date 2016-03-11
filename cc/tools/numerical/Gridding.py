# -*- coding: utf-8 -*-

"""
Tools for making coordinate grids.

Author: R. Lombaert

"""

import numpy as np


def makeGrid(minval,maxval,gridpoints=0,log=0,make_int=0):
    
    """
    Make grid between max and min value.
    
    @param minval: lower boundary of grid range
    @type minval: float
    @param maxval: upper boundary of grid range
    @type maxval: float
    
    @keyword gridpoints: number of grid points, including boundaries. If 0 it 
                         is replaced by 2, if 1 then minval is gridpoint
                                
                         (default: 0)
    @type gridpoints: int
    @keyword log: if grid should be calculated in logspace

                  (default: 0)
    @type log: bool
    @keyword make_int: 0 if final gridpoints should be rounded to nearest 
                       integer
                                
                       (default: 0)
    @type make_int: bool
    
    @return: the grid points including the boundaries
    @rtype: array

    """
    
    #-- Make floats of in and out
    minval = float(minval)
    maxval = float(maxval)
    
    #-- Check validity of number of grid points.
    if int(gridpoints) < 2:
        gridpoints = 2
    else:
        gridpoints = int(gridpoints)

    #-- In case of logspace    
    if log:
        grid = np.logspace(np.log10(minval),np.log10(maxval),gridpoints)
    #-- In case of linear space
    else:
        grid = np.linspace(minval,maxval,gridpoints)
    
    #-- Round to the nearest integer if requested
    if bool(make_int): grid = np.around(grid)
        
    return grid
    