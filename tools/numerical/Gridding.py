# -*- coding: utf-8 -*-

"""
Module for calculating all sorts of grids.

Author: R. Lombaert

"""

from scipy import linspace
from scipy import log10



def makeGrid(minval,maxval,gridpoints=0,log=0,floatornot=0):
    
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
    @keyword floatornot: 0 if final gridpoints should be rounded to nearest 
                         integer
                                
                         (default: 0)
    @type floatornot: bool
    
    @return: the grid points including the boundaries
    @rtype: array

    """
    
    if gridpoints == 0:
        gridpoints = 2
    else:
        gridpoints = int(gridpoints)
    final_grid = []
    
    if log:
        minval = float(log10(minval))
        maxval = float(log10(maxval))
        grid = linspace(minval,maxval,gridpoints)
        if floatornot:
            for i in xrange(len(grid)):
                final_grid.append(int(round(10**(grid[i]))))
        else:
            for i in xrange(len(grid)):
                final_grid.append(10**(grid[i]))
        
    else:
        minval = float(minval)
        maxval = float(maxval)
        grid = linspace(minval,maxval,gridpoints)
        if floatornot:
            for i in xrange(len(grid)):
                final_grid.append(int(round(grid[i])))
        else:
            for i in xrange(len(grid)):
                final_grid.append(grid[i])
    return final_grid
    
