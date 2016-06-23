# -*- coding: utf-8 -*-

"""
A tool set for reading and using molecular level populations.

Author: R. Lombaert

"""

import numpy as np
from cc.tools.io import DataIO



class PopReader(Reader):
    
    ''' 
    The PopReader class.
    
    Provides methods to read and manage level populations calculated with 
    radiative transfer such as through GASTRoNOoM and MCP/ALI.
    
    Subclasses provide the read and parse methods that are code-specific.
    
    '''
    
    def __init__(self,fn):
        
        ''' 
        Initialize an PopReader object by setting the contents dict.
        
        @param fn: The filename of the level populations. 
        @type fn: str
        
        '''
        super(Reader,self).__init__()
        self.fn = fn
        
    
    
    def readPop():
    
        '''
        Read the level populations as a function of radius. 
        
        Each level is stored
        as an index according to the column from which it was read.
    
        '''
        
        data = np.loadtxt(self.fn)
        self.contents['r'] = data[:,0]
        self.contents['pop'] = data[:,1:]
        
    
    
    
    def getRadius():
    
        '''
        Return the radial grid.
        
        @return: The radius grid read from the file
        @rtype: array
        
        '''
        
        return self.contents['r']
        
        
        
    def getPop(index=None):
    
        '''
        Return the level populations for a given index or irrespective of index.
        
        @keyword index: The index of the level, if None all levels are returned
                        as an array. Each level is then accessed through 
                        pop[:,index].
        @type index: int
        
        @return: The level populations in the form of a 1d or 2d array, 
                 depending if a single level or multiple levels are requested.
        @rtype: array
        
        '''
        
        if index is None: 
            return self.contents['pop']
        
        return self.contents['pop'][:,index]
            