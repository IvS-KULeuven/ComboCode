# -*- coding: utf-8 -*-

"""
A tool set for reading and using molecular level populations.

Author: R. Lombaert and M. van de Sande

"""

import os, collections
import numpy as np
from cc.tools.io import DataIO
from cc.tools.readers.Reader import Reader

import matplotlib.pyplot as p

from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d



class PopReader(Reader):
    
    ''' 
    The PopReader class.
    
    Inherits from the Reader class.
    
    Provides methods to manage level populations calculated with 
    radiative transfer such as through GASTRoNOoM and MCP/ALI. Provides basic 
    methods for reading level populations.
    
    Subclasses provide the read and parse methods that are code-specific, in 
    case the basic methods do not apply:
        - MlineReader which reads mline output and combines spectroscopic 
          molecular data and level populations. 
    
    Typically spectroscopic information is not available here but is given in 
    MolReader objects. Collision rates are given in CollisReader objects. For
    MCP/ALI those are combined in LamdaReader.
    
    '''
    
    def __init__(self,fn,*args,**kwargs):
        
        ''' 
        Initialize an PopReader object. The filename and additional args/kwargs 
        are passed to the Reader parent class, where the dict is made and the 
        filename is stored. 
        
        Additional args/kwargs are used for the dict creation of the parent of 
        Reader.
        
        Note that the filename can be given as None, or may not be an existing
        file. In this case, nothing is read, and the user can set the impact 
        parameter and the level populations manually with setP and setPop.
        
        @param fn: The filename of the level populations. 
        @type fn: str
        
        '''
        
        super(PopReader,self).__init__(fn,*args,**kwargs)
        
        #-- Initialise the pars/pop dict entries
        self['pars'] = {}
        self['pop'] = {}
        
        #-- Only call the read method if the class is not inheriting PopReader
        if self.__class__ == PopReader and os.path.isfile(self.fn):
            self.read()
            
    
    
    def read(self):
    
        '''
        Read the level populations as a function of impact parameter. 
        
        Each level is stored
        as an index according to the column from which it was read.
    
        '''
        
        #-- Basic level populations are read. Applies to MCP/ALI output. Not to
        #   GASTRoNOoM, for which the read function is redefined in MlineReader
        data = np.loadtxt(self.fn)
        
        #-- Save impact parameter grid (already increasing), and create pop dict
        self['p'] = data[:,0]
        
        #-- Loop over the level indices. Note 1-based indexing! Column with 
        #   index 0 is impact parameter.
        self['pars']['ny'] = len(data[0])-1
        for i in range(1,self['pars']['ny']+1):
            self['pop'][i] = data[:,i]
    
    
    
    def setNY(self,ny): 
    
        '''
        Set the number of levels. Usually the highest level index.
        
        @param ny: The number of levels
        @type ny: int
        
        '''
        
        self['pars']['ny'] = ny
        
    
    
    def setP(self,p):
        
        '''
        Instead of reading the impact parameters, set them explicitly. 
        
        @param p: the impact parameters (cm)
        @type p: array
        
        '''
    
        self['p'] = p
        
    
    
    def setPop(self,index,n):
        
        '''
        Instead of reading the populations, set them here per level.
        
        @param index: The level index
        @type index: int
        @param n: The level populations as a function of impact parameter for 
                  level with index
        @type n: array
        
        '''
        
        self['pop'][index] = n
        
        
    
    def getP(self):
    
        '''
        Return the impact parameter grid.
        
        @return: The impact parameter grid in cm
        @rtype: array
        
        '''
        
        return self['p']
        
        
        
    def getPop(self,index=None):
    
        '''
        Return the level populations for a set of level indices.
        
        Note that this is the level index, not lower J. For CO, the J quantum 
        number would be index-1. 
        
        Note that the array approach is not used because that makes indexing 
        (0-based for python, 1-based for fortran) more confusing and prone to
        mistakes. Hence, a dictionary explicity index-key approach is preferred.
        This is internal only. 
        
        If an index is given as an iterable, an array is still returened
        with shape = (len(index),len(self['p'])) for ease of use. This includes
        when all level populations are requested, i.e. when index is None. 
        
        If one wants to use the dictionary itself, self['pop'] is of course 
        available.
        
        @keyword index: The index of the level, if None all levels are returned
                        as a dict. Each level is then accessed through 
                        pop[index].
                        
                        (default: None)
        @type index: int
        
        @return: The level populations in the form of a 1d or 2d array, 
                 depending if a single level or multiple levels are requested.
        @rtype: array
        
        '''
        
        if index is None: 
            return np.array([self['pop'][i] for i in self.getLI()])
        elif isinstance(index,collections.Iterable) \
                and not isinstance(index,str):
            return np.array([self['pop'][i] for i in index])
        else:
            return self['pop'][int(index)]
    
    
    
    def getLI(self):
        
        '''
        Return the indices of the excitation levels included in the level pops.
        
        @return: The level indices
        @rtype: array
        
        '''
        
        return range(1,self['pars']['ny']+1)
        
        
    
    def setInterp(self,itype='spline',*args,**kwargs):
    
        '''
        Set the interpolator for the level populations.
        
        Additional arguments can be passed to the interpolator object.
        
        @keyword itype: The type of interpolator. Either spline or linear.
        
                       (default: 'spline)
        @type itype: str
    
        '''
        
        #-- Select the interpolation type
        itype = itype.lower()
        if itype == 'linear':
            interp = interp1d
        else:
            interp = spline1d
        
        self['ipop'] = dict()
        
        #-- Set the interpolators
        for i in self.getLI():
            self['ipop'][i] = interp(x=self['p'],y=self['pop'][i],\
                                     *args,**kwargs)
            
            
    
    def getInterp(self,index):
    
        '''
        Get the interpolator for a given level index. 
        
        @param index: The level index
        @type index: int
        
        @return: The level population interpolator as a function of impact 
                 parameter in cm. 
        @rtype: interpolator
                 
        '''
        
        return self['ipop'][index]
        
        
    
    def plotPop(self,fn=None):
    
        '''
        Plot the level populations for all included levels.
        
        @keyword fn: Filename of the plot, including path and extension. If None
                     the plot is shown in the python session.

                     (default: None)
        @type fn: str
        
        '''
        
        p.clf()
        fig1 = p.figure(1)
        for ii in range(1,self['pars']['ny']+1):
            p.loglog(self.getP(),self.getPop(ii))
        p.xlabel('P (cm)', fontsize = 14)
        p.ylabel('Level Populations', fontsize = 14)
        p.ylim(ymin=1e-12)
        p.xlim(xmin=self.getP()[0])
        p.xlim(xmax=self.getP()[-1])
        if not fn is None:
            if not os.path.splitext(fn)[1]: fn += '.pdf'
            p.savefig(fn,bbox_inches='tight')
        else:
            p.show()