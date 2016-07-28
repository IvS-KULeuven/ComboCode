# -*- coding: utf-8 -*-

"""
A class for reading and managing collisional rate data.

Author: R. Lombaert

"""

import os, collections
import numpy as np

from cc.tools.readers.SpectroscopyReader import SpectroscopyReader
from cc.tools.io import DataIO

import matplotlib.pyplot as p

from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d



class CollisReader(SpectroscopyReader):
    
    '''
    A Reader for collisional rate data.
    
    This class inherits from SpectroscopyReader. 
    
    Methods for reading and managing collisional rate data are provided. In case
    the basic methods are not sufficient, subclasses provide code-specific 
    methods. The basic methods apply to the collis files for GASTRoNOoM. 
    
    Inheriting from CollisReader are:
        - LamdaReader for reading Lamda molecular data files, used in MCP/ALI.
    
    Other molecular spectroscopic data are typically handled by MolReader. Level
    populations are handled by PopReader. For GASTRoNOoM this is combined in 
    MlineReader for rad trans output, and available through RadiatReader for 
    spectroscopic input.    
    
    '''
    
    def __init__(self,fn,*args,**kwargs):
        
        ''' 
        Initialize an CollisReader object. The filename and additional 
        args/kwargs are passed to the Reader parent class, where the dict is 
        made and the filename is stored. 
        
        Additional args/kwargs are used for the dict creation.
        
        @param fn: The filename of collisional rate data. 
        @type fn: str
        
        '''
        
        super(CollisReader,self).__init__(fn,*args,**kwargs)
        #-- Only call the read method if the class is not inheriting 
        #   CollisReader
        if self.__class__ == CollisReader: 
            self.read()
    
    
    
    def read(self): 
    
        '''
        Read the collision rates file. Assumes GASTRoNOoM format.
        
        To read ALI/MCP collision rates (which are in the Lamda format), make
        use of the LamdaReader, which redefines this method. The other retrieval
        methods remain valid.        
        
        Each transition is stored as an index, and gives upper and lower level
        index as well as the collision rate.
        
        '''
        
        #-- Read the collis file which is just one long column. Assumes there is
        #   at least one zero value for padding! 
        collis = np.loadtxt(self.fn)
        self['pars'] = dict()
        
        #-- The number of transitions is given by the number of non-zero values
        #   in the beginning of the file. 
        ntrans = DataIO.findZero(0,collis)
        n0 = DataIO.findNumber(ntrans,collis)-ntrans
        
        #-- The number of temperatures in the file can be calculated now
        ntemp = (len(collis)-2*(ntrans+n0))/(ntrans+n0+1)
        
        #-- Prep the coll_trans array and add indices for coll_trans
        dtype = [('index',int),('lup',int),('llow',int),('rates',np.ndarray)]
        self['coll_trans'] = np.empty(shape=(ntrans,),dtype=dtype)
        self['coll_trans']['index'] = range(1,ntrans+1)
        
        #-- Add the level indices
        self['coll_trans']['lup'] = collis[:ntrans]
        self['coll_trans']['llow'] = collis[ntrans+n0:2*ntrans+n0]
        
        #-- Retrieve the temperatures. 
        start_i = 2*(ntrans+n0)
        Tgrid = [collis[start_i+i*(ntrans + n0 + 1)] for i in range(ntemp)]
        
        #-- Check if any of them is zero, and 
        #   readjust ntemp (sometimes a 0 temp with 0 rates is present in file)
        Tgrid_real = [Ti for Ti in Tgrid if Ti != 0.]
        ntemp_real = len(Tgrid_real)
        
        #-- Retrieve rates and insert into array. Loop for ntemp, and ignore
        #   rates if T is 0 there.
        rates = np.empty(shape=(ntrans,ntemp_real))
        start_i = start_i + 1
        for i,Ti in enumerate(Tgrid):
            if Ti == 0.0: continue
            this_i = start_i+i*(ntrans + n0 + 1)
            rates[:,i] = collis[this_i:this_i+ntrans]
        
        #-- Save into coll_trans array
        for i in range(ntrans):
            self['coll_trans']['rates'][i] = rates[i,:]
            
        #-- Save additional information
        self['pars']['ncoll_trans'] = ntrans
        self['pars']['ncoll_temp'] = ntemp_real
        self['coll_temp'] = np.array(Tgrid_real)
        
    
    
    def getTI(self,itype='coll_trans',lup=None,llow=None):
        
        '''
        Return the indices of the transitions read from the collisional rate 
        data.
        
        A specific index (or array of indices) can be requested by specifying 
        the lower and upper level index.
        
        LamdaReader first inherits from MolReader, so the default will be itype
        'trans' for LamdaReader.
        
        @keyword itype: The type of index. Default is for CollisReader. Needs to
                        be here in case LamdaReader calls this method, where the
                        type has to be specified. It will call both this method
                        after the version in MolReader and pass on the call.
                        
                        (default: 'coll_trans')
        @type itype: str 
        @keyword lup: The index of the upper energy level in the transition. If
                      both this and llow are None, all transition indices are 
                      returned.
        
                      (default: None)
        @type lup: int
        @keyword llow: The index of the lower energy level in the transition. If
                       both this and lup are None, all transition indices are 
                       returned.
        
                       (default: None)
        @type llow: int
        
        @return: The transition indices
        @rtype: array
        
        '''
        
        return super(CollisReader,self).getTI(itype=itype,lup=lup,llow=llow)
        


    def getTUpper(self,index=None,itype='coll_trans'):
        
        '''
        Return the indices of the upper states of the <nline> included 
        transitions.
        
        These are NOT the quantum numbers! Especially for CO-type molecules, 
        the J-level is equal to level_index-1.
        
        In case a single upper level is requested via index, the level index is
        extracted from the array.
                
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        @keyword itype: The type of index. Default is for CollisReader. 
                        MolReader needs this to be 'trans'.
                        
                        (default: 'coll_trans')
        @type itype: str
                
        @return: The upper level indices/x ordered by transition index.
        @rtype: array
        
        '''
        
        return super(CollisReader,self).getTUpper(itype=itype,index=index)
        
    
    
    def getTLower(self,index=None,itype='coll_trans'):
    
        '''
        Return the indices of the lower states of the <nline> included 
        transitions.
        
        These are NOT the quantum numbers! Especially for CO-type molecules, 
        the J-level is equal to level_index-1.
        
        In case a single lower level is requested via index, the level index is
        extracted from the array.
                
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        @keyword itype: The type of index. Default is for CollisReader. 
                        MolReader needs this to be 'trans'.
                        
                        (default: 'coll_trans')
        @type itype: str
                
        @return: The lower level indices/x ordered by transition index.
        @rtype: array
        
        '''
        
        return super(CollisReader,self).getTLower(itype=itype,index=index)
        
                

    def getRates(self,index=None,lup=None,llow=None):
    
        '''
        Retrieve the collision rates, sorted by collisional transition index.
        
        An index can be specified to retrieve a specific set of collision rates. 
        
        Alternatively, the lup and llow can be specified. index takes 
        precedence
        
        The rates are given as a function of the temperature, accessible by 
        getTemp. 

        In case a single set of rates is requested via index, the rates array is
        extracted from the encompassing array.
                
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        @keyword lup: The index of the upper energy level in the transition. If
                      both this and llow are None, all transition indices are 
                      returned.
        
                      (default: None)
        @type lup: int
        @keyword llow: The index of the lower energy level in the transition. If
                       both this and lup are None, all transition indices are 
                       returned.
        
                       (default: None)
        @type llow: int
                
        @return: The collision rates (in cm^3 s^-1) are returned as arrays 
                 as a function of temperature for given coll_trans indices.
        @rtype: array
        
        '''
        
        #-- In case index is given or max one of llow or lup are given, just 
        #   retrieve
        if not index is None or llow is None or lup is None: 
            return self.get('coll_trans','rates',index)
        
        #-- Else, figure out index first, then get the rates.
        index = self.getTI(llow=llow,lup=lup,itype='coll_trans')
        return self.get('coll_trans','rates',index)
        
    
    
    def getTemp(self): 
    
        '''
        Retrieve the temperature grid for the collision rates. 
        
        @return: The temperature grid in K
        @rtype: array
        
        '''
        
        return self['coll_temp']
    
    
    
    def setInterp(self,index=None,itype='spline',*args,**kwargs):
    
        '''
        Set the interpolator for the collision rates versus temperature.
        
        Additional arguments can be passed to the interpolator object.
        
        @keyword index: The transition index. If default, all collision rates 
                        are interpolated. If an iterable object (such as a list)
                        the method iterates over the indices. If a single value,
                        only one set of rates is interpolated.
                        
                        (default: None)
        @type index: int/list
        @keyword itype: The type of interpolator. Either spline or linear.
        
                       (default: 'spline')
        @type itype: str
    
        '''
        
        #-- Select the interpolation type
        itype = itype.lower()
        if itype == 'linear':
            interp = interp1d
        else:
            interp = spline1d
        
        #-- Set the indices: 
        if index is None:
            index = range(1,self['pars']['ncoll_trans']+1)
        elif not isinstance(indices,collections.Iterable) \
                or isinstance(indices,str):
            index = [index]
                
        self['icoll'] = dict()
        
        #-- Set the interpolators
        for i in index:
            self['icoll'][i] = interp(x=self['coll_temp'],\
                                      y=self.get('coll_trans','rates',i),\
                                      *args,**kwargs)
            
            
    
    def getInterp(self,index):
    
        '''
        Get the interpolator for a given transition index. 
        
        @param index: The level index
        @type index: int
        
        @return: The colision rates (in cm^3 s^-1) interpolator as a function of
                 temperature in K. 
        @rtype: interpolator
                 
        '''
        
        return self['icoll'][index]
        
        
        
    def plotCollis(self,fn=None,indices=None):
    
        '''
        Plot the collision rates for given transitions as a function of 
        temperature.
        
        @keyword fn: Filename of the plot, including path and extension. If None
                     the plot is shown in the python session.

                     (default: None)
        @type fn: str
        @keyword indices: The coll_trans indices to be plotted. Default in case
                          all transitions are to be plotted.
        
                          (default: None)
        @type indices: array
        
        '''
        
        p.clf()
        fig1 = p.figure(1)
        if indices is None: 
            indices = range(1,self['pars']['ncoll_trans']+1)
        for ii in indices:
            p.semilogy(self.getTemp(),self.getRates(index=ii))
        p.xlabel('T (K)', fontsize = 14)
        p.ylabel('Collision Rates (cm$^3$ s$^{-1}$)', fontsize = 14)
        p.xlim(xmin=min(self.getTemp()))
        p.xlim(xmax=max(self.getTemp()))
        if not fn is None:
            p.savefig(fn,bbox_inches='tight')
        else:
            p.show()
    
