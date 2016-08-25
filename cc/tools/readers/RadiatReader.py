# -*- coding: utf-8 -*-

"""
A GASTRoNOoM radiat file parser. 

Reads the GASTRoNOoM file syntax and provides the information in the form 
consistent with the MolReader, the class from which this class inherits.

MolReader provides the functionality to return all info from Einstein 
coefficients and frequencies, to weights, energy levels and transition levels.

Author: R. Lombaert

"""

import os 
import numpy as np

import cc.path
from cc.tools.io import DataIO
from cc.tools.readers.MolReader import MolReader


class RadiatReader(MolReader):
    
    '''
    Class for working with radiat files from GASTRoNOoM's input. 
    
    Inherits from MolReader, which provides the basic methods for managing 
    spectroscopic data. This class replaces some of those methods to be applied
    to spectroscopic data given in GASTRoNOoM's input format. This is the former
    Radiat class, integrated with the new Reader structure. 
    
    Note that RadiatReader is independent from MlineReader. MlineReader reads 
    the mline spectroscopic output data, which should in principle be identical
    to the input data read with RadiatReader from MOLECULE_radiat.dat and 
    MOLECULE_indices.dat.
    
    Collision rate input data are read with CollisReader.
    
    '''
    
    def __init__(self,fn,ny,nline,*args,**kwargs):
        
        '''
        Initializing an instance of the RadiatReader class.
        
        Additional args and kwargs are passed to the dict creation (parent of 
        Reader)
        
        The number of levels and transitions are saved as instance properties.
        
        @param fn: The radiat filename, including filepath.
        @type fn: string    
        @param ny: The number of levels included in the spectroscopy
        @type ny: int
        @param nline: The number of transitions included in the spectroscopy
        @type nline: int
        
        '''
        
        #-- Create parent instance, and set properties.
        super(RadiatReader,self).__init__(fn,*args,**kwargs)
        self.ny = ny
        self.nline = nline
        
        #-- Add the pars dict to match with the mlinereader.
        self['pars'] = {'nline':nline,'ny':ny}
        
        #-- Read the radiat file. This is a very nonstandard format, so the 
        #   parent read method is not appropriate for this
        self.read()
        
        
        
    def read(self):
         
        '''
        Read the MOLECULE_radiat file and set a dictionary for the instance
        with all the information. 
        
        Done on creation of an instance of the class.
        
        '''
        
        #-- Read the radiat file which is just one long column
        radiat = np.loadtxt(self.fn)
        
        #-- Prep the transition array
        dtype = [('index',int),('lup',int),('llow',int),('frequency',float),\
                 ('einsteinA',float)]
        self['trans'] = np.empty(shape=(self.nline,),dtype=dtype)
        self['trans']['index'] = range(1,self.nline+1)
        
        #-- Prep the level array
        dtype = [('index',int),('weight',float),('energy',float)]
        self['level'] = np.empty(shape=(self.ny,),dtype=dtype)
        self['level']['index'] = range(1,self.ny+1)
        
        #-- The number of values per property is either nline or ny
        nvals = [self.nline,self.nline,self.ny,self.ny,self.nline,self.nline]
        
        #-- Numbers are padded with zeroes at end. The amount of padding is 
        #   unknown. if total length is equal to sum of all step sizes, reading 
        #   radiat is very simple.
        if sum(nvals) == len(radiat):
            n0_nline = 0
            n0_ny = 0
            
        #-- Otherwise a more elaborate method must be used to determine the 
        #   amount of zero-padding that is used for each nline and ny.
        else:
            #-- Einstein A is never 0, so this is easy.
            n0_nline = DataIO.findNumber(self.nline,radiat)-self.nline
            
            #-- Two nline-sized blocks, one ny-sized block. Because the first
            #   energy level might be 0 cm^-1, this is not the real n0_ny
            start_i = 2*(self.nline+n0_nline)+self.ny
            n0_ny1 = DataIO.findNumber(start_i,radiat)-start_i
            
            #-- Second ny-sized block gives the real n0_ny, since llow can never
            #   be zero:
            start_i = DataIO.findZero(start_i+n0_ny1,radiat)
            n0_ny = DataIO.findNumber(start_i,radiat)-start_i
            
        #-- Collect property names, types and number of zero-padding
        props = ['einsteinA','frequency','weight','energy','llow','lup']
        ptypes = ['trans','trans','level','level','trans','trans']
        n0s = [n0_nline,n0_nline,n0_ny,n0_ny,n0_nline,n0_nline]
        
        #-- Save the data into the arrays
        for i,(nval,ptype,prop) in enumerate(zip(nvals,ptypes,props)):
            #-- Determine starting index for this block
            start_i = sum(nvals[:i])+sum(n0s[:i])
            
            #-- Retrieve information
            self[ptype][prop] = radiat[start_i:start_i+nval]
        