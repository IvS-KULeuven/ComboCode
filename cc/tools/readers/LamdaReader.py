# -*- coding: utf-8 -*-

"""
A class for reading Lamda molecular data files. See:
http://home.strw.leidenuniv.nl/~moldata/

Author: R. Lombaert

"""

import os
import numpy as np

from cc.tools.readers.CollisReader import CollisReader
from cc.tools.readers.MolReader import MolReader
from cc.tools.io import DataIO



class LamdaReader(MolReader,CollisReader):
    
    '''
    A Reader for Lamda molecular data files. 
    
    For more information on the Lamda molecular data files, see:
    http://home.strw.leidenuniv.nl/~moldata/
    Lamda files are the input spectroscopic and collisional data files for 
    MCP/ALI.

    Provides methods to manage molecular spectroscopic data (excitation 
    levels, Einstein A coefficients, transitions, etc) through inheritance from 
    the MolReader class, and methods to manage collisional rate data through 
    inheritance from the CollisReader class. 
    
    Note that both parents have methods getTI, getTUpper, getTLower. It 
    functions the same, and the default will be MolReader's method. If you want 
    to retrieve trans indices for the collision rates, specify the itype as 
    'coll_trans'.
    
    '''
    
    def __init__(self,fn,*args,**kwargs):
        
        '''
        Creating a LamdaReader object ready for reading Lamda molecular data.
        
        A full filename is required.
        
        Additional args/kwargs are used for the dict creation of the parent of 
        Reader.
        
        @param fn: The lamda data filename, including filepath.
        @type fn: string        
        
        '''

        super(LamdaReader, self).__init__(fn=fn,*args,**kwargs)
        self.read()
        


    def read(self):
    
        '''
        Read the Lamda file, including molecular spectroscopy and collision 
        rates.
        
        Each level and transition is stored as an index, following the 
        prescription of MolReader and CollisReader.
        
        '''
        
        #-- Read the files lines, and set the single number properties.
        d = DataIO.readFile(self.fn)
        self['pars'] = dict()
        pars = ['molecule','mol_weight','ny']
        dtypes = [str,float,int]
        indices = [1,3,5]
        for i,p,dt in zip(indices,pars,dtypes):
            self['pars'][p] = dt(d[i].split()[0])
        
        #-- Read the levels and their properties. Swap columns 1 and 2 to match
        #   MlineReader's format.
        d1 = np.genfromtxt(self.fn,skip_header=7,max_rows=self['pars']['ny'],\
                           dtype='i8,f8,f8',usecols=(0,2,1),\
                           names=['index','weight','energy'])
        self['level'] = d1
        
        #-- Determine the index from ny and preset lines
        index = 7 + self['pars']['ny'] + 1
        self['pars']['nline'] = int(d[index].split()[0])
        
        #-- Read the transitions and their properties. Swap columns 3 and 4 to
        #   match MlineReader's format.
        d2 = np.genfromtxt(self.fn,usecols=[0,1,2,4,3],dtype='i8,i8,i8,f8,f8',\
                           skip_header=index+2,max_rows=self['pars']['nline'],\
                           names=['index','lup','llow','frequency','einsteinA'])
        self['trans'] = d2
        
        #-- Continue single line parameters.
        index += 1 + self['pars']['nline'] + 2
        self['pars']['npartners'] = int(d[index].split()[0])
        index += 2
        self['pars']['partners'] = d[index:index+self['pars']['npartners']]
        index += self['pars']['npartners'] + 1
        self['pars']['ncoll_trans'] = int(d[index].split()[0])
        index += 2
        self['pars']['ncoll_temp'] = int(d[index].split()[0])
        index += 2
        self['coll_temp'] = np.array(d[index].split(),dtype=float)
        
        #-- Read the collision rates. 
        d3 = np.genfromtxt(self.fn,usecols=[0,1,2],dtype='i8,i8,i8',\
                           skip_header=index+2,names=['index','lup','llow'])
        d4 = np.genfromtxt(self.fn,skip_header=index+2,\
                           usecols=range(3,3+self['pars']['ncoll_temp']))        
        dtype = [('index',int),('lup',int),('llow',int),('rates',np.ndarray)]
        self['coll_trans'] = np.empty(shape=(self['pars']['ncoll_trans'],),\
                                      dtype=dtype)
        self['coll_trans']['index'] = d3['index']
        self['coll_trans']['lup'] = d3['lup']
        self['coll_trans']['llow'] = d3['llow']
        for i in range(self['pars']['ncoll_trans']):
            self['coll_trans']['rates'][i] = d4[i,:]

