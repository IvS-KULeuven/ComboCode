# -*- coding: utf-8 -*-

"""
A class for reading and managing FITS files.

Author: R. Lombaert

"""

import os
from scipy import array 
import pyfits

from cc.tools.io.LPDataReader import LPDataReader
from cc.tools.io import DataIO



class FitsReader(LPDataReader):
    
    '''
    A FITS file reader for line profiles.
    
    '''
    
    def __init__(self,filename):
        
        '''
        A FITS file reader for line profiles.
        
        Filename of the FITS file is passed to the object upon creation.
        
        @param filename: The FITS filename, including filepath.
        @type filename: string
        
        '''
        
        super(FitsReader, self).__init__(filename)
        self.readFits()
        
        
    
    def readFits(self):
        
        ''' 
        Read the FITS file. 
        
        Assumes Tmb flux values in K, with respect to velocity. 
        
        '''
        
        hdr = pyfits.getheader(self.filename)
        n_points = hdr.get('naxis1')
        if n_points is None or n_points == 1: raw_input('Could not find array size of fits file.')
        vlsr_keys = ['VELO-LSR','VCORR','VLSR']
        vlsr = None
        while vlsr is None:
            vlsr = hdr.get(vlsr_keys.pop(0))
        self.contents['vlsr'] = vlsr
        try:
            lp = pyfits.getdata(self.filename)
        except MemoryError:
            print 'WARNING! Reading %s results in a '%self.filename + \
                  'MemoryError. Ignoring datafile for now.'
            self.contents['flux'] = []
            self.contents['velocity'] = []
            return
        self.contents['flux'] = lp[0][0][0]
        deltav = hdr.get('deltav')
        vel_grid = [(vlsr + (i-(n_points+1)/2.)*deltav)/1000. 
                    for i in xrange(n_points)]
        if vel_grid[0] > vel_grid [-1]: vel_grid.reverse()
        self.contents['velocity'] = vel_grid
        
