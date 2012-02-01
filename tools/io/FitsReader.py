# -*- coding: utf-8 -*-

"""
A class for reading and managing FITS files.

Author: R. Lombaert

"""

import os
from scipy import array,arange
import pyfits

from cc.tools.io.LPDataReader import LPDataReader
from cc.tools.io import DataIO



class FitsReader(LPDataReader):
    
    '''
    A FITS file reader for line profiles.
    
    '''
    
    def __init__(self,filename,vlsr=None):
        
        '''
        A FITS file reader for line profiles.
        
        Filename of the FITS file is passed to the object upon creation.
        
        @param filename: The FITS filename, including filepath.
        @type filename: string
        
        @keyword vlsr: The system velocity with respect to local standard of 
                       rest. Only used when this is not found in the fits file.
                       In km/s.
                       
                       (default: None)
        @type vlsr: float
        
        '''
        
        super(FitsReader, self).__init__(filename)
        self.readFits(vlsr)
        
        
    
    def readFits(self,vlsr=None):
        
        ''' 
        Read the FITS file. 
        
        Assumes Tmb flux values in K, with respect to velocity. 
        
        @keyword vlsr: The system velocity with respect to local standard of 
                       rest. Only used when this is not found in the fits file.
                       In km/s.
                       
                       (default: None)
        @type vlsr: float
        
        '''
        
        #- get header, number of grid points and stepsize
        hdr = pyfits.getheader(self.filename)
        n_points = hdr.get('naxis1')
        if n_points is None or n_points == 1: 
            raw_input('Could not find array size of fits file.')
        crpix1 = hdr.get('CRPIX1')
        
        #- load the data
        try:
            lp = pyfits.getdata(self.filename)
        except MemoryError:
            print 'WARNING! Reading %s results in a '%self.filename + \
                  'MemoryError. Ignoring datafile for now.'
            self.contents['flux'] = []
            self.contents['velocity'] = []
            return
        self.contents['flux'] = lp[0][0][0]
        
        #- Check vlsr
        vlsr_keys = ['VELO-LSR','VCORR','VLSR']
        vlsr_fits = None
        while vlsr_fits is None:
            vlsr_fits = hdr.get(vlsr_keys.pop(0))
        
        #- Create velocity OR frequency grid depending on fits file vlsr
        if vlsr_fits is None:
            #- use input vlsr
            self.contents['vlsr'] = vlsr
            crval1 = hdr.get('CRVAL1')
            cdelt1 = hdr.get('CDELT1')
            restfreq = hdr.get('RESTFREQ')
            freq_grid = restfreq + crval1 + (arange(n_points) - crpix1)*cdelt1
            #- c in cm/s, vlsr in km/s, conversion factor of 10**5, final: km/s
            vel_grid = vlsr - (freq_grid-restfreq)/restfreq*self.c/100000.
        else:           
            #- use vlsr from fits file
            self.contents['vlsr'] = vlsr_fits/1000.
            deltav = hdr.get('deltav')
            #- vlsr_fits in m/s, we want km/s, conversion factor of 1000
            vel_grid = (vlsr_fits + (arange(n_points) - crpix1)*deltav)/1000.
        
        #- Make sure the velocity grid is ascending.
        if vel_grid[0] > vel_grid [-1]: 
            vel_grid = vel_grid[::-1]
            self.contents['flux'] = self.contents['flux'][::-1]
        self.contents['velocity'] = vel_grid
        


