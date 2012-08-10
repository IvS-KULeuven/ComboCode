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
    
    def __init__(self,filename,info_path=os.path.join(os.path.expanduser('~'),\
                                                      'ComboCode','Data')):
        
        '''
        A FITS file reader for line profiles.
        
        Filename of the FITS file is passed to the object upon creation.
        
        @param filename: The FITS filename, including filepath.
        @type filename: string
        
        @keyword info_path: The path to the folder containing the info file on
                            stars, called Star.dat. 
                            
                            (default: ~/ComboCode/Data)
        @type info_path: string

        '''
        
        super(FitsReader, self).__init__(filename,info_path)
        self.readFits()
        
    
    
    def readFits(self):
        
        ''' 
        Read the FITS file. 
        
        Assumes Tmb flux values in K, with respect to velocity. 
        
        The source velocity is taken from the fits file if available. Otherwise
        it is taken from Star.dat.
        
        '''
        
        #- get header, number of grid points and stepsize
        hdr = pyfits.getheader(self.filename)
        self.contents['hdr'] = hdr
        n_points = hdr.get('naxis1')
        if n_points is None or n_points == 1: 
            raw_input('Could not find array size of fits file.')
        crpix1 = hdr.get('CRPIX1')
        
        #- load the data
        #try:
        lp = pyfits.getdata(self.filename)
        #except MemoryError:
            #print 'WARNING! Reading %s results in a '%self.filename + \
                  #'MemoryError. Ignoring datafile for now.'
            #self.contents['flux'] = []
            #self.contents['velocity'] = []
            #return
        if len(lp.shape) == 4 and lp.shape[0] == lp.shape[1] == lp.shape[2] == 1:
            self.contents['flux'] = lp[0][0][0]
        elif lp.shape[0] != 1:
            self.contents['flux'] = lp
        else:
            raise IOError('Unknown fits file formate for the FitsReader. Contact Robin.')

        #- Check vlsr
        vlsr_keys = ['VELO-LSR','VCORR','VLSR']
        vlsr_fits = None
        while vlsr_fits is None and not vlsr_keys == []:
            vlsr_fits = hdr.get(vlsr_keys.pop(0))
        
        #- Create velocity OR frequency grid depending on fits file vlsr
        #- If vlsr_fits if 0.0, it is possible the velocity grid is already
        #- scaled to the right vlsr. In that case, take vlsr from an inputfile.
        if vlsr_fits is None:
            #- use input vlsr
            self.checkVlsr()
            #- if None, vlsr not in fits file, so take guess from Star.dat
            #- if 0.0, it IS in fits file, and needs to be used as such
            crval1 = hdr.get('CRVAL1')
            cdelt1 = hdr.get('CDELT1')
            ctype1 = hdr.get('CTYPE1')
            if ctype1.strip() == 'VRAD':
                #-- If VRAD type on axis 1 use this: Does not need vlsr!
                vel_grid = crval1 + (arange(n_points) - crpix1)*cdelt1
            else:            
                restfreq = hdr.get('RESTFREQ')
                freq_grid = restfreq + crval1 + (arange(n_points) - crpix1)*cdelt1
                #- c in cm/s, vlsr in km/s, conversion factor of 10**5, final: km/s
                vel_grid = self.getVlsr() - (freq_grid-restfreq)/restfreq*self.c/100000.
        else:           
            #- use vlsr from fits file
            #- if it is 0.0, save the initial guess from Star.dat, 
            #if abs(vlsr_fits) < 1e03: self.checkVlsr()
            #- else save value from fits file.
            self.contents['vlsr'] = vlsr_fits/1000.
            self.checkVlsr()
            deltav = hdr.get('deltav')
            #- vlsr_fits in m/s, we want km/s, conversion factor of 1000
            #- use vlsr_fits, regardless of it being 0 or not.
            vel_grid = (vlsr_fits + (arange(n_points) - crpix1)*deltav)/1000.
        
        #- Make sure the velocity grid is ascending.
        if vel_grid[0] > vel_grid [-1]: 
            vel_grid = vel_grid[::-1]
            self.contents['flux'] = self.contents['flux'][::-1]
        self.contents['velocity'] = vel_grid
        
        #- Get the date of observation
        self.contents['date_obs'] = hdr.get('DATE-OBS')
        