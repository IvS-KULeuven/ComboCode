# -*- coding: utf-8 -*-

"""
A module for reading and managing fits-based line profile data files. 

Author: R. Lombaert

"""

import os
from scipy import array,arange
#import pyfits
from astropy.io import fits as pyfits

import cc.path
from cc.tools.readers.LPDataReader import LPDataReader
from cc.tools.io import DataIO



def fixFits(fn):

    '''
    Fix fits files in case opening them with FitsReader complains about trailing
    spaces and nonstandard values. 
    
    @param fn: The filename and path to the fitsfile
    @type fn: str
    
    '''
    
    hdulist = pyfits.open(fn,mode='update')
    hdulist.flush()
    hdulist.close()
    


def changeFitsHeader(fn_old,fn_new,key,value,comment=None,**kwargs):

    '''
    Change a key in the header of a fits file and save to a new file.

    Additionali just keywords for the writeto pyfits function can be passed along
    with kwargs. eg output_verify='ignore'.

    @param fn_old: The old filename of the fits file
    @type fn_old: str
    @param fn_new: The new filename of the fits file
    @type fn_new: str
    @param key: The header key to be changed. Has to be present already
                (for now)
    @type key: str
    @param value: The new value
    @type value: str, int, float

    @keyword comment: If desired, the comment for the keyword can be changed
                      here. By default, no change is made.

                      (default: None)
    @type comment: str

    '''
    
    hdulist = pyfits.open(fn_old, mode = 'update')
    hdulist[0].header.update(key,value,comment=comment)
    hdulist.flush()
    hdulist.close()


class FitsReader(LPDataReader):

    '''
    A FITS file reader for line profile data.

    '''

    def __init__(self,fn,star_name=None,*args,**kwargs):

        '''
        A FITS file reader for line profile data.

        Filename of the FITS file is passed to the object upon creation.
        
        Additional args/kwargs are used for the dict creation of the parent of 
        Reader.
        
        @param fn: The FITS filename, including filepath.
        @type fn: string

        @keyword star_name: The star name if the filename doesn't follow naming
                            conventions. None otherwise.

                            (default: None)
        @type star_name: str

        '''

        super(FitsReader, self).__init__(fn=fn,star_name=star_name,\
                                         *args,**kwargs)
        self.type = 'fits' 
        self.readFits()



    def readFits(self):

        '''
        Read the FITS file.

        Assumes Tmb flux values in K, with respect to velocity.

        The source velocity is taken from the fits file if available. Otherwise
        it is taken from Star.dat.

        '''

        #- get header, number of grid points and stepsize
        hdr = pyfits.getheader(self.fn)
        self.contents['hdr'] = hdr
        n_points = hdr.get('naxis1')
        if n_points is None or n_points == 1:
            raw_input('Could not find array size of fits file.')
        crpix1 = hdr.get('CRPIX1')

        #- load the data
        #try:
        lp = pyfits.getdata(self.fn)
        #except MemoryError:
            #print 'WARNING! Reading %s results in a '%self.fn + \
                  #'MemoryError. Ignoring datafile for now.'
            #self.contents['flux'] = []
            #self.contents['velocity'] = []
            #return
        if len(lp.shape) == 4 and lp.shape[0] == lp.shape[1] == lp.shape[2] == 1:
            self.contents['flux'] = lp[0][0][0]
        elif len(lp.shape) == 3 and lp.shape[0] == lp.shape[1] == 1:
            self.contents['flux'] = lp[0][0]
        elif lp.shape[0] != 1:
            self.contents['flux'] = lp
        else:
            raise IOError('Unknown fits file formate for the FitsReader.'+\
                          'Check dimensions fits data. (Then contact Robin)')

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
