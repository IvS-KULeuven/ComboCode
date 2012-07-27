# -*- coding: utf-8 -*-

"""
A class for reading and managing line profile data files.

Author: R. Lombaert

"""

from scipy import std

from cc.tools.io.Reader import Reader
from cc.data import Data


class LPDataReader(Reader):
    
    '''
    A data reader for line profiles.
    
    '''
    
    def __init__(self,filename):
        
        '''
        A data reader for line profiles.
        
        Filename of the data file is passed to the object upon creation.
        
        @param filename: The data filename, including filepath.
        @type filename: string
        
        '''
        
        super(LPDataReader, self).__init__()
        self.filename = filename
        self.c = 2.99792458e10          #in cm/s
    
    
    
    def getVelocity(self):
        
        '''
        Return the velocity grid of this data file.
        
        @return: The velocity grid
        @rtype: array/list
        
        '''
        
        return self.contents['velocity']
    
    
    
    def getFlux(self):
        
        '''
        Return the flux of the line profile. Usually Tmb.
        
        @return: The flux 
        @rtype: array/list
        
        '''
        
        return self.contents['flux']
    
    
    def getNoise(self,vexp=None):
        
        '''
        Return the noise value for this dataset.
        
        @keyword vexp: The terminal gas expansion velocity. Used to determine the
                     velocity range for noise calculation. Only relevant if 
                     noise was not yet set.
                     
                     (default: None)
        @type vexp: float
        
        @return: The std of the dataset
        @rtype: float
        
        '''
        
        if not self.contents.has_key('noise'):
            self.setNoise(vexp)
        return self.contents['noise']
        
        
   
    def getVlsr(self):
        
        '''
        Return the source velocity of the observed object in this dataset.
        
        This is the value read from the fits file, and if not available or 
        applicable, the initial guess provided by the user is returned.
        
        @return: The source velocity of the observed object in the dataset
        @rtype: float
        
        '''

        return self.contents['vlsr']
    
    
    
    def getDateObs(self):
    
        """
        If available, get the date of observation.
        
        """
        
        return self.contents['date_obs']
        
    
    
    def setNoise(self,vexp):
        
        """
        Calculate the noise of the dataset, by taking the standard deviation
        of the velocity range lower than v_lsr - 2*vexp. 
        
        If this is not available, the velocity range is taken to be where it is
        smaller than 1.2*vexp.
        
        If the size of the above selection still is too small, None is 
        returned.
        
        @param vexp: The terminal gas expansion velocity. Used to determine the
                     velocity range for noise calculation.
        @type vexp: float
        
        @return: The noise of the dataset
        @rtype: float
        
        """
        
        vel = self.contents['velocity']
        flux = self.contents['flux']
        vlsr = self.contents['vlsr']
        vexp = abs(float(vexp))
        noise = Data.getStd(vel,flux,wmin=vlsr-5*vexp,wmax=vlsr-2*vexp)
        if noise is None:
            noise = Data.getStd(vel,flux,wmin=vlsr-5*vexp,wmax=vlsr-1.2*vexp)
        if noise is None:
            print 'WARNING! Noise cannot be determined, not enough data ' + \
                  'points outside the emission line.'
        self.contents['noise'] = noise
        