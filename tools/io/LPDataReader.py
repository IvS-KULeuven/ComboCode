# -*- coding: utf-8 -*-

"""
A class for reading and managing line profile data files.

Author: R. Lombaert

"""

import os
from scipy import std

from cc.tools.io.Reader import Reader
from cc.data import Data
from cc.tools.io import DataIO


class LPDataReader(Reader):
    
    '''
    A data reader for line profiles.
    
    '''
    
    def __init__(self,filename,info_path=os.path.join(os.path.expanduser('~'),\
                                                      'ComboCode','Data')):
        
        '''
        A data reader for line profiles.
        
        Filename of the data file is passed to the object upon creation.
        
        @param filename: The data filename, including filepath.
        @type filename: string
        
        @keyword info_path: The path to the folder containing the info file on
                            stars, called Star.dat. 
                            
                            (default: ~/ComboCode/Data)
        @type info_path: string
        
        '''
                
        super(LPDataReader, self).__init__()
        self.filename = filename
        self.star_name_gastronoom = os.path.split(self.filename)[1].split('_')[0]
        self.info_path = info_path
        self.c = 2.99792458e10          #in cm/s
        self.contents['vlsr'] = None
    
    
    
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
    
    
    def getVexp(self):
        
        '''
        Return the fitted value of the gas terminal velocity. 
        
        @return: The gas terminal velocity determined by the auto fit of the LP
        @rtype: float
        
        '''
        
        return self.contents['vexp']
        
    
    
    def setVexp(self,vexp):
        
        '''
        Set a value for the gas terminal velocity. 
        
        @param vexp: The gas terminal velocity, determined through whichever 
                     way you want. This value is only used to set the noise 
                     parameter in the object. 
        @type vexp: float
        
        '''
        
        self.contents['vexp'] = float(vexp)
        
    
    def getNoise(self,vexp=None):
        
        '''
        Return the noise value for this dataset.
        
        @keyword vexp: The terminal gas velocity, only required if vexp is not 
                       set yet for this LPDataReader object!
                       
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
        applicable, the initial guess provided by Star.dat is returned.
        
        If the star is not included in Star.dat, the v_lsr is set to 0.0.
        
        @return: The source velocity of the observed object in the dataset
        @rtype: float
        
        '''

        return self.contents['vlsr']
    
    
    
    def checkVlsr(self):
        
        """
        Check if the Vlsr was set correctly. If not, it is taken from Star.dat.
        
        It is assumed the vlsr in the fits file is correct within 25%. If it is
        not, the value from Star.dat is used.
        
        """
        
        try:
            star_index = DataIO.getInputData(path=self.info_path,\
                                            keyword='STAR_NAME_GASTRONOOM')\
                                            .index(self.star_name_gastronoom)
            vlsr = DataIO.getInputData(path=self.info_path,\
                                       keyword='V_LSR')[star_index]
        except KeyError,ValueError: 
            print 'Star not found in Star.dat for %s. '%(self.filename) + \
                  'Add star to Star.dat!'
            raise IOError()
        if self.getVlsr() is None: 
            self.contents['vlsr'] = vlsr
        elif vlsr == 0.0:
            if not abs(self.getVlsr()) < 0.25:
                self.contents['vlsr'] = vlsr
        elif abs(self.getVlsr()/vlsr) > 1.3 or abs(self.getVlsr()/vlsr) < 0.75:
            self.contents['vlsr'] = vlsr
            
        
    
    def getDateObs(self):
    
        """
        If available, get the date of observation.
        
        """
        
        return self.contents['date_obs']
        
    
    
    def setNoise(self,vexp=None):
        
        """
        Calculate the noise of the dataset, by taking the standard deviation
        of the velocity range lower than v_lsr - 2*vexp. 
        
        If this is not available, the velocity range is taken to be where it is
        smaller than 1.1*vexp.
        
        If the size of the above selection still is too small, None is 
        returned.
        
        The gas terminal velocity has to have been set previously or passed to
        this method as a keyword.
        
        @keyword vexp: The terminal gas velocity, only required if vexp is not 
                       set yet for this LPDataReader object!
                       
                       (default: None)
        @type vexp: float
        
        @return: The noise of the dataset
        @rtype: float
        
        """
        
        if vexp <> None and not self.contents.has_key('vexp'):
            self.contents['vexp'] = float(vexp)
        if not self.contents.has_key('vexp'):
            print 'Cannot set noise without an estimate for the terminal gas'+\
                  ' velocity. Pass as keyword to this method!'
            return None
        vel = self.contents['velocity']
        flux = self.contents['flux']
        vlsr = self.contents['vlsr']
        vexp = self.contents['vexp']
        noise = Data.getStd(wave=vel,flux=flux,wmin=vlsr-5*vexp,\
                            wmax=vlsr-2*vexp,minsize=10)
        if noise is None:
            noise = Data.getStd(wave=vel,flux=flux,wmin=vlsr-3*vexp,\
                                wmax=vlsr-1.8*vexp,minsize=10)
        if noise is None:
            noise = Data.getStd(wave=vel,flux=flux,wmin=vlsr-3*vexp,\
                                wmax=vlsr-1.1*vexp,minsize=10)
        if noise is None:
            noise = Data.getStd(wave=vel,flux=flux,wmin=vlsr-3*vexp,\
                                wmax=vlsr-1.1*vexp,minsize=5)
        if noise is None:
            print 'WARNING! Noise cannot be determined, not enough data ' + \
                  'points outside the emission line.'
        self.contents['noise'] = noise
        