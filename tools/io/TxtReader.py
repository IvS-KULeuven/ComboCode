# -*- coding: utf-8 -*-

"""
A class for reading and managing txt files. 

Author: R. Lombaert

"""

import os
from scipy import array 
import pyfits

from cc.tools.io import DataIO
from cc.tools.io.LPDataReader import LPDataReader



class TxtReader(LPDataReader):
    
    '''
    A FITS file reader for line profiles.
    
    '''
    
    def __init__(self,filename):
        
        '''
        A txt file reader for line profiles.
        
        Filename of the txt file is passed to the object upon creation.
        
        @param filename: The txt filename, including filepath.
        @type filename: string
        
        '''
        
        super(TxtReader, self).__init__(filename)
        self.readTxt()
        
        
    
    def readTxt(self,vlsr=None):
        
        ''' 
        Read the txt file. 
        
        Assumes Tmb flux values in K, with respect to velocity. 
        
        @keyword vlsr: The system velocity with respect to local standard of 
                       rest. Included in the contents of the data object.
                       In km/s. Used when calculating eg noise values. 
                       
                       (default: None)
        @type vlsr: float
        
        '''
        
        start_i = self.filename[-6:] == '.table' and 1 or 0
        data = DataIO.readCols(filename=self.filename,start_row=start_i,nans=1)
        self.contents['velocity'] = data[0]
        self.contents['flux'] = data[1]
        if self.contents['velocity'][0] > self.contents['velocity'][-1]: 
            self.contents['velocity'] = self.contents['velocity'][::-1]
            self.contents['flux'] = self.contents['flux'][::-1]
        self.contents['date_obs'] = 'N.A.'
        self.contents['vlsr'] = vlsr