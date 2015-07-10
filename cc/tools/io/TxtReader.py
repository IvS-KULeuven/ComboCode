# -*- coding: utf-8 -*-

"""
A class for reading and managing txt files. 

Author: R. Lombaert

"""

import os
from scipy import array 
import pyfits

import cc.path
from cc.tools.io import DataIO
from cc.tools.io.LPDataReader import LPDataReader



class TxtReader(LPDataReader):
    
    '''
    A FITS file reader for line profiles.
    
    '''
    
    def __init__(self,filename,star_name=None):
        
        '''
        A txt file reader for line profiles.
        
        Filename of the txt file is passed to the object upon creation.
        
        Initialisation of the object will also check Vlsr in Star.dat and add 
        it as a keyword in the contents dictionary.
        
        @param filename: The txt filename, including filepath.
        @type filename: string
               
        @keyword star_name: The star name if the filename doesn't follow naming
                            conventions. None otherwise.
                            
                            (default: None)
        @type star_name: str
        
        '''
        
        super(TxtReader, self).__init__(filename=filename,star_name=star_name)
        self.readTxt()
        self.checkVlsr()
        
        
    
    def readTxt(self):
        
        ''' 
        Read the txt file. 
        
        Assumes Tmb flux values in K, with respect to velocity. 
                
        '''
        
        data = DataIO.readCols(filename=self.filename,start_row=0,nans=1)
        if self.filename[-6:] == '.ISPEC':
            del data[0]
            data[0] = data[0]/1000.
        self.contents['velocity'] = data[0]
        self.contents['flux'] = data[1]
        if self.contents['velocity'][0] > self.contents['velocity'][-1]: 
            self.contents['velocity'] = self.contents['velocity'][::-1]
            self.contents['flux'] = self.contents['flux'][::-1]
        self.contents['date_obs'] = 'N.A.'
        self.contents['vlsr'] = None