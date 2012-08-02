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
    
    def __init__(self,filename,info_path=os.path.join(os.path.expanduser('~'),\
                                                      'ComboCode','Data')):
        
        '''
        A txt file reader for line profiles.
        
        Filename of the txt file is passed to the object upon creation.
        
        Initialisation of the object will also check Vlsr in Star.dat and add 
        it as a keyword in the contents dictionary.
        
        @param filename: The txt filename, including filepath.
        @type filename: string
               
        @keyword info_path: The path to the folder containing the info file on
                            stars, called Star.dat. 
                            
                            (default: ~/ComboCode/Data)
        @type info_path: string       
        
        '''
        
        super(TxtReader, self).__init__(filename,info_path)
        self.readTxt()
        self.checkVlsr()
        
        
    
    def readTxt(self):
        
        ''' 
        Read the txt file. 
        
        Assumes Tmb flux values in K, with respect to velocity. 
                
        '''
        
        start_i = self.filename[-6:] == '.table' and 1 or 0
        data = DataIO.readCols(filename=self.filename,start_row=start_i,nans=1)
        self.contents['velocity'] = data[0]
        self.contents['flux'] = data[1]
        if self.contents['velocity'][0] > self.contents['velocity'][-1]: 
            self.contents['velocity'] = self.contents['velocity'][::-1]
            self.contents['flux'] = self.contents['flux'][::-1]
        self.contents['date_obs'] = 'N.A.'
        self.contents['vlsr'] = None