# -*- coding: utf-8 -*-

"""
A module for reading and managing text-based line profile data files. 

Author: R. Lombaert

"""

import os
from scipy import array 

import cc.path
from cc.tools.io import DataIO
from cc.tools.readers.LPDataReader import LPDataReader



class TxtReader(LPDataReader):
    
    '''
    A text file reader for line profile data.
    
    '''
    
    def __init__(self,fn,star_name=None,*args,**kwargs):
        
        '''
        A txt file reader for line profile data.
        
        Filename of the txt file is passed to the object upon creation.
        
        Initialisation of the object will also check Vlsr in Star.dat and add 
        it as a keyword in the contents dictionary.
        
        Text files must consist of two columns: velocity (km/s), flux (unit)
        
        Additional args/kwargs are used for the dict creation of the parent of 
        Reader.
        
        @param fn: The txt filename, including filepath.
        @type fn: string
               
        @keyword star_name: The star name if the filename doesn't follow naming
                            conventions. None otherwise.
                            
                            (default: None)
        @type star_name: str
        
        '''
        
        super(TxtReader, self).__init__(fn=fn,star_name=star_name,\
                                        *args,**kwargs)
        self.type = 'txt'
        self.readTxt()
        self.checkVlsr()
        
        
    
    def readTxt(self):
        
        ''' 
        Read the txt file. 
        
        Assumes Tmb flux values in K, with respect to velocity. 
                
        '''
        
        data = DataIO.readCols(filename=self.fn,start_row=0,nans=1)
        if self.fn[-6:] == '.ISPEC':
            del data[0]
            data[0] = data[0]/1000.
        self.contents['velocity'] = data[0]
        self.contents['flux'] = data[1]
        if self.contents['velocity'][0] > self.contents['velocity'][-1]: 
            self.contents['velocity'] = self.contents['velocity'][::-1]
            self.contents['flux'] = self.contents['flux'][::-1]
        self.contents['date_obs'] = 'N.A.'
        self.contents['vlsr'] = None