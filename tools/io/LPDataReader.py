# -*- coding: utf-8 -*-

"""
A class for reading and managing line profile data files.

Author: R. Lombaert

"""

from cc.tools.io.Reader import Reader



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
    