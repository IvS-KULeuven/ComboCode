# -*- coding: utf-8 -*-

"""
A tool set for providing molecular data. Subclasses 
provide code-specific read and partse methods.

Author: R. Lombaert

"""

import numpy as np
from cc.tools.io import DataIO



class MolReader(Reader):
    
    ''' 
    The MolReader class.
    
    Provides methods to manage molecule data (collision rates, excitation 
    levels). 
    
    Subclasses provide the read and parse methods that are code-specific.
    
    '''
    
    def __init__(self):
        
        ''' 
        Initialize an MolReader object by setting the contents dict.
        
        '''
        super(Reader,self).__init__()
        
        
    
    