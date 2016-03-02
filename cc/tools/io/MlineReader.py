# -*- coding: utf-8 -*-

"""
A class for reading and managing Mline output.

Author: M. van de Sande and R. Lombaert

"""

import os
import numpy as np

from cc.tools.io.Reader import Reader
from cc.tools.io import DataIO



class MlineReader(Reader):
    
    '''
    A Reader for Mline output files. 
    
    NYI - Work in progress.
    
    '''
    
    def __init__(self,filename,ny_up=None,ny_low=None,nline=None,n_impact=None,\
                 n_impact_extra=None):
        
        '''
        Creating an Mline object ready for reading Mline output.
        
        Reading mline output through full filename, including filepath.
        
        The mlinefile number is given as *, which is then replaced as needed.
        
        @param filename: The mline filename, including filepath. The mline
                         file number is given by *.
        @type filename: string
        
        '''
        
        super(SphinxReader, self).__init__()
        self.molecule = os.path.splitext(filename)[0].split('_')[-1]
        self.filename = filename.replace('ml1','ml*').replace('ml2','ml*')\
                                .replace('ml3','ml*')
        


