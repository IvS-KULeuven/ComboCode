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
    
    def __init__(self,filename):
        
        '''
        Creating an Mline object ready for reading Mline output.
        
        Reading mline output through full filename, including filepath.
        
        The mlinefile number is given as *, which is then replaced as needed.
        
        @param filename: The mline filename, including filepath. The mline
                         file number is given by *.
        @type filename: string
        
        '''
        
        super(MlineReader, self).__init__()
        self.molecule = os.path.splitext(filename)[0].split('_')[-1]
        self.filename = filename.replace('ml1','ml*').replace('ml2','ml*')\
                                .replace('ml3','ml*')
        self.readMoleculeSettings()
        self.readMlineParameters()
        


    def readMoleculeSettings(self):
    
        '''
        Read molecule spectroscopy settings from the ml1 file.
        
        This includes ny_up, ny_low, nline, n_impact and n_impact_extra.
        
        '''



    def readMlineParameters(self):
        
        '''
        Read the extra mline parameters from the mline log file written by CC.
        
        Only when this file is available!
        
        Nothing is read from the database. This method is standalone, similar to
        SphinxReader. 
        
        '''
        
