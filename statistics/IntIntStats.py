# -*- coding: utf-8 -*-

"""
Performing statistics on the integrated intensity of resolved molecular lines.

Author: R. Lombaert

"""

import os
import scipy
from scipy import argmin,array
import operator

from cc.tools.io import DataIO
from cc.statistics.Statistics import Statistics



class IntIntStats(Statistics):
    
    """
    Environment with several tools to perform statistics on integrated 
    intensities of resolved molecular lines.
    
    """
        
    def __init__(self,star_name,path_code='codeSep2010',\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):        
        
        """ 
        Initializing an instance of IntIntStats.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2010')
        @type path_code: string
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        
        """
        
        super(IntIntStats,self).__init__(star_name=star_name,\
                                         path_combocode=path_combocode,\
                                         code='GASTRoNOoM',path_code=path_code)
        #- Dict keeping integrated intensities of datafiles. 
        #- key: Transition(), value: list of values for every datafile included
        #- ~ self.transitions[INSTRUMENT] where INSTRUMENT usually is FREQ_RESO
        self.dintint = dict()
        #- Dict keeping integrated intensities of sphinx models
        #- key: Transition(), value: list of values for every model included
        #- ~ self.star_grid
        self.mintint = dict()
        
    
    
    def calcIntInts(self,instrument='FREQ_RESO'):
        
        """
        Calculate the integrated intensity of datafiles as well as associated 
        sphinx files based on self.instruments[instrument] and self.star_grid
        respectively.
        
        The data intensities are stored in the dictionary self.dintint, the 
        model intensities in self.mintint. 
        
        They keys in the dictionaries are the Transition to which the sphinx 
        or datafiles belong. 
        
        The values are lists in accordance with self.instruments[instrument] 
        and self.star_grid for data and sphinx files respectively.
        
        @keyword instrument: The instrument for which this is done. For now 
                             only 'FREQ_RESO', ie spectrally resolved lines.
                             
                             (default: 'FREQ_RESO')
        @type instrument: string
        
        """
        
        instrument = instrument.upper()
        if instrument != 'FREQ_RESO':
            print 'Integrated line intensities cannot yet be calculated for '+\
                  'data other than spectrally resolved line such as APEX, ' + \
                  'JCMT, HIFI, etc.'
            return
        
        