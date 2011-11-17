# -*- coding: utf-8 -*-

"""
A radiat file parser. 

Returns all info from Einstein coefficients and frequencies, 
to weights, energy levels and transition levels.

Author: R. Lombaert

"""

import os 

from cc.tools.io import DataIO
from cc.modeling.objects import Star



class Radiat():
    
    '''
    Class for working with radiat files from GASTRoNOoM's input. 
    
    '''
    
    def __init__(self,molecule,use_indices_dat=0,\
                 path=os.path.join(os.path.expanduser('~'),'GASTRoNOoM','src',\
                                   'data')):
        
        '''
        Initializing an instance of the Radiat class.
        
        @param molecule: The molecule (short hand name) for which radiat info is 
                         loaded.
        @type molecule: Molecule()
        @keyword use_indices_dat: If True, get the radiat filename from 
                                  ~/ComboCode/Data/Indices.dat
                                        
                                  (default: 0)
        @type use_indices_dat: bool
        @keyword path: The path where the radiat file can be found. Uses 
                       radiat_backup appended to the path in case use_indices_dat
                       is True.
                            
                       (default: ~/GASTRoNOoM/src/data/)
        @type path: string
        
        '''
        
        self.molecule = molecule
        self.path = path
        if use_indices_dat:
            f = Star.getInputData(path=os.path.join(\
                                        self.molecule.path_combocode,'Data'),\
                                  keyword='RADIAT',filename='Indices.dat',\
                                  start_index=4)[self.molecule.indices_index]
            self.filename = os.path.join(self.path,'radiat_backup',f)
        else:
            self.filename = os.path.join(self.path,\
                                        '%s_radiat.dat'%self.molecule.molecule)
        self.c = 2.99792458e10          #in cm/s
        self.__read()
        
        
        
    def __read(self):
         
        '''
        Read the MOLECULE_radiat file and set a dictionary for the instance
        with all the information. 
        
        Done on creation of an instance of the class.
        
        '''
        
        radiat = DataIO.readCols(filename=self.filename,make_array=0)[0]
        self.dict = {}
        #- The starting indices for the next step are determined incrementally, 
        #- with a right term added + arbitrary number of zeroes at the end of 
        #- each step in the loop
        step_sizes = [self.molecule.nline,self.molecule.nline,
                      self.molecule.ny_up+self.molecule.ny_low,
                      self.molecule.ny_up+self.molecule.ny_low,
                      self.molecule.nline,self.molecule.nline]
        if sum(step_sizes) == len(radiat):
            nozeroes = True
        else:
            nozeroes = False
        
        #- The ending indices are always the starting index + this index
        pars = ['EINSTEIN','FREQUENCY','WEIGHT','ENERGY','LOWER','UPPER']
        i = 0
        for delta,par in zip(step_sizes,pars):
            #- grab list from input
            self.dict[par] = radiat[i:i+delta]
            #- Determine starting index of the next list
            try:
                #- if no zeroes present: don't look for zeroes 
                #- (important cuz an energy level may be zero)
                i = nozeroes \
                        and (i + delta) \
                        or DataIO.findNumber(i+delta,radiat)     
            except IndexError:
                if par == 'UPPER':
                    #- This is the last step, you expect an index error here, 
                    #- since it's the end of the file
                    if not(len(radiat) == i+delta or radiat[i+delta] == 0):
                        #- If the next number is zero, or if the file ends here, 
                        #- everything is OK: Move on, otherwise raise error.
                        raise IndexError("When reading %s,"%self.filename + \
                                         "the level indices were " + \
                                         "formatted incorrectly. Aborting...")
                else:
                    #- Too early to get this error! Aborting!
                    raise IndexError('When reading %s,'%self.filename + \
                                     'file ended too soon before'+\
                                     ' definition of all parameters. '+\
                                     'Aborting...')
        
        #- if zeroes are present, one ofthe energy levels may also be really zero
        #- check this and correct for it here.
        if self.dict['ENERGY'][-1] == 0.0:    
            self.dict['ENERGY'] = self.dict['ENERGY'][0:-1]
            self.dict['ENERGY'][0:0] = [0.0]
        #Another check up...
        if sum(step_sizes) != sum([len(self.dict[par]) for par in pars]):
            raise IndexError('Fewer or more entries found for all the ' + \
                             'parameters in %s than expected. Aborting...'\
                             %self.filename)
        self.dict['LOWER'] = [int(x) for x in self.dict['LOWER']]
        self.dict['UPPER'] = [int(x) for x in self.dict['UPPER']]
             
             
             
    def getFrequency(self,unit='GHz'):
        
        '''
        Return a list of the frequencies for the <nline> included transitions. 
        
        @keyword unit: The unit of the returned values after reading
                       Can be: GHz, MHz, Hz, MICRON, MM, CM, M
                        
                       (default: 'GHz')
        @type unit: string
        
        '''
        
        if unit.upper() == 'GHZ':
            return self.dict['FREQUENCY']
        elif unit.upper() == 'MHZ':
            return [w*10.0**3 for w in self.dict['FREQUENCY']]
        elif unit.upper() == 'KHZ':
            return [w*10.0**6 for w in self.dict['FREQUENCY']]
        elif unit.upper() == 'HZ':
            return [w*10.0**9 for w in self.dict['FREQUENCY']]
        elif unit.upper() == 'MICRON':
            #- e9*e-4=e5 (GHz -> Hz, cm -> micron)
            return [self.c/(w*10.0**5) for w in self.dict['FREQUENCY']]                         
        elif unit.upper() == 'MM':
            #- e9*e-1=e8 (GHz -> Hz, cm -> mm)
            return [self.c/(w*10.0**8) for w in self.dict['FREQUENCY']]                         
        elif unit.upper() == 'CM':
            #- e9 (GHz -> Hz)
            return [self.c/(w*10.0**9) for w in self.dict['FREQUENCY']]                         
        elif unit.upper() == 'M':
            return [self.c/(w*10.0**11) for w in self.dict['FREQUENCY']]         
        
        
        
    def getEnergyLevels(self):
        
        '''
        Return a list of all energy levels for the <ny_up+ny_low> states 
        in the default unit of the radiat files.
        
        '''
        
        return self.dict['ENERGY']
    
    
    
    def getEinsteinCoefficients(self):
        
        '''
        Return a list of the <nline> Einstein coefficients.
        
        '''
        
        return self.dict['EINSTEIN']
    
    
    
    def getWeights(self):
        
        '''
        Return a list of the weights of the <ny_up+ny_low> states.
        
        '''
        
        return self.dict['WEIGHT']
    
    
    
    def getLowerStates(self):
        
        '''
        Return a list of the lower states of the <nline> included transitions.
        
        '''
        
        return self.dict['LOWER']
    
    
    
    def getUpperStates(self):
        
        '''
        Return a list of the upper states of the <nline> included transitions.
        
        '''
        
        return self.dict['UPPER']