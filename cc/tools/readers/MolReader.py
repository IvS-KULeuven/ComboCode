# -*- coding: utf-8 -*-

"""
A tool set for providing molecular data. Subclasses 
provide code-specific read and partse methods.

Author: R. Lombaert

"""

import numpy as np
from cc.tools.io import DataIO
from cc.tools.readers.SpectroscopyReader import SpectroscopyReader

from astropy import units as u


class MolReader(SpectroscopyReader):
    
    ''' 
    The MolReader class. 
    
    Inherits from the Reader class.
    
    Provides methods to manage molecular spectroscopic data (excitation 
    levels, Einstein A coefficients, transitions, etc). The only exception to 
    this is the collisional rate data. 
    
    Subclasses provide the read and parse methods that are code-specific: 
        - MlineReader for GASTRoNOoM output
        - RadiatReader for GASTRoNOoM input
        - LamdaReader for ALI/MCP input and Lamda data files in general (see
          http://home.strw.leidenuniv.nl/~moldata/).
        
    The level populations are read with PopReader for MCP/ALI, and MlineReader
    for GASTRoNOoM. The collision rates are read with CollisReader for 
    GASTRoNOoM, and LamdaReader for MCP/ALI.
    
    '''
    
    def __init__(self,fn,*args,**kwargs):
        
        ''' 
        Initialize an MolReader object. The filename and additional args/kwargs 
        are passed to the Reader parent class, where the dict is made and the 
        filename is stored. 
        
        Additional args/kwargs are used for the dict creation.
        
        Note that this is the only class directly inheriting from Reader that 
        does not have its own read method (unlike PopReader and CollisReader). 
        This is because it is always combined with a different type of class
        in multiple inheritance, depending on the input type of the spectroscopy
        (GASTRoNOoM versus Lamda)
        
        @param fn: The filename of molecular spectroscopic data. 
        @type fn: str
        
        '''
        
        super(MolReader,self).__init__(fn,*args,**kwargs)
        
        
    
    def getTFrequency(self,index=None,unit=u.GHz):
        
        '''
        Return the frequencies for the <nline> included transitions. 

        In case a single frequency is requested via index, the frequency is
        extracted from the array.
                
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        @keyword unit: The unit of the returned values. Can be any valid units 
                       str from the astropy units module (freq or wavelength), 
                       or the unit itself. 
                       In case of wavelength, the values are *not* reversed.
                        
                       (default: u.GHz)
        @type unit: string/unit
        
        @return: The frequencies/y ordered by transition index.
        @rtype: float/array
        
        '''
        
        #-- Get the frequencies in GHz
        freqs = self.get('trans','frequency',index)
        
        #-- Convert the units. Grab the unit first
        if isinstance(unit,str): unit = getattr(u,unit)
        return (freqs*u.GHz).to(unit,equivalencies=u.spectral()).value
    
    
    
    def getTEinsteinA(self,index=None):
    
        '''
        Return the Einstein A coefficients for the <nline> included transitions. 

        In case a single Einstein coefficient is requested via index, the value
        is extracted from the array.
                
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        
        @return: The Einstein A coefficient(s) (s^-1) ordered by transition 
                 index.
        @rtype: float/array
        
        '''
    
        return self.get('trans','einsteinA',index)
    
    
    
    def getTI(self,itype='trans',lup=None,llow=None):
        
        '''
        Return the indices of the transitions read from the molecular 
        spectroscopy data.
        
        A specific index (or array of indices) can be requested by specifying 
        the lower and upper level index.
        
        LamdaReader first inherits from MolReader, so the default will be itype
        'trans' for LamdaReader.
        
        @keyword itype: The type of index. Default is for MolReader. Needs to be
                        here in case LamdaReader calls this method, where the
                        type has to be specified. It will call both this method
                        and the version in CollisReader and pass on the call.
                        
                        (default: 'trans')
        @type itype: str 
        @keyword lup: The index of the upper energy level in the transition. If
                      both this and llow are None, all transition indices are 
                      returned.
        
                      (default: None)
        @type lup: int
        @keyword llow: The index of the lower energy level in the transition. If
                       both this and lup are None, all transition indices are 
                       returned.
        
                       (default: None)
        @type llow: int
        
        @return: The transition indices
        @rtype: array
        
        '''
        
        return super(MolReader,self).getTI(itype=itype,lup=lup,llow=llow)
        
        

    def getTUpper(self,index=None,itype='trans'):
        
        '''
        Return the indices of the upper states of the <nline> included 
        transitions.
        
        These are NOT the quantum numbers! Especially for CO-type molecules, 
        the J-level is equal to level_index-1.

        In case a single upper level is requested via index, the level index is
        extracted from the array.
        
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        @keyword itype: The type of index. Default is for MolReader. 
                        CollisReader needs this to be 'coll_trans'.
                        
                        (default: 'trans')
        @type itype: str
                
        @return: The upper level indices/x ordered by transition index.
        @rtype: int/array
        
        '''
        
        return super(MolReader,self).getTUpper(itype=itype,index=index)
        
    
    
    def getTLower(self,index=None,itype='trans'):
    
        '''
        Return the indices of the lower states of the <nline> included 
        transitions.
        
        These are NOT the quantum numbers! Especially for CO-type molecules, 
        the J-level is equal to level_index-1.

        In case a single lower level is requested via index, the level index is
        extracted from the array.
                
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        @keyword itype: The type of index. Default is for MolReader. 
                        CollisReader needs this to be 'coll_trans'.
                        
                        (default: 'trans')
        @type itype: str
                
        @return: The lower level indices/x ordered by transition index.
        @rtype: int/array
        
        '''
        
        return super(MolReader,self).getTLower(itype=itype,index=index)


    
    def getLI(self):
        
        '''
        Return the indices of the excitation levels read from the molecular 
        spectroscopy data.
        
        @return: The level indices
        @rtype: array
        
        '''
        
        return self['level']['index']
        
    
    
    def getLWeight(self,index=None):
        
        '''
        Return the weights of the <ny_up+ny_low> excitation levels.

        In case a single weight is requested via index, the weight value is
        extracted from the array.
                        
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        
        @return: The weights ordered by level index.
        @rtype: float/array
        
        '''
        
        return self.get('level','weight',index)
        
        
    
    def getLEnergy(self,index=None,unit=1./u.cm):

        '''
        Return the energies of the <ny_up+ny_low> excitation levels.

        In case a single energy level is requested via index the energy value is
        extracted from the array.
                        
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        @keyword unit: The unit of the returned values. Can be any valid units 
                       str from the astropy units module (energy), or the unit 
                       itself. 'cm-1' and 'cm^-1' are accepted as well.
                        
                       (default: 1./u.cm)
        @type unit: string/unit
        
        @return: The energies ordered by level index.
        @rtype: float/array
        
        '''
        
        #-- Get the energies in cm^-1
        energies = self.get('level','energy',index)
        
        #-- Convert the units. Grab the unit first
        if isinstance(unit,str) and unit.lower() in ['cm-1','cm^-1']: 
            unit = 1./u.cm 
        elif isinstance(unit,str): 
            unit = getattr(u,unit)
        
        #-- In case of temperature, and extra step is needed
        if (isinstance(unit,u.Quantity) and unit.unit.is_equivalent(u.K)) \
                or (isinstance(unit,u.UnitBase) and unit.is_equivalent(u.K)):
            energies = (energies/u.cm).to(u.erg,equivalencies=u.spectral())
            return energies.to(unit,equivalencies=u.temperature_energy()).value
        else: 
            return (energies/u.cm).to(unit,equivalencies=u.spectral()).value
        