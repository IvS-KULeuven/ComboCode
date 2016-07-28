# -*- coding: utf-8 -*-

"""
A tool set for reading spectroscopy transition-based data.

Author: R. Lombaert

"""

import numpy as np, collections

from cc.tools.io import DataIO
from cc.tools.readers.Reader import Reader


class SpectroscopyReader(Reader):
    
    ''' 
    The SpectroscopyReader class.
    
    Inherits from the Reader class, and thus functions much like a 
    dictionary. 
    
    SpectroscopyReader functions primarily as an intermediary class between 
    Reader and MolReader/CollisReader, both classes that make heavy use of 
    transition based indexing, and does relate to upper and lower level indices.
        
    Typically not used stand-alone. Refer to MolReader and CollisReader for
    practical uses. Both classes make use of this parent class' methods.
    
    '''
    
    def __init__(self,fn,*args,**kwargs):
        
        ''' 
        Initialize an SpectroscopyReader object.
        
        Additional args & kwargs passed to __init__ are passed to dict __init__. 
        
        @param fn: The filename of the file that is being read. 
        @type fn: str
        
        '''
        
        #-- Create the dictionary instance 
        super(SpectroscopyReader,self).__init__(fn,*args,**kwargs)
        
    
    
    def get(self,ptype,prop,index=None):
    
        '''
        Return a property for a given type of property possibly for a given 
        index. 
        
        In case a single value is requested via index, the property value is
        extracted from the array.
        
        @param ptype: The type of property. 'coll_trans', 'trans' or 'level'.
        @type ptype: str
        @param prop: The property itself. eg 'energy' for level, or 'lup' for 
                     trans, or 'rates' for coll_trans.
        @type prop: str
        
        @keyword index: The index. In case of default, all are returned. Can be
                        any array-like object that includes indices
                        
                        (default: None)
        @type index: int/array
        
        @return: The property sorted by property type index, or a single element
        @rtype: float/int/array
        
        '''
        
        #-- Return all if no index specified, otherwise set an iterable.
        if index is None:
            return self[ptype][prop]
            
        #-- Prefer explicit selection on indexing. Doesn't assume indexing in 
        #   files goes 1 -> i_max 
        selection = self[ptype][prop][np.in1d(self[ptype]['index'],index)]
        
        #-- If a non-iterable object was passed as index, return just one value
        #   if only one value was indeed found. Otherwise, just return as is.
        if selection.shape != (1,) or isinstance(index,collections.Iterable):
            return selection
        else:
            return selection[0]
        
    
    
    def getTI(self,itype,lup=None,llow=None):
        
        '''
        Return the indices of the transitions read from the molecular 
        spectroscopy data or from the collisional rate data.
        
        A specific index (or array of indices) can be requested by specifying 
        the lower and upper level index.

        @param itype: The type of index. MolReader needs this to be 'trans'.
                      CollisReader needs this to be 'coll_trans'.
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
        
        #-- Return all if no index given
        if lup is None and llow is None: 
            return self[itype]['index']
        
        #-- Check for matches for given lup and llow. Will be empty array if no
        #   match found.
        bools = np.ones(shape=self[itype]['index'].shape,dtype=bool)
        if not lup is None: 
            bools *= self[itype]['lup'] == lup
        if not llow is None:
            bools *= self[itype]['llow'] == llow
        
        return self[itype]['index'][bools]
    
    
    
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
        @rtype: array
        
        '''
        
        return self.get(itype,'lup',index)
        
    
    
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
        @rtype: array
        
        '''
        
        return self.get(itype,'llow',index)