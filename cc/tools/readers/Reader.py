# -*- coding: utf-8 -*-

"""
A tool set for reading all kinds of output/input/data.

Author: R. Lombaert

"""

from cc.tools.io import DataIO


class Reader(dict):
    
    ''' 
    The Reader class.
    
    Inherits from the builtin dictionary class, and thus functions much like a 
    dictionary. For now, no base methods from dict are overwritten.
    
    Reader functions primarily as a final stop before the dict class is called. 
    Any Reader objects must inherit from this (directly or indirectly). 
    
    Can read in the data/input/outputfiles that are given to this class and 
    store them.
    
    Classes inheriting from Reader usually implement their own read and write
    methods, but some basic capabilities are available here.
    
    '''
    
    def __init__(self,fn,*args,**kwargs):
        
        ''' 
        Initialize an Reader object by setting the contents dict.
        
        Additional args & kwargs passed to __init__ are passed to dict __init__. 
        
        @param fn: The filename of the file that is being read. 
        @type fn: str
        
        '''
        
        #-- Create the dictionary instance 
        super(Reader,self).__init__(*args,**kwargs)
        
        #-- Remember the filename (possibly including a wildcard character)
        self['fn'] = fn
        self.fn = self['fn']
        
        #-- The contents contains the raw files read with DataIO.readFile.
        self['contents'] = dict()
        
    
    
    def readFile(self,wildcard='*',*args,**kwargs):
        
        '''
        Read a filename and store its contents in the Reader object.
        
        The contents are stored in a dictionary, with the filename as key. Note 
        that one filename can refer to multiple files with the use of a wildcard
        character. This character can be replaced upon calling this method. 
        
        The contents can be returned with the method getFile.
        
        If the file is not found, an empty list is stored instead.
        
        Additional args/kwargs are passed to DataIO.readFile (such as delimiter
        and replace_spaces)
        
        @keyword wildcard: if a wildcard character is present in the filename, 
                           it can be replaced here.
                           
                           (default: '*')
        @type wildcard: string
    
        '''
        
        #-- Replace the wildcard if it is present.
        fn = self.fn.replace('*',wildcard)
        
        #-- Read the file
        self['contents'][fn] = DataIO.readFile(fn,*args,**kwargs)
    
    
    
    def getFile(self,wildcard='*',*args,**kwargs):
        
        '''
        Return the contents of the file of this Reader instance.
        
        Wildcard characters can be replaced. 
        
        Additional keywords can be passed to the readFile method in case the 
        file was not read yet.
        
        @keyword wildcard: if a wildcard character is present in the filename, 
                           it can be replaced here.
                           
                           (default: '*')
        @type wildcard: string
        
        @return: a list of lines in the file, each line represented by a list 
                 of strings, which were seperated by the chosen delimiter, if 
                 applicable. 
        @rtype: list
        
        '''
        
        #-- Replace the wildcard if it is present.
        fn = self.fn.replace('*',wildcard)
        
        #-- Check if the file was already read. 
        if not self['contents'].has_key(fn):
            self.readFile(wildcard,*args,**kwargs)
            
        return self['contents'][fn]
    