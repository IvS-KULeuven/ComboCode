# -*- coding: utf-8 -*-

"""
A class for reading and managing Mline output.

Author: R. Lombaert and M. van de Sande

"""

import os
import numpy as np

from astropy import constants as cst

from cc.tools.readers.PopReader import PopReader
from cc.tools.readers.MolReader import MolReader
from cc.tools.io import DataIO


class MlineReader(MolReader,PopReader):
    
    '''
    A Reader for Mline output files. 
    
    This class inherits from both PopReader and MolReader, and thus provides 
    spectroscopic information of the molecule as given in the mline output. In
    principle this information should be identical to what is read by 
    RadiatReader from the GASTRoNOoM spectroscopic input files. This information
    is read from ml1.
    
    In addition,
    mline output gives other information calculated by the radiative transfer, 
    and that is available here as well. The level populations specifically are
    managed through methods in PopReader. This class then also provides 
    additional data such as source function, optical depth, etc. This 
    information is read from ml3. 
    
    MlineReader only depends on GASTRoNOoM output and does not require 
    information from the databases or input files.
    
    '''
    
    def __init__(self,fn,*args,**kwargs):
        
        '''
        Creating an Mline object ready for reading Mline output.
        
        Reading mline output through full filename, including filepath.
        
        The mlinefile number is given as *, which is then replaced as needed.
        
        The wildcard character does not have to be present, and is inserted 
        automatically. It is handled down the line by the parent classes.
        
        Additional args and kwargs are passed to the dict creation (parent of 
        Reader)
        
        Note that properties based on the filename are saved as instance 
        properties (such as id, molecule name, filename).
        
        @param fn: The mline filename, including filepath. The mline
                   file number is given by *.
        @type fn: string        
        
        '''
        
        #-- Insert wildcard character to encompass all mline output files
        fn = fn.replace('ml1','ml*').replace('ml2','ml*').replace('ml3','ml*')
        super(MlineReader, self).__init__(fn=fn,*args,**kwargs)
        
        #-- Save molecule and id info.
        self.molecule = os.path.splitext(self.fn)[0].split('_')[-1]
        self.id = os.path.splitext(self.fn)[0].split('_')[-2].replace('ml*','')
                                                 
        #-- Read the files. Note that parent classes do this also, but only if 
        #   they are not used as parent class through inheritance, and thus work
        #   standalone
        self.read()
    
    
    
    def read(self): 
    
        '''
        Read the mline output files ml1 and ml3. 
        
        The information is stored in the instance dictionary. 
        
        '''
        
        #-- Read the main output files.
        self.__readML1()
        self.__readML3()

        
        
    def __readML1(self):
    
        '''
        Read molecule spectroscopic properties and basic circumstellar radial 
        profiles from the ml1 file.
        
        This is the information that is returned by MolReader methods.
        
        Note that the level indexing is J + 1 for simple molecules like CO. For
        any molecule, the 0-energy level has index 1!
        
        '''
        
        #-- Molecule settings are in ml1
        fn = self.fn.replace('ml*','ml1')
        
        #-- First block contains keywords, so use readDict
        d1 = DataIO.readDict(filename=fn,start_row=8,end_row=38,\
                             convert_floats=1,convert_ints=1,\
                             key_modifier='lower')
        self['pars'] = d1
        
        #-- Number of header lines: number of pars, 9 lines in first header, 1
        #   line with column names.
        nhdr = len(d1) + 9 + 1
        
        #-- Second block contains transitions
        d2 = np.genfromtxt(fn,usecols=[0,1,2,3,4],skip_header=nhdr,\
                           dtype='i8,i8,i8,f8,f8',max_rows=d1['nline'],\
                           names=['index','lup','llow','frequency','einsteinA'])
        
        #-- Add 1 to the lup and llow columns so they match up with the indices
        #   of the levels (for some reason j = 0 level is given as index 0 for 
        #   llow in the transition definition, while the level index is actually
        #   1). While for J-numbers this works, this approach breaks down for eg
        #   H2O.
        d2['lup'] += 1 
        d2['llow'] += 1 
        self['trans'] = d2
        
        #-- Third block contains levels (note hdr contains line with col names)
        nhdr += d1['nline'] + 1
        d3 = np.genfromtxt(fn,skip_header=nhdr,max_rows=d1['ny'],\
                           dtype='i8,f8,f8',names=['index','weight','energy'])
        
        self['level'] = d3
        
        #-- Fourth block contains circumstellar properties as a function of 
        #   impact parameter (note that P1 = R_outer, Pmax = R_star, so reverse)
        nhdr += d1['ny']
        colnames = ['p_rstar','vel','nmol','nh2','amol','Tg','Td']
        d4 = np.genfromtxt(fn,names=colnames,\
                           skip_header=nhdr,max_rows=d1['n_impact'])[::-1]
        self['props'] = d4
        
        #-- Add a convenient reference to the instance dict so PopReader and
        #   MolReader parent methods can work with the impact parameter grid.
        #   R_STAR is given in cm.
        self['p'] = d4['p_rstar'] * self['pars']['r_star']
        
        
        
    def __readML3(self):
    
        '''
        Read the information from the ML3 file. 
        
        This includes in six blocks, with subblocks for each impact parameter:
            - Scattering integral (si)
            - Source function (sf)
            - Line opacities (lo): if negative (or 1e-60 when 
              use_no_maser_option=1) the line is masing. 
            - Number densities (pop)
            - Derivative of scat int wrt line opacity times line opacity at line
              center (DsiDloXlo)
              (see GASTRoNOoM src code source_mline/lamit.f for details: D1)
            - Derivative of scat int wrt src function times src function at line
              center (DsiDsfXsf)
              (see GASTRoNOoM src code source_mline/lamit.f for details: D2)
        
        '''
                
        #-- 1) First three blocks list values for each transition and each 
        #      impact parameter. No white lines. If a white line occurs, it is 
        #      because nline%8 == 0. Hence number of text lines is 
        #      int(nline)/8+1 per impact parameter. Number of sub-blocks is then
        #      the number of impact parameters. Total length of the block is 
        #      nline_l*n_impact. Data starts at previous n_block + 3*N_block, 
        #      with N_block in [0->5] for [scat int, sf, line opac, pop, DsiDlo,
        #      DsiDsf]. 3 stands for 2 whitespace lines and one line with dashes
        #
        #   2) Data is fixed format, with 14 characters per entry (cannot use
        #      spaces as delimiters, because they might be filled with - sign)
        #
        #   3) Data is stored based on level index (for level pops), or trans
        #      index (for everything else), and each index gives an array as a 
        #      function of impact parameter, which is in units of R_STAR, and is
        #      given in self['props']['P']. Converted to cgs it is 
        #      given in self['props']['P_cm']. 
        #
        #-- Define the function that reads the fixed-format file per block
        def readBlock(N_block):
            
            '''
            Read a block of data from ml3. 
            
            @param N_block: The index of the block (0 => 5)
            @type N_block: int
            
            @return: The data 
            @rtype: array
            
            '''
            
            #-- Number of lines for this block
            nblock = nblock_ny if N_block == 3 else nblock_nline
            
            #-- Set starting index, based on number of block and block length
            ind0 = (N_block-1)*nblock_nline + 3*N_block
            
            #-- If N_block > 3, one block has only ny=ny_up+ny_low values, 
            #   otherwise add another block with nline values per subblock
            if N_block > 3: 
                ind0 += nblock_ny
            else: 
                ind0 += nblock_nline

            #-- Read data and return. fn defined in mother function
            return np.genfromtxt(fn,skip_header=ind0,max_rows=nblock,\
                                 delimiter=[14]*8)
        
        #-- Grab filename and set the contents values for the six blocks.
        fn = self.fn.replace('ml*','ml3')
        props = ['si','sf','lo','pop','DsiDloXlo','DsiDsfXsf']
        for pr in props: 
            self[pr] = dict()
        
        #-- Gather some relevant parameters. Length of block based on n_impact 
        #   and number of values (nline or ny). Assuming 8 columns.
        n_impact = self['pars']['n_impact']
        nline = self['pars']['nline']
        nline_l = int(nline)/8+1
        nblock_nline = n_impact * nline_l
        ny = self['pars']['ny']
        ny_l = int(ny)/8+1
        nblock_ny = n_impact * ny_l
        
        #-- Read data, looping over the 6 blocks.         
        dds = [readBlock(N) for N in xrange(6)]
        
        #-- Set all the transition-specific information. Note indexing (0-based
        #   in python, 1-based in fortran). Note that the arrays are reversed to
        #   match the increasing impact parameter grid.
        for dd,pr in zip(dds,props):
            if pr == 'pop': continue
            for i in xrange(nline): 
                self[pr][i+1] = dd[i/8::nline_l,i%8][::-1]
                
        #-- Set the level populations. Note indexing (0-based in python, 1-based
        #   in fortran). Note that the arrays are reversed to match the 
        #   increasing impact parameter grid.
        for i in xrange(ny): 
            self['pop'][i+1] = dds[3][i/8::ny_l,i%8][::-1]

    
    
    def getProp(self,prop):
    
        '''
        Return a radial circumstellar profile associated with the impact 
        parameter grid for the molecule of this mline model. 
        
        The impact parameter grid in cm is available through getP().
        
        @param prop: The requested property. One of p_rstar, vel (km/s), nmol
                     (cm^-3), nh2 (cm^-3), amol, Tg (K), Td (K).
        @type prop: str
        
        @return: The abundance profile
        @rtype: array
        
        '''
        
        return self['props'][prop]
        