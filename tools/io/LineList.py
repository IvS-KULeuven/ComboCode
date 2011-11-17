# -*- coding: utf-8 -*-

"""
Toolbox for reading, parsing, producing line lists for numerous molecules, and
from different databases.

Author: R. Lombaert

"""

import os 
import re
import string

from cc.tools.io import DataIO
from cc.modeling.objects import Transition



class LineList():
    
    '''
    Reading, managing and parsing of line lists taken from online databases.
    
    '''
    
    def __init__(self,molecule,x_min,x_max,x_unit='micron',cdms=0,jpl=0,\
                 lamda=0,min_strength=None,max_exc=None,include_extra=0,\
                 path=os.path.join(os.path.expanduser('~'),'LineLists')):
        '''
        Creating a LineList object. 
        
        Any database from which a linelist is requested requires an entry for 
        given molecule and database in path/MOLECULE_DATABASE.dat .
        
        @param molecule: The molecule for which LineList is requested
        @type molecule: Molecule()
        @param x_min: Minimum value for reading data
        @type x_min: float
        @param x_max: Maximum value for reading data
        @type x_max: float
        
        @keyword x_unit: Unit of wave/freq (micron,ghz,mhz,hz)
                              
                         (default: 'micron')
        @type x_unit: string                      
        @keyword cdms: read molecule data from CDMS catalog
        
                       (default: 0)
        @type cdms: bool
        @keyword jpl: read molecule data from JPL catalog
        
                      (default: 0)
        @type jpl: bool
        @keyword lamda: read molecule data from LAMDA catalog
        
                        (default: 0)
        @type lamda: bool
        @keyword min_strength: if None all are included, else only lines with 
                               strengths above this value are included 
                               (log scale, nm2*Mhz)
                               
                               (default: None)
        @type min_strength: float
        @keyword max_exc: if None all are included, else only lines with 
                          excitation energies below this value are included 
                          (cm-1)
                          
                          (default: None)
        @type max_exc: float
        @keyword path: path where the line lists have been stored, names of 
                       files are molecule_DATABASE.dat 
                       
                       (default: '~/LineLists/')
        @type path: string
        @keyword include_extra: include extra information such as integrated 
                                line intensity and excitation energy                                
                                
                                (default: 0)
        @type include_extra: bool
        
        '''
        self.molecule = molecule
        self.c = 2.99792458e10          #in cm/s
        self.cdms = cdms
        self.jpl = jpl
        self.lamda = lamda
        self.x_unit = x_unit.lower()
        self.min_strength = min_strength
        self.max_exc = max_exc
        if self.x_unit == 'micron':
            self.x_min = self.c/(x_max*10**(-4))*10**(-6)
            self.x_max = self.c/(x_min*10**(-4))*10**(-6)
        elif self.x_unit == 'ghz':
            self.x_min = x_min*10**3
            self.x_max = x_max*10**3
        elif self.x_unit == 'mhz':
            self.x_min = x_min
            self.x_max = x_max
        elif self.x_unit == 'hz':
            self.x_min = x_min*10**-6
            self.x_max = x_min*10**-6
        else:
            raise IOError('Unit %s not yet supported, or wrongly spelled.' + \
                          'Check LineList.__init__.__doc__ for more info.')
        self.path = path
        self.input_unit = 'mhz'
        self.line_list = []
        self.include_extra = include_extra 



    def makeCatInt(self,numeral):
        
        '''
        Return an integer matching the given string from the catalog 
        (to deal with 'A#' numerals).
        
        @param numeral: The numeral to be converted to integer
        @type numeral: string
        
        '''
        
        try:
            return int(numeral)
        except ValueError:
            alpha = string.letters[:26]
            alphanumerals = [100+10*i for i,letter in enumerate(alpha)]
            return alphanumerals[alpha.index(numeral[0])] + int(numeral[1])



    def parseStandardCatalog(self,data,catstring):
        
        '''
        Parse Line Lists of standard format, such as for JPL or CDMS.
        
        @param data: The content of the inputfile, no replaced spaces or 
                       delimiter used when reading with DataIO.readFile!
        @type data: list[string]
        @param catstring: the name of the catalog to be included in the output,
                          for instance 'CDMS', or 'JPL'
        @type catstring: string
        
        '''
        
        data = [line 
                for line in data 
                if float(line[0:13]) <= self.x_max \
                    and float(line[0:13]) >= self.x_min]
        if self.min_strength <> None:
            data = [line 
                    for line in data 
                    if float(line[21:29]) >= self.min_strength]
        if self.max_exc <> None:
            data = [line 
                    for line in data 
                    if float(line[31:41]) <= self.max_exc]
        vib_pattern = re.compile(r'(v\d?=\d)')
        if self.include_extra:
            data = [[float(line[0:13]),\
                     line[61:63].strip() \
                        and self.makeCatInt(line[61:63]) \
                        or 0,self.makeCatInt(line[55:57]),\
                     line[57:59].strip() \
                        and self.makeCatInt(line[57:59]) \
                        or 0,\
                     line[59:61].strip() \
                        and self.makeCatInt(line[59:61]) \
                        or 0,\
                     line[73:75].strip() \
                        and self.makeCatInt(line[73:75]) \
                        or 0,\
                     self.makeCatInt(line[67:69]),\
                     line[69:71].strip() \
                        and self.makeCatInt(line[69:71]) \
                        or 0,\
                     line[71:73].strip() \
                        and self.makeCatInt(line[71:73]) \
                        or 0,\
                     vib_pattern.search(line[81:len(line)]) \
                        and vib_pattern.search(line[81:len(line)]).groups()[0]\
                        or '',\
                     catstring,\
                     float(line[21:29]),\
                     float(line[31:41])] 
                    for line in data]
        else:
            data = [[float(line[0:13]),\
                     line[61:63].strip() \
                        and self.makeCatInt(line[61:63]) \
                        or 0,\
                     self.makeCatInt(line[55:57]),\
                     line[57:59].strip() \
                        and self.makeCatInt(line[57:59]) \
                        or 0,\
                     line[59:61].strip() \
                        and self.makeCatInt(line[59:61]) \
                        or 0,\
                     line[73:75].strip() \
                        and self.makeCatInt(line[73:75]) \
                        or 0,\
                     self.makeCatInt(line[67:69]),line[69:71].strip() \
                        and self.makeCatInt(line[69:71]) \
                        or 0,\
                     line[71:73].strip() \
                        and self.makeCatInt(line[71:73]) \
                        or 0,\
                     vib_pattern.search(line[81:len(line)]) \
                        and vib_pattern.search(line[81:len(line)]).groups()[0]\
                        or '',\
                     catstring] 
                    for line in data]
        return data
        


    def readCDMS(self):
        
        '''
        Read data from CDMS line list catalogs for a specific molecule.
        
        '''
        
        data = DataIO.readFile(os.path.join(self.path,\
                                          self.molecule.molecule+'_CDMS.dat'),\
                               replace_spaces=0)
        print 'Reading data from CDMS database for %s.'%self.molecule.molecule
        uncertainties = [float(line[13:21]) for line in data]
        if min(uncertainties) < 0 and max(uncertainties) == 0:
            input_xmin = self.x_min
            input_xmax = self.x_max
            self.x_min = (self.c/(input_xmin*10**6))**-1
            self.x_max = (self.c/(input_xmax*10**6))**-1
            self.input_unit = 'cm-1'
        elif min(uncertainties) < 0 and max(uncertainties) > 0:
            raise ValueError('Uncertainties in CDMS input file for ' + \
                             'molecule %s are ambiguous.'\
                             %self.molecule.molecule)
        data = self.parseStandardCatalog(data,'CDMS')
        if self.input_unit == 'cm-1':
            data = sorted([[i == 0 \
                                and self.c*entry*10**-6 \
                                or entry \
                            for i,entry in enumerate(line)] 
                           for line in data])
        self.line_list.extend(data)
        


    def readJPL(self):
        
        '''
        Read data from JPL line list catalogs for a specific molecule.
        
        '''
        
        data = DataIO.readFile(os.path.join(self.path,\
                                           self.molecule.molecule+'_JPL.dat'),\
                               replace_spaces=0)
        print 'Reading data from JPL database for %s.'%self.molecule.molecule
        data = self.parseStandardCatalog(data,'JPL')
        self.line_list.extend(data)



    def readLAMDA(self):
        
        '''
        Read data from LAMDA line list catalogs for a specific molecule.
        
        NOT YET IMPLEMENTED!
        
        '''
        print 'LAMDA catalog not yet supported.'
    


    def loadLineList(self,include_extra=0):
        
        '''
        Start reading all requested line list data.
        
        '''
        
        if self.cdms:
            try:
                self.readCDMS()
            except IOError:
                print 'CDMS line list not found for %s in %s.' \
                      %(self.molecule.molecule,self.path)
        if self.jpl:
            try:
                self.readJPL()
            except IOError:
                print 'JPL line list not found for %s in %s.' \
                      %(self.molecule.molecule,self.path)
        if self.lamda:
            try:
                self.readLAMDA()
            except IOError:
                print 'LAMDA line list not found for %s in %s.' \
                      %(self.molecule.molecule,self.path)
        self.line_list = sorted(self.line_list)
    


    def getLineList(self):
        
        '''
        Return the loaded line list, with frequencies in MHz.
        
        '''
        
        if not self.line_list:
            self.loadLineList()
        return self.line_list
        


    def makeTransitions(self,telescope=None,offset=None,n_quad=None,\
                        use_maser_in_sphinx=None):
        
        '''
        Make a list of transitions from the line list that was read.
        
        If any of the extra transition parameters is None, the default in 
        Transition.Transition.__init__() is used.
        
        @keyword telescope: The telescope that observed given line 
        
                            (default: None)
        @type telescope: string
        @keyword offset: The offset from center pixel of observed transition
                         
                         (default: None)
        @type offset: float
        @keyword n_quad: The number of grid points in the formal integration 
                         when calculating the line profile in sphinx.
                         
                         (default: None)
        @type n_quad: int
        @keyword use_maser_in_sphinx: Allow masering in sphinx calculations
        
                                      (default: None)
        @type use_maser_in_sphinx: bool
        @return: The transitions
        @rtype: list[Transition]
        
        '''
        pars = dict([(parname,par) 
                     for par,parname in zip([telescope,offset,\
                                             n_quad,use_maser_in_sphinx],\
                                            ['telescope','offset',\
                                             'n_quad','use_maser_in_sphinx']) 
                     if par <> None])
        return [Transition.Transition(molecule=self.molecule,\
                                      frequency=float(trans[0])*10**6,\
                                      exc_energy=float(trans[12]),\
                                      int_intensity_log=float(trans[11]),\
                                      vup=int(trans[1]),\
                                      jup=int(trans[2]),\
                                      kaup=int(trans[3]),\
                                      kcup=int(trans[4]),\
                                      vlow=int(trans[5]),\
                                      jlow=int(trans[6]),\
                                      kalow=int(trans[7]),\
                                      kclow=int(trans[8]),\
                                      vibrational=trans[9],\
                                      **pars) 
                for trans in self.getLineList()]
                
