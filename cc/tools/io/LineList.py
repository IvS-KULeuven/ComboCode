# -*- coding: utf-8 -*-

"""
Toolbox for reading, parsing, producing line lists for numerous molecules, and
from different databases.

Author: R. Lombaert

"""

import os 
import re
import string

import cc.path
from cc.tools.io import DataIO
from cc.modeling.objects import Transition



class LineList():
    
    '''
    Reading, managing and parsing of line lists taken from online databases.
    
    '''
    
    def __init__(self,molecule,x_min,x_max,x_unit='micron',cdms=0,jpl=0,\
                 lamda=0,min_strength=None,max_exc=None,include_extra=0):
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
        
        data = DataIO.readFile(os.path.join(cc.path.ll,\
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
        
        data = DataIO.readFile(os.path.join(cc.path.ll,\
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
                self.cdms = 0
                print 'CDMS line list not found for %s in %s.' \
                      %(self.molecule.molecule,cc.path.ll)
        if self.jpl:
            try:
                self.readJPL()
            except IOError:
                print 'JPL line list not found for %s in %s.' \
                      %(self.molecule.molecule,cc.path.ll)
                self.jpl = 0
        if self.lamda:
            try:
                self.readLAMDA()
            except IOError:
                print 'LAMDA line list not found for %s in %s.' \
                      %(self.molecule.molecule,cc.path.ll)
        self.line_list = sorted(self.line_list)
    


    def getLineList(self):
        
        '''
        Return the loaded line list, with frequencies in MHz.
        
        '''
        
        if not self.line_list:
            self.loadLineList()
        return self.line_list
        


    def makeTransitions(self,telescope=None,offset=0.0,n_quad=100,\
                        fraction_tau_step=1e-2,min_tau_step=1e-4,\
                        write_intensities=0,tau_max=12.,tau_min=-6.,\
                        check_tau_step=1e-2,):
        
        '''
        Make a list of transitions from the line list that was read.
        
        Default transitions parameters are the same as those used in
        Transition.Transition.__init__().
        
        @keyword telescope: The telescope that observed given line 
        
                            (default: None)
        @type telescope: string
        @keyword offset: The offset from center pixel of observed transition
                         
                         (default: None)
        @type offset: float
        @keyword n_quad: The number of grid points in the formal integration 
                         when calculating the line profile in sphinx.
                         
                         (default: 100)
        @type n_quad: int
        @keyword fraction_tau_step: tau_total*fraction_tau_step gives min. 
                                    delta_tau in strahl.f. If too low, 
                                    min_tau_step will be used.
        
                                    (default: 1e-2)
        @type fraction_tau_step: float
        @keyword min_tau_step: minimum of delta_tau in strahl.f
        
                               (default: 1e-4)
        @type min_tau_step: float
        @keyword write_intensities: set to 1 to write the intensities of first 
                                    50 impact-parameters at the end of sphinx
        
                                    (default: 0)
        @type write_intensities: bool
        @keyword tau_max: maximum optical depth used for the calculation of the 
                          formal integral
        
                          (default: 12.)
        @type tau_max: float
        @keyword tau_min: maximum optical depth used for the calculation of the
                          formal integral
        
                          (default: -6.)
        @type tau_min: float
        @keyword check_tau_step: check.par.in sphinx if step in tau not too 
                                 large
        
                                 (default: 0.01)
        @type check_tau_step: float
        
        @return: The transitions
        @rtype: list[Transition]
        
        '''

        pars = {'fraction_tau_step':fraction_tau_step,\
                'min_tau_step':min_tau_step,'tau_max':tau_max,\
                'write_intensities':write_intensities,'tau_min':tau_min,\
                'check_tau_step':check_tau_step,'n_quad':n_quad,\
                "offset":offset,'telescope':telescope}
        if self.molecule.spec_indices == 0 and self.cdms == 1:
            trans_list = [Transition.Transition(molecule=self.molecule,\
                                        frequency=float(trans[0])*10**6,\
                                        exc_energy=self.include_extra \
                                                        and float(trans[12]) \
                                                        or None,\
                                        int_intensity_log=self.include_extra \
                                                           and float(trans[11])\
                                                           or None,\
                                        vup=int(trans[3]),\
                                        jup=int(trans[2]),\
                                        vlow=int(trans[7]),\
                                        jlow=int(trans[6]),\
                                        vibrational=trans[9],\
                                        **pars) 
                          for trans in self.getLineList()]        
        else:
            trans_list = [Transition.Transition(molecule=self.molecule,\
                                        frequency=float(trans[0])*10**6,\
                                        exc_energy=self.include_extra \
                                                        and float(trans[12]) \
                                                        or None,\
                                        int_intensity_log=self.include_extra \
                                                           and float(trans[11])\
                                                           or None,\
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
                    
        return trans_list