# -*- coding: utf-8 -*-

"""
Toolbox for reading, parsing, producing line lists for numerous molecules, and
from different databases.

Author: R. Lombaert

"""

import os 
import re
import string
from astropy import units as u

import cc.path
from cc.tools.io import DataIO
from cc.modeling.objects import Transition



class LineList():
    
    '''
    Reading, managing and parsing of line lists taken from online databases.
    
    '''
    
    def __init__(self,fn,x_min=1e-10,x_max=1e10,unit='micron',\
                 min_strength=None,max_exc=None):
        '''
        Creating a LineList object. 
        
        @param fn: The full path and filename to the spectroscopy file. Must 
                   contain either the CDMS or the JPL string.
        @type fn: str
        
        @keyword x_min: Minimum value for reading data in frequency/wavelength
                        If default, the lowest frequency is the one in the file.
                        
                        (default: 1e-10)
        @type x_min: float
        @keyword x_max: Maximum value for reading data in frequency/wavelength
                        If default, the highest frequency is the one in the file
                      
                        (default: 1e10)
        @type x_max: float
        @keyword unit: Unit of x_min/x_max (micron,GHz,MHz,Hz). Any 
                       astropy.constants unit works. This can also be a Unit()
                       class object.
                              
                       (default: 'micron')
        @type unit: string                      
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
        
        '''
        
        self.fn = fn
        self.min_strength = min_strength
        self.max_exc = max_exc
        
        #-- Figure out if it's a JPL or a CDMS catalog
        if self.fn.upper().find('JPL') != -1: 
            self.catstring = 'JPL'
        else:
            self.catstring = 'CDMS'
        
        #-- Grab the unit
        if isinstance(unit,str) and unit.lower() in ['1 / cm','cm-1','cm^-1']: 
            unit = u.Unit("1 / cm")
        elif isinstance(unit,str): 
            unit = getattr(u,unit)
            
        #-- Convert the units to the default MHz (unit of the files), but 
        #   reverse min/max in case a wavelength is given. Wave number/Frequency
        #   are OK
        if unit.is_equivalent(u.micron):
            self.x_min = (x_max*unit).to(u.MHz,equivalencies=u.spectral())
            self.x_max = (x_min*unit).to(u.MHz,equivalencies=u.spectral())
        else:
            self.x_max = (x_max*unit).to(u.MHz,equivalencies=u.spectral())
            self.x_min = (x_min*unit).to(u.MHz,equivalencies=u.spectral())
        
        #-- Initialise the internal line list and read.
        self.line_list = []
        self.read()
        


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



    def __parseCatalog(self,data):
        
        '''
        Parse Line Lists of standard format, such as for JPL or CDMS.
        
        @param data: The content of the inputfile, no replaced spaces or 
                     delimiter used when reading with DataIO.readFile!
        @type data: list[string]
        
        '''
        
        data = [line 
                for line in data 
                if float(line[0:13]) <= self.x_max.value \
                    and float(line[0:13]) >= self.x_min.value]
        if self.min_strength:
            data = [line 
                    for line in data 
                    if float(line[21:29]) >= self.min_strength]
        if self.max_exc:
            data = [line 
                    for line in data 
                    if float(line[31:41]) <= self.max_exc]
        vib_pattern = re.compile(r'(v\d?=\d)')
    
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
                 self.catstring,\
                 float(line[21:29]),\
                 float(line[31:41])] 
                for line in data]
        return data
        


    def __readCDMS(self):
        
        '''
        Read data from CDMS line list catalogs for a specific molecule.
        
        '''
        
        data = DataIO.readFile(self.fn,\
                               replace_spaces=0)
        print 'Reading data from CDMS database for'
        print self.fn
        
        #-- If the uncertainties are negative, change the unit of min/max to 
        #   cm-1
        uncertainties = [float(line[13:21]) for line in data]
        if min(uncertainties) < 0 and max(uncertainties) == 0:
            self.x_min = self.x_min.to(1./u.cm,equivalencies=u.spectral())
            self.x_max = self.x_max.to(1./u.cm,equivalencies=u.spectral())
        elif min(uncertainties) < 0 and max(uncertainties) > 0:
            raise ValueError('Uncertainties in CDMS input file for ' + \
                             'file %s are ambiguous.'\
                             %self.fn)

        data = self.__parseCatalog(data)

        #-- If unit was changed, change the f values to MHz, the default unit
        rcm = u.Unit("1 / cm")
        if self.x_min.unit == rcm:
            data = sorted([[(entry*rcm).to(u.MHz,equivalencies=u.spectral()) 
                                if not i else entry
                            for i,entry in enumerate(line)] 
                           for line in data])
        self.line_list = data
        


    def __readJPL(self):
        
        '''
        Read data from JPL line list catalogs for a specific molecule.
        
        '''
        
        data = DataIO.readFile(self.fn,\
                               replace_spaces=0)
        print 'Reading data from JPL database for'
        print self.fn
        data = self.__parseCatalog(data)
        self.line_list = data



    def read(self):
        
        '''
        Start reading all requested line list data.
        
        '''
        
        self.line_list = []
        if self.catstring == 'CDMS':
            self.__readCDMS()
        elif self.catstring == 'JPL':
            self.__readJPL()



    def getLineList(self):
        
        '''
        Return the loaded line list, with frequencies in MHz.
        
        '''
        
        return self.line_list
        


    def makeTransitions(self,molecule,telescope=None,offset=0.0,n_quad=100,\
                        fraction_tau_step=1e-2,min_tau_step=1e-4,\
                        write_intensities=0,tau_max=12.,tau_min=-6.,\
                        check_tau_step=1e-2,):
        
        '''
        Make a list of transitions from the line list that was read.
        
        Default transitions parameters are the same as those used in
        Transition.Transition.__init__().
        
        Requires a Molecule() object to be passed to the method call. 
        
        @param molecule: The molecule for these transitions
        @type molecule: Molecule()
        
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
        if molecule.spec_indices == 0:
            trans_list = [Transition.Transition(molecule=molecule,\
                                        frequency=float(trans[0])*10**6,\
                                        exc_energy=float(trans[12]),\
                                        int_intensity_log=float(trans[11]),\
                                        vup=int(trans[3]),\
                                        jup=int(trans[2]),\
                                        vlow=int(trans[7]),\
                                        jlow=int(trans[6]),\
                                        vibrational=trans[9],\
                                        **pars) 
                          for trans in self.getLineList()]        
        else:
            trans_list = [Transition.Transition(molecule=molecule,\
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
                    
        return trans_list