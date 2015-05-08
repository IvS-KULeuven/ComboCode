# -*- coding: utf-8 -*-

"""
A toolbox for reading dust opacities in every shape and form.

Author: R. Lombaert

"""

import os
from scipy import array

from cc.tools.io import DataIO


class KappaReader(object):
    
    """
    An interface for reading opacities of multiple dust species.
    
    """
    
    def __init__(self,path_opac=os.path.join(os.path.expanduser('~'),\
                                             'MCMax','Opacities'),\
                 path_cc=os.path.join(os.path.expanduser('~'),'ComboCode')):
        
        """
        Initiating an instance of the KappaReader.
        
        @keyword path_opac: The path to the opacity home folder
        
                               (default: ~/MCMax/Opacities/)
        @type path_opac: str
        
        @keyword path_cc: Location of the ComboCode folder.
    
                          (default: ~/ComboCode/)
        @type path_cc: str
        """
        
        path_cc = os.path.join(path_cc,'Data')
        self.lspecies = DataIO.getInputData(path=path_cc,\
                                            keyword='SPECIES_SHORT',\
                                            filename='Dust.dat')
        self.lfilenames = DataIO.getInputData(path=path_cc,\
                                              keyword='PART_FILE',\
                                              filename='Dust.dat')
        self.lspec_dens = DataIO.getInputData(path=path_cc,\
                                              keyword='SPEC_DENS',\
                                              filename='Dust.dat')
        self.path_opac = path_opac
        self.kappas = dict()
        self.qext_a = dict()
        self.waves = dict()
    


    def readKappas(self,species):
        
        """
        Read kappas (cm2/g) and Q_ext/a (cm-1) for a dust species from the 
        MCMax INPUT files. 
        
        This also reads the absorption and scattering kappas separately. 
        
        @param species: The dust species (from Dust.dat)
        @type species: string
                        
        """
        
        if self.waves.has_key(species):
            return
        try:
            ispecies = self.lspecies.index(species)
        except ValueError:
            print 'Species not found in Dust.dat.'
            return
        fn = os.path.join(self.path_opac,self.lfilenames[ispecies])
        sd = self.lspec_dens[ispecies]
        if fn[-9:] == '.particle':
            part_file = DataIO.readFile(filename=fn,delimiter=' ') 
            wav = array([float(q[0]) 
                         for q in part_file if len(q) == 4]) 
            kappa = [array([float(q[1]) 
                            for q in part_file if len(q) == 4]),
                     array([float(q[2]) 
                            for q in part_file if len(q) == 4]),
                     array([float(q[3]) 
                            for q in part_file if len(q) == 4])]
        else: 
            part_file = DataIO.readCols(filename=fn)
            wav = part_file[0]
            kappa = part_file[1:]
        self.waves[species] = wav
        self.kappas[species] = kappa
        self.qext_a[species] = array(kappa) * 4/3. * sd
        
    
    
    def getKappas(self,species,index=0):
        
        """
        Return the kappas for given species.
        
        The index determines if you want extinction, scattering or absorption.
        
        @keyword index: The index of the kappas in the .opacity/.particle file. 
                        0: extinction, 1: absorption, 2: scattering
                        
                        (default: 0)
        @type index: int
        
        @return: (wavelength, kappas) [micron,cm2/g]
        @rtype: (array,array)
        
        """
        
        self.readKappas(species)
        if self.waves.has_key(species):
            return (self.waves[species],self.kappas[species][index])
        else:
            return (None,None)
        
    
    
    def getExtEff(self,species,index=0):
        
        """
        Return the extinction efficiencies per grain size for given species.
        
        The index determines if you want extinction, scattering or absorption.
        
        @keyword index: The index of the kappas in the .opacity/.particle file. 
                        0: extinction, 1: absorption, 2: scattering
                        
                        (default: 0)
        @type index: int
        
        @return: (wavelength, q_ext/a) [micron,cm-1]
        @rtype: (array,array)
        
        """
        
        self.readKappas(species)
        if self.waves.has_key(species):
            return (self.waves[species],self.qext_a[species][index])
        else:
            return (None,None)
        