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
    
    def __init__(self):
        
        """
        Initiating an instance of the KappaReader.
        
        """
        
        dust_path = os.path.join(os.path.expanduser('~'),'ComboCode','Data')
        self.lspecies = DataIO.getInputData(path=dust_path,\
                                            keyword='SPECIES_SHORT',\
                                            filename='Dust.dat')
        self.lfilenames = DataIO.getInputData(path=dust_path,\
                                              keyword='PART_FILE',\
                                              filename='Dust.dat')
        self.lspec_dens = DataIO.getInputData(path=dust_path,\
                                              keyword='SPEC_DENS',\
                                              filename='Dust.dat')
        self.kappas = dict()
        self.qext_a = dict()
        self.waves = dict()
    


    def readKappas(self,species):
        
        """
        Read kappas (cm2/g) and Q_ext/a (cm-1) for a dust species from the 
        MCMax INPUT files.
        
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
        fn = os.path.join(os.path.expanduser('~'),'MCMax','src',\
                          self.lfilenames[ispecies])
        sd = self.lspec_dens[ispecies]
        
        part_file = DataIO.readFile(filename=fn,delimiter=' ') 
        wav = array([float(wl[0]) 
                     for wl in part_file if len(wl) == 4]) 
        kappa = array([float(q[1]) 
                       for q in part_file if len(q) == 4]) 
        self.waves[species] = wav
        self.kappas[species] = kappa
        self.qext_a[species] = kappa * 4/3. * sd
        
    
    
    def getKappas(self,species):
        
        """
        Return the kappas for given species.
        
        @return: (wavelength, kappas) [micron,cm2/g]
        @rtype: (array,array)
        
        """
        
        self.readKappas(species)
        if self.waves.has_key(species):
            return (self.waves[species],self.kappas[species])
        else:
            return (None,None)
        
    
    
    def getExtEff(self,species):
        
        """
        Return the extinction efficiencies per grain size for given species.
        
        @return: (wavelength, q_ext/a) [micron,cm-1]
        @rtype: (array,array)
        
        """
        
        self.readKappas(species)
        if self.waves.has_key(species):
            return (self.waves[species],self.qext_a[species])
        else:
            return (None,None)
        