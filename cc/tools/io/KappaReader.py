# -*- coding: utf-8 -*-

"""
A toolbox for reading dust opacities in every shape and form.

Author: R. Lombaert

"""

import os
import numpy as np
from numpy import array

import cc.path
from cc.tools.io import DataIO


class KappaReader(object):
    
    """
    An interface for reading opacities of multiple dust species.
    
    """
    
    def __init__(self):
        
        """
        Initiating an instance of the KappaReader.
        
        """
        
        self.lspecies = DataIO.getInputData(path=cc.path.usr,\
                                            keyword='SPECIES_SHORT',\
                                            filename='Dust.dat')
        self.lfilenames = DataIO.getInputData(path=cc.path.usr,\
                                              keyword='PART_FILE',\
                                              filename='Dust.dat')
        self.lspec_dens = DataIO.getInputData(path=cc.path.usr,\
                                              keyword='SPEC_DENS',\
                                              filename='Dust.dat')
        self.kappas = dict()
        self.qext_a = dict()
        self.waves = dict()
        self.fns = dict()
        


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
        fn = os.path.join(cc.path.mopac,self.lfilenames[ispecies])
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
        self.fns[species] = fn
        self.waves[species] = wav
        self.kappas[species] = kappa
        self.qext_a[species] = array(kappa) * 4/3. * sd
        
    
    
    def getWavelength(self,species):
    
        """
        Return the wavelength grid for a given species. 
        
        @param species: The dust species (from Dust.dat)
        @type species: string
        
        @return: wavelength (micron)
        @rtype: array
        
        
        """
        
        self.readKappas(species)
        if self.waves.has_key(species):
            return self.waves[species]
        else:
            return np.empty(0)
        
        
        
    def getKappas(self,species,index=0):
        
        """
        Return the kappas for given species.
        
        The index determines if you want extinction, absorption or scattering.
        
        @param species: The dust species (from Dust.dat)
        @type species: string
        
        @keyword index: The index of the kappas in the .opacity/.particle file. 
                        0: extinction, 1: absorption, 2: scattering
                        
                        (default: 0)
        @type index: int
        
        @return: kappas (cm2/g)
        @rtype: array
        
        """
        
        self.readKappas(species)
        if self.kappas.has_key(species):
            return self.kappas[species][index]
        else:
            return np.empty(0)
        
    
    
    def getExtEff(self,species,index=0):
        
        """
        Return the extinction efficiencies per grain size for given species.
        
        The index determines if you want extinction, scattering or absorption.
        
        @keyword index: The index of the kappas in the .opacity/.particle file. 
                        0: extinction, 1: absorption, 2: scattering
                        
                        (default: 0)
        @type index: int
        
        @return: q_ext/a [micron,cm-1]
        @rtype: array
        
        """
        
        self.readKappas(species)
        if self.qext_a.has_key(species):
            return self.qext_a[species][index]
        else:
            return np.empty(0)
        