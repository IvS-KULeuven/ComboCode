# -*- coding: utf-8 -*-

"""
Methods for column density calculation.

Author: R. Lombaert

"""

import os
from scipy.integrate import trapz
from scipy import average, argmax
import math

from cc.tools.io import DataIO
from cc.data import Data


class ColumnDensity(object):
    
    """
    Environment to calculate column densities and similar properties for a 
    Star() object.
    
    """
    
    def __init__(self,star):
        
        """
        Initializing an instance of the ColumnDensity class. 
        
        To be calculated:
         - r_des: The minimum radius, ie where the dust species starts to exist.
                  This is usually where the density of the species rises to 
                  1% of the maximum density of the shell (in cm).
         - t_des: The temperature at the destruction radius r_des of the species
         - r_max: The maximum radius at which species exists, taken to be 
                  outer radius, or where the density of species drops below 
                  10**(-10) times the maximum density in the shell (in cm).
         - t_min: The temperature at the maximum radius r_max of the species
         - coldens: The column density in g/cm2 of the dust species between
                    r_des and r_max
        
        @param star: The model parameter set
        @type star: Star()
        
        """
        
        self.star = star
        dust_path = os.path.join(star.path_combocode,'Data')
        self.dustlist = DataIO.getInputData(path=dust_path,\
                                            keyword='SPECIES_SHORT',\
                                            filename='Dust.dat')
        self.dustmolar = DataIO.getInputData(path=dust_path,\
                                            keyword='MOLAR_WEIGHT',\
                                            filename='Dust.dat',make_float=1)
        self.r_solar = star.r_solar      #in cm
        self.m_solar = star.m_solar      #in g
        self.year = star.year            #in s 
        self.avogadro = 6.022e23         
        self.au = star.au
        self.r_des = dict()
        self.r_max = dict()
        self.t_min = dict()
        self.t_des = dict()
        self.r_min_cd = dict()
        self.r_max_cd = dict()
        self.compd = dict()
        self.rad = []
        self.coldens = dict()
        self.fullcoldens = dict()
        self.dustfractions = dict()
        self.readDustInfo()
        
    

    def readDustInfo(self):
        
        """
        Read all column densities, min/max temperatures and min/max radii for 
        the species involved in the MCMax model.
        
        """
        
        dens = self.star.getMCMaxOutput(keyword='DENSITY',\
                                        incr=self.star['NRAD']*\
                                             self.star['NTHETA'])
        dens = Data.reduceArray(dens,self.star['NTHETA'])
        temp = self.star.getMCMaxOutput(keyword='TEMPERATURE',\
                                        incr=self.star['NRAD']*\
                                             self.star['NTHETA'])
        temp = Data.reduceArray(temp,self.star['NTHETA'])
        compf = os.path.join(os.path.expanduser('~'),'MCMax',\
                             self.star.path_mcmax,'models',\
                             self.star['LAST_MCMAX_MODEL'],'composition.dat')
        comp = DataIO.readCols(compf)
        self.rad = comp.pop(0)*self.au
        self.r_outer = self.rad[-1]
        for species in self.star['DUST_LIST']:
            #- Save the actual density profile for this dust species, as well
            #- as calculating the full column density of a dust species.
            self.dustfractions[species] = comp.pop(0)
            self.compd[species] = self.dustfractions[species]*dens
            self.fullcoldens[species] = trapz(x=self.rad,y=self.compd[species])      
            
            #- Determine the column density from 90% of the dust species formed
            #- onward, based on the mass fractions!
            #- Not before, because the comparison with H2 must be made,
            #- and this will skew the result if not solely looking at where the
            #- dust has (almost) all been formed.
            #- We also save min amd max radii, for use with the H2 calculation
            a_species = self.star['A_%s'%species]
            maxdens = max(self.compd[species])
            mindens = maxdens*10**(-10)
            radsel = self.rad[(self.dustfractions[species]>0.9*a_species)*\
                              (self.compd[species]>mindens)]
            denssel = self.compd[species]\
                                [(self.dustfractions[species]>0.9*a_species)*\
                                 (self.compd[species]>mindens)]
            self.coldens[species] = trapz(x=radsel,y=denssel)        
            self.r_min_cd[species] = radsel[0]
            self.r_max_cd[species] = radsel[-1]
            
            #- Determine the actual destruction radius and temperature.
            #- Taken where the density reaches 1% of the maximum density
            #- (not mass fraction).
            self.r_des[species] = self.rad[self.compd[species]>(maxdens*0.01)][0]
            self.t_des[species] = temp[self.compd[species]>(maxdens*0.01)][0]
            
            #- e-10 as limit for minimum is ok, because if shell is 100000 R*
            #- the mass conservation dictates ~ (10^5)^2 = 10^10 (r^2 law) 
            #- decrease in density. Shells this big dont occur anyway.
            self.r_max[species] = self.rad[self.compd[species]>mindens][-1]
            self.t_min[species] = temp[self.compd[species]>mindens][-1]
    
    
    
    def dustFullColDens(self,species):
        
        """
        Calculate the full column density of a dust species in the shell. 
        
        This is NOT the value used for determining the abundance of the species
        
        @param species: The dust species
        @type species: string
        
        @return: The column density of the dust species in the full envelope
                 in g/cm2
        @rtype: float
        
        """
        
        if not float(self.star['A_%s'%species]):
            return 0
        try:
            ispecies = self.dustlist.index(species)
        except ValueError:
            print 'Dust species %s not found in Dust.dat.'%species
            return 0
        if not self.dustmolar[ispecies]:
            print 'No molar weight given for dust species %s in Dust.dat.'\
                  %species
            return 0
        return self.fullcoldens[species]
        
        
        
    def dustMolecAbun(self,species):
        
        """ 
        Calculate the molecular abundance of a dust species with respect to H2.
        
        @param species: the dust species (from Dust.dat)
        @type species: string
        
        @return: The molecular abundance with respect to H2 of the dust species
        @rtype: float
        
        """

        if not float(self.star['A_%s'%species]):
            return 0
        try:
            ispecies = self.dustlist.index(species)
        except ValueError:
            print 'Dust species %s not found in Dust.dat.'%species
            return 0
        if not self.dustmolar[ispecies]:
            print 'No molar weight given for dust species %s in Dust.dat.'\
                  %species
            return 0
        nspecies = self.coldens[species]*self.avogadro/self.dustmolar[ispecies]
        nh2 = self.hydrogenColDens(species)
        return nspecies/nh2
        
        
    
    def hydrogenColDens(self,species):
        
        """
        Calculate the column density of molecular hydrogen between two radial
        distances. 
        
        The value is derived from the total mass-loss rate, assuming everything
        is H2.
        
        @param species: the dust species for which the H2 col dens is 
                        calculated. Required for the radial information. 
        @type species: string
        
        @return: The molecular hydrogen column density between the given radii
                 (cm-2)
        @rtype: float
        
        """

        mdot_gas = float(self.star['MDOT_GAS'])*self.m_solar/self.year
        vexp_gas = float(self.star['VEL_INFINITY_GAS']) * 100000
        rin = self.r_min_cd[species]
        rout = self.r_max_cd[species]
        mh2 = 2.
        sigma = (1./rin-1./rout)*mdot_gas/vexp_gas/4./math.pi
        nh2 = sigma * self.avogadro / mh2
        return nh2    
    
    
        