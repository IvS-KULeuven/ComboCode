# -*- coding: utf-8 -*-

"""
Module for calculating the energy balance. 

Author: R. Lombaert, H. Olofsson & M. Maercker (Chalmers, Sweden)

For the use of the EnergyBalance module, please contact one of the three authors
listed here.

"""

import os, collections, functools, copy
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d
from scipy.interpolate import interp1d
from scipy.integrate import odeint, trapz, cumtrapz
from astropy import constants as cst
from astropy import units as u

import matplotlib.pyplot as plt

from cc.data import Data
from cc.plotting import Plotting2
from cc.tools.io import DataIO
from cc.tools.numerical import Operators as op
from cc.tools.numerical import Gridding
from cc.modeling.profilers import Velocity, Density, Profiler, Grainsize
from cc.modeling.profilers import Mdot, Radiance, Opacity, Temperature
from cc.tools.readers import CollisReader, PopReader, LamdaReader, MlineReader
from cc.tools.readers import RadiatReader
import cc.path

#-- Define a few global constants. Calling cst inside the functions slows them
#   down too much

#-- The Boltzmann constant
global k_b
k_b = cst.k_B.cgs.value # erg/K

#-- Hydrogen mass
global mh 
mh = cst.m_p.cgs.value # g



def dTdr(T,r,v,gamma,rates=None,warn=1):

    '''
    The differential equation for the kinetic temperature profile. 
    
    Currently only takes into account adiabatic expansion. 
    
    This function is used by EnergyBalance to iterate and solve the ODE.
    
    @param T: The temperature at which to evaluate the differential equation
    @type T: float
    @param r: The radial point(s) (cm)
    @type r: array/float
    @param v: The velocity profile object
    @type v: Velocity() 
    @param gamma: The adiabatic coefficient profile as function of T
    @type gamma: Profiler()
    
    @keyword rates: The heating and cooling rates, summed up (erg/s/cm3, H-C). 
                    Includes the density factor, without the velocity component.
                    This speeds up iteration, limiting the radial grid 
                    evaluation where possible. Default if only adiabatic cooling
                    is taken into account.
              
                    (default: None)
    @type rates: Profiler()
    @keyword warn: Warn when extrapolation occurs in an interpolation object. 
                   Not applicable to functional evaluation, or gamma in which 
                   case the extrapolation is deemed safe (constant value at 
                   low and high end temperatures consistent with diatomic gases)
    
                   (default: 1)
    @type warn: bool    
    
    @return: The derivative with respect to radius (K/cm)
    @rtype: array/float
    
    '''
    
    #-- Set the gamma-based factor. gamma is assumed to be a Profiler object.
    tg = gamma.eval(T,warn=0)
    factor1 = 2-2*tg
    factor2 = tg-1
    
    #-- Determine the velocity and its derivative from the Velocity() object
    vi = v.eval(r,warn=warn)
    dvi = v.diff(r,warn=warn)
    
    #-- Calculate dT/dr
    #   1) adiabatic cooling term
    term1 = factor1 * (1./r + 0.5*1./vi*dvi) * T
    
    #   2) other heating and cooling terms (can be none at all). This includes
    #      the density factor, except the velocity.
    rr = rates.eval(r,warn=warn) if rates else 0.
    term2 = factor2 * rr / vi
    
    #-- And the final summation of the two terms
    Tprime = term1 + term2
   
    return Tprime



class EnergyBalance(object):
    
    '''
    The energy balance class.
    
    Author: R. Lombaert, H. Olofsson & M. Maercker (Chalmers, Sweden)

    For the use of the EnergyBalance module, please contact one of the three 
    authors listed here.
    
    Calculates energy balance and provides tools for reading and writing I/O, 
    and setting input physics.
    
    For now limited to the energy balance from the dust formation radius up to
    the outer radius.
        
    An example for calculating the energy balance (see the respective functions
    for in-depth information): 
    >>> from cc.modeling.physics import EnergyBalance as EB
    >>> eb = EB.EnergyBalance()
    >>> eb.iterT()
    >>> eb.plotT()
    >>> eb.plotRates()
    This calculates and plots the temperature profiles and heating/cooling
    rates across all iterations using the default settings given in 
    cc.path.aux/inputEnergyBalance_standard.dat
    
    Default in the inputEnergyBalance_standard file can be adapted through 
    keywords passed to EnergyBalance. A few examples: 
    >>> eb = EB.EnergyBalance(a=[0.005,0.25,300,1])
    which creates a grain-size grid from 0.005 to 0.25 micron with 300 
    points distributed in log scale. 
    
    >>> eb = EB.EnergyBalance(hterms=['dg','dt'])
    adds dust-gas collisional heating and heat exchange to the adiabatic 
    cooling as heat/cool terms.
    
    >>> eb = EB.EnergyBalance(template='gastronoom')
    uses the input template to reproduce the GASTRoNOoM settings. 
    
    >>> eb = EB.EnergyBalance(gamma='h2')
    calculates the adiabatic coefficient from the temperature-dependent measured
    coefficient of molecular hydrogen.
    
    >>> eb = EB.EnergyBalance(heatmode='gs2014')
    calculates the heating and cooling terms from dust-gas interaction based on 
    the formalism of Gail & Sedlmayr 2014 as opposed to the classical 
    implementations of Decin et al. 2006 and Schoïer et al. 2001. 
    
    Furthermore, the iterative procedure can be adapted to user-specific 
    needs:
    >>> eb.iterT(conv=0.005,imax=200,step_size=0.05)
    sets the relative convergence criterion to 0.005 from default 0.01, the
    maximum iterations to 200, and the step_size of the gradual maximum 
    allowed temperature change between iterations to 5%.   
    
    '''
    
    def __init__(self,fn=None,template='standard',**kwargs):
        
        '''
        Creating an EnergyBalance object used to calculate the energy balance.
        
        Can work in conjunction with iteration with a RT code. 
        
        Takes into account several contributors to the heating and cooling of
        an outflow:
            - Adiabatic cooling
            - dust-gas collisional indirect heating (elastic collisions)
            - dust-gas direct heat exchange (accommodation)
            
        Not yet implemented:
            - Radiative line cooling and heating (for H2, CO, H2O, ...)
            - Photoelectric heating 
            - Heating by cosmic rays
        
        An instance of EnergyBalance contains three types of properties: 
            - basic input parameters given by the inputfile
            - independent variables and basic profiles set before calculations
            - profiles calculated from the basic parameters, variables, and 
              profiles, reset whenever any of the previous two properties are
              changed.
        
        First set is given in aux/inputEnergyBalance.dat.
        
        Second set includes: radius, grain size, gas velocity, 
        gas/dust mass-loss rates, dust opacities, dust temperature, stellar 
        radiance. 
        
        Third set includes: gas temperature, gas/dust density, drift velocity.
        
        Note that for now line cooling parameters are fixed. Hence after a new 
        RT model was calculated, a new energy balance must be calculated. In the 
        future, the level pops will be updateable. 
        
        @keyword fn: The parameter inputfile filename. If not given, a default
                     inputfile is read from aux/inputEnergyBalance.dat.
        
                     (default: aux/inputEnergyBalance.dat)
        @type fn: str
        @keyword template: The inputfile template used for initializing the 
                           energy balance. Default is the standard, with as much
                           consistency as possible. Alternatives are 'mcp' and
                           'gastronoom' to reproduce the respective code's 
                           results. Templates cannot be changed, but inputfiles
                           given by fn overwrite any input given there.
                           
                           (default: 'standard')
        @type template: str
        
        @keyword kwargs: Input parameters listed in the default 
                         inputEnergyBalance file can also be passed through the
                         initialisation of this instance. They will replace the 
                         values given in the inputfile.
                       
                         (default: {})
        @type kwargs: dict
                
        '''
        
        #-- A welcome message, and credentials
        m = "Welcome to the EnergyBalance module of ComboCode! \n"+\
            "This module was written at Chalmers University of Technology. "+\
            "The main author is R. Lombaert, with support from H. Olofsson "+\
            "and M. Maercker. Please contact one of these people if you wish "+\
            "to use the EnergyBalance with the purpose of publishing results."
        print m
        
        #-- Initialise all variables
        self.__reset()
        
        #-- Set default input template, and the custom parameter file
        self.template = template.lower()
        if self.template not in ['mcp','gastronoom']: self.template = 'standard'
        self.fn = fn
        self.include_lc = 0
            
        #-- Read input parameters, set constants, coordinate grids, stellar
        #   properties, basic profiles, adiabatic coefficient, line cooling,
        #   velocity 
        self.readInputParameters(**copy.deepcopy(kwargs))
        self.setGrids()
        self.setStar()
        self.setProfiles()
        self.setGamma()
        self.setVelocity()
        
        #-- Everything has been initialised. Set the current Temperature profile
        self.__setT()
        
        #-- Set the abundance for each molecule. The initial T has been set, so
        #   also calculate the line cooling terms if requested.
        for m in self.molecules: 
            self.setAbundance(m)
            if self.include_lc: self.setLineCooling(m)
    
    
    
    def formatInput(self,ilst):
    
        '''
        Format an input keyword list for a variable. 
        
        Can be given as a string, split by spaces, or as the list output from 
        str.split().
        
        Assumes the first element is the requested function or a constant, 
        followed by key=value pairs that serve as the input keyword arguments 
        for the function if the first element is not a constant. The first
        element thus does not contain a '=' sign.
        
        The key=value pairs are formatted into a dictionary.
        
        Example str: 
        'read_opacity species=AMCDHSPREI index=1 order=5'
        Example list:
        ['read_opacity','species=AMCDHSPREI','index=1','order=5']
        
        @param ilst: The input line, arguments split by spaces
        @type ilst: list[str]/str
        
        @return:
        @rtype: list[str,dict]
        
        '''
        
        #-- Check if a list or a str is given. Apply the ' ' split.
        if not isinstance(ilst,collections.Iterable) or isinstance(ilst,str):
            ilst = str(ilst).split()
        
        #-- If the list is empty, return an empty entry
        if not ilst: 
            return ['',{}]
            
        #-- Read the first element and parse the rest a dictionary. Either the
        #   first element gives a constant, or it gives a read method for files
        #   The rest are the function input arguments in key=value pairs
        ddk = {'convert_lists':1,'convert_floats':1,'convert_ints':1}
        return [ilst[0],DataIO.readDict(lines=ilst[1:],**ddk)]
        
        
        
    def readInputParameters(self,**kwargs):
    
        '''
        Read the input parameters for the energy balance. 
        
        If a filename is not given, the one set upon initialisation is used. If
        none was given upon initialisation the default one from aux/ is read.

        Resets any variables that depend on independent variables.
                
        @keyword fn: The parameter inputfile filename. If None, a default
                     inputfile is read from aux/
        
                     (default: None)
        @type fn: str

        @keyword kwargs: Input parameters listed in the default 
                         inputEnergyBalance file can also be passed through the
                         initialisation of this instance. They will replace the 
                         values given in the inputfile.
                       
                         (default: {})
        @type kwargs: dict
        
        '''
        
        #-- Only reset if self.T is already yet defined.
        if self.T: self.__reset()
        
        #-- Read the template
        fnbase = os.path.join(cc.path.aux,'inputEnergyBalance_')
        templates = {k: '{}{}.dat'.format(fnbase,k)
                     for k in ['standard','mcp','gastronoom']}
        self.pars = DataIO.readDict(templates[self.template],convert_lists=1,\
                                    convert_floats=1,multi_keys=['molecule'])
        if not self.pars.has_key('molecule') or self.pars['molecule'] == ['']:
            self.pars['molecule'] = []
            
        #-- Read from the filename if it is given. 
        if not self.fn is None: 
            dd = DataIO.readDict(self.fn,convert_lists=1,convert_floats=1,\
                                 multi_keys=['molecule'])
            m = dd.pop('molecule',[])
            self.pars.update(dd)
            self.pars['molecule'] += m
        
        #-- Insert any extra arguments given in function call
        m = kwargs.pop('molecule',[])
        if isinstance(m, str): m = [m] 
        self.pars.update(kwargs)
        self.pars['molecule'] += m    
        
        #-- Format the molecule information, and extract molecule names. 
        self.pars['molecule'] = [mlst.split() for mlst in self.pars['molecule']]
        abuns = []
        self.molecules = []
        for mlst in self.pars['molecule']:
            m = mlst.pop(0)
            if not m in self.molecules: 
                self.molecules.append(m)
                abuns.append(mlst)
        self.pars['molecule'] = abuns
        
        #-- Read left over information and parse it as a dictionary. This 
        #   contains the abundance information. Only used when no pops are given
        #   and for Hpe in case of 12C16O
        self.pars['molecule'] = [self.formatInput(mlst)
                                 for mlst in self.pars['molecule']]
        
        #-- Then initialise the heating and cooling terms requested
        self.H = {term: {} for term in self.pars['hterms']}
        self.C = {term: {} for term in self.pars['cterms'] if term != 'lc'}
        
        #-- Adiabatic cooling is always calculated
        self.C['ad'] = {}
        
        #-- Prepare the settings for the molecular abundances. Note that these
        #   abundance are taken from the model output when possible in case 
        #   level populations are used. If those are not available, the info
        #   given in the molecule input is used. 
        for m in self.molecules:
            self.abun = {}
        
        #-- Add the molecules for line cooling if applicable, and set the 
        #   containers for collision rates, level pops, spectroscopy, Texc, 
        #   abundances. Spectroscopy dict will refer to either collis or pop
        #   depending on the code (LamdaReader or MlineReader, respectively)
        if 'lc' in self.pars['cterms']:
            self.pars['cterms'] = [p for p in self.pars['cterms'] if p != 'lc']
            self.include_lc = 1
            self.collis = {}
            self.pop = {}
            self.Texc = {}
            self.ipop = {}
            self.mol = {}
            
            #-- Set the dict entries for each molecule. Also add a decorated 
            #   function for each molecule that adds the molecule as an 
            #   argument by default, so Clc can be called without arguments.
            #   (note that 'lc' does not occur in the self.C dict keys)
            for m in self.molecules: 
                fname = 'lc_{}'.format(m)
                self.C[fname] = {}
                self.pars['cterms'].append(fname)
                setattr(self,'C'+fname,functools.partial(self.Clc,m))
                
            #-- And also set the proper readers depending on the code. 
            if self.pars['rtcode'].lower() == 'gastronoom':
                self.colread = CollisReader.CollisReader
                self.popread = MlineReader.MlineReader
            else:
                self.colread = LamdaReader.LamdaReader
                self.popread = PopReader.PopReader
        
        #-- If no line cooling is requested, remove any and all population files
        else: 
            self.pars['pop'] = []
            
    
    
    def setGrids(self,r=None,a=None,l=None):
    
        '''
        Initialise the coordinate grids: radius (cm), grain size (cm), 
        wavelength (cm). 
        
        The grids can be given as lists, which are used as input for 
        Gridding.makeGrid. If not given, they are taken from the inputfile.

        @keyword r: The radius grid (cm). [rmin,rmax,nsteps,log?]
                    
                    (default: None)
        @type r: list        
        @keyword a: The grain size grid (cm). [amin,amax,nsteps,log?] Can be 
                    single value in case of one average grain size.
                    
                    (default: None)
        @type a: list/float
        @keyword l: The wavelength grid (cm). [lmin,lmax,nsteps,log?]
                    
                    (default: None)
        @type l: list      
          
        '''
        
        #-- Only reset if self.T is already yet defined.
        if self.T: self.__reset()
        
        #-- Take the grids from the inputfile if not given here
        if r is None: r = self.pars['r']
        if a is None: a = self.pars['a']
        if l is None: l = self.pars['l']
        
        #-- Check if a is a single value or a grid
        a = Data.arrayify(a)
        if a.size < 2:
            self.a = a
        else:
            self.a = Gridding.makeGrid(*a)
            
        #-- Set radius and wavelength
        self.r = Gridding.makeGrid(*r)
        self.l = Gridding.makeGrid(*l)
        
        

    def setStar(self,rstar=None,Tstar=None,Ltype=None):
        
        '''
        Set the stellar properties. 
        
        If any of these are None, they are taken from the inputEnergyBalance.dat
        file.
        
        Defines the Radiance profile based on these properties.
        
        @keyword rstar: The stellar radius (cm)
                        
                        (default: None)
        @type rstar: float
        @keyword Tstar: The stellar effective temperature (K)
                        
                        (default: None)
        @type Tstar: float
        @keyword Ltype: The type of stellar radiance profile (blackbody, ...)
                        
                        (default: None)
        @type Ltype: float
        
        '''
        
        #-- Only reset if self.T is already yet defined.
        if self.T: self.__reset()
        
        #-- Take the profile definitions from the inputfile if not given here
        if rstar is None: rstar = self.pars['rstar']
        if Tstar is None: Tstar = self.pars['Tstar']
        if Ltype is None: Ltype = self.pars['Ltype']
        
        self.rstar = rstar
        self.Tstar = Tstar
        
        #-- Set the stellar radiance profile and stellar surface area
        self.rad = Radiance.Radiance(l=self.l,func=Ltype,T=Tstar)
        self.rad.setSurface(rstar)
        


    def setProfiles(self,opac=None,mdot=None,mdot_dust=None,T=None,Td=None,
                    mu=None):
    
        '''
        Initialise the independent profiles:opacity, mass-loss rate, 
        initial temperature, dust temperature 
        
        The grids can be given as a 1- or 2-item list, which are used as input 
        for the respective Profiler classes. If not given, they are taken from 
        the inputfile.
        
        Format of the input lists: 
            - [constant value]
            - [func/str,{kwargs}]

        @keyword opac: The opacity profiles (l, cm^2/g). [func,{pars}]
                    
                       (default: None)
        @type opac: [func,dict]
        @keyword mdot: The mass-loss rate profile (r, Msun/yr) [func,{pars}] or 
                       constant.
                    
                       (default: None)
        @type mdot: [func,dict]
        @keyword mdot_dust: The dust mass-loss rate profile (r, Msun/yr) 
                            [func,{pars}] or constant.
                    
                            (default: None)
        @type mdot_dust: [func,dict]
        @keyword T: The initial temperature profile (r, K). [func,{pars}]
                    The first arg is the function, followed by T0 and r0. Can be
                    string Td as well to set it to the dust temperature.
                    
                    (default: None)
        @type T: [func,dict]
        @keyword Td: The dust temperature profile (r, K). [func,{pars}]
                     The first arg is the function.
                    
                     (default: None)
        @type Td: [func,dict]
        @keyword mu: The mean molecular weight if nonstandard. Otherwise 
                     calculated from fH and fHe. Could be specified in case, 
                     e.g., some sort of dissociation is ongoing in the inner 
                     wind.
                     
                     (default: None)
        @type mu: float
        
        '''

        #-- Only reset if self.T is already yet defined.
        if self.T: self.__reset()
        
        #-- Take the profile definitions from the inputfile if not given here
        if opac is None: opac = self.formatInput(self.pars['opac'])
        if mdot is None: mdot = self.formatInput(self.pars['mdot'])
        if mdot_dust is None: 
            mdot_dust = self.formatInput(self.pars['mdot_dust'])
        if T is None: T = self.formatInput(self.pars['Tinit'])
        if Td is None: Td = self.formatInput(self.pars['Td'])
        if mu is None: mu = self.pars.get('mu',None)
        
        #-- Set the opacity profile: l, func, pars
        self.opac = Opacity.Opacity(self.l,opac[0],**opac[1])
        
        #-- Set the mass-loss-rate profiles. Profiler checks itself for constant
        #   If a constant, then mdot[1] is an empty dict
        self.mdot = Mdot.Mdot(self.r,mdot[0],**mdot[1])
        
        #-- Set the dust temperature profile
        self.Td = Temperature.Temperature(self.r,Td[0],**Td[1])

        #-- Use a power law for T in the inner wind?
        self.inner = T[1]['inner']

        #-- Set the initial temperature profile
        if T[0].lower() == 'td':
            self.T_iter[0] = Temperature.Temperature(self.r,Td[0],\
                                                     inner=self.inner,**Td[1])
        else:
            self.T_iter[0] = Temperature.Temperature(self.r,T[0],**T[1])
        
        #-- Determine the T epsilon going from rstar to r0. This is used for the
        #   inner wind. This can be turned off by setting inner=0 in Tinit.
        #   Reset the inner wind power law for the T profiler as well. 
        self.T0 = self.T_iter[0].T0
        self.r0 = self.T_iter[0].r0
        self.inner_eps = -np.log(self.Tstar/self.T0)/np.log(self.rstar/self.r0)
        self.T_iter[0].setInnerEps(self.inner_eps)
        
        #-- Set the dust mass-loss rate, and add the inner radius if needed.
        if not mdot_dust[1].has_key('r0'): mdot_dust[1]['r0'] = self.r0
        self.mdot_dust = Mdot.Mdot(self.r,mdot_dust[0],**mdot_dust[1])
        
        #-- Set the mean molecular weight
        if mu is None: 
            mu = Density.mean_molecular_weight(self.pars['fH'],self.pars['fHe'])
        self.mu = mu
        
        
        
    def setGamma(self,gamma=None):
    
        '''
        Set the adiabatic coefficient. 
        
        Three options:
            - Constant value
            - Step function at 350K (going 7/5 => 5/3 at lower T, MCP model)
            - T-dependent gamma for H2 gas
            
        @keyword gamma: The adiabatic coefficient profile (K, /). Is either a
                        constant (float), or 'h2' in which case the adiabatic
                        coefficient is read from a file as a function of T, or
                        h2_mcp in which case gamma is a step function: 7./5. for 
                        T>350 K, 5./3. otherwise. 
                        
                        (default: None)
        @type gamma: float/str
        
        '''
        
        #-- When gamma == 'h2' or 'h2_mcp', vary gamma as function of T, 
        #   appropriate for H2. 
        if isinstance(self.pars['gamma'],str): 
            gamma = self.pars['gamma'].lower()
            #-- The MCP implementation
            if gamma == 'h2_mcp':
                #-- At low T the diatomic molecule H2 (gamma~1.4) loses internal 
                #   degrees of freedom, and becomes like a mono-atomic gas 
                #   (gamma~1.67).
                kwargs = {'x': self.T_iter[0].eval(inner_eps=self.inner_eps,\
                                                   warn=not self.inner), 
                          'func': Profiler.step,
                          'ylow': 5./3., 'yhigh': 7./5., 'xstep': 350. }
                gamma = Profiler.Profiler(**kwargs)
            #-- Interpolate the adiabatic coefficient taken from NIST database 
            else: 
                fn = os.path.join(cc.path.aux,'h2_physical_properties.dat')
                kwargs = {'x': self.T_iter[0].eval(inner_eps=self.inner_eps,\
                                                   warn=not self.inner),
                          'func': Profiler.interp_file,
                          'ikwargs': {'ext': 3, 'k': 3}, 
                          'filename': fn, 'ycol': -1}
                gamma = Profiler.Profiler(**kwargs)
        
        #-- Alternatively, gamma is constant, and gives the value of the 
        #   coefficient
        else: 
            kwargs = {'x': self.T_iter[0].eval(inner_eps=self.inner_eps,\
                                               warn=not self.inner), \
                      'func': Profiler.constant, 'c': self.pars['gamma']}
            gamma = Profiler.Profiler(**kwargs)
        
        self.gamma = gamma



    def setVelocity(self,v=None):
    
        '''
        Set the velocity profile, based on the requested function.
        
        Done separately from other profiles, because this can depend on the 
        temperature at the inner radius (independent), and on the adiabatic 
        coefficient (for the sound velocity). 
        
        This is done before the line cooling term is set.
        
        @keyword v: The velocity profile (r, cm/s). [func,{pars}]
                    
                    (default: None)
        @type v: [func,dict] or [cst]

        '''

        #-- Take v from the inputfile if not given
        if v is None: v, vdd = self.formatInput(self.pars['v'])
        else: v, vdd = v[0], v[1]
        
        #-- if sound velocity requested for v0, and not enough pars given, 
        #   set automatically
        if v == 'vbeta' and vdd.get('v0',0) == 'vs' and not vdd.has_key('T'):
            vdd['T'], vdd['mu'] = self.T0, self.mu
            vdd['gamma'] = self.gamma.eval(self.T0)
            
        #-- Set the velocity profile: r, func, dfunc, order interp dfunc, pars
        self.v = Velocity.Velocity(self.r,v,**vdd)



    def setLineCooling(self,m):
    
        '''
        Set the line cooling parameters.    
        
        This includes reading the collision rates and level populations.
        
        This is done per molecule, and upon initialisation of the EnergyBalance.
        
        A distinction is made between the source of the spectroscopy/level pops:
          - GASTRoNOoM: The info is taken from the MlineReader and CollisReader
          - ALI: The info is taken from the LamdaReader, PopReader and .par 
                 output file
        
        @param m: The molecule name from the input molecules list.
        @type m: str
        
        '''
        
        #-- Pops were read when the abundances were set. Read the collision 
        #   rates
        imol = self.molecules.index(m)
        self.collis[m] = self.colread(self.pars['collis'][imol])
    
        #-- ipop remembers the iteration number for which the level 
        #   populations were set, so they can be updated later.
        self.ipop[m] = self.i
        
        #-- If no level populations are given, define some default ones. 
        if not self.pars['pop'][imol]: 
            self.setPopInitial(m)
            
        #-- Pops are code-dependent. Spectroscopy as well. Set that here
        elif self.pars['rtcode'].lower() == 'gastronoom': 
            #-- Note that the abundance is in the mline output for GASTRoNOoM
            #   Hence, the populations have already been read.
            #   Refer the spectroscopy to the mline file for later use
            self.mol[m] = self.pop[m]
        else:
            self.pop[m] = self.popread(self.pars['pop'][imol])
        
            #-- Refer the spectroscopy to the lamda file for later use
            self.mol[m] = self.collis[m]
        
        #-- Set the coll rate and level pop interpolators
        self.pop[m].setInterp(itype='spline',k=3,ext=3)
        self.collis[m].setInterp(itype='spline',k=1,ext=0)

        #-- Calculate the cooling term for the initial temperature
        self.Clc(m,update_pop=0)
        
        #-- Calculate the excitation temperature
        self.setTexc(m)
        


    def setPopInitial(self,m):
    
        '''
        Set an initial set of level populations. 
        
        Based on the Boltzmann distribution for a given excitation temperature.
        The kinetic temperature is taken for now. Should likely be scalable in 
        the future. 
        
        Solves the set of equations:
            - Sum(n_i) = 1
            - n_i = n_0 * (g_l/g_0) exp((E_0 - E_l)/Tkin)        

        @param m: The molecule name from the input molecules list.
        @type m: str
        
        '''
        
        #-- In case of GASTRoNOoM, need the radiat file to be read
        if self.pars['rtcode'] == 'gastronoom':
            imol = self.molecules.index(m)
            fn = self.pars['collis'][imol].replace('collis','radiat')
            ny = max(self.collis[m]['coll_trans']['lup'])
            self.mol[m] = RadiatReader.RadiatReader(fn,ny=ny)
        
        #-- Get T profile and other information
        T = self.T.eval(inner_eps=self.inner_eps,warn=not self.inner)
        g0 = self.mol[m].getLWeight(1)
        E0 = self.mol[m].getLEnergy(1,unit='erg')
        ny = self.mol[m]['pars']['ny']
        g = self.mol[m].getLWeight(range(2,ny+1))
        E = self.mol[m].getLEnergy(range(2,ny+1),unit='erg')        

        #-- Force gups/Eups into a column vector for array multiplication
        g.shape = (g.size,1)
        E.shape = (E.size,1)
        
        #-- The equation that calculates n_0
        Ediff = E0 - E
        nfactor = g/g0*np.exp(np.outer(Ediff,1./(k_b*0.999*T)))
        denominator = 1 + sum(nfactor)
        n0 = 1./denominator
        
        #-- Calculate all the other levels
        n = n0*nfactor
        
        #-- Create a PopReader object with filename "None"
        self.pop[m] = PopReader.PopReader(None)
        self.pop[m].setP(self.r)
        self.pop[m].setNY(ny)
        self.pop[m].setPop(index=1,n=n0)
        for i in range(2,ny+1):
            self.pop[m].setPop(index=i,n=n[i-2,:])
    
    
    
    def setAbundance(self,m):
    
        '''
        Set the abundance profile for a molecule.
        
        Two possibilities: 
            - Level populations are given, which include molecular abundance
              profiles as well, and are RT code dependent
            - No level pops are given, so use the default abundance profile

        The second possibility requires a read method and a filename to be given
        in the molecule definition, e.g.
        molecule=12C16O np.loadtxt fname=waql.par usecols=[1,4] skiprows=9 unpack=1
        Alternatively, a constant value can be given as well (not advised).
        
        @param m: The molecule name from the input molecules list.
        @type m: str
        
        '''
        
        #-- Which molecule?
        imol = self.molecules.index(m)
        
        #-- Check if level pops are given (those files contain an abundance 
        #   profile as well). If not, use the molecule definition. 
        if not self.pars['pop'] or len(self.pars['pop']) <= imol \
                or not self.pars['pop'][imol]: 
            mlst = self.pars['molecule'][imol]
            #-- Run a check if an abundance is actually given. 
            if not mlst[0]: 
                q = "No abundance input given in molecule={} ".format(m)+\
                    "definition. Cannot calculate default level populations "+\
                    "or set photoelectric heating."
                raise IOError(q)
        
            #-- A constant value might be requested, otherwise read the file
            #   If constant value, the kwargs dict will be empty.
            if len(mlst[1]) == 0: 
                p, abun = self.r, mlst[0]
            else:
                dd = DataIO.read(func=mlst[0],**mlst[1])

                #-- Check whether the values are given in rstar or in cm
                if dd[0][0]/self.pars['rstar']<0.1: 
                    p = dd[0]*self.pars['rstar']
                else:
                    p = dd[0]

                #-- In case three columns are read, nmol and nh2 are read 
                #   separately
                if len(dd) == 2:
                    abun = dd[1]
                elif len(dd) == 3: 
                    abun = dd[1]/dd[2]
        
        #-- If pops are given, extract the abundance profile, depending on the
        #   rt code. Molecule information in inputfile is then not used.
        elif self.pars['rtcode'].lower() == 'gastronoom': 
            #-- Abundances for mline in GASTRoNOoM are in the mline output, 
            #   which contains the level populations. So read that here already
            self.pop[m] = self.popread(self.pars['pop'][imol])
            
            #-- Extract the abundance profile as a function of impact parameter
            p = self.pop[m].getP()
            abun = self.pop[m].getProp('amol')
        else:
            #-- Abundance profile for MCP/ALI is in a separate parameter output
            #   file
            fn = self.pars['pop'][imol].replace('pop','log')                
            
            #-- Find out number of cells. Note vols 0 and 1 are not used, so -2.
            nvol = DataIO.getKeyData(filename=fn,incr=0,keyword='nvol')[0]
            nvol = int([line.split('=')[1] 
                        for line in nvol 
                        if 'nvol' in line.lower()][0])
            
            #-- Find out where the data start + 2 lines for units and separator
            data = DataIO.readFile(fn,' ')
            ixmol = DataIO.findKey(i=0,key='xmol',data=data) + 2
            p, abun = np.genfromtxt(fn,usecols=[1,4],skip_header=ixmol+1,\
                                    unpack=1,max_rows=nvol-2)
        
        #-- Constant value is requested
        if len(p) != len(abun):
            abun_interp = mlst[0]

        #-- We need flexible extrapolation, use interp1d and linear interp
        #   Linear because cubic is too unstable for the last radial step in
        #   a steep abundance decline. 
        else:
            abun_interp = interp1d(x=p,y=abun,kind='linear',bounds_error=0,\
                                   assume_sorted=1,fill_value=(abun[0],0.))
        self.abun[m] = Profiler.Profiler(x=self.r,func=abun_interp)
        
        
        
    def setTexc(self,m):
    
        '''
        Calculate the excitation temperature for the level populations given in 
        self.pop. This method is called when the pops are read for the first 
        time (or set by setInitialPop), and when the level populations are 
        updated.
        
        The excitation temperature is set as an array of dimensions 
        (number of transition indices, number of impact parameter), hence for a
        given transition index, you can retrieve the excitation temperature as
        self.Texc[m][i-1,:] as a function of impact parameter. Alternatively, 
        you can access Texc by calling getTexc and passing either the indices or
        the lup/llow.
        
        The impact parameter grid can be retrieved from self.pop[m].getP().
        
        Note that the excitation temperature is calculated for all radiative
        transitions included in the molecular spectroscopy; not for all 
        collisional transitions. 
        
        @param m: The molecule tag
        @type m: str
        
        '''
        
        #-- Retrieve the transition indices, and then extract relevant 
        #   information for all transitions.
        indices = self.mol[m].getTI(itype='trans')
        lup = self.mol[m].getTUpper(indices)
        llow = self.mol[m].getTLower(indices)
        gl = self.mol[m].getLWeight(llow)
        gu = self.mol[m].getLWeight(lup)
        El = self.mol[m].getLEnergy(llow,unit='erg')
        Eu = self.mol[m].getLEnergy(lup,unit='erg')
        
        #-- Retrieve the associated populations
        nl = self.pop[m].getPop(llow)
        nu = self.pop[m].getPop(lup)
        
        #-- Force g/E into column vectors for array multiplication
        gl.shape = (gl.size,1)
        gu.shape = (gu.size,1)
        El.shape = (El.size,1)
        Eu.shape = (Eu.size,1)
        
        #-- The equation that calculates Texc
        self.Texc[m] = -1./np.log(nu/nl*gl/gu)*(Eu-El)/k_b
        


    def getTexc(self,m,indices=None,lup=None,llow=None):
    
        """
        Return the excitation temperature profile as a function of impact 
        parameter for given transition indices of a given molecule. 
        
        Can also take upper and/or lower level indices to determine the 
        transition indices. 

        The excitation temperature is returned as an array of dimensions 
        (number of transition indices, number of impact parameters).

        @param m: The molecule tag
        @type m: str
        
        @keyword indices: The transition indices. Can be indirectly given 
                          through upper and/or lower level indices. If default
                          and no llow/lup are given, all Texc are returned for
                          this molecule.
                          
                          (default: None)
        @type indices: list/array
        @keyword lup: The upper level indices used to extract the transition
                      indices. If this or llow are defined, the keyword indices 
                      is overwritten.
                      
                      (default: None)
        @type lup: list/array
        @keyword llow: The upper level indices used to extract the transition
                       indices. If this or lup are defined, the keyword indices 
                       is overwritten.
                      
                       (default: None)
        @type llow: list/array
        
        @return: The excitation temperature, with array shape (n_index,n_p)
        @rtype: array
        
        """
        
        #-- If either lup or llow are given, get transition indices from them
        if not (lup is None and llow is None):
            if not lup is None: lup = Data.arrayify(lup)
            if not llow is None: llow = Data.arrayify(llow)
            indices = self.mol[m].getTI(lup=lup,llow=llow)
        
        #-- If indices is still None, all transitions are returned
        if indices is None: return self.Texc[m]
        
        #-- Make sure the indices are given as an array, and select. Note 
        #   python's 0-based indexing vs the fortran 1-based indexing.
        indices = Data.arrayify(indices)
        return self.Texc[m][indices-1,:]
    
    
    
    def __reset(self):
    
        '''
        Reset all variables that depend on the radius, velocity, etc.
        
        '''
        
        #-- Density profiles
        self.gdens = None
        self.ddens = None
        self.nd = None
        
        #-- Dust velocity profile
        self.vd = None
        self.w = None
        
        #-- Temperature profile and iterations
        self.T = None
        self.T_iter = dict()
        self.i = 0
        
        #-- And the heating and cooling terms
        self.H = {}
        self.C = {}
        self.rates = {}
    
    
    
    def __next_iter(self):
    
        '''
        Reset some of the properties that depend on temperature so they can be 
        calculated anew.
        
        '''
        
        #-- Dust and drift velocity profiles depend on T if w_thermal != none
        self.vd = None
        self.w = None
        
        #-- Set the current temperature profile to the new calculation
        self.__setT()
               
        #-- Recalculate the dust velocity. Drift will be done as well.
        self.setVdust() 
    
    
    
    def setDrift(self):
    
        '''
        Calculate the drift velocity profile. 
        
        This assumes a stellar radiance profile is given, of which the 
        wavelength grid is used to integrate Q_l * L_l, thereby interpolating
        the opacity. 
        
        The two available modes are standard and beta: 
          - standard: Calculates the drift velocity based on the balance between
          the radiation pressure and the drag force. 
          - vbeta: Follows MCP, where the terminal drift velocity is calculated 
          from the terminal gas velocity, based on the balance between radiation
          pressure and drag force, with a beta law going to that max velocity.
        
        Note that the vbeta mode follows the F mode for dust velocity in MCP. 
        
        
        '''

        if self.w is None:
            #-- Set density for the mean molecular weight.
            self.setDensity('gas')

            #-- Define shared input parameters for driftRPDF function
            kwargs = {'a': self.a, 'l': self.l, 'P': self.pars['P'], 
                      'sd': self.pars['sd'], 'opac': self.opac, 
                      'radiance': self.rad, 'alpha': self.pars['alpha'],
                      'mu': self.gdens.getMeanMolecularWeight()}

            #-- Check if the MCP method is requested, otherwise do 'standard'
            method = self.formatInput(self.pars['w_mode'])
            if method[0] == 'vbeta2D':
                #-- Calculates the terminal drift velocity from the terminal gas
                #   velocity, and sets a beta law for the inner wind up to it
                #-- For this the terminal gas velocity is needed, and the mdot
                #   at the inner radius (no variable mass loss or thermal comp)
                kwargs['v'] = max(self.v.eval())
                kwargs['mdot'] = self.mdot.eval(self.r0)
                kwargs['w_thermal'] = 'none'
                wmax = Velocity.driftRPDF(r=self.r0,**kwargs)[0]

                #-- Set the drift as a beta law with sensible values, unless 
                #   they were defined in the w_mode key. Even the max velocity
                #   can be fixed in the input through vinf.
                v0 = method[1].get('v0',self.v.eval(self.r0))
                r0 = method[1].get('r0',self.r0)
                beta = method[1].get('beta',1.0)
                wmax = method[1].get('vinf',wmax)
                self.w = Velocity.Drift(r=self.r,a=self.a,\
                                        func=Velocity.vbeta2D,\
                                        r0=r0,v0=v0,vinf=wmax,beta=beta)
            
            else:
                #-- If standard, add the relevant driftRPDF keywords and set the
                #   profile
                kwargs['w_thermal'] = self.pars['w_thermal']
                kwargs['v'] = self.v
                kwargs['r'] = self.r
                kwargs['mdot'] = self.mdot 
                kwargs['func'] = Velocity.driftRPDF
                
                #-- Include T if a thermal velocity term is needed
                vtherm_types = ['kwok','mean','rms','prob','epstein']
                if self.pars['w_thermal'].lower() in vtherm_types:
                    #-- T is always set upon initialisation
                    kwargs['T'] = self.T
                
                self.w = Velocity.Drift(**kwargs)
        
        
    
    def setVdust(self):
    
        '''
        Set the dust velocity profile. The drift is averaged over grain size.
        
        '''
        
        if self.vd is None:
            self.setDrift()
            kwargs = {'r': self.r, 'func': Velocity.vdust,
                      'w': self.w, 'v': self.v, 'norm_type': 'standard'}
            self.vd = Velocity.Velocity(**kwargs)
        
        
        
    def setDensity(self,dtype): 
    
        '''
        Calculate the density profile for gas or dust. 
        
        The dust density profile is not dependent on grain size, since the drift
        averaged over grain size.
        
        @param dtype: The density profile type, being 'gas' or 'dust'.
        @type dtype: str
        
        '''
        
        #-- Set the density profile for either dust or gas. Nothing is done if 
        #   dtype is not recognized.
        if dtype.lower() == 'gas':
            if self.gdens is None:
                kwargs = {'r': self.r,'func':Density.dens_mdot,
                          'mdot': self.mdot, 'v': self.v}
                self.gdens = Density.Density(**kwargs)
                kwargs = {'fH': self.pars['fH'], 'fHe': self.pars['fHe']}
                self.gdens.calcNumberDensity(**kwargs)

                #-- Replace mean molecular weight in case a nonstandard value is
                #   requested
                self.gdens.mu = self.mu
        
        if dtype.lower() == 'dust':
            #-- Calculate the mass density
            if self.ddens is None:
                #-- Calculate the dust velocity profile. The drift is averaged
                #   over the grain size.
                self.setVdust()
                kwargs = {'r': self.r,'func':Density.dens_mdot, 'dust': 1,
                          'mdot': self.mdot_dust, 'v': self.vd}
                self.ddens = Density.Density(**kwargs)

            #-- Calculate the number density
            if self.nd is None and self.a.size > 1:
                #-- Set the grain size distribution. self.nd refers to a
                #   Grainsize.Distribution object!
                kwargs = {'r': self.r, 'a': self.a, 'sd': self.pars['sd'],
                          'func': getattr(Grainsize,self.pars['GSD']), 
                          'w': self.w,
                          'v': self.v, 'mdot_dust': self.mdot_dust}
                self.nd = Grainsize.Distribution(**kwargs)
                
            elif self.nd is None:
                #-- Working with average grain size. Calculate number density
                #   from the mass density: rho/(4/3 pi rho_specific a^3)
                self.ddens.calcNumberDensity(sd=self.pars['sd'],a=self.a)
                #-- In this case, self.nd just refers to the Density() object
                self.nd = self.ddens
        
    
    
    def updatePop(self,m,fn='',updateLC=0):
    
        '''
        Update the level populations.
        
        This reads the level populations from the same file as listed upon
        creation of the EnergyBalance and assumes those level populations are
        appropriate for the temperature profile of the iteration number self.i.
        
        This is done per molecule.
        
        @param m: The molecule name from the input molecules list.
        @type m: str
        
        @keyword fn: The filename of a populations file, in case it changes or 
                     is added anew after default populations have been used. 
                     Ignored if the file does not exist. 
                     
                     (default: '')
        @type fn: str
        @keyword updateLC: The flag that moves the self.ipop to the current (or
                           next) iteration, so that calcClc knows a new line 
                           cooling term must be updated. This is always set to 1
                           if the populations are updated. Setting this to 1 
                           upon function call, tells calcClc the LC term must be
                           calculated anew regardless of reading new pops (and 
                           thus uses the old pops with the new T profile of this
                           iteration). This key can be set in the inputfile.
                           
                           (default: 0)
        @type updateLC: bool
        
        '''
        
        #-- First check if populations are given. If not, don't update them.
        imol = self.molecules.index(m)
        
        #-- Maybe a new file is added for pops
        if os.path.isfile(str(fn)): 
            self.pars['pop'][imol] = fn
        
        #-- If a file name is available, check if updating is needed
        if self.pars['pop'][imol]:
            #-- Select the files and read them with the Reader objects.
            npop = self.popread(self.pars['pop'][imol])
            if False in  [np.array_equal(npop.getPop(i),self.pop[m].getPop(i)) 
                          for i in npop.getLI()]:
                self.pop[m] = npop

                #-- Set the level pop interpolators and calc Texc
                self.setTexc(m)
                self.pop[m].setInterp(itype='spline',k=3,ext=3)
                
                #-- Make sure to tell calcClc to recalculate the LC term
                updateLC = 1
                
        #-- Tell calcClc to recalculate the cooling term for this molecule
        if updateLC or self.pars.get('updateLC',0):
            #-- If there's a key for this iteration index already, the cooling
            #   term was already calculated, so set the index for the update to
            #   the next iteration.
            if self.C['lc_{}'.format(m)].has_key(self.i):
                self.ipop[m] = self.i + 1
            else:
                self.ipop[m] = self.i


            
    def __setT(self):
    
        '''
        Set the current temperature profile for this iteration. When i == 0,
        this is the initial guess.
        
        This is done explicitly, because upon initialisation the state of the
        object depends on self.T being None or not.
        
        '''
        
        self.T = self.T_iter[self.i]
        
        
    
    def iterT(self,conv=0.01,imax=50,step_size=0.,dTmax=0.10,warn=1,\
              *args,**kwargs):
    
        '''
        Iterate the temperature profile until convergence criterion is reached.
        
        Extra arguments are passed on to the spline1d interpolation of the 
        total cooling and heating terms through the calcT() call.
        
        @keyword conv: The maximum relative allowed change in T for convergence.
        
                       (default: 0.01)
        @type conv: float
        @keyword imax: Maximum number of allowed iterations. Code stops when 
                       this number is reached, but not before dTmax = 1.
        
                       (default: 50)
        @type imax: int
        @keyword step_size: The minimum step size in maximum allowed relative
                            temperature change. The code decides dynamically 
                            based on how close the iteration is to convergence
                            which multiple of this step size the next iteration
                            allows. If T change percentage grows to fast, 
                            decrease this number. Default is 0, in which case 
                            the maximum allowed T change is kept constant at
                            dTmax.
                            
                            (default: 0)
        @type step_size: float
        @keyword dTmax: The starting value for the maximum allowed relative
                        temperature change between iterations. This maximum 
                        increases as the code reaches convergences, through a 
                        set multiple of step_size. If step_size is 0, dTmax 
                        stays constant
                            
                        (default: 0.10)
        @type dTmax: float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        
        '''
        
        print('-- Iterating T(r) now.')
        dTnsteps = (1.-dTmax)/step_size if step_size else 0.
        steps = 0.
        
        #-- If first iteration in this call of iterT, calculate it in any case
        #   This is not necessarily iteration 0! The pops may have been updated.
        print('Iteration {} for T(r)...'.format(self.i+1))
        self.calcT(dTmax,warn=warn,*args,**kwargs)
        
        #-- Then continue iterating until the relative difference between this 
        #   and the previous result is smaller than convergence criterion in all
        #   radial points, or when the maximum number of iterations is reached
        while True:
            dTdiff = np.abs(1.-self.T.eval(warn=not self.inner)\
                          /self.T_iter[self.i-1].eval(warn=not self.inner))
            if (np.all(dTdiff<conv) and steps == dTnsteps) or self.i == imax: 
                break
            if np.all(dTdiff<(dTnsteps-steps)*conv):
                if dTnsteps-steps > 8.: steps += 4.
                elif dTnsteps-steps > 4.: steps += 2.
                else: steps += 1.
            print('Iteration {} for T(r)...'.format(self.i+1))
            self.calcT(dTmax+step_size*steps,warn=warn,*args,**kwargs)
            if not np.all(np.isfinite(self.T.eval(inner_eps=self.inner_eps,\
                                                  warn=not self.inner))):
                print('NaNs found in T-profile. Breaking off iteration.')
                break
    
    
    
    def calcT(self,dTmax=1,warn=1,ode_kwargs={},*args,**kwargs):
    
        '''
        Calculate the temperature profile based on the differential equation 
        given by Goldreich & Scoville (1976). The function is defined in dTdr.
        
        The initial temperature is taken to be the second argument of the 
        initial temperature profile given by T in inputEnergyBalance.dat.  
        This is typically the condensation temperature.
        
        Additionals arguments are passed on to the spline1d interpolation of the
        total cooling and heating terms, e.g. k=3, ext=0 are defaults.
        
        @keyword dTmax: The maximum allowed relative temperature change for this
                        T calculation. Set to 100% by default.
                            
                        (default: 1)
        @type dTmax: float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        @keyword ode_kwargs: Extra arguments for the ODE solver. They are added
                             to the dict ode_args made by the function, and 
                             hence overwrites any defaults. In principle, the 
                             defaults are fine, but can be overwritten if needed
        
                             (default: {})
        @type ode_kwargs: dict
        
        '''
        
        #-- Note that T is always set upon initialisation and at the end of a 
        #   new T calculation.
        
        #-- Move to the next iteration. Cooling and heating terms belong to this
        #   next iteration (even tho they are based on initial guess T profile,
        #   except for the adiabatic cooling term)
        self.i += 1
        
        #-- Calculate the gas density if it hasn't been calculated yet
        self.setDensity('gas')
        
        #-- Calculate the heating and cooling terms
        #   Note that the adiabatic term is used explicitly in the dTdr 
        #   function. The rate is calculated with the new T profile at the end
        #   for plotting.
        for term in self.pars['hterms']:
            getattr(self,'H'+term)()
        for term in self.pars['cterms']:
            getattr(self,'C'+term)()            
        
        #-- The abundance term can be included here as well to reduce the number
        #   of profile evaluations in dTdr.
        nh2 = self.gdens.eval(dtype='nh2')
        
        #-- The ascale factor is the correction for non-zero He and atomic H 
        #   abundances (Decin et al. 2006). Reduces to original if fH = fHe = 0.
        ascale = self.gdens.fH + 1. + self.gdens.fHe * (self.gdens.fH + 2.)
        
        #-- The density term does not yet include the velocity factor. This is
        #   included when dTdr is calculated.
        dens_term = nh2 * k_b * ascale

        #-- Then set the heating and cooling terms as profiler objects
        yrates = sum([self.H[k][self.i] for k in self.H.keys()]) \
               - sum([self.C[k][self.i] for k in self.C.keys() if k != 'ad'])

        #-- Remember the rates profiler for each iteration
        self.rates[self.i] = yrates/dens_term
        
        #-- Note that if no rates are requested, this will give an empty array
        if self.r.size == Data.arrayify(yrates).size:
            rp = Profiler.Profiler(self.r,spline1d(self.r,self.rates[self.i],\
                                                   *args,**kwargs))
        else:
            rp = None 
        
        #-- Calculate the next iteration of the temperature profile
        ode_args = {'func': dTdr, 'y0': self.T0, 't': self.r,
                    'args': (self.v,self.gamma,rp,warn)}
        ode_args.update(ode_kwargs)
        Tr = odeint(**ode_args)[:,0]
        
        #-- Check the temperature change. Limit it to the requested maximum 
        #   allowed change. Only do this from the second iteration, to allow 
        #   a big jump from the initial condition.
        if self.i < 0: 
            dTmax = self.pars['dTmax']
        print('Changing T by {:.1f}%.'.format(dTmax*100.))
        Ti = self.T.eval(inner_eps=self.inner_eps,warn=not self.inner)        
        Tr = np.where(Tr<=0.,np.ones_like(Tr),Tr)
        Tr = Ti + dTmax*(Tr-Ti)
        
        #-- Create a profiler for the new temperature structure. Extrapolation
        #   done by returning boundary values, ie T0 and the temp at outer 
        #   boundary
        Tinterp = spline1d(self.r,Tr,k=3,ext=3)
        keys = {'inner':self.inner,'inner_eps':self.inner_eps,'r0':self.r0,
                'T0':self.T0}
        self.T_iter[self.i] = Temperature.Temperature(self.r,Tinterp,**keys)
        
        #-- Initialize some T-dependent properties for the next iteration. Also
        #   sets the new temperature profile as the current one.
        self.__next_iter()
    
        #-- Calculate the adiabatic term for plotting. Also sets the current
        #   temperature profile as the new calculation.
        self.Cad()

    
    
    def Cad(self):
        
        '''
        Calculate the adiabatic cooling rate. This is not used by the solver of 
        the dTdr differential equation, and is only for plotting purposes. 
        
        Uses the latest calculation of the T-profile. 
        
        Equation derived from Decin et al. 2006 and Danilovich 2016.
        
        ''' 

        #-- if cooling rate has already been calculated: don't do anything
        if self.C['ad'].has_key(self.i): return 
        
        #-- Set the gamma-based factor evaluated at the new T profile
        #   calculated.
        gamma = self.gamma.eval(self.T.eval(inner_eps=self.inner_eps,\
                                            warn=not self.inner),warn=0)
        
        #-- Makes the rate negative, in accordance with -C in the second term of
        #   the diff eq.
        factor1 = 2-2*gamma
        factor2 = gamma-1
        
        #-- Determine the velocity and its derivative from the Velocity() object
        #   Note that these will never complain bout extrapolation unless they
        #   are profiles read beforehand instead of functions such as vbeta
        vi = self.v.eval()
        dvi = self.v.diff()
        
        #-- Calculate the conversion term K/cm => erg/s/cm3: allows to compare 
        #   C terms with adiabatic term.
        self.setDensity('gas')
        nh2 = self.gdens.eval(dtype='nh2')
        ascale = self.gdens.fH + 1. + self.gdens.fHe * (self.gdens.fH + 2.)
        conversion = nh2 * k_b * vi * ascale / factor2
        
        #-- Calculate the temperature dependent term
        main_term = (1./self.r + 0.5*1./vi*dvi) \
                        * self.T.eval(inner_eps=self.inner_eps,\
                                      warn=not self.inner)
        
        #-- Calculate the rate. Multiply by -1 to have positive cooling rate.
        #   [erg/s/cm3]
        self.C['ad'][self.i] = -1. * conversion * factor1 * main_term
    
    
    
    def Hdg(self):
    
        '''
        Calculate the heating rate by dust-gas collisions. 
        
        Note that this term evaluates to zero for the inner wind at r<r0, if the
        dust density is 0 there. This can be enforced by choosing an mdot_step
        function.
        
        Sputtering is applied if a grain size distribution is used. For a 
        constant grain size, this is not done, since it's not realistic to 
        remove all dust if the drift becomes too large.
        
        The equation is derived from Goldreich & Scoville 1976, based on Schoier
        et al 2001, and Decin et al. 2006: 
        
        General case: 
        Hdg = n_d sigma_d w  x  1/2 rho_g w^2
        =>
        Hdg = pi/2 n_H2 m_H (fH+2) (1+4fHe) (1-P)^(2/3) Int(n_d w^3 a^2 da)
        
        MCP/ALI case (average grain size a): 
        Hdg = pi n_H2 m_H n_d w^3 a^2
        
        GASTRoNOoM case (MRN): 
        Hdg = pi/2 A n_H2^2 m_H (fH+2)^2 (1+4fHe) Int(w^3 a^-1.5 da)
        
        '''
        
        #-- if heating rate has already been calculated: don't do anything
        if self.H['dg'].has_key(self.i): return 

        #-- Initialize the current temperature profile as the initial guess
        #   This is not done if i > 0, since the last solution is set as the 
        #   current T profile after a calculation
        if not self.T: 
            self.__setT()

        #-- Initialise density and a few parameters
        self.setDrift()
        self.setDensity('gas')
        self.setDensity('dust')
        nh2 = self.gdens.eval(dtype='nh2')
        fH = self.pars['fH']
        fHe = self.pars['fHe']
        P = self.pars['P']
        mode = self.pars['heatmode'].lower()
        
        #-- Set the constant factor. Same for both modes
        cfac = np.pi*nh2*mh*(fH+2.)*(1.+4.*fHe)*(1.-P)**(-2./3.)
        
        #-- Hdg according to Gail & Sedlmayr
        if mode == 'gs2014':
            #-- The thermal velocity and factor 4/3. for specular reflection: 
            vT = 4./3.*np.sqrt(8.*k_b*self.T.eval(inner_eps=self.inner_eps,\
                                                  warn=not self.inner)\
                             /(self.gdens.mu*np.pi*mh))
            
        #-- Hdg according to Decin et al. 2006: extra factor 0.5 for Ekin, no
        #   thermal component.
        else:
            cfac *= 0.5
            vT = 0.

        #-- Differ between having a grain size distribution or not. Note that 
        #   this factor remains constant for the same gas density, regardless of
        #   fH, fHe. It is a different matter for gsfac where drift enters.
        if self.a.size > 1: 
            #-- GSD: integrate drift, dust number density, a^2 over grain size
            #   self.nd is a Grainsize.Distribution object! Apply sputtering.
            vT = np.transpose([vT])
            nd = self.nd.eval(w=self.w.eval(),w_sputter=self.pars['w_sputter'])
            w = self.w.eval()
            afac = nd*self.a**2.*w**2.*(vT**2.+w**2.)**0.5
            gsfac = trapz(x=self.a,y=afac,axis=1)
        else: 
            #-- Avg grain size: no integration
            #   self.nd refers to the self.ddust, or a dust Density() object!
            #   no sputtering for a single grain size
            nd = self.nd.eval(dtype='ndust')
            w = np.squeeze(self.w.eval())
            gsfac = nd*self.a**2*w**2.*(vT**2.+w**2.)**0.5

        #-- Set dust-gas collisional heating [erg/s/cm3]
        self.H['dg'][self.i] = cfac*gsfac
        


    def Hdt(self):
    
        '''
        Calculate the heating rate by dust-gas heat exchange. 
        
        The equation is derived from Burke & Hollenbach 1983, based on 
        Groenewegen 1994, Schoier et al. 2001, and Decin et al. 2006: 
        
        Note that this term evaluates to zero for the inner wind at r<r0, if the
        dust density is 0 there. This can be enforced by choosing an mdot_step
        function.
        
        Sputtering is applied if a grain size distribution is used. For a 
        constant grain size, this is not done, since it's not realistic to 
        remove all dust if the drift becomes too large.
        
        The older implementations used by MCP/ALI and GASTRoNOoM follow Burke & 
        Hollenbach 1983 and Groenwegen 1994. The general-case implementation 
        follows the more recent work of Gail & Sedlmayr, adding in proper vT and
        drift terms. The GS2014 is not yet fully implemented, but the Hdt term 
        is already available here.
        
        General case (following Gail & Sedlmayr 2014, see Eq 15.19): 
        Hdt = alpha pi k_b n_h2 (fH+2.)(1.-P)^(-2./3.) Int(n_d a^2 da) (T_d - T)
        sqrt(8*vT^2+drift^2) (1/(gamma-1))
        
        MCP/ALI case (n_d: dust number density, average grain size a): 
        Hdt = 2 alpha pi k_b n_h2 (n_d a^2) (T_d - T) vT
        
        GASTRoNOoM case (n_d: distribution following MRN): 
        Hdt = 2 alpha pi k_b n_h2 (fH+2.) Int(n_d a^2 da) (T_d - T) vT

        In all cases, alpha is the accommodation coefficient (from Groenewegen
        1994): 
        alpha = 0.35 exp(-sqrt(0.002*(T_d + T)))+0.1
        
        and vT is the thermal velocity:
        vT = sqrt(8 k_b t/(mu pi m_H))
        
        '''
        
        #-- if heating rate has already been calculated: don't do anything
        if self.H['dt'].has_key(self.i): return 
        
        #-- Initialise density and a few parameters
        self.setDrift()
        self.setDensity('gas')
        self.setDensity('dust')
        nh2 = self.gdens.eval(dtype='nh2')
        fH = self.pars['fH']
        P = self.pars['P']
        mode = self.pars['heatmode']
        
        #-- Extract the dust and gas temperatures
        td = self.Td.eval()
        t = self.T.eval(inner_eps=self.inner_eps,warn=not self.inner)
        
        #-- Calculate the average accommodation coefficient
        accom = 0.35*np.exp(-np.sqrt((td+t)*0.002))+0.1
        
        #-- Calculate the constant factor as a function of r 
        cfac = np.pi*k_b*nh2*(fH+2.)*(1.-P)**(-2./3.)
        
        #-- The thermal velocity: 
        vT = np.sqrt(8*k_b*t/(self.gdens.mu*np.pi*mh))
        
        #-- Calculate the temp factor
        tdiff = (td-t)
        
        #-- Differ between having a grain size distribution or not.
        if self.a.size > 1: 
            #-- GSD: integrate dust number density, a^2 over grain size
            #   self.nd is a Grainsize.Distribution object! Apply sputtering.
            nd = self.nd.eval(w=self.w.eval(),w_sputter=self.pars['w_sputter'])
            gsfac = trapz(x=self.a,y=nd*self.a**2,axis=1)
        else: 
            #-- Avg grain size: no integration
            #   self.nd refers to the self.ddust, or a dust Density() object!
            #   A single grain size: no sputtering.
            nd = self.nd.eval(dtype='ndust')
            gsfac = nd*self.a**2

        #-- Set the dust-gas collisional heating term. First mode is Gail & 
        #   Sedlmayr's equation from 2014, based on Draine 1980. Default is the
        #   accommodation heat exchange by Burke & Hollenbach 1983
        if mode.lower() == 'gs2014':
            #-- Avg over a^2. Normalisation disappears when combined with gsfac
            #   Sputtering is not included for the average of w, since it would
            #   skew the distribution of drift velocities to a lower value. We 
            #   cannot approximate what happens to destroyed grains at this 
            #   point, hence we don't know how this would affect the drift 
            drift = self.w.avgDrift(norm_type='collisional',nd=self.nd)
            
            #-- fg/2. based on gamma, with fg = 2/(gamma-1)
            fg = 1./(self.gamma.eval(self.T.eval(inner_eps=self.inner_eps,\
                                                 warn=not self.inner),\
                                     warn=0)-1.)
            
            #-- Velocity factor 
            vfac = (8.*vT**2+drift**2)**0.5
            self.H['dt'][self.i] = accom * cfac * fg * tdiff * gsfac * vfac
            
        else: 
            #-- Extra fac 2 enters through 2kT (BH1983)
            self.H['dt'][self.i] = 2. * accom * cfac * tdiff * vT * gsfac
            


    def Hpe(self):
    
        '''
        Calculate the heating rate by the photoelectric effect. 
        
        Two methods are available at present: 
        1) Following Draine 1978 and Huggins et al. 1988, as used by MCP
        (see also Crosas & Menten 1997)
        2) Following Bakes & Tielens 1994., as implemented by Decin et al.
        2006 in GASTRoNOoM. 
               
        No options to tweak these methods has been implemented yet, but a lot of
        consistency checks should be done, and some of the assumptions can be 
        improved on with a consistent calculation (e.g. Av in method 2).
        
        Photoelectric heating as derived by Bakes & Tielens is maybe a good 
        basis for further development of the heating term, but see also Woitke
        2015 (2015EPJWC.10200011W) for recent developments.
        
        The derivation in method 2 makes a lot of assumptions and are applied
        specifically to the case of GASTRoNOoM (e.g. fixed grain size 
        distribution). Hence, the implementation is considered appropriate for 
        the default settings in the inputEnergyBalance_gastronoom.dat file. Any
        deviations from that must be treated carefully.
        
        No sputtering is applied. Small grains are the most relevant for Hpe, 
        and any sputtering would reduce the size of existing dust grains, thus
        increasing the amount of small grains. Since we cannot have Hpe directly
        depend on the grain size distribution (and instead work with a scaling 
        factor amin_scale), we do not apply sputtering to the photoelectric 
        heating.
        
        '''
        
        #-- if heating rate has already been calculated: don't do anything
        if self.H['pe'].has_key(self.i): return 
        
        #-- Initialise density
        self.setDensity('gas')
        self.setDensity('dust')
        
        #-- Method 1: Scale the heating rate with the gas density and a constant
        #   factor (Draine 1978) -- Maybe ngas is better here
        method = self.formatInput(self.pars['pe_method'])
        if method[0].lower() == 'draine':
            nh2 = self.gdens.eval(dtype='nh2')
            #-- Default Kpe 1e-26 erg/s, following Huggins et al. 1988 and 
            #   Schoier & Olofsson 2001 
            Kpe = method[1].get('Kpe',1e-26)
            self.H['pe'][self.i] = Kpe * nh2
        
        #-- Method 2: More elaborate calculation of the heating rate, though
        #   caution is advised due to the uncertainties in quite a few of the
        #   assumptions.
        elif method[0].lower() == 'bakes':
            T = self.T.eval(inner_eps=self.inner_eps,warn=not self.inner)
            G0 = method[1].get('G0',1.)
            amin_scale = method[1].get('amin_scale',0.2)
            d2g = self.ddens.eval()/self.gdens.eval()
            nhtot = self.gdens.eval(dtype='nhtot')

            #-- Calculate H_pe-hat from Decin 2006 EQ. 9 (see Sect 2.2.2)
            #   First: electron density. Assumed to come from photodissociation
            #   of CO into C and O, whereby C is almost immediately ionised. 
            #   Hence, difference between C/O > 1 and C/O < 1, since the former 
            #   has free C available for ionisation. 
            #   Note that we want the relative CO abundance compared to the max
            #   at R_STAR or R_INNER.
            aco = self.abun['12C16O'].eval()
            xco = aco/aco[0]
            abun_c = method[1].get('abun_c')
            abun_o = method[1].get('abun_o')
            if abun_c > abun_o:
                ne = nhtot * abun_c * (1-xco*abun_o/abun_c)

            #-- If C/O < 1, all C is in CO: ne = 0, in the inner wind. More e in
            #   outer wind where CO gets photodissociated.  
            #   Note that a small bump is formed possibly in the inner wind, but
            #   this is due to the GASTRoNOoM interpolation, and having values 
            #   very close to zero. It's not important
            #   considering how small the rate is in the inner wind.
            else:
                ne = nhtot * abun_c * (1-xco)
            
            #-- Calculate the factor itself, noting the a_min scaling and G0
            Hpehat = amin_scale*1e-24*nhtot
            Hpehat[ne>0] *= 0.03/(1+2e-4*G0*np.sqrt(T[ne>0])/ne[ne>0])
            Hpehat[ne==0] = 0.

            #-- Calculate the correction on H_pe-hat based on the optical depth
            #   Note the use of G0 again. tau_uv = 1.8 * A_v = 1.8 * 1.6e-22 N_H
            #   H Column density from R_outer integrated inwards, as a function
            #   of self.r. Note r_outer-r to have a positive column density.
            Nh = cumtrapz(nhtot[::-1],x=(self.r[-1]-self.r)[::-1],\
                              initial=0.)[::-1]

            #-- Extra scaling factor for the dust-to-gas ratio in the extinction
            Av =  1.6e-22*Nh/0.01*d2g
            self.H['pe'][self.i] = Hpehat*G0*np.exp(-1.8*Av)



    def Hcr(self):
    
        '''
        Calculate the heating rate by cosmic rays. 
        
        The rate is taken from Goldsmith & Langer 1978, and is also used by 
        Groenwegen 1994 and Decin et al 2006. This is a very approximate formula
        and will need updating in the future. 
        
        Two modes are available: The standard one, using n(H2), and the one used
        by Groenewegen 1994 and Decin et al. 2006 (where all gas mass is placed 
        into H2, even if fH or fHe are non-zero). Controlled through the keyword
        cr_method == 'groenewegen' or cr_method == 'standard' in the inputfile.
        
        '''
        
        #-- if heating rate has already been calculated: don't do anything
        if self.H['cr'].has_key(self.i): return 
        
        #-- Initialise density
        self.setDensity('gas')
        
        #-- Goldsmith & Langer 1978 give the rate for ionisation of H2 by cosmic
        #   rays. However Groenewegen 1994 and Decin 2006 use this rate by 
        #   placing all gas mass into H2 (even if fH and fHe are not zero), ie
        #   by dividing the density rho by 2*m_H
        #   Technically this should only be nh2, but we follow Groenewegen/Decin
        #   for now. Why?! Both options are implemented.
        nh2 = self.gdens.eval(dtype='nh2')
        fH = self.gdens.fH
        fHe = self.gdens.fHe
        
        #-- Decin et al 2006/Groenewegen 1994
        if self.pars['cr_method'].lower() == 'groenewegen':
            self.H['cr'][self.i] = 6.4e-28 * nh2 * (fH/2.+1.)*(1.+4.*fHe)
        else: 
            self.H['cr'][self.i] = 6.4e-28 * nh2
            
        
        
    def Ch2(self):
    
        '''
        Calculate the line cooling rate by vibrational excitation of H_2.
        
        Two modes are available (used by MCP/ALI and GASTRoNOoM, respectively):
        1) Following Groenewegen 1994, based on Hartquist et al 1980. Based
        on fitting of tabulated data of H2 cooling under LTE conditions
        2) Following Decin et al 2006, based on GS1976, Hollenbach & McKee
        1979 and Hollenbach & McKee 1989.
        
        The keyword h2_method (groenewegen, or decin) determines which of the 
        two is used. 
        
        '''
        
        #-- if cooling rate has already been calculated: don't do anything
        if self.C['h2'].has_key(self.i): return 
        
        #-- Initialise density
        self.setDensity('gas')
        nh2 = self.gdens.eval(dtype='nh2')
        T = self.T.eval(inner_eps=self.inner_eps,warn=not self.inner)
        
        #-- Method by Groenewegen 1994
        if self.pars['h2_method'].lower() == 'groenewegen':
            self.C['h2'][self.i] = 2.6111e-21*nh2*(T/1000.)**4.74
        
        #-- Method by Decin et al 2006
        else:
            #-- Based on deexcitation collisions between H2-H2, and H2-H, so 
            #   need H-number density
            nh = self.gdens.eval(dtype='nh')
            
            #-- Spontaneous emission rate from first vibrational state of H2 
            #   (in s^-)1
            A10 = 3e-7
            
            #-- Energy (in erg) of the emitted photon (0.6 eV converted)
            hnu10 = 9.61305939e-13 
            
            #-- The collisional deexcitation rate factor alpha
            coll_deex_h = 1.0e-12*np.sqrt(T)*np.exp(-1000./T)
            coll_deex_h2 = 1.4e-12*np.sqrt(T)*np.exp(-18100./(T+1200.))
            alpha = nh*coll_deex_h + nh2*coll_deex_h2
            
            #-- Nominator and denominator for n1
            nomin = alpha*np.exp(-hnu10/k_b/T)
            denomin = alpha*(1+np.exp(-hnu10/k_b/T))+A10
            
            #-- Number density of vibrationally excited molecules:
            n1 = nh2*nomin/denomin
            self.C['h2'][self.i] = A10 * hnu10 * n1
        
        

    def Clc(self,m,update_pop=1):
    
        '''
        Calculate the radiative line cooling rate for a molecule.
        
        For GS2014, no correction term yet for the excitation temperature per 
        level. 
        
        Follows Sahai 1990.
        
        Cubic spline interpolation for the level populations. Linear 
        interpolation and extrapolation for the collision rates (as for 
        GASTRoNOoM). 
        
        Currently not yet implemented to use a sqrt(T/T0) extrapolation at lower
        boundary as is done by MCP/ALI.
        
        Note that EnergyBalance will internally call this function with the 
        molecule tag tacked onto the Clc function name. 
        
        @param m: The molecule name from the input molecules list.
        @type m: str
        @keyword update_pop: Check if the level populations have been updated
        
                             (default: 1)
        @type update_pop: bool        
        
        '''
        
        def calcLevelLC(llow,r,T):
        
            '''
            Calculate the line cooling contribution for a single transition 
            index, with a fixed lower level i and for which j > i, for the 
            entire radial grid.
            
            Keep in mind, the goal is to include all transitions from every 
            level to every level. It doesn't actually matter if the energy is 
            lower or higher: The Sahai + Einstein equations change the sign if 
            the lower level is really the upper level in terms of energy. 
            
            We can use the structure of the collision rate files, where all 
            possible collisional transitions are included. By returning all 
            transitions to a given lower level, we already have j > i. 
            
            A note must be made here. Normally one would want to work with 
            all levels that have higher energy than Elow. However, for CO 
            this leads to issues because some v=1 levels have lower energy
            than some v=0 levels, while they are still sorted going v=0 to 
            jmax, then v=1 to jmax, ie not sorted by energy. The collision
            rates however assume that they are sorted by energy. Meaning 
            that for some collisional transitions Eup-Elow becomes < 0 
            because of how the CO spectroscopy is sorted. This is not 
            necessarily a problem, hence why we assume Eup-Elow must be > 0
            in what follows, and force it to be through abs(Eup-Elow). 
            We assume the collision rate files are sorted properly, thus 
            retrieve all "higher energy levels" by simply passing the llow
            index to the getTI method. This leads to results that are 
            identical with GASTRoNOoM CO cooling rates.
            
            
            @param llow: index of the transition lower level.
            @type llow: int
            @param r: The radial grid in cm
            @type r: array
            @param T: The temperature in K for which to calculate the cooling
            @type T: array
            
            @return: The contribution to the line cooling for a given lower 
                     level i (in ergs * cm^3 / s). 
            @rtype: float
            
            '''
            
            #-- Energy, weight and population of the lower level
            Elow = self.mol[m].getLEnergy(index=llow,unit='erg')
            glow = self.mol[m].getLWeight(index=llow)
            poplow = self.pop[m].getInterp(llow)(r)
            
            #-- Get the transition indices that go to the level with index llow
            #   Enforce indices to be an array.
            indices = self.collis[m].getTI(itype='coll_trans',llow=(llow,))
            
            #-- In case no upper levels were found, llow is the highest level
            #   available, and there are no collisions to be taken into account
            #   return 0
            if not indices.size: 
                return 0.
            
            #-- Retrieve the level indices, energies, weights and populations of
            #   upper levels
            lups = self.collis[m].getTUpper(index=indices,itype='coll_trans') 
            Eups = self.mol[m].getLEnergy(lups,unit='erg')
            popups = np.array([self.pop[m].getInterp(lup)(r) for lup in lups])
            gups = self.mol[m].getLWeight(index=lups)
            
            #-- Force gups/Eups into a column vector for array multiplication
            gups.shape = (gups.size,1)
            Eups.shape = (Eups.size,1)
            
            #-- Retrieve the rates between llow and all lups
            #   Also get the respective weights and energies
            Culs = np.array([self.collis[m].getInterp(i)(T) for i in indices])
            
            #-- Calculate the reversed rate for lower to upper level.
            #   Based on the Einstein relation: Cul/Clu = gl/gu exp(Eul/kT)
            expfac = np.exp(np.outer(-1*abs(Elow-Eups),1./(k_b*T)))
            Clus = Culs*np.multiply(expfac,gups/glow) 
            
            #-- Calculate the sum of the cooling contribution across all j>i 
            #   transitions
            return np.sum((Clus*poplow-Culs*popups)*abs(Eups-Elow),axis=0)    
        
        
        #-- if cooling rate has already been calculated: don't do anything
        if self.C['lc_{}'.format(m)].has_key(self.i): 
            return 
        
        #-- Initialise the gas density
        self.setDensity('gas')
        
        #-- Check if the level populations have to be updated
        if update_pop: self.updatePop(m)
        
        #-- Else if the level populations have not yet been updated, use the 
        #   cooling term calculated for the temperature for which the level pops
        #   were calculated. Note that self.pop[m][0] always exists. It is 
        #   calculated upon initialisation.
        print "Passing through Clc at iteration %i"%self.i
        if self.i > self.ipop[m]:
            print "Taking previous Clc for iteration %i"%self.ipop[m] 
            j = self.ipop[m]
            self.C['lc_{}'.format(m)][self.i] = self.C['lc_{}'.format(m)][j]
            return
        
        print "Calculating new Clc for iteration %i"%self.i
        #-- Extract the relevant density profiles
        nh2 = self.gdens.eval(dtype='nh2')
        amol = self.abun[m].eval(warn=0)

        #-- Calculate the line cooling term for this molecule.
        llows = self.pop[m].getLI()
        LCterm = np.empty(shape=(len(self.r),len(llows)))
        for i in llows: 
            LCterm[:,i-1] = calcLevelLC(i,self.r,\
                                        self.T.eval(inner_eps=self.inner_eps,\
                                                    warn=not self.inner))
        LCtotal = nh2*nh2*amol*np.sum(LCterm,axis=1)
        
        #-- Do NOT Multiply by -1. This already gives the net energy lost
        self.C['lc_{}'.format(m)][self.i] = LCtotal



    def plotRateIterations(self,iterations=[],dTsign='C',mechanism='ad',\
                           scale=1,fn=None,cfg=None,**kwargs):

        '''
        Plot the heating or cooling rates for several iterations of a given
        mechanism and heat exchange sign (cooling or heating).
        
        @keyword iterations: The iteration indices to be plotted. If default, 
                             all iterations are plotted
                             
                             (default: [])
        @type iterations: list
        @keyword dTsign: The dT type: heating (H) or cooling (C).
        
                         (default: 'C')
        @type dTsign: str
        @keyword mechanism: The cooling/heating mechanism to be plotted.
                            
                            (default: 'ad')
        @type mechanism: type
        @keyword scale: Scale the heating and cooling rates with r^4/1e14cm.
        
                        (default: 1)
        @type scale: bool
        @keyword fn: The filename and path of the plot (without extension). If
                     default, a filename can be given in a cfg file, and if that
                     is not the case, the plot is simply shown but not saved.
                     
                     (default: None)
        @type fn: str
        @keyword cfg: The filename to the cfg file for this plot with extra 
                      settings. 
                      
                      (default: None)
        @type cfg: str
        
        @keyword kwargs: Additional keywords passed to the plotting method. 
                         Overwrites any keys given in the cfg file.
                          
                         (default: {})
        @type kwargs: dict
        
        '''
        
        #-- Read the dict if available, and extract the filename. 
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn = cfg_dict.pop('filename')
        
        #-- If not iterations are specified, plot all of them.
        if not iterations:
            iterations = range(1,self.i)
            
        #-- Make sure we have a list, and proper input format. Adapt filename
        iterations = Data.arrayify(iterations)
        #if 0 in iterations: iterations = [i for i in iterations if i!=0]
        dTsign = dTsign[0].upper()
        if not fn is None:
            fn += '_{}{}'.format(dTsign,mechanism)
            fn += '_iter{}-{}'.format(iterations[0],iterations[-1])
            ddict['filename'] = fn
        
        #-- Set parameters. kwargs > cfg > locally defined
        pars = {'xlogscale': 1, 'ylogscale': 1, 'xaxis': 'r (cm)',
                'extension': 'pdf'}
        pars.update(cfg_dict)
        pars.update(kwargs)
        
        #-- Additional settings for plot
        pmech = mechanism.replace('_','\_')
        ddict = {'yaxis': dTsign+pmech+'-rate (ergs s$^{-1}$ cm$^{-3}$)',
                 'keytags': [str(i) for i in iterations], 'x': [], 'y': []}
        if scale:
            ddict['yaxis'] += r' $\times$ (r/10$^{14}$ cm)$^4$'

        #-- Where are the mechanisms? Only select first two characters in mech
        mech_types = {'ad': 'C', 'dg': 'H', 'dt': 'H', 'lc': 'C',
                      'pe': 'H', 'h2': 'C', 'cr': 'H' }
        mtype = mech_types[mechanism[:2]]
        
        #-- Set the sign. Only select first two characters in mech
        Csigns = {'ad': 1, 'dg': -1, 'dt': -1, 'lc': 1,
                  'pe': -1, 'h2': 1, 'cr': -1 }
        sign = Csigns[mechanism[:2]]
        if dTsign == 'H': 
            sign = -1 * sign
            
        #-- The heating/cooling terms for several iterations
        for i in iterations:
            isel = getattr(self,mtype)[mechanism][i]*sign > 0 
            ix = self.r[isel]
            ddict['x'].append(ix)
            iy = getattr(self,mtype)[mechanism][i][isel]*sign
            if scale: iy = iy*(ix*1e-14)**4
            ddict['y'].append(iy)

        #-- Add in other parameters
        if abs(min([min(yi) for yi in ddict['y']])) < 1e-20: 
            ddict['ymin'] = 1e-20
        ddict.update(pars)
        fn = Plotting2.plotCols(**ddict)
        
        if not fn is None:
            print('Heating and cooling rates plotted at:')
            print(fn+'.pdf')
        
        
        
    def plotRates(self,scale=1,fn=None,iteration=None,cfg=None,join=0,**kwargs):
    
        '''
        Plot the heating and cooling rates for a single iteration.
        
        @keyword scale: Scale the heating and cooling rates with r^4/1e14cm.
        
                        (default: 1)
        @type scale: bool
        @keyword fn: The filename and path of the plot (without extension). If
                     default, a filename can be given in a cfg file, and if that
                     is not the case, the plot is simply shown but not saved.
                     
                     (default: None)
        @type fn: str
        @keyword iteration: The iteration for which to plot the rates. Default 
                            is the last iteration. 
                            
                            (default: None)
        @type iteration: int
        @keyword cfg: The filename to the cfg file for this plot with extra 
                      settings. 
                      
                      (default: None)
        @type cfg: str
        @keyword join: Join together the plot for heating and cooling terms.
                       
                       (default: 0)
        @type join: bool
        
        @keyword kwargs: Additional keywords passed to the plotting method. 
                         Overwrites any keys given in the cfg file.
                          
                         (default: {})
        @type kwargs: dict
        
        '''
        
        #-- Read the dict if available, and extract the filename. 
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn = cfg_dict.pop('filename')
            
        #-- Check which iteration to plot, use default if not give or invalid. 
        if iteration is None or iteration > self.i: iteration = self.i
        
        #-- Add the iteration if any but the last one is requested
        if fn and iteration != self.i: fn += '_iter{}'.format(self.i)
        
        #-- Default line_types
        lts = {'ad': '-.g', 'dg': '-.b', 'dt': '--r',
               'pe': '--g', 'h2': '-.k', 'cr': '-k' }
        mkeys = [k for k in self.C.keys() if k[:2] == 'lc']
        lts_extra = ['-r', '--b', '--k', '-m', '--m', '-y', '--y', '-g']
        for mk in mkeys:
            lts[mk] = lts_extra.pop(0)
        
        #-- Set parameters. kwargs > cfg > locally defined
        pars = {'xlogscale': 1, 'ylogscale': 1, 'xaxis': 'r (cm)',
                'extension': 'pdf'}
        pars.update(cfg_dict)
        pars.update(kwargs)
        
        #-- Iterate over heating and cooling type, two plots total
        fns = []
        for sign in ['H','C']:
            ddict = {'yaxis': sign+'-rate (ergs s$^{-1}$ cm$^{-3}$)',
                     'keytags': [], 'x': [], 'y': [], 'line_types': []}
            
            if scale:
                ddict['yaxis'] += r' $\times$ (r/10$^{14}$ cm)$^4$'
            
            #-- Set the filename if one is given with a suffix
            if fn: 
                ddict['filename'] = '{}_{}'.format(fn,sign)
            
            #-- The heating (alt. cooling) terms
            for term in sorted(getattr(self,sign).keys()):
                sterm = '$_\\mathrm{{{s}}}$'.format(s=term.replace('_','\_'))
                ddict['keytags'].append(sign+sterm)
                isel = getattr(self,sign)[term][iteration] > 0 
                ix = self.r[isel]
                ddict['x'].append(ix)
                iy = getattr(self,sign)[term][iteration][isel]
                if scale: iy = iy*(ix*1e-14)**4
                ddict['y'].append(iy)
                ddict['line_types'].append(lts[term])
            
            #-- The negative cooling (alt. heating) terms
            negsign = 'C' if sign == 'H' else 'H'
            for term in sorted(getattr(self,negsign).keys()):
                isel = getattr(self,negsign)[term][iteration] < 0 
                if not isel.any(): continue
                sterm = '$_\\mathrm{{{s}}}$'.format(s=term.replace('_','\_'))
                ddict['keytags'].append(sign+sterm)
                ix = self.r[isel]
                ddict['x'].append(ix)
                iy = -1.*getattr(self,negsign)[term][iteration][isel]
                if scale: iy = iy*(ix*1e-14)**4
                ddict['y'].append(iy)
                ddict['line_types'].append(lts[term])

            #-- Add in other parameters
            if abs(min([min(yi) for yi in ddict['y'] if yi.size])) < 1e-20: 
                ddict['ymin'] = 1e-20
            ddict.update(pars)
            ifn = Plotting2.plotCols(**ddict)
            fns.append(ifn)
        
        #-- Print the plotted filenames, if any
        if join and fns[-1] and fns[-1][-3:] == 'pdf': 
            DataIO.joinPdf(old=fns,new=fn+'.pdf',del_old=1)
            print('Heating and cooling rates plotted at:')
            print(fn+'.pdf')
        elif fns[-1]:
            print('Heating and cooling rates plotted at:')
            print('\n'.join(fns))
    
    
    
    def plotT(self,fn=None,iterations=[],cfg=None,**kwargs):
    
        '''
        Plot the temperature profiles of the different iterations. 
        
        @keyword fn: The filename and path of the plot (without extension). If
                     default, a filename can be given in a cfg file, and if that
                     is not the case, the plot is simply shown but not saved.
                     
                     (default: None)
        @type fn: str
        @keyword iterations: The iteration indices to be plotted. If default, 
                             all iterations are plotted. The first and last 
                             iterations are always plotted.
                             
                             (default: [])
        @type iterations: list
        @keyword cfg: The filename to the cfg file for this plot with extra 
                      settings. 
                      
                      (default: None)
        @type cfg: str
        @keyword kwargs: Additional keywords passed to the plotting method. 
                         Overwrites any keys given in the cfg file.
                         
                         (default: {})
        @type kwargs: dict
        
        '''
        
        #-- Read the dict if available, and extract the filename. 
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn = cfg_dict.pop('filename')
        
        #-- If no iterations are specified, plot all of them.
        if iterations: 
            iterations = [i for i in iterations if i < self.i and i > 0]
        elif not iterations:
            iterations = range(1,self.i)
        
        #-- Add first and last iteration
        iterations = [0] + iterations + [self.i]
        
        #-- Set parameters. kwargs > cfg > locally defined
        pars = {'xlogscale': 1, 'ylogscale': 1, 'yaxis': 'T (K)',
                'xaxis': 'r (cm)', 'filename': fn, 
                'keytags': [r'T$_\mathrm{d}$']
                          +['Iteration {}'.format(i) for i in iterations],
                'extension': 'pdf'}
        pars.update(cfg_dict)
        pars.update(kwargs)
        
        #-- Set the T-profiles and plot
        y = [self.Td.eval()]+[self.T_iter[i].eval(inner_eps=self.inner_eps,\
                                                  warn=not self.inner) 
                              for i in iterations]
        pfn = Plotting2.plotCols(x=self.r,y=y,**pars)
        
        #-- Print the plotted filename
        if pfn: 
            print('The temperature profiles are plotted at:')
            print(pfn)
        
    
    
    def plotHCTerm(self,fn=None,iterations=[],dTsign='C',cfg=None,**kwargs):
        
        '''
        Plot the total heating and cooling term (excluding adiabatic cooling)
        calculated by calcT before the differential equation is solved. 
        
        This includes the density factor that enters, and so is essentially 
        (H - C)/rho. The velocity does not enter here, since that is calculated
        explicitly in dTdr. Hence the y-axis units K/s.
        
        @keyword fn: The filename and path of the plot (without extension). If
                     default, a filename can be given in a cfg file, and if that
                     is not the case, the plot is simply shown but not saved.
                     
                     (default: None)
        @type fn: str
        @keyword iterations: The iteration indices to be plotted. If default, 
                             all iterations are plotted. The first and last 
                             iteration is always plotted.
                             
                             (default: [])
        @type iterations: list
        @keyword dTsign: The dT type: heating (H) or cooling (C). Either the 
                         positive (C) or the negative (H) sum is plotted. 
        
                         (default: 'C')
        @type dTsign: str
        @keyword cfg: The filename to the cfg file for this plot with extra 
                      settings. 
                      
                      (default: None)
        @type cfg: str
        @keyword kwargs: Additional keywords passed to the plotting method. 
                         Overwrites any keys given in the cfg file.
                         
                         (default: {})
        @type kwargs: dict
        
        '''
        
        #-- Read the dict if available, and extract the filename. 
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn = cfg_dict.pop('filename')
        
        #-- If no iterations are specified, plot all of them.
        if iterations: 
            iterations = [i for i in iterations if i < self.i and i > 1]
        elif not iterations:
            iterations = range(2,self.i)
        
        #-- Add first and last iteration
        iterations = [1] + iterations + [self.i]
        
        #-- Set parameters. kwargs > cfg > locally defined
        pars = {'xlogscale': 1, 'ylogscale': 1, 
                'yaxis': 'H - C (K s$^{-1}$)',
                'xaxis': 'r (cm)', 'filename': fn, 
                'keytags': ['Iteration {}'.format(i) for i in iterations],
                'extension': 'pdf'}
        pars.update(cfg_dict)
        pars.update(kwargs)
        
        #-- Set the H - C terms and plot
        dTsign = -1 if dTsign.upper() == 'H' else 1
        y = [self.rates[i]*dTsign for i in iterations]
        pfn = Plotting2.plotCols(x=self.r,y=y,**pars)
        
        #-- Print the plotted filename
        if pfn: 
            print('The H-C terms are plotted at:')
            print(pfn)
            
            
    
    def plotTexc(self,m,indices=None,llow=None,lup=None,fn=None,cfg=None,\
                 **kwargs):
        
        '''
        Plot the excitation temperature as a function of impact parameter for a
        selection of radiative transitions, given by either the transition 
        indices or the upper and/or lower level indices. 
        
        Retrieves the excitation temperature from self.Texc. They are 
        appropriate for the level populations used by the current iteration, and
        are updated when the level populations are updated.
        
        At this time, plotting Texc for other iterations is not possible.
        
        @param m: The molecule tag
        @type m: str
        
        @keyword indices: The transition indices. Can be indirectly given 
                          through upper and/or lower level indices. If default
                          and no llow/lup are given, all Texc are returned for
                          this molecule.
                          
                          (default: None)
        @type indices: list/array
        @keyword lup: The upper level indices used to extract the transition
                      indices. If this or llow are defined, the keyword indices 
                      is overwritten.
                      
                      (default: None)
        @type lup: list/array
        @keyword llow: The upper level indices used to extract the transition
                       indices. If this or lup are defined, the keyword indices 
                       is overwritten.
                      
                       (default: None)
        @type llow: list/array
        @keyword fn: The filename and path of the plot (without extension). If
                     default, a filename can be given in a cfg file, and if that
                     is not the case, the plot is simply shown but not saved.
                     
                     (default: None)
        @type fn: str
        @keyword cfg: The filename to the cfg file for this plot with extra 
                      settings. 
                      
                      (default: None)
        @type cfg: str
        @keyword kwargs: Additional keywords passed to the plotting method. 
                         Overwrites any keys given in the cfg file.
                         
                         (default: {})
        @type kwargs: dict
        
        
        
        '''
        
        #-- Read the dict if available, and extract the filename. 
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn = cfg_dict.pop('filename')
            
        #-- Check which iteration to plot, use default if not give or invalid. 
        if iteration is None or iteration > self.i: iteration = self.i
        
        #-- Add the iteration if any but the last one is requested
        if fn and iteration != self.i: fn += '_iter{}'.format(self.i)
        
        #-- If either lup or llow are given, get transition indices from them
        if not (lup is None and llow is None):
            if not lup is None: lup = Data.arrayify(lup)
            if not llow is None: llow = Data.arrayify(llow)
            indices = self.mol[m].getTI(lup=lup,llow=llow)
        
        #-- Redefine lup and llow to  match the given transitions. 
        #   Do it here because lup/llow are needed for plot keytags. 
        lup = self.mol[m].getTUpper(indices)
        llow = self.mol[m].getTLower(indices)

        #-- Get the impact parameter grid
        p = self.pop[m].getP()
        
        #-- Get kinetic T profile and the excitation temperatures
        T = self.T.eval(p,inner_eps=self.inner_eps,warn=not self.inner)
        Texc = self.getTexc(m=m,indices=indices)
                
        #-- Set parameters. kwargs > cfg > locally defined
        pars = {'xlogscale': 1, 'ylogscale': 0, 
                'yaxis': 'T$_\mathrm{exc}$ (K)',
                'xaxis': 'p (cm)', 'filename': fn, 
                'extension': 'pdf',
                'keytags': ['T$_\mathrm{kin}$']+
                           ['llup = {}, llow = {}'.format(nui,nli)
                            for nui,nli in zip(lup,llow)]}
        pars.update(cfg_dict)
        pars.update(kwargs)
        
        #-- Make the plot
        y = [T]+[Texc[i,:] for i in range(Texc.shape[0])]
        pfn = Plotting2.plotCols(x=p,y=y,**pars)
