# -*- coding: utf-8 -*-

"""
Module for calculating the energy balance. 

Author: R. Lombaert

"""

import os, collections
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline1d
from scipy.integrate import odeint, trapz 
from astropy import constants as cst
from astropy import units as u

from cc.data import Data
from cc.plotting import Plotting2
from cc.tools.io import DataIO
from cc.tools.numerical import Operators as op
from cc.tools.numerical import Gridding
from cc.modeling.profilers import Velocity, Density, Profiler, Grainsize
from cc.modeling.profilers import Mdot, Radiance, Opacity, Temperature
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
                    Default if only adiabatic cooling is taken into account
              
                    (default: None)
    @type rates: Profiler()
    @keyword warn: Warn when extrapolation occurs.
    
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
    
    Calculates energy balance and provides tools for reading and writing I/O, 
    and setting input physics.
    
    For now limited to the energy balance from the dust formation radius up to
    the outer radius.
    
    NYI: line cooling/heating, comisc ray heating, photoelectric heating.
    
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
        
        #-- Initialise all variables
        self.__reset()
        
        #-- Set default input template, and the custom parameter file
        self.template = template.lower()
        if self.template not in ['mcp','gastronoom']: self.template = 'standard'
        self.fn = fn
            
        #-- Read input parameters, set constants, coordinate grids, stellar
        #   properties, basic profiles, and the adiabatic coefficient
        self.readInputParameters(**kwargs)
        self.setGrids()
        self.setStar()
        self.setProfiles()
        self.setGamma()
        
        #-- Everything has been initialised. Set the current Temperature profile
        self.__setT()
    
    
    
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
                                    convert_floats=1)
        
        #-- Read from the filename if it is given. 
        if not self.fn is None: 
            dd = DataIO.readDict(self.fn,convert_lists=1,convert_floats=1)
            self.pars.update(dd)
        
        #-- Insert any extra arguments given in function call
        self.pars.update(kwargs)
        
        #-- Then initialise the heating and cooling terms requested
        self.H = {term: {} for term in self.pars['hterms']}
        self.C = {term: {} for term in self.pars['cterms']}
        
        #-- Adiabatic cooling is always calculated
        self.C['ad'] = {}
        
    
    
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
        


    def setProfiles(self,v=None,opac=None,mdot=None,mdot_dust=None,T=None,\
                    Td=None):
    
        '''
        Initialise the basic profiles: gas velocity, opacity, mass-loss rate, 
        initial temperature, dust temperature 
        
        The grids can be given as lists, which are used as input for 
        the respective Profiler classes. If not given, they are taken from the 
        inputfile.

        @keyword v: The velocity profile (r, cm/s). [func,*pars]
                    
                    (default: None)
        @type v: list        
        @keyword opac: The opacity profiles (l, cm^2/g). [func,*pars]
                    
                       (default: None)
        @type opac: list
        @keyword mdot: The mass-loss rate profile (r, Msun/yr) [func,*pars] or 
                       constant.
                    
                       (default: None)
        @type mdot: list
        @keyword mdot_dust: The dust mass-loss rate profile (r, Msun/yr) 
                            [func,*pars] or constant.
                    
                            (default: None)
        @type mdot_dust: list
        @keyword T: The initial temperature profile (r, K). [func,*pars]
                    The first arg is the function, second argument always the 
                    initial temperature. Can be string Td as well to set it to
                    the dust temperature.
                    
                    (default: None)
        @type T: list/str
        @keyword Td: The dust temperature profile (r, K). [func,*pars]
                     The first arg is the function.
                    
                     (default: None)
        @type Td: list
        
        '''

        #-- Only reset if self.T is already yet defined.
        if self.T: self.__reset()
        
        #-- Take the profile definitions from the inputfile if not given here
        if v is None: v = self.pars['v']
        if opac is None: opac = self.pars['opac']
        if mdot is None: mdot = self.pars['mdot']
        if mdot_dust is None: mdot_dust = self.pars['mdot_dust']
        if T is None: T = self.pars['Tinit']
        if Td is None: Td = self.pars['Td']
                
        #-- Set the velocity profile: r, func, dfunc, order interp dfunc, pars
        self.v = Velocity.Velocity(self.r,v[0],None,3,*v[1:])
        
        #-- Set the opacity profile:
        #   Check if KappaReader is requested for opacities. opac then contains
        #   the parameters for the KappaReader interpolator.
        if opac.pop(0) == 'KappaReader':
            #-- Reads the requested opacity with index + order, and interpolates
            opac = [Opacity.read_opacity(*opac)]
        
        #-- func will be the only element in opac if KappaReader is requested
        self.opac = Opacity.Opacity(self.l,opac[0],None,3,*opac[1:])
        
        #-- Set the mass-loss-rate profiles. Check if these are constant.
        mdot, mdot_dust = Data.arrayify(mdot), Data.arrayify(mdot_dust)
        if mdot.size < 2:
            self.mdot = mdot
        else:
            self.mdot = Mdot.Mdot(self.r,mdot[0],None,3.,*mdot[1:])
        if mdot_dust.size < 2:
            self.mdot_dust = mdot_dust
        else:
            self.mdot_dust = Mdot.Mdot(self.r,mdot_dust[0],None,3.,\
                                       *mdot_dust[1:])
        
        #-- Set the dust temperature profile
        self.Td = Temperature.Temperature(self.r,Td[0],None,3,*Td[1:])

        #-- Set the initial temperature profile and the temperature at the inner
        #   boundary.
        if isinstance(T,str) and T.lower() == 'td':
            self.T_iter[0] = self.Td
        else:
            self.T_iter[0] = Temperature.Temperature(self.r,T[0],None,3,*T[1:])
        self.T0 = T[0]
        
        
        
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
                kwargs = {'x': self.T_iter[0].eval(), 
                          'func': Profiler.step,
                          'ylow': 5./3., 'yhigh': 7./5., 'xstep': 350. }
                gamma = Profiler.Profiler(**kwargs)
            #-- Interpolate the adiabatic coefficient taken from NIST database 
            else: 
                fn = os.path.join(cc.path.aux,'h2_physical_properties.dat')
                kwargs = {'x': self.T_iter[0].eval(),
                          'func': Profiler.interp_file,
                          'ikwargs': {'ext': 3, 'k': 3}, 
                          'filename': fn, 'ycol': -1}
                gamma = Profiler.Profiler(**kwargs)
        
        #-- Alternatively, gamma is constant, and gives the value of the 
        #   coefficient
        else: 
            kwargs = {'x': self.T_iter[0].eval(), \
                      'func': Profiler.constant, 'c': self.pars['gamma']}
            gamma = Profiler.Profiler(**kwargs)
        
        self.gamma = gamma



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
        
        '''

        if self.w is None:
            kwargs = {'r': self.r, 'a': self.a, 'func': Velocity.driftRPDF,
                      'l': self.l, 'v': self.v, 'mdot': self.mdot,
                      'P': self.pars['P'], 'sd': self.pars['sd'],
                      'opac': self.opac, 'radiance': self.rad,
                      'w_thermal': self.pars['w_thermal'],
                      'alpha': self.pars['alpha'],
                      'mu': self.gdens.getMeanMolecularWeight()}
            
            #-- Include T if a thermal velocity term is needed
            if self.pars['w_thermal'].lower() in ['kwok','mean','rms','prob','epstein']:
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
        
        
        
    def setDensity(self,dtype,custom_mu=None): 
    
        '''
        Calculate the density profile for gas or dust. 
        
        The dust density profile is not dependent on grain size, since the drift
        averaged over grain size.
        
        @param dtype: The density profile type, being 'gas' or 'dust'.
        @type dtype: str
        
        @keyword custom_mu: Set a custom mean molecular weight. Mainly for 
                            testing purposes. If default, calculated from fH and
                            fHe.
                            
                            (default: None)
        @type custom_mu: float.
        
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
                if custom_mu:
                    self.gdens.mu = custom_mu
        
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
        
    
    
    def __setT(self):
    
        '''
        Set the current temperature profile for this iteration. When i == 0,
        this is the initial guess.
        
        This is done explicitly, because upon initialisation the state of the
        object depends on self.T being None or not.
        
        '''
        
        self.T = self.T_iter[self.i]
        
        
    
    def iterT(self,conv=0.01,imax=50,step_size=0.05,dTmax=0.20,warn=1,\
              *args,**kwargs):
    
        '''
        Iterate the temperature profile until convergence criterion is reached.
        
        Extra arguments are passed on to calcT().
        
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
                            allows.
                            
                            (default: 0.05)
        @type step_size: float
        @keyword dTmax: The starting value for the maximum allowed relative
                        temperature change between iterations. This maximum 
                        increases as the code reaches convergences, through a 
                        set multiple of step_size.
                            
                        (default: 0.20)
        @type dTmax: float
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        
        '''
        
        print('-- Iterating T(r) now.')
        dTnsteps = (1.-dTmax)/step_size
        steps = 0.
        
        #-- If first iteration, calculate one step in any case
        if self.i == 0:
            print('Iteration 1 for T(r)...')
            self.calcT(dTmax)
        
        #-- Then continue iterating until the relative difference between this 
        #   and the previous result is smaller than convergence criterion in all
        #   radial points, or when the maximum number of iterations is reached
        while True:
            dTdiff = np.abs(1.-self.T.eval()/self.T_iter[self.i-1].eval())
            if (np.all(dTdiff<conv) and steps == dTnsteps) or self.i == imax: 
                break
            if np.all(dTdiff<(dTnsteps-steps)*conv):
                if dTnsteps-steps > 8.: steps += 4.
                elif dTnsteps-steps > 4.: steps += 2.
                else: steps += 1.
            print('Iteration {} for T(r)...'.format(self.i+1))
            self.calcT(dTmax+step_size*steps,warn=warn)
            if not np.all(np.isfinite(self.T.eval())):
                print('nans found in T-profile. Breaking off iteration.')
                break
    
    
    
    def calcT(self,dTmax=1,order=3,warn=1):
    
        '''
        Calculate the temperature profile based on the differential equation 
        given by Goldreich & Scoville (1976). The function is defined in dTdr.
        
        The initial temperature is taken to be the second argument of the 
        initial temperature profile given by T in inputEnergyBalance.dat.  
        This is typically the condensation temperature.
        
        @keyword dTmax: The maximum allowed relative temperature change for this
                        T calculation. Set to 100% by default.
                            
                        (default: 1)
        @type dTmax: float
        @keyword order: The spline1d interpolation order for the rates and 
                        temperatures.
                            
                        (default: 3)
        @type order: int
        @keyword warn: Warn when extrapolation occurs.
        
                       (default: 1)
        @type warn: bool
        
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
        self.rates[self.i] = yrates
        
        #-- Add the density term
        yrates = yrates/dens_term
        
        #-- Note that if no rates are requested, this will give an empty array
        if self.r.size == Data.arrayify(yrates).size:
            rp = Profiler.Profiler(self.r,spline1d(self.r,yrates,k=order,ext=1))
        else:
            rp = None 
        
        #-- Calculate the next iteration of the temperature profile
        kwargs = {'func': dTdr, 'y0': self.T.eval()[0], 't': self.r,
                  'args': (self.v,self.gamma,rp,warn)}
        Tr = odeint(**kwargs)[:,0]
        
        #-- Check the temperature change. Limit it to the requested maximum 
        #   allowed change. Only do this from the second iteration, to allow 
        #   a big jump from the initial condition.
        if self.i < 0: 
            dTmax = self.pars['dTmax']
        print('Changing T by {:.1f}%.'.format(dTmax*100.))
        Ti = self.T.eval()        
        Tr = np.where(Tr<=0.,np.ones_like(Tr),Tr)
        Tr = Ti + dTmax*(Tr-Ti)
        
        #-- Create a profiler for the new temperature structure. Extrapolation
        #   done by returning boundary values, ie T0 and the temp at outer 
        #   boundary
        Tinterp = spline1d(self.r,Tr,k=order,ext=3)
        self.T_iter[self.i] = Temperature.Temperature(self.r,Tinterp)
        
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
        gamma = self.gamma.eval(self.T.eval(),warn=0)
        
        #-- Makes the rate negative, in accordance with -C in the second term of
        #   the diff eq.
        factor1 = 2-2*gamma
        factor2 = gamma-1
        
        #-- Determine the velocity and its derivative from the Velocity() object
        vi = self.v.eval()
        dvi = self.v.diff()
        
        #-- Calculate the conversion term K/cm => erg/s/cm3: allows to compare 
        #   C terms with adiabatic term.
        self.setDensity('gas')
        nh2 = self.gdens.eval(dtype='nh2')
        ascale = self.gdens.fH + 1. + self.gdens.fHe * (self.gdens.fH + 2.)
        conversion = nh2 * k_b * vi * ascale / factor2
        
        #-- Calculate the temperature dependent term
        main_term = (1./self.r + 0.5*1./vi*dvi) * self.T.eval()
        
        #-- Calculate the rate. Multiply by -1 to have positive cooling rate.
        self.C['ad'][self.i] = -1. * conversion * factor1 * main_term
    
    
    
    def Hdg(self):
    
        '''
        Calculate the heating rate by dust-gas collisions. 
        
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
            vT = 4./3.*np.sqrt(8.*k_b*self.T.eval()/(self.gdens.mu*np.pi*mh))
            
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
            #   self.nd is a Grainsize.Distribution object!
            vT = np.transpose([vT])
            nd = self.nd.eval()
            w = self.w.eval()
            afac = nd*self.a**2.*w**2.*(vT**2.+w**2.)**0.5
            gsfac = trapz(x=self.a,y=afac,axis=1)
        else: 
            #-- Avg grain size: no integration
            #   self.nd refers to the self.ddust, or a dust Density() object!
            nd = self.nd.eval(dtype='ndust')
            w = np.squeeze(self.w.eval())
            gsfac = nd*self.a**2*w**2.*(vT**2.+w**2.)**0.5

        #-- Set dust-gas collisional heating
        self.H['dg'][self.i] = cfac*gsfac
        


    def Hdt(self):
    
        '''
        Calculate the heating rate by dust-gas heat exchange. 
        
        The equation is derived from Burke & Hollenbach 1983, based on 
        Groenewegen 1994, Schoier et al. 2001, and Decin et al. 2006: 
        
        General case: 
        Hdt = 
        =>
        Hdt = 
        
        MCP/ALI case (average grain size a): 
        Hdt = 
        
        GASTRoNOoM case (MRN): 
        Hdt = 
        
        '''
        
        #-- if heating rate has already been calculated: don't do anything
        if self.H['dt'].has_key(self.i): return 
        
        #-- Initialise density and a few parameters
        self.setDensity('gas')
        self.setDensity('dust')
        nh2 = self.gdens.eval(dtype='nh2')
        fH = self.pars['fH']
        P = self.pars['P']
        mode = self.pars['heatmode']
        
        #-- Extract the dust and gas temperatures
        td = self.Td.eval()
        t = self.T.eval()
        
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
            #   self.nd is a Grainsize.Distribution object!
            nd = self.nd.eval()
            gsfac = trapz(x=self.a,y=nd*self.a**2,axis=1)
        else: 
            #-- Avg grain size: no integration
            #   self.nd refers to the self.ddust, or a dust Density() object!
            nd = self.nd.eval(dtype='ndust')
            gsfac = nd*self.a**2

        #-- Set the dust-gas collisional heating term. First mode is Gail & 
        #   Sedlmayr's equation from 2014, based on Draine 1980. Default is the
        #   accommodation heat exchange by Burke & Hollenbach 1983
        if mode.lower() == 'gs2014':
            #-- Avg over a^2. Normalisation disappears when combined with gsfac
            drift = self.w.avgDrift(norm_type='collisional',nd=self.nd)
            
            #-- fg/2. based on gamma, with fg = 2/(gamma-1)
            fg = 1./(self.gamma.eval(self.T.eval(),warn=0)-1.)
            
            #-- Velocity factor 
            vfac = (8.*vT**2+drift**2)**0.5
            self.H['dt'][self.i] = accom * cfac * fg * tdiff * gsfac * vfac
            
        else: 
            #-- Extra fac 2 enters through 2kT (BH1983)
            self.H['dt'][self.i] = 2. * accom * cfac * tdiff * vT * gsfac
            


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
        
        @keywords kwargs: Additional keywords passed to the plotting method. 
                          Overwrites any keys given in the cfg file.
                          
                          (default: {})
        @type kwargs: dict
                
        @return: The filename of the plot is returned, including extension.
        @rtype: str
        
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
        if 0 in iterations: iterations = [i for i in iterations if i!=0]
        dTsign = dTsign[0].upper()
        if not fn is None:
            fn += '_{}{}'.format(dTsign,mechanism)
            fn += '_iter{}{}'.format(iterations[0],iterations[-1])
            ddict['filename'] = fn
        
        #-- Set parameters. kwargs > cfg > locally defined
        pars = {'xlogscale': 1, 'ylogscale': 1, 'xaxis': 'r (cm)',
                'extension': 'pdf'}
        pars.update(cfg_dict)
        pars.update(kwargs)
        
        #-- Additional settings for plot
        ddict = {'yaxis': dTsign+mechanism+'-rate (ergs s$^{-1}$ cm$^{-3}$)',
                 'keytags': [str(i) for i in iterations], 'x': [], 'y': []}
        if scale:
            ddict['yaxis'] += r' $\times$ (r/10$^{14}$ cm)$^4$'

        #-- Where are the mechanisms?
        mech_types = {'ad': 'C', 'dg': 'H', 'dt': 'H', 'co': 'C', 'h2o': 'C',
                      'pe': 'H', 'h2': 'C', 'cr': 'H' }
        mtype = mech_types[mechanism]
        
        #-- Set the sign
        Csigns = {'ad': 1, 'dg': -1, 'dt': -1, 'co': 1, 'h2o': 1,
                  'pe': -1, 'h2': 1, 'cr': -1 }
        sign = Csigns[mechanism]
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
        ddict.update(pars)
        fn = Plotting2.plotCols(**ddict)
        
        if not fn is None:
            print('Heating and cooling rates plotted at:')
            print(fn+'.pdf')
        
        
        
    def plotRates(self,scale=1,fn=None,iteration=None,cfg=None,join=0,**kwargs):
    
        '''
        Plot the heating and cooling rates for this iteration.
        
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
        
        @keywords kwargs: Additional keywords passed to the plotting method. 
                          Overwrites any keys given in the cfg file.
                          
                          (default: {})
        @type kwargs: dict
                
        @return: The filename of the plot is returned, including extension.
        @rtype: str
        
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
        lts = {'ad': '-.g', 'dg': '-.b', 'dt': '--r', 'co': '-r', 'h2o': '--b',
               'pe': '--g', 'h2': '-.k', 'cr': '-k' }
        
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
                sterm = '$_\\mathrm{{{s}}}$'.format(s=term)
                ddict['keytags'].append(sign+sterm)
                isel = getattr(self,sign)[term][self.i] > 0 
                ix = self.r[isel]
                ddict['x'].append(ix)
                iy = getattr(self,sign)[term][self.i][isel]
                if scale: iy = iy*(ix*1e-14)**4
                ddict['y'].append(iy)
                ddict['line_types'].append(lts[term])
            
            #-- The negative cooling (alt. heating) terms
            negsign = 'C' if sign == 'H' else 'H'
            for term in sorted(getattr(self,negsign).keys()):
                isel = getattr(self,negsign)[term][self.i] < 0 
                if not isel.any(): continue
                sterm = '$_\\mathrm{{{s}}}$'.format(s=term)
                ddict['keytags'].append(sign+sterm)
                ix = self.r[isel]
                ddict['x'].append(ix)
                iy = -1.*getattr(self,negsign)[term][self.i][isel]
                if scale: iy = iy*(ix*1e-14)**4
                ddict['y'].append(iy)
                ddict['line_types'].append(lts[term])

            #-- Add in other parameters
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
                             all iterations are plotted. The last iteration is
                             always plotted.
                             
                             (default: [])
        @type iterations: list
        @keyword cfg: The filename to the cfg file for this plot with extra 
                      settings. 
                      
                      (default: None)
        @type cfg: str
        @keywords kwargs: Additional keywords passed to the plotting method. 
                          Overwrites any keys given in the cfg file.
                          
                          (default: {})
        @type kwargs: dict
        
        @return: The filename of the plot is returned, including extension.
        @rtype: str
        
        '''
        
        #-- Read the dict if available, and extract the filename. 
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn = cfg_dict.pop('filename')
        
        #-- If not iterations are specified, plot all of them.
        if iterations: 
            iterations = [i for i in iterations if i < self.i-1]
        elif not iterations:
            iterations = range(0,self.i-1)
        
        #-- Make sure we have a list, and proper input format.
        iterations = Data.arrayify(iterations)
        
        #-- Set parameters. kwargs > cfg > locally defined
        pars = {'xlogscale': 1, 'ylogscale': 1, 'yaxis': 'T (K)',
                'xaxis': 'r (cm)', 'filename': fn, 
                'keytags': [r'T$_\mathrm{d}$']\
                          +['Iteration {}'.format(i) for i in iterations]\
                          +['Iteration {}'.format(self.i)],
                'extension': 'pdf'}
        pars.update(cfg_dict)
        pars.update(kwargs)
        
        #-- Set the T-profiles and plot
        y = [self.Td.eval()]+[self.T_iter[i].eval() for i in iterations]\
           +[self.T_iter[self.i].eval()]
        pfn = Plotting2.plotCols(x=self.r,y=y,**pars)
        
        #-- Print the plotted filename
        if pfn: 
            print('The temperature profiles are plotted at:')
            print(pfn)
        
