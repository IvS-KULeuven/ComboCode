# -*- coding: utf-8 -*-

"""
Module including functions for stellar parameters, and the STAR class and 
its methods and attributes.

Author: R. Lombaert

"""

import types
from glob import glob
import os
from scipy import pi, log, sqrt
from scipy import array, exp, zeros
from scipy import integrate, linspace
from scipy import argmin,argmax, empty
from scipy.interpolate import interp1d
import operator
from numpy import savetxt

import cc.path
from cc.data import Data
from cc.tools.io import Database
from cc.tools.io import DataIO, Atmosphere
from cc.tools.numerical import Interpol
from cc.modeling.objects import Molecule
from cc.modeling.objects import Transition
from cc.modeling.tools import ColumnDensity
from cc.modeling.codes import MCMax


def getStar(star_grid,modelid,idtype='GASTRONOOM'):
    
    '''
    Grab a Star() object from a list of such objects, given a model id.
    
    If no modelid is found, an empty list is returned. If Star() objects are 
    found (even only one), a list of them is returned.
    
    Based on the cooling modelid.
    
    @param star_grid: the Star() objects
    @type star_grid: list[Star()]
    @param modelid: the given modelid for which the selection is made.
    @type modelid: string
    
    @keyword idtype: The type of model id
                    
                     (default: GASTRONOOM)    
    @type idtype: string

    @return: The models matching the modelid
    @rtype: list[Star()]
    
    '''
    
    modelid, idtype = str(modelid), str(idtype)
    return [s for s in star_grid if s['LAST_%s_MODEL'%idtype] == modelid]
    
    
    
def makeStars(models,id_type,path,code):
    
    '''
    Make a list of dummy Star() objects.

    @param models: model_ids for the new models
    @type models: list[string]
    @param id_type: The type of id (PACS, GASTRONOOM, MCMAX)
    @type id_type: string
    @param path: Output folder in the code's home folder
    @type path: string
    @param code: The code (which is not necessarily equal to id_type, such as 
                 for id_type == PACS)
    @type code: string

    @return: The parameter sets, mostly still empty!
    @rtype: list[Star()]
    
    '''
    
    extra_pars = dict([('path_'+code.lower(),path)])
    star_grid = [Star(example_star={'LAST_%s_MODEL'%id_type.upper():model},\
                      **extra_pars) 
                 for model in models]
    return star_grid
      

    
class Star(dict):
    
    """
    Star class maintains information about a stellar model and its properties.

    Inherits from dict.
    
    """



    def __init__(self,path_gastronoom='',path_mcmax='',extra_input=None,\
                 example_star=dict()):
        
        """
        Initiate an instance of the STAR class.
        
        @keyword path_gastronoom: path in ~/GASTRoNOoM/ for modeling out/input
                                  
                                  (default: None)
        @type path_gastronoom: string
        @keyword path_mcmax: the folder in ~/MCMax/ for modeling out/input
        
                             (default: None)
        @type path_mcmax: string
        @keyword example_star: if not None the STAR object is exact duplicate 
                               of example_star. Can be a normal dictionary as 
                               well. Paths are not copied and need to be given 
                               explicitly.
                                    
                               (default: None)
        @type example_star: dict or Star()                                  
        @keyword extra_input: extra input that you wish to add to the dict
        
                              (default: None)
        @type extra_input: dict or Star()
        
        @return: STAR object in the shape of a dictionary which includes all 
                 stellar data available, if None are passed for both options it
                 is an empty dictionary; if an example star is passed it has a 
                 dict that is an exact duplicate of the example star's dict
        @rtype: Star()
        
        """    
            
        super(Star, self).__init__(example_star)
        if extra_input <> None: self.update(extra_input)
        self.Rsun = 6.95508e10         #in cm  Harmanec & Prsa 2011
        self.Msun = 1.98547e33      #in g   Harmanec & Prsa 2011
        self.Mearth = 5.97237e27        # in g
        self.Tsun = 5779.5747            #in K   Harmanec & Psra 2011
        self.Lsun = 3.846e33           #in erg/s
        self.year = 31557600.            #julian year in seconds
        self.au = 149598.0e8             #in cm
        self.c = 2.99792458e10          #in cm/s
        self.h = 6.62606957e-27         #in erg*s, Planck constant
        self.k = 1.3806488e-16          #in erg/K, Boltzmann constant
        self.sigma = 5.67040040e-5         #in erg/cm^2/s/K^4   Harmanec & Psra 2011
        self.pc = 3.08568025e16         #in cm
        self.mh = 1.672661e-24           #in g, mass hydrogen atom
        self.G = 6.674e-8               # in cm^3 g^-1 s^-2
        
        self.path_gastronoom = path_gastronoom        
        self.path_mcmax = path_mcmax
        
        #-- Convenience paths
        cc.path.mout = os.path.join(cc.path.mcmax,self.path_mcmax)
        cc.path.gout = os.path.join(cc.path.gastronoom,self.path_gastronoom)
        
        self.dust_list = None
        
        

    def __getitem__(self,key):

        """
        Overriding the standard dictionary __getitem__ method.
        
        @param key: Star()[key] where key is a string for which a corresponding
                    dictionary value is searched. If the key is not present in 
                    the dictionary, an attempt is made to calculate it from 
                    already present data; if it fails a KeyError is still 
                    raised. 
        @type key: string            
        
        @return: The value from the Star() dict for key
        @rtype: any
        
        """
        
        if not self.has_key(key):
            self.missingInput(key)
            return super(Star,self).__getitem__(key)
        elif super(Star,self).__getitem__(key) == '%':
            del self[key]
            self.missingInput(key)
            value = super(Star,self).__getitem__(key)
            self[key] = '%'
            return value 
        else:
            return super(Star,self).__getitem__(key)



    def __cmp__(self,star):
        
        """
        Overriding the standard dictionary __cmp__ method.
        
        A parameter set (dictionary of any type) is compared with this instance
        of Star(). 
        
        An attempt is made to create keys with values in each dict, if the 
        other has keys that are not present in the first. If this fails, False
        is returned.     
        
        @param star: A different parameter set. 
        @type star: dict or Star()             
        
        @return: The comparison between this object and star
        @rtype: bool
        
        """
        
        try:
            all_keys = set(self.keys() + star.keys())
            for k in all_keys:
                if not self.has_key(): 
                    self[k]
                if not star.has_key():
                    star[k]
            #if len(self) > len(star):
                #changed_keys = [star[k] 
                                #for k in self.keys() 
                                #if not star.has_key(k)]
                #print "Initialized keys in STAR2 from STAR1 == STAR2 comparison : \n" + str(changed_keys)
            #if len(self) < len(star):
                #changed_keys = [self[k] for k in star.keys() if not self.has_key(k)]
                #print "Initialized keys in STAR1 from STAR1 == STAR2 comparison : \n" + str(changed_keys)
        except KeyError:
            print 'Comparison error: Either STAR1 or STAR2 contains a key ' + \
                  'that cannot be initialized for the other.'
            print 'Both STAR instances are considered to be unequal.'
        finally:
            if isinstance(star, super(Star)):
                return cmp(super(Star,self), super(Star,star))
            else:
                return cmp(super(Star,self), star)                 

                      
                            
    def readDustProperties(self):
        
        '''
        When requesting the list of dust species, all dust species with nonzero
        abundance are returned. The list is created here, and remembered. 
        
        When this is done, the dust properties are also read for those species.
        
        The info is saved in self.dust.
        
        The dust list itself can be accessed through getDustList(). This method
        also calls the readDustProperties method to ensure all info is always 
        available. self.dust_list has a fixed order of appearance, according to
        Dust.dat. Important for, e.g., MCMax.py. 
        
        '''
        
        #-- If already read, don't read again
        if self.dust_list <> None: return
    
        #-- Read the info
        dust_info = dict()
        for k,s in zip(['SPECIES_SHORT','SPEC_DENS','MOLAR_WEIGHT','T_DES',\
                        'T_DESA','T_DESB','PART_FILE'],\
                        ['species','sd','molar','tdes','tdesa','tdesb','fn']):
            dust_info[s] = DataIO.getInputData(keyword=k,filename='Dust.dat')
        
        #-- Check which dust species are requested and exist in Dust.dat
        dust_list = [species 
                     for species in dust_info['species']
                     if self.has_key('A_' + species)]
        dust_list = [species 
                     for species in dust_list
                     if float(self['A_' + species]) != 0]
        
        #-- order in which the species appear is fixed according to the 
        #   order of the species in the Dust.dat input file. Note that this
        #   order does not matter for the database. 
        #   rgrain species are put first, following the dust_list order.
        #   Then the rest. This is the same as what is used in MCMax.py.
        dust_list = sorted(dust_list,\
                           key=lambda x: (not self.has_key('RGRAIN_%s'%x),\
                                          dust_list.index(x)))
        self.dust_list = tuple(dust_list)
        
        print '=========='
        print 'Dust species taken into account during modeling: %s'\
              %(', '.join(self.dust_list))
          
        #-- Save all relevant info based on which species are requested here
        self.dust = dict()
        for species in self.dust_list:
            index = dust_info['species'].index(species)
            self.dust[species] = dict()
            for key in ['sd','molar','tdes','tdesa','tdesb','fn']:
                self.dust[species][key] = dust_info[key][index]        
                
        
    
    def getDustList(self):
        
        '''
        Return (and initialize) the list of nonzero abundance dust species in 
        the model. 
        
        @return: The dust species
        @rtype: tuple(str)
        
        '''
        
        if self.dust_list is None: self.readDustProperties()
        return self.dust_list
    
    
    
    def addCoolingPars(self):
        
        '''
        Add Star parameters from the cooling database. 
        
        Any existing parameters are overwritten!
        
        '''
        
        cooling_path = os.path.join(cc.path.gout,\
                                    'GASTRoNOoM_cooling_models.db')
        cool_db = Database.Database(cooling_path)
        cooling_dict = cool_db[self['LAST_GASTRONOOM_MODEL']].copy()
        cooling_keys = ['T_STAR','R_STAR','TEMDUST_FILENAME','MDOT_GAS']
        for k in cooling_dict.keys():
            if k not in cooling_keys: del cooling_dict[k]            
        cooling_dict['R_STAR'] = float(cooling_dict['R_STAR'])/self.Rsun
        self.update(cooling_dict)



    def writeDensity(self):
        
        '''
        Write dust mass density and n(h2) profile (in Rstar).
        
        Only if MCMax or GASTRoNOoM model_id is available! 
        
        '''
        
        mcmid = self['LAST_MCMAX_MODEL']
        gasid = self['LAST_GASTRONOOM_MODEL']
        if mcmid:    
            rad = self.getDustRad(unit='rstar')
            dens = self.getDustDensity()
            fn = os.path.join(cc.path.mout,'models',mcmid,\
                              'density_profile_%s.dat'%mcmid)
            DataIO.writeCols(fn,[rad,dens])
        if gasid:
            rad = self.getGasRad(unit='rstar')
            nh2 = self.getGasNumberDensity()
            fn = os.path.join(cc.path.gout,'models',gasid,\
                              'nh2_density_profile_%s.dat'%gasid)
            DataIO.writeCols(fn,[rad,nh2])


    def readKappas(self):
        
        '''
        Read the kappas.dat file of an MCMax model.
    
        '''
        
        opas = DataIO.readCols(os.path.join(cc.path.mout,'models',\
                                            self['LAST_MCMAX_MODEL'],\
                                            'kappas.dat'))
        return opas.pop(0),opas
                            

    #OBSOLETE. NOT USED. PERHAPS IN THE FUTURE. DONT CALL THIS UPON INITIALIZATION OF STAR()
    #def convertRadialUnit(self):
        
        #'''
        #Convert any radial unit for shell parameters to R_STAR.
        
        #'''
        
        #if self['R_SHELL_UNIT'] != 'R_STAR':
            #shell_units = ['CM','M','KM','AU']
            #unit_index = shell_units.index(self['R_SHELL_UNIT'].upper())
            #unit_conversion = [1./(self.Rsun*self['R_STAR']),\
                               #10.**2/(self.Rsun*self['R_STAR']),\
                               #10.**5/(self.Rsun*self['R_STAR']),\
                               #self.au/(self.Rsun*self['R_STAR'])]
            #for par in ['R_INNER_GAS','R_INNER_DUST','R_OUTER_GAS',\
                        #'R_OUTER_DUST'] \
                    #+ ['R_MAX_' + species for species in self.getDustList()] \
                    #+ ['R_MIN_' + species for species in self.getDustList()]:
                #if self.has_key(par):
                    #self[par] = self[par]*unit_conversion[unit_index]
        #else:
            #pass
                                             

    
    def removeMutableMCMax(self,mutable,var_pars):
        
        """
        Remove mutable parameters after an MCMax run.
    
        @param mutable: mutable keywords
        @type mutable: list[string]
        @param var_pars: parameters that are varied during gridding, these will
                         not be removed and are kept constant throughout the 
                         iteration
        @type var_pars: list[string]
        
        """
        
        #- remove keys which should be changed by output of new mcmax model, 
        #- but only the mutable input!!! 
        for key in self.keys():
            if key in mutable \
                      + ['R_MAX_' + species for species in self.getDustList()]\
                      + ['T_DES_' + species for species in self.getDustList()]\
                      + ['R_DES_' + species for species in self.getDustList()]\
                    and key not in var_pars:
                del self[key]
        
        #- Check the effective destruction temperature of every species, and 
        #- see if max and min T's are as requested.
        self.checkT()
        
        #- No point in keeping zeroes around for T_DES or T_MIN
        for species in self.getDustList():
            for par in ('T_DES_' + species, 'T_MIN_' + species):
                if self.has_key(par):
                    if not float(self[par]):
                        del self[par]
                    
  

    def removeMutableGastronoom(self,mutable,var_pars):
        
        """
        Remove mutable parameters after a GASTRoNOoM run.
    
        @param mutable: mutable parameters
        @type mutable: list[string]
        @param var_pars: parameters that are varied during gridding, these will
                         not be removed and are kept constant throughout the 
                         iteration
        @type var_pars: list[string]
        
        """
        
        for key in self.keys():
            if key in mutable and key not in var_pars:
                del self[key]
                    
 
    
    def updateMolecules(self,parlist):
        
        '''
        Update variable information in the molecule instances of this star.
        
        @param parlist: parameters that have to be updated.
        @type parlist: list[string]
        
        '''
        
        for molec in self['GAS_LIST']:
            molec.updateParameters(pardict=dict([(k,self[k]) 
                                                 for k in parlist]))
    


    def normalizeDustAbundances(self):
        
        """
        Normalize the dust abundances such that they add up to a total of 1.
        
        If the MRN_DUST keyword for MCMax is nonzero, all nonzero abundances 
        are set to 1. The abundance given in the inputfile does not matter in 
        this case.
        
        """
        
        abun_ori = [self['A_%s'%sp] for sp in self.getDustList()]
        self['A_DUST_ORIGINAL'] = abun_ori
        if int(self['MRN_DUST']): 
            self['A_NO_NORM'] = 1
            print 'WARNING! Take care when extracting output in MCMax using '+\
                  'these scripts, if MRN_DUST == 1! Especially if some ' + \
                  'abundances are set manually and some according to MRN: ' + \
                  'these are not normalized to 1, since this depends on the '+\
                  'MRN distributed dust species.'
            for sp in self.getDustList():
                mrn_count = 0
                if self['MRN_NGRAINS'] != len(self.getDustList()) \
                        and self.has_key('RGRAIN_%s'%sp):
                    self.__setitem__('A_%s'%sp,2)
                    mrn_count += 1
                elif self['MRN_NGRAINS'] == len(self.getDustList()):
                    self.__setitem__('A_%s'%sp,2)        
                    mrn_count += 1
                if mrn_count != self['MRN_NGRAINS']:
                    raise IOError('MRN_NGRAINS not equal to amount of RGRAIN_sp keywords.')
        total = sum(abun_ori)
        if not int(self['A_NO_NORM']) and '%.3f'%total != '1.000':
            print 'Normalizing dust abundances to 1, from a total of %f.'%total
            abun_new = [a/total for a in abun_ori]
            print ', '.join(['%.2f'%a for a in abun_ori]), ' is changed to ', \
                  ', '.join(['%.2f'%a for a in abun_new]), ' for ', \
                  ', '.join(self.getDustList()), '.'
            [self.__setitem__('A_%s'%sp,a) for a,sp in zip(abun_new,\
                                                           self.getDustList())]



    def calcA_NO_NORM(self):
        
        """
        Set the default value of A_NO_NORM to 0.
        
        """
        
        if not self.has_key('A_NO_NORM'):
            self['A_NO_NORM'] = 0
        else:
            pass
    
        

    def __addLineList(self):
        
        """ 
        Take molecular transitions from a line list and wavelength range.
        
        Based on the GASTRoNOoM radiat and indices data files. See Molecule.py
        for more info.
        
        """
        
        gas_list = []
        if type(self['LS_TELESCOPE']) is types.StringType:
            self['LS_TELESCOPE'] = [self['LS_TELESCOPE']]
        if not self['LS_NO_VIB']:
            self['LS_NO_VIB'] = []
        elif type(self['LS_NO_VIB']) is types.StringType:
            self['LS_NO_VIB'] = [self['LS_NO_VIB']]
        ctrl = [tr.getInputString(include_nquad=0) for tr in self['GAS_LINES']]
        for molec in self['GAS_LIST']:
            for telescope in self['LS_TELESCOPE']:
                if telescope == 'PACS':
                    ls_min = 50
                    ls_max = 200
                elif telescope == 'SPIRE':
                    ls_min = 180
                    ls_max = 700
                nl = Transition.makeTransitionsFromRadiat(molec=molec,\
                            telescope=telescope,ls_min=ls_min,ls_max=ls_max,\
                            n_quad=self['N_QUAD'],offset=self['LS_OFFSET'],\
                            use_maser_in_sphinx=self['USE_MASER_IN_SPHINX'],\
                            path_gastronoom=self.path_gastronoom,\
                            ls_unit='micron',\
                            no_vib=molec.molecule in self['LS_NO_VIB'])
                #-- exclude transitions if they are already included. 
                #   This is to avoid adding a double with a different n_quad, in 
                #   case it was included manually. Can use GAS_LINES key, as it 
                #   contains both manual and other lines, but the latter are 
                #   also assigned self['N_QUAD'] anyway
                nl = [tr for tr in nl
                         if tr.getInputString(include_nquad=0) not in ctrl]
                gas_list.extend(nl)
        self['GAS_LINES'].extend(gas_list)


    
    def calcLS_NO_VIB(self):
        
        """
        Set the default value of LS_NO_VIB (remove vibrational states from the
        calculation and plots) to zero.        
        
        """
        
        if not self.has_key('LS_NO_VIB'):
            self['LS_NO_VIB'] = []
        else:
            pass
    
    
        
    def calcN_QUAD(self):
        
        """
        Set the default value of N_QUAD to 100. 
        
        Only used when auto selecting transition based on a wavelength range.
        
        """
        
        if not self.has_key('N_QUAD'):
            self['N_QUAD'] = 100
        else:
            pass
    
    
    
    def calcLS_OFFSET(self):
        
        """
        Set the default value of LS_OFFSET to 0.0. 
        
        Only used when auto selecting transitions based on a wavelength range.
        
        """
        
        if not self.has_key('LS_OFFSET'):
            self['LS_OFFSET'] = 0.0
        else:
            pass
    
        

    def checkT(self):
        
        """
        Search input list for minimum temperature.
    
        Method prints the actual minimum T for which the model was calculated.
        
        Note that the density drops to zero gradually and that the criterium
        has to be sudden change of slope. Check criterium if the printed T is 
        not good enough as determination of max radius IS correct.
        
        """
        
        coldens = ColumnDensity.ColumnDensity(self)
        self.calcT_INNER_DUST()
        for index,species in enumerate(self.getDustList()):
            self['T_DES_%s'%species] = coldens.t_des[species]
            self['R_DES_%s'%species] = coldens.r_des[species]\
                                        /self.Rsun/self['R_STAR']
            print 'The EFFECTIVE maximum temperature for species %s '%species+\
                  'is %.2f K, at radius %.2f R_STAR.'\
                  %(self['T_DES_%s'%species],self['R_DES_%s'%species])
        
        species_list_min = [species 
                            for species in self.getDustList() 
                            if self.has_key('T_MIN_%s'%species) \
                                or self.has_key('R_MAX_%s'%species)]
        for species in species_list_min:
            print 'The EFFECTIVE minimum temperature for species'+\
                  ' %s is %.2f K at maximum radius %.2f R_STAR.'\
                  %(species,coldens.t_min[species],\
                    coldens.r_max[species]/self.Rsun/self['R_STAR'])
            if self.has_key('T_MIN_%s'%species):
                print 'The REQUESTED minimum temperature for species '+\
                      '%s is %.2f K.'%(species,self['T_MIN_%s'%species])
            if self.has_key('R_MAX_%s'%species):
                print 'The REQUESTED maximum radius for species'+\
                      '%s is %.2f R_STAR.'%(species,self['R_MAX_%s'%species])
            print 'The EFFECTIVE outer radius of the shell is %.2f R_STAR.'\
                  %(coldens.r_outer/self.Rsun/self['R_STAR'])
            print 'Note that if R_MAX is ~ the effective outer radius, the ' +\
                  'requested minimum temperature may not be reached.'
        return
    

    
    def getFullNameMolecule(self,short_name):
        
        '''
        Get the full name of a molecule, based on it's short name, 
        if it is present in the GAS_LIST.
                
        @return: Name the name. None if not available.
        @rtype: string
        
        '''
        
        molecule = [molec.molecule_full 
                    for molec in self['GAS_LIST'] 
                    if molec.molecule == short_name]
        if not molecule:
            return None
        #- Return first only, if multiple there's multiple requested molecules 
        #- of the same type (fi different abundance)
        else:
            return molecule[0]     



    def getShortNameMolecule(self,full_name):
        
        '''
        Get the short name of a molecule, based on it's full name, 
        if it is present in the GAS_LIST.
        
        @param full_name: The full name of a molecule from Molecule.dat
        @type full_name: string
        
        @return: None if not available, otherwise the short hand name.
        @rtype: string
        
        '''
        
        molecule = [molec.molecule 
                    for molec in self['GAS_LIST'] 
                    if molec.molecule_full == full_name]
        if not molecule:
            return None
        #- Return first only, if multiple there's multiple requested molecules 
        #- of the same type (fi different abundance)
        else:
            return molecule[0]     
      

     
    def getMolecule(self,molec_name):
        
        '''
        Get a Molecule() object based on the molecule name. 
        
        A Star() object always has only one version of one molecule.
        
        @param molec_name: short name of the molecule
        @type molec_name: string
        
        @return: The molecule
        @rtype: Molecule()
        
        '''
        
        try:
            return [molec 
                    for molec in self['GAS_LIST'] 
                    if molec.molecule == molec_name][0]
        except IndexError:
            return None
    
    
    
    def getTransition(self,sample):
        
        '''
        Return a Transition() object that has the same parameters as sample. 
        
        The actual model ids are not included in this comparison! 
        
        None is returned if no match is found. 
        
        @param sample: A sample transition to be cross referenced with the 
                       transitions in this Star() object. If a match is found, 
                       it is returned.
        @type sample: Transition()
        @return: If a match is found, this transition is returned.
        @rtype: Transition()
        
        '''
         
        i = 0
        while i < len(self['GAS_LINES']) and sample != self['GAS_LINES'][i]:
            i += 1
        if i == len(self['GAS_LINES']):
            return None
        else:
            return self['GAS_LINES'][i]
        
    
     
    def getTransList(self,**kwargs):
        
        '''
        Return a list of (transmodelid, molecmodelid, dictionary) for every 
        transition in the Star model.
        
        '''
        
        trl = Transition.extractTransFromStars([self],**kwargs)
        trl_info = [(trans.getModelId(),\
                     trans.molecule.getModelId(),\
                     trans.makeDict())
                    for trans in trl]
        return trl_info



    def getTransitions(self,molec):
        
        '''
        Return a list of all transitions associated with a single molecule.
        
        @param molec: the shorthand notation of the molecule
        @type molec: string
        
        @return: All transitions for this molecule
        @rtype: list[Transition()]
        
        '''
        
        return [trans 
                for trans in self['GAS_LINES'] 
                if trans.molecule.molecule==molec]



    def getDustFn(self,species=''):
        
        '''
        Return the dust output filename for density and temperature.
        
        You can choose the species for which the output file is retrieved. 
        Unless thermal contact for this model is on. 
        
        @keyword ftype: The species for which dust info is read. Default if the
                        average profiles are requested. Can be any species in 
                        Dust.dat as long as the model calculation included it.
                      
                        (default: '')
        @type ftype: str
        
        @return: The filename of the requested dust output file.
        @rtype: str
        
        '''
        
        fp = os.path.join(cc.path.mout,'models',self['LAST_MCMAX_MODEL'])
        
        #- if T_CONTACT: no specific species denstemp files, 
        #- so denstemp.dat is taken
        if species and not self['T_CONTACT']:
            ispecies = self.getDustList().index(species)
            fn = os.path.join(fp,'denstempP%.2i.dat'%(ispecies+1))
        else:
            fn = os.path.join(fp,'denstemp.dat')
        
        return fn
        
        
        
    def getDustRad(self,species='',unit='cm'):
        
        '''
        Return the dust radial grid in cm, au or Rstar. 
              
        @keyword species: The dust species for which to return the grid. If 
                          default or if T_CONTACT is on, the file is denstemp
                          otherwise it is the species specific file.
        
                          (default: '')
        @type species: str
        @keyword unit: The unit of the returned grid. Can be 'cm','rstar','au',
                       'm'
        
                       (default: 'cm')
        @type unit: str
        
        @return: array giving the radial grid.
        @rtype: array
        
        '''

        if not self['LAST_MCMAX_MODEL']: return empty(0)
    
        fn = self.getDustFn(species)
        rad = array(DataIO.getMCMaxOutput(incr=int(self['NRAD']),\
                                          keyword='RADIUS',filename=fn))
        
        unit = str(unit).lower()
        if unit == 'au':
            rad = rad/self.au
        elif unit == 'rstar':
            rad = rad/self.Rsun/self['R_STAR']
        elif unit == 'm':
            rad = rad*10**-2
            
        return rad   
    
    
    
    def getDustDensity(self,species=''):
        
        '''
        Return the dust density profile from the file made for MCMax.
        
        self.getDustRad returns the radial grid associated with this profile. 
        
        An empty array is returned if the model does not exist.
        
        @keyword species: The dust species for which to return the grid. If 
                          default or if T_CONTACT is on, the file is denstemp
                          otherwise it is the species specific file.
        
                          (default: '')
        @type species: str
        @return: the density profile (g/cm3)
        @rtype: array
        
        '''
        
        if not self['LAST_MCMAX_MODEL']: return empty(0)
        
        #-- Read the dust density profile and reduce the array by averaging
        #   over the azimuthal coordinate.
        fn = self.getDustFn(species)
        incr = int(self['NRAD'])*int(self['NTHETA'])
        dens_ori = DataIO.getMCMaxOutput(filename=fn,incr=incr,\
                                         keyword='DENSITY')
        dens = Data.reduceArray(dens_ori,self['NTHETA'])
        
        return dens         
         
         

    def getDustTemperature(self,add_key=0,species=''):
         
        '''
        Return the dust temperature profile from MCMax model.
        
        self.getDustRad returns the radial grid associated with this profile. 
        
        An empty array is returned if the model does not exist.
        
        The total dust temperature without separate components for the 
        different dust species is returned if no species is requested or if 
        thermal contact is on.
        
        @keyword species: The dust species for which to return the grid. If 
                          default or if T_CONTACT is on, the file is denstemp
                          otherwise it is the species specific file.
        
                          (default: '')
        @type species: str
        @keyword add_key: Add a key for a legend to the ouput as third tuple
                          element.
                          
                          (default: 0)
        @type add_key: bool
        
        @return: The temperature profile (K) as well as a key, if requested.
        @rtype: (array,string)
        
        '''
        
        if not self['LAST_MCMAX_MODEL']: 
            return add_key and (empty(0),'') or empty(0)

        #-- Read the dust temperature profile and reduce the array by averaging
        #   over the azimuthal coordinate.
        fn = self.getDustFn(species)
        incr = int(self['NRAD'])*int(self['NTHETA'])
        temp_ori = DataIO.getMCMaxOutput(incr=incr,keyword='TEMPERATURE',\
                                         filename=fn)
        temp = Data.reduceArray(temp_ori,self['NTHETA'])

        if add_key:
            key = '$T_{\mathrm{d, avg}}$ for %s'\
                  %self['LAST_MCMAX_MODEL'].replace('_','\_')
            return temp,key
        else: 
            return temp
        
    
    
    def getGasNumberDensity(self,molecule=False,**kwargs):
        
        '''
        Give the n(h2) number density profile of the gas read from a GASTRoNOoM
        model.
        
        Additional input keywords for self.getCoolFn() can be passed along.
        
        An empty is returned in case no model is available. 
        
        @keyword molecule: Return the molecular number density instead of the
                           H_2 number density.
                           
                           (default: False)
        @type molecule: bool 
        
        @return: The number density (in cm-3) profile of n(H_2)
        @rtype: array
        
        '''
        
        if not self['LAST_GASTRONOOM_MODEL']: return empty(0)
        fgr_file = self.getCoolFn(**kwargs)
        
        kws = dict()
        if molecule:
            kws['keyword'] = 'N(MOLEC)'
            kws['key_index'] = 8
        else:
            kws['keyword'] = 'N(H2)'
        
        nmol = DataIO.getGastronoomOutput(filename=fgr_file,return_array=1,\
                                          **kws)
                                              
        return nmol
    
    
    
    
    
    def getGasVelocity(self,**kwargs):
        
        '''
        Give the velocity profile of the gas read from a GASTRoNOoM model.
        
        Additional input keywords for self.getCoolFn() can be passed along.
        
        An empty is returned in case no model is available. 
        
        @return: The velocity (in cm/s) profile
        @rtype: array
        
        '''
        
        if not self['LAST_GASTRONOOM_MODEL']: return empty(0)
        fgr_file = self.getCoolFn(**kwargs)
        vel = DataIO.getGastronoomOutput(filename=fgr_file,keyword='VEL',\
                                         return_array=1)
        return vel
        
        
    
    def getGasTemperature(self,**kwargs):
        
        '''
        Give the temperature profile of the gas read from a GASTRoNOoM model.
        
        Additional input keywords for self.getCoolFn() can be passed along.
        
        An empty is returned in case no model is available. 
        
        @return: The temperature (in K) profile
        @rtype: array
        
        '''
        
        if not self['LAST_GASTRONOOM_MODEL']: return empty(0)
        fgr_file = self.getCoolFn(**kwargs)
        temp = DataIO.getGastronoomOutput(filename=fgr_file,keyword='TEMP',\
                                          return_array=1)
        return temp
    


    def getGasRad(self,unit='cm',ftype='fgr_all',**kwargs):
        
        '''
        Return the gas radial grid in cm, AU or Rstar. 
        
        Additional input keywords for self.getCoolFn() can be passed along.
        
        An empty is returned in case no model is available. 
        
        @keyword unit: The unit of the returned grid. Can be 'cm','rstar','au',
                       'm'
        
                       (default: 'cm')
        @type unit: str
        @keyword ftype: The cooling output file type. Either '1', '2', '3', 
                        'fgr', 'fgr_all', or 'rate'.
                      
                        (default: 'fgr_all')
        @type ftype: str
        
        @return: the radial grid
        @rtype: array
        
        '''
        
        if not self['LAST_GASTRONOOM_MODEL']: return empty(0)
        
        #-- coolrate***.dat does not give radii.
        if ftype == 'rate': return empty(0)
        
        unit = str(unit).lower()
        fgr_file = self.getCoolFn(ftype=ftype,**kwargs)
        rad = DataIO.getGastronoomOutput(filename=fgr_file,keyword='RADIUS',\
                                         return_array=1)
        #-- fgr_all gives radius in cm. Others in rstar. Convert others to cm
        if ftype != 'fgr_all':
            rad = rad*self['R_STAR']*self.Rsun
        
        if unit == 'au':
            rad = rad/self.au
        elif unit == 'rstar':
            rad = rad/self.Rsun/self['R_STAR']
        elif unit == 'm':
            rad = rad*10**-2
        return rad
        
        
    
    def getCoolFn(self,ftype='fgr_all',mstr='',modelid=''):
        
        '''
        Return the cooling output filename.
        
        You can define the type of cooling file you want, as well as an 
        additional identification string for the molecule/sampling.
        
        @keyword ftype: The cooling output file type. Either '1', '2', '3', 
                        'fgr', 'fgr_all', or 'rate'.
                      
                        (default: 'fgr_all')
        @type ftype: str
        @keyword mstr: The additional identication string. Not applicable to 
                       'fgr' or 'fgr_all'. Can be any molecule, or 'sampling'.
                       File must exist to be used further!
                       
                       (default: '')
        @type mstr: str
        @keyword modelid: The model id to be used. If default, the cooling id 
                          is used. If defined, it could refer to an id 
                          different from the cooling id such as an mline id.
                          
                          (default: '')
        @type modelid: str
        
        @return: The filename of the requested cooling output file.
        @rtype: str
        
        '''
        
        if not modelid:
            modelid = self['LAST_GASTRONOOM_MODEL']
        ftype = str(ftype)
        if ftype == 'fgr' or ftype == 'fgr_all':
            mstr = ''
        if mstr:
            mstr = '_' + mstr
        fn = os.path.join(cc.path.gout,'models',modelid,\
                          'cool%s%s%s.dat'%(ftype,modelid,mstr))
        return fn
        
        

    def calcTLR(self):  
        
        """
        Stefan-Boltzmann's law.
            
        Star() object needs to have at least 2 out of 3 parameters (T,L,R), 
        with in L and R in solar values and T in K.
    
        The one missing parameter is calculated. 
    
        This method does nothing if all three are present.
        
        """
        
        if not self.has_key('T_STAR'):
            self['T_STAR']=(float(self['L_STAR'])/float(self['R_STAR'])**2.)\
                                **(1/4.)*self.Tsun
        elif not self.has_key('L_STAR'):
            self['L_STAR']=(float(self['R_STAR']))**2.*\
                                (float(self['T_STAR'])/self.Tsun)**4.
        elif not self.has_key('R_STAR'):
            self['R_STAR']=(float(self['L_STAR'])*\
                                (self.Tsun/float(self['T_STAR']))**4)**(1/2.)
        else:
            pass 



    def calcLL_GAS_LIST(self):
        
        '''
        Define Molecule() objects for the molecules requested in the LineList
        mode.
        
        '''
        
        if not self.has_key('LL_GAS_LIST'):
            if type(self['LL_MOLECULES']) is types.StringType:
                self['LL_GAS_LIST'] = [Molecule.Molecule(linelist=1,\
                                                molecule=self['LL_MOLECULES'])]
            else:
                self['LL_GAS_LIST'] = [Molecule.Molecule(molecule=molec,
                                                         linelist=1) 
                                       for molec in self['LL_MOLECULES']]
        else:
            pass
        
   
    
    def calcUSE_MASER_IN_SPHINX(self):
        
        '''
        Set the default value of USE_MASER_IN_SPHINX parameter.
        
        '''
        
        if not self.has_key('USE_MASER_IN_SPHINX'):
            self['USE_MASER_IN_SPHINX']=0
        else:
            pass


    
    def calcLOGG(self):
        
        """
        Set the default value of LOGG to 0.
        
        """
        
        if not self.has_key('LOGG'):
            self['LOGG']=0
        else:
            pass


    def calcF_CONT_TYPE(self):
        
        """
        Set the default value of F_CONT_TYPE to MCMax. This is the type of 
        derivation of measured continuum fluxes. Can be: ISO, MSX, PHOT, MCMax
        
        """
        
        if not self.has_key('F_CONT_TYPE'):
            self['F_CONT_TYPE'] = 'MCMax'
        else:
            pass
    
    
    def calcAH2O_RATE(self):
        
        '''
        Calculate the outflow rate of H2O, by multiplying the H2O abundance 
        with the mass-loss rate. 
        
        Value is set in units of Msun/yr
        
        '''
        
        if not self.has_key('AH2O_RATE'):
            self['AH2O_RATE'] = self['F_H2O'] * self['MDOT_GAS']
        else:
            pass
    


    def calcAH2O_RATE(self):
        
        '''
        Calculate the outflow rate of H2O, by multiplying the H2O abundance 
        with the mass-loss rate. 
        
        Value is set in units of Msun/yr
        
        '''
        
        if not self.has_key('AH2O_RATE'):
            self['AH2O_RATE'] = self['F_H2O'] * self['MDOT_GAS']
        else:
            pass


        
    def calcT_INNER_DUST(self):
        
        """
        Find the dust temperature at the inner radius in Kelvin.
        
        Taken from last mcmax model, and defined by the dust species able to 
        exist at the highest temperature; if no mcmax model is present, the 
        temperature is taken to be zero, indicating no inner radius T is 
        available.
        
        """        
        
        if not self.has_key('T_INNER_DUST'):
            rad = self.getDustRad(unit='rstar')
            temp = self.getDustTemperature()
            if not rad.size:
                self['T_INNER_DUST'] = 0
                return
            rin = self['R_INNER_DUST']
            self['T_INNER_DUST'] = temp[argmin(abs(rad-rin))]
        else:
            pass
        


    def calcTEMPERATURE_EPSILON_GAS(self):
        
        """
        If not present in input, TEMPERATURE_EPSILON_GAS is equal to 0.5.
        
        """
        
        if not self.has_key('TEMPERATURE_EPSILON_GAS'):
            self['TEMPERATURE_EPSILON_GAS'] = 0.5
        else:
            pass



    def calcTEMPERATURE_EPSILON2_GAS(self):
        
        """
        If not present in input, TEMPERATURE_EPSILON2_GAS is equal to 0, 
        in which case it will be ignored when making input file.
        
        """
        
        if not self.has_key('TEMPERATURE_EPSILON2_GAS'):
            self['TEMPERATURE_EPSILON2_GAS'] = 0
        else:
            pass



    def calcRADIUS_EPSILON2_GAS(self):
        
        """
        If not present in input, RADIUS_EPSILON2_GAS is equal to 0, \
        in which case it will be ignored when making input file.
        
        """
        
        if not self.has_key('RADIUS_EPSILON2_GAS'):
            self['RADIUS_EPSILON2_GAS'] = 0
        else:
            pass



    def calcTEMPERATURE_EPSILON3_GAS(self):
        
        """
        If not present in input, TEMPERATURE_EPSILON3_GAS is equal to 0, 
        in which case it will be ignored when making input file.
        
        """
        
        if not self.has_key('TEMPERATURE_EPSILON3_GAS'):
            self['TEMPERATURE_EPSILON3_GAS'] = 0
        else:
            pass



    def calcRADIUS_EPSILON3_GAS(self):
        
        """
        If not present in input, RADIUS_EPSILON3_GAS is equal to 0, 
        in which case it will be ignored when making input file.
        
        """
        
        if not self.has_key('RADIUS_EPSILON3_GAS'):
            self['RADIUS_EPSILON3_GAS'] = 0
        else:
            pass



    def calcDUST_TO_GAS_CHANGE_ML_SP(self):
        
        """
        Set default value of sphinx/mline specific d2g ratio to the 
        semi-empirical d2g ratio, ie based on MDOT_DUST and MDOT_GAS. 
        
        In order to turn this off, set this parameter to 0 in the input file, 
        in which case the iterated acceleration d2g ratio is used.
        
        Both MDOT_GAS and MDOT_DUST have to be defined explicitly if this 
        parameter is not. 
        
        This parameter has to be defined explicitly if one of MDOT_GAS and 
        MDOT_DUST is not defined explicitly. 
        
        Note that the DUST_TO_GAS keyword is the internal representation of the
        dust_to_gas ratio and should never be explicitly defined. For all 
        practical purposes, use DUST_TO_GAS_CHANGE_ML_SP.
        
        """
        
        if not self.has_key('DUST_TO_GAS_CHANGE_ML_SP'):
            if not self.has_key('MDOT_DUST'):
                raise IOError('Both MDOT_DUST and DUST_TO_GAS_CHANGE_ML_SP '+\
                              'are undefined.')
            if not self.has_key('MDOT_GAS'):
                raise IOError('Both MDOT_GAS and DUST_TO_GAS_CHANGE_ML_SP '+\
                              'are undefined.')
            self['DUST_TO_GAS_CHANGE_ML_SP'] = self['DUST_TO_GAS']
        else:
            pass



    def calcR_INNER_GAS(self):
        
        """
        If not present in input, R_INNER_GAS is equal to R_INNER_DUST
    
        """
        
        if not self.has_key('R_INNER_GAS'):
            self['R_INNER_GAS'] = self['R_INNER_DUST']
        else:
            pass



    def calcUSE_DENSITY_NON_CONSISTENT(self):
        
        """
        Set USE_DENSITY_NON_CONSISTENT off by default.
        
        """
        
        if not self.has_key('USE_DENSITY_NON_CONSISTENT'):
            self['USE_DENSITY_NON_CONSISTENT'] = 0
        else:
            pass        



    def calcR_OUTER_DUST(self):
        
        """
        If not present in input, R_OUTER_DUST is calculated from 
        R_OUTER_DUST_AU.
            
        """
        
        if not self.has_key('R_OUTER_DUST'):
            if self.has_key('R_OUTER_DUST_AU'):
                self['R_OUTER_DUST'] = self['R_OUTER_DUST_AU']*self.au\
                                            /self['R_STAR']/self.Rsun
            elif self.has_key('R_OUTER_MULTIPLY'):
                self['R_OUTER_DUST'] = self['R_INNER_DUST']\
                                            *self['R_OUTER_MULTIPLY']
        else:
            pass
        


    def calcR_INNER_DUST(self):
    
        """
        Calculate the inner radius from MCMax output in stellar radii.
        
        If no MCMax model is calculated yet, R_{i,d} is the stellar radius.
        
        Else, the inner dust radius is taken where the density reaches a 
        threshold, defined by R_INNER_DUST_MODE:
        
            - MAX: Density reaches a maximum value, depends on the different 
                   condensation temperatures of the dust species taken into 
                   account 
            - ABSOLUTE: Density becomes larger than 10**(-30)
            - RELATIVE: Density becomes larger than 1% of maximum density
        
        Unless defined in the CC input file, the dust radius is updated every 
        time a new iteration starts.
        
        If no MCMax model is known, and destruction temperature iteration is 
        off, the inner radius is 2 stellar radii for calculation time reasons.
        
        """
    
        if not self.has_key('R_INNER_DUST'):
            if self.has_key('R_INNER_DUST_AU'):
                self['R_INNER_DUST'] = self['R_INNER_DUST_AU']*self.au\
                                                /self['R_STAR']/self.Rsun
            else:
                rad = self.getDustRad()
                if rad.size == 0: 
                    self['R_INNER_DUST'] = 1.0
                    return
                dens = self.getDustDensity()
                if self['R_INNER_DUST_MODE'] == 'MAX':
                    ri_cm = rad[argmax(dens)]
                elif self['R_INNER_DUST_MODE'] == 'ABSOLUTE':
                    ri_cm = rad[dens>10**(-30)][0]
                else:
                    ri_cm = rad[dens>0.01*max(dens)][0]
                self['R_INNER_DUST'] = ri_cm/self.Rsun/self['R_STAR']
        else:
            pass



    def calcR_INNER_DUST_MODE(self):
         
        """
        The mode of calculating the inner radius from MCMax output, if present.
        
        Can be ABSOLUTE (dens>10**-50) or RELATIVE (dens[i]>dens[i+1]*0.01).
        
        CLASSIC reproduces the old method. 
        
        Set here to the default value of ABSOLUTE.
        
        """
        
        if not self.has_key('R_INNER_DUST_MODE'):
            self['R_INNER_DUST_MODE'] = 'ABSOLUTE'
        else:
            pass        



    def calcRID_TEST(self):
         
        """
        The mode of determining the dust temp profile. 
        
        Only for testing purposes.
        
            - Default is 'R_STAR', ie the temperature is taken from the stellar
              radius onward, regardless of what the inner radius is. 
        
            - 'R_INNER_GAS' is used for taking the dust temperature from the 
              inner gas radius onward. 
        
            - 'BUGGED_CASE' is the old version where r [R*] > R_STAR [Rsun]. 
        
        """
        
        if not self.has_key('RID_TEST'):
            self['RID_TEST'] = 'R_STAR'
        else:
            pass 



    def calcSPEC_DENS_DUST(self):
        
        """
        Calculating average specific density of dust shell.
        
        This is based on the input dust species abundances and their specific 
        densities.
    
        """
    
        if not self.has_key('SPEC_DENS_DUST'):
            if int(self['MRN_DUST']):
                these_sd = [self.dust[sp]['sd']
                            for sp in self.getDustList()
                            if self.has_key('RGRAIN_%s'%sp)]
                sd_mrn = sum(these_sd)/len(these_sd)
                a_sd = sum([self['A_' + sp]*self.dust[sp]['sd']
                            for sp in self.getDustList()
                            if self['A_%s'%sp] != 2.])
                self['SPEC_DENS_DUST'] = (sd_mrn + a_sd)/2.
            else:
                these_sd = [self['A_' + sp]*self.dust[sp]['sd']
                            for sp in self.getDustList()]
                self['SPEC_DENS_DUST'] = sum(these_sd)
        else:
            pass
            


    def calcLAST_MCMAX_MODEL(self):
        
        """
        Creates empty string if not present yet.
    
        """
        
        if not self.has_key('LAST_MCMAX_MODEL'):
            self['LAST_MCMAX_MODEL'] = ''
        else:
            pass
        


    def calcLAST_PACS_MODEL(self):
        
        """
        Sets to None if not present yet.
    
        """
        
        if not self.has_key('LAST_PACS_MODEL'):
            self['LAST_PACS_MODEL'] = None
        else:
            pass
        


    def calcLAST_SPIRE_MODEL(self):
        
        """
        Sets to None if not present yet.
        
        Note that this is an index if it IS present, and can be zero. Always 
        check with None instead of boolean. 
    
        """
        
        if not self.has_key('LAST_SPIRE_MODEL'):
            self['LAST_SPIRE_MODEL'] = None
        else:
            pass
        


    def calcLAST_GASTRONOOM_MODEL(self):
        
        """
        Creates empty string if not present yet.
    
        """
        
        if not self.has_key('LAST_GASTRONOOM_MODEL'):
            self['LAST_GASTRONOOM_MODEL'] = ''
        else:
            pass
    
    
    def calcDRIFT_TYPE(self):
        
        """
        Set the type of drift between dust and gas taken into account. 
        
        Is either consistent (from momentum transfer calculation or zero). 
        
        """
        
        if not self.has_key('DRIFT_TYPE'):
            if self['V_EXP_DUST'] == self['VEL_INFINITY_GAS']:
                self['DRIFT_TYPE'] = 'ZERO'
            else: 
                self['DRIFT_TYPE'] = 'CONSISTENT'
        else:
            pass
        
    
    
    def calcDRIFT(self):
        
        """
        Find terminal drift velocity from last calculated GASTRoNOoM model.
        
        Units are km/s and is given for grain size 0.25 micron.
        
        If no GASTRoNOoM model exists, the drift is taken to be 0.
        
        """
        
        if not self.has_key('DRIFT'):
            try:
                #- The last 10 entries should be roughly constant anyway, 
                #- ~terminal values
                self['DRIFT'] = self.getAverageDrift()[-5]/10.**5    
            except IOError:
                self['DRIFT'] = 0
                print 'No GASTRoNOoM model has been calculated yet, drift ' + \
                      'is unknown and set to the default of 0.'
        else:
            pass
        


    def calcM_DUST(self):
        
        """
        Find total dust mass, based on sigma_0 in the case of a power law.
    
        """
        
        if not self.has_key('M_DUST'):
            if self['DENSTYPE'] == 'POW':
                if self['DENSPOW'] == 2:
                    self['M_DUST']  \
                        = 2*pi*self['DENSSIGMA_0']\
                        *(self['R_INNER_DUST']*self['R_STAR']*self.Rsun)**2\
                        *log(self['R_OUTER_DUST']/float(self['R_INNER_DUST']))\
                        /self.Msun
                else:
                    self['M_DUST'] \
                        = 2*pi*self['DENSSIGMA_0']\
                        *(self['R_INNER_DUST']*self['R_STAR']*self.Rsun)**2\
                        /(2.-self['DENSPOW'])/self.Msun\
                        *((self['R_OUTER_DUST']/float(self['R_INNER_DUST']))\
                        **(2.-self['DENSPOW'])-1.)
            else:
                pass
        else:
            pass
        

    
    def calcMDOT_DUST(self): 
        
        '''
        Calculate the value of MDOT_DUST from the DUST_TO_GAS_RATIO_ML_SP. 
        
        Requires MDOT_GAS and VEL_INFINITY_GAS to be defined. 
        
        This parameter is recalculated after every iteration and updates
        V_EXP_DUST in the equation.
        
        MDOT_DUST can be given explicitly in the inputfile in which case it 
        remains unchanged.
        
        MDOT_DUST is used to calculate the real DUST_TO_GAS ratio parameter. So
        through explicit definition of 2 parameters out of MDOT_GAS, MDOT_DUST  
        and DUST_TO_GAS_CHANGE_ML_SP you can control what the internal 
        dust-to-gas ratio should be.
        
        If DUST_TO_GAS_CHANGE_ML_SP is not given, MDOT_DUST and MDOT_GAS have 
        to be defined explicitly.
        
        '''
        
        if not self.has_key('MDOT_DUST'):
            if not self.has_key('DUST_TO_GAS_CHANGE_ML_SP'):
                raise IOError('Both MDOT_DUST and DUST_TO_GAS_CHANGE_ML_SP '+\
                              'are undefined.')
            if not self.has_key('MDOT_GAS'):
                raise IOError('Both MDOT_DUST and MDOT_GAS are undefined.')
            self['MDOT_DUST'] = float(self['DUST_TO_GAS_CHANGE_ML_SP'])\
                                /float(self['VEL_INFINITY_GAS'])\
                                *float(self['V_EXP_DUST'])\
                                *float(self['MDOT_GAS'])
        else:
            pass


    def calcMDOT_GAS(self): 
        
        '''
        Calculate the value of MDOT_GAS from the DUST_TO_GAS_RATIO_ML_SP.
        
        Requires MDOT_DUST and VEL_INFINITY_GAS to be defined. 
        
        This parameter is recalculated after every iteration and updates
        V_EXP_DUST in the equation.
        
        MDOT_GAS can be given explicitly in the inputfile in which case it 
        remains unchanged.
        
        MDOT_GAS is used to calculate the real DUST_TO_GAS ratio parameter. So
        through explicit definition of 2 parameters out of MDOT_GAS, MDOT_DUST  
        and DUST_TO_GAS_CHANGE_ML_SP you can control what the internal 
        dust-to-gas ratio should be.
        
        If DUST_TO_GAS_CHANGE_ML_SP is not given, MDOT_GAS has to be defined 
        explicitly.
        
        '''
        
        if not self.has_key('MDOT_GAS'):
            if not self.has_key('DUST_TO_GAS_CHANGE_ML_SP'):
                raise IOError('Both MDOT_GAS and DUST_TO_GAS_CHANGE_ML_SP '+\
                              'are undefined.')
            if not self.has_key('MDOT_DUST'):
                raise IOError('Both MDOT_DUST and MDOT_GAS are undefined.')
            self['MDOT_GAS'] = float(self['VEL_INFINITY_GAS'])\
                               /float(self['V_EXP_DUST'])\
                               *float(self['MDOT_DUST'])\
                               /float(self['DUST_TO_GAS_CHANGE_ML_SP'])
        else:
            pass

        
    def calcMDOT_MODE(self):
        
        '''
        Set the default value of MDOT_MODE to constant.
        
        '''
        
        if not self.has_key('MDOT_MODE'):
            self['MDOT_MODE'] = 'CONSTANT'
        else:
            pass



    def calcMDOT_GAS_START(self):
        
        '''
        Set the default value of MDOT_GAS_START equal to MDOT_GAS.
        
        '''
        
        if not self.has_key('MDOT_GAS_START'):
            self['MDOT_GAS_START'] = self['MDOT_GAS']
        else:
            pass


    def calcSHELLMASS_CLASS(self):
        
        '''
        Set the order of magnitude of SHELLMASS = Mdot/v_inf. 
        0: Mdot/v_inf < 5e-8
        1: 5e-8 <= Mdot/v_inf < 2e-7
        2: 2e-7 <= Mdot/v_inf < 5e-7
        3: 5e-7 <= Mdot/v_inf
        
        '''
        
        if not self.has_key('SHELLMASS_CLASS'):
            if self['SHELLMASS'] < 5e-8: 
                #self['SHELLMASS_CLASS'] = (0,r'$\dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} < 5 \times 10^{-8}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}\ \mathrm{km}^{-1}\ \mathrm{s}$')
                self['SHELLMASS_CLASS'] = (0,r'$\dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} < 5 \times 10^{-8}$')
            elif self['SHELLMASS'] >= 3e-7: 
                #self['SHELLMASS_CLASS'] = (3,r'$\dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} \geq 3 \times 10^{-7}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}\ \mathrm{km}^{-1}\ \mathrm{s}$')
                self['SHELLMASS_CLASS'] = (3,r'$\dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} \geq 3 \times 10^{-7}$')
            elif self['SHELLMASS'] >= 5e-8 and self['SHELLMASS'] < 1.0e-7: 
                self['SHELLMASS_CLASS'] = (1,r'$5 \times 10^{-8}$ $\leq \dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} < 1 \times 10^{-7}$') 
                #self['SHELLMASS_CLASS'] = (1,r'$5 \times 10^{-8}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}\ \mathrm{km}^{-1}\ \mathrm{s}$ $\leq \dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} < 1 \times 10^{-7}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}\ \mathrm{km}^{-1}\ \mathrm{s}$') 
            else: 
                #self['SHELLMASS_CLASS'] = (2,r'$1 \times 10^{-7}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}\ \mathrm{km}^{-1}\ \mathrm{s}$ $\leq \dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} < 3 \times 10^{-7}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}\ \mathrm{km}^{-1}\ \mathrm{s}$')
                self['SHELLMASS_CLASS'] = (2,r'$1 \times 10^{-7}$ $\leq \dot{M}_\mathrm{g}/v_{\infty\mathrm{,g}} < 3 \times 10^{-7}$') 
        else:
            pass
        
        

    def calcMDOT_CLASS(self):
        
        '''
        Set the order of magnitude of MDOT. 
        0: Mdot < 3e-7
        1: 3e-7 <= Mdot < 3e-6
        2: 3e-6 <= Mdot < 1e-5
        3: 1e-5 <= Mdot
        
        '''
        
        if not self.has_key('MDOT_CLASS'):
            if self['MDOT_GAS'] < 3e-7: 
                self['MDOT_CLASS'] = (0,r'$\dot{M}_\mathrm{g} < 3 \times 10^{-7}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}$')
            elif self['MDOT_GAS'] >= 1e-5: 
                self['MDOT_CLASS'] = (3,r'$\dot{M}_\mathrm{g} \geq 1 \times 10^{-5}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}$')
            elif self['MDOT_GAS'] >= 3e-7 and self['MDOT_GAS'] < 3e-6: 
                self['MDOT_CLASS'] = (1,r'$3 \times 10^{-7}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}$ $\leq \dot{M}_\mathrm{g} < 3 \times 10^{-6}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}$') 
            else: 
                self['MDOT_CLASS'] = (2,r'$3 \times 10^{-6}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}$ $\leq \dot{M}_\mathrm{g} < 1 \times 10^{-5}\ \mathrm{M}_\odot\ \mathrm{yr}^{-1}$') 
        else:
            pass
        
        
    def calcQ_STAR(self):
        
        ''' 
        Set the stellar pulsation constant (che1992).
        
        '''
        
        if not self.has_key('Q_STAR'):
            self['Q_STAR'] = self['P_STAR']*self['M_STAR']**0.5\
                                    *self['R_STAR']**(-3/2.)
        else:
            pass
        
    
    def calcSCD_CLASS(self):
        
        '''
        Set the order of magnitude of SHELLCOLDENS. 
        0: scd < 0.06
        1: 0.06 <= scd < 0.15
        2: 0.15 <= scd < 0.4
        3: 0.4 <= scd
        
        '''
        
        if not self.has_key('SCD_CLASS'):
            if self['SHELLCOLDENS'] < 0.07: 
                self['SCD_CLASS'] = (0,r'$\bar{m} < 0.07\ \mathrm{g\;cm}^{-2}$')
            elif self['SHELLCOLDENS'] >= 0.07 and self['SHELLCOLDENS'] < 0.15: 
                self['SCD_CLASS'] = (1,r'$0.07\ \mathrm{g\;cm}^{-2}$ $\leq \bar{m} < 0.15\ \mathrm{g\;cm}^{-2}$')
            elif self['SHELLCOLDENS'] >=0.4: 
                self['SCD_CLASS'] = (3,r'$\bar{m} \geq 0.4\ \mathrm{g\;cm}^{-2}$')
            else: 
                self['SCD_CLASS'] = (2,r'$0.15\ \mathrm{g\;cm}^{-2}$ $\leq \bar{m} < 0.4\ \mathrm{g\;cm}^{-2}$')
        else:
            pass
            
    
    def calcL_CLASS(self):
        
        '''
        Set the order of magnitude of L_STAR.
        
        0: lstar < 6000
        1: 6000 <= lstar < 8000
        2: 8000 <= lstar < 10000
        3: 10000 <= lstar
        
        '''
        
        if not self.has_key('L_CLASS'):
            if self['L_STAR'] < 6000: 
                self['L_CLASS'] = (0,r'$L_\star < 6000$ $\mathrm{L}_\odot$')
            elif self['L_STAR'] >= 10000: 
                self['L_CLASS'] = (3,r'$L_\star \geq 10000$ $\mathrm{L}_\odot$')
            elif self['L_STAR'] >= 8000 and self['L_STAR'] < 10000: 
                self['L_CLASS'] = (2,r'$8000$ $\mathrm{L}_\odot$ $\leq L_\star < 10000$ $\mathrm{L}_\odot$')
            else: 
                self['L_CLASS'] = (1,r'$6000$ $\mathrm{L}_\odot$ $\leq L_\star < 8000$ $\mathrm{L}_\odot$')
        else:
            pass
        
    
    
    def calcSCATTYPE(self):
        
        '''
        Set the default scattering type for MCMax to 'ISOTROPIC'. 
        
        Can also be 'NONE', or 'FULL'. 
        
        'FULL' requires dust opacity .particle files with full scattering 
        matrices.
        
        '''
        
        if not self.has_key('SCATTYPE'):
            self['SCATTYPE'] = 'ISOTROPIC'
        else:
            pass
    
    
    
    def calcT_CONTACT(self):
        
        '''
        Set thermal contact off. 
        
        '''
        
        if not self.has_key('T_CONTACT'):
            self['T_CONTACT'] = 0
        else:
            pass
    
        
        
    def calcT_CLASS(self):
        
        '''
        Set the order of magnitude of T_STAR.
        
        0: tstar < 2000
        1: 2000 <= tstar < 2200
        2: 2200 <= tstar < 2500
        3: 2500 <= tstar
        
        '''
        
        if not self.has_key('T_CLASS'):
            if self['T_STAR'] < 2000: 
                self['T_CLASS'] = (0,r'$T_\star < 2000\ \mathrm{K}$')
            elif self['T_STAR'] >= 2500: 
                self['T_CLASS'] = (3,r'$T_\star \geq 2500\ \mathrm{K}$')
            elif self['T_STAR'] >= 2250 and self['T_STAR'] < 2500: 
                self['T_CLASS'] = (2,r'$2250\ \mathrm{K}$ $\leq T_\star < 2500\ \mathrm{K}$') 
            else: 
                self['T_CLASS'] = (1,r'$2000\ \mathrm{K}$ $\leq T_\star < 2250\ \mathrm{K}$')
        else:
            pass
    
    
    def calcVG_CLASS(self):
        
        '''
        Set the order of magnitude of VEL_INFINITY_GAS
        
        0: vg < 10
        1: 10 <= vg < 15
        2: 15 <= vg < 20
        3: 20 <= vg
        
        '''
        
        if not self.has_key('VG_CLASS'):
            if self['VEL_INFINITY_GAS'] < 10.: 
                self['VG_CLASS'] = (0,r'$v_{\infty\mathrm{,g}} < 10\ \mathrm{km\;s}^{-1}$')
            elif self['VEL_INFINITY_GAS'] >= 20.: 
                self['VG_CLASS'] = (3,r'$v_{\infty\mathrm{,g}} \geq 20\ \mathrm{km\;s}^{-1}$')
            elif self['VEL_INFINITY_GAS'] >= 15. and self['VEL_INFINITY_GAS'] < 20.: 
                self['VG_CLASS'] = (2,r'$15\ \mathrm{km\;s}^{-1}$ $\leq v_{\infty\mathrm{,g}} < 20\ \mathrm{km\;s}^{-1}$') 
            else: 
                self['VG_CLASS'] = (1,r'$10\ \mathrm{km\;s}^{-1}$ $\leq v_{\infty\mathrm{,g}} < 15\ \mathrm{km\;s}^{-1}$') 
        else:
            pass
        
        
    
    def calcV_EXP_DUST(self):
        
        """
        Calculate dust terminal velocity from gas terminal velocity and drift.
        
        Given in km/s.
        
        """

        if not self.has_key('V_EXP_DUST'):
            self['V_EXP_DUST']= float(self['VEL_INFINITY_GAS']) \
                                    + float(self['DRIFT'])
        else:
            pass    



    def calcREDDENING(self):
    
        '''
        A boolean flag for applying interstellar reddening or not. This is 
        model (read: distance) dependent, hence belongs in Star() objects.
        
        Having this available here makes it possible to compare using reddening 
        or not. 
        
        Default value is set to 0.
        
        '''
        
        if not self.has_key('REDDENING'):
            self['REDDENING']= 0
        else:
            pass    



    def calcREDDENING_MAP(self):
    
        '''
        The interstellar extinction map used for determining the interstellar
        extinction in K-band at a given distance, in the direction of given 
        longitude and latitude (set in Star.dat). 
        
        Default is Marshall et al. 2006 (marshall), but is replaced by Drimmel 
        et al. 2003 (drimmel) in case ll and bb are outside the range of 
        availability in Marshall. 
        
        Alternatives are Arenou et al. 1992 (arenou) and Schlegel et al. 1998 
        (schlegel).
        
        '''
        
        if not self.has_key('REDDENING_MAP'):
            self['REDDENING_MAP']= 'marshall'
        else:
            pass    



    def calcREDDENING_LAW(self):
    
        '''
        The extinction law used to redden model spectra. 
        
        Default is the combination of the laws by Fitzpatrick et al. 2004 
        (Optical) and Chiar & Tielens 2006 (IR), see IvS repo for more details
        at ivs.sed.reddening. 
        
        Alternatives include cardelli1989, donnell1994, fitzpatrick1999,
        fitzpatrick2004, chiar2006.
        
        '''
        
        if not self.has_key('REDDENING_LAW'):
            self['REDDENING_LAW']= 'fitz2004chiar2006'
        else:
            pass    



    def calcRT_SED(self):
        
        '''
        Set the default value of MCMax ray-tracing of the SED to False.
        
        '''
        
        if not self.has_key('RT_SED'):
            self['RT_SED']= 0
        else:
            pass         



    def calcIMAGE(self):
        
        '''
        Set the default value of MCMax image to False.
        
        '''
        
        if not self.has_key('IMAGE'):
            self['IMAGE']= 0
        else:
            pass    


    def calcVISIBILITIES(self):
        
        '''
        Set the default value of MCMax visibilities to False.
        
        '''
        
        if not self.has_key('VISIBILITIES'):
            self['VISIBILITIES']= 0
        else:
            pass    

    def calcREDO_OBS(self):
        
        '''
        Set the default value of redoing observation files to False.
        
        '''
        
        if not self.has_key('REDO_OBS'):
            self['REDO_OBS']= 0
        else:
            pass    
        
 
    
    def calcR_MAX(self,missing_key):
        
        """
        Calculate the maximum existence radii for dust species.
        
        Based on T_MIN_SPECIES for the species, and derived from mcmax output.
        
        If not MCMax model available, a power law is assumed. If T_MIN is not 
        given, no boundaries are assumed. 
        
        Is given in solar radii.
        
        @param missing_key: the missing max radius for a species that is needed
                            Of the format R_MAX_SPECIES.
        @type missing_key: string
        
        """
        
        if not self.has_key(missing_key):
            #- R_MAX is for T_MIN
            try: 
                tmin = float(self[missing_key.replace('R_MAX','T_MIN',1)])
                if self['LAST_MCMAX_MODEL']:
                    species = missing_key[6:]
                    rad = self.getDustRad(species=species,unit='rstar')
                    temp = self.getDustTemperature(species=species)
                    temp_interp = interp1d(rad,temp)
                    try: 
                        #-- interp1d raises a ValueError when outside of bounds
                        rmax = temp_interp(tmin)
                    except ValueError:
                        #- if T_MIN (for R_MAX) > temp[0] then no dust can
                        #- be present of the species
                        #- RMAX should be made very small, ie R*
                        if tmin > temp[0]:
                            rmax = self['R_STAR']

                        #- on the other hand, if TMIN < temp[-1]
                        #- all dust is allowed, so raise KeyError, such that
                        #- R_MAX is set to None
                        elif tmin < temp[-1]:
                            raise KeyError
                        
                        #- In the alternative cases, something is seriously 
                        #  wrong!
                        else:
                            raise IOError('Something went wrong when '+\
                                          'searching for R_MAX corresponding'+\
                                          ' to T_MIN... Debug!')
                else:
                    #-- Grab the power typical for optically thin temperature
                    #   profiles in the Rayleigh limit, i.e. s=1
                    power = -2./(4+1)
                    #-- Take the reciprocal relation of the dust temperature
                    #   profile (such as in Profiler.dustTemperaturePowerLaw())
                    #   rmax is in Rstar!
                    rmax = (tmin/self['T_STAR'])**(1/power)/2.
                self[missing_key] = rmax
            except KeyError:
                self[missing_key] = None
        else:
            pass        


                            
    def calcT_DES(self,sp):
        
        """
        Find the max temperature at which a dust species can exist.
        
        First, the CC inputfile is searched for T_MAX_SPECIES, in which case
        the sublimation temperature is constant. T_MAX is never made by Star()!
        
        If not present, Dust.dat info is taken, being either a sublimation 
        temperature, or the coefficients to calculate a pressure dependent
        sublimation temperature. These are set using T_DESA_ and T_DESB_SPECIES
        
        Note that tdesa and tdesb from Dust.dat are the coefficients given in 
        Kama et al 2009. MCMax uses a slightly different definition, and the 
        notation has crossed that of the paper. In any case, MCMax's definition
        of tdesa and tdesb is defined as shown here.
        
        This assumes TDESITER to be on.
        
        @param sp: The dust species
        @type sp: string
        
        """
        
        if not self.has_key('T_DESA_' + sp) \
                or not self.has_key('T_DESB_' + sp):
            try:
                #-- Check if a max temperature was given for the dust species
                if not self['T_MAX_' + sp]:
                    del self['T_MAX_' + sp]
                #-- If yes, then set the coefficients so that a constant 
                #   condensation temperature is used
                self['T_DESA_'+sp] = 10.0**(4)/self['T_MAX_' + sp]
                self['T_DESB_'+sp] = 10.0**(-4)
            except KeyError:
                if self.dust[sp]['tdesa']:
                    #-- Set the MCMax coefficients based on the Kama 2009 vals
                    self['T_DESA_'+sp] = 1e4*self.dust[sp]['tdesb']\
                                                /self.dust[sp]['tdesa']
                    self['T_DESB_'+sp] = 1e4/self.dust[sp]['tdesa']
                else:
                    #-- If no coefficients are available, use a more or less 
                    #   constant condensation temperature based on TDES in 
                    #   Dust.dat.
                    self['T_DESA_'+sp] = 1e4/self.dust[sp]['tdes']
                    self['T_DESB_'+sp] = 1e-4
        else:
            pass




    def calcR_SHELL_UNIT(self):
        
        '''
        Set default value of R_SHELL_UNIT to R_STAR.
        
        '''
        
        if not self.has_key('R_SHELL_UNIT'):
            self['R_SHELL_UNIT'] = 'R_STAR'
        else:
            pass



    def getAverageDrift(self):
        
        '''
        Return an array with the average drift velocity as a function of 
        radius, from coolfgr_all, in cm/s.
        
        '''
        
        inputfile = self.getCoolFn(ftype='fgr_all')
        drift = DataIO.getGastronoomOutput(inputfile,keyword='VDRIFT')  
        opa_gs_max = 2.5e-1
        opa_gs_min = 5.0e-3
        return array(drift)/sqrt(0.25)*1.25\
                            *(opa_gs_max**(-2.)-opa_gs_min**(-2.))\
                            /(opa_gs_max**(-2.5)-opa_gs_min**(-2.5))
        


    def calcDENSTYPE(self):
        
        """
        Define the type of density distribution.
        
        Default is 'MASSLOSS' for first iteration, otherwise SHELLFILE.
        
        If second iteration, a DENSFILE is created taking into account the 
        acceleration zone. This file is only created if not already present. 
        
        The dust density profile is calculated from the h2 number density, 
        after scaling to the dust mass-loss rate and correcting for the dust
        velocity profile. 
        
        """
        
        if not self.has_key('DENSTYPE') or not self.has_key('DENSFILE'):
            if self['MDOT_MODE'] != 'CONSTANT':
                exstr = '_var'
            else:
                exstr = ''
            filename = os.path.join(cc.path.gout,'data_for_mcmax',\
                                    '_'.join(['dens',\
                                              self['LAST_GASTRONOOM_MODEL'],\
                                    'mdotd%s%.2e.dat'%(exstr,\
                                                       self['MDOT_DUST'])]))
            if os.path.isfile(filename):
                self['DENSFILE'] = filename
                self['DENSTYPE'] = "SHELLFILE"
            else:
                if self['LAST_GASTRONOOM_MODEL']:
                    if self.has_key('DENSTYPE'):
                        if self['DENSTYPE'] == "MASSLOSS": 
                            raise IOError
                    #-- Grab the velocity profile so the gas velocity can be 
                    #   converted to the dust velocity.
                    rad = self.getGasRad()
                    vg = self.getGasVelocity()
                    #-- Use the H2 density profile to take into account any 
                    #   type of variable mass loss (including exponents in 
                    #   r_points_mass_loss.
                    nh2 = self.getGasNumberDensity()
                    #-- Get the drift profile, corrected for the average grain 
                    #   size
                    drift = self.getAverageDrift()     
                    self['DENSTYPE'] = "SHELLFILE"
                    #-- Calc dust density based on md/mg instead of d2g to take
                    #   into account velocity profiles instead of terminal vels
                    dens = nh2*self.mh*2.*self['MDOT_DUST']/self['MDOT_GAS']\
                            *vg/(vg+drift)
                    #-- GASTRoNOoM calculates smoother density profiles than 
                    #   this formula ever can accomplish
                    #dens = float(self['MDOT_DUST'])*self.Msun\
                    #        /((vg+drift)*4.*pi*rad**2.*self.year)
                    self['DENSFILE'] = filename
                    DataIO.writeCols(filename,[rad/self.au,dens])        
                    print '** Made MCMax density input file at %s.'%filename
                else:
                    print '** Writing and/or reading DENSFILE output and/or '+\
                          'input failed. Assuming standard mass-loss density'+\
                          ' distribution.'
                    self['DENSTYPE'] = "MASSLOSS"
                    self['DENSFILE'] = ''
        else:
            pass



    def calcSHELLMASS(self):
        
        """
        Calculate the average mass per shell of the circumstellar envelope. 

        Calculated by Mdot_gas/vexp.
    
        """
        
        if not self.has_key('SHELLMASS'):
            #self['SHELLMASS'] = float(self['MDOT_GAS'])*self.Msun\
                                  #/((self['VEL_INFINITY_GAS']*10**5)*self.year)
            self['SHELLMASS'] = self['MDOT_GAS']/self['VEL_INFINITY_GAS']
        else:
            pass


    def calcSHELLDENS(self):
        
        """
        Calculate the average density of the circumstellar envelope. 
    
        """
        
        if not self.has_key('SHELLDENS'):
            self['SHELLDENS'] = float(self['MDOT_GAS'])*self.Msun\
                                  /((self['VEL_INFINITY_GAS']*10**5)*self.year\
                                    *(self['R_STAR']*self.Rsun)**2*4.*pi)
        else:
            pass
        
        
        
    def calcSHELLCOLDENS(self):
        
        """
        Calculate a proxy for the average column density of the circumstellar 
        shell. 
        
        This is (intuitively) rho * R_STAR, which is important for radiative 
        excitation (density tracing the source of the radiation, R_STAR setting
        the scale of the envelope). Two types of radiative excitation can be 
        related to this: direct stellar light, and thermal dust emission.
        
        Especially important for water, but in a balance with other excitation
        mechanisms.
        
        Note that this quantity is also related to the optical depth through
        tau = kappa*coldens.
    
        """
        
        if not self.has_key('SHELLCOLDENS'):
            self['SHELLCOLDENS'] = self['SHELLDENS']\
                                    *self['R_STAR']*self.Rsun
        else:
            pass
        
        
    def calcSHELLDENS2(self):
        
        """
        Calculate a proxy for the average degree of collisional excitation in 
        the circumstellar shell.
        
        This is (intuitively) sqrt(rho * rho * R_STAR): two density factors 
        tracing both collisional partners, and R_STAR setting the scale of the 
        envelope.
        
        Sqrt is taken for easy comparison between this and the mass-loss rate
        to which it is directly proportional.
        
        Especially important for CO, but also to some degree for water where it 
        is in balance with other excitation mechanisms.
        
        Calculated by taking SHELLDENS**2*R_STAR ~ R_STAR^3/2.
        
        """
        
        if not self.has_key('SHELLDENS2'):
            self['SHELLDENS2'] = sqrt(self['SHELLDENS']**2\
                                         *self['R_STAR']*self.Rsun)
        else:
            pass

        
    def calcDENSFILE(self):
        
        """
        Pointer to the calcDENSTYPE method in case DENSFILE is missing.
    
        """
        
        self.calcDENSTYPE()
        
                            
                                                            
    def calcMRN_DUST(self):
        
        '''
        Set the default value for MRN_DUST to 0.
        
        '''
        
        if not self.has_key('MRN_DUST'):
            self['MRN_DUST'] = 0
        else:
            pass
    
    
    
    def calcMRN_INDEX(self):
    
        '''
        Set the default value for MRN_INDEX to 3.5 (standard power law in ISM
        dust grain size distribution).
        
        '''
        
        if not self.has_key('MRN_INDEX'):
            self['MRN_INDEX'] = 3.5
        else:
            pass
        
        
    
    def calcMRN_NGRAINS(self):
        
        '''
        Set the default balue for MRN_NGRAINS to the max number of dust species
        involved. 
        
        This means that all dust species are treated in the mrn treatment of 
        MCMax. 
        
        If the max is set to less species, then the extra species are treated 
        as normal, with manually set abundances. 
        
        '''
        
        if not self.has_key('MRN_NGRAINS'):
            self['MRN_NGRAINS'] = len(self.getDustList())
        else:
            pass
        
        
    
    def calcMRN_RMAX(self):
        
        '''
        Set the default value for the maximum grain size in micron. 
        
        Abundances of bigger grains will be set to 0.
        
        '''
        
        if not self.has_key('MRN_RMAX'):
            self['MRN_RMAX'] = 1000.
        else:
            pass
        
        
    
    def calcMRN_RMIN(self):
        
        '''
        Set the default value for the minimum grain size in micron. 
        
        Abundances of smaller grains will be set to 0.
        
        '''
        
        if not self.has_key('MRN_RMIN'):
            self['MRN_RMIN'] = 0.01
        else:
            pass
        
    
    def calcSCSET(self):
        
        '''
        Set default of self-consistent settling to False.
        
        '''
        
        if not self.has_key('SCSET'):
            self['SCSET'] = 0
        else:
            pass
        
    
    def calcSCSETEQ(self):
        
        '''
        Set default of self-consistent settling to True.
        
        Only relevant if SCSET == 1.
        
        '''
        
        if not self.has_key('SCSETEQ'):
            self['SCSETEQ'] = 1
        else:
            pass
        
        
        
    def calcALPHATURB(self):
        
        '''
        Set default of the turbulent mixing strenght to 1e-4.
                
        '''
        
        if not self.has_key('ALPHATURB'):
            self['ALPHATURB'] = 1e-4
        else:
            pass
        
        
    
    def calcMOLECULE(self):
        
        '''
        Set the MOLECULE keyword to empty list if not given in the input.
        
        '''
        
        if not self.has_key('MOLECULE'):
            self['MOLECULE'] = []
        else:
            pass
            


    def calcGAS_LIST(self):
        
        """
        Set the GAS_LIST keyword based on the MOLECULE keyword. 
        
        The input MOLECULE format from the CC input is converted into 
        Molecule() objects.
        
        """
        
        if not self.has_key('GAS_LIST') and self['MOLECULE']:
            if len(self['MOLECULE'][0]) == 2:
                #- First convert the long GASTRoNOoM input molecule names to 
                #- the short names, since if len() is 2, it comes from 
                #- PlottingSession.setPacsFromDb
                molec_indices \
                    = [DataIO.getInputData(keyword='MOLEC_TYPE',make_float=0,\
                                           filename='Molecule.dat')\
                                          .index(molec[0]) 
                       for molec in self['MOLECULE']]
                molecules_long = [molec[0] for molec in self['MOLECULE']]
                self['MOLECULE'] \
                    = [[DataIO.getInputData(keyword='TYPE_SHORT',\
                                            filename='Molecule.dat')[index]] \
                        + [molec[1]] 
                       for molec,index in zip(self['MOLECULE'],molec_indices)]
                self['TRANSITION'] \
                    = [[DataIO.getInputData(keyword='TYPE_SHORT',\
                                            filename='Molecule.dat')\
                            [molec_indices[molecules_long.index(trans[0])]]] \
                        + trans[1:]
                       for trans in self['TRANSITION']]
                #- Pull the info from the db
                self['GAS_LIST'] = []
                for molec,model_id in self['MOLECULE']:
                    self['GAS_LIST'].append(Molecule.makeMoleculeFromDb(\
                                        model_id=model_id,molecule=molec,\
                                        path_gastronoom=self.path_gastronoom))
            else:
                for key,index in zip(['R_OUTER','CHANGE_FRACTION_FILENAME',\
                                             'SET_KEYWORD_CHANGE_ABUNDANCE',\
                                             'NEW_TEMPERATURE_FILENAME',\
                                             'SET_KEYWORD_CHANGE_TEMPERATURE',\
                                             'ENHANCE_ABUNDANCE_FACTOR',\
                                             'ABUNDANCE_FILENAME'],\
                                            [13,16,17,18,19,15,14]):
                    if self['%s_H2O'%key]:
                        self['MOLECULE'] = \
                            [[(i==index and molec[0] in ['1H1H16O','p1H1H16O',\
                                                         '1H1H17O','p1H1H17O',\
                                                         '1H1H18O','p1H1H18O']) 
                                    and self['%s_H2O'%key] 
                                    or str(entry) 
                              for i,entry in enumerate(molec)]
                             for molec in self['MOLECULE']]                        
                #-- Check if startype is not BB, because if starfile is given 
                #   while BB is requested, the starfile cannot be given to the
                #   Molecule class. 
                starfile = self['STARTYPE'] != 'BB' and self['STARFILE'] or ''
                self['GAS_LIST'] = \
                    [Molecule.Molecule(\
                        molecule=molec[0],ny_low=int(molec[1]),\
                        ny_up=int(molec[2]),nline=int(molec[3]),\
                        n_impact=int(molec[4]),n_impact_extra=int(molec[5]),\
                        abun_molec=float(molec[6]),\
                        abun_molec_rinner=float(molec[7]),\
                        abun_molec_re=float(molec[8]),\
                        rmax_molec=float(molec[9]),itera=int(molec[10]),\
                        lte_request=int(molec[11]),\
                        use_collis_radiat_switch=int(molec[12]),\
                        abundance_filename=molec[14],\
                        enhance_abundance_factor=float(molec[15]),\
                        opr=self['OPR'],\
                        ratio_12c_to_13c=self['RATIO_12C_TO_13C'],\
                        ratio_16o_to_18o=self['RATIO_16O_TO_18O'],\
                        ratio_16o_to_17o=self['RATIO_16O_TO_17O'],\
                        r_outer=float(molec[13]) \
                                    and float(molec[13]) \
                                    or self['R_OUTER_GAS'],\
                        outer_r_mode=float(molec[13]) \
                                        and 'FIXED' \
                                        or self['OUTER_R_MODE'],\
                        dust_to_gas_change_ml_sp=self\
                                                 ['DUST_TO_GAS_CHANGE_ML_SP'],\
                        set_keyword_change_abundance=int(molec[17]),\
                        change_fraction_filename=molec[16],\
                        set_keyword_change_temperature=int(molec[19]),\
                        new_temperature_filename=molec[18],\
                        starfile=starfile)
                     
                     for molec in self['MOLECULE']]
            
            #- safety check
            requested_molecules = set([molec.molecule 
                                      for molec in self['GAS_LIST']])
            if not len(self['GAS_LIST']) == len(requested_molecules): 
                raise IOError('Multiple parameter sets for a single molecule'+\
                              ' passed. This is impossible! Contact Robin...')     
            print 'Gas molecules that are taken into account are ' + \
                  ', '.join(sorted([molec[0] for molec in self['MOLECULE']]))+\
                  '.'        
        elif not self.has_key('GAS_LIST') and not self['MOLECULE']:
            self['GAS_LIST'] = []
        else:
            pass



    def calcR_OUTER_H2O(self):
        
        '''
        Set default value of R_OUTER_H2O to 0.
        
        '''
        
        if not self.has_key('R_OUTER_H2O'):
            self['R_OUTER_H2O'] = 0
        else:
            pass 
      


    def calcNEW_TEMPERATURE_FILENAME_H2O(self):
        
        '''
        Set default value of NEW_TEMPERATURE_FILENAME_H2O to ''.
        
        '''
        
        if not self.has_key('NEW_TEMPERATURE_FILENAME_H2O'):
            self['NEW_TEMPERATURE_FILENAME_H2O'] = ''
        else:
            pass 



    def calcCHANGE_FRACTION_FILENAME_H2O(self):
        
        '''
        Set default value of CHANGE_FRACTION_FILENAME_H2O to ''.
        
        '''
        
        if not self.has_key('CHANGE_FRACTION_FILENAME_H2O'):
            self['CHANGE_FRACTION_FILENAME_H2O'] = ''
        else:
            pass 
      


    def calcSET_KEYWORD_CHANGE_TEMPERATURE_H2O(self):
        
        '''
        Set default value of SET_KEYWORD_CHANGE_TEMPERATURE_H2O to ''.
        
        '''
        
        if not self.has_key('SET_KEYWORD_CHANGE_TEMPERATURE_H2O'):
            self['SET_KEYWORD_CHANGE_TEMPERATURE_H2O'] = 0
        else:
            pass 



    def calcSET_KEYWORD_CHANGE_ABUNDANCE_H2O(self):
        
        '''
        Set default value of SET_KEYWORD_CHANGE_ABUNDANCE_H2O to ''.
        
        '''
        
        if not self.has_key('SET_KEYWORD_CHANGE_ABUNDANCE_H2O'):
            self['SET_KEYWORD_CHANGE_ABUNDANCE_H2O'] = 0
        else:
            pass 
      


    def calcENHANCE_ABUNDANCE_FACTOR_H2O(self):
        
        '''
        Set default value of ENHANCE_ABUNDANCE_FACTOR_H2O to ''.
        
        '''
        
        if not self.has_key('ENHANCE_ABUNDANCE_FACTOR_H2O'):
            self['ENHANCE_ABUNDANCE_FACTOR_H2O'] = 0
        else:
            pass 
      


    def calcABUNDANCE_FILENAME_H2O(self):
        
        '''
        Set default value of ABUNDANCE_FILENAME_H2O to ''.
        
        '''
        
        if not self.has_key('ABUNDANCE_FILENAME_H2O'):
            self['ABUNDANCE_FILENAME_H2O'] = ''
        else:
            pass 
      

    
    def calcR_OUTER_EFFECTIVE(self):
        
        '''
        Get the effective outer radius (either from Mamon, or a fixed value).
        
        '''
        
        filename = os.path.join(cc.path.gout,'models',\
                                self['LAST_GASTRONOOM_MODEL'],\
                                'input%s.dat'%self['LAST_GASTRONOOM_MODEL'])
        
        if not self.has_key('R_OUTER_EFFECTIVE'):
            self['R_OUTER_EFFECTIVE'] \
                = float(DataIO.readFile(filename=filename,delimiter=' ')[0][4])
        else:
            pass



    def calcKEYWORD_DUST_TEMPERATURE_TABLE(self):
        
        '''
        Set KEYWORD_DUST_TEMPERATURE_TABLE to False for now. 
        
        If it was not yet defined, there is not ftemperature file anyway.
        
        '''
        
        if not self.has_key('KEYWORD_DUST_TEMPERATURE_TABLE'):
            if self['DUST_TEMPERATURE_FILENAME']:
                self['KEYWORD_DUST_TEMPERATURE_TABLE'] = 1
            else:
                self['KEYWORD_DUST_TEMPERATURE_TABLE'] = 0
        else:
            pass
    

    
    def calcNUMBER_INPUT_DUST_TEMP_VALUES(self):
        
        '''
        Set NUMBER_INPUT_DUST_TEMP_VALUES to length of input file for dust temp 
        
        If it does not exist set to 0.
        
        '''
        
        if not self.has_key('NUMBER_INPUT_DUST_TEMP_VALUES'):
            if self['DUST_TEMPERATURE_FILENAME']:
                self['NUMBER_INPUT_DUST_TEMP_VALUES'] \
                    = len([1 
                           for line in DataIO.readFile(\
                                            self['DUST_TEMPERATURE_FILENAME']) 
                           if line])
            else:
                self['NUMBER_INPUT_DUST_TEMP_VALUES'] = 0
        else:
            pass  

    
                            
    def calcDUST_TEMPERATURE_FILENAME(self):
        
        """
        Look for the temperature stratification of the star.
    
        If a last mcmax model is available, the filename is given, (for now 2d). 
        
        Else an empty string is given, and a power law is used in GASTRoNOoM.
        
        """
        
        if not self.has_key('DUST_TEMPERATURE_FILENAME'):
            filename = self['RID_TEST'] != 'R_STAR' \
                            and os.path.join(cc.path.mout,\
                                             'data_for_gastronoom',\
                                             '_'.join(['Td',\
                                                     self['LAST_MCMAX_MODEL'],\
                                                     self['RID_TEST']\
                                                      +'.dat']))\
                            or os.path.join(cc.path.mout,\
                                            'data_for_gastronoom',\
                                            '_'.join(['Td',\
                                                      self['LAST_MCMAX_MODEL']\
                                                      + '.dat']))
            if os.path.isfile(filename):
                self['DUST_TEMPERATURE_FILENAME'] = filename
            else:
                if self['LAST_MCMAX_MODEL']:
                    rad = self.getDustRad(unit='rstar')
                    temp = self.getDustTemperature()
                    if self['RID_TEST'] == 'R_STAR':
                         temp = temp[rad > 1]
                         rad = rad[rad > 1]
                    elif self['RID_TEST'] == 'R_INNER_GAS':
                         temp = temp[rad > self['R_INNER_GAS']]
                         rad = rad[rad > self['R_INNER_GAS']]
                    elif self['RID_TEST'] == 'BUGGED_CASE':
                         temp = temp[rad > self['R_STAR']]
                         rad = rad[rad > self['R_STAR']]
                    self['DUST_TEMPERATURE_FILENAME'] = filename
                    DataIO.writeCols(filename,[rad,temp])
                    self['KEYWORD_DUST_TEMPERATURE_TABLE'] = 1
                    self['NUMBER_INPUT_DUST_TEMP_VALUES'] = len(rad)
                    print '** Made dust temperature stratifaction file at %s.'\
                          %filename
                    if self['NUMBER_INPUT_DUST_TEMP_VALUES'] > 999:
                        ss = '** WARNING! The dust temperature file contains'+\
                             ' more than 999 grid points. GASTRoNOoM will '+\
                             'fail because it requires less than 1000 points.'
                        print ss
                else:
                    self['DUST_TEMPERATURE_FILENAME'] = ''
        else:
            pass                        
                            
                            
  
    def calcGAS_LINES(self):
        
        """
        Making transition line input for gas data (auto search) 
        and additional no-data lines.
        
        The Transition() objects are created then for these lines and added
        to the GAS_LINES list.
    
        """
        
        if not self.has_key('GAS_LINES'):
            self['GAS_LINES'] = list()
            #-- To make sure the GAS_LIST is done, and the conversion of 
            #   TRANSITION to the right molecule names is done 
            #   (in case of PlottingSession.setPacsFromDb is used)
            self.calcGAS_LIST()     
            
            #-- Check if specific transition were requested in addition to data            
            #   Note that these include autosearch transitions if requested
            #   (See ComboCode.py)
            if self.has_key('TRANSITION'):
                #-- Keep if molecule is available in this model
                molecules = [m.molecule for m in self['GAS_LIST']]
                self['TRANSITION'] = [trans 
                                      for trans in self['TRANSITION'] 
                                      if trans[0] in molecules]
                
                #-- Check if identical transitions both with and without n_quad 
                #   are requested: Only keep the one with. (ie in the case of 
                #   radio_autosearch returns a manually requested transition)
                trl_ras = [tr for tr in self['TRANSITION'] if len(tr) == 11]
                trl_man = [tr for tr in self['TRANSITION'] if len(tr) == 12]
                trl_filt = [tr for tr in trl_ras 
                               if tr not in [tt[:-1] for tt in trl_man]]
                self['TRANSITION'] = trl_man + trl_filt
                
                #-- transitions taken from radio db have only 11 keys, and thus 
                #   have no n_quad. This is added from the Star() object. 
                nl = [Transition.makeTransition(star=self,trans=trans) 
                      for trans in self['TRANSITION']]
                nl = [trans for trans in nl if trans]
                self['GAS_LINES'].extend(nl)
                
            #- Check if molecular line catalogues have to be browsed to create 
            #- line lists in addition to the data
            if self['LINE_SELECT']:
                if self['LINE_SELECT'] == 1: 
                    self.__addLineList()
                elif self['LINE_SELECT'] == 2: 
                    trl = DataIO.readDict(self['LS_FILE'],\
                                          multi_keys=['TRANSITION'])
                    trl_sorted = DataIO.checkEntryInfo(trl['TRANSITION'],11,\
                                                       'TRANSITION')
                    nl = [Transition.makeTransition(trans=trans,star=self) 
                          for trans in trl_sorted]
                    nl = [trans for trans in nl if trans]
                    
                    #-- Kick out those transitions that were requested manually
                    #   Can use GAS_LINES key, as it contains both manual and 
                    #   other lines, but the latter are also assigned 
                    #   self['N_QUAD'] anyway
                    ctrl = [tr.getInputString(include_nquad=0) 
                            for tr in self['GAS_LINES']]
                    nl = [tr for tr in nl
                             if tr.getInputString(include_nquad=0) not in ctrl]
                    self['GAS_LINES'].extend(nl)    
    
            #-- Sort the transitions.
            self['GAS_LINES'] = sorted(list(self['GAS_LINES']),\
                                       key=lambda x: str(x))
            #-- Check uniqueness. Same N_QUAD double transitions can still occur
            self['GAS_LINES'] = Transition.checkUniqueness(self['GAS_LINES'])
            
            #-- Is this still needed? 
            requested_transitions = set([str(trans) 
                                         for trans in self['GAS_LINES']]) 
            if not len(self['GAS_LINES']) == len(requested_transitions):
                print 'Length of the requested transition list: %i'\
                      %len(self['GAS_LINES'])
                print 'Length of the requested transition list with only ' + \
                      'the "transition string" parameters: %i'\
                      %len(requested_transitions)
                print 'Guilty transitions:'
                trans_strings = [str(trans) for trans in self['GAS_LINES']]
                print '\n'.join([str(trans) 
                                 for trans in self['GAS_LINES'] 
                                 if trans_strings.count(str(trans))>1])
                raise IOError('Multiple parameter sets for a single ' + \
                              'transition requested. This is impossible! '+ \
                              'Check if N_QUAD and TRANSITION definitions ' + \
                              'have the same value in your inputfile.' + \
                              'If yes, check code/contact Robin.')
        else:
            pass
    


    def calcSTARTYPE(self):
        
        """
        Set the default value for STARTYPE, which is the blackbody 
        assumption (BB). 
        
        """
        
        if not self.has_key('STARTYPE'):
            self['STARTYPE'] = 'BB'
        else:
            pass



    def calcSTARFILE(self):
        
        """
        Set the default value for STARFILE, which is an empty string 
        (ie STARTYPE is BB, no inputfile). 
        
        """
        
        if not self.has_key('STARFILE'):
            if self['STARTYPE'] == 'BB':
                self['STARFILE'] = ''
            elif self['STARTYPE'] == 'ATMOSPHERE':
                modeltypes = ['comarcs','marcs','kurucz']
                modeltype = None
                for mt in modeltypes:
                    if mt in self['ATM_FILENAME']:
                        modeltype = mt
                        continue
                if modeltype is None: 
                    raise IOError('Atmosphere model type is unknown.')
                DataIO.testFolderExistence(cc.path.starf)
                atmfile = self['ATM_FILENAME']
                atmos = Atmosphere.Atmosphere(modeltype,filename=atmfile)
                atmosmodel = atmos.getModel(teff=self['T_STAR'],\
                                            logg=self['LOGG'])
                starfile = os.path.join(cc.path.starf,'%s_teff%s_logg%s.dat'\
                                        %(os.path.splitext(atmos.filename)[0],\
                                        str(atmos.teff_actual),\
                                        str(atmos.logg_actual)))
                if not os.path.isfile(starfile):
                    savetxt(starfile,atmosmodel,fmt=('%.8e'))
                print 'Using input model atmosphere at '
                print starfile
                self['STARFILE'] = starfile
            elif self['STARTYPE'] == 'TABLE':
                self['STARTABLE'] = self['STARTABLE'].strip('"').strip("'")
                if not (os.isfile(self['STAR_TABLE']) \
                        and os.path.split(self['STARTABLE'])[0]):
                    self['STARFILE'] = os.path.join(cc.path.starf,\
                                                    self['STARTABLE'])
                else:
                    self['STARFILE'] = self['STARTABLE']
                print 'Using input star spectrum at '
                print self['STARFILE']
        else:
            pass



    def calcLINE_SELECT(self):
        
        ''' 
        If the LINE_SELECT keyword is not present, set to False.
        
        '''
        
        if not self.has_key('LINE_SELECT'):
            self['LINE_SELECT'] = 0
        else:
            pass



    def calcDUST_TO_GAS(self):
        
        '''
        Calculate the empirical value oft he dust to gas ratio.
        
        '''
        
        if not self.has_key('DUST_TO_GAS'):
            self['DUST_TO_GAS'] = float(self['MDOT_DUST'])\
                                    *float(self['VEL_INFINITY_GAS'])\
                                    /float(self['MDOT_GAS'])\
                                    /float(self['V_EXP_DUST'])
        else:
            pass


  
    def calcDUST_TO_GAS_INITIAL(self):
        
        '''
        Set a default value for the initial dust-to-gas ratio at 0.002.
        
        '''
        
        if not self.has_key('DUST_TO_GAS_INITIAL'):
            self['DUST_TO_GAS_INITIAL'] = 0.002
        else:
            pass  



    def calcDUST_TO_GAS_ITERATED(self):
        
        '''
        Fetch the iterated result of the dust-to-gas ratio from cooling.
        
        '''
        
        if not self.has_key('DUST_TO_GAS_ITERATED'):
            try:
                filename = os.path.join(cc.path.gout,'models',\
                                        self['LAST_GASTRONOOM_MODEL'],\
                                        'input%s.dat'\
                                        %self['LAST_GASTRONOOM_MODEL'])
                self['DUST_TO_GAS_ITERATED'] = float(DataIO.readFile(\
                                                        filename=filename,\
                                                        delimiter=' ')[0][6])
            except IOError:
                self['DUST_TO_GAS_ITERATED'] = None
        else:
            pass  



    def getOpticalDepth(self,wavelength=0):
        
        '''
        Calculate the optical depth.
        
        If wavelength keyword is given, tau at wavelength is returned. 
        
        Otherwise, the full wavelength array is returned.
        
        @keyword wavelength: the wavelength in micron. If 0, the whole 
                             wavelength array is returned.
                             
                             (default: 0)
        @type wavelength: float
        
        @return: The optical depth at requested wavelength or the full
                 wavelength and optical depth arrays
        @rtype: float or (array,array)
        
        '''
        
        wavelength = float(wavelength)
        rad = self.getDustRad()
        dens = self.getDustDensity()
        wave_list,kappas = self.readWeightedKappas()
        if wavelength:
            wave_index = argmin(abs(wave_list-wavelength))
            return integrate.trapz(y=dens*kappas[wave_index],x=rad)
        else:
            return (wave_list,array([integrate.trapz(y=dens*kappas[i],x=rad)
                                     for i in xrange(len(wave_list))]))
        
    
    def calcINCLUDE_SCAT_GAS(self):
        
        '''
        Set the keyword INCLUDE_SCAT_GAS to 0.
        
        The keyword decides whether to take into account the scattering 
        coefficients in GASTRoNOoM as if they contributed to the absorption
        coefficients. 
        
        '''
        
        if not self.has_key('INCLUDE_SCAT_GAS'):  
            self['INCLUDE_SCAT_GAS'] = 0
        else:
            pass
        
    
    def readWeightedKappas(self):
        
        '''
        Return the wavelength and kappas weighted with their respective dust 
        mass fractions.
        
        Typically you only want the absorption coefficients because GASTRoNOoM
        does not take into account scattering. You could try approximating 
        the effect of scattering on the acceleration, but at this point this is 
        not taken into account accurately.
        
        @return: The wavelength and weighted kappas grid
        @rtype: (array,array)
        
        '''
        
        wave_list,kappas = self.readKappas()
        if self['INCLUDE_SCAT_GAS']:
            #-- First the absorption coefficients of all dust species are given 
            #   Then the scattering coefficients. So iterate twice over the 
            #   dust list.
            wkappas = [sum([float(self['A_%s'%(species)])*float(kappas[i][j])
                            for i,species in enumerate(self.getDustList()*2)])
                       for j in xrange(len(kappas[0]))]
        else: 
            #-- Only iterate once over the dust list to just take the 
            #   absorption coefficients.
            wkappas = [sum([float(self['A_%s'%(species)])*float(kappas[i][j])
                            for i,species in enumerate(self.getDustList())])
                       for j in xrange(len(kappas[0]))]
        return array(wave_list),array(wkappas)
        
        
    
    def calcRATIO_12C_TO_13C(self):
        
        '''
        Set default value for ratio_12c_to_13c to 0.
        
        '''
        
        if not self.has_key('RATIO_12C_TO_13C'):  
            self['RATIO_12C_TO_13C'] = 0
        else:
            pass
    


    def calcRATIO_16O_TO_17O(self):
        
        '''
        Set default value for ratio_16o_to_17o to 0.
        
        '''
        
        if not self.has_key('RATIO_16O_TO_17O'):  
            self['RATIO_16O_TO_17O'] = 0
        else:
            pass
            


    def calcRATIO_16O_TO_18O(self):
        
        '''
        Set default value for ratio_16o_to_18o to 0.
        
        '''
        
        if not self.has_key('RATIO_16O_TO_18O'):  
            self['RATIO_16O_TO_18O'] = 0
        else:
            pass    
        


    def calcOPR(self):
        
        '''
        Set default value for opr to 0.
        
        '''
        
        if not self.has_key('OPR'):  
            self['OPR'] = 0
        else:
            pass                
        


    def calcUSE_NEW_DUST_KAPPA_FILES(self):
        
        '''
        Set the default value of USE_NEW_DUST_KAPPA_FILES to 1.
        
        '''
        
        if not self.has_key('USE_NEW_DUST_KAPPA_FILES'):  
            self['USE_NEW_DUST_KAPPA_FILES'] = 1
        else:
            pass         



    def calcTEMDUST_FILENAME(self):
        
        """
        Making extinction efficiency input files for GASTRoNOoM from MCMax 
        output mass extinction coefficients.
        
        If no MCMax output available, this file is temdust.kappa, the standard.
        
        In units of cm^-1, Q_ext/a.
        
        """
        
        if not self.has_key('TEMDUST_FILENAME'):
            if self['NLAM'] > 2000:
                #- For now not supported, GASTRoNOoM cannot take more than 2000
                #- wavelength opacity points
                raise IOError('NLAM > 2000 not supported due to GASTRoNOoM '+\
                              'opacities!')    
            filename = os.path.join(cc.path.gdata,'temdust_%s.dat'\
                                    %self['LAST_MCMAX_MODEL'])
            if not int(self['USE_NEW_DUST_KAPPA_FILES']) \
                    or not self['LAST_MCMAX_MODEL']:
                self['TEMDUST_FILENAME'] = 'temdust.kappa'
            elif os.path.isfile(filename):
                self['TEMDUST_FILENAME'] = os.path.split(filename)[1]
            else:
                try:
                    wavelength,q_ext = self.readWeightedKappas()
                    q_ext *= self['SPEC_DENS_DUST']*(4.0/3)
                    wavelength = list(wavelength)
                    wavelength.reverse()
                    q_ext = list(q_ext)
                    q_ext.reverse()
                    self['TEMDUST_FILENAME'] = os.path.split(filename)[1]
                    DataIO.writeCols(filename,[wavelength,q_ext])
                    print '** Made opacity file at ' + filename + '.'
                except IOError:
                    self['TEMDUST_FILENAME'] = 'temdust.kappa'
        else:
            pass    
     
    
    
    def calcR_OH1612_AS(self):
         
        '''
        Set the R_OH1612_AS to the default value of 0 as.
        
        '''
        
        if not self.has_key('R_OH1612_AS'):  
            self['R_OH1612_AS'] = 0
        else:
            pass


    
    def calcR_OH1612(self):
         
        '''
        Calculate the R_OH1612 in R_STAR.
        
        '''
        
        if not self.has_key('R_OH1612'):  
            self['R_OH1612'] = Data.convertAngular(self['R_OH1612_AS'],\
                                                   self['DISTANCE'])\
                                    /self['R_STAR']/self.Rsun
        else:
            pass         


    
    def calcR_OH1612_NETZER(self):
         
        '''
        Calculate the radial OH maser peak distance in cm.
        
        Taken from Netzer & Knapp 1987, eq. 29. 
        
        The interstellar radiation factor is taken as A = 5.4 
        (avg Habing Field)

        '''
        
        if not self.has_key('R_OH1612_NETZER'):  
            mg = self['MDOT_GAS']/1e-5
            vg = self['VEL_INFINITY_GAS']
            self['R_OH1612_NETZER'] = ((5.4*mg**0.7/vg**0.4)**-4.8\
                                        + (74.*mg/vg)**-4.8)**(-1/4.8)*1e16\
                                        /self['R_STAR']/self.Rsun
        else:
            pass         
    
    
    
    def getBlackBody(self):
        
        '''
        Calculate the black body intensity profile.
        
        @return: The wavelength and black body intensity grid
        @rtype: (array,array)
        
        '''
        
        #- Define wavelength grid in cm
        w = 10**(linspace(-9,2,5000))
        freq = self.c/w
        #- Calculate the blackbody
        bb = 2*self.h*freq**3/self.c**2 * \
             (1/(exp((self.h*freq)/(self.k*self['T_STAR']))-1))
        return w*10**(4),bb*10**(23)
    
    
    
    def getObservedBlackBody(self):
        
        '''
        Scale the blackbody intensity following the distance and stellar radius.
        
        This is not the flux!
        
        @return: The wavelength grid and rescaled blackbody intensity
        @rtype: (array,array)
        
        '''
        
        w,bb = self.getBlackBody()
        return w,bb*(self['R_STAR']*self.Rsun)**2\
                   /(self['DISTANCE']*self.pc)**2
        


    def missingInput(self,missing_key):
        
        """
        Try to resolve a missing key.
        
        @param missing_key: the missing key for which an attempt will be made 
                            to calculate its value based on already present 
                            parameters
        @type missing_key: string
        
        """
        
        if missing_key in ('T_STAR','L_STAR','R_STAR'):
            self.calcTLR()
        elif missing_key in ['R_MAX_' + species 
                             for species in self.getDustList()]:
            self.calcR_MAX(missing_key)
        elif missing_key in ['R_DES_' + species 
                             for species in self.getDustList()]:
            self.checkT()
        elif missing_key in ['T_DESA_' + species 
                             for species in self.getDustList()] + \
                            ['T_DESB_' + species 
                             for species in self.getDustList()]:
            self.calcT_DES(missing_key[7:])
        elif missing_key in ['T_DES_' + species 
                             for species in self.getDustList()]:
            self.checkT()
        elif hasattr(self,'calc' + missing_key):
            getattr(self,'calc' + missing_key)()
        else:
            pass

 