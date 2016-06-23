# -*- coding: utf-8 -*-

"""
Toolbox for Molecules, used in various applications for GASTRoNOoM.

Author: R. Lombaert

"""

import os 
import re

import cc.path
from cc.tools.io import DataIO
from cc.tools.io.Database import Database
from cc.tools.io import Radiat



def makeMoleculeFromDb(molec_id,molecule,path_gastronoom='codeSep2010',\
                       mline_db=None):
    
    '''
    Make a Molecule() from a database, based on model_id and molec_id.
    
    Returns None if the molecule is not available in the database.
    
    @param molec_id: the model_id of the molecule
    @type molec_id: string
    @param molecule: the short hand name of the molecule        
    @type molecule: string
    
    @keyword path_gastronoom: the output path in the ~/GASTRoNOoM/. directory
    
                              (default: codeSep2010)
    @type path_gastronoom: string
    @keyword mline_db: The mline database, which can be passed in case one 
                       wants to reduce overhead. Not required though.
                       
                       (default: None)
    @type mline_db: Database()
    
    @return: the molecule with all its information embedded
    @rtype: Molecule()
    
    '''
    
    #-- Convenience path
    cc.path.gout = os.path.join(cc.path.gastronoom,path_gastronoom)
    
    #-- Retrieve cooling id log
    cooling_log = os.path.join(cc.path.gout,'models',molec_id,'cooling_id.log')
    if os.path.isfile(cooling_log):
        model_id = DataIO.readFile(cooling_log)[0]
    #- ie mline id is the same as model id, the first calced for this id
    else: 
        model_id = molec_id
        
    if mline_db is None:
        molec_db = Database(os.path.join(cc.path.gout,\
                                         'GASTRoNOoM_mline_models.db'))
    else:
        molec_db = mline_db
    
    if not molec_db.has_key(model_id) \
            or not molec_db[model_id].has_key(molec_id) \
            or not molec_db[model_id][molec_id].has_key(molecule):
        return None
    molec_dict = molec_db[model_id][molec_id][molecule].copy()
    extra_pars = molec_dict['MOLECULE'].split()[1:]
    for k,v in zip(['ny_low','ny_up','nline','n_impact','n_impact_extra'],\
                   extra_pars):
        molec_dict[k] = int(v)
    for key in ['MOLECULE','CHANGE_DUST_TO_GAS_FOR_ML_SP',\
                'NUMBER_INPUT_ABUNDANCE_VALUES','KEYWORD_TABLE',\
                'MOLECULE_TABLE','ISOTOPE_TABLE']:
        if molec_dict.has_key(key):
            del molec_dict[key]
    molec_dict = dict([(k.lower(),v) for k,v in molec_dict.items()])
    molec = Molecule(molecule=molecule,**molec_dict)
    molec.setModelId(molec_id)
    return molec
    


class Molecule():
    
    '''
    A class to deal with molecules in GASTRoNOoM.
    
    '''
    
    def __init__(self,molecule,ny_up=0,ny_low=0,nline=0,n_impact=0,\
                 n_impact_extra=0,abun_molec=1.0e-10,abun_molec_rinner=1.0e-10,\
                 abun_molec_re=1.0e-10,rmax_molec=1.,itera=0,lte_request=None,\
                 use_collis_radiat_switch=0,dust_to_gas_change_ml_sp=0,\
                 use_no_maser_option=0,use_maser_in_sphinx=1,\
                 fehler=1e-4,xdex=2.,n_freq=30,start_approx=0,\
                 use_fraction_level_corr=1,fraction_level_corr=0.8,\
                 number_level_max_corr=1e-12,\
                 ratio_12c_to_13c=0,ratio_16o_to_17o=0,ratio_16o_to_18o=0,\
                 opr=0,r_outer=0,outer_r_mode='MAMON',abundance_filename=None,\
                 change_fraction_filename=None,set_keyword_change_abundance=0,\
                 set_keyword_change_temperature=0,enhance_abundance_factor=0,\
                 new_temperature_filename=None,linelist=0,starfile=''):
        
        '''
        Initiate a Molecule class, setting all values for the allowed 
        transition parameters to zero (by default).
        
        @param molecule: shorthand name of the molecule
        @type molecule: string
        
        @keyword ny_up: number of levels in first vibration state v=1
        
                        (default: 0)
        @type ny_up: int
        @keyword ny_low: number of levels in the ground vibration state v=0
        
                         (default: 0)
        @type ny_low: int
        @keyword nline: number of allowed transitions in molecule
        
                        (default: 0)
        @type nline: int  
        @keyword n_impact: number of depth points in radius mesh
        
                           (default: 0)
        @type n_impact: int       
        @keyword n_impact_extra: number of depth points in radius_mesh 
                                 (< n_impact) used for variable mass-loss 
                                 (0 if constant mdot)
        
                                 (default: 0)
        @type n_impact_extra: int            
        @keyword itera: number of iterations in mline for molecule, LTE 
                        approximation if zero
          
                        (default: 0)
        @type itera: string    
        @keyword abun_molec: molecular abundance at the stellar radius.
                             Default is arbitrary, and used if molecule is co 
                             or h2o.
          
                             (default: 1.0e-10)
        @type abun_molec: float    
        @keyword abun_molec_rinner: molecular abundance at inner shell radius.
                                    Default is arbitrary, and used if molecule
                                    is co or h2o.
          
                                    (default: 1.0e-10)
        @type abun_molec_rinner: float            
        @keyword abun_molec_re: molecular abundance at rmax_molec. Default is 
                                arbitrary, and used if molecule is co or h2o.
          
                                (default: 1.0e-10)
        @type abun_molec_re: float 
        @keyword rmax_molec: The radius from which the Willacy abundance
                             profiles are used. They are rescaled to 
                             abun_molec_re. Default is arbitrary, and used if 
                             molecule is co or h2o.
                             
                             (default: 1.)
        @type rmax_molec: float
        @keyword use_collis_radiat_switch: in case of unstable mline, such as 
                                           for para h2o sometimes
                                           
                                           (default: 0)
        @type use_collis_radiat_switch: bool
        @keyword fehler: convergence criterium in mline
        
                         (default: 1e-4)
        @type fehler: float
        @keyword xdex: Controls the distribution of the impact parameters in the
                       interval between R_STAR and R_OUTER. 
                       
                       (default: 2.)
        @type xdex: float 
        @keyword n_freq: Number of frequency points in line profile
        
                         (default: 30)
        @type n_freq: int
        @keyword start_approx: set to 0 when one wants to start with LTE-approx 
                               as starting n(NY,N_IMPACT); set to 1 when 
                               starting from another model - with same NY, 
                               N_IMPACT, ...
        
                               (default: 0)
        @type start_approx: bool
        @keyword use_fraction_level_corr: set to 1 if one wants to put a limit 
                                          on the level-population correction 
                                          (BES3).
        
                                          (default: 1)
        @type use_fraction_level_corr: bool
        @keyword fraction_level_corr: user-defined fraction for maximum change 
                                      in level-population correction; useful in 
                                      case of H2O
        
                                      (default: 0.8)
        @type fraction_level_corr: float
        @keyword number_level_max_corr: user-defined level population. Only the 
                                        level corrections for levels with a 
                                        higher level population will be used to 
                                        determine convergence criterion
        
                                        (default: 1e-12)
        @type number_level_max_corr: float
        @keyword ratio_12c_to_13c: 12c/13c ratio, only relevant for 13co and 
                                   other molecules with that isotope
                                   
                                   (default: 0)
        @type ratio_12c_to_13c: int
        @keyword ratio_16o_to_17o: 16o/17o ratio, only relevant for h2_17o and 
                                   other molecules with that isotope
                                   
                                   (default: 0)
        @type ratio_16o_to_17o: int
        @keyword ratio_16o_to_18o: 16o/18o ratio, only relevant for h2_18o and 
                                   other molecules with that isotope
                                   
                                   (default: 0)
        @type ratio_16o_to_18o: int
        @keyword opr: ortho-to-para water ratio, only relevant for ph2o,
                      ph2_17o,ph2_18o and other molecules with para h2o
                      
                      (default: 0)
        @type opr: int    
        @keyword r_outer: the outer radius of the shell for this molecule, 
                          0 if MAMON
                          
                          (default: 0)
        @type r_outer: float
        @keyword lte_request: using LTE in mline only (with or without 
                              collision rates: determined by itera), if default
                              lte_request is 0 if itera != 0 and 1 if itera ==0
        
                              (default: 0)
        @type lte_request: bool
        @keyword outer_r_mode: the mode used for calculating r_outer 
                               (FIXED or MAMON)
                               
                               (default: 'MAMON')
        @type outer_r_mode: string
        @keyword dust_to_gas_change_ml_sp: if 0 not used, otherwise this is an 
                                           alternative value for the
                                           dust-to-gas ratio in mline/sphinx 
                                           for this molecule and its 
                                           transitions.
                                           
                                           (default: 0)
        @type dust_to_gas_change_ml_sp: float
        @keyword abundance_filename: if enhance_abundance_factor is not zero, 
                                        this includes the filename and/or path 
                                        to the file that includes the profile. 
                                        
                                        (default: None)
        @type abundance_filename: string
        @keyword enhance_abundance_factor: if 0 the Willacy abundance profiles 
                                           are uses, if not zero the 
                                           abundance_filename is used and 
                                           scaled with the factor given here. 
                                           THIS DOES NOT RESCALE ABUNDANCES BY 
                                           WILLACY! Only used for filename 
                                           abundances, hence why this parameter
                                           also turns this feature on/off
                                           
                                           (default: 0)
        @type enhance_abundance_factor: float
        @keyword set_keyword_change_abundance: Change the abundance calculated 
                                               in cooling by a radius dependent
                                               factor
                                               
                                               (default: 0)
        @type set_keyword_change_abundance: bool
        @keyword change_fraction_filename: the filename of the enhancement 
                                           factors if 
                                           set_keyword_change_abundance != 0
                                           
                                           (default: None)
        @type change_fraction_filename: string
        @keyword set_keyword_change_temperature: Use a different temperature 
                                                 structure in mline and sphinx
            
                                                 (default: 0)
        @type set_keyword_change_temperature: bool
        @keyword new_temperature_filename: the filename for the temperature 
                                           structure if 
                                           set_keyword_change_temperature != 0
        
                                           (default: None)
        @type new_temperature_filename: string
        @keyword use_no_maser_option: Do not allow masers (neg opacs) in mline
                                      RT by setting negative line opacs to 1e-60
                                      If use_maser_in_sphinx is on, mline 
                                      will do a final run including masers 
                                      anyway to see what would happen if they 
                                      were inluded by allowing negative opacs
                                      for the line profile calculations in 
                                      sphinx (but not for the convergence in 
                                      mline).
                                      
                                      (default: 0)
        @type use_no_maser_option: bool
        @keyword use_maser_in_sphinx: When on, does a final mline run including 
                                      masers, allowing negative opacities. When 
                                      off, sets the masing line opacities to 
                                      1e-60 when writing out the ml3 file.
                                     
                                      (default: 1)
        @type use_maser_in_sphinx: bool
        @keyword linelist: The molecule is created for the LineList module. No
                           radiative information is read from GASTRoNOoM input
                           files.
        
                           (default: 0)
        @type linelist: bool
        @keyword starfile: input filename for a stellar input spectrum (either
                           user defined or from a model atmosphere spectrum)
                           
                           (default: '')
        @type starfile: str
                           
        
        '''
 
        self.molecule = str(molecule)
        self.ny_up = int(ny_up)
        self.ny_low = int(ny_low)
        self.nline = int(nline)
        self.n_impact = int(n_impact)
        self.n_impact_extra = int(n_impact_extra)
        self.molecule_index = DataIO.getInputData(keyword='TYPE_SHORT',\
                                                  filename='Molecule.dat')\
                                                 .index(self.molecule)
        mdata = ['MOLEC_TYPE','NAME_SHORT','NAME_PLOT',\
                 'SPEC_INDICES','USE_INDICES_DAT']
        attrs = ['molecule_full','molecule_short','molecule_plot',\
                 'spec_indices','use_indices_dat']
        mfloat = [0,0,0,1,1]
        for k,a,mf in zip(mdata,attrs,mfloat):
            setattr(self,a,DataIO.getInputData(keyword=k,make_float=mf,\
                                               filename='Molecule.dat',\
                                               rindex=self.molecule_index,))
        
        self.itera = int(itera)
        #- lte_request may be undefined, but then it would be faulty input, 
        #- where we do want the code to crash... 
        if self.itera==0 and lte_request is None:
            self.lte_request = 1
        #- Normally u never use lte if taking into account collision rates
        elif self.itera != 0:      
            self.lte_request = 0
        elif self.itera==0 and lte_request is not None:
            self.lte_request = lte_request
        self.abun_molec = abun_molec
        self.abun_molec_rinner = abun_molec_rinner
        self.abun_molec_re = abun_molec_re
        self.rmax_molec = rmax_molec
        self.use_collis_radiat_switch = int(use_collis_radiat_switch)
        self.ratio_12c_to_13c = ratio_12c_to_13c
        self.ratio_16o_to_17o = ratio_16o_to_17o
        self.ratio_16o_to_18o = ratio_16o_to_18o
        self.opr = opr
        self.dust_to_gas_change_ml_sp = float(dust_to_gas_change_ml_sp)
        self.use_no_maser_option = int(use_no_maser_option)
        self.use_maser_in_sphinx = int(use_maser_in_sphinx)
        self.fehler = fehler
        self.xdex = xdex
        self.n_freq = int(n_freq)
        self.start_approx = int(start_approx)
        self.use_fraction_level_corr = int(use_fraction_level_corr)
        self.fraction_level_corr = fraction_level_corr
        self.number_level_max_corr = number_level_max_corr
        
        #-- Set the molecule inputfiles for abundance and temperature, applying
        #   the path from Path.dat if the file does not exist or the path to the
        #   file is not given. (eg a subfolder might be given, but that works)
        for k in ['abundance_filename','change_fraction_filename',\
                  'new_temperature_filename']:
            fn = locals()[k]
            if fn and not (os.path.isfile(fn) and os.path.split(fn)[0]):
                fn = os.path.join(cc.path.molf,fn)
                setattr(self,k,fn)
            else:
                setattr(self,k,fn)
        self.enhance_abundance_factor = float(enhance_abundance_factor)
        self.set_keyword_change_abundance = int(set_keyword_change_abundance)
        self.set_keyword_change_temperature = \
                    int(set_keyword_change_temperature)

        #-- Mainly for plotting purposes: The relative, multiplicative abundance 
        #   factor with respect to main isotope (and OPR) is calculated
        #   This does not take into account enhance_abundance_factor!
        self.abun_factor = self.getAbunFactor() 
        self.outer_r_mode = outer_r_mode
        if self.outer_r_mode == 'MAMON': self.r_outer = 0
        else: self.r_outer = float(r_outer)
        self.__model_id = None
        if not linelist:
            if self.use_indices_dat:
                tag = '_'.join([self.molecule,str(self.ny_low),\
                                str(self.ny_up),str(self.nline)])
                i = DataIO.getInputData(start_index=4,keyword='MOLECULE',\
                                        filename='Indices.dat').index(tag)
                self.indices_index = i
            self.radiat = Radiat.Radiat(molecule=self)
            if self.spec_indices:
                if self.use_indices_dat:
                    f = DataIO.getInputData(path=cc.path.usr,start_index=4,\
                                            keyword='INDICES',rindex=i,\
                                            filename='Indices.dat')
                    filename = os.path.join(cc.path.gdata,'indices_backup',f)
                else:
                    filename = os.path.join(cc.path.gdata,\
                                         '{}_indices.dat'.format(self.molecule))
                rf = DataIO.readFile(filename,' ')
                self.radiat_indices = [[int(i) for i in line] for line in rf]
        else:
            self.radiat = None
            self.radiat_indices = None
        self.starfile = starfile


    def __str__(self):
        
        '''
        Printing a molecule as it should appear in the GASTRoNOoM input file.
        
        '''
        
        ll = [self.molecule_full,self.ny_low,self.ny_up,self.nline,\
              self.n_impact,self.n_impact_extra]
        return 'MOLECULE={} {:d} {:d} {:d} {:d} {:d}'.format(*ll)



    def updateParameters(self,pardict):
        
        '''
        Update parameters.
        
        @param pardict: the parameters with respective values for the update
        @type pardict: dict()
        
        '''
        
        for k,v in pardict.items():
            if hasattr(self,k.lower()): setattr(self,k.lower(),v)



    def makeDict(self,path=None,in_progress=0):
        
        '''
        Return a dict with molecule string, and other relevant parameters.
        
        @keyword path: If a different path is needed, it can be passed here, 
                       for files. For instance, when making dictionaries for 
                       Molecule() objects in the case of supercomputer copies.
                      
                       (default: None)
        @type path: string
        @keyword in_progress: add an extra dict entry "IN_PROGRESS" if the 
                              molecule is still being calculated somewhere.
                              
                              (default: 0)
        @type in_progress: bool
        
        @return: The molecule dictionary including all relevant, defining 
                 information
        @rtype: dict()
        
        '''
        
        dd = dict([('MOLECULE',str(self).replace('MOLECULE=','')),\
                   ('ITERA',self.itera),\
                   ('ABUN_MOLEC',self.abun_molec),\
                   ('ABUN_MOLEC_RINNER',self.abun_molec_rinner),\
                   ('ABUN_MOLEC_RE', self.abun_molec_re),\
                   ('RMAX_MOLEC',self.rmax_molec),\
                   ('LTE_REQUEST',self.lte_request),\
                   ('OUTER_R_MODE',self.outer_r_mode),\
                   ('USE_COLLIS_RADIAT_SWITCH',self.use_collis_radiat_switch),\
                   ('R_OUTER',self.r_outer),\
                   ('USE_NO_MASER_OPTION',self.use_no_maser_option),\
                   ('USE_MASER_IN_SPHINX',self.use_maser_in_sphinx),\
                   ('FEHLER',self.fehler),\
                   ('XDEX',self.xdex),\
                   ('N_FREQ',self.n_freq),\
                   ('START_APPROX',self.start_approx),\
                   ('USE_FRACTION_LEVEL_CORR',self.use_fraction_level_corr),\
                   ('FRACTION_LEVEL_CORR',self.fraction_level_corr),\
                   ('NUMBER_LEVEL_MAX_CORR',self.number_level_max_corr)
                   ])
        
        for par,isot in [('RATIO_16O_TO_18O','18O'),\
                         ('RATIO_16O_TO_17O','17O'),\
                         ('RATIO_12C_TO_13C','13C'),\
                         ('OPR','p1H')]:
            if isot in self.molecule:
                dd[par] = getattr(self,par.lower())
        if self.dust_to_gas_change_ml_sp:
            dd['CHANGE_DUST_TO_GAS_FOR_ML_SP'] = 1
            dd['DUST_TO_GAS_CHANGE_ML_SP'] = self.dust_to_gas_change_ml_sp

        if self.enhance_abundance_factor:
            dd['ENHANCE_ABUNDANCE_FACTOR'] = self.enhance_abundance_factor
            niav = len(DataIO.readCols(filename=self.abundance_filename)[0])
            dd['NUMBER_INPUT_ABUNDANCE_VALUES'] = niav
            dd['KEYWORD_TABLE'] = 1
            moltab = self.molecule_full[:self.molecule_full.index('.')]
            dd['MOLECULE_TABLE'] = moltab 
            dd['ISOTOPE_TABLE'] = self.molecule_full
            if path <> None:
                afn = os.path.join(path,\
                                   os.path.split(self.abundance_filename)[1])
            else:
                afn = self.abundance_filename
            dd['ABUNDANCE_FILENAME'] = '"{}"'.format(afn)

        if self.set_keyword_change_abundance:
            dd['SET_KEYWORD_CHANGE_ABUNDANCE'] \
                            = self.set_keyword_change_abundance
            if path <> None:
                cffn = os.path.join(path,\
                                os.path.split(self.change_fraction_filename)[1])
            else:
                cffn = self.change_fraction_filename
            dd['CHANGE_FRACTION_FILENAME'] = '"{}"'.format(cffn)

        if self.set_keyword_change_temperature:
            dd['SET_KEYWORD_CHANGE_TEMPERATURE'] \
                            = self.set_keyword_change_temperature
            if path <> None:
                tfn = os.path.join(path,\
                                os.path.split(self.new_temperature_filename)[1])
            else:
                tfn = self.new_temperature_filename
            dd['NEW_TEMPERATURE_FILENAME'] = '"{}"'.format(tfn)

        if self.starfile:
            dd['USE_STARFILE'] = 1
            if path <> None:
                sfn = os.path.join(path,os.path.split(self.starfile)[1])
                
            else:
                sfn = self.starfile
            dd['STARFILE'] = '"{}"'.format(sfn)
        
        if int(in_progress):
            dd['IN_PROGRESS'] = 1   

        return dd        
                         


    def setModelId(self,model_id):
        
        '''
        Set a model_id for the molecule, which identifies model_id for MLINE!
        
        @param model_id: The model_d to be associated with this molecule
        @type model_id: string
        
        '''
        
        self.__model_id = model_id
        


    def makeLabel(self):
        
        '''
        Return a short-hand label for this particular molecule, 
        taken from the molecule.dat file.
        
        @return: The label associated with this molecule
        @rtype: string
        
        '''
        
        return '\ {}\ '.format(self.molecule_plot)    
            


    def getModelId(self):
        
        '''
        Return the model_id associated with this molecule. None if not yet set.
        
        @return: the model id is returned.
        @rtype: string
        
        '''
        
        return self.__model_id
    


    def isMolecule(self):
        
        '''
        Return True to help the codes know that this is a molecule 
        and not a transition.
        
        @return: True if molecule or False if transition
        @rtype: bool
        
        '''
        
        return True
    


    def getAbunFactor(self):
        
        '''
        Return the abundance factor of the molecule, with respect to its main 
        isotope/ortho version.
        
        '''
        
        raw_factors = DataIO.getInputData(keyword='ABUN_FACTOR',\
                                          filename='Molecule.dat',\
                                          rindex=self.molecule_index)
        factors = [factor == '1' and 1 or float(getattr(self,factor.lower())) 
                   for factor in raw_factors.split('*')]
        #-- abun_factor should only take into account isotope/OPR factors
        #if float(self.enhance_abundance_factor): 
        #    factors.append(float(self.enhance_abundance_factor))
        total_factor = 1
        while factors:
            total_factor *= factors.pop()
        return total_factor



    def isWater(self):
         
        '''
        Is this molecule a water molecule?
        
        (ortho, para, isotopologue thereof included)
        
        @return: True or False
        @rtype: bool
        
        '''
        
        return self.molecule in ['1H1H16O','p1H1H16O','1H1H17O',\
                                 'p1H1H17O','1H1H18O','p1H1H18O']
        
        
    
    def readMline(self):
    
        """
        Read the mline output for this molecule, given that a model_id is 
        available. 
        
        The mline output is available in the MlineReader object, as a property
        of Molecule().
        
        [NYI]
        
        
        """