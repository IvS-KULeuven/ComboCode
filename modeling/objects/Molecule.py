# -*- coding: utf-8 -*-

"""
Toolbox for Molecules, used in various applications for GASTRoNOoM.

Author: R. Lombaert

"""

import os 
import re

from cc.tools.io import DataIO
from cc.tools.io import Database
from cc.tools.io import Radiat



def makeMoleculeFromDb(molec_id,molecule,path_gastronoom='codeSep2010',\
                       path_combocode=os.path.join(os.path.expanduser('~'),\
                                                   'ComboCode')):
    
    '''
    Make a Molecule() from a database, based on model_id and molec_id.
    
    @param molec_id: the model_id of the molecule
    @type molec_id: string
    @param molecule: the short hand name of the molecule        
    @type molecule: string
    @keyword path_gastronoom: the output path in the ~/GASTRoNOoM/. directory
    
                              (default: codeSep2010)
    @type path_gastronoom: string
    @keyword path_combocode: the path to the ComboCode home directory
    
                             (default: ~/ComboCode/)
    @type path_combocode: 
    @return: the molecule with all its information inbedded
    @rtype: Molecule()
    
    '''
    
    filepath = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                            path_gastronoom)
    if os.path.isfile(os.path.join(filepath,'models',molec_id,\
                                   'cooling_id.log')):
        model_id = DataIO.readFile(os.path.join(filepath,'models',molec_id,\
                                                'cooling_id.log'))[0]
    #- ie mline id is the same as model id, the first calced for this id
    else: 
        model_id = molec_id
    molec_db = Database.Database(os.path.join(filepath,\
                                              'GASTRoNOoM_mline_models.db'))
    molec_dict = molec_db[model_id][molec_id][molecule].copy()
    extra_pars = molec_dict['MOLECULE'].split()[1:]
    for k,v in zip(['ny_low','ny_up','nline','n_impact','n_impact_extra'],\
                   extra_pars):
        molec_dict[k] = int(v)
    for key in ['MOLECULE','CHANGE_DUST_TO_GAS_FOR_ML_SP',\
                'NUMBER_INPUT_ABUNDANCE_VALUES','KEYWORD_TABLE',\
                'MOLECULE_TABLE']:
        if molec_dict.has_key(key):
            del molec_dict[key]
    molec_dict = dict([(k.lower(),v) for k,v in molec_dict.items()])
    molec = Molecule(molecule=molecule,path_combocode=path_combocode,\
                     **molec_dict)
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
                 ratio_12c_to_13c=0,ratio_16o_to_17o=0,ratio_16o_to_18o=0,\
                 opr=0,r_outer=0,outer_r_mode='MAMON',abundance_filename=None,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),\
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
        @keyword path_combocode: CC home folder
          
                                 (default: '/home/robinl/ComboCode')
        @type path_combocode: string       
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
        @keyword linelist: The molecule is created for the LineList module. No
                           radiative information is read from GASTRoNOoM input
                           files.
        
                           (default: 0)
        @type linelist: bool
        @keyword starfile: input filename for a stellar input spectrum (either
                           user defined or from a model atmosphere spectrum)
                           Note that the file must be linked to
                           ~/ComboCode/StarFiles/starfile_tablestar.dat
                           for it to be used by GASTRoNOoM, and that the code 
                           itself must have been compiled with the correct 
                           b.f and table_star.f files. 
                           
                           (default: '')
        @type starfile: str
                           
        
        '''
 
        self.molecule = str(molecule)
        self.ny_up = int(ny_up)
        self.ny_low = int(ny_low)
        self.nline = int(nline)
        self.n_impact = int(n_impact)
        self.n_impact_extra = int(n_impact_extra)
        self.path_combocode = path_combocode
        self.molecule_index = DataIO.getInputData(keyword='TYPE_SHORT',\
                                path=os.path.join(self.path_combocode,'Data'),\
                                filename='Molecule.dat').index(self.molecule)
        self.molecule_full = DataIO.getInputData(keyword='MOLEC_TYPE',\
                                path=os.path.join(self.path_combocode,'Data'),\
                                filename='Molecule.dat',make_float=0)\
                             [self.molecule_index]
        self.molecule_short = DataIO.getInputData(keyword='NAME_SHORT',\
                                path=os.path.join(self.path_combocode,'Data'),\
                                filename='Molecule.dat')[self.molecule_index]
        self.molecule_plot = DataIO.getInputData(keyword='NAME_PLOT',\
                                path=os.path.join(self.path_combocode,'Data'),\
                                filename='Molecule.dat')[self.molecule_index]
        self.spec_indices = DataIO.getInputData(keyword='SPEC_INDICES',\
                                path=os.path.join(self.path_combocode,'Data'),\
                                filename='Molecule.dat')[self.molecule_index]
        self.use_indices_dat = DataIO.getInputData(keyword='USE_INDICES_DAT',\
                                path=os.path.join(self.path_combocode,'Data'),\
                                filename='Molecule.dat')[self.molecule_index]
        if self.use_indices_dat and not linelist:
             self.indices_index = DataIO.getInputData(start_index=4,\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='MOLECULE',filename='Indices.dat')\
                                .index('_'.join([self.molecule,\
                                                 str(self.ny_low),\
                                                 str(self.ny_up),\
                                                 str(self.nline)]))
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
        self.use_collis_radiat_switch = use_collis_radiat_switch
        self.ratio_12c_to_13c = ratio_12c_to_13c
        self.ratio_16o_to_17o = ratio_16o_to_17o
        self.ratio_16o_to_18o = ratio_16o_to_18o
        self.opr = opr
        self.dust_to_gas_change_ml_sp = float(dust_to_gas_change_ml_sp)
        self.enhance_abundance_factor = float(enhance_abundance_factor)
        if abundance_filename <> None:
            self.abundance_filename = os.path.split(abundance_filename)[0] \
                                        and abundance_filename \
                                        or os.path.join(os.path.expanduser('~'),\
                                                        'GASTRoNOoM',\
                                                        'CustomAbundances',\
                                                        abundance_filename)
        #- Mainly for plotting purposes: The relative, multiplicative abundance 
        #- factor with respect to main isotope is calculated for isotopologue
        self.abun_factor = self.setAbunFactor() 
        self.outer_r_mode = outer_r_mode
        if self.outer_r_mode == 'MAMON': self.r_outer = 0
        else: self.r_outer = float(r_outer)
        self.set_keyword_change_abundance = int(set_keyword_change_abundance)
        self.change_fraction_filename = change_fraction_filename
        self.set_keyword_change_temperature = \
                    int(set_keyword_change_temperature)
        self.new_temperature_filename = new_temperature_filename
        self.__model_id = None
        if not linelist:
            self.radiat = Radiat.Radiat(molecule=self,\
                                    use_indices_dat=self.use_indices_dat)
            if self.spec_indices:
                path = os.path.join(os.path.expanduser('~'),'GASTRoNOoM','src',\
                                    'data')
                if self.use_indices_dat:
                    f = DataIO.getInputData(path=os.path.join(self.path_combocode,\
                                                            'Data'),\
                                        keyword='INDICES',\
                                        filename='Indices.dat',\
                                        start_index=4)[self.indices_index]
                    filename = os.path.join(path,'indices_backup',f)
                else:
                    filename = os.path.join(path,'%s_indices.dat'%self.molecule)
                self.radiat_indices = [[int(i) for i in line] 
                                     for line in DataIO.readFile(filename,' ')]
        else:
            self.radiat = None
            self.radiat_indices = None
        self.starfile = starfile


    def __str__(self):
        
        '''
        Printing a molecule as it should appear in the GASTRoNOoM input file.
        
        '''
        
        return 'MOLECULE=%s %i %i %i %i %i' %(self.molecule_full,self.ny_low,\
                                             self.ny_up,self.nline,\
                                             self.n_impact,self.n_impact_extra)



    def updateParameters(self,pardict):
        
        '''
        Update parameters.
        
        @param pardict: the parameters with respective values for the update
        @type pardict: dict()
        
        '''
        
        for k,v in pardict.items():
            if hasattr(self,k.lower()): setattr(self,k.lower(),v)



    def makeDict(self,path=None):
        
        '''
        Return a dict with molecule string, and other relevant parameters.
        
        @keyword path: If a different path is needed, it can be passed here, 
                            for files
                      
                            (default: None)
        @type path: string
        
        '''
        
        new_dict = dict([('MOLECULE',str(self).replace('MOLECULE=','')),\
                         ('ITERA',self.itera),\
                         ('ABUN_MOLEC',self.abun_molec),\
                         ('ABUN_MOLEC_RINNER',self.abun_molec_rinner),\
                         ('ABUN_MOLEC_RE', self.abun_molec_re),\
                         ('RMAX_MOLEC',self.rmax_molec),\
                         ('LTE_REQUEST',self.lte_request),\
                         ('OUTER_R_MODE',self.outer_r_mode),\
                         ('USE_COLLIS_RADIAT_SWITCH',self.use_collis_radiat_switch),\
                         ('R_OUTER',self.r_outer)])
        for par,isot in [('RATIO_16O_TO_18O','18O'),\
                         ('RATIO_16O_TO_17O','17O'),\
                         ('RATIO_12C_TO_13C','13C'),\
                         ('OPR','p1H')]:
            if isot in self.molecule:
                new_dict[par] = getattr(self,par.lower())
        if self.dust_to_gas_change_ml_sp:
            new_dict['CHANGE_DUST_TO_GAS_FOR_ML_SP'] \
                    = self.dust_to_gas_change_ml_sp and 1 or 0
            new_dict['DUST_TO_GAS_CHANGE_ML_SP'] \
                    = self.dust_to_gas_change_ml_sp
        if self.enhance_abundance_factor:
            new_dict['ENHANCE_ABUNDANCE_FACTOR'] \
                    = self.enhance_abundance_factor
            new_dict['NUMBER_INPUT_ABUNDANCE_VALUES'] \
                    = len(DataIO.readCols(filename=self.abundance_filename)[0])
            new_dict['KEYWORD_TABLE'] = 1
            new_dict['MOLECULE_TABLE'] \
                    = self.molecule_full[:self.molecule_full.index('.')]
            new_dict['ISOTOPE_TABLE'] = self.molecule_full
            if path <> None:
                 new_dict['ABUNDANCE_FILENAME'] \
                    = '"%s"'%os.path.join(path,\
                             os.path.split(self.abundance_filename)[1])
            else:
                 new_dict['ABUNDANCE_FILENAME'] \
                    = '"%s"'%self.abundance_filename
        if self.set_keyword_change_abundance:
            new_dict['SET_KEYWORD_CHANGE_ABUNDANCE'] \
                    = self.set_keyword_change_abundance
            if path <> None:
                 new_dict['CHANGE_FRACTION_FILENAME'] \
                    = '"%s"'%os.path.join(path,\
                             os.path.split(self.change_fraction_filename)[1])
            else:
                 new_dict['CHANGE_FRACTION_FILENAME'] \
                    = '"%s"'%self.change_fraction_filename
        if self.set_keyword_change_temperature:
            new_dict['SET_KEYWORD_CHANGE_TEMPERATURE'] \
                    = self.set_keyword_change_temperature
            if path <> None:
                 new_dict['NEW_TEMPERATURE_FILENAME'] \
                    = '"%s"'%os.path.join(path,\
                             os.path.split(self.new_temperature_filename)[1])
            else:
                 new_dict['NEW_TEMPERATURE_FILENAME'] \
                    = '"%s"'%self.new_temperature_filename
        if self.starfile:
            new_dict['STARFILE'] = self.starfile
        return new_dict        
                         


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
        
        return '\ %s\ '%self.molecule_plot    
            


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
    


    def setAbunFactor(self):
        
        '''
        Return the abundance factor of the molecule, with respect to its main 
        isotope/ortho version.
        
        '''
        
        factors = [factor == '1' and 1 or float(getattr(self,factor.lower())) 
                   for factor in DataIO.getInputData(\
                                        path=os.path.join(self.path_combocode,\
                                                          'Data'),\
                                        keyword='ABUN_FACTOR',\
                                        filename='Molecule.dat')\
                                 [self.molecule_index].split('*')]
        if float(self.enhance_abundance_factor): 
            factors.append(float(self.enhance_abundance_factor))
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
        