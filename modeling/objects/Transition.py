# -*- coding: utf-8 -*-

"""
Toolbox for Transitions, used in various applications concerning GASTRoNOoM.

Author: R. Lombaert

"""

import os 
import re
import scipy
from scipy import pi, exp, linspace, argmin, array, diff, mean
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import types

from ivs.sigproc import filtering

from cc.data import LPTools
from cc.modeling.objects import Molecule 
from cc.tools.io import Database, DataIO
from cc.tools.io import SphinxReader
from cc.tools.io import FitsReader, TxtReader
from cc.tools.numerical import Interpol
from cc.statistics import BasicStats as bs


def extractTransFromStars(star_grid,sort_freq=1,sort_molec=1,pacs=0):
    
    '''
    Extract a list a of unique transitions included in a list of Star() objects
    and sort them.
    
    An extra flag is added to the input keyword to determine inclusion of PACS
    lines or not.
    
    The list of transitions is copies to make sure no funky references mess 
    with the original transitions.
    
    @param star_grid: The list of Star() objects from which all 'GAS_LINES' 
                       keys are taken and collected in a set. 
    @type star_grid: list[Star()]
    
    @keyword sort_freq: Sort the transitions according to their frequencies. 
                        Otherwise, they are sorted by their wavelengths.
    
                        (default: 1)
    @type sort_freq: bool
    @keyword sort_molec: Sort the transitions by molecule first.
                      
                         (default: 1) 
    @type sort_molec: bool
    @keyword pacs: 1: Only select PACS lines. 0: Don't select PACS lines. 2: Do
                      not make a distinction between the telescope type.
                      
                      (default: 0)
    @type pacs: int
    
    @return: a list of unique transitions included in all Star() objects in
             star_grid
    @rtype: list[Transition()]
    
    '''
    
    pacs = int(pacs)
    selection = list(sorted(set([trans 
                                 for star in star_grid
                                 for trans in star['GAS_LINES']]),\
                            key=lambda x:sort_freq \
                                and (sort_molec and x.molecule or '',\
                                     x.frequency) \
                                or  (sort_molec and x.molecule or '',\
                                     x.wavelength)))
    if not pacs:
        return [trans for trans in selection if 'PACS' not in trans.telescope]
    if pacs == 1:
        return [trans for trans in selection if 'PACS' in trans.telescope]
    else:
        return selection



def updateLineSpec(trans_list):
    
    '''
    Update telescope.spec file.
    
    This method checks for requested transitions that are not yet present in 
    .spec file of the telescope, and adds them.
    
    @param trans_list: The transitions in the CC input
    @type trans_list: list[Transition()]
    
    '''
    
    telescopes = list(set([trans.telescope for trans in trans_list]))
    for telescope in telescopes:
        if not ('HIFI' in telescope or 'PACS' in telescope):
            print 'Warning! Telescope beam efficiencies for %s'%telescope + \
                  ' are added arbitrarily and thus are not the correct values.'
        old_spec = DataIO.readFile(os.path.join(os.path.expanduser('~'),\
                                  'GASTRoNOoM','src','data',telescope+'.spec'))
        line_spec_list = [line 
                          for line in old_spec 
                          if line.find('LINE_SPEC') >-1 and line[0] != '#']
        these_trans = [trans 
                       for trans in trans_list 
                       if trans.telescope == telescope]
        these_trans = [trans.getLineSpec() 
                       for trans in these_trans 
                       if trans.getLineSpec() not in line_spec_list]
        line_spec_list = sorted(line_spec_list + these_trans,\
                              key=lambda x: [x.split()[0],float(x.split()[9])])
        line_spec_list = [line.replace(' ','\t') for line in line_spec_list]
        try:
            line_spec_index = [line[0:9] 
                               for line in old_spec].index('LINE_SPEC')
            new_spec = old_spec[0:line_spec_index] + line_spec_list
        except ValueError:
            new_spec = old_spec + line_spec_list
        DataIO.writeFile(os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                      'src','data',telescope+'.spec'),\
                         new_spec+['\n######################################'])

 

def makeTransition(trans,star=None,def_molecs=None,\
                   path_combocode=os.path.join(os.path.expanduser('~'),\
                                               'ComboCode')):
    
    '''
    Create a Transition instance based on a Star object and a standard CC input 
    line, with 12 entries in a list.
    
    If a star object is not given, the method creates Molecule() objects itself
    For now only 12C16O, 1H1H16O and p1H1H16O are possible.
    
    @param trans: the input line with 12 entries, the first being the molecule,
                  followed by all 11 CC input parameters for a transition
    @type trans: list[string]
    
    @keyword star: The star object providing basic info
    
                   (default: None)
    @type star: Star()
    @keyword def_molecs: Default molecules needed for the requested transitions
                         None is returned if the requested molecule is not 
                         present. If both this and star are not given, a few
                         default molecules are loaded. 
                         
                         (default: None)
    @type def_molecs: dict(string: Molecule())
    @keyword path_combocode: CC home folder
          
                             (default: '/home/robinl/ComboCode')
    @type path_combocode: string       
        
    @return: The transition object is returned with all info included
    @rtype: Transition()
    
    '''
    
    if star is None and def_molecs is None: 
        def_molecs = {'12C16O':Molecule.Molecule('12C16O',61,61,240,\
                                               path_combocode=path_combocode),\
                      '1H1H16O':Molecule.Molecule('1H1H16O',39,90,1157,\
                                               path_combocode=path_combocode),\
                      'p1H1H16O':Molecule.Molecule('p1H1H16O',32,90,1029,\
                                               path_combocode=path_combocode)}
        molec = def_molecs.get(trans[0].replace('TRANSITION=',''),None)
        path_gastronoom = None
        umis = 0
    elif star <> None:
        molec = star.getMolecule(trans[0].replace('TRANSITION=',''))
        umis = star['USE_MASER_IN_SPHINX']
        path_combocode = star.path_combocode
        path_gastronoom = star.path_gastronoom
    else:
        molec = def_molecs.get(trans[0].replace('TRANSITION=',''),None)
        path_gastronoom = None
        umis = 0
    
    if molec <> None:
        return Transition(molecule=molec,\
                          vup=int(trans[1]),jup=int(trans[2]),\
                          kaup=int(trans[3]),kcup=int(trans[4]),\
                          vlow=int(trans[5]),jlow=int(trans[6]),\
                          kalow=int(trans[7]),kclow=int(trans[8]),\
                          telescope=trans[9],offset=float(trans[10]),\
                          n_quad=int(trans[11]),use_maser_in_sphinx=umis,\
                          path_combocode=path_combocode,\
                          path_gastronoom=path_gastronoom)
    else:
        return None    



def makeTransitionsFromTransList(filename,path_combocode=os.path.join(os.path\
                                            .expanduser('~'),'ComboCode')):
    
    '''
    Make Transition objects for the transitions listed in a line list file,
    used by CC. 
    
    The syntax is the same as the input Transition lines in a CC inputfile.
    
    @param filename: The filename to the linelist
    @type filename: string
    
    @keyword path_combocode: The CC home folder
    
                             (default: ~/ComboCode/)
    @type path_combocode: string
    
    @return: The Transitions in the file are returned in object form.
    @rtype: list[Transition]
    
    '''
    
    def_molecs = {'12C16O':Molecule.Molecule('12C16O',61,61,240,\
                                             path_combocode=path_combocode),\
                  '1H1H16O':Molecule.Molecule('1H1H16O',39,90,1157,\
                                              path_combocode=path_combocode),\
                  'p1H1H16O':Molecule.Molecule('p1H1H16O',32,90,1029,\
                                               path_combocode=path_combocode)}
    trl = DataIO.readDict(filename,multi_keys=['TRANSITION'])
    trl_sorted = DataIO.checkEntryInfo(trl['TRANSITION'],14,'TRANSITION')
    trans = [makeTransition(trans=t,def_molecs=def_molecs,\
                            path_combocode=path_combocode) 
             for t in trl_sorted]
    return trans

    

def makeTransitionFromSphinx(filename,path_combocode=os.path.join(os.path\
                                            .expanduser('~'),'ComboCode')):
    
    '''
    Make a Transition() based on the filename of a Sphinx file.
    
    For this information is taken from the Molecule() model database.
    
    @param filename: The sphinx file name, including path
    @type filename: string
    
    @keyword path_combocode: The CC home folder
    
                             (default: ~/ComboCode/)
    @type path_combocode: string
    
    @return: The transition
    @rtype: Transition()
    
    '''
    
    filepath,filename = os.path.split(filename)
    if not filepath:
        raise IOError('Please include the full filepath of the Sphinx ' + \
                      'filename, needed to determine the name of the database.')
    
    #- clip model_id and save it into trans_id
    filepath,trans_id = os.path.split(filepath)
    #- clip 'models'
    filepath = os.path.split(filepath)[0]
    #- clip path_gastronoom and save it in path_gastronoom
    filepath,path_gastronoom = os.path.split(filepath)
    
    if os.path.isfile(os.path.join(filepath,path_gastronoom,'models',trans_id,\
                                   'cooling_id.log')):
        model_id = DataIO.readFile(os.path.join(filepath,path_gastronoom,\
                                                'models',trans_id,\
                                                'cooling_id.log'))[0]
        #- If an mline_id.log exists, a cooling_id.log will always exist also
        if os.path.isfile(os.path.join(filepath,path_gastronoom,'models',\
                                       trans_id,'mline_id.log')):
            molec_id = DataIO.readFile(os.path.join(filepath,path_gastronoom,\
                                                    'models',trans_id,\
                                                    'mline_id.log'))[0]
        else: #- ie trans id is the same as molec id, first trans calced for id
            molec_id = trans_id
    #- ie mline and trans id are same as model id, first calced for id
    else: 
        model_id = trans_id
        molec_id = trans_id
        
    filename = filename.rstrip('.dat').split('_')[2:]
    molec_name = filename.pop(0)
    molec = Molecule.makeMoleculeFromDb(molecule=molec_name,molec_id=molec_id,\
                                        path_gastronoom=path_gastronoom,\
                                        path_combocode = path_combocode)
    telescope = filename.pop(-2)
    pattern = re.compile(r'^(\D+)(\d*.?\d*)$')
    numbers = [pattern.search(string).groups() for string in filename]
    numbers = dict([(s.lower(),n) for s,n in numbers])
    trans = Transition(molecule=molec,telescope=telescope,\
                       path_combocode=path_combocode,\
                       path_gastronoom=path_gastronoom,**numbers)
    trans.setModelId(trans_id)
    trans_db = Database.Database(os.path.join(filepath,path_gastronoom,\
                                              'GASTRoNOoM_sphinx_models.db'))
    trans_dict = trans_db[model_id][molec_id][trans_id][str(trans)].copy()
    [setattr(trans,k.lower(),v)
            for k,v in trans_dict.items()
            if k != 'TRANSITION']
    return trans



def checkUniqueness(trans_list):
    
    '''
    Check uniqueness of a list of transitions.
    
    Same transitions are replaced by a single transition with all datafiles 
    from the originals. 
    
    Based on the parameters of the transitions. If they are the same, only one
    transition is added to the output list, but the datafiles are all included.
    
    Datafiles do not identify a transition!
    
    @param trans_list: The transitions to be checked for uniqueness
    @type trans_list: list[Transition()]
    
    @return: The merged transitions with all datafiles included.
    @rtype: tuple[Transition()]
    
    '''
    
    merged = []
    for trans in trans_list:
        if trans not in merged: 
            merged.append(trans)
        else:
            merged[merged.index(trans)].addDatafile(trans.datafiles)
    return merged
    
    

class Transition():
    
    '''
    A class to deal with transitions in GASTRoNOoM.
    
    '''
    
    def __init__(self,molecule,telescope=None,vup=0,jup=0,kaup=0,kcup=0,\
                 nup=None,vlow=0,jlow=0,kalow=0,kclow=0,nlow=None,offset=0.0,\
                 frequency=None,exc_energy=None,int_intensity_log=None,\
                 n_quad=100,use_maser_in_sphinx=0,\
                 vibrational='',path_gastronoom=None,datafiles=None,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):
        
        '''
        
        Initiate a Transition instance, setting all values for the 
        quantummechanical parameters to zero (by default).
        
        @param molecule: The molecule to which this transition belongs
        @type molecule: Molecule()
        
        @keyword telescope: The telescope with which the transition is observed
        
                            (default: None)
        @type telescope: string
        @keyword vup: The upper vibrational level
        
                      (default: 0)
        @type vup: int
        @keyword jup: The upper rotational level
        
                      (default: 0)
        @type jup: int
        @keyword kaup: The upper level of the first projection of the ang mom
        
                       (default: 0)
        @type kaup: int       
        @keyword kcup: The upper level of the second projection of the ang mom
        
                       (default: 0)
        @type kcup: int
        @keyword nup: if not None it is equal to Kaup and only relevant for 
                      SO/hcn-type molecules
        
                      (default: 0)
        @type nup: int
        @keyword vlow: The lower vibrational level
        
                       (default: 0)
        @type vlow: int
        @keyword jlow: The lower rotational level
        
                       (default: 0)
        @type jlow: int
        @keyword kalow: The lower level of the first projection of the ang mom
        
                        (default: 0)
        @type kalow: int
        @keyword kclow: The lower level of the second projection of the ang mom
        
                        (default: 0)
        @type kclow: int
        @keyword nlow: if not None it is equal to Kalow, and only relevant for 
                       SO-type molecules
        
                       (default: None)
        @type nlow: int
        @keyword offset: The offset of the radiation peak with respect to the 
                         central pixel of the observing instrument. Only 
                         relevant for non-point corrected PACS or non Pacs data
                         (so not relevant when INSTRUMENT_INTRINSIC=1)
                         
                         (default: 0.0)
        @type offset: float
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        @keyword n_quad: Number of impact par. in quadrature int(Inu pdp).
                         Relevant for the ray-tracing of the line profiles
        
                         (default: 100)
        @type n_quad: int
        @keyword use_maser_in_spinx: using maser calc in sphinx code
                                     
                                     (default: 0)
        @type use_maser_in_spinx: bool
        @keyword frequency: if not None the frequency of the transition is 
                            taken to be this parameter in Hz, if None the 
                            frequency of this transition is derived from the 
                            GASTRoNOoM radiative input file
                            
                            No indices are searched for if frequency is given 
                            here. Usually used only for line listing.
                            
                            (default: None)
        @type frequency: float
        @keyword exc_energy: the excitation energy (cm-1, from CDMS/JPL) if you
                             want it included in the Transition() object, 
                             not mandatory!
                             
                             (default: None)
        @type exc_energy: float                     
        @keyword int_intensity_log: the integrated intensity of the line in 
                                    logscale from CDMS/JPL, if you want it 
                                    included in the Transition() object, 
                                    not mandatory!
                                    
                                    (default: None)
        @type int_intensity_log: float
        @keyword vibrational: In the case of line lists, this keyword indicates
                              the type of vibrational excitation, relevant fi 
                              for H2O. Empty string if vibrational groundstate
                              
                              (default: '')
        @type vibrational: string
        @keyword path_gastronoom: model output folder in the GASTRoNOoM home
        
                                  (default: None)
        @type path_gastronoom: string
        @keyword datafiles: (multiple) filename(s) and path(s) to a datafile 
                            for transition, specific for the telescope used for 
                            this dataset. Only applicable for single transition 
                            files, such as ground-based data. None if no files 
                            available.
                                    
                            (default: None)
        @type datafiles: string/list
        
        '''
        
        self.molecule = molecule
        if telescope is None:
            self.telescope = 'N.A.'
        else:
            telescope = telescope.upper()
            if telescope.find('H2O') != -1 and telescope.find('PACS') != -1 \
                    and not self.molecule.isWater():
                telescope = telescope.replace('-H2O','')
            elif telescope.find('H2O') == -1 and telescope.find('PACS') != -1 \
                    and self.molecule.isWater():
                print 'WARNING! Water lines should not be included in the ' + \
                      '%s.spec file. Create a file with the '%telescope + \
                      'same name, appending -H2O to the telescope name and ' +\
                      'removing all LINE_SPEC lines.'
                telescope = '%s-H2O'%telescope
            self.telescope = telescope
        self.vup = int(vup)
        self.jup = int(jup)
        self.kaup = int(kaup)
        self.kcup = int(kcup)
        self.vlow = int(vlow)
        self.jlow = int(jlow)
        self.kalow = int(kalow)
        self.kclow = int(kclow)
        self.offset = float(offset)
        self.path_combocode = path_combocode
        self.n_quad = int(n_quad)
        self.use_maser_in_sphinx = int(use_maser_in_sphinx)
        self.__model_id = None
        if nup is None or nlow is None:
            self.nup = self.kaup            
            self.nlow = self.kalow
            #- In case of SO, the quantum number is treated as Kaup/low but is
            #- in fact Nup/low
        self.exc_energy = exc_energy
        self.int_intensity_log = int_intensity_log
        self.vibrational = vibrational
        self.sphinx = None
        self.path_gastronoom = path_gastronoom
        if type(datafiles) is types.StringType:
            self.datafiles = [datafiles]
        else:
            self.datafiles = datafiles
        self.lpdata = None 
        self.radiat_trans = None
        if frequency is None:
             #-- sets frequency from GASTRoNOoM input in s^-1
            self.__setIndices()  
        else:
            self.frequency = frequency
        self.c = 2.99792458e10 #in cm
        self.wavelength = self.c/self.frequency #in cm
        self.best_vlsr = None
        self.fittedlprof = None
        self.best_mfilter = None
        self.intintpacs = dict()
        self.intinterrpacs = dict()
        self.intintpacs_blends = dict()
        
    def __str__(self):
        
        '''
        Printing a transition as it should appear in the GASTRoNOoM input file.
        
        @return: The Transition string
        @rtype: string
        
        '''
        
        return 'TRANSITION=%s %i %i %i %i %i %i %i %i %s %.2f' \
               %(self.molecule.molecule_full,self.vup,self.jup,self.kaup,\
                 self.kcup,self.vlow,self.jlow,self.kalow,self.kclow,\
                 self.telescope,self.offset)
                    


    def __eq__(self,other):
        
        '''
        Compare two transitions and return true if equal.
        
        The condition is the dictionary of this Transition() returned by the 
        makeDict() method.
        
        @return: The comparison
        @rtype: bool
        
        '''
        
        try:        
            if self.makeDict() == other.makeDict():
                return True
            else:
                return False
        except AttributeError:
            return False
                


    def __ne__(self,other):
        
        '''
        Compare two transitions and return true if not equal.
        
        The condition is the dictionary of this Transition() returned by the 
        makeDict() method.
        
        @return: The negative comparison
        @rtype: bool
        
        '''
        
        try:
            if self.makeDict() != other.makeDict():
                return True
            else:
                return False
        except AttributeError:
            return True


 
    def __hash__(self):
        
        '''
        Return a hash number based on the string of the transition, but make 
        sure that every comparison includes the isTransition method of the 
        trans.
        
        @return: The hash number:
        @rtype: int
        
        '''
        
        return hash(str(self.makeDict()))



    def getInputString(self):
         
         '''
         Return a string giving the CC input line for this transition.
         
         This includes N_QUAD and has the shorthand naming convention of the 
         molecule as opposed to the str() method.
         
         @return: the input line for this transition
         @rtype: string
         
         '''
         
         return 'TRANSITION=%s %i %i %i %i %i %i %i %i %s %.2f %i' \
                %(self.molecule.molecule,self.vup,self.jup,self.kaup,\
                  self.kcup,self.vlow,self.jlow,self.kalow,self.kclow,\
                  self.telescope,self.offset,self.n_quad)
          

    
    def getLineSpec(self):
        
        '''
        Make a string that suits telescope.spec inputfile.
        
        @return: The Line spec string
        @rtype: string
        
        '''
        
        return 'LINE_SPEC=%s %i %i %i %i %i %i %i %i %.2f %.2f ! %s'\
               %(self.molecule.molecule_full,self.vup,self.jup,self.kaup,\
                 self.kcup,self.vlow,self.jlow,self.kalow,self.kclow,\
                self.getBeamwidth(),self.getEfficiency(),\
                self.molecule.molecule)



    def getBeamwidth(self):
        
        '''
        Calculate the beamwidth specific for the telescope.
        
        @return: The beamwidth
        @rtype: float
        
        '''
        
        #- get telescope diameter in cm
        if 'PACS' in self.telescope.upper():
            wav = array([55.,68,73,75,84,94,110,136,145,150,168,187])/10000.
            beam = array([8.4875,8.19,8.19,8.4525,8.2875,8.45,9.0625,9.6875,\
                          9.96,10.5425,11.0125,11.96])
            interp = interp1d(wav,beam)
            try:
                return interp(self.wavelength)
            except ValueError:
                print 'WARNING! The wavelength of %s falls '%str(self) +\
                      'outside of the interpolation range when determining '+\
                      'the beamwidth. Extrapolating linearly...'
                if self.wavelength < wav[0]:
                    return Interpol.linInterpol(wav[:2],beam[:2],\
                                                self.wavelength)
                else:
                    return Interpol.linInterpol(wav[-2:],beam[-2:],\
                                                self.wavelength)
        filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM','src',\
                                'data',self.telescope+'.spec')
        telescope_diameter = [float(line.split('=')[1][0:line.split('=')[1]\
                                    .index('!')].strip())
                              for line in DataIO.readFile(filename)
                              if line.find('TELESCOPE_DIAM') == 0][0] * 100.
        
        #- 1.22 is diffraction limited specification, 
        #- last factor is conversion to arcseconds
        return 1.22*self.wavelength/telescope_diameter*60.*60.*360./(2.*pi)  
        


    def getEfficiency(self):
        
        '''
        Calculate telescope beam efficiency. 
        
        Telescope specific!
        
        This beam efficiency is included in the .spec files for GASTRoNOoM.
        
        This number however is not used in any calculations of the line profile 
        and is included here for reference only. 
        
        The numbers currently are correct only for HIFI.
        
        @return: The beam efficiency
        @rtype: float
        
        '''
        
        if self.telescope.find('PACS') > -1 \
                or self.telescope.find('SPIRE') > -1:
            #- some random number, it is NOT relevant for PACS or SPIRE...
            return 0.60                 
        if self.telescope.find('HIFI') > -1:
            if self.frequency > 1116.2e9 and self.frequency < 1400.e9: 
                #- band 5a and 5b
                eff_mb0 = 0.66
            else:
                eff_mb0 = 0.76
            eff_f = 0.96
            sigma = 3.8e-4  #in cm
            eff_mb = eff_mb0 * exp(-1*(4*pi*sigma/self.wavelength)**2)
            return eff_mb/eff_f
        else:
            return 1.00



    def makeDict(self,in_progress=0):
        
        '''
        Return a dict with transition string, and other relevant parameters.
        
        @keyword in_progress: add an extra dict entry "IN_PROGRESS" if the 
                              transition is still being calculated on Vic.
                              
                              (default: 0)
        @type in_progress: bool
        
        @return: The transition dictionary including all relevant, defining 
                 information
        @rtype: dict()
        
        '''
        
        if int(in_progress):
            return dict([('TRANSITION',str(self).replace('TRANSITION=','')),\
                         ('N_QUAD',self.n_quad),\
                         ('USE_MASER_IN_SPHINX',self.use_maser_in_sphinx),\
                         ('IN_PROGRESS',1)])
        else:
            return dict([('TRANSITION',str(self).replace('TRANSITION=','')),\
                         ('N_QUAD',self.n_quad),\
                         ('USE_MASER_IN_SPHINX',self.use_maser_in_sphinx)])
 


    def makeSphinxFilename(self,number='*'):
        
        '''
        Return a string in the sphinx filename format.
         
        @keyword number: the number in the filename (sph*, sph1 or sph2, hence 
                         can be *, 1, 2)
                         
                         (default: '*')
        @type number: string
        
        @return: The sphinx filename for this transition
        @rtype: string
        
        '''
        
        try:
            number = str(int(number))
        except ValueError:
            number = str(number)
        quantum_dict = dict()
        if not self.molecule.spec_indices:
            for i,attr in enumerate(['vup','jup','vlow','jlow']):
                quantum_dict[i] = (attr,getattr(self,attr))
        elif self.molecule.spec_indices == 2:
            for i,attr in enumerate(['vup','Jup','Nup','vlow','jlow','Nlow']):
                quantum_dict[i] = (attr,getattr(self,attr.lower()))
        elif self.molecule.spec_indices == 3:
            for i,attr in enumerate(['vup','jup','vlow','jlow']):
                quantum_dict[i] = (attr,getattr(self,attr))            
        else:
            for i,attr in enumerate(['vup','jup','Kaup','Kcup',\
                                     'vlow','jlow','Kalow','Kclow']):
                quantum_dict[i] = (attr,getattr(self,attr.lower()))
        return '_'.join(['sph%s%s'%(number,self.getModelId()),\
                                    self.molecule.molecule] + \
                        ['%s%i'%(quantum_dict[k][0],quantum_dict[k][1])
                           for k in sorted(quantum_dict.keys())] + \
                        [self.telescope,'OFFSET%.2f.dat'%(self.offset)])
        


    def getFrequency(self):
        
        '''
        Get frequency of the transition from the radiat.dat files. 
        
        @return: frequency in Hz
        @rtype: float
        
        '''

        return self.frequency
        


    def getEnergyUpper(self):
         
        '''
        Return the energy level of the upper state of this transition.

        @return: energy level in cm^-1
        @rtype: float
    
        '''
        
        if self.molecule.radiat is None:
            print '%s_radiat.dat not found. Cannot find energy levels.'\
                  %self.molecule.molecule
            return
        energy = self.molecule.radiat.getEnergyLevels()
        return  float(energy[self.up_i-1])



    def getEnergyLower(self):

        '''
        Return the energy level of the lower state of this transition.

        @return: energy level in cm^-1
        @rtype: float
    
        '''
        
        if self.molecule.radiat is None:
            print '%s_radiat.dat not found. Cannot find energy levels.'\
                  %self.molecule.molecule
            return
        energy = self.molecule.radiat.getEnergyLevels()
        return  float(energy[self.low_i-1])



    def __setIndices(self):
         
        '''
        Set the index of this transition in the radiat file of GASTRoNOoM.
        
        The index from the indices file for lower and upper state are set.
        
        '''
        #- For 12C16O and 13C16O:
        #- indices = [i<60 and [i+1,0,i] or [i+1,1,i-60] for i in range(120)]
        #- As shown in above line, the first 60 (0-59) j's are associated with
        #- index (1-60) and v=0, the next 60 (60-119) j's are associated with 
        #- index (61-120) and v=1 

        if not self.molecule.spec_indices:
            self.up_i = self.jup + self.vup*self.molecule.ny_low + 1
            self.low_i = self.jlow + self.vlow*self.molecule.ny_low + 1
        else:
            indices = self.molecule.radiat_indices
            #- some molecules have only 2 or 3 relevant quantum numbers
            quantum_up = [q 
                            for i,q in enumerate([self.vup,self.jup,\
                                                    self.kaup,self.kcup]) 
                            if i<len(indices[0])-1]  
            quantum_low = [q 
                                for i,q in enumerate([self.vlow,self.jlow,\
                                                    self.kalow,self.kclow]) 
                                if i<len(indices[0])-1]
            #- Get index of the transition quantum numbers in the indices list
            #- If not present in list, ValueError is raised: probably caused 
            #- by using a linelist that doesn't include this transition.
            self.up_i  = indices[[i[1:] 
                                for i in indices].index(quantum_up)][0]
            self.low_i = indices[[i[1:] 
                                for i in indices].index(quantum_low)][0]
        self.radiat_trans = self.molecule.radiat.getTransInfo(low_i=self.low_i,\
                                                              up_i=self.up_i)
        if self.radiat_trans is False:
            raw_input('Something fishy is going on in Transition.py... '+\
                    'non-unique transition indices! Abort')
        self.frequency = float(self.radiat_trans['frequency'])
         

    def makeLabel(self):
        
        '''
        Return a short-hand label for this particular transition. 
        
        These labels can be used for plot line identifications for instance.
        
        If vibrational is not None, it always concerns a line list and is 
        included as well.
        
        @return: The transition label
        @rtype: string
        
        '''
        
        if self.vibrational:
            if not self.molecule.spec_indices:
                return '%s: %i,%i - %i,%i' \
                       %(self.vibrational,self.vup,self.jup,self.vlow,\
                         self.jlow)
            elif self.molecule.spec_indices == 2:
                return '%s: %i,%i$_{%i}$ - %i,%i$_{%i}$'\
                       %(self.vibrational,self.vup,self.jup,self.nup,\
                         self.vlow,self.jlow,self.nlow)
            elif self.molecule.spec_indices == 3:
                return '%s: %i$_{%i}$ - %i$_{%i}$' \
                       %(self.vibrational,self.jup,self.nup,self.jlow,\
                         self.nlow)
            else:
                return '%s: %i,%i$_{%i,%i}$ - %i,%i$_{%i,%i}$'\
                       %(self.vibrational,self.vup,self.jup,self.kaup,\
                         self.kcup,self.vlow,self.jlow,self.kalow,self.kclow)
        else:
            if not self.molecule.spec_indices:
                if self.vup == 0 and self.vlow ==0:
                    return r'$J=%i-%i$' %(self.jup,self.jlow)
                else:
                    return r'$\nu=%i$, $J=%i-%i' \
                           %(self.vup,self.jup,self.jlow)
            elif self.molecule.isWater():
                ugly = r'J_{\mathrm{K}_\mathrm{a}, \mathrm{K}_\mathrm{c}}'
                if self.vup == 0 and self.vlow ==0:
                    return r'$%s=%i_{%i,%i}-%i_{%i,%i}$'\
                            %(ugly,self.jup,self.kaup,self.kcup,\
                              self.jlow,self.kalow,self.kclow)
                else:
                    dvup = {1:2, 2:3}
                    label = r'$\nu_%i=1$, $%s=%i_{%i,%i}-%i_{%i,%i}$'\
                             %(dvup[self.vup],ugly,self.jup,self.kaup,\
                               self.kcup,self.jlow,self.kalow,self.kclow)
                    return label
            elif not self.molecule.isWater() and \
                    (self.molecule.spec_indices == 2 \
                     or self.molecule.spec_indices == 3):
                return '%i,%i$_{%i}$ - %i,%i$_{%i}$'\
                       %(self.vup,self.jup,self.nup,self.vlow,self.jlow,\
                         self.nlow)
            else:
                return '%i,%i$_{%i,%i}$ - %i,%i$_{%i,%i}$'\
                       %(self.vup,self.jup,self.kaup,self.kcup,\
                         self.vlow,self.jlow,self.kalow,self.kclow)
            


    def setModelId(self,model_id):
        
        '''
        Set a model_id for the transition, identifying the model for SPHINX! 
        
        May or may not be equal to the molecule's id, depending on the input 
        parameters.
        
        @param model_id: The model id
        @type model_id: string
        
        '''
        
        self.__model_id = model_id
        


    def getModelId(self):
        
        '''
        Return the model_id associated with this transition. 
        
        None if not yet set.
        
        Empty string if the model calculation failed.
        
        @return: The model id of this transition
        @rtype: string
        
        '''
        
        return self.__model_id
        

    
    def isMolecule(self):
        
        '''
        Is this a Molecule() object?
        
        @return: Molecule()?
        @rtype: bool
        
        '''
        
        return False



    def isTransition(self):
        
        '''
        Is this a Transition() object?
        
        @return: Transition()?
        @rtype: bool
        
        '''
        
        return True
    


    def isDone(self):
        
        '''
        Return True if successfully calculated by sphinx, False otherwise.
        
        @return: Finished calculation?
        @rtype: bool
        
        '''
        
        return self.__model_id
        
    
    
    def readSphinx(self):
         
        '''
        Read the sphinx output if the model id is not None or ''.
        
        '''
         
        if self.sphinx is None and self.getModelId() <> None \
                and self.getModelId() != '':
            filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                    self.path_gastronoom,'models',\
                                    self.getModelId(),\
                                    self.makeSphinxFilename())
            self.sphinx = SphinxReader.SphinxReader(filename)
     
     
     
    def addDatafile(self,datafile):
        
        '''
        Add a datafile name/multiple filenames for this transition. 
        
        If datafile has been given a valid value, the self.lpdata list is set
        back to None, so the datafiles can all be read again. 
        
        @param datafile: the full filename, or multiple filenames
        @type datafile: string/list
        
        '''
        
        if type(datafile) is types.StringType:
            if self.datafiles is None:
                self.datafiles = [datafile]
            else: 
                self.datafiles.append(datafile)
        else:
            if self.datafiles is None:
                self.datafiles = datafile
            else: 
                self.datafiles.extend(datafile)
        if datafile: self.lpdata = None
            

    def readData(self):
         
        '''
        Read the datafiles associated with this transition if available.
        
        '''
         
        if self.lpdata is None:
            self.lpdata = []
            if self.datafiles <> None:
                for df in self.datafiles:
                    info_path = os.path.join(self.path_combocode,'Data')
                    if df[-5:] == '.fits':
                        lprof = FitsReader.FitsReader(df,info_path)
                    else:
                        lprof = TxtReader.TxtReader(df,info_path)
                    self.lpdata.append(lprof)
            else:
                print 'No data found for %s. Setting v_lsr to 0.0'%str(self)+\
                      ' and cannot estimate noise or vexp values.'
                    
    
    
    def setData(self,trans):
        
        """
        Set data equal to the data in a different Transition object, but for 
        the same transition. Does no erase data if they have already been set 
        in this object. 
        
        A check is ran to see if the transition in both objects really is the 
        same. Sphinx model id may differ! 
        
        Avoids too much overhead when reading data from files. 
        
        The same is done for the fitresults from LPTools.fitLP().
        
        @param trans: The other Transition() object, assumes the transition 
                      is the same, but possibly different sphinx models.
        @type trans: Transition()        
        
        """
        
        if self.lpdata is None and self == trans: 
            self.lpdata = trans.lpdata
        if self.fittedlprof is None and self == trans: 
            self.fittedlprof = trans.fittedlprof
    
    
    def getVlsr(self):
        
        """
        Return the vlsr read from the fits file, or taken from the Star.dat file.
        
        This is different from the getBestVlsr() method, which determines the 
        best matching vlsr between data and sphinx, if both are available. 
        
        @return: the source velocity taken from the fits file OR Star.dat. 0
                 if data are not available.
        @rtype: float
        
        """
        
        self.readData()
        if self.lpdata: 
            return self.lpdata[0].getVlsr()
        else:
            return 0.0
            
    
    def fitLP(self):
        
        '''
        Run the autofit routine for line profiles. 
        
        The gas terminal velocity, its error, the soft parabola function and 
        possibly the extra gaussian will be set. It is possible the SP is a 
        gaussian instead, if a soft parabola did not work well. 
        
        The method makes sure the data have been read. The fitting routine is 
        only done for the first dataset in the list!
        
        Lastly, this method makes sure the noise is calculated for the first
        data object in this class.
        
        '''
        
        self.readData()
        if self.lpdata and self.fittedlprof is None:
            self.fittedlprof = \
              LPTools.fitLP(lprof=self.lpdata[0],\
                            info_path=os.path.join(self.path_combocode,'Data'))
            vexp = self.fittedlprof['vexp']
            self.lpdata[0].setNoise(vexp)
            
    
    
    def getNoise(self):
        
        '''
        Return the noise value of the FIRST data object in this transition, if
        available. 
        
        Note that None is returned if fitLP has not yet been ran, or if no data
        are available.
        
        @return: The noise value
        @rtype: float
        
        '''
        
        if self.lpdata:
            if self.lpdata[0].getNoise() is None:
                self.fitLP()
            return self.lpdata[0].getNoise()
        else:
            return None
            

    def getVexp(self):
        
        '''
        Get the gas terminal velocity as estimated from a line profile fitting
        routine. 
        
        @return: vexp
        @rtype: float
        
        '''
        
        if self.lpdata:
            if self.lpdata[0].getNoise() is None:
                self.fitLP()
            return self.fittedlprof['vexp']
        else:

            return None
  
            
    def getBestVlsr(self):
        
        """ 
        If self.best_vlsr is None, the best source velocity will be guessed by
        comparing chi^2 values between different values for the source velocity
        
        May be different from input value vlsr, the original expected 
        source velocity (from Star.dat)!
        
        Based on the first dataset in self.lpdata. Note that multiple datasets
        might be included, but the scaling will be done for only the first one. 
        Since different datasets are for the same telescope, the same line
        and (usually) the same conditions, the source velocity is not expected
        to be different for the different datasets anyway. 
        
        Method will attempt to call self.readData and readSphinx if either are 
        not available. If data are still not available, the initial guess is 
        returned. If data are available, but sphinx isn't, vlsr from the fits
        files is returned, and the initial guess is returned in case of txt 
        files.
       
        @return: the best guess vlsr, or the initial guess if no sphinx or data
                 are available [will return vlsr included in fitsfiles if 
                 available].
        @rtype: float
        
        """
        
        #-- check if vlsr was already calculated
        if self.best_vlsr <> None:
            return self.best_vlsr

        #-- attempt to read data and find the initial vlsr guess
        vlsr = self.getVlsr()
        #-- check if sphinx is finished, if not return vlsr from data container
        self.readSphinx()
        if vlsr == 0.0 or not self.sphinx: 
            return vlsr

        #-- Auto fit the line profile with a soft parabola and/or gaussian.
        #   This will set the vexp, evexp, soft parabola and gaussian profiles
        self.fitLP()
        vexp = self.fittedlprof['vexp']
        
        #-- get all the profiles and noise values
        noise = self.lpdata[0].getNoise()
        dvel = self.lpdata[0].getVelocity()
        dtmb = self.lpdata[0].getFlux()
        mvel = self.sphinx.getVelocity()
        mtmb = self.sphinx.getLPTmb()
        
        #-- Finding the best vlsr:
        #   1) filter the sphinx model onto the data grid, after rescaling the
        #      sphinx velocity grid to the given vlsr of the data.
        #   2) Check in the interval [vlsr-0.5vexp:vlsr+0.5*vexp] with steps 
        #      equial to the data bin size if there is a better match between
        #      model and data. This gives the 'best_vlsr'
        res = dvel[1]-dvel[0]
        mtmb_filter = filtering.filter_signal(x=mvel+vlsr,y=mtmb,ftype='box',\
                                              x_template=dvel,window_width=res)
        #-- Number of values tested is int(0.5*vexp/res+1),0.5*vexp on one side 
        #   and on the other side
        nstep = int(0.5*vexp/res+1)
        
        #-- Check if there are enough zeroes in the model flux grid
        #   ie if either of the following 2 statements evaluate to True, a non-
        #   zero element is found, and the technique used here for matching 
        #   different vlsr cannot be used
        if mtmb_filter[1][:nstep].any() or mtmb_filter[1][-nstep:].any():
            raise ValueError('Warning! Not enough zeroes in the grid! ' + \
                             'Talk to Robin!')
        #-- Since usually the data velocity grid extends far beyond what's 
        #   given by the sphinx velocity grid, we can just shift by adding and
        #   removing elements at the start and end of the list.
        mtmb_grid = [mtmb_filter[1]]
        
        for i in range(1,nstep):
            mtmbi_left = mtmb_filter[1][i:]
            while len(mtmbi_left) < len(dtmb):
                mtmbi_left = scipy.append(mtmbi_left,0)
            mtmbi_right = mtmb_filter[1][:-i]
            while len(mtmbi_right) < len(dtmb):
                mtmbi_right = scipy.insert(mtmbi_right,0,0)
            mtmb_grid.extend([mtmbi_left,mtmbi_right])
            
        #-- Calculate the chi squared for every filtered model
        chisquared = [bs.calcChiSquared(data=dtmb[dtmb>=-3*noise],\
                                        model=mtmbi[dtmb>=-3*noise],\
                                        noise=noise)
                      for mtmbi in mtmb_grid]
                      
        #-- Get the minimum chi squared, set the best_vlsr and set the best 
        #   filtered model profile
        self.chi2_best_vlsr = min(chisquared)
        #-- Tracing back the step associated with the index (uses int division)
        #   0 gives 0 again. 1 gives 1, 2 gives 1. 3 gives 2, 4 gives 2. etc.
        imin = argmin(chisquared)
        best_step = (imin-1)/2+1
        #-- Determining whether to add or subtract best_step*res from the vlsr
        modifier = imin%2 == 0 and 1 or -1
        self.best_vlsr = vlsr + modifier*best_step*res
        #-- Note that the velocity grid of best_mfilter is the data velocity
        self.best_mfilter = mtmb_grid[imin]
        
        #print "Best V_lsr: %f km/s, "%self.best_vlsr + \
              #"original V_lsr: %f km/s for transition %s, %s."\
              #%(vlsr,str(self),self.getModelId())
        return self.best_vlsr
        
        
   
    def getIntIntIntSphinx(self,units='si'):
        
        """
        Calculate the integrated intrinsic intensity of the sphinx line profile
        in SI or cgs units. Velocity is converted to frequency before 
        integration.
        
        Returns None if no sphinx profile is available yet!

        @keyword units: The unit system in which the integrated intensity is 
                        returned. Can be 'si' or 'cgs'.
                        
                        (default: 'si')
        @param units: string
        
        @return: The integrated intrinsic intensity of the line profile in 
                 SI or cgs units. (W/m2 or erg/s/cm2)
        @rtype: float
        
        """
        
        units = units.lower()
        if self.sphinx is None:
            self.readSphinx()
        if self.sphinx is None:
            return
        #-- Get the velocity grid of the line (without vlsr), convert to cm/s 
        mvel = self.sphinx.getVelocityIntrinsic()*10**5
        #-- Get the intrinsic intensity of the line in erg/s/cm2/Hz
        mint = self.sphinx.getLPIntrinsic()
        
        #-- Convert velocity grid to frequency grid, with self.frequency as 
        #   zero point (the rest frequency of the line, without vlsr) in Hz
        #   Negative velocities increase the frequency (blueshift), positive 
        #   velocities decrease the frequency (redshift).
        freqgrid = self.frequency*(1-mvel/self.c)
        
        #-- Integrate Fnu over frequency to get integrated intensity and 
        #   multiply by -1 due to a descending frequency grid rather than 
        #   ascending (causing the integrated value to be negative).
        intint_cgs = -1*trapz(x=freqgrid,y=mint)
        intint_si = intint_cgs*10**-3
        if intint_cgs < 0:
            raise IOError('Negative integrated flux found! Double check what is happening!')
        if units == 'si':
            return intint_si
        else:
            return intint_cgs
            
    
    
    def getIntConIntSphinx(self):
        
        """
        Calculate the integrated convolved intensity of the sphinx line profile
        over velocity.
        
        Returns None if no sphinx profile is available yet!

        @return: The integrated convolved intensity of the line profile in
                 erg km/s/s/cm2
        @rtype: float
        
        """
        
        if self.sphinx is None:
            self.readSphinx()
        if self.sphinx is None:
            return
        mvel = self.sphinx.getVelocity()
        mcon = self.sphinx.getLPConvolved()
        return trapz(x=mvel,y=mcon)
    
    
    
    def getIntTmbSphinx(self):
        
        """
        Calculate the integrated Tmb of the sphinx line profile over velocity.
        
        Returns None if no sphinx profile is available yet!

        @return: The integrated model Tmb profile in K km/s
        @rtype: float
        
        """
        
        if self.sphinx is None:
            self.readSphinx()
        if self.sphinx is None:
            return
        mvel = self.sphinx.getVelocity()
        mtmb = self.sphinx.getLPTmb()
        return trapz(x=mvel,y=mtmb)
    
    
    
    def getPeakTmbSphinx(self):
    
        """
        Get the peak Tmb of the sphinx line profile. 
        
        Returns None if no sphinx profile is available yet.
        
        Is equal to the mean of the 5 points around the center of the profile.
        
        @return: The peak Tmb of the sphinx line profile
        @rtype: float        
        
        """
        
        if self.sphinx is None:
            self.readSphinx()
        if self.sphinx is None:
            return
        mtmb = self.sphinx.getLPTmb()
        imid = len(mtmb)/2
        return mean(mtmb[imid-2:imid+3])
         
    
    
    def getIntTmbData(self):
        
        """
        Calculate the integrated Tmb of the data line profile over velocity.
        
        Note that only the first of data profiles is used for this, if there 
        are multiple profiles available for this transition. (i.e. multiple
        observations of the same transition with the same telescope)
        
        Makes use of the results from the fitLP method. If no extra gaussian is
        used, the integrated data profile is returned. Otherwise, the soft 
        parabola fit is integrated instead to avoid taking into account an 
        absorption feature in the profile.
        
        Returns None if no data are available. 
        
        This does not work for PACS or SPIRE data.
        
        @return: The integrated data Tmb in K km/s
        @rtype: float
        
        """
        
        if self.fittedlprof is None:
            self.fitLP()
        if self.fittedlprof is None:
            return
        if self.fittedlprof['fitgauss'] <> None:
            #-- Integrating the fitted SoftPar, rather than data
            #   due to detected absorption component in the line profile.
            return self.fittedlprof['fgintint']
        else:
            #-- Using broader integration window for data
            #   due to a Gaussian-like profile, rather than a SP.
            return self.fittedlprof['dintint']
        
        
    
    def getPeakTmbData(self):
    
        """
        Get the peak Tmb of the data line profile. 
        
        Returns None if no sphinx profile or data are available yet.
        
        Is equal to the mean of the 5 points in the data profile around the 
        center of the sphinx profile (ie in the same velocity bin as 
        getPeakTmbSphinx). Makes use of the best_vlsr, so that is ran first.
        
        This does not work for PACS or SPIRE data.

        @return: The peak Tmb of the sphinx line profile
        @rtype: float        
        
        """

        if self.fittedlprof is None:
            self.fitLP()
        if self.fittedlprof is None:
            return None
        
        info_path = os.path.join(self.path_combocode,'Data')
        #-- Do not use the best vlsr for the data peak determination. This 
        #   should be model independent.
        return LPTools.getPeakLPData(lprof=self.lpdata[0],\
                                     info_path=info_path)
    
    
    def setIntIntPacs(self,fn,dint,dint_err,st_blends=None):
        
        """
        Set the integrated intensity for this transition measured in given 
        filename. (in SI units of W/m2)
        
        A negative value is given for those lines suspected of being in a blend
        both by having at least 2 model transitions in a fitted line's breadth
        or by having a fitted_fwhm/pacs_fwhm of ~ 1.4 or more.
        
        @param fn: The data filename of PACS that contains the measured 
                   integrated intensity.
        @type fn: string
        @param dint: The value for the integrated intensity in W/m2.
        @type dint: float
        @param dint_err: The fitting uncertainty on the intensity  + absolute 
                         flux calibration uncertainty of 20%.
        @type dint_err: float
        
        @keyword st_blends: List of sample transitions involved in a line blend
                            detected by finding multiple central wavs in a PACS
                            wavelength resolution bin
        
                            (default: None)
        @type st_blends: list[Transition()]
        """
        
        self.intintpacs[fn] = float(dint)
        self.intinterrpacs[fn] = float(dint_err)
        self.intintpacs_blends[fn] = st_blends
        
    
    def getIntIntPacs(self,fn):
    
        '''
        If already set, the integrated intensity can be accessed here based on
        filename (multiple measurements can be available for a single 
        transition). 
        
        Always set and returned in W/m2!
        
        Must have been set through setIntIntPacs()! Otherwise returns None.
        
        A negative value is given for those lines suspected of being in a blend
        both by having at least 2 model transitions in a fitted line's breadth
        or by having a fitted_fwhm/pacs_fwhm of ~ 1.4 or more.
        
        @param fn: The data filename of PACS that contains the measured 
                   integrated intensity.
        @type fn: string
        
        @return: The integrated intensity measured in PACS for this filename, 
                 in SI units of W/m2, and the fitting uncertainty + absolute 
                 flux calibration uncertainty of 20%.
        @rtype: (float,float)
        
        '''
        
        if self.intintpacs.has_key(fn):
            return (self.intintpacs[fn],\
                    self.intinterrpacs[fn],\
                    self.intintpacs_blends[fn])
        else:
            return (None,None,None)
        
    
    def getLoglikelihood(self):
        
        """
        Calculate the loglikelihood of comparison between sphinx and dataset.
        
        Gives a measure for the goodness of the fit of the SHAPE of the 
        profiles.
        
        Done only for the first dataset! Makes use of the filtered sphinx 
        profile for the best vlsr, see self.getBestVlsr()
        
        Returns None if sphinx or data profile are not available. 
        
        Rescales the sphinx profile according to the difference in integrated
        Tmb between dataset and sphinx profile.

        @return: The loglikelihood
        @rtype: float
        
        """
        
        if self.best_vlsr is None:
            self.getBestVlsr()
        if self.best_vlsr is None:
            return None
        vel = self.lpdata[0].getVelocity()
        noise = self.lpdata[0].getNoise()
        window = self.fittedlprof['intwindow']
        vexp = self.fittedlprof['vexp']
        vlsr = self.getVlsr()
        if self.fittedlprof['fitgauss'] <> None:
            dsel = self.fittedlprof['fitprof'].evaluate(vel)
            dsel = dsel[abs(vel-vlsr)<=window*vexp]
        else:
            dsel = self.lpdata[0].getFlux()[abs(vel-vlsr)<=window*vexp]
        if self.getPeakTmbData() <= 5.*self.getNoise():
            #-- If the data are very noisy, use the fitted line profile to 
            #   determine the shift_factor, instead of the data themself.
            shift_factor = self.fittedlprof['fgintint']/self.getIntTmbSphinx()
        else:
            #-- Note that even if data are not noisy, the fitted lprof is still
            #   used here, in case an absorption is detected. 
            shift_factor = self.getIntTmbData()/self.getIntTmbSphinx()
        msel = self.best_mfilter[abs(vel-vlsr)<=window*vexp]
        msel = msel*shift_factor
        return bs.calcLoglikelihood(data=dsel,model=msel,noise=noise)