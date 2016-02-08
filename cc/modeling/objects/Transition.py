# -*- coding: utf-8 -*-

"""
Toolbox for Transitions, used in various applications concerning GASTRoNOoM.

Author: R. Lombaert

"""

import os 
import re
import scipy
import subprocess
from glob import glob
from scipy import pi, exp, linspace, argmin, array, diff, mean, isfinite
from scipy.interpolate import interp1d
from scipy.integrate import trapz
import types

from ivs.sigproc import funclib

import cc.path
from cc.modeling.objects import Molecule 
from cc.tools.io import Database, DataIO
from cc.tools.io import SphinxReader
from cc.tools.io import FitsReader, TxtReader
from cc.tools.numerical import Interpol
from cc.statistics import BasicStats as bs
from cc.data import LPTools


def getLineStrengths(trl,mode='dint',nans=1,n_data=0,scale=0,**kwargs):

    
    '''
    Get the line strengths from Transition() objects, either from data or from
    models defined by the mode. 
    
    Additional keywords for the functions that return LSs from the Transition()
    object can be passed along as **kwargs.
    
    Modes:
        - dint: Integrated line strength of an observed line in spectral mode
        - mint: Intrinsic integrated line strength of a modeled line
        - dtmb: Integrated line strength in main-beam temperature (Kelvin*km/s) 
                of a heterodyne instrument
        - mtmb: Modeled line strength after convolution with beam profile and 
                conversion to main-beam temperature.
        - cint: A combination of dint and mint. Data objects are assumed to be 
                listed first.
        - ctmb: A combination of dtmb and mtmb. Data objects are assumed to be 
                listed first.
    For data: Blended lines are returned as negative. If a line is in a blend,
    but not attributed a line strength, the line strength of the 'main 
    component' is returned, also as a negative value.
    
    For now no errors for any mode, except mode==dint or cint.
    
    @param trl: The Transition() objects.
    @type trl: list[Transition()]
    
    @keyword mode: The mode in which the method is called. Either 'dint', 
                   'mint', 'mtmb', 'dtmb', 'cint' or 'ctmb' values.
                   
                   (default: 'dint')
    @type mode: str
    @keyword nans: Set undefined line strengths as nans. Errors are set as a 
                   nan if it concerns mode==dint. Otherwise, they are not set.
                   
                   (default: 1)
    @type nans: bool
    @keyword n_data: The number of data Star() objects, assuming they are the 
                     first in the star_grid. Only required if mode == 'combo'.
                     
                     (default: 0)
    @type n_data: int
    @keyword scale: Scale data to antenna of 1 m**2, necessary if comparing data
                    from different telescopes
                  
                    (default: 0)
    @type scale: bool
    
    @return: The requested line strengths in W/m2, or K*km/s, as well as errors
             if applicable. If a combo mode is requested, errors are given when
             available, and listed as None if not available (Plotting2 module 
             knows how to deal with this).
    @rtype: (list[float],list[float])
    
    '''
    
    modes = {'dint': 'getIntIntUnresolved',\
             'mint': 'getIntIntIntSphinx',\
             'dtmb': 'getIntTmbData',\
             'mtmb': 'getIntTmbSphinx'}
    
    n_data = int(n_data)
    if mode[0] == 'c':
        mode = 'd%s'%mode[1:]
    elif mode[0] == 'd':
        n_data = len(trl)
    else: 
        n_data = 0
    mode = mode.lower()
    if mode not in modes.keys():
        print 'Mode unrecognized. Can be dint, mint, dtmb, mtmb, cint or ctmb.'
        return
    
    allints = []
    allerrs = []
    for it,t in enumerate(trl):
        scaling = t.telescope_size**2
        if n_data and it == n_data:
            mode = 'm%s'%mode[1:]
        if t is None:
            allints.append(nans and float('nan') or None)
            continue

        nls = getattr(t,modes[mode])(**kwargs)

        if mode == 'dint':
            if nls[0] == 'inblend':
                for tb in nls[2]:
                    bls,ebls,blends = getattr(tb,modes[mode])(**kwargs)
                    if bls != 'inblend': 
                        nls = (abs(bls)*(-1),ebls,blends)
                        break
            nint = (nans and nls[0] is None) and float('nan') or nls[0]
            nerr = (nans and nls[0] is None) and float('nan') or nls[1]
            allints.append(nint)
            allerrs.append(nerr)
            
        #-- To be done: Make sure scaling is not applied when nls[0] is None and nans is off.
        elif mode == 'dtmb':
            if scale == 1:
                allints.append(((nans and nls is None) and float('nan') or nls[0])/scaling)
                allerrs.append((nans and nls is None) and float('nan') or nls[1])
            else:
                allints.append((nans and nls is None) and float('nan') or nls[0])
                allerrs.append((nans and nls is None) and float('nan') or nls[1])
            
        else:
            if scale == 1:
                allints.append(((nans and nls is None) and float('nan') or nls)/scaling)
                allerrs.append(None)                           
            else:    
                allints.append((nans and nls is None) and float('nan') or nls)
                allerrs.append(None)            
    
    
    allints, allerrs = array(allints), array(allerrs)
    return (allints,allerrs)


def getTransFromStarGrid(sg,criterion,mode='index'):
    
    '''
    Select a transition from a list of Star() objects based on a given 
    criterion.
    
    Different modes are possible: 'index' and 'sample'. The former returns a 
    transition with given list-index in the Star()['GAS_LINES'] which should
    always return the same transition for every object. The latter returns all 
    transitions that are equal to the a given sample transition (ie following 
    the equality rules of a Transition() object).
    
    @param sg: The grid of Star() objects.
    @type sg: list[Star()]
    @param criterion: The selection criterion, based on the mode.
    @type criterion: int/Transition()
    
    @keyword mode: The selection mode. For now either 'index' or 'sample'.
    
                   (default: 'index')
    @type mode: string
    
    @return: The selected Transition() objects are returned. 
    @rtype: list[Transition()]
    
    '''
    
    if mode.lower() == 'index':
        criterion = int(criterion)
    elif mode.lower() == 'sample':
        pass
    else:
        print 'Keyword mode was not recognized. Set to "index" or "sample".'
        return
    
    if mode == 'index':
        trl = [s['GAS_LINES'][criterion] for s in sg]
    elif mode == 'sample':
        trl = [s.getTransition(criterion) for s in sg]
    return trl



def extractTransFromStars(star_grid,sort_freq=1,sort_molec=1,dtype='all'):
    
    '''
    Extract a list a of unique transitions included in a list of Star() objects
    and sort them.
    
    A selection can be made on data type, which is either a telescope or 
    instrument, or all, or all unresolved, or all resolved.
    
    The list of transitions is copied to make sure no funky references mess 
    with the original transitions. list() and set() make sure of this.
    
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
    @keyword dtype: 'all': return all lines, 
                    'resolved': return all resolved lines,
                    'unresolved': return all unresolved lines, 
                    'pacs'/'spire'/'apex'/'jcmt'/...: return all telescope or 
                    instrument specific lines selection is based on presence of
                    string in telescope name.
                    Invalid definitions return an empty list.
                    
                    (default: 'all')
    @type dtype: str
    
    @return: a list of unique transitions included in all Star() objects in
             star_grid
    @rtype: list[Transition()]
    
    '''
    
    dtype = dtype.upper()
    selection = list(sorted(set([trans 
                                 for star in star_grid
                                 for trans in star['GAS_LINES']]),\
                            key=lambda x:sort_freq \
                                and (sort_molec and x.molecule or '',\
                                     x.frequency) \
                                or  (sort_molec and x.molecule or '',\
                                     x.wavelength)))
    if dtype == 'UNRESOLVED':
        return [trans for trans in selection if trans.unresolved]
    elif dtype == 'RESOLVED':
        return [trans for trans in selection if not trans.unresolved]
    elif dtype == 'ALL':
        return selection
    else:
        return [trans for trans in selection if dtype in trans.telescope]



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
        if not ('HIFI' in telescope):
            print 'Warning! Telescope beam efficiencies for %s'%telescope + \
                  ' are added arbitrarily and thus are not the correct values.'
        old_spec = DataIO.readFile(os.path.join(cc.path.gdata,telescope+'.spec'))
        line_spec_list = [line 
                          for line in old_spec 
                          if line.find('LINE_SPEC') >-1 and line[0] != '#']
        line_spec_quant = [line.split()[0:9]
                           for line in line_spec_list]
        
        #-- Check for doubles (left over issue from when the check was done 
        #   based on the full line spec including the beamwidth, which has 
        #   changed slightly since then. The entry with the real efficiency 
        #   (not 1) is chosen over the other one for reference. The beamwidths
        #   are relatively close.
        new_lsl = []
        new_lsq = []
        for lsq,lsl in zip(line_spec_quant,line_spec_list):
            if line_spec_quant.count(lsq) == 1: 
                new_lsq.append(lsq)
                new_lsl.append(lsl)
            else:
                if lsq not in new_lsq:
                    if float(lsl.split()[10]) != 1.:
                        new_lsq.append(lsq)
                        new_lsl.append(lsl)
                        
        these_trans = [tr 
                       for tr in trans_list 
                       if tr.telescope == telescope]
        these_trans = [tr.getLineSpec()
                       for tr in these_trans
                       if tr.getLineSpec().split()[0:9] not in new_lsq]
        
        new_lsl = sorted(new_lsl + these_trans,\
                         key=lambda x: [x.split()[0],float(x.split()[9])])
        new_lsl = [line.replace(' ','\t') for line in new_lsl]
        try:
            line_spec_index = [line[0:9] 
                               for line in old_spec].index('LINE_SPEC')
            new_spec = old_spec[0:line_spec_index] + new_lsl
        except ValueError:
            new_spec = old_spec + new_lsl

        DataIO.writeFile(os.path.join(cc.path.gdata,telescope+'.spec'),\
                         new_spec+['\n######################################'])
 

def makeTransition(trans,star=None,def_molecs=None):
    
    '''
    Create a Transition instance based on a Star object and a standard CC input 
    line, with 11 or 12 entries in a list. This method assumes the Transition 
    definition has been screened by DataIO.checkEntryInfo. If the 12th entry is
    not present, n_quad is taken from the Star() object. It is set to 100 if no
    Star() object is given.
    
    If a star object is not given, the method creates Molecule() objects itself
    For now only 12C16O, 13C16O, 1H1H16O and p1H1H16O are possible.
    
    @param trans: the input line with 11 or 12 entries, the first being the 
                  molecule, followed by all 10 or 11 CC input parameters for a 
                  transition. It is possible also to give the GASTRoNOoM syntax 
                  for a transition, without n_quad in the end. In that case, 
                  n_quad is taken from Star() or equal to 100. Any parameters 
                  given after the 11th, respectively the 10th, parameter are 
                  ignored.
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

    @return: The transition object is returned with all info included
    @rtype: Transition()
    
    '''
    
    #-- The GASTRoNOoM syntax, with long molecule name (always includes a '.')
    if trans[0].find('.') != -1:
        imolec = DataIO.getInputData(keyword='MOLEC_TYPE',\
                                     filename='Molecule.dat',\
                                     make_float=0).index(trans[0])
        molec_short = DataIO.getInputData(keyword='TYPE_SHORT',\
                                          filename='Molecule.dat')[imolec]
    else:
        molec_short = trans[0]

    #-- If not Star() or default molecules given, define a few:
    if star is None and def_molecs is None: 
        def_molecs = {'12C16O':Molecule.Molecule('12C16O',61,61,240),\
                      '13C16O':Molecule.Molecule('13C16O',61,61,240),\
                      '1H1H16O':Molecule.Molecule('1H1H16O',39,90,1157),\
                      'p1H1H16O':Molecule.Molecule('p1H1H16O',32,90,1029)}
    
    #-- Put in default info if Star() is not given, otherwise take from Star()
    if star is None:
        molec = def_molecs.get(molec_short,None)
        path_gastronoom = None
        umis = 0
        n_quad = 100
    else:
        molec = star.getMolecule(molec_short)
        umis = star['USE_MASER_IN_SPHINX']
        n_quad = star['N_QUAD']
        path_gastronoom = star.path_gastronoom
    
    #-- If more than 11 entries are present, n_quad may be given. Try it!
    #   This means manual definitions of n_quad may be different for transitions
    #   in the same Star() object (and sphinx id). Typically not the case, but 
    #   can help in case a single transition fails numerically.
    if len(trans) > 11:
        try:
            n_quad = int(trans[11])
        except ValueError:
            pass
    
    if molec <> None:
        return Transition(molecule=molec,\
                          vup=int(trans[1]),jup=int(trans[2]),\
                          kaup=int(trans[3]),kcup=int(trans[4]),\
                          vlow=int(trans[5]),jlow=int(trans[6]),\
                          kalow=int(trans[7]),kclow=int(trans[8]),\
                          telescope=trans[9],offset=float(trans[10]),\
                          n_quad=n_quad,use_maser_in_sphinx=umis,\
                          path_gastronoom=path_gastronoom)
    else:
        return None    



def getModelIds(filepath):
    
    '''
    Return the modelids for a given model folder, as well as path_gastronoom.    
    
    @param filepath: The path to the model folder.
    @type filepath: str
    
    @return: the 3 ids and the path_gastronoom are returned
    @rtype: (model_id,molec_id,trans_id,path_gastronoom)
    
    '''
    
    #-- clip model_id and save it into trans_id
    filepath,trans_id = os.path.split(filepath)
    #-- clip 'models'
    filepath = os.path.split(filepath)[0]
    #-- clip path_gastronoom and save it
    path_gastronoom = os.path.split(filepath)[1]
    
    #-- Convenience path
    cc.path.gout = os.path.join(cc.path.gastronoom,path_gastronoom)
    
    cooling_log = os.path.join(cc.path.gout,'models',trans_id,'cooling_id.log')
    if os.path.isfile(cooling_log):
        model_id = DataIO.readFile(cooling_log)[0]
        #-- If an mline_id.log exists, a cooling_id.log will always exist also
        mline_log = os.path.join(cc.path.gout,'models',trans_id,'mline_id.log')
        if os.path.isfile(mline_log):
            molec_id = DataIO.readFile(mline_log)[0]
        #-- ie trans id is the same as molec id, first trans calced for id
        else: 
            molec_id = trans_id
    #-- ie mline and trans id are same as model id, first calced for id
    else: 
        model_id = trans_id
        molec_id = trans_id
    
    return (model_id,molec_id,trans_id,path_gastronoom)



def makeTransitionFromSphinx(filename,use_maser_in_sphinx=0,\
                             pull_keys_from_sphdb=1,mline_db=None):
    '''
    Make a Transition() based on the filename of a Sphinx file.
    
    For this, information is taken from the Molecule() model database and the
    sphinx file has to be associated with a molecule available there. 
    
    Information is also be pulled from the sph database, but can be turned off.
    
    @param filename: The sphinx file name, including path
    @type filename: string
    
    @keyword pull_keys_from_sphdb: Pull keywords from the sphinx database
    
                                   (default: 1)
    @type pull_keys_from_sphdb: bool
    @keyword use_maser_in_sphinx: The value that is given to the keyword 
                                  USE_MASER_IN_SPHINX for the the Transiton()
                                  object. If different values are needed, you 
                                  will have to change them manually. Only 
                                  required if no information is taken from the 
                                  sph db
                                  
                                  (default: 0)
    @type use_maser_in_sphinx: bool
    @keyword mline_db: The mline database, which can be passed in case one 
                       wants to reduce overhead. Not required though.
                       
                       (default: None)
    @type mline_db: Database()
    
    @return: The transition
    @rtype: Transition()
    
    '''
    
    if pull_keys_from_sphdb and save_in_sphdb:
        raise IOError('Cannot pull keys from database as well as save to ' + \
                       'the same database. Adapt input.')
    
    filepath,fn = os.path.split(filename)
    if not filepath:
        raise IOError('Please include the full filepath of the Sphinx ' + \
                      'filename, needed to determine the name of the ' + \
                      'database.')
    
    model_id,molec_id,trans_id,path_gastronoom = getModelIds(filepath)
    
    #-- Split the filename in bits that can be used as input for Transition()
    file_components = fn.rstrip('.dat').split('_')[2:]
    molec_name = file_components.pop(0)
    
    #-- Make the molecule
    molec = Molecule.makeMoleculeFromDb(molecule=molec_name,molec_id=molec_id,\
                                        path_gastronoom=path_gastronoom,\
                                        mline_db=mline_db)
    
    #-- If no entry available in the mline database, return None
    if molec is None:
        return None
    
    telescope = file_components.pop(-2)
    
    #-- Create a string pattern for the quantum numbers
    pattern = re.compile(r'^(\D+)(\d*.?\d*)$')
    numbers = [pattern.search(string).groups() for string in file_components]
    numbers = dict([(s.lower(),n) for s,n in numbers])

    #-- If no sphinx database is available, get info from files where possible
    #   or from input to this method
    if not pull_keys_from_sphdb:
        #-- Extract n_quad from sph2 file
        sph2file = DataIO.readFile(filename=filename,delimiter=None)
        i1 = sph2file[3].find('frequency points and ')
        i2 = sph2file[3].find(' impact parameters')
        numbers['n_quad'] = int(sph2file[3][i1+len('frequency points and '):i2])
        numbers['use_maser_in_sphinx'] = use_maser_in_sphinx
        
    #-- Make the transition object
    trans = Transition(molecule=molec,telescope=telescope,\
                       path_gastronoom=path_gastronoom,**numbers)
    
    #-- Set the model id and return
    trans.setModelId(trans_id)
    
    #-- If keys are required to be pulled from sph database, do that now
    if pull_keys_from_sphdb:
        trans_db = Database.Database(os.path.join(cc.path.gastronoom,\
                                                  path_gastronoom,\
                                                'GASTRoNOoM_sphinx_models.db'))
        trans_dict = trans_db[model_id][molec_id][trans_id][str(trans)].copy()
        [setattr(trans,k.lower(),v) for k,v in trans_dict.items()
                                    if k != 'TRANSITION']
    
    return trans

    
    
def sphinxDbRecovery(path_gastronoom,use_maser_in_sphinx=0):

    '''    
    Reconstruct a sphinx database based on existing model_ids, presence of sph2
    files and inclusion in the mline database. 
    
    Not based at first on the cooling database because it is possible different
    mline or sphinx models with different ids exist for the same cooling id.
    
    @param path_gastronoom: The path_gastronoom to the output folder
    @type path_gastronoom: string
    
    @keyword use_maser_in_sphinx: The value that is given to the keyword 
                                  USE_MASER_IN_SPHINX for the database entry. 
                                  If different values are needed, you will have
                                  to change them manually.
                                  
                                  (default: 0)
    @type use_maser_in_sphinx: bool
    
    '''

    #-- Convenience path
    cc.path.gout = os.path.join(cc.path.gastronoom,path_gastronoom)
    umis = int(use_maser_in_sphinx)
    
    #-- Make a backup of the existing sphinx db if it exists.
    sph_db_path = os.path.join(cc.path.gout,'GASTRoNOoM_sphinx_models.db')
    if os.path.isfile(sph_db_path):
        i = 0
        backup_file =  '%s_backupSphDbRetrieval_%i'%(sph_db_path,i)
        while os.path.isfile(backup_file):
            i += 1
            backup_file = '%s_backupSphDbRetrieval_%i'%(sph_db_path,i)
        subprocess.call(['mv %s %s'%(sph_db_path,backup_file)],shell=True)
    
    #-- Make a list of all sph2 files that are present in path_gastronoom
    all_sph2 = glob(os.path.join(cc.path.gout,'models','*','sph2*'))
    
    #-- Read the mline_db to decrease overhead
    mline_db = Database.Database(os.path.join(cc.path.gout,\
                                              'GASTRoNOoM_mline_models.db'))
    trans_db = Database.Database(os.path.join(cc.path.gout,\
                                              'GASTRoNOoM_sphinx_models.db'))
    #-- for all sph2 files, extract all id's, make a transition and add to db
    for i,sph2 in enumerate(all_sph2):
        if i%1000 == 0:
            print('Saving sphinx result %i, with filename %s.'%(i,sph2))
        fp,filename = os.path.split(sph2)
        model_id,molec_id,trans_id,pg = getModelIds(fp)
        trans = makeTransitionFromSphinx(filename=sph2,\
                                         use_maser_in_sphinx=umis,\
                                         pull_keys_from_sphdb=0,\
                                         mline_db=mline_db)
        #-- if transition is returned as None, the associated mline model was 
        #   not found in the mline database. ie don't add this to sphinx db 
        if trans is None:
            continue
        if not trans_db.has_key(model_id):
            trans_db[model_id] = dict()
        if not trans_db[model_id].has_key(molec_id):
            trans_db[model_id][molec_id] = dict()
        if not trans_db[model_id][molec_id].has_key(trans_id):
            trans_db[model_id][molec_id][trans_id] = dict()
        if not trans_db[model_id][molec_id][trans_id].has_key(str(trans)):
            trans_db[model_id][molec_id][trans_id][str(trans)] \
                = trans.makeDict()
        else: 
            print "There's already an entry in the sphinx db for"
            print "model_id: %s"%model_id
            print "mline_id: %s"%molec_id
            print "sphinx_id: %s"%trans_id
            print "transition: %s"%str(trans)
            print "Check what is going on!"
    
    trans_db.sync()
    
    
    
def makeTransitionsFromRadiat(molec,telescope,ls_min,ls_max,ls_unit='GHz',\
                              n_quad=100,offset=0.0,use_maser_in_sphinx=0,\
                              path_gastronoom=None,no_vib=0):
    
    '''
    Make Transition() objects from a Radiat file of a molecule, within a given
    wavelength/frequency range. 
    
    Requires a Molecule() object to work!
    
    @param molec: The molecule for which the line list is made.
    @type molec: Molecule()
    @param telescope: The telescope for which the Transition() list is made.
    @type telescope: string
    @param ls_min: The minimum allowed wavelength/frequency for the transitions
    @type ls_min: float
    @param ls_max: The maximum allowed wavelength/frequency for the transitions
    @type ls_max: float
    
    @keyword ls_unit: The unit of the wavelength/frequency range. Can be: GHz, 
                      MHz, Hz, MICRON, MM, CM, M
    
                      (default: 'GHz')
    @type ls_unit: string
    @keyword n_quad: The N_QUAD value for GASTRoNOoM sphinx calculations. 
    
                     (default: 100)
    @type n_quad: int
    @keyword offset: The offset from center position for calculations in sphinx
    
                     (default: 0.0)
    @type offset: float
    @keyword use_maser_in_spinx: using maser calc in sphinx code
                                     
                                 (default: 0)
    @type use_maser_in_spinx: bool
    @keyword path_gastronoom: model output folder in the GASTRoNOoM home
        
                              (default: None)
    @type path_gastronoom: string
    @keyword no_vib: Do not include vibrational states in the output list.
                     
                     (default: 0)
    @type no_vib: bool
    
    @return: The newly made Transition() objects for all transitions in range
    @rtype: list[Transition()]    
    
    '''
    
    radiat = molec.radiat
    wave = radiat.getFrequency(unit=ls_unit)
    low = radiat.getLowerStates()
    up = radiat.getUpperStates()
    if not molec.spec_indices:
        #- molec.ny_low is the number of levels in gs vib state
        #- molec.ny_up is the number of levels above gs vib state
        #- generally ny_up/ny_low +1 is the number of vib states
        ny_low = molec.ny_low
        n_vib = int(molec.ny_up/ny_low) + 1 
        indices = [[i+1,int(i/ny_low),i-ny_low*int(i/ny_low)]
                    for i in xrange(molec.ny_low*n_vib)]
        nl = [Transition(molecule=molec,telescope=telescope,\
                         vup=int(indices[u-1][1]),\
                         jup=int(indices[u-1][2]),\
                         vlow=int(indices[l-1][1]),\
                         jlow=int(indices[l-1][2]),\
                         offset=offset,n_quad=n_quad,\
                         use_maser_in_sphinx=use_maser_in_sphinx,\
                         path_gastronoom=path_gastronoom)
              for l,u,w in zip(low,up,wave)
              if w > ls_min and w < ls_max]
    else:
        indices = molec.radiat_indices
        quantum = ['v','j','ka','kc']
        nl = [] 
        for l,u,w in zip(low,up,wave):
            if w > ls_min and w < ls_max:
                quantum_dict = dict()
                #- some molecs only have 2 or 3 quantum numbers
                for i in xrange(1,len(indices[0])):    
                    quantum_dict[quantum[i-1]+'up'] = \
                                            int(indices[u-1][i])
                    quantum_dict[quantum[i-1]+'low'] = \
                                            int(indices[l-1][i])
                nl.append(Transition(molecule=molec,n_quad=n_quad,\
                                     telescope=telescope,offset=offset,\
                                     use_maser_in_sphinx=use_maser_in_sphinx,\
                                     path_gastronoom=path_gastronoom,\
                                     **quantum_dict))
    if no_vib:
        nl = [line for line in nl if line.vup == 0]
    return nl



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
            #-- Only add data files if there are any to begin with.
            if trans.datafiles <> None:
                ddict = dict(zip(trans.datafiles,trans.fittedlprof))
                merged[merged.index(trans)].addDatafile(ddict)
    return merged
    
    

class Transition():
    
    '''
    A class to deal with transitions in GASTRoNOoM.
    
    '''
    
    def __init__(self,molecule,telescope=None,vup=0,jup=0,kaup=0,kcup=0,\
                 nup=None,vlow=0,jlow=0,kalow=0,kclow=0,nlow=None,offset=0.0,\
                 frequency=None,exc_energy=None,int_intensity_log=None,\
                 n_quad=100,use_maser_in_sphinx=0,vibrational='',\
                 path_gastronoom=None):
        
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
                         relevant for non-point corrected PACS/SPIRE or 
                         resolved data (so not relevant when the intrinsic line
                         strengths are used from the models)
                         
                         (default: 0.0)
        @type offset: float
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
                #print 'WARNING! Water lines should not be included in the ' +\
                      #'%s.spec file. Create a file with the '%telescope + \
                      #'same name, appending -H2O to the telescope name and '+\
                      #'removing all LINE_SPEC lines.'
                telescope = '%s-H2O'%telescope
            self.telescope = telescope
            self.readTelescopeProperties()
        self.vup = int(vup)
        self.jup = int(jup)
        self.kaup = int(kaup)
        self.kcup = int(kcup)
        self.vlow = int(vlow)
        self.jlow = int(jlow)
        self.kalow = int(kalow)
        self.kclow = int(kclow)
        self.offset = float(offset)
        self.n_quad = int(n_quad)
        self.use_maser_in_sphinx = int(use_maser_in_sphinx)
        self.__model_id = None
        if nup is None or nlow is None:
            self.nup = self.kaup            
            self.nlow = self.kalow
            #- In case of SO/PO/OH, the quantum number is treated as Kaup/low but is
            #- in fact Nup/low
        self.exc_energy = exc_energy
        self.int_intensity_log = int_intensity_log
        self.vibrational = vibrational
        self.sphinx = None
        self.path_gastronoom = path_gastronoom
        self.unresolved = 'PACS' in self.telescope or 'SPIRE' in self.telescope
        self.datafiles = None
        self.lpdata = None 
        self.fittedlprof = None
        self.radiat_trans = None
        if frequency is None:
             #-- sets frequency from GASTRoNOoM input in s^-1
            self.__setIndices()  
        else:
            self.frequency = frequency
        self.c = 2.99792458e10          #in cm
        self.h = 6.62606957e-27         #in erg*s Planck constant
        self.k = 1.3806488e-16          #in erg/K Boltzmann constant
        self.wavelength = self.c/self.frequency #in cm
        #-- The vlsr from Star.dat (set by the unresolved-data objects such as
        #   Spire or Pacs)
        self.vlsr = None
        self.best_vlsr = None
        self.best_mtmb = None
        #-- PACS and SPIRE spectra are handled differently (setIntIntUnresolved)
        if self.unresolved: 
            self.unreso = dict()
            self.unreso_err = dict()
            self.unreso_blends = dict()
        else:
            self.unreso = None
            self.unreso_err = None
            self.unreso_blends = None
        #

        
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



    def readTelescopeProperties(self):
    
        """
        Read the telescope properties from Telescope.dat. 
        
        This currently includes the telescope size in m, and the default 
        absolute flux calibration uncertainty. 
        
        """
        
        all_telescopes = DataIO.getInputData(keyword='TELESCOPE',start_index=5,\
                                             filename='Telescope.dat')
        if 'PACS' in self.telescope: 
            telescope = 'PACS'
        else:
            telescope = self.telescope
        try:
            tel_index = all_telescopes.index(telescope)
        except ValueError:
            raise ValueError('Requested telescope not found in Telescope.dat.')
        
        self.telescope_size = DataIO.getInputData(keyword='SIZE',start_index=5,\
                                                  filename='Telescope.dat',\
                                                  rindex=tel_index)
        self.tel_abs_err = DataIO.getInputData(keyword='ABS_ERR',start_index=5,\
                                               filename='Telescope.dat',\
                                               rindex=tel_index)


    def getInputString(self,include_nquad=1):
         
        '''
        Return a string giving the CC input line for this transition.
        
        This includes N_QUAD by default, but can be excluded if needed. The 
        shorthand naming convention of the molecule is used, as opposed to the
        str() method.
        
        @keyword include_nquad: Include the nquad number at the end of the str
                                
                                (default: 1)
        @type include_nquad: bool
        
        @return: the input line for this transition
        @rtype: string
        
        '''
        
        if include_nquad:         
            return 'TRANSITION=%s %i %i %i %i %i %i %i %i %s %.1f %i' \
                   %(self.molecule.molecule,self.vup,self.jup,self.kaup,\
                     self.kcup,self.vlow,self.jlow,self.kalow,self.kclow,\
                     self.telescope,self.offset,self.n_quad)
        else:
            return 'TRANSITION=%s %i %i %i %i %i %i %i %i %s %.1f' \
                   %(self.molecule.molecule,self.vup,self.jup,self.kaup,\
                     self.kcup,self.vlow,self.jlow,self.kalow,self.kclow,\
                     self.telescope,self.offset)

    
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
            data = DataIO.readCols(os.path.join(cc.path.aux,\
                                                'Pacs_beamsize_v4.dat'))
            wav = data[0]/10000.
            beam = data[1]
            
            #-- Go from micron to cm. self.wavelength is in cm!
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
        #filename = os.path.join(cc.path.gdata,self.telescope+'.spec')
        #telescope_diameter = [float(line.split('=')[1][0:line.split('=')[1]\
        #                            .index('!')].strip())
        #                      for line in DataIO.readFile(filename)
        #                      if line.find('TELESCOPE_DIAM') == 0][0] * 100.
        
        #- 1.22 is diffraction limited specification, 
        #- last factor is conversion to arcseconds
        telescope_diameter = self.telescope_size * 100
        return 1.22*self.wavelength/telescope_diameter*60.*60.*360./(2.*pi)  
        


    def getEfficiency(self):
        
        '''
        Calculate telescope beam efficiency. 
        
        Telescope specific!
        
        The beam efficiency is included in the .spec files for GASTRoNOoM.
        
        This number however is not used in any calculations of the line profile 
        and is included here for reference only. 
        
        The numbers currently are correct only for HIFI.
        
        @return: The beam efficiency
        @rtype: float
        
        '''
        
        if self.unresolved:
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
        


    def getEnergyUpper(self,unit='cm-1'):
         
        '''
        Return the energy level of the upper state of this transition.
        
        You can choose the unit. Options: ['cm-1','K','erg'] Default is taken
        if the unit is not recognized.
    
        @keyword unit: The unit of the returned value. ['cm-1','K','erg']
        
                       (default: cm-1)
        @type unit: str

        @return: energy level in cm^-1
        @rtype: float
    
        '''
        
        if self.molecule.radiat is None:
            print '%s_radiat.dat not found. Cannot find energy levels.'\
                  %self.molecule.molecule
            return
        if unit.lower() not in ['cm-1','k','erg']:
            print('WARNING: unit not recognized. Returning cm-1.')
        
        energy = self.molecule.radiat.getEnergyLevels()
        
        if unit.upper() == 'K':
            return float(energy[self.up_i-1])*self.c*self.h/self.k
        elif unit.lower() == 'erg':
            return float(energy[self.up_i-1])*self.c*self.h
        else:
            return float(energy[self.up_i-1])

    
    
    def getEnergyLower(self,unit='cm-1'):

        '''
        Return the energy level of the lower state of this transition.
        
        You can choose the unit. Options: ['cm-1','K','erg']
    
        @keyword unit: The unit of the returned value. ['cm-1','K','erg']
        
                       (default: cm-1)
        @type unit: str

        @return: energy level in cm^-1
        @rtype: float
    
        '''
        
        if self.molecule.radiat is None:
            print '%s_radiat.dat not found. Cannot find energy levels.'\
                  %self.molecule.molecule
            return
        if unit.lower() not in ['cm-1','k','erg']:
            print('WARNING: unit not recognized. Returning cm-1.')
            
        energy = self.molecule.radiat.getEnergyLevels()
        
        if unit.upper() == 'K':
            return  float(energy[self.low_i-1])*self.c*self.h/self.k
        elif unit.lower() == 'erg':
            return  float(energy[self.low_i-1])*self.c*self.h
        else:
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
                      'non-unique transition indices for %s! Abort!'\
                      %self.getInputString(include_nquad=0))
        self.frequency = float(self.radiat_trans['frequency'])


    def makeAxisLabel(self,include_mol=1):
        
        '''
        Make a label suitable to be used on an axis in a plot.
        
        At the moment always returns a label typical for integrated line 
        strength.
        
        @keyword include_mol: Include the molecule name in the label.
        
                              (default: 1)
        @type include_mol: bool
        
        @return: The label
        @rtype: str
        
        '''
        
        t = self.makeLabel(inc_vib=0).replace('$','',2)
        if include_mol:
            m = self.molecule.molecule_plot.replace('$','',4)
            axislabel = 'I_{\mathrm{%s} (%s)}'%(m,t)
        else:
            axislabel = 'I_{%s}'%(t)
        return axislabel
    
        

    def makeLabel(self,inc_vib=1,return_vib=0):
        
        '''
        Return a short-hand label for this particular transition. 
        
        These labels can be used for plot line identifications for instance.
        
        If self.vibrational is not None, it always concerns a line list and is 
        included as well.
        
        @keyword inc_vib: Include the vibrational state in the label. Is always
                          True is self.vibrational is not None.
        
                          (default: 1)
        @type inc_vib: bool
        @keyword return_vib: Only return a label for the vibrational state. 
                             Does not work if vibrational is not None.
        
                             (default: 0)
        @type return_vib: bool
        
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
                if ((self.vup == 0 and self.vlow ==0) or not inc_vib)\
                        and not return_vib:
                    return r'$J=%i - %i$' %(self.jup,self.jlow)
                elif return_vib:
                    return r'$\nu=%i$' %(self.vup)
                else:
                    return r'$\nu=%i$, $J=%i-%i$' \
                           %(self.vup,self.jup,self.jlow)
            elif self.molecule.spec_indices == 1:
                ugly = r'J_{\mathrm{K}_\mathrm{a}, \mathrm{K}_\mathrm{c}}'
                if ((self.vup == 0 and self.vlow ==0) or not inc_vib) \
                        and not return_vib:
                    return r'$%s=%i_{%i,%i} - %i_{%i,%i}$'\
                            %(ugly,self.jup,self.kaup,self.kcup,\
                              self.jlow,self.kalow,self.kclow)
                elif return_vib:
                    if self.vup == 0:
                        return r'$\nu=0$'
                    else:
                        dvup = {1:2, 2:3}
                        return r'$\nu_%i=1$'%(dvup[self.vup])
                else:
                    dvup = {1:2, 2:3}
                    return r'$\nu_%i=1$, $%s=%i_{%i,%i} - %i_{%i,%i}$'\
                             %(dvup[self.vup],ugly,self.jup,self.kaup,\
                               self.kcup,self.jlow,self.kalow,self.kclow)
            elif (self.molecule.spec_indices == 2 \
                     or self.molecule.spec_indices == 3):
                ugly = r'J_{\mathrm{N}}'
                if ((self.vup == 0 and self.vlow ==0) or not inc_vib)\
                        and not return_vib:
                    return r'$%s=%i_{%i} - %i_{%i}$'\
                           %(ugly,self.jup,self.nup,self.jlow,self.nlow)
                elif return_vib:
                    return r'$\nu=%i$' %(self.vup)
                else:
                    return r'$\nu=%i$, $%s=%i_{%i} - %i_{%i}$'\
                           %(self.vup,ugly,self.jup,self.nup,\
                             self.jlow,self.nlow)

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
            gout = os.path.join(cc.path.gastronoom,self.path_gastronoom)
            filename = os.path.join(gout,'models',self.getModelId(),\
                                    self.makeSphinxFilename())
            self.sphinx = SphinxReader.SphinxReader(filename)
     
     
     
    def addDatafile(self,datadict):
        
        '''
        Add a datadict for this transition. Includes filenames as keys, and 
        fit results as values (can be None, in which case the filename is 
        excluded)
        
        The filenames are saved in self.datafiles, the fitresults in 
        self.fittedlprof.
        
        If datadict has been given a valid value, the self.lpdata list is set
        back to None, so the datafiles can all be read again. 
        
        When called for a SPIRE or PACS transition, no datafiles are added.
        
        The fit results include: 
        The gas terminal velocity, its error, the soft parabola function and 
        possibly the extra gaussian will be set. It is possible the SP is a 
        gaussian instead, if a soft parabola did not work well (the no gamma
        parameter is included in the dictionary). 
        
        The fit results are taken from the radio database which has the option 
        to (re-)run the LPTools.fitLP method. If no fit results are available
        there, no data can be used for this instance of Transition().
        
        Typically, the entry db[star_name][transition] is what is given here.
        
        @param datadict: the filenames and associated fit results
        @type datadict: dict()
        
        '''
        
        #-- If datafile is nto a dict or an empty dict(), no data available, so
        #    leave things as is 
        if not type(datadict) is types.DictType or not datadict:
            return
        
        #-- PACS and SPIRE spectra handled differently (setIntIntUnresolved)
        if self.unresolved: 
            return
        
        #-- Check if fit results are available for the file. If not, remove it. 
        #   Check if the file path is defined in all cases. If not, add radio
        #   data folder to it. Save filenames in datafiles.
        if self.datafiles is None:
            self.datafiles = []
            self.fittedlprof = []
            
        for k in sorted(datadict.keys()):
            if datadict[k] <> None:
                if not os.path.split(k)[0]: 
                    self.datafiles.append(os.path.join(cc.path.dradio,k))
                else: 
                    self.datafiles.append(k)
                self.fittedlprof.append(datadict[k])
        
        #-- If lists are empty, no valid data were found. 
        if not self.datafiles:
            self.datafiles, self.fittedlprof = None,None 
            print 'No data found for %s. If there should be data: '%str(self)+\
                  'Did you run fitLP() in the Radio db?'
        #-- Datafiles have been updated, so reset the lpdata and fittedlprof 
        #   properties. Data will be read anew when next they are requested.
        self.lpdata = None
        
        

    def readData(self):
         
        '''
        Read the datafiles associated with this transition if available.
        
        '''
        
        if self.unresolved:
            return
        if self.lpdata is None:
            if self.datafiles <> None:
                self.lpdata = []
                for idf,df in enumerate(self.datafiles):
                    if df[-5:] == '.fits':
                        lprof = FitsReader.FitsReader(filename=df)
                    else:
                        lprof = TxtReader.TxtReader(filename=df)
                    lprof.setNoise(self.getVexp(idf))
                    self.lpdata.append(lprof)            
                    
    
    
    def setData(self,trans):
        
        """
        Set data equal to the data in a different Transition object, but for 
        the same transition. Does not erase data if they have already been set 
        in this object. 
        
        A check is ran to see if the transition in both objects really is the 
        same. Sphinx model id may differ! 
        
        Avoids too much overhead when reading data from files. 
        
        The same is done for the fitresults from LPTools.fitLP() taken out of 
        the Radio database.
        
        @param trans: The other Transition() object, assumes the transition 
                      is the same, but possibly different sphinx models.
        @type trans: Transition()        
        
        """
        
        if self == trans: 
            if self.datafiles is None: 
                self.datafiles = trans.datafiles
                self.fittedlprof = trans.fittedlprof
                self.lpdata = trans.lpdata
        
        
    
    def getVlsr(self,index=0):
        
        """
        Return the vlsr read from the fits file of a resolved-data object, or 
        taken from the Star.dat file in case of unresolved data. Note that 
        resolved data may also return vlsr from Star.dat if the vlsr in the 
        data file is significantly different from the value in Star.dat.
        
        This is taken from the first lpdata object available by default, but 
        can be chosen through the index keyword..
        
        Returns 0.0 if not an unresolved line, and there are no data available.
        
        This is different from the getBestVlsr() method, which determines the 
        best matching vlsr between data and sphinx, if both are available. 
        
        @keyword index: The data list index of the requested noise value
        
                        (default: 0)
        @type index: int
        
        @return: the source velocity taken from the fits file OR Star.dat. 0
                 if data are not available. In km/s!
        @rtype: float
        
        """
        
        self.readData()
        if self.lpdata: 
            return self.lpdata[index].getVlsr()
        elif not self.unresolved and not self.lpdata:
            return 0.0
        else:
            return self.vlsr
            
            
    
    def getNoise(self,index=0):
        
        '''
        Return the noise value of the FIRST data object in this transition, if
        available. 
        
        Note that None is returned if no data are available.
        
        A different index than the default allows access to the other data 
        objects.
        
        @keyword index: The data list index of the requested noise value
        
                        (default: 0)
        @type index: int
        
        @return: The noise value
        @rtype: float
        
        '''
        
        self.readData()
        if self.lpdata:
            return self.lpdata[index].getNoise(vexp=self.getVexp(index=index))
        else:
            return None
            


    def getVexp(self,index=0):
        
        '''
        Get the gas terminal velocity as estimated from a line profile fitting
        routine for the FIRST data object. 
        A different index than the default allows access to the other data 
        objects.
        
        @keyword index: The data list index of the requested noise value
        
                        (default: 0)
        @type index: int
        
        @return: vexp
        @rtype: float
        
        '''
        
        if self.fittedlprof:
            return self.fittedlprof[index]['vexp']
        else:
            return None
  
            
    def getBestVlsr(self,index=0):
        
        """ 
        If self.best_vlsr is None, the best source velocity will be guessed by
        comparing chi^2 values between different values for the source velocity
        
        May be different from input value vlsr, the original expected 
        source velocity (from Star.dat)!
        
        By default, based on the first dataset in self.lpdata. Note that
        multiple datasets might be included. If so, a different index can be 
        given to do the scaling based on a different data object. Scaling is 
        done for only one datase. Since different datasets are for the same 
        telescope, the same line and (usually) the same conditions, the source 
        velocity is not expected to be different for the different datasets 
        anyway. 
        
        Method will attempt to call self.readData and readSphinx if either are 
        not available. If data are still not available, the initial guess is 
        returned. If data are available, but sphinx isn't, vlsr from the fits
        files is returned, and the initial guess is returned in case of txt 
        files.
       
        @keyword index: The data list index of the requested noise value
        
                        (default: 0)
        @type index: int
        
        @return: the best guess vlsr, or the initial guess if no sphinx or data
                 are available [will return vlsr included in fitsfiles if 
                 available].
        @rtype: float
        
        """
        
        
        
        #-- check if vlsr was already calculated
        if self.best_vlsr <> None:
            return self.best_vlsr

        #-- Read the data.
        #   This will set the vexp, soft parabola and gaussian profiles
        self.readData()
        
        #-- Cannot be done for unresolved lines. 
        #   Cannot be done if sphinx has not been calculated.
        #   Then, return the vlsr from Star.dat
        self.readSphinx()
        if self.unresolved or not self.lpdata or not self.sphinx:
            return self.getVlsr(index=index)
        
        #-- get all the profiles and noise values
        noise = self.getNoise(index=index)
        dvel = self.lpdata[index].getVelocity()
        dtmb = self.lpdata[index].getFlux()
        mvel = self.sphinx.getVelocity()
        mtmb = self.sphinx.getLPTmb()
        
        #-- Finding the best vlsr:
        #   1) interpolate the sphinx model, after rescaling
        #      the sphinx velocity grid to the given vlsr of the data. The flux 
        #      is assumed to be 0.0 outside the sphinx profile. 
        interpolator = interp1d(x=mvel+self.getVlsr(index=index),\
                               y=mtmb,fill_value=0.0,bounds_error=False)

        #   2) Check in the interval [vlsr-0.5vexp:vlsr+0.5*vexp] with steps 
        #      equal to the data bin size if there is a better match between
        #      model and data. This gives the 'best_vlsr'
        #-- Number of values tested is int(0.5*vexp/res+1),0.5*vexp on one side 
        #   and on the other side
        res = dvel[1]-dvel[0]
        nstep = int(0.5*self.getVexp(index=index)/res+1)
        
        #-- Use the interpolation for the nstep*2+1 shifted velocity grids
        mtmb_grid = [interpolator(dvel+i*res) for i in range(-nstep,nstep+1)]
        
        #-- Note that we shift the data velocity grid, while we should be 
        #   shifting the model velocity grid instead with several vlsr values. 
        #   Therefore perform the inverse operation to determine the actual vlsr
        vlsr_grid = [self.getVlsr(index)-i*res for i in range(-nstep,nstep+1)]

        #-- Calculate the chi squared for every shifted model
        chisquared = [bs.calcChiSquared(data=dtmb[dtmb>=-3*noise],\
                                        model=mtmbi[dtmb>=-3*noise],\
                                        noise=noise)
                      for mtmbi in mtmb_grid]
                      
        #-- Get the minimum chi squared, set the best_vlsr and set the best 
        #   shifted model profile
        imin = argmin(chisquared)
        self.chi2_best_vlsr = chisquared[imin]
        self.best_vlsr = vlsr_grid[imin]

        #-- Note that the velocity grid of best_mtmb is the data velocity
        self.best_mtmb = mtmb_grid[imin]
        
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
        @type units: string
        
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
        #-- Get the cont_subtracted intrinsic intensity in erg/s/cm2/Hz
        mint = self.sphinx.getLPIntrinsic(cont_subtract=1)
        
        #-- Convert velocity grid to frequency grid, with self.frequency as 
        #   zero point (the rest frequency of the line, without vlsr) in Hz
        #   Negative velocities increase the frequency (blueshift), positive 
        #   velocities decrease the frequency (redshift).
        freqgrid = self.frequency*(1-mvel/self.c)
        
        #-- Integrate Fnu over frequency to get integrated intensity and 
        #   multiply by -1 due to a descending frequency grid rather than 
        #   ascending (causing the integrated value to be negative).
        intint_cgs = -1*trapz(x=freqgrid[isfinite(mint)],\
                              y=mint[isfinite(mint)])
        intint_si = intint_cgs*10**-3
        if intint_cgs < 0:
            print 'WARNING! Negative integrated flux found for %s with id %s!'\
                  %(str(self),self.getModelId())
            #raise IOError('Negative integrated flux found! Double check what is happening!')
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
        mcon = self.sphinx.getLPConvolved(cont_subtract=1)
        
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
        mtmb = self.sphinx.getLPTmb(cont_subtract=1)
        
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
        mtmb = self.sphinx.getLPTmb(cont_subtract=1)
        imid = len(mtmb)/2
        return mean(mtmb[imid-2:imid+3])
         
    
    
    def getIntTmbData(self,index=0,use_fit=0):
        
        """
        Calculate the integrated Tmb of the data line profile over velocity.
        
        Note that by default only the first of data profiles is used for this, 
        if there are multiple profiles available for this transition. (i.e. 
        multiple observations of the same transition with the same telescope)
        
        A different index than the default allows access to the other data 
        objects.
        
        Makes use of the results from the LPTools.fitLP method taken from the 
        Radio database upon adding datafiles. If no extra 
        gaussian is used, the integrated data profile is returned. Otherwise, 
        the soft parabola fit is integrated instead to avoid taking into 
        account an absorption feature in the profile.
        
        The fitted line profile can be forced to be used for the integrated line
        strength.
        
        The uncertainty is also returned. Three options: 
        - The default absolute flux calibration uncertainty from Telescope.dat
        - The above + the fitting uncertainty [TBI]
        - The abs flux cal uncertainty set in Radio.py + the fitting uncertainty [TBI]
        The fitting uncertainty is currently not yet implemented, nor the option
        to add the the flux calibration uncertainty to Radio.py. [TBI]
        
        Returns None if no data are available. 
        
        This does not work for PACS or SPIRE data.
        
        @keyword index: The data list index of the requested noise value
        
                        (default: 0)
        @type index: int
        
        @return: The integrated data Tmb in K km/s and the relative uncertainty
        @rtype: (float,float)
        
        """
        
        self.readData()
        if self.fittedlprof is None:
            return
        if self.fittedlprof[index].has_key('abs_err'): 
            abs_err = self.fittedlprof[index]['abs_err']
        else:
            abs_err = self.tel_abs_err 
        if use_fit or self.fittedlprof[index]['fitabs'] <> None:
            #-- Integrating the fitted SoftPar, rather than data
            #   due to detected absorption component in the line profile.
            return (self.fittedlprof[index]['fgintint'],abs_err)
        else:
            #-- Using broader integration window for data
            #   due to a Gaussian-like profile, rather than a SP.
            return (self.fittedlprof[index]['dintint'],abs_err)
        
        
    
    def getPeakTmbData(self,index=0):
    
        """
        Get the peak Tmb of the data line profile. 
        
        Returns None if no sphinx profile or data are available yet.
        
        Is equal to the mean of the 5 points in the data profile around the 
        center of the sphinx profile (ie in the same velocity bin as 
        getPeakTmbSphinx). Makes use of the best_vlsr, so that is ran first.
        
        Note that by default only the first of data profiles is used for this, 
        if there are multiple profiles available for this transition. (i.e. 
        multiple observations of the same transition with the same telescope)
        
        A different index than the default allows access to the other data 
        objects.
        
        This does not work for PACS or SPIRE data.
        
        @keyword index: The data list index of the requested noise value
        
                        (default: 0)
        @type index: int
        
        @return: The peak Tmb of the sphinx line profile
        @rtype: float        
        
        """
        
        self.readData()
        if self.fittedlprof is None:
            return None
        
        #-- Do not use the best vlsr for the data peak determination. This 
        #   should be model independent.
        return LPTools.getPeakLPData(lprof=self.lpdata[index])
    
    
    
    def setIntIntUnresolved(self,fn,dint,dint_err,vlsr,st_blends=None,):
        
        """
        Set the integrated intensity for this transition measured in given 
        filename. (in SI units of W/m2)
        
        A negative value is given for those lines suspected of being in a blend
        both by having at least 2 model transitions in a fitted line's breadth
        or by having a fitted_fwhm/intrinsic_fwhm of ~ 1.2 or more.
        
        The vlsr is also passed through this function as that is only available
        from the data objects (in this case the Spire or Pacs class). For 
        unresolved lines, vlsr is read from Star.dat.
        
        @param fn: The data filename that contains the measured 
                   integrated intensity.
        @type fn: string
        @param dint: The value for the integrated intensity in W/m2. If the 
                     line is part of a blend that has already been added, this
                     may also say 'inblend'. All transitions involved in the 
                     blend are then given by st_blends. 
        @type dint: float or string
        @param dint_err: The fitting uncertainty on the intensity  + absolute 
                         flux calibration uncertainty of 20%. Relative value!
        @type dint_err: float
        @param vlsr: The source velocity with respect to the local standard of 
                     rest in cm/s.
        @type vlsr: float
        
        @keyword st_blends: List of sample transitions involved in a line blend
                            detected by finding multiple central wavs in a
                            wavelength resolution bin
        
                            (default: None)
        @type st_blends: list[Transition()]
        
        """
        
        if dint == 'inblend':
            self.unreso[fn] = dint
        else:
            self.unreso[fn] = float(dint)
        self.unreso_err[fn] = dint_err <> None and float(dint_err) or None
        self.unreso_blends[fn] = st_blends
        #-- Set the vlsr, but convert to km/s for plotting purposes.
        self.vlsr = vlsr * 10**(-5)
        
        
    
    def getIntIntUnresolved(self,fn=''):
    
        '''
        If already set, the integrated intensity can be accessed here based on
        filename (multiple measurements can be available for a single 
        transition). 
        
        Always set and returned in W/m2!
        
        Must have been set through setIntIntUnresolved()! Otherwise returns 
        None.
        
        A negative value is given for those lines suspected of being in a blend
        both by having at least 2 model transitions in a fitted line's breadth
        or by having a fitted_fwhm/intrinsic_fwhm of ~ 1.2 or more.
        
        @keyword fn: The data filename that contains the measured integrated
                     intensity. Can be set to '' or None if simply the
                     first entry in the keys() list is to be used. Mostly only 
                     one line strength is associated with the object anyway.
                     
                     (default: '')
        @type fn: string
        
        @return: The integrated intensity measured in unresolved data for this 
                 filename, in SI units of W/m2, and the fitting uncertainty + 
                 absolute flux calibration uncertainty (relative value!), and 
                 the blends if applicable.
        @rtype: (float,float,list)
        
        '''
        
        if self.unreso.has_key(fn):
            return (self.unreso[fn],\
                    self.unreso_err[fn],\
                    self.unreso_blends[fn])
        elif not fn and self.unreso.keys():
            k = self.unreso.keys()[0]
            return (self.unreso[k],\
                    self.unreso_err[k],\
                    self.unreso_blends[k])
        else:
            return (None,None,None)
        
    
    
    def getLoglikelihood(self,use_bestvlsr=1,index=0,partial=0,vcut=0.0):
        
        """
        Calculate the loglikelihood of comparison between sphinx and dataset.
        
        Gives a measure for the goodness of the fit of the SHAPE of the 
        profiles.
        
        Note that by default only the first of data profiles is used for this, 
        if there are multiple profiles available for this transition. (i.e. 
        multiple observations of the same transition with the same telescope)
        
        A different index than the default allows access to the other data 
        objects.
        
        Done for the dataset with given index! Makes use of the interpolated 
        sphinx profile for the best vlsr, see self.getBestVlsr() if use_bestvlsr 
        is True. If this keyword is False, interpolates the sphinx model for the
        vlsr from Star.dat or the fits file.
        
        Returns None if sphinx or data profile are not available. 
        
        Rescales the sphinx profile according to the difference in integrated
        Tmb between dataset and sphinx profile.

        @keyword use_bestvlsr: Use the fitted best-guess for the v_lsr when 
                               determining the velocity grid for the model. If 
                               not, the vlsr from the Star.dat file or the fits
                               file is used. 
                               
                               (default: 1)
        @type use_bestvlsr: bool
        @keyword partial: Use a partial profile rather than the entire profile. 
                          Set to 1 for larger than the cutoff value (vcut), set 
                          to -1 for smaller than the cutoff value.
                          
                          (default: 0)
        @type partial: int
        @keyword vcut: The cut off value in km/s used by partial.
                       
                       (default: 0.0)
        @type vcut: float        
        @keyword index: The data list index of the requested noise value
        
                        (default: 0)
        @type index: int
        
        @return: The loglikelihood
        @rtype: float
        
        """
        
        if use_bestvlsr and self.best_vlsr is None:
            self.getBestVlsr(index=index)
        if use_bestvlsr and self.best_vlsr is None:
            print 'Using standard v_lsr from Star.dat or fits file for LLL.'
            use_bestvlsr = 0
            
        vel = self.lpdata[index].getVelocity()
        noise = self.getNoise(index=index)
        window = self.fittedlprof[index]['intwindow']
        vexp = self.getVexp(index=index)
        vlsr = self.getVlsr(index=index)
        
        #-- Select the line profile within the relevant window, and cut off 
        #   part in case partial lll is requested
        if partial > 0:
            selection = (abs(vel-vlsr)<=window*vexp)*(vel>vcut)
        elif partial < 0:
            selection = (abs(vel-vlsr)<=window*vexp)*(vel<vcut)
        else:
            selection = abs(vel-vlsr)<=window*vexp

        if self.fittedlprof[index]['fitabs'] <> None:
            pars = array(self.fittedlprof[index]['fitprof'][1])
            functype = self.fittedlprof[index]['fitprof'][0]
            dsel = funclib.evaluate(functype,vel[selection],pars)
        else:
            dsel = self.lpdata[index].getFlux()[selection]
        if self.getPeakTmbData(index=index) <= 5.*self.getNoise(index=index):
            #-- If the data are very noisy, use the fitted line profile to 
            #   determine the shift_factor, instead of the data themself.
            num = self.fittedlprof[index]['fgintint']
            shift_factor = num/self.getIntTmbSphinx()
        else:
            #-- Note that even if data are not noisy, the fitted lprof is still
            #   used here, in case an absorption is detected. 
            num,abs_err = self.getIntTmbData(index=index)
            shift_factor = num/self.getIntTmbSphinx()
        
        if use_bestvlsr:
            msel = self.best_mtmb[selection]
            msel = msel*shift_factor
        else:
            mvel = self.sphinx.getVelocity()
            mtmb = self.sphinx.getLPTmb()
            interpolator = interp1d(x=mvel+self.getVlsr(index=index),y=mtmb,\
                                    fill_value=0.0,bounds_error=False)
            mtmb_interp = interpolator(vel[selection])
            msel = mtmb_interp*shift_factor
        
        return bs.calcLoglikelihood(data=dsel,model=msel,noise=noise)