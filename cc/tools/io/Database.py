# -*- coding: utf-8 -*-

"""
A simple interface to work with a database saved on the hard disk. 

Author: Robin Lombaert

"""

import os
import cPickle
import time
import subprocess
import portalocker
from glob import glob

import cc.path
from cc.tools.io import DataIO



def updateAllDbs(func,db_name,*args,**kwargs):

    ''' 
    Update all databases located in any of your modeling folders in 
    cc.path.gastronoom. 
    
    Choose the method to run as well as the database to update, and add the args
    and kwargs for that method. 
    
    This method synchronizes the databases, so run with care. 
    
    @param func: The method to use for updating a database
    @type func: function
    @param db_name: Name of the database to be updated
    @type db_name: str
    
    '''

    #-- Select all model folders in cc.path.gastronoom
    gpaths = sorted(glob(os.path.join(cc.path.gastronoom,'*','GASTRoNOoM*.db')))
    gpaths = list(set([os.path.split(gp)[0] for gp in gpaths]))
    
    #-- Run db update method for all databases
    for gp in gpaths: 
    
        #-- Set filename
        fn = os.path.join(gp,db_name)
        print "******************************"
        print "Now converting database at:"
        print fn
        
        #-- Call requested function
        db = func(db_fn=fn,*args,**kwargs)
        db.sync()



def convertDbKeys(path_input=''):

    '''
    Convert your local databases such that some keywords are shifted around 
    according to their relevance for the different subcodes of GASTRoNOoM.
    
    Moves USE_MASER_IN_SPHINX from sphinx to mline. Sets it to the default value
    of 1, since it was never changed for mline. 
    
    Moves USE_NO_MASER_OPTION, N_FREQ, START_APPROX, USE_FRACTION_LEVEL_CORR, 
    FRACTION_LEVEL_CORR, NUMBER_LEVEL_MAX_CORR from cooling to mline.
    
    Moves WRITE_INTENSITIES, TAU_MAX, TAU_MIN, CHECK_TAU_STEP from cooling to 
    sphinx. 
    
    Adds the new keywords N_IMPACT_EXTRA[_RIN/ROUT] to cooling.
    
    Adds the new keywords FRACTION_TAU_STEP, MIN_TAU_STEP to sphinx. 
    
    Adds the new keyword FEHLER to mline.
    
    Converts any 'double' notation in str format to the floats, e.g. 
    TAU_MIN='-6d0' becomes TAU_MIN=-6.
    
    Converts all cooling, mline, sphinx and pacs databases found in the 
    GASTRoNOoM home folder. 
    
    Can Update ComboCode input files for the maser keywords.
    
    @keyword path_input: The location of your ComboCode inputfiles. Use empty
                         string or None if you do not want to update your input-
                         files automatically. Is directly inserted into glob so
                         takes wildcards, eg /Users/robinl/ComboCode/input/*.dat
    
                         (default: '') 
    @type path_input: str
    
    '''
    
    #-- Keys to be moved around or added as new. Given as (key,default) pairs
    keys_cool_ml = [('USE_NO_MASER_OPTION',0,int),('N_FREQ',30,int),\
                    ('START_APPROX',0,int),('USE_FRACTION_LEVEL_CORR',1,int),\
                    ('FRACTION_LEVEL_CORR',0.8,float),\
                    ('NUMBER_LEVEL_MAX_CORR',1e-12,float)]
    keys_cool_sph = [('WRITE_INTENSITIES',0,int),('TAU_MAX',12,float),\
                     ('TAU_MIN',-6,float),('CHECK_TAU_STEP',0.01,float)]

    #-- USE_MASER_IN_SPHINX belongs here as it was never properly used in sphinx
    #   If this was already corrected earlier, adding the key here will do 
    #   nothing.
    #-- Not adding USE_STARFILE. It's now in Input_Keywords_Mline for checking 
    #   when it is 1 or 2. When it's 0, it's not even added (since STARFILE 
    #   is not relevant either.)
    new_keys_ml = [('FEHLER',1e-4,float),('USE_MASER_IN_SPHINX',1,int)]    
    new_keys_sph = [('FRACTION_TAU_STEP',1e-2,float),\
                    ('MIN_TAU_STEP',1e-4,float)]
    new_keys_cool = [('N_IMPACT_EXTRA',0,int),\
                     ('N_IMPACT_EXTRA_RIN',100.,float),\
                     ('N_IMPACT_EXTRA_ROUT',150.,float)]
    convert_double = ['TAU_MAX','TAU_MIN','NUMBER_LEVEL_MAX_CORR',\
                      'STEP_RS_RIN','STEP_RIN_ROUT','FRACTION_TAU_STEP',\
                      'MIN_TAU_STEP']

    gpaths = sorted(glob(os.path.join(cc.path.gastronoom,'*','GASTRoNOoM*.db')))
    gpaths = list(set([os.path.split(gp)[0] for gp in gpaths]))
    print "New databases will be located at filename_old.db_dbConversion"
    #-- Add to Mline and rm from cooling databases
    for gp in gpaths: 
        cfn = os.path.join(gp,'GASTRoNOoM_cooling_models.db')
        mfn = os.path.join(gp,'GASTRoNOoM_mline_models.db')
        sfn = os.path.join(gp,'GASTRoNOoM_sphinx_models.db')
        pfns = sorted(glob(os.path.join(gp,'stars','*','GASTRoNOoM*.db')))

        print "******************************"
        print "Now converting databases from:"
        print "\n".join([cfn,mfn,sfn]+pfns)
        
        #-- Make copies of the dbs, which will be changed. Avoids erroneous
        #   changes and data loss.
        for fn in [cfn,mfn,sfn]+pfns: 
            os.system('cp %s %s'%(fn,fn+'_dbConversion'))
        cdb = Database(cfn+'_dbConversion')
        mdb = Database(mfn+'_dbConversion')
        sdb = Database(sfn+'_dbConversion')
        
        #-- Add new cooling keys
        for key,defval,valtype in new_keys_cool:
            #-- If key was already moved previously, nothing will change
            #   Method doesn't add anything if already present. 
            cdb = addKeyCooling(key=key,val=defval,db=cdb)
        
        #-- Convert the "double" notations in the cooling db:
        for cmid in cdb.keys():
            for k,v in cdb[cmid].items():
                if k in convert_double and isinstance(v,str):
                    cdb[cmid][k] = v.replace('d','e')
                    cdb.addChangedKey(cmid)
                    
        #-- Move keywords from cooling to mline. Include the new keywords in 
        #   case they were already used by someone. If not, they won't be in the
        #   cooling db either, and they'll just be added to sph db.
        for key,defval,valtype in keys_cool_ml + new_keys_ml:
            for cmid in cdb.keys():
                if not key in cdb[cmid].keys():
                    #-- If key was already moved previously, nothing will change
                    #   Method doesn't add anything if already present. 
                    mdb = addKeyMline(key=key,val=defval,db=mdb,id=cmid)
                    continue
                val = cdb[cmid].pop(key)
                cdb.addChangedKey(cmid)
                if key in convert_double:
                    val = float(val.replace('d','e'))
                mdb = addKeyMline(key=key,val=valtype(val),db=mdb,id=cmid)
        
        #-- Move keywords from cooling to sphinx. Include the new keywords in 
        #   case they were already used by someone. If not, they won't be in the
        #   cooling db either, and they'll just be added to sph db.
        for key,defval,valtype in keys_cool_sph + new_keys_sph:
            for cmid in cdb.keys():
                if not key in cdb[cmid].keys():
                    #-- If key was already moved previously, nothing will change
                    #   Method doesn't add anything if already present. 
                    sdb = addKeySphinx(key=key,val=defval,db=sdb,id=cmid)
                    #-- Then add the key to all PACS dbs!
                    for pp in pfns:
                        pdb = addKeyPacs(key=key,val=defval,\
                                         db_fn=pp+'_dbConversion',id=cmid)
                        pdb.sync()
                    continue
                val = cdb[cmid].pop(key)
                cdb.addChangedKey(cmid)
                if key in convert_double:
                    val = float(val.replace('d','e'))
                sdb = addKeySphinx(key=key,val=valtype(val),db=sdb,id=cmid)                            
                #-- Then add the key to all PACS dbs!
                for pp in pfns:
                    pdb = addKeyPacs(key=key,val=valtype(val),\
                                     db_fn=pp+'_dbConversion',id=cmid)
                    pdb.sync()
                    
        #-- Remove USE_MASER_IN_SPHINX from sphinx databases
        sdb = rmKeySphinx(key='USE_MASER_IN_SPHINX',db=sdb)
        #-- Remove from pacs databases
        for pp in pfns:
            pdb = rmKeyPacs(key='USE_MASER_IN_SPHINX',db_fn=pp+'_dbConversion')
            pdb.sync()

        cdb.sync()
        mdb.sync()
        sdb.sync()
    
    #-- Obsolete:
    # -- Add old default value of USE_FRACTION_LEVEL_CORR
#     if path_input:
#         print "******************************"
#         print "Now converting adding old defaults to inputfiles at:"
#         print path_input
#         ifiles = glob(path_input)
#         comment = '# set to 1 if one wants to put a limit on the ' + \
#                   'level-population correction (BES3).\n'
#         for ff in ifiles:
#             lines = DataIO.readFile(ff,None,replace_spaces=0)
#             ldict = DataIO.readDict(ff)
#             -- First find a good place to add the keyword. Just after 
#               CHECK_TAU_STEP is good as it is present in all inputfiles.
#             for i,l in enumerate(lines):
#                 if l.find('CHECK_TAU_STEP') != -1:
#                     break
#             -- Then add the old default value.
#             if not ldict.has_key('USE_FRACTION_LEVEL_CORR'):
#                 lines[i:i] = ['USE_FRACTION_LEVEL_CORR=1              '+comment]
#             DataIO.writeFile(ff,lines,mode='w',delimiter='')
    
    #-- Add the proper default maser keys to inputfiles.
    if not path_input: return
    ifiles = glob(path_input)
    comment = '# set to 1 if one wants to omit masers occuring when solving '+\
              'the radiative transfer equation\n'
    for ff in ifiles:
        lines = DataIO.readFile(ff,None,replace_spaces=0)
        ldict = DataIO.readDict(ff)
        for i,l in enumerate(lines):
            if l.find('USE_MASER_IN_SPHINX') != -1:
                k = 'USE_MASER_IN_SPHINX=1                 '+l[l.find('#'):]
                lines[i] = k
                break
        if not ldict.has_key('USE_NO_MASER_OPTION'):
            lines[i:i] = ['USE_NO_MASER_OPTION=1                 '+comment]
        DataIO.writeFile(ff,lines,mode='w',delimiter='')



def updateDustMCMaxDatabase(filename):

    '''
    Update dust filenames in MCMax database with the new OPAC_PATH system. 
    
    @param filename: The file and path to the MCMax database. 
    @type filename: str
    
    '''
    
    i = 0
    new_filename = '%s_new'%(filename)
    
    db_old = Database(filename)
    db_new = Database(new_filename)
    
    path = os.path.join(cc.path.usr,'Dust_updatefile.dat')
    dustfiles = DataIO.readCols(path)
    pfn_old = list(dustfiles[0])
    pfn_new = list(dustfiles[1])
    
    for k,v in db_old.items():
        dd = v['dust_species']
        dd_new = dict()
        for pfn,cont in dd.items():
            try:
                new_key = pfn_new[pfn_old.index(pfn)]
                dd_new[new_key] = cont
            except ValueError:
                dd_new[pfn] = cont
        v['dust_species'] = dd_new
        db_new[k] = v
    db_new.sync()
        


def convertMCMaxDatabase(path_mcmax):
    
    '''
    Convert MCMax database to the dict format.
    
    This change was made to speed up the use of the database and makes use of 
    the Database() class.
    
    @param path_mcmax: the name of the MCMac subfolder.
    @type path_mcmax: string

    '''
    
    print '** Converting MCMax database to dictionary format...'
    #-- Convenience path
    cc.path.mout = os.path.join(cc.path.mcmax,path_mcmax)
    db_path = os.path.join(cc.path.mout,'MCMax_models.db')
    i = 0
    backup_file = '%s_backup%i'%(db_path,i)
    while os.path.isfile(backup_file):
        i += 1
        backup_file = '%s_backup%i'%(db_path,i)
    subprocess.call(['mv %s %s'%(db_path,backup_file)],\
                    shell=True,stdout=subprocess.PIPE)
    mcmax_db = open(backup_file,'r')
    old_db = []
    try:
        while True:
            model = cPickle.load(mcmax_db)
            old_db.append(model)
    except EOFError:
        print '** End of old database reached.'
    finally:
        mcmax_db.close()
    db = Database(db_path)
    print '** Inserting entries into the new database.'
    for commands,model_id in old_db:
        photon_count = DataIO.convertFloat(commands.pop(0),convert_int=1)
        commands = [c.split('=') for c in commands]
        commanddict = dict([(k,DataIO.convertFloat(v,convert_int=1)) 
                            for k,v in commands 
                            if k[0:4] not in ['opac','part','abun',\
                                              'Tdes','minr','maxr']])
        #- Correcting erroneous input entries in the old version
        if commanddict.has_key('densfile'):
            commanddict['densfile'] = "'%s'"%commanddict['densfile']
        if not commanddict.has_key('FLD'):
            commanddict['FLD'] = '.false.'
        if not commanddict.has_key('randomwalk'):
            commanddict['randomwalk'] = '.true.'
        
        commanddict['photon_count'] = photon_count
        speciespars = [(k,DataIO.convertFloat(v,convert_int=1)) 
                       for k,v in commands 
                       if k[0:4] in ['opac','part','abun',\
                                     'Tdes','minr','maxr']]       
        commanddict['dust_species'] = dict()
        i = 0 
        while speciespars:
            i += 1
            this_species = dict()
            for k,v in speciespars:
                if int(k[-2:]) == i:
                    if k[:4] == 'part' or k[:4] == 'opac':
                        speciesfile = v.strip("'")
                    else:
                        this_species[k.strip('%.2i'%i)] = v
            commanddict['dust_species'][os.path.split(speciesfile)[1]] \
                = this_species
            speciespars = [(k,v) for k,v in speciespars 
                           if int(k[-2:]) != i]
        db[model_id] = commanddict
    db.sync()   
    print '** Done!'
    
    
def cleanDatabase(db_path):
    
    '''
    Remove any db entries with a dictionary that includes the IN_PROGRESS key.
    
    Works for cooling, mline and sphinx databases.
    
    @param db_path: full path to the database.
    @type db_path: string
    
    '''
    
    code = os.path.split(db_path)[1].split('_')[-2]
    if code not in ['cooling','mline','sphinx','MCMax']:
        raise IOError('Database path is not related to a GASTRoNOoM or MCMax '+\
                      'database.')
    db = Database(db_path)
    #-- Locking the database while this is done to avoid issues.
    dbfile = db._open('r')
    print '****************************************************************'
    print '** Checking {} database for in-progress models now...'.format(code)
    for cool_id,vcool in db.items():
        #-- For cooling IN PROGRESS entry is found in vcool
        if code in ['cooling','MCMax']:
            if vcool.has_key('IN_PROGRESS'):
                del db[cool_id]
                print 'Removed in-progress model with id {}.'.format(cool_id)
            continue
        for ml_id,vml in vcool.items():
            for key,val in vml.items():
                #-- For mline IN PROGRESS entry is found in the molecule dict.
                if code == 'mline':
                    if val.has_key('IN_PROGRESS'):
                        del db[cool_id][ml_id][key]
                        db.addChangedKey(cool_id)
                        print 'Removed in-progress molecule {} '.format(key)+\
                              'id {}.'.format(ml_id)
                    continue

                for trans,vsph in val.items():
                    #-- For sphinx IN PROGRESS entry is found in the trans dict.
                    #   No need to check code, it's the last possibility.
                    if vsph.has_key('IN_PROGRESS'):
                        del db[cool_id][ml_id][key][trans]
                        db.addChangedKey(cool_id)
                        print 'Removed in-progress transition '+ \
                              '{} with id {}.'.format(trans,key)
                #-- Remove trans_ids that don't contain transitions. Run .keys()
                #   again because maybe in-progress deletion removed all entries
                #   from id
                if code == 'sphinx' and not db[cool_id][ml_id][key].keys():
                    del db[cool_id][ml_id][key]
                    db.addChangedKey(cool_id)
                    print 'Removed empty trans id {}.'.format(key)                    

            #-- Remove molec_ids that don't contain molecules. Run .keys() again
            #   because maybe in-progress deletion removed all entries from id
            if code == 'mline' and not db[cool_id][ml_id].keys():
                del db[cool_id][ml_id]
                db.addChangedKey(cool_id)
                print 'Removed empty trans id {}.'.format(ml_id)   
                                    
    print '** Unlocking and synchronizing the database...'
    dbfile.close()
    db.sync()
    print '** Done!'
    print '****************************************************************'



def coolingDbRetrieval(path_gastronoom,r_outer=None):
    
    '''    
    Reconstruct a cooling database based on the mline database and the
    GASTRoNOoM inputfiles.
    
    Only works if the water MOLECULE convenience keywords, the MOLECULE R_OUTER
    and/or the MOLECULE ENHANCE_ABUNDANCE_FACTOR keywords were not adapted!
    
    @param path_gastronoom: The path_gastronoom to the output folder
    @type path_gastronoom: string
    
    @keyword r_outer: The outer radius used for the cooling model, regardless
                      of the outer_r_mode parameter.
                      
                      (default: None)
    @type r_outer: float
    
    '''
    
    #-- Convenience path
    cc.path.gout = os.path.join(cc.path.gastronoom,path_gastronoom)
                                
    coolkeys_path = os.path.join(cc.path.aux,'Input_Keywords_Cooling.dat')
    coolkeys = DataIO.readCols(coolkeys_path,make_float=0,make_array=0)[0]
    extra_keys = ['ENHANCE_ABUNDANCE_FACTOR','MOLECULE_TABLE','ISOTOPE_TABLE',\
                  'ABUNDANCE_FILENAME','NUMBER_INPUT_ABUNDANCE_VALUES',\
                  'KEYWORD_TABLE']
    coolkeys = [k for k in coolkeys if k not in extra_keys]
    cool_db_path = os.path.join(cc.path.gout,'GASTRoNOoM_cooling_models.db')
    ml_db_path = os.path.join(cc.path.gout,'GASTRoNOoM_mline_models.db')
    subprocess.call(['mv %s %s_backupCoolDbRetrieval'\
                     %(cool_db_path,cool_db_path)],shell=True)
    cool_db = Database(db_path=cool_db_path)
    ml_db = Database(db_path=ml_db_path)
    for ml_id in ml_db.keys():
        file_path = os.path.join(cc.path.gout,'models',\
                                 'gastronoom_%s.inp'%ml_id)
        input_dict = DataIO.readDict(file_path)
        input_dict = dict([(k,v) for k,v in input_dict.items() 
                                 if k in coolkeys])
        cool_db[ml_id] = input_dict
        if not r_outer is None:
            cool_db[ml_id]['R_OUTER'] = r_outer
    cool_db.sync()
    
    

def addKeyCooling(key,val,db_fn='',db=None):
        
    '''
    Add a (key,value) pair to every entry in the cooling database. 
    
    Not added to a particular entry if already present.
    
    @param key: The name of the keyword to be added.
    @type key: str
    @param val: The default value for the keyword.
    @type val: any
    
    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)
    for k in db.keys():
        if not key in db[k].keys():
            db[k][key] = val
            db.addChangedKey(k)
    return db



def rmKeyCooling(key,val,db_fn='',db=None):
        
    '''
    Remove a key from every entry in the cooling database. 
    
    @param key: The name of the keyword to be removed.
    @type key: str
    
    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)
    for k in db.keys():
        if key in db[k].keys():
            del db[k][key]
            db.addChangedKey(k)
    return db
    


def addKeyMline(key,val,db_fn='',db=None,id=''):
        
    '''
    Add a (key,value) pair to every entry in the mline database.
    
    Not added to a particular entry if already present.
        
    @param key: The keyword to be added
    @type key: str
    @param val: The default value of the keyword
    @type val: any
    
    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    @keyword id: If the (key,val) pair is only to be added to one cooling id, 
                 give that id here. If not given, the pair is added to all ids
                 If id not in db, nothing is done.
                 
                 (default: '')
    @type id: str
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)
    
    if id and not db.has_key(id):
        return db
    
    cids = db.keys() if not id else [id]
    for k in cids:
        for l in db[k].keys():
            for mol in db[k][l].keys():
                if not key in db[k][l][mol].keys():
                    db[k][l][mol][key] = val
                    db.addChangedKey(k)
    return db
    
    

def rmKeyMline(key,db_fn='',db=None,id=''):
        
    '''
    Remove a key from every entry in the mline database.    
    
    @param key: The keyword to be removed
    @type key: str
    
    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    @keyword id: If the key is only to be removed from one cooling id, give that
                 id here. If not given, the key is removed from all ids
                 If id not in db, nothing is done.
                 
                 (default: '')
    @type id: str
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)

    if id and not db.has_key(id):
        return db
    
    cids = db.keys() if not id else [id]
    for k in cids:
        for l in db[k].keys():
            for mol in db[k][l].keys():
                if key in db[k][l][mol].keys():
                    del db[k][l][mol][key]
                    db.addChangedKey(k)
    return db
    


def rmKeySphinx(key,db_fn='',db=None,id=''):

    '''
    Remove a key from the sphinx database entries. 
    
    @param key: They keyword to be removed
    @type key: str
    
    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    @keyword id: If the key is only to be removed from one cooling id, give that
                 id here. If not given, the key is removed from all ids
                 If id not in db, nothing is done.
                 
                 (default: '')
    @type id: str
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)    

    if id and not db.has_key(id):
        return db
    
    cids = db.keys() if not id else [id]
    for k in cids:
        for l in db[k].keys():
            for o in db[k][l].keys():
                for trans in db[k][l][o].keys():
                    if key in db[k][l][o][trans].keys():
                        del db[k][l][o][trans][key]
                        db.addChangedKey(k)
    return db
    
    

def addKeySphinx(key,val,db_fn='',db=None,id=''):

    '''
    Add a (key,value) pair to every entry of the sphinx database entries. 
    
    @param key: They keyword to be added
    @type key: str
    @param val: The default value of the keyword
    @type val: any

    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    @keyword id: If the (key,val) pair is only to be added to one cooling id, 
                 give that id here. If not given, the pair is added to all ids
                 If id not in db, nothing is done.
                 
                 (default: '')
    @type id: str
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)    

    if id and not db.has_key(id):
        return db
    
    cids = db.keys() if not id else [id]
    for k in cids:
        for l in db[k].keys():
            for o in db[k][l].keys():
                for trans in db[k][l][o].keys():
                    if not key in db[k][l][o][trans].keys():
                        db[k][l][o][trans][key] = val
                        db.addChangedKey(k)
    return db
       

    
def rmKeyPacs(key,db_fn='',db=None,id=''):

    '''
    Remove a key from the PACS database entries. 
    
    @param key: They keyword to be removed
    @type key: str
    
    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    @keyword id: If the key is only to be removed from one cooling id, give that
                 id here. Note this is the COOLING ID, not the pacs id. 
                 Iteration over all pacs ids is always done. However, if id is 
                 given it is cross checked with the cooling id in the entry and 
                 only then added. If not given, the key is removed from all 
                 pacs ids as well as cooling ids.
                 If cooling id not in db, nothing is done.
                 
                 (default: '')
    @type id: str
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)

    for k in db.keys(): #- pacs id: dict
        #-- If cooling id does not match, do nothing and move on.
        if id and db[k]['cooling_id'] != id: continue 
        for l in db[k]['trans_list']: #- list of tuples
            if key in l[2].keys():
                del l[2][key]
                db.addChangedKey(k)
    return db
    
    
    
def addKeyPacs(key,val,db_fn='',db=None,id=''):

    '''
    Add a (key,value) pair to every entry of the PACS database entries. 
    
    @param key: They keyword to be added
    @type key: str
    @param val: The default value of the keyword
    @type val: any
    
    @keyword db_fn: The filename and path of the database. Only required if db 
                    is not given.
                    
                    (default: '')
    @type db_fn: string
    @keyword db: The database. Is updated and returned. If not given, a filename
                 is required.
                 
                 (default: None)
    @type db: Database()
    @keyword id: If the (key,val) pair is only to be added to one pacs id, 
                 give that id here. If not given, the pair is added to all ids.
                 If id not in db, nothing is done.
                 
                 (default: '')
    @type id: str
    
    @return: The new database, not yet synchronized.
    @rtype: Database()
    
    '''
    
    if db is None and not db_fn:
        return
        
    if db is None:
        db = Database(db_fn)

    for k in db.keys(): #- pacs id: dict
        #-- If cooling id does not match, do nothing and move on.
        if id and db[k]['cooling_id'] != id: continue 
        for l in db[k]['trans_list']: #- list of tuples
            if not key in l[2].keys():
                l[2][key] = val
                db.addChangedKey(k)
    return db    


    
def replaceSubstring(db,oldss,newss):

    '''
    Replace a substring of values in database entries, such as the home 
    directory of filenames. This is applied to all strings in the database vals!
    
    This applies to the value of a dictionary (key,value) pair on the second 
    level of the database (the first level being (model_id, dictionary). Does
    not work for embedded dictionaries in a dictionary, such as the dust_species
    embedded dictionary in an MCMax database!
    
    The database is not synchronized in this method.
    
    @param db: The database (for now only MCMax or cooling) 
    @type db: Database()
    @param oldss: The old home name (eg /home/mariev/ or /home/)
    @type oldss: str
    @param newss: The new home name (efg /Users/robinl/ or /Users/)
    @type newss: str
    
    '''
    
    for m,dd in db.items(): 
        for k,v in dd.items():
            if isinstance(v,str) and oldss in v:
                dd[k] = v.replace(oldss,newss)
                db.addChangedKey(m)
                
    

class Database(dict):
    
    '''
    A database class.
    
    The class creates and manages a dictionary saved to the hard disk. 
    
    It functions as a python dictionary with the extra option of synchronizing
    the database instance with the dictionary saved on the hard disk. 
    
    No changes will be made to the hard disk copy, unless Database.sync() is 
    called.
    
    Note that changes made on a deeper level than the (key,value) pairs of the 
    Database (for instance in the case where value is a dict() type itself) 
    will not be automatically taken into account when calling the sync() 
    method. The key for which the value has been changed on a deeper level has 
    to be added to the Database.__changed list by calling addChangedKey(key)
    manually.
    
    Running the Database.sync() method will not read the database from the hard
    disk if no changes were made or if changes were made on a deeper level 
    only. In order to get the most recent version of the Database, without 
    having made any changes, use the .read() method. Note that if changes were 
    made on a deeper level, they will be lost.
    
    Example:
    
    >>> import os
    >>> import Database
    >>> filename = 'mytest.db'
    >>> db = Database.Database(filename)
    No database present at mytest.db. Creating a new one.
    >>> db['test'] = 1
    >>> db['test2'] = 'robin'
    >>> db.sync()
    >>> db2 = Database.Database(filename)
    >>> print db2['test']
    1
    >>> print db2['test2'] 
    robin
    >>> db2['test'] = 2
    >>> db2.sync()
    >>> db.sync()
    >>> print db['test']
    1
    >>> db.read()
    >>> print db['test']
    2
    >>> del db2['test2']
    >>> db2.sync()
    >>> print db['test2']
    robin
    >>> db.read()
    >>> print db['test2']
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    KeyError: 'test2'
    >>> test_dict = dict()
    >>> db['test'] = test_dict
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test']
    {}
    >>> db['test']['test'] = 1
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test']
    {}
    >>> db.addChangedKey('test')
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test']
    {'test': 1}
    >>> db.setdefault('test','defkey')
    {'test': 1}
    >>> db.setdefault('test3','defval')
    'defval'
    >>> db.sync()
    >>> db2.read()
    >>> print db2['test3']
    defval
    >>> os.system('rm %s'%filename)
    0
    '''
    
    
    def __init__(self,db_path):
        
        '''
        Initializing a Database class.
        
        Upon initialization, the class will read the dictionary saved at the 
        db_path given as a dictionary.
        
        Note that cPickle is used to write and read these dictionaries.
        
        If no database exists at db_path, a new dictionary will be created.
        
        @param db_path: The path to the database on the hard disk.
        @type db_path: string
  
        '''
        
        super(Database, self).__init__()
        self.path = db_path
        self.folder = os.path.split(self.path)[0]
        self.read()
        self.__changed = []
        self.__deleted = []
      
      
      
    def __delitem__(self,key):
        
        '''
        Delete a key from the database.
        
        This deletion is also done in the hard disk version of the database 
        when the sync() method is called. 
        
        This method can be called by using syntax:
        del db[key]
        
        @param key: a dict key that will be deleted from the Database in memory
        @type key: a type valid for a dict key
        
        '''
        
        self.__deleted.append(key)
        return super(Database,self).__delitem__(key)
        
    
    
    
    def __setitem__(self,key,value):
        
        '''
        Set a dict key with value. 
        
        This change is only added to the database saved on the hard disk when
        the sync() method is called. 
        
        The key is added to the Database.__changed list.
        
        This method can be called by using syntax:
        db[key] = value
        
        @param key: a dict key that will be added to the Database in memory
        @type key: a type valid for a dict key
        @param key: value of the key to be added
        @type value: any
        
        '''
        
        self.__changed.append(key)
        return super(Database,self).__setitem__(key,value)
        
        
        
    def setdefault(self,key,*args):
        
        '''
        Return key's value, if present. Otherwise add key with value default
        and return. 
        
        Database.__changed is updated with the key if it is not present yet.
        
        @param key: the key to be returned and/or added.
        @type key: any valid dict() key
        @param args: A default value added to the dict() if the key is not 
        present. If not specified, default defaults to None.
        @type args: any type
        @return: key's value or default
                
        '''
        
        if not self.has_key(key):
            self.__changed.append(key)
        return super(Database,self).setdefault(key,*args)
        
        
            
    def pop(self,key,*args):
        
        '''
        If database has key, remove it from the database and return it, else 
        return default. 
    
        If both default is not given and key is not in the database, a KeyError
        is raised. 
        
        If deletion is successful, this change is only added to the database 
        saved on the hard disk when the sync() method is called. 
        
        The key is added to the Database.__deleted list, if present originally.
                
        @param key: a dict key that will be removed from the Database in memory
        @type key: a type valid for a dict key
        @param args: value of the key to be returned if key not in Database
        @type args: any
        @return: value for key, or default
        
        '''
        
        if self.has_key(key):
            self.__deleted.append(key)
        return super(Database,self).pop(key,*args)
        
        
    def popitem(self):
        
        '''
        Remove and return an arbitrary (key, value) pair from the database.
        
        A KeyError is raised if the database has an empty dictionary.
        
        If removal is successful, this change is only added to the database 
        saved on the hard disk when the sync() method is called. 
        
        The removed key is added to the Database.__deleted list.
                
        @return: (key, value) pair from Database
        
        '''
        
        (key,value) = super(Database,self).popitem()
        self.__deleted.append(key)
        return (key,value)
            
            
        
    def update(self,*args,**kwargs):
        
        '''
        Update the database with new entries, as with a dictionary. 
        
        This update is not synched to the hard disk! Instead Database.__changed
        includes the changed keys so that the next sync will save these changes
        to the hard disk.
        
        @param args: A dictionary type object to update the Database.
        @type args: dict()
        @keyword kwargs: Any extra keywords are added as keys with their values.
        @type kwargs: any type that is allowed as a dict key type.
        
        '''
        
        self.__changed.extend(kwargs.keys())
        self.__changed.extend(args[0].keys())
        return super(Database,self).update(*args,**kwargs)
               
               
        
    def read(self):
        
        '''
        Read the database from the hard disk.
        
        Whenever called, the database in memory is updated with the version 
        saved on the hard disk.
        
        Any changes made outside the session of this Database() instance will
        be applied to the database in memory! 
        
        Any changes made to existing keys in current memory before calling 
        read() will be undone! Use sync() instead of read if you want to keep
        current changes inside the session. 
        
        If no database is present at the path given to Database() upon 
        initialisation, a new Database is made by saving an empty dict() at the
        requested location.        
        
        Reading and saving of the database is done by cPickle-ing the dict(). 
        
        '''
        
        try:
            while True:
                #dbfile = open(self.path,'r')
                #portalocker.lock(dbfile, portalocker.LOCK_EX)
                dbfile = self._open('r')
                try:
                    try:
                        db = cPickle.load(dbfile)
                        dbfile.close()
                        break
                    except ValueError:
                        print 'Loading database failed: ValueError ~ ' + \
                              'insecure string pickle. Waiting 5 seconds ' + \
                              'and trying again.' 
                        dbfile.close()
                        time.sleep(5)
                except EOFError:
                    print 'Loading database failed: EOFError. Waiting 5 ' + \
                          'seconds and trying again.'
                    dbfile.close()
                    time.sleep(5)
            self.clear()
            super(Database,self).update(db)
        except IOError:
            print 'No database present at %s. Creating a new one.'%self.path
            self.__save()
                
                
                
    def sync(self):
        
        ''' 
        Update the database on the harddisk and in the memory.
         
        The database is read anew, ie updated with the hard disk version to 
        account for any changes made by a different program. Next, the changes
        made to the database in memory are applied, before saving the database
        to the hard disk again.
        
        Any items deleted from the database in memory will also be deleted from
        the version saved on the hard disk!
        
        The keys that are changed explicitly are all listed in self.__changed,
        to which entries can be added manually using the addChangedKey method, 
        or automatically by calling .update(), .__setitem__() or .setdefault().
        
        '''
        
        if self.__changed or self.__deleted:
            current_db = dict([(k,v) 
                               for k,v in self.items() 
                               if k in set(self.__changed)])
            while True:    
                self.read()
                self.__deleted = list(set(self.__deleted))
                for key in self.__deleted:
                    try:
                        super(Database,self).__delitem__(key)
                    except KeyError:
                        pass
                super(Database,self).update(current_db)
                backup_file = self.__save()
                try:
                    #-- Read the object, if TypeError, catch and repeat (which  
                    #   can happen if db written into by two instances of 
                    #   Database at the same time)
                    testread = Database(self.path)
                    #-- If the read object is not the same as the one in memory, 
                    #   repeat writing as well. 
                    if testread != self:
                        raise TypeError
                    #-- Remove backup if all is fine. If not, it won't be 
                    #   removed: tracer for issues if they occur.
                    if backup_file and os.path.isfile(backup_file):
                        subprocess.call(['rm %s'%(backup_file)],shell=True)
                    break
                except TypeError: 
                    #-- Just wait a few seconds to allow other instances to 
                    #   finish writing
                    time.sleep(2)
            self.__deleted = []
            self.__changed = []
        
        #-- Nothing changed in this instance of the db. Just read the db saved
        #   to hard disk to update this instance to the real-time version. 
        else:
            self.read()
    
    
    def __save(self):
        
        '''
        Save a database. 
        
        Only called by Database() internally. Use sync() to save the Database
        to the hard disk.
        
        Reading and saving of the database is done by cPickle-ing the dict(). 
        
        @return: the filename of the backup database is returned
        @rtype: string
        
        '''
        
        backup_file = ''
        if os.path.isfile(self.path):
            i = 0
            backup_file =  '%s_backup%i'%(self.path,i)
            while os.path.isfile(backup_file):
                i += 1
                backup_file = '%s_backup%i'%(self.path,i)
            subprocess.call(['mv %s %s'%(self.path,backup_file)],\
                            shell=True)
        #-- Write the file, dump the object
        #dbfile = open(self.path,'w')
        #portalocker.lock(dbfile, portalocker.LOCK_EX)
        dbfile = self._open('w')
        cPickle.dump(self,dbfile)
        dbfile.close()
        return backup_file
            
    
    
    def _open(self,mode):
    
        '''
        Open the database on the disk for writing, reading or appending access.
        
        A lock is added to the database, which remains in place until the file 
        object is closed again. 
        
        @return: The opened file 
        @rtype: file()
        
        '''
        
        dbfile = open(self.path,mode)
        portalocker.lock(dbfile, portalocker.LOCK_EX)
        return dbfile
        
    
    
    def addChangedKey(self,key):
        
        '''
        Add a key to the list of changed keys in the database.
        
        This is useful if a change was made to an entry on a deeper level, 
        meaning that the __set__() method of Database() is not called directly.
        
        If the key is not added to this list manually, it will not make it into
        the database on the hard disk when calling the sync() method.
        
        @param key: the key you want to include in the next sync() call.
        @type key: string
        
        '''
        
        if key not in self.__changed: self.__changed.append(key)
    
    
    
    def getDeletedKeys(self):
        
        '''
        Return a list of all keys that have been deleted from the database in
        memory.
        
        @return: list of keys
        @rtype: list
        '''
        
        return self.__deleted
    
    
    
    def getChangedKeys(self):
        
        '''
        Return a list of all keys that have been changed in the database in
        memory.
        
        @return: list of keys
        @rtype: list
        '''
        
        return self.__changed



if __name__ == "__main__":
    import doctest
    doctest.testmod()        


"""
def getDbStructure(db_path,id_index=None,code=None):
    
    '''
    Return id_index and code based on db_path, or id_index, or code.
    
    @param db_path: the path + filename of the database
    @type db_path: string
    
    @keyword id_index: index of model id in the database entries, specify only 
                       if you know what you're doing! default None is used if 
                       the code keyword is used or the code is taken from the 
                       database filename. 
                       
                       (default: None)
    @type id_index: int
    @keyword code: name of the (sub-) code for which the deletion is done 
                   (pacs, mcmax, cooling, mline, sphinx), default None can be 
                   used when id_index is given, or if the database filename 
                   includes the codename you want to use deletion for
                   
                   (default: None)
    @type code: string
    
    @return: The id_index and code are returned. The code can be None if 
             id_index was already defined beforehand, in which case the precise
             code doesn't matter, however the method will try to determine the 
             code based on the filename
    @rtype: (int,string)
    
    '''
    
    if id_index <> None and code <> None:
        raise IOError('Either specify the code, or the id_index or none of ' +\
                      'both.')
    code_indices = dict([('cooling',1),('mcmax',1),('sphinx',2),('mline',2),\
                         ('pacs',2)])
    if id_index is None:
        if code <> None:
            id_index = code_indices[code.lower()]
        else:
            for k,v in code_indices.items():
                if k.lower() in os.path.split(db_path)[1].lower():
                    if id_index <> None:
                        this_str = 'There is an ambiguity in the filename of'+\
                                   ' the database. At least two of the codes'+\
                                   ' are in the name.'
                        raise ValueError(this_str)
                    id_index = v
                    code = k
    if id_index is None:
        this_str2 = 'Cannot figure out which code the database is used for. '+\
                    'Please specify the "code" keyword.'
        raise IOError(this_str2)
    if code is None:
        for k in code_indices.keys():
                if k.lower() in os.path.split(db_path)[1].lower():
                    if code <> None:
                        print 'WARNING! There is an ambiguity in the ' + \
                              'filename of the database. At least two of the'+\
                              ' codes are in the name.'
                    else:
                        code = k
    return id_index,code                             
                             


def deleteModel(model_id,db_path,id_index=None,code=None):
    
    '''
    Delete a model_id from a database. A back-up is created!
    
    This method is only used by MCMax.py, which uses a database in the old 
    format.
    
    @param model_id: the model_id
    @type model_id: string
    @param db_path: the path + filename of the database
    @type db_path: string
    
    @keyword id_index: index of model id in the database entries, specify only 
                       if you know what you're doing! default None is used if 
                       the code keyword is used or if the code is taken from 
                       the database filename.
                       
                       (default: None)
    @type id_index: int
    @keyword code: name of the (sub-) code for which the deletion is done 
                   (pacs, mcmax, cooling, mline, sphinx), default None can be 
                   used when id_index is given, or if the database filename 
                   includes the codename you want to use deletion for
                   
                   (default: None)
    @type code: string
    
    '''
    
    id_index,code = getDbStructure(db_path=db_path,id_index=id_index,code=code)
    subprocess.call(['mv ' + db_path + ' ' + db_path+'old'],\
                    shell=True,stdout=subprocess.PIPE)
    gastronoom_db_old = open(db_path+'old','r')
    gastronoom_db = open(db_path,'w')
    print "Making the following change(s) to database at %s: \n"%db_path + \
          "(Note that if nothing is printed, nothing happened and the model"+\
          "id wasn't found.)"
    i = 0
    while True:
        try:
            model = cPickle.load(gastronoom_db_old)
            if model_id == model[id_index]:
                print 'Deleting model id %s from database at %s.'\
                      %(model_id,db_path)
                i += 1
            else:
                cPickle.dump(model,gastronoom_db)  
        except EOFError:
            print 'Done! %i models were deleted from the database.'%i
            break
    gastronoom_db_old.close()
    gastronoom_db.close()


def convertPacsDbToDictDb(db_path):
    
    '''
    ***OBSOLETE***
    
    Convert an old PACS database to the new dictionary based format. 
    
    Keeps the old db in existence!
    
    This change was made to speed up the use of the database and makes use of 
    the Database() class.
    
    @param db_path: Path to the PACS database
    @type db_path: string
    
    '''
    subprocess.call([' '.join(['mv',db_path,db_path+'_old'])],shell=True)
    pacs_db_old = open(db_path+'_old','r')      
    pacs_db_new = open(db_path,'w')      
    old_db = []
    while True:
        try:
            model = cPickle.load(pacs_db_old)
            old_db.append(model)
        except EOFError:
            break  
    new_db = dict()
    for entry in old_db:
        trans_list = entry[0]
        filename = entry[1]
        pacs_id = entry[2]
        model_id = entry[3]
        if not new_db.has_key(pacs_id): 
            new_db[pacs_id] = dict([('filenames',[])])
        new_db[pacs_id]['trans_list'] = trans_list
        new_db[pacs_id]['filenames'].append(filename)
        new_db[pacs_id]['cooling_id'] = model_id
    cPickle.dump(new_db,pacs_db_new)
    pacs_db_old.close()
    pacs_db_new.close()
    

    
def convertGastronoomDatabases(path_gastronoom):
    
    '''
    ***OBSOLETE***
    
    Convert all GASTRoNOoM databases (cooling, mline, sphinx) to the dict
    format.
    
    This change was made to speed up the use of the database and makes use of 
    the Database() class.
       
    @param path_gastronoom: the name of the GASTRoNOoM subfolder.
    @type path_gastronoom: string
    
    '''
    
    cooling_path = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                path_gastronoom,'GASTRoNOoM_cooling_models.db')
    convertCoolingDb(path=cooling_path)
    mline_path = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                              path_gastronoom,'GASTRoNOoM_mline_models.db')
    convertMlineDb(path=mline_path)
    sphinx_paths = glob(os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                    path_gastronoom,'models','GASTRoNOoM*.db'))
    convertSphinxDb(paths=sphinx_paths,path_gastronoom=path_gastronoom)
    print '** Done!'
    print '*******************************'


    
def convertCoolingDb(path):
    
    '''
    ***OBSOLETE***
    
    Convert the cooling db from list to dictionary format.
    
    This change was made to speed up the use of the database and makes use of 
    the Database() class.
    
    @param path: the full path to the db
    @type path: string
    
    '''
    
    print '** Converting Cooling database to dictionary format...'    
    subprocess.call([' '.join(['mv',path,path+'_old'])],shell=True)
    cool_db_old = open(path+'_old','r')      
    old_db = []
    while True:
        try:
            model = cPickle.load(cool_db_old)
            old_db.append(model)
        except EOFError:
            break  
    cool_db_old.close()
    new_db = dict()
    for entry in old_db:
        command_list = entry[0]
        model_id = entry[1]
        new_db[model_id] = command_list
    saveDatabase(db_path=path,db=new_db)
    

    
def convertMlineDb(path):
    
    '''
    ***OBSOLETE***
    
    Convert the mline db from list to dictionary format.
        
    This change was made to speed up the use of the database and makes use of 
    the Database() class.
    
    @param path: the full path to the db
    @type path: string
    
    '''
    
    print '** Converting Mline database to dictionary format...'
    subprocess.call([' '.join(['mv',path,path+'_old'])],shell=True)
    ml_db_old = open(path+'_old','r')      
    old_db = []
    while True:
        try:
            model = cPickle.load(ml_db_old)
            old_db.append(model)
        except EOFError:
            break  
    ml_db_old.close()
    new_db = dict()
    for entry in old_db:
        molecule = entry[0]
        molec_dict = entry[1]
        molec_id = entry[2]
        model_id = entry[3]
        if not new_db.has_key(model_id): 
            new_db[model_id] = dict()
        if not new_db[model_id].has_key(molec_id): 
            new_db[model_id][molec_id] = dict()
        new_db[model_id][molec_id][molecule] = molec_dict
    saveDatabase(db_path=path,db=new_db)
    
    
    
def convertSphinxDb(paths,path_gastronoom):
    
    '''
    ***OBSOLETE***
    
    Convert the sphinx db from list to dictionary format.
    
    This change was made to speed up the use of the database and makes use of 
    the Database() class.
    
    @param paths: list of all full paths to the dbs
    @type paths: list[string]
    @param path_gastronoom: the name of the GASTRoNOoM subfolder.
    @type path_gastronoom: string
    
    '''
    
    #[subprocess.call([' '.join(['mv',path,path+'_old'])],shell=True)
    #        for path in paths]if model[0] == self.command_list:
                if self.replace_db_entry and 0:
                    mcmax_db.close()
                    Database.deleteModel(model_id=str(model[1]),db_path=os.path.join(os.path.expanduser('~'),'MCMax',self.path,'MCMax_models.db'))
                    raise EOFError
                else:
                    print '** MCMax model has been calculated before with ID ' + str(model[1]) + '.'
                    self.model_id = model[1]
                break
    print '** Converting Sphinx database to dictionary format...'
    old_db = []
    for path in paths:
        sph_db_old = open(path,'r')
        while True:
            try:
                model = cPickle.load(sph_db_old)
                old_db.append(model)
            except EOFError:
                break 
        sph_db_old.close()
    new_db = dict()
    for entry in old_db:
        transition = entry[0]
        trans_dict = entry[1]
        trans_id = entry[2]
        molec_id = entry[3]
        model_id = entry[4]
        if not new_db.has_key(model_id): 
            new_db[model_id] = dict()
        if not new_db[model_id].has_key(molec_id): 
            new_db[model_id][molec_id] = dict()
        if not new_db[model_id][molec_id].has_key(trans_id): 
            new_db[model_id][molec_id][trans_id] = dict()
        new_db[model_id][molec_id][trans_id][transition] = trans_dict
    new_path = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                            path_gastronoom,'GASTRoNOoM_sphinx_models.db')
    saveDatabase(db_path=new_path,db=new_db)



def updateCoolingOpas(path_gastronoom):
    
    '''
    ***OBSOLETE***
    
    Update cooling database at path_gastronoom with the new opacity keys.
    
    This method was used when the conversion from a single input opacity file in
    GASTRoNOoM to an inputparameter filename for these opacities was made.
    
    @param path_gastronoom: the name of the GASTRoNOoM subfolder.
    @type path_gastronoom: string
    
    '''
    
    print '** Converting cooling database opacity entries.'
    cooling_path = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                path_gastronoom,'GASTRoNOoM_cooling_models.db')
    subprocess.call([' '.join(['mv',cooling_path,cooling_path+\
                               '_old_opakeys'])],shell=True)
    old_db = getDatabase(db_path=cooling_path+'_old_opakeys')
    new_db = dict()
    new_db.update(old_db)
    for model_id in new_db.keys():    
        if not old_db[model_id]['#OPA_FILE'] \
                and not old_db[model_id]['#TEMDUST_FILE']:
            new_db[model_id]['TEMDUST_FILENAME'] = '"temdust.kappa"'
            new_db[model_id]['USE_NEW_DUST_KAPPA_FILES'] = 0
        elif old_db[model_id]['#OPA_FILE'] \
                and old_db[model_id]['#TEMDUST_FILE']:
            new_db[model_id]['TEMDUST_FILENAME'] = '"%s"'\
                    %os.path.split(old_db[model_id]['#TEMDUST_FILE'])[1]
            new_db[model_id]['USE_NEW_DUST_KAPPA_FILES'] = 1        
        elif not old_db[model_id]['#OPA_FILE'] \
                and old_db[model_id]['#TEMDUST_FILE'] \
                    == '/home/robinl/GASTRoNOoM/src/data/qpr_files/'+\
                       'temdust_silicates.kappa':
            new_db[model_id]['TEMDUST_FILENAME'] = '"temdust.kappa"'
            new_db[model_id]['USE_NEW_DUST_KAPPA_FILES'] = 0        
        else:
            subprocess.call([' '.join(['mv',cooling_path+'_old_opakeys',\
                    cooling_path])],shell=True)        
            raise ValueError('Something fishy is going on... Discrepant OPA' +\
                             ' and TEMDUST files for %s!'%model_id)
        del new_db[model_id]['#OPA_FILE']
        del new_db[model_id]['#TEMDUST_FILE']
    saveDatabase(db_path=cooling_path,db=new_db)          
    gast_data = os.path.join(os.path.expanduser('~'),'GASTRoNOoM','src','data')
    qpr_ori = os.path.join(gast_data,'qpr_files','qpr_silicates_jus1992.dat') 
    if os.path.isfile(qpr_ori):
        print '** Copying original opacity files... Remember to copy qpr.dat'+\
              ' and temdust.kappa to your VIC account as well!'
        qpr = os.path.join(gast_data,'qpr.dat')
        temdust = os.path.join(gast_data,'temdust.kappa')
        temd_ori = os.path.join(gast_data,'qpr_files',\
                                'temdust_silicates.kappa')
        subprocess.call([' '.join(['cp',qpr_ori,qpr])],shell=True)
        subprocess.call([' '.join(['cp',temd_ori,temdust])],shell=True)
    print '** Done!'
    print '*******************************'
         

##############################
## FIND FIRST COOLING MODEL ##
############################## 
   
def getStartPosition(gastronoom_db,model_id):
   '''Get the start position in the db based on the cooling model_id, to speed up the cross-checking procedure.
   
   Input:   gastronoom_db = File object, the opened database
      model_id=string, the cooling model_id that you are searching for.
   Output:  start_position, the db will HAVE TO BE reset to this position, the position is -1 if the id hasn't been found
   '''
   try:
      while True:
         last_position = gastronoom_db.tell()
         model = cPickle.load(gastronoom_db)
         if str(model[-1]) == model_id:             #Find the first position of the cooling model_id occurrence, and remember it for repositioning the file object once a check has been done.
            start_position = last_position
            break
   except EOFError:
      start_position = -1
   return start_position

##############################
#CONVERT OLD SPHINX DATABASE##
##############################

def convertOldSphinxDbToNew(path):
   '''Convert an old Sphinx database to the new cooling model_id specific database format. Keeps the old db in existence!'''
   gastronoom_db_old = open(os.path.join(os.path.expanduser('~'),'GASTRoNOoM',path,\
                                 'GASTRoNOoM_sphinx_models.db'),'r')     
   while True:   
      try:
         model = cPickle.load(gastronoom_db_old)
         this_cooling_id = model[-1]
         gastronoom_db = open(os.path.join(os.path.expanduser('~'),'GASTRoNOoM',path,'models',\
                                          'GASTRoNOoM_sphinx_%s.db'%this_cooling_id),'a')
         cPickle.dump(model,gastronoom_db)  
         gastronoom_db.close()
      except EOFError:
         print 'Database conversion finished.'
         gastronoom_db_old.close()
         break

##############################
#REMOVE OBS INFO FROM  MCMAX DB
##############################     

def removeObsInfoMCMaxDb(path):
   '''Convert an old MCMax database that includes observation information (ie ray tracing )to a new database without it.
   Keeps the old db in existence!
   
   Input:   path=string, path where the database is saved (ie path_mcmax). Db is assumed to be called MCMax_models.db!
   '''
   subprocess.call([' '.join(['mv',os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models.db'),\
                              os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models_old_obsinfo.db')])],shell=True)
   mcmax_db_old = open(os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models_old_obsinfo.db'),'r')     
   mcmax_db_new = open(os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models.db'),'a')     
   while True:   
      try:
         model = cPickle.load(mcmax_db_old)
         model[0] = [m for m in model[0] if m != os.path.join(os.path.expanduser('~'),'MCMax','src','Spec.out')]
         cPickle.dump(model,mcmax_db_new)  
      except EOFError:
         print 'Database conversion finished.'
         mcmax_db_old.close()
         mcmax_db_new.close()
         break

##############################
#CONVERT OLD MCMAX DATABASE ##
##############################     
      
def convertOldMCMaxDbToNew(path):
   '''Convert an old MCMax database to the new database based on actual inputfiles instead of command line options.
   Keeps the old db in existence!
   
   Input:   path=string, path where the database is saved. Db is assumed to be called MCMax_models.db!'''
   subprocess.call([' '.join(['mv',os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models.db'),\
                              os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models_old.db')])],shell=True)
   mcmax_db_old = open(os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models_old.db'),'r')     
   mcmax_db_new = open(os.path.join(os.path.expanduser('~'),'MCMax',path,'MCMax_models.db'),'a')     
   while True:   
      try:
         model = cPickle.load(mcmax_db_old)
         model[0] = [m.replace('-s ','') for m in model[0]]
         abun_lines = [(i,m.split()) for i,m in enumerate(model[0]) if m[:4] == 'abun']
         n = len(abun_lines)
         for i,m in abun_lines:
            model[0][i] = m[0]
            model[0][i+n:i+n] = [m[1]]
         cPickle.dump(model,mcmax_db_new)  
      except EOFError:
         print 'Database conversion finished.'
         mcmax_db_old.close()
         mcmax_db_new.close()
         break
         

def updateDatabase(db_path,db):
    
    ''' 
    Update the database saved on the harddisk. The database is read again, 
    updated and then saved to ensure no data recently added to the database are
    lost.
        
    Input:    db_path=string, the path to the database.
        db=dict, the new version of the database.
    Output:  current_db=dict(), the current version of the db is returned.
    
    '''
    
    current_db = getDatabase(db_path=db_path)
    current_db.update(db)
    saveDatabase(db_path=db_path,db=current_db)
    return current_db

##############################
## SAVE GASTRONOOM DATABASE ##
##############################

def saveDatabase(db_path,db):
    '''Save a database.
    
    Input:    db_path=string, the path to the database.
        db=dict, the new version of the database.
    '''
    dbfile = open(db_path,'w')
    cPickle.dump(db,dbfile)
    dbfile.close()

##############################
### GET GASTRONOOM DATABASE###
##############################

def getDatabase(db_path):
    
    '''
    Return a database.
    
    Input:    db_path=string, the path to the database.
    '''
    try:
        dbfile = open(db_path,'r')
        while True:
            try:
                db = cPickle.load(dbfile)
                break
            except ValueError:
                print 'Loading database failed: ValueError ~ insecure string pickle. Waiting 10 seconds and trying again.' 
                time.sleep(10)
        dbfile.close()
        return db
    except IOError:
        return dict()
        



def browseDatabase(db_path,model_id,code=None,id_index=None,conditional=None,cond_index=None):
    '''Browse a database and return the entry from it based on the id.
    
    Input:    db_path=str, the filepath to the database to browse
        model_id=str, the model_id (molec specific in case of sphinx and mline, general if cooling)
        OPTIONAL id_index=int, default=None, index of model id in the database entries, specify only if you know what you're doing! 
            default None is used if the code keyword is used or the code is taken from the database filename. 
        OPTIONAL code=string, default=None, name of the (sub-) code for which the db is browsed (pacs, mcmax, cooling, mline, sphinx),
            default None can be used when id_index is given, or if the database filename includes the codename
        OPTIONAL conditional=anything, default=None, if an extra condition is imposed for one part of the db entry, cond_id has to be
            specified.
        OPTIONAL cond_index=int, default=None, index of the conditional term in the db entry, must be specified if conditional is not None
    Output:  dict giving the input parameters
    '''
    if (conditional is None and cond_index <> None) or (conditional <> None and cond_index is None):
        raise IOError('If a conditional term is added, its database index must be given as well; and vice versa. Aborting.')
    id_index,code = getDbStructure(db_path=db_path,code=code,id_index=id_index)
    star_db = open(db_path,'r')
    try:
        while True:
            entry = cPickle.load(star_db)
            if model_id == entry[id_index]:
                if conditional <> None:
                    if entry[cond_index] == conditional:
                        star_db.close()
                        break
                    else:
                        pass
                else:
                    star_db.close()
                    break
    except EOFError:
        star_db.close()
        raise EOFError('The %s you requested, is not available in the %s database.'%(model_id,code is None and 'requested' or code))
    return entry

##############################
#ADD PARAMETER TO MODEL IN DB#
##############################

def addParameterToModel(model_id,database_path,parameter,value):
    '''Add a single parameter to a model's parameter dictionary in a database. A backup of the old database is made!
    
    Input:    model_id=string, the model_id you're referring to
        database_path=string, the database
        parameter=string, name of the parameter as it appears in the inputfile for the respective code
        value=int/float/..., the value of the parameter
    '''
    subprocess.call(['mv ' + database_path + ' ' + database_path+'old'],shell=True,stdout=subprocess.PIPE)
    gastronoom_db_old = open(database_path+'old','r')
    gastronoom_db = open(database_path,'w')
    print 'Making a change to database at %s'%database_path
    while True:
        try:
            model = cPickle.load(gastronoom_db_old)
            if model_id == model[1]:
                print 'Adding the %s parameter with value %s to the dictionary of the model with id %s.'%(parameter,str(value),model_id)
                new_dict = model[0]
                new_dict[parameter] = value
                model = [new_dict,model_id]
            cPickle.dump(model,gastronoom_db)  
        except EOFError:
            break
    gastronoom_db_old.close()
    gastronoom_db.close()
    


def addParameterToMCMaxDb(database_path,parameter,value,last_par):
    
    '''
    Add a single parameter to all models' parameter dictionaries in an MCMax database, if it is not yet present. A backup of the old database is made!
    
    Input:    database_path=string, the database
        parameter=string, name of the parameter as it appears in the inputfile for the respective code
        value=int/float/..., the value of the parameter
        last_par=string, the parameter is inserted after this parameter
    '''
    subprocess.call(['mv ' + database_path + ' ' + database_path+'old'],shell=True,stdout=subprocess.PIPE)
    gastronoom_db_old = open(database_path+'old','r')
    gastronoom_db = open(database_path,'w')
    print 'Making a change to database at %s'%database_path
    while True:
        try:
            model = cPickle.load(gastronoom_db_old)
            model_id = model[1]
            new_dict = model[0]
            keys = [entry.split('=')[0] for entry in new_dict]
            index = keys.index(last_par)+1
            if keys[index] != parameter:    
                print 'Adding the %s parameter with value %s to the dictionary of the model with id %s.'%(parameter,str(value),model_id)
                new_dict[index:index] = ['%s=%s'%(parameter,str(value))]
            model = [new_dict,model_id]
            cPickle.dump(model,gastronoom_db)
        except EOFError:
            break
    gastronoom_db_old.close()
    gastronoom_db.close()
    

"""