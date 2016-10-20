# -*- coding: utf-8 -*-

"""
Running GASTRoNOoM and managing output from GASTRoNOoM.

Author: R. Lombaert

"""

import os
import cPickle 
from glob import glob
import subprocess      
from scipy import array

import cc.path
from cc.tools.io import DataIO
from cc.tools.io import Atmosphere
from cc.modeling.codes.ModelingSession import ModelingSession
from cc.modeling.objects.Molecule import Molecule



class Gastronoom(ModelingSession):
    
    """ 
    Class that includes all methods required for creating a GATRoNOoM model.
    
    """
        
    def __init__(self,path_gastronoom='runTest',vic=None,sphinx=0,\
                 replace_db_entry=0,cool_db=None,ml_db=None,sph_db=None,\
                 skip_cooling=0,recover_sphinxfiles=0,\
                 new_entries=[],single_session=0):
    
        """ 
        Initializing an instance of a GASTRoNOoM modeling session.
        
        A single cooling model and a given set of molecules and transitions are
        calculated here or are retrieved from the databases.
        
        @keyword path_gastronoom: modeling folder in GASTRoNOoM home
        
                                  (default: 'runTest')
        @type path_gastronoom: string
        @keyword vic: the vic manager for running sphinx models on VIC3 
        
                      (default: None)
        @type vic: Vic()
        @keyword sphinx: Running Sphinx?
        
                         (default: 0)
        @type sphinx: bool
        @keyword replace_db_entry: replace an entry in the databases with a 
                                   newly calculated model with a new model id 
                                   (e.g. if some general data not included in 
                                   the inputfiles is changed)
                                   
                                   (default: 0)
        @type replace_db_entry: bool
        @keyword new_entries: The new model_ids when replace_db_entry is 1
                                   of other models in the grid. These are not 
                                   replaced!
                                   
                                   (default: [])
        @type new_entries: list[str]     
        @keyword skip_cooling: Skip running cooling in case a model is not 
                               found in the database, for instance if it is 
                               already known that the model will fail
        
                               (default: 0)
        @type skip_cooling: bool
        @keyword recover_sphinxfiles: Try to recover sphinx files from the disk
                                      in case they were correctly calculated, 
                                      but not saved to the database for one 
                                      reason or another. 
                                      
                                      (default: 0) 
        @type recover_sphinxfiles: bool
        @keyword cool_db: the cooling database
        
                          (default: None)
        @type cool_db: Database()
        @keyword ml_db: the mline database
        
                        (default: None)
        @type ml_db: Database()
        @keyword sph_db: the sphinx database
        
                         (default: None)
        @type sph_db: Database()
        @keyword single_session: If this is the only CC session. Speeds up db
                                 check.
                                 
                                 (default: 0)
        @type single_session: bool
                
        """
        
        super(Gastronoom,self).__init__(code='GASTRoNOoM',\
                                        path=path_gastronoom,\
                                        replace_db_entry=replace_db_entry,\
                                        new_entries=new_entries,\
                                        single_session=single_session)
        #-- Convenience path
        cc.path.gout = os.path.join(cc.path.gastronoom,self.path)
        self.vic = vic
        self.sphinx = sphinx
        cool_keys = os.path.join(cc.path.aux,'Input_Keywords_Cooling.dat')
        ml_keys = os.path.join(cc.path.aux,'Input_Keywords_Mline.dat')
        sph_keys = os.path.join(cc.path.aux,'Input_Keywords_Sphinx.dat')
        self.cooling_keywords = [line.strip() 
                                 for line in DataIO.readFile(cool_keys) 
                                 if line]
        self.mline_keywords = [line.strip() 
                               for line in DataIO.readFile(ml_keys) 
                               if line]
        self.sphinx_keywords = [line.strip() 
                                for line in DataIO.readFile(sph_keys) 
                                if line]
        DataIO.testFolderExistence(os.path.join(cc.path.gout,'data_for_mcmax'))
        self.trans_bools = []
        self.mline_done = False
        self.cool_done = False
        #-- If a cooling model is in progress, the model manager will hold until
        #   the other cc session is finished. 
        self.in_progress = False
        #-- Same for molecules.
        self.molec_in_progress = []
        #-- Transitions in progress are simply ignored.
        self.trans_in_progress = []
        self.cooling_molec_keys = ['ENHANCE_ABUNDANCE_FACTOR',\
                                   'ABUNDANCE_FILENAME',\
                                   'NUMBER_INPUT_ABUNDANCE_VALUES',\
                                   'KEYWORD_TABLE','MOLECULE_TABLE',\
                                   'ISOTOPE_TABLE']
        self.no_ab_molecs = ['12C16O','13C16O','1H1H16O','p1H1H16O',\
                             '1H1H17O','p1H1H17O','1H1H18O','p1H1H18O']
        
        #- Read standard input file with all parameters that should be included
        filename = os.path.join(cc.path.aux,'inputGASTRoNOoM.dat')
        self.standard_inputfile = DataIO.readDict(filename,\
                                                  comment_chars=['#','!'])
        self.skip_cooling = skip_cooling
        self.recover_sphinxfiles = recover_sphinxfiles
        self.cool_db = cool_db
        self.ml_db = ml_db
        self.sph_db = sph_db
        #self.pacs_db = pacs_db
        


    def addTransInProgress(self,trans):
        
        '''
        Remember a transition in progress. They will be checked at the 
        end of a VIC run to see if they have been correctly calculated. If being
        calculated in another CC session, they will not be checked. But the user
        can always reload.
        
        @param trans: The transition in progress
        @type trans: Transition()
        
        '''
        
        self.trans_in_progress.append(trans)



    def addMolecInProgress(self,molec):
        
        '''
        Remember a molecule in progress. The model manager will wait until it 
        is finished.
        
        @param molec: The molecule in progress
        @type molec: Molecule()
        
        '''
        
        self.molec_in_progress.append(molec)
        
        
        
    def execGastronoom(self,filename,subcode):
        
        '''
        Execute a subprocess of GASTRoNOoM: cooling, mline or sphinx
        
        @param filename: the full path+filename of the inputfile
        @type filename: string
        @param subcode: one of ['cooling','mline','sphinx'], the code to run
        @type subcode: string
        
        '''
  
        print '** Running %s...'%subcode
        if not subcode.lower() in ['cooling','mline','sphinx']:
            raise IOError('Subcode of GASTRoNOoM wrongly specified.')
        subprocess.call(['echo %s | %s'%(filename,subcode.lower())],shell=True)
        print '** DONE!'
        print '***********************************'


    
    def makeIdLog(self,new_id,molec_id=None):
        
        '''
        Make a text file with the original cooling id in it.
        
        This is used when creating a Transition() from a sphinx outputfilename.
        
        @param new_id: the new id for either mline or sphinx
        @type new_id: string
        @keyword molec_id: if given, an mline_id.log will be made including the
                           mline id to which molecules are linked
                           
                           (default: None)
        @type molec_id: string
        
        '''
        
        DataIO.writeFile(filename=os.path.join(cc.path.gout,'models',new_id,\
                         'cooling_id.log'),input_lines=[self.model_id])
        if not molec_id is None:
            DataIO.writeFile(filename=os.path.join(cc.path.gout,'models',\
                             new_id,'mline_id.log'),input_lines=[molec_id])



    def updateModel(self,model_id=None):

        """
        Updating model name in command list, eg in case mline or sphinx require
        a new model_id, different from the cooling model_id.
        
        @keyword model_id: if None the instance's model_id is taken, otherwise 
                           this one is used
                           
                           (default: None)
        @type model_id: string
        
        """
        
        for par in ['OUTPUT_DIRECTORY','PARAMETER_FILE','OUTPUT_SUFFIX']:
            self.command_list[par] = self.command_list[par].\
                    replace(self.command_list[par][self.command_list[par].\
                    find('model_'):self.command_list[par].find('model_')+25],\
                    not model_id is None and model_id or self.model_id,2)
        


    def deleteCoolingId(self,model_id):
        
        '''
        Delete everything from all databases that has anything to do with 
        this cooling id. 
        
        The databases on the hard disk are updated as well!
        
        @param model_id: the cooling model_id
        @type model_id: string
        
        '''
        
        print 'Replacing database entries for old cooling id %s.'%model_id
        del self.cool_db[model_id]
        try:
            del self.ml_db[model_id]
        except KeyError:
            pass
        try:
            del self.sph_db[model_id]
        except KeyError:
            pass



    def checkCoolingDatabase(self,molec_dict):

        """
        Checking cooling database.
        
        @param molec_dict: molecule info for this cooling model, ie CO and H2O
        @type molec_dict: dict()
        
        @return: The presence of the cooling model in the database
        @rtype: bool
        
        """
        
        #-- Lock the cooling database by opening it in read mode. It's closed
        #   once the database check is finalised. Note that in a case of a crash
        #   during the for loop, the python shell must be exited to unlock the 
        #   sphinx database again. The sync() is now done only once at the very
        #   end since the file on the disk will not change.
        if not self.single_session: self.cool_db.sync()
        cool_dbfile = self.cool_db._open('r')        
        for i,model_id in enumerate(sorted(self.cool_db.keys())):
            cool_dict = self.cool_db[model_id]
            model_bool = self.cCL(self.command_list.copy(),cool_dict,'cooling',\
                                  extra_dict=molec_dict)
            if model_bool:
                if cool_dict.has_key('IN_PROGRESS'):
                    self.in_progress = True
                    print 'Cooling model is currently being calculated in a ' +\
                          'different CC modeling session with ID %s.'\
                          %(model_id)
                    self.model_id = model_id
                    self.updateModel()
                    finished = 1
                    break
                elif self.replace_db_entry and model_id not in self.new_entries: 
                    self.deleteCoolingId(model_id)
                    print 'Deleting cooling model with ID %s from database.'\
                          %(model_id) + 'Calculating anew with id %s.'\
                          %self.model_id
                    finished = 0
                    break
                else:
                    print 'GASTRoNOoM cooling model has been calculated ' + \
                          'before with ID %s.'%model_id
                    self.model_id = model_id
                    self.updateModel()
                    finished = 1
                    break
            
            #-- Reached the end of db without match. Make new entry in db, in
            #   progress. Cant combine this with next line in case the last
            #   model gives a match.
            if i == len(self.cool_db)-1:
                print 'No match found in GASTRoNOoM cooling database. ' + \
                      'Calculating new model.'
                finished = 0
        
        #-- In case of an empty db, the above loop is not accessed.
        if not self.cool_db.keys():
            print 'No match found in GASTRoNOoM cooling database. ' + \
                  'Calculating new model.'
            finished = 0
        
        #-- Add the model in progress to the cooling db
        if finished == 0:    
            #-- OK to update. molec_dict is a copy. Because have to
            #   add the other input keywords for cooling to the H2O info. 
            #   This is saved to the db
            molec_dict.update(self.command_list)
            molec_dict['IN_PROGRESS'] = 1
            self.cool_db[self.model_id] = molec_dict
        
        #-- Synchronize and unlock db.
        cool_dbfile.close()
        if not self.single_session: self.cool_db.sync()
        return finished
                


    def checkMlineDatabase(self):
        
        """
        Check mline database.
        
        @return: The presence of mline models in the database, equivalent with 
                 self.molec_list
        @rtype: list[bool]
        
        """
        
        #-- Lock the mline database by opening it in read mode. It's closed
        #   once the database check is finalised. Note that in a case of a crash
        #   during the for loop, the python shell must be exited to unlock the 
        #   sphinx database again. The sync() is now done only once at the very
        #   end since the file on the disk will not change.
        if not self.single_session: self.ml_db.sync()
        ml_dbfile = self.ml_db._open('r')
        if not self.ml_db.has_key(self.model_id):
            self.ml_db[self.model_id] = dict([(self.model_id,dict())])
            for molec in self.molec_list:
                molec.setModelId(self.model_id)
                #-- Add an in-progress entry to the db
                md = molec.makeDict(in_progress=1)
                self.ml_db[self.model_id][self.model_id][molec.molecule] = md
            self.ml_db.addChangedKey(self.model_id)
            ml_dbfile.close()
            self.ml_db.sync()
            return [False]*len(self.molec_list)

        model_bools = []
        for molec in self.molec_list:
            for molec_id in [k for k,v in sorted(self.ml_db[self.model_id].items())
                               if molec.molecule in v.keys()]:
                db_molec_dict = self.ml_db[self.model_id][molec_id]\
                                          [molec.molecule]
                if self.cCL(this_list=molec.makeDict(),\
                            modellist=db_molec_dict,\
                            code='mline',\
                            ignoreAbun=molec.molecule in self.no_ab_molecs):
                    molec.setModelId(molec_id)
                    model_bools.append(True)
                    if db_molec_dict.has_key('IN_PROGRESS'):
                        self.addMolecInProgress(molec)
                        print 'Mline model is currently being ' + \
                              'calculated in a different CC modeling '+\
                              'session for %s with ID %s.'\
                              %(molec.molecule,molec.getModelId())                     
                    else:
                        print 'Mline model has been calculated before for '\
                              '%s with ID %s.'\
                              %(molec.molecule,molec.getModelId())         
                    break
                    
            #-- No match found. Calculate molecule anew. Now give it an id
            #   matching existing ids in terms of pars, or make a new one 
            if molec.getModelId() is None:
                model_bools.append(False)
                cool_dbd = self.ml_db[self.model_id]
                for i,(k,v) in enumerate(sorted(cool_dbd.items())):
                    #-- No molecule in this id yet. 
                    #   This should NEVER occur! Raise error.
                    #   Such empty dicts are removed upon fail/success check
                    #   in doMline.
                    if not v.keys(): 
                        raise KeyError('Empty molec id found in the database.'+\
                                       ' This should not be possible.')
                    
                    #-- molecule already in this id. Keep searching, unless it 
                    #   is the last entry. Then we have to make a new id
                    if molec.molecule in v.keys() \
                                and i != len(cool_dbd.keys())-1: 
                            continue
                        
                    #-- Now check if par dict is same except the molecule
                    #   definition and some molecule specific parameters. Only
                    #   check this specific keys list. If match, choose this ID
                    #   Make sure molec is not in v.keys() already (to avoid
                    #   the case of this being the last entry of the db)
                    if not molec.molecule in v.keys():
                        mm = v.keys()[0]
                        cks = ['OUTER_R_MODE','CHANGE_DUST_TO_GAS_FOR_ML_SP',\
                               'DUST_TO_GAS_CHANGE_ML_SP','STARFILE',\
                               'USE_STARFILE','USE_NO_MASER_OPTION',\
                               'USE_MASER_IN_SPHINX','FEHLER','N_FREQ',\
                               'START_APPROX','USE_FRACTION_LEVEL_CORR',\
                               'FRACTION_LEVEL_CORR','NUMBER_LEVEL_MAX_CORR']
                        if self.cCL(this_list=molec.makeDict(),\
                                    modellist=v[mm],code='mline',\
                                    check_keys=cks):
                            break
                    
                    #-- If end of dict is reached, no proper match was found
                    #   Create new id and put it in k. 
                    if i == len(cool_dbd.keys())-1: 
                        #-- First check if the original model id is already
                        #   in use (it may have been deleted if first use 
                        #   failed). Never copy output for this though, 
                        #   since it is the original 
                        if not cool_dbd.has_key(self.model_id):
                            k = self.model_id
                        else:                         
                            k = self.makeNewId()
                            self.makeIdLog(new_id=k)
                            self.copyOutput(molec,self.model_id,k)
                        self.ml_db[self.model_id][k] = dict()
                
                #-- It is possible cool_dbd is actually empty. Then use the
                #   model id as molec id. No need for an id log in this case
                #   nor any copying of output.
                if not cool_dbd.keys():
                    k = self.model_id
                    self.ml_db[self.model_id][self.model_id] = dict()
                
                #-- The last k value (or the new one in case end of dict was
                #   reached) is set as the model id.
                molec.setModelId(k)
                
                #-- Add an in-progress entry to the db
                md = molec.makeDict(in_progress=1)
                self.ml_db[self.model_id][k][molec.molecule] = md
                self.ml_db.addChangedKey(self.model_id)
                
                #-- Inform the user
                print 'Mline model for %s '%molec.molecule + \
                      'has not been calculated before. Calculate anew with '+ \
                      'ID %s.'%(molec.getModelId())
        ml_dbfile.close()
        if not self.single_session: self.ml_db.sync()
        return model_bools            



    def checkSphinxDatabase(self):
        
        """
        Check Sphinx database.
        
        The presence of the sphinx models in the database is saved in 
        self.trans_bools, equivalent to self.trans_list
        
        """
        
        #-- Remember which molecules have been added to new id, if applicable
        copied_molecs = []    

        #-- Lock the sphinx database by opening it in read mode. It's closed
        #   once the database check is finalised. Note that in a case of a crash
        #   during the for loop, the python shell must be exited to unlock the 
        #   sphinx database again. The sync() is now done only once at the very
        #   end since the file on the disk will not change.
        if not self.single_session: self.sph_db.sync()
        sph_dbfile = self.sph_db._open('r')
        if not self.sph_db.has_key(self.model_id):
            self.sph_db[self.model_id] = dict()
        for trans in self.trans_list:
            molec_id = trans.molecule.getModelId()
            molec = trans.molecule
            if not molec_id:
                self.trans_bools.append(False)
                trans.setModelId('')
            elif not self.sph_db[self.model_id].has_key(molec_id):
                trans.setModelId(molec_id)
                nd = dict([(str(trans),trans.makeDict(1))])
                self.sph_db[self.model_id][molec_id] = dict([(molec_id,nd)])
                self.sph_db.addChangedKey(self.model_id)
                self.trans_bools.append(False)
            else:    
                tdb = self.sph_db[self.model_id][molec_id]
                gd_ids = [k for k,v in sorted(tdb.items())
                            if str(trans) in v.keys()]
                for trans_id in gd_ids:
                    db_trans_dict = self.sph_db[self.model_id][molec_id]\
                                               [trans_id][str(trans)].copy()
                    if self.cCL(this_list=trans.makeDict(),\
                                modellist=db_trans_dict,code='sphinx'):
                        trans.setModelId(trans_id)
                        self.trans_bools.append(True)
                        if not self.vic is None \
                                and db_trans_dict.has_key('IN_PROGRESS'):
                            self.vic.addTransInProgress(trans)
                            print 'Sphinx model is currently being '+\
                                  'calculated for %s of %s with ID %s.'\
                                  %(str(trans),molec.molecule,\
                                    trans.getModelId())                                 
                        elif self.vic is None \
                              and db_trans_dict.has_key('IN_PROGRESS'):
                            self.addTransInProgress(trans)
                            print 'Sphinx model is currently being ' + \
                                  'calculated in a different CC modeling '+\
                                  'session for %s of %s with ID %s.'\
                                  %(str(trans),molec.molecule,\
                                    trans.getModelId())                     
                        else:
                            print 'Sphinx model has been calculated before '+ \
                                  'for %s of %s with ID %s.'\
                                  %(str(trans),molec.molecule,\
                                    trans.getModelId())
                        break
                        
                #-- No match found. Calculate trans anew. Now give it an id
                #   matching existing ids in terms of pars, or make a new one 
                if trans.getModelId() is None:
                    self.trans_bools.append(False)
                    molec_dbd = self.sph_db[self.model_id][molec_id]

                    for i,(k,v) in enumerate(sorted(molec_dbd.items())):
                        #-- No trans in this id yet.
                        #   This should NEVER occur! Raise error.
                        #   Such empty dicts are removed upon fail/success check
                        #   in doSphinx/checkSphinxOutput.
                        if not v.keys(): 
                            raise KeyError('Empty trans id found in the ' + \
                                           'database. This should not be ' + \
                                           'possible.')

                        #-- Trans already in this id. Keep searching, unless it 
                        #   is the last entry. Then we have to make a new id
                        if str(trans) in v.keys() \
                                and i != len(molec_dbd.keys())-1: 
                            continue         

                        #-- Now check if par dict is same except the transition
                        #   definition. Trick cCL by adding the db trans def to
                        #   this trans' dict. If match, choose this ID
                        #   Make sure trans is not in v.keys() already (to avoid
                        #   the case of this being the last entry of the db)
                        if not str(trans) in v.keys():
                            tt = v.keys()[0]
                            extra_dict={'TRANSITION':v[tt]['TRANSITION']}
                            if self.cCL(this_list=trans.makeDict(),\
                                        modellist=v[tt],code='sphinx',\
                                        extra_dict=extra_dict):
                                break
                                
                        #-- If end of dict is reached, no proper match was found
                        #   Create new id and put it in k. 
                        if i == len(molec_dbd.keys())-1: 
                            #-- First check if the original molec id is already
                            #   in use (it may have been deleted if first use 
                            #   failed). Never copy output for this though, 
                            #   since it is the original 
                            if not molec_dbd.has_key(molec_id):
                                k = molec_id
                                copied_molecs.append((molec.molecule,k))
                            else: 
                                k = self.makeNewId()
                                self.makeIdLog(new_id=k,molec_id=molec_id)
                            self.sph_db[self.model_id][molec_id][k] = dict()
                    
                    #-- It is possible molec_dbd is actually empty. Then use the
                    #   molec id as trans id. No need for an id log in this case
                    #   nor any copying of output.
                    if not molec_dbd.keys():
                        k = molec_id
                        copied_molecs.append((molec.molecule,k))
                        self.sph_db[self.model_id][molec_id][molec_id] = dict()
                        
                    #-- The last k value (or the new one in case end of dict was
                    #   reached) is set as the model id.
                    trans.setModelId(k)
                    
                    #-- When the trans is not yet calculated, either a new id
                    #   was made or an existing one is used. It still needs to 
                    #   be checked if all mline and cooling info is available. 
                    #   You only want to do this once per session for each 
                    #   molecule, because ls/ln checks add a lot of overhead. 
                    #   copyOutput double checks if links already exist
                    if (molec.molecule,k) not in copied_molecs:
                        self.copyOutput(trans,molec_id,k)
                        copied_molecs.append((molec.molecule,k)) 
                    td = trans.makeDict(1)
                    self.sph_db[self.model_id][molec_id][k][str(trans)] = td
                    self.sph_db.addChangedKey(self.model_id)

        sph_dbfile.close()
        if not self.single_session: self.sph_db.sync()
        


    def copyOutput(self,entry,old_id,new_id):
        
        '''
        Copy modelling output based on model_id.
        
        @param entry: the modeling object for which output is copied
        @type entry: Molecule() or Transition()
        @param old_id: The old model_id
        @type old_id: string
        @param new_id: the new_model_id
        @type new_id: string
        
        '''
        
        folder_old = os.path.join(cc.path.gout,'models',old_id)
        folder_new = os.path.join(cc.path.gout,'models',new_id)
        lsprocess = subprocess.Popen('ls %s'%folder_old,shell=True,\
                                     stdout=subprocess.PIPE)
        lsfile = lsprocess.communicate()[0].split('\n')
        lsfile = [os.path.split(line)[1] 
                     for line in lsfile 
                     if ((line[0:2] == 'ml' or line[0:4] == 'cool') \
                            and not entry.isMolecule()) \
                         or line[0:7] == 'coolfgr' \
                         or line[0:4] == 'para' \
                         or line[0:5] == 'input']
        if not entry.isMolecule():
            lsfile = [line 
                         for line in lsfile 
                         if not ((line[0:2] == 'ml' \
                                  and os.path.splitext(line)[0].split('_')[-1]\
                                      != entry.molecule.molecule))]
                                            
            lsfile = [line 
                         for line in lsfile 
                         if not (line[0:4] == 'cool' \
                            and (line.split('_')[-1].replace('.dat','') \
                                            != entry.molecule.molecule \
                            or line.split('_')[-1].replace('.dat','')=='sampling'\
                            or line[0:7] == 'coolfgr'))]
                             
        new_lsfile = [line.replace(old_id,new_id) for line in lsfile]
        DataIO.testFolderExistence(folder_new)
        lsprocess = subprocess.Popen('ls %s'%folder_new,shell=True,\
                                     stdout=subprocess.PIPE)
        already_done = lsprocess.communicate()[0].split('\n')
        for ls,nls in zip(lsfile,new_lsfile):
            if not nls in already_done:
                subprocess.call(['ln -s %s %s'%(os.path.join(folder_old,ls),\
                                                os.path.join(folder_new,nls))],\
                                shell=True)



    def doCooling(self,star):
        
        """
        Run Cooling.

        First, database is checked for retrieval of old model. 

        @param star: The parameter set for this session
        @type star: Star()
        
        """
        
        #-- Collect H2O and CO molecule definitions for inclusion in the 
        #   cooling inputfile. Also includes abundance_filename info for H2O if
        #   requested
        if not star.getMolecule('1H1H16O') is None:
            h2o_dict = star.getMolecule('1H1H16O').makeDict()
        else:            
            h2o_dict = Molecule('1H1H16O',45,45,648,50).makeDict()
        if not star.getMolecule('12C16O') is None:
            co_dict = star.getMolecule('12C16O').makeDict()
        else:
            co_dict = Molecule('12C16O',61,61,240,50).makeDict()
        
        #-- no abundance profiles should be possible for CO. 
        if co_dict.has_key('MOLECULE_TABLE'):
            raise IOError('CO cannot be attributed a custom abundance ' + \
                          'profile at this time.')
        
        #-- F_H2O is irrelevant if an abundance file is passed for oH2O
        if h2o_dict.has_key('MOLECULE_TABLE'):
            del self.command_list['F_H2O']
            
        #-- Collect all H2O molecular information important for cooling
        molec_dict = dict([(k,h2o_dict[k]) 
                            for k in self.cooling_molec_keys 
                            if h2o_dict.has_key(k)])

        #-- Check database: only include H2O extra keywords if 
        #   abundance_filename is present. CO can't have this anyway.
        model_bool = self.checkCoolingDatabase(molec_dict=molec_dict.copy())    
        
        #- Run cooling if above is False
        if not model_bool:
            DataIO.testFolderExistence(os.path.join(cc.path.gout,'models',\
                                                    self.model_id))
            commandfile = ['%s=%s'%(k,v) 
                           for k,v in sorted(self.command_list.items())
                           if k != 'R_POINTS_MASS_LOSS'] + \
                          ['####'] + \
                          ['%s=%s'%('MOLECULE',co_dict['MOLECULE'])] + \
                          ['%s=%s'%(k,h2o_dict[k]) 
                           for k in self.cooling_molec_keys + ['MOLECULE']
                           if h2o_dict.has_key(k)] + ['####']
            if self.command_list.has_key('R_POINTS_MASS_LOSS'):
                commandfile.extend(['%s=%s'%('R_POINTS_MASS_LOSS',v) 
                                    for v in self.command_list\
                                                    ['R_POINTS_MASS_LOSS']] + \
                                   ['####'])
            filename = os.path.join(cc.path.gout,'models',\
                                    'gastronoom_' + self.model_id + '.inp')
            DataIO.writeFile(filename,commandfile)
            if not self.skip_cooling:
                self.execGastronoom(subcode='cooling',filename=filename)
                self.cool_done = True
            if os.path.isfile(os.path.join(cc.path.gout,'models',\
                                           self.model_id,'coolfgr_all%s.dat'\
                                           %self.model_id)):
                #-- Note that there is no need to create a log parameter file
                #   since the inputfiles are still available, and they always
                #   contain all cooling keywords. (maybe change for h2o cooling)
                if self.cool_db[self.model_id].has_key('IN_PROGRESS'):
                    del self.cool_db[self.model_id]['IN_PROGRESS']
                    self.cool_db.addChangedKey(self.model_id)
            else:
                print 'Cooling model calculation failed. No entry is added '+ \
                      'to the database.'
                del self.cool_db[self.model_id]    
                self.model_id = ''
            if not self.single_session: self.cool_db.sync()                    



    def doMline(self,star):
        
        """
        Run mline.
        
        First, database is checked for retrieval of old models. 

        @param star: The parameter set for this session
        @type star: Star()
        
        """
        
        #-- Make sure to reset this in case an iteration between cooling and 
        #   mline is happening
        self.mline_done = False
        model_bools = self.checkMlineDatabase()
        del self.command_list['R_OUTER']
        del self.command_list['OUTER_R_MODE']
        for molec,model_bool in zip(self.molec_list,model_bools):
            if not model_bool:
                self.updateModel(molec.getModelId())
                commandfile = ['%s=%s'%(k,v) 
                               for k,v in sorted(self.command_list.items())
                               if k != 'R_POINTS_MASS_LOSS'] +\
                              ['####'] + \
                              ['%s=%s'%(k,v) 
                               for k,v in sorted(molec.makeDict().items())] +\
                              ['####']
                if self.command_list.has_key('R_POINTS_MASS_LOSS'):
                    commandfile.extend(['%s=%s'%('R_POINTS_MASS_LOSS',v) 
                                        for v in self.command_list\
                                                    ['R_POINTS_MASS_LOSS']] +\
                                       ['####'])
                filename = os.path.join(cc.path.gout,'models',\
                                        'gastronoom_%s.inp'%molec.getModelId())
                DataIO.writeFile(filename,commandfile)                
                self.execGastronoom(subcode='mline',filename=filename)
                self.mline_done=True
                path = os.path.join(cc.path.gout,'models',molec.getModelId())
                fns = 'ml*{}_{}.dat'.format(molec.getModelId(),molec.molecule)
                if len(glob(os.path.join(path,fns))) == 3:
                    #-- Remove in-progress entry.
                    if self.ml_db[self.model_id][molec.getModelId()]\
                            [molec.molecule].has_key('IN_PROGRESS'):
                        del self.ml_db[self.model_id][molec.getModelId()]\
                                      [molec.molecule]['IN_PROGRESS']
                    
                    #-- Write mline keywords not included in sph files but used 
                    #   in the database in an extra log file. Only do this if it
                    #   doesn't already exist. The parameters should be the same
                    #   for all transitions with this model id.
                    mlfn = 'mline_parameters_{}.log'.format(molec.molecule)
                    mlfn = os.path.join(path,mlfn)
                    if not os.path.isfile(mlfn):
                        #-- Add MOLECULE too. Cuz, why not. For TRANSITION, that
                        #   info is recreated from sph files. Not so for mline.
                        mlfile = ['{}={}'.format(k,v)
                                  for k,v in sorted(molec.makeDict().items())]
                        DataIO.writeFile(mlfn,mlfile)
                else:
                    del self.ml_db[self.model_id][molec.getModelId()]\
                                  [molec.molecule] 
                    #-- Remove the molecule id if it does not contain molecules
                    #   anymore. The id is thus unused.
                    if not self.ml_db[self.model_id][molec.getModelId()].keys():
                        del self.ml_db[self.model_id][molec.getModelId()]
                    print 'Mline model calculation failed for'\
                          '%s. No entry is added to the database.'\
                          %(molec.molecule)
                    molec.setModelId('')
                    
                #-- Synchronize db: Both when successful or failure. 
                self.ml_db.addChangedKey(self.model_id)
                if not self.single_session: self.ml_db.sync()
                
                
        if set([molec.getModelId() for molec in self.molec_list]) == set(['']):  
            #- no mline models calculated: stop GASTRoNOoM here
            self.model_id = ''
            print 'Mline model calculation failed for all requested ' + \
                  'molecules. Stopping GASTRoNOoM here!'
        else:        
            #- at least one molecule was successfully calculated, so start  
            #- Sphinx, hence if vic is requested, the cooling model_id can now  
            #- be added to the models list
            if not self.vic is None and self.sphinx: 
                #- add the command list to the vic models list
                self.vic.addModel(self.model_id,self.command_list)
            
   

    def doSphinx(self,star):
        
        """
        Run Sphinx.
        
        First, database is checked for retrieval of old models. 

        @param star: The parameter set for this session
        @type star: Star()
        
        """
        
        self.checkSphinxDatabase()
        print '%i transitions out of %i not yet calculated.'\
              %(len([boolean for boolean in self.trans_bools if not boolean]),\
                len(self.trans_bools))
        for i,(trans_bool,trans) in enumerate(zip(self.trans_bools,\
                                                  self.trans_list)):
            if not trans_bool and trans.getModelId():
                if not self.sphinx:
                    #- Only transitions with no db entry will get empty model id
                    del self.sph_db[self.model_id][trans.molecule.getModelId()]\
                                   [trans.getModelId()][str(trans)]
                    self.sph_db.addChangedKey(self.model_id)
                    trans.setModelId('')
                elif not self.vic is None:
                    #- add transition to the vic translist for this cooling id
                    self.vic.addTrans(trans)
                elif self.recover_sphinxfiles: 
                    self.checkSphinxOutput(trans)
                else:
                    self.updateModel(trans.getModelId())
                    commandfile = ['%s=%s'%(k,v) 
                                   for k,v in sorted(self.command_list.items()) 
                                   if k != 'R_POINTS_MASS_LOSS'] + ['####'] + \
                                  ['%s=%s'%(k,v) 
                                   for k,v in sorted(trans.molecule.makeDict()\
                                                                .items())] + \
                                  ['####'] + \
                                  ['%s=%s'%(k,v) 
                                   for k,v in sorted(trans.makeDict()\
                                                                .items())] + \
                                  ['######']
                    if self.command_list.has_key('R_POINTS_MASS_LOSS'):
                        commandfile.extend(['%s=%s'%('R_POINTS_MASS_LOSS',v) 
                                            for v in self.command_list\
                                                    ['R_POINTS_MASS_LOSS']] + \
                                           ['####'])
                    filename = os.path.join(cc.path.gout,'models',\
                                            'gastronoom_%s.inp'\
                                            %trans.getModelId())
                    DataIO.writeFile(filename,commandfile)                
                    print 'Starting calculation for transition %i out of %i.'\
                          %(i+1,len(self.trans_bools))
                    self.execGastronoom(subcode='sphinx',filename=filename)
                    self.checkSphinxOutput(trans)
                    if not self.single_session: self.sph_db.sync()
                    
        #- check if at least one of the transitions was calculated: then 
        #- self.model_id doesnt have to be changed
        self.finalizeSphinx() 
        
        #-- Sync the sphinx db in case self.sphinx is False or 
        #   self.recover_sphinxfiles is True, to make sure sph_db is up-to-date.
        #   In case models are calculated or vic is ran this is done in other
        #   places in the code.
        if not self.sphinx or self.recover_sphinxfiles:
            self.sph_db.sync()
            
        mline_not_available = set([trans.molecule.getModelId() 
                                   for boolean,trans in zip(self.trans_bools,\
                                                            self.trans_list) 
                                   if not boolean])  \
                              == set([''])
        if not self.vic is None and self.sphinx and (False in self.trans_bools \
              and not mline_not_available):
            self.vic.queueModel()
        elif not self.vic is None and self.sphinx and (False not in self.trans_bools\
              or mline_not_available):
            self.vic.reset()
            
 
 
    def finalizeSphinx(self):
        
        '''
        Check if at least one of the transitions has been calculated correctly. 
        
        if not, self.model_id is set to "". 
        
        '''
        
        if set([trans.getModelId() for trans in self.trans_list]) == set(['']):
            self.model_id = ''
        for trans in self.trans_in_progress:
            if not self.sph_db[self.model_id][trans.molecule.getModelId()]\
                              [trans.getModelId()].has_key(str(trans)) \
                    or self.sph_db[self.model_id][trans.molecule.getModelId()]\
                                  [trans.getModelId()][str(trans)]\
                                  .has_key('IN_PROGRESS'):
                trans.setModelId('')  
                


    def checkSphinxOutput(self,trans):
        
        '''
        Check if sphinx output is complete and update the database with 
        calculated models. 
        
        Requires model_id and path defined in Gastronoom instance.
        
        @param trans: the transition that is being checked
        @type trans: Transition()
        
        '''
        
        filename = trans.makeSphinxFilename(number='*')
        path = os.path.join(cc.path.gout,'models',trans.getModelId())
        #- Sphinx puts out 2 files per transition
        if len(glob(os.path.join(path,filename))) == 2:                    
            if self.sph_db[self.model_id][trans.molecule.getModelId()]\
                          [trans.getModelId()][str(trans)]\
                          .has_key('IN_PROGRESS'):
                del self.sph_db[self.model_id][trans.molecule.getModelId()]\
                               [trans.getModelId()][str(trans)]['IN_PROGRESS']

            #-- Write sphinx keywords not included in sph filenames but used in
            #   the database in an extra log file. Only do this if it doesn't 
            #   already exist. The parameters should be the same for all 
            #   transitions with this model id.
            sphfn = os.path.join(path,'sphinx_parameters.log')
            if not os.path.isfile(sphfn):
                sphfile = ['{}={}'.format(k,v)
                           for k,v in sorted(trans.makeDict().items())
                           if k != 'TRANSITION']
                DataIO.writeFile(sphfn,sphfile)
            
            #-- Print to the shell that the line is finished.
            print 'Sphinx model calculated successfully for '+\
                  '%s of %s with id %s.'%(str(trans),trans.molecule.molecule,\
                                          trans.getModelId())                                                             
        else:
            del self.sph_db[self.model_id][trans.molecule.getModelId()]\
                           [trans.getModelId()][str(trans)]
            #-- Remove the transition id if it does not contain transitions
            #   anymore. The id is thus unused.
            if not self.sph_db[self.model_id][trans.molecule.getModelId()]\
                              [trans.getModelId()].keys():
                del self.sph_db[self.model_id][trans.molecule.getModelId()]\
                               [trans.getModelId()]
            print 'Sphinx model calculation failed for %s of %s with id %s.'\
                  %(str(trans),trans.molecule.molecule,trans.getModelId())
            print 'No entry is added to the Sphinx database.'
            trans.setModelId('')
        
        self.sph_db.addChangedKey(self.model_id)
        
        
    def setCommandKey(self,comm_key,star,star_key=None,alternative=None):
        
        '''
        Try setting a key in the command_list from a star instance. 
        
        If the key is unknown, it is left open and will be filled in from the 
        standard gastronoom inputfile.
        
        @param comm_key: the name of the keyword in the command list
        @type comm_key: string
        @param star: Parameter set for this session
        @type star: Star()
        
        @keyword star_key: the name of the keyword in Star() (minus '_%s'
                           %key_type (DUST or GAS), which is added as well in a 
                           second attempt if the first without the addition is 
                           not found). If None, it's equal to comm_key
                           
                           (default: None)
        @type star_key: string
        @keyword alternative: a default value passed from the standard 
                              inputfile that is used if the keyword or the 
                              keyword + '_%s'%key_type is not found in Star().
                              
                              (default: None)
        @type alternative: string
        
        @return: Success? 
        @rtype: bool
        
        '''    
        
        keyword_int_list = ['ITERA_COOLING','LOG_DEPTH_STEP_POWER',\
                            'USE_MLINE_COOLING_RATE_CO',\
                            'USE_NEW_DUST_KAPPA_FILES','STEP_RIN_ROUT',\
                            'STEP_RS_RIN']
        exp_not_list = ['STEP_RIN_ROUT','STEP_RS_RIN']
        make_int = comm_key in keyword_int_list
        
        #- Make sure large integers are given in exponential notation
        #- Fortran cannot work with very large integers.
        exp_not = comm_key in exp_not_list
        
        return super(Gastronoom, self).setCommandKey(comm_key,star,'GAS',\
                                                     star_key,alternative,\
                                                     make_int,exp_not)
 


    def doGastronoom(self,star):
        
        """
        Run GASTRoNOoM-cooling. 
        
        The input parameter list is prepared here.
        
        @param star: Parameter set for this session
        @type star: Star()
        
        """

        print '***********************************'
        print '** Making input file for GASTRoNOoM'
        
        #-- Add the previous cooling model_id to the list of new entries, so it
        #   does not get deleted if replace_db_entry == 1. 
        #   This id is not known in new_entries, as the new_entries are passed
        #   for the previous models, not the current one.(ie when iterations>1)
        if self.model_id: 
            self.new_entries.append(self.model_id)
        
        #-- Make sure to reset this in case an iteration between cooling and 
        #   mline, cooling and mcmax is happening
        self.cool_done = False
        self.model_id = self.makeNewId()
        self.trans_list=star['GAS_LINES']    
        self.molec_list=star['GAS_LIST']
        self.command_list = dict()
        self.command_list['DATA_DIRECTORY'] = '"' + cc.path.gdata + '"'
        self.command_list['OUTPUT_DIRECTORY'] \
                = '"' + os.path.join(cc.path.gout,'models',self.model_id) +'/"'
        self.command_list['PARAMETER_FILE'] \
                = '"' + os.path.join(cc.path.gout,'models',self.model_id,\
                                     'parameter_file_%s.dat'%self.model_id)+'"'
        self.command_list['OUTPUT_SUFFIX'] = self.model_id
        
        if star['TEMDUST_FILENAME'] == 'temdust.kappa':
            self.command_list['SPEC_DENS_DUST'] = 3.3
        else:
            self.command_list['SPEC_DENS_DUST'] = star['SPEC_DENS_DUST']        
        
        self.command_list['DUST_TEMPERATURE_FILENAME'] \
                = '"%s"'%star['DUST_TEMPERATURE_FILENAME']
        self.command_list['TEMDUST_FILENAME'] = '"%s"'%star['TEMDUST_FILENAME']
        self.command_list['R_STAR'] = float(star['R_STAR'])*star.Rsun
        
        if star['MDOT_MODE'].upper() != 'CONSTANT':
            self.command_list['USE_DENSITY_NON_CONSISTENT'] \
                    = str(int(star['USE_DENSITY_NON_CONSISTENT']))
            if star['MDOT_MODE'].upper() == 'OTHER':
                try:
                    self.command_list['R_POINTS_MASS_LOSS'] \
                            = star['R_POINTS_MASS_LOSS']
                except KeyError:
                    raise KeyError('R_POINTS_MASS_LOSS not present. Include '+\
                                   'this parameter in your inputfile or set '+\
                                   'MDOT_MODE to something different than OTHER.')
        
        #- always has to be included, 0.5 if not present in the inputfile, 
        #- also sets power law innermost region
        self.command_list['TEMPERATURE_EPSILON'] \
                = star['TEMPERATURE_EPSILON_GAS']    
                
        #- the next few keywords only if the mode is not cooling, but epsilon
        if star['TEMPERATURE_MODE_GAS'] != 'cooling':        
            if float(star['TEMPERATURE_EPSILON2_GAS']):
                self.command_list['TEMPERATURE_EPSILON2'] \
                        = star['TEMPERATURE_EPSILON2_GAS']
                self.command_list['RADIUS_EPSILON2'] \
                        = star['RADIUS_EPSILON2_GAS']
                #- Only add third power law if the second power law is non-zero
                if float(star['TEMPERATURE_EPSILON3_GAS']):     
                    self.command_list['TEMPERATURE_EPSILON3'] \
                            = star['TEMPERATURE_EPSILON3_GAS']
                    self.command_list['RADIUS_EPSILON3'] \
                            = star['RADIUS_EPSILON3_GAS']    
        
        #-- Add the stochastic velocity filename if mode is filename, or power 
        #   and inner value if mode is powerlaw. Mode and Vel are always added.
        if star['STOCHASTIC_VEL_MODE'].upper() == 'POWERLAW':
            self.command_list['STOCHASTIC_VEL_POWER'] \
                    = star['STOCHASTIC_VEL_POWER']
            self.command_list['STOCHASTIC_VEL_INNER'] \
                    = star['STOCHASTIC_VEL_INNER']
        elif star['STOCHASTIC_VEL_MODE'].upper() == 'FILENAME':
            self.command_list['STOCHASTIC_VEL_FILENAME'] \
                    = star['STOCHASTIC_VEL_FILENAME']
        
        self.setCommandKey('DUST_TO_GAS',star,'DUST_TO_GAS_INITIAL',\
                           self.standard_inputfile['DUST_TO_GAS'])
        
        #- Both pars set explicitly because they are also present in 
        #- mline_keywords
        self.setCommandKey('R_OUTER',star,'R_OUTER',\
                           self.standard_inputfile['R_OUTER'])     
        self.setCommandKey('OUTER_R_MODE',star,'OUTER_R_MODE',\
                           self.standard_inputfile['OUTER_R_MODE'])
        add_keys = [k  
                    for k in self.standard_inputfile.keys() 
                    if not (self.command_list.has_key(k) \
                        or k in self.mline_keywords + self.sphinx_keywords)]
        [self.setCommandKey(k,star,alternative=self.standard_inputfile[k]) 
         for k in add_keys]
        print '** DONE!'
        print '***********************************'
        
        #- model/output naming entries are excluded when comparing new models
        #- with database models
        #- Start the model calculation
        self.doCooling(star)

        #- Removing mutable input is done in ModelingManager now, as well as 
        #- starting up sphinx and mline...
        
