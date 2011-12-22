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

from cc.tools.io import DataIO
from cc.modeling.ModelingSession import ModelingSession
from cc.modeling.objects.Molecule import Molecule



class Gastronoom(ModelingSession):
    
    """ 
    Class that includes all methods required for creating a GATRoNOoM model.
    
    """
        
    def __init__(self,path_combocode=os.path.join(os.path.expanduser('~'),\
                                                  'ComboCode'),\
                 path_gastronoom='runTest',vic=None,sphinx=0,\
                 replace_db_entry=0,cool_db=None,ml_db=None,sph_db=None,\
                 pacs_db=None):
    
        """ 
        Initializing an instance of a GASTRoNOoM modeling session.
        
        A single cooling model and a given set of molecules and transitions are
        calculated here or are retrieved from the databases.
        
        @keyword path_gastronoom: modeling folder in GASTRoNOoM home
        
                                  (default: 'runTest')
        @type path_gastronoom: string
        @keyword path_combocode: CC home folder
        
                                 (default: '/home/robinl/ComboCode')
        @type path_combocode: string
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
        @keyword cool_db: the cooling database
        
                          (default: None)
        @type cool_db: Database()
        @keyword ml_db: the mline database
        
                        (default: None)
        @type ml_db: Database()
        @keyword sph_db: the sphinx database
        
                         (default: None)
        @type sph_db: Database()
        @keyword pacs_db: the pacs database
        
                          (default: None)
        @type pacs_db: Database()
        
        """
        
        super(Gastronoom,self).__init__(code='GASTRoNOoM',\
                                        path=path_gastronoom,\
                                        path_combocode=path_combocode,\
                                        replace_db_entry=replace_db_entry)
        self.vic = vic
        self.trans_in_progress = []
        self.sphinx = sphinx
        cool_keys = os.path.join(self.path_combocode,'CC',\
                                 'Input_Keywords_Cooling.dat')
        ml_keys = os.path.join(self.path_combocode,'CC',\
                               'Input_Keywords_Mline.dat')
        sph_keys = os.path.join(self.path_combocode,'CC',\
                                'Input_Keywords_Sphinx.dat')
        self.cooling_keywords = [line.strip() 
                                 for line in DataIO.readFile(cool_keys) 
                                 if line]
        self.mline_keywords = [line.strip() 
                               for line in DataIO.readFile(ml_keys) 
                               if line]
        self.sphinx_keywords = [line.strip() 
                                for line in DataIO.readFile(sph_keys) 
                                if line]
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                   'GASTRoNOoM',self.path,'data_for_mcmax'))
        self.trans_bools = []
        self.done_mline = False
        self.cooling_molec_keys = ['ENHANCE_ABUNDANCE_FACTOR',\
                                   'ABUNDANCE_FILENAME',\
                                   'NUMBER_INPUT_ABUNDANCE_VALUES',\
                                   'KEYWORD_TABLE','MOLECULE_TABLE',\
                                   'ISOTOPE_TABLE']
        self.no_ab_molecs = ['12C16O','13C16O','1H1H16O','p1H1H16O',\
                             '1H1H17O','p1H1H17O','1H1H18O','p1H1H18O']
        
        #- Read standard input file with all parameters that should be included
        filename = os.path.join(os.path.expanduser("~"),\
                                'GASTRoNOoM','inputGASTRoNOoM.dat')
        self.standard_inputfile = DataIO.readDict(filename,\
                                                  comment_chars=['#','!'])
        self.cool_db = cool_db
        self.ml_db = ml_db
        self.sph_db = sph_db
        self.pacs_db = pacs_db



    def addTransInProgress(self,trans):
        
        '''
        Remember a transition in progress on VIC. They will be checked at the 
        end of the VIC run to see if they have been correctly calculated.
        
        @param trans: The transition in progress
        @type trans: Transition()
        
        '''
        
        self.trans_in_progress.append(trans)



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
        exec_path = os.path.join(os.path.expanduser('~'),'GASTRoNOoM','src',\
                                 'exec',subcode.lower())
        subprocess.call(['echo %s | %s'%(filename,exec_path)],shell=True)
        print '** DONE!'
        print '***********************************'


    
    def makeIdLog(self, new_id,molec_id=None):
        
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
        
        DataIO.writeFile(filename=os.path.join(os.path.expanduser('~'),\
                         'GASTRoNOoM',self.path,'models',new_id,\
                         'cooling_id.log'),input_lines=[self.model_id])
        if molec_id <> None:
            DataIO.writeFile(filename=os.path.join(os.path.expanduser('~'),\
                             'GASTRoNOoM',self.path,'models',new_id,\
                             'mline_id.log'),input_lines=[molec_id])



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
                    model_id <> None and model_id or self.model_id,2)
        


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
        for pacs_id in [k for k,v in self.pacs_db.items() 
                          if v['cooling_id'] == model_id]:
            del self.pacs_db[pacs_id]  



    def checkCoolingDatabase(self,molec_dict):

        """
        Checking cooling database.
        
        @param molec_dict: molecule info for this cooling model, ie CO and H2O
        @type molec_dict: dict()
        
        @return: The presence of the cooling model in the database
        @rtype: bool
        
        """
        
        for model_id,cool_dict in self.cool_db.items():
            model_bool = self.compareCommandLists(self.command_list.copy(),\
                                                  cool_dict,'cooling',\
                                                  extra_dict=molec_dict)
            if model_bool:
                if self.replace_db_entry: 
                    self.deleteCoolingId(model_id)
                    return False
                else:
                    print 'GASTRoNOoM cooling model has been calculated ' + \
                          'before with ID %s.'%model_id
                    self.model_id = model_id
                    self.updateModel()
                    return True
        print 'No match found in GASTRoNOoM cooling database. Calculating ' + \
              'new model.'
        return False



    def checkMlineDatabase(self):
        
        """
        Check mline database.
        
        @return: The presence of mline models in the database, equivalent with 
                 self.molec_list
        @rtype: list[bool]
        
        """
        
        if not self.ml_db.has_key(self.model_id):
            [molec.setModelId(self.model_id) for molec in self.molec_list]
            self.ml_db[self.model_id] = dict([(self.model_id,dict())])
            return [False]*len(self.molec_list)
        model_bools = []
        new_molec_id = ''
        for molec in self.molec_list:
            for molec_id in [k for k,v in self.ml_db[self.model_id].items() 
                               if molec.molecule in v.keys()]:
                if self.compareCommandLists(this_list=molec.makeDict(),\
                                         modellist=self.ml_db[self.model_id]\
                                                             [molec_id]\
                                                             [molec.molecule],\
                                            code='mline',\
                                            ignoreAbun=molec.molecule \
                                                        in self.no_ab_molecs):
                    molec.setModelId(molec_id)
                    model_bools.append(True)
                    print 'Mline model has been calculated before for %s'\
                           %molec.molecule + \
                           ' with ID %s.'%(molec.getModelId())         
                    break
            if molec.getModelId() is None:
                model_bools.append(False)
                [molec.setModelId(k) 
                 for k,v in sorted(self.ml_db[self.model_id].items())
                 if molec.getModelId() is None \
                    and molec.molecule not in v.keys()]
                #- still not defined, so make new model_id
                if molec.getModelId() is None:    
                    if not new_molec_id:
                        new_molec_id = self.makeNewId()
                        self.makeIdLog(new_id=new_molec_id)
                        self.copyOutput(molec,self.model_id,new_molec_id)
                        self.ml_db[self.model_id][new_molec_id] = dict()
                        self.ml_db.addChangedKey(self.model_id)
                    molec.setModelId(new_molec_id)    
                print 'Mline model for %s '%molec.molecule + \
                      'has not been calculated before. Calculate anew with '+ \
                      'ID %s.'%(molec.getModelId())
        return model_bools            



    def checkSphinxDatabase(self):
        
        """
        Check Sphinx database.
        
        The presence of the sphinx models in the database is saved in 
        self.trans_bools, equivalent to self.trans_list
        
        """
        
        if not self.sph_db.has_key(self.model_id):
            self.sph_db[self.model_id] = dict()
        new_trans_id = ''
        #- Remember which molecules have been added to new id, if applicable
        molecules_copied_to_new_id = []    
        for trans in self.trans_list:
            molec_id = trans.molecule.getModelId()
            if not molec_id:
                self.trans_bools.append(False)
                trans.setModelId('')
            elif not self.sph_db[self.model_id].has_key(molec_id):
                trans.setModelId(molec_id)
                self.sph_db[self.model_id][molec_id] = \
                    dict([(molec_id,dict([(str(trans),trans.makeDict(1))]))])
                self.sph_db.addChangedKey(self.model_id)
                self.trans_bools.append(False)
            else:    
                for trans_id in [k for k,v in self.sph_db[self.model_id]\
                                                         [molec_id].items() 
                                   if str(trans) in v.keys()]:
                    db_trans_dict = self.sph_db[self.model_id][molec_id]\
                                               [trans_id][str(trans)].copy()
                    if self.compareCommandLists(this_list=trans.makeDict(),\
                                                modellist=db_trans_dict,\
                                                code='sphinx'):
                        trans.setModelId(trans_id)
                        self.trans_bools.append(True)
                        if self.vic <> None \
                                and db_trans_dict.has_key('IN_PROGRESS'):
                            self.vic.addTransInProgress(trans)
                            print 'Sphinx model is currently being '+\
                                  'calculated for %s of %s with ID %s.'\
                                  %(str(trans),trans.molecule.molecule,\
                                    trans.getModelId())                                 
                        elif self.vic is None \
                              and db_trans_dict.has_key('IN_PROGRESS'):
                            self.addTransInProgress(trans)
                            print 'Sphinx model is currently being ' + \
                                  'calculated in a different CC modeling '+\
                                  'session for %s of %s with ID %s.'\
                                  %(str(trans),trans.molecule.molecule,\
                                    trans.getModelId())                     
                        else:
                            print 'Sphinx model has been calculated before '+ \
                                  'for %s of %s with ID %s.'\
                                  %(str(trans),trans.molecule.molecule,\
                                    trans.getModelId())
                        break
                if trans.getModelId() is None:
                    self.trans_bools.append(False)
                    [trans.setModelId(k)
                     for k,v in sorted(self.sph_db[self.model_id][molec_id]\
                                          .items())
                     if trans.getModelId() is None \
                        and str(trans) not in v.keys()]
                    if trans.getModelId() is None:
                        if not new_trans_id:
                            new_trans_id = self.makeNewId()
                            self.makeIdLog(new_id=new_trans_id,\
                                           molec_id=molec_id)
                            self.sph_db[self.model_id][molec_id][new_trans_id]\
                                    = dict()
                        if trans.molecule not in molecules_copied_to_new_id:
                            self.copyOutput(trans,molec_id,new_trans_id)
                            molecules_copied_to_new_id.append(trans.molecule)
                        trans.setModelId(new_trans_id)
                    self.sph_db[self.model_id][molec_id][trans.getModelId()]\
                               [str(trans)] = trans.makeDict(1)
                    self.sph_db.addChangedKey(self.model_id)
                


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
        
        folder_old = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                  self.path,'models',old_id)
        folder_new = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                  self.path,'models',new_id)
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
                         if not (line[0:2] == 'ml' \
                            and line.split('_')[-1].replace('.dat','') \
                                            != entry.molecule.molecule)]
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
        
        if star.getMolecule('1H1H16O') <> None:
            h2o_dict = star.getMolecule('1H1H16O').makeDict()
        else:            
            h2o_dict = Molecule('1H1H16O',45,45,648,50,\
                    path_combocode=self.path_combocode).makeDict()
        if star.getMolecule('12C16O') <> None:
            co_dict = star.getMolecule('12C16O').makeDict()
        else:
            co_dict = Molecule('12C16O',61,61,240,50,\
                    path_combocode=self.path_combocode).makeDict()
        
        #- The check for abundance profiles for both h2o & co (it is impossible 
        #- for cooling to have both at the same time) is done in ComboCode.py
        molec_dict = dict([(k,h2o_dict[k]) 
                                    for k in self.cooling_molec_keys 
                                    if h2o_dict.has_key(k)] \
                            + [(k,co_dict[k]) 
                                    for k in self.cooling_molec_keys 
                                    if co_dict.has_key(k)])
        if h2o_dict.has_key('MOLECULE_TABLE') \
                and co_dict.has_key('MOLECULE_TABLE'):
            raise IOError('WARNING! CO/H2O both have customized abundances. '+\
                          'Only the CO customized abundance will be used ' + \
                          'for cooling.')
        
        #- Check database
        model_bool = self.checkCoolingDatabase(molec_dict=molec_dict)    
        
        #- Run cooling if above is False
        if not model_bool:
            DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                       self.code,self.path,'models',\
                                       self.model_id))
            commandfile = ['%s=%s'%(k,v) 
                           for k,v in sorted(self.command_list.items())
                           if k != 'R_POINTS_MASS_LOSS'] + ['####'] + \
                          ['%s=%s'%(k,co_dict[k]) 
                           for k in self.cooling_molec_keys + ['MOLECULE'] 
                           if co_dict.has_key(k)] + \
                          ['%s=%s'%(k,h2o_dict[k]) 
                           for k in self.cooling_molec_keys + ['MOLECULE'] 
                           if h2o_dict.has_key(k)] + ['####']
            if self.command_list.has_key('R_POINTS_MASS_LOSS'):
                commandfile.extend(['%s=%s'%('R_POINTS_MASS_LOSS',v) 
                                    for v in self.command_list\
                                                    ['R_POINTS_MASS_LOSS']] + \
                                   ['####'])
            filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                    self.path,\
                                    'gastronoom_' + self.model_id + '.inp')
            DataIO.writeFile(filename,commandfile)                
            self.execGastronoom(subcode='cooling',filename=filename)
            if os.path.isfile(os.path.join(os.path.expanduser("~"),\
                                           'GASTRoNOoM',self.path,'models',\
                                           self.model_id,'coolfgr_all%s.dat'\
                                           %self.model_id)):
                molec_dict.update(self.command_list)                  
                self.cool_db[self.model_id] = molec_dict
            else:
                print 'Cooling model calculation failed. No entry is added '+ \
                      'to the database.'
                self.model_id = ''
                


    def doMline(self,star):
        
        """
        Run mline.
        
        First, database is checked for retrieval of old models. 

        @param star: The parameter set for this session
        @type star: Star()
        
        """
        
        model_bools = self.checkMlineDatabase()
        del self.command_list['R_OUTER']
        del self.command_list['OUTER_R_MODE']
        for molec,model_bool in zip(self.molec_list,model_bools):
            if not model_bool:
                self.updateModel(molec.getModelId())
                commandfile = ['%s=%s'%(k,v) 
                               for k,v in sorted(self.command_list.items())
                               if k != 'R_POINTS_MASS_LOSS'] + ['####'] + \
                              ['%s=%s'%(k,v) 
                               for k,v in sorted(molec.makeDict().items())] +\
                              ['####']
                if self.command_list.has_key('R_POINTS_MASS_LOSS'):
                    commandfile.extend(['%s=%s'%('R_POINTS_MASS_LOSS',v) 
                                        for v in self.command_list\
                                                    ['R_POINTS_MASS_LOSS']] +\
                                       ['####'])
                filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                        self.path,'gastronoom_%s.inp'\
                                        %molec.getModelId())
                DataIO.writeFile(filename,commandfile)                
                self.execGastronoom(subcode='mline',filename=filename)
                self.done_mline=True
                if len([f for f in glob(os.path.join(os.path.expanduser("~"),\
                                        'GASTRoNOoM',self.path,'models',\
                                        molec.getModelId(),'ml*%s_%s.dat'\
                                        %(molec.getModelId(),molec.molecule)))])\
                        == 3:
                    self.ml_db[self.model_id][molec.getModelId()]\
                              [molec.molecule] = molec.makeDict()
                    self.ml_db.addChangedKey(self.model_id)
                    self.ml_db.sync()
                else:
                    print 'Mline model calculation failed for'\
                          '%s. No entry is added to the database.'\
                          %(molec.molecule)
                    molec.setModelId('')
        if set([molec.getModelId() for molec in self.molec_list]) == set(['']):  
            #- no mline models calculated: stop GASTRoNOoM here
            self.model_id = ''
            print 'Mline model calculation failed for all requested ' + \
                  'molecules. Stopping GASTRoNOoM here!'
        else:        
            #- at least one molecule was successfully calculated, so start  
            #- Sphinx, hence if vic is requested, the cooling model_id can now  
            #- be added to the models list
            if self.vic <> None and self.sphinx: 
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
        for trans_bool,trans in zip(self.trans_bools,self.trans_list):
            if not trans_bool and trans.getModelId():
                if not self.sphinx:
                    #- Only transitions with no db entry will get empty model id
                    del self.sph_db[self.model_id][trans.molecule.getModelId()]\
                                   [trans.getModelId()][str(trans)]
                    self.sph_db.addChangedKey(self.model_id)
                    trans.setModelId('')
                elif self.vic <> None:
                    #- add transition to the vic translist for this cooling id
                    self.vic.addTrans(trans)
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
                    filename = os.path.join(os.path.expanduser('~'),\
                                            'GASTRoNOoM',self.path,\
                                            'gastronoom_%s.inp'\
                                            %trans.getModelId())
                    DataIO.writeFile(filename,commandfile)                
                    self.execGastronoom(subcode='sphinx',filename=filename)
                    self.checkSphinxOutput(trans)
                    self.sph_db.sync()
                    
        #- check if at least one of the transitions was calculated: then 
        #- self.model_id doesnt have to be changed
        self.finalizeSphinx() 
        mline_not_available = set([trans.molecule.getModelId() 
                                   for boolean,trans in zip(self.trans_bools,\
                                                            self.trans_list) 
                                   if not boolean])  \
                              == set([''])
        if self.vic <> None and self.sphinx and (False in self.trans_bools \
              and not mline_not_available):
            self.vic.queueModel()
        elif self.vic <> None and self.sphinx and (False not in self.trans_bools\
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
        #- Sphinx puts out 2 files per transition
        if len(glob(os.path.join(os.path.expanduser("~"),'GASTRoNOoM',\
                self.path,'models',trans.getModelId(),filename))) == 2:                    
            if self.sph_db[self.model_id][trans.molecule.getModelId()]\
                          [trans.getModelId()][str(trans)]\
                          .has_key('IN_PROGRESS'):
                del self.sph_db[self.model_id][trans.molecule.getModelId()]\
                               [trans.getModelId()][str(trans)]['IN_PROGRESS']
            self.sph_db.addChangedKey(self.model_id)
            print 'Sphinx model calculated successfully for '+\
                  '%s of %s with id %s.'%(str(trans),trans.molecule.molecule,\
                                          trans.getModelId())                                                             
        else:
            del self.sph_db[self.model_id][trans.molecule.getModelId()]\
                           [trans.getModelId()][str(trans)]
            self.sph_db.addChangedKey(self.model_id)            
            print 'Sphinx model calculation failed for %s of %s with id %s.'\
                  %(str(trans),trans.molecule.molecule,trans.getModelId())
            print 'No entry is added to the Sphinx database.'
            trans.setModelId('')
        
        
        
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
                            'USE_NO_MASER_OPTION','USE_MLINE_COOLING_RATE_CO',\
                            'USE_NEW_DUST_KAPPA_FILES']
        if comm_key in keyword_int_list: make_int = 1
        else: make_int = 0     
        return super(Gastronoom, self).setCommandKey(comm_key,star,'GAS',\
                                                     star_key,alternative,\
                                                     make_int)
 


    def doGastronoom(self,star):
        
        """
        Run GASTRoNOoM-cooling. 
        
        The input parameter list is prepared here.
        
        @param star: Parameter set for this session
        @type star: Star()
        
        """

        print '***********************************'
        print '** Making input file for GASTRoNOoM'
        self.model_id = self.makeNewId()
        self.trans_list=star['GAS_LINES']    
        self.molec_list=star['GAS_LIST']
        self.command_list = dict()
        self.command_list['DATA_DIRECTORY'] \
                = '"' + os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                     'src','data/') + '"'
        self.command_list['OUTPUT_DIRECTORY'] \
                = '"' + os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                     self.path,'models',self.model_id) + '/"'
        self.command_list['PARAMETER_FILE'] \
                = '"' + os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                     self.path,'models',self.model_id,\
                                     'parameter_file_%s.dat'%self.model_id)+'"'
        self.command_list['OUTPUT_SUFFIX'] = self.model_id
        
        if star['TEMDUST_FILENAME'] == 'temdust.kappa':
            self.command_list['SPEC_DENS_DUST'] = 3.3
        else:
            self.command_list['SPEC_DENS_DUST'] = star['SPEC_DENS_DUST']        
        
        self.command_list['DUST_TEMPERATURE_FILENAME'] \
                = '"%s"'%star['DUST_TEMPERATURE_FILENAME']
        self.command_list['TEMDUST_FILENAME'] = '"%s"'%star['TEMDUST_FILENAME']
        self.command_list['R_STAR'] = float(star['R_STAR'])*star.r_solar
        
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
                        or k in self.mline_keywords + self.sphinx_keywords + \
                                ['TEMPERATURE_EPSILON2','RADIUS_EPSILON2'])]
        [self.setCommandKey(k,star,alternative=self.standard_inputfile[k]) 
         for k in add_keys]
        print '** DONE!'
        print '***********************************'
        
        #- model/output naming entries are excluded when comparing new models
        #- with database models
        #- Start the model calculation
        self.doCooling(star)
        star['LAST_GASTRONOOM_MODEL'] = self.model_id
        #- Removing mutable input is done in ModelingManager now, as well as 
        #- starting up sphinx and mline...
        
