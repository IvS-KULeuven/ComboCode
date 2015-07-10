# -*- coding: utf-8 -*-

"""
Interface for modeling work.

Author: R. Lombaert

"""

import os

import cc.path
from cc.modeling.codes.MCMax import MCMax
from cc.modeling.codes.Gastronoom import Gastronoom
from cc.tools.io import Database



class ModelingManager():
    
    """ 
    A modeling manager which maintains information for all modeling sessions.
    
    """
    
    def __init__(self,var_pars,processed_input,iterations=1,\
                 mcmax=0,gastronoom=0,sphinx=0,iterative=0,\
                 num_model_sessions=1,vic_manager=None,replace_db_entry=0,\
                 path_gastronoom='runTest',path_mcmax='runTest',\
                 skip_cooling=0,recover_sphinxfiles=0):
        
        """ 
        Initializing a ModelingManager instance.
        
        From this class, the MCMax and GASTRoNOoM codes are ran for each 
        iteration, as well as the database retrieval of older models. 
        
        @param var_pars: gridded parameters in the CC session
        @type var_pars: list[string]
        @param processed_input: The processed paramaters from the CC inputfile
        @type processed_input: dict
        @keyword iterations: Number of iterations
                             
                             (default: 1)
        @type iterations: int
        @keyword mcmax: Running MCMax?
        
                        (default: 0)
        @type mcmax: bool
        @keyword gastronoom: Running GASTRoNOoM?
        
                             (default: 0
        @type gastronoom: bool
        @keyword sphinx: Running Sphinx?
        
                         (default: 0)
        @type sphinx: bool
        @keyword iterative: Ray-trace MCMax models on every iteration
            
                            (default: 0)
        @type iterative: bool
        @keyword num_model_sessions: number of sessions, i.e. len(star_grid)
        
                                     (default: 1)
        @type num_model_sessions: int
        @keyword vic_manager: the vic manager to run models (sphinx) on VIC3
        
                              (default: None)
        @type vic_manager: Vic()
        @keyword replace_db_entry: replace an entry in the database with a 
                                   newly calculated model with a new model id 
                                   (eg if some general data not included in 
                                   the inputfiles is changed)
                                   
                                   (default: 0)
        @type replace_db_entry: bool
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
        @keyword path_mcmax: modeling folder in MCMax home
        
                             (default: 'runTest')
        @type path_mcmax: string
        @keyword path_gastronoom: modeling folder in GASTRoNOoM home
        
                                  (default: 'runTest')
        @type path_gastronoom: string
        
        """
        
        self.var_pars, self.iterations = var_pars, int(iterations)
        self.mcmax, self.gastronoom = int(mcmax), int(gastronoom)
        self.sphinx = int(sphinx)
        self.skip_cooling = skip_cooling
        self.input_dict = processed_input
        self.iterative = int(iterative)
        self.star_grid_old = [[] for i in xrange(num_model_sessions)]
        self.vic = vic_manager
        self.replace_db_entry = replace_db_entry
        self.new_entries_mcmax = []
        self.new_entries_cooling = []
        self.trans_bool_list = []
        self.mline_done_list = []
        self.mcmax_done_list = []
        self.path_mcmax = path_mcmax
        self.path_gastronoom = path_gastronoom
        self.recover_sphinxfiles = recover_sphinxfiles
        
        #-- Convenience paths
        cc.path.gout = os.path.join(cc.path.gastronoom,self.path_gastronoom)
        cc.path.mout = os.path.join(cc.path.mcmax,self.path_mcmax)
        self.setDatabases()
        
        #-- Initialize code booleans for checking if a code was ran
        self.mcmax_done = False
        self.mline_done = False
        
        
    def setDatabases(self):
         
        '''
        Initialize all databases relevant for this grid.
        
        '''
        
        if self.gastronoom:
            cool_db_path = os.path.join(cc.path.gout,\
                                        'GASTRoNOoM_cooling_models.db')
            ml_db_path = os.path.join(cc.path.gout,\
                                      'GASTRoNOoM_mline_models.db')
            sph_db_path = os.path.join(cc.path.gout,\
                                       'GASTRoNOoM_sphinx_models.db')
            self.cool_db = Database.Database(db_path=cool_db_path)
            self.ml_db = Database.Database(db_path=ml_db_path)
            self.sph_db = Database.Database(db_path=sph_db_path)       
        if self.mcmax:
            mcmax_db_path = os.path.join(cc.path.mout,'MCMax_models.db')
            self.mcmax_db = Database.Database(db_path=mcmax_db_path)
        if self.vic <> None:
            self.vic.setSphinxDb(self.sph_db)
        
        
        
    def startModeling(self,star,star_index):
        
        """ 
        Start the modeling process on a model star.
        
        @param star: The parameter set for this session
        @type star: Star()
        @param star_index: The index of the Star() object in the full list in 
                           CC. Only used to track earlier iterations if 
                           iterative==1
        @type star_index: int
        
        """
        
        #-- Note that db synchronisation is done every time a model is 
        #   successfully calculated, every time a previous code has been ran,
        #   and every time a change is made to the database, in terms of model
        #   id reservation. 
        
        for i in range(self.iterations):
            if self.mcmax: 
                print '***********************************'
                print '** Starting MCMax calculation.'
                print '** Iteration # %i for Model %i.'%(i+1,star_index+1)
                #-- Check if the previously calculated GASTRoNOoM model was
                #   calculated. This means a lot of time has passed since
                #   the last synchronisation of the mcmax db, so sync it
                if self.mline_done:
                    self.mcmax_db.sync()
                #-- Initiate a dust session which is used for every iteration
                if i == 0: 
                    dust_session = MCMax(path_mcmax=self.path_mcmax,\
                                         db=self.mcmax_db,\
                                         new_entries=self.new_entries_mcmax,\
                                         replace_db_entry=self.replace_db_entry)
                self.mcmax_done = False
                dust_session.doMCMax(star)
                if dust_session.mcmax_done: 
                    self.mcmax_done = True
                
                #- If last iteration, ray trace.
                if i+1 == self.iterations or self.iterative:
                    dust_session.rayTrace(star)
                     
                #- Remember every iteration if iterative==True.
                if self.iterative:
                    self.star_grid_old[star_index].append(star.copy())
                
                #- add/change MUTABLE input keys, which MAY be overwritten by
                #- the input file inputComboCode.dat
                #- Only relevant if a model match was found/calculated.
                if dust_session.model_id:
                    star.removeMutableMCMax(dust_session.mutable,self.var_pars)
                    star.update(self.input_dict)
                    if self.replace_db_entry:
                        self.new_entries_mcmax.append(dust_session.model_id)
                
            if self.gastronoom:    
                #- Initiate a gas session which is used for every iteration
                if i == 0: 
                    gas_session = Gastronoom(vic=self.vic,\
                                        path_gastronoom=self.path_gastronoom,\
                                        cool_db=self.cool_db,\
                                        ml_db=self.ml_db,\
                                        sph_db=self.sph_db,\
                                        sphinx=self.sphinx,\
                                        skip_cooling=self.skip_cooling,\
                                        replace_db_entry=self.replace_db_entry,\
                                        new_entries=self.new_entries_cooling,\
                                        recover_sphinxfiles=self.recover_sphinxfiles)
                    self.mline_done = False
                if self.mcmax_done:
                    #-- MCMax was ran successfully, in other words, quite a bit 
                    #   of time has passed, so update the cool database.
                    self.cool_db.sync()
                
                #-- If last iteration and no mline is requested, don't run 
                #   cooling. There is no point anyway (and means you don't have 
                #   to wait for your MCMax model. 
                if (i+1 == self.iterations) and not star['GAS_LIST']:
                    continue
                
                print '***********************************'
                print '** Starting GASTRoNOoM calculation.'
                print '** Iteration # %i for Model %i.'%(i+1,star_index+1)
                                
                #-- Otherwise, run cooling and do the rest of the loop
                gas_session.doGastronoom(star)
                
                #- add/change MUTABLE input keys, which MAY be overwritten by 
                #- the input file inputComboCode.dat
                star.removeMutableGastronoom(gas_session.mutable,self.var_pars)
                star.update(self.input_dict)
                star.updateMolecules(parlist=gas_session.mutable)
                
                #-- Remember the new model_ids if the replace_db_entry option 
                #   is on. 
                if self.replace_db_entry:
                    self.new_entries_cooling.append(gas_session.model_id)
                
                #- If last iteration, run mline and sphinx. Note that the 
                #- model_id cannot have changed: it's either the original 
                #- cooling model_id or no model_id at all in case of total fail
                if (i+1 == self.iterations) and gas_session.model_id: 
                    if gas_session.cool_done:
                        self.ml_db.sync()
                    gas_session.doMline(star)
                    if gas_session.mline_done:
                        self.mline_done = True
                    #- Check if the model id is still valid after the mline run
                    if gas_session.model_id:
                        gas_session.doSphinx(star)
                print '***********************************'
        
        #- remember trans bools if sphinx is enabled, so you can trace which 
        #- models have been calculated in this session and which were retrieved 
        #- from database. Also remember if mline was required to be ran this 
        #- session, or rather if everything was taken out of the db
        if self.gastronoom:
            self.mline_done_list.append(self.mline_done)
            if self.sphinx: 
                self.trans_bool_list.append(gas_session.trans_bools)
        
        #- remember if mcmax was ran, or rather the model taken from the db.
        #- In first or second iteration, doesn't matter.
        if self.mcmax: 
            self.mcmax_done_list.append(self.mcmax_done)
        
        #- old MCMax models are kept if iterative is requested. Before every
        #- change the model is saved to the logging list of models. In practice, 
        #- there will be lists in star_grid_old of Star models, in which every 
        #- model is the next step in the iteration, regardless of a mcmax or 
        #- gastronoom run. When plotting, each of these lists will be plotted 
        #- together showing the evolution of the modeling session for this 
        #- parameter set.
                
