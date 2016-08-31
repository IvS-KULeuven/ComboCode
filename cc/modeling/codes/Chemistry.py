# -*- coding: utf-8 -*-

"""
Running Chemistry and managing output from Chemistry.

Author: M. Van de Sande (based on MCMax.py of R. Lombaert)

"""

import os
import subprocess
from glob import glob

import cc.path
from cc.tools.io import DataIO, Database
from cc.modeling.codes.ModelingSession import ModelingSession


     

class Chemistry(ModelingSession):
    
    """ 
    Class that includes all methods required for creating an Chemistry model. 
    
    """
    
    def __init__(self,path_chemistry='runTest',replace_db_entry=0,db=None,\
                 single_session=0):
        
        """ 
        Initializing an instance of ModelingSession.
        
        @keyword db: the Chemistry database
        
                          (default: None)
        @type db: Database()
        @keyword replace_db_entry: replace an entry in the Chemistry database with 
                                   a newly calculated model with a new model id 
                                   (for instance if some general data not 
                                   included in the inputfiles is changed)
                                   
                                   (default: 0)
        @type replace_db_entry: bool
        @keyword path_chemistry: modeling folder in Chemistry home
        
                             (default: 'runTest')
        @type path_chemistry: string
        @keyword new_entries: The new model_ids when replace_db_entry is 1
                                   of other models in the grid. These are not 
                                   replaced!
                                   
                                   (default: [])
        @type new_entries: list[str]     
        @keyword single_session: If this is the only CC session. Speeds up db
                                 check.
                                 
                                 (default: 0)
        @type single_session: bool
                
        """
        
        super(Chemistry, self).__init__(code='Chemistry',path=path_chemistry,\
                                    replace_db_entry=replace_db_entry,\
                                    single_session=single_session)
        #-- Convenience path
        cc.path.cout = os.path.join(cc.path.chemistry,self.path)
        #DataIO.testFolderExistence(os.path.join(cc.path.mout,\
                                                #'data_for_gastronoom'))
        self.db = db
        
        #-- If an chemistry model is in progress, the model manager will hold until
        #   the other cc session is finished. 
        self.in_progress = False
        
        #- Read standard input file with all parameters that should be included
        #- as well as some dust specific information
        self.inputfilename = os.path.join(cc.path.aux,'inputChemistry.dat')
        self.standard_inputfile = DataIO.readDict(self.inputfilename,\
                                                  convert_floats=1,\
                                                  convert_ints=1,\
                                                  comment_chars=['#','*'])
        chemistry_keys = os.path.join(cc.path.aux,'Input_Keywords_Chemistry.dat')
        self.chemistry_keywords = [line.strip() 
                                   for line in DataIO.readFile(chemistry_keys) 
                                   if line]      



            
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
        
        keyword_int_list = ['PERFORM_ROUTINE']
        if comm_key.lower() in keyword_int_list: make_int = 1
        else: make_int = 0     
        return super(Chemistry, self).setCommandKey(comm_key,star,'DUST',\
                                                star_key,alternative,\
                                                make_int)


    
    
    
    def cCL(self,*args,**kwargs):
        
        ''' 
        Short-hand helper function for compareCommandLists.
        
        '''
        
        return self.compareCommandLists(*args,**kwargs) 
    
    
        
    def checkDatabase(self):

        """
        Checking Chemistry database.
        
        @return: The presence of the Chemistry model in the database
        @rtype: bool
        
        """
        
        #-- Lock the Chemistry database by opening it in read mode. It's closed
        #   once the database check is finalised. Note that in a case of a crash
        #   during the for loop, the python shell must be exited to unlock the 
        #   sphinx database again. The sync() is now done only once at the very
        #   end since the file on the disk will not change.
        if not self.single_session: self.db.sync()
        chem_dbfile = self.db._open('r')
        db_ids = sorted(self.db.keys())
        for i,model_id in enumerate(db_ids):
            chem_dict = self.db[model_id]
            model_bool = self.compareCommandLists(self.command_list.copy(),\
                                                  chem_dict, 'chemistry')
            if model_bool:
                if chem_dict.has_key('IN_PROGRESS'):
                    self.in_progress = True
                    print 'Chemistry model is currently being calculated in a ' +\
                          'different CC modeling session with ID %s'\
                          %(model_id)
                    self.model_id = model_id
                    finished = 1
                    break
                elif self.replace_db_entry \
                        and model_id not in self.new_entries: 
                    print 'Replacing Chemistry database entry for old ID %s'\
                          %model_id
                    del self.db[model_id]
                    finished = 0
                    break
                else:
                    print 'Chemistry model has been calculated ' + \
                          'before with ID %s'%model_id
                    self.model_id = model_id
                    finished = 1
                    break
        
            #-- Reached the end of db without match. Make new entry in db, in
            #   progress. Cant combine this with next line in case the last
            #   model gives a match.
            if i == len(self.db)-1:
                print 'No match found in Chemistry database. ' + \
                      'Calculating new model.'
                finished = 0
        
        #-- In case of an empty db, the above loop is not accessed.
        if not self.db.keys():
            print 'No match found in Chemistry database. Calculating new model.'
            finished = 0
        
        #-- Add the model in progress to the Chemistry db
        if finished == 0:    
            self.model_id = self.makeNewId()
            self.db[self.model_id] = self.command_list.copy()
            self.db[self.model_id]['IN_PROGRESS'] = 1
            
        #-- In case of an empty db, the above loop is not accessed.
        if not self.db.keys():
            print 'No match found in Chemistry database. Calculating new model.'
            finished = 0
        
        #-- Synchronize and unlock db.
        chem_dbfile.close()
        if not self.single_session: self.db.sync()
        return finished
        
        
            
    def doChemistry(self,star):
        
        """
        Running Chemistry.
        
        @param star: The parameter set for this session
        @type star: Star()
        
        """

        print '***********************************'                                       
        #- Create the input dictionary for this Chemistry run
        print '** Making input file for Chemistry'
        #-- Add the previous model_id to the list of new entries, so it does 
        #   not get deleted if replace_db_entry == 1. 
        if self.model_id: 
            self.new_entries.append(self.model_id)
        self.model_id = ''
        self.command_list = dict()
        # Bijzonder gevallen hier ipv in loop over orig
        # Dan is input_lines/command_list wat er in de database staat
        # input_lines = d
        
        if star['PERFORM_ROUTINE'] == 0:
            self.command_list['ROUTINE_RADIUS'] = 0
        else:
            self.command_list['ROUTINE_RADIUS'] = star['ROUTINE_RADIUS']
                #*star['R_STAR']*star.Rsun
        self.command_list['R_STAR'] = star['R_STAR']*star.Rsun
        #self.command_list['R_INNER_CHEM'] = star['R_INNER_CHEM']*\
            #star['R_STAR']*star.Rsun
        #self.command_list['R_OUTER_CHEM'] = star['R_OUTER_CHEM']*\
            #star['R_STAR']*star.Rsun
        
        self.command_list['REACTIONS_FILE'] = '"'+os.path.join(cc.path.csource,\
            'rates',star['REACTIONS_FILE'])+'"'
        self.command_list['SPECIES_FILE'] = '"'+os.path.join(cc.path.csource,\
            'specs',star['SPECIES_FILE'])+'"'
        self.command_list['FILECO'] = '"'+os.path.join(cc.path.csource,\
            'shielding',star['FILECO'])+'"'
        self.command_list['FILEN2'] = '"'+os.path.join(cc.path.csource,\
            star['FILEN2'])+'"'
        
        add_keys = [k  
                    for k in self.standard_inputfile.keys() 
                    if not self.command_list.has_key(k)]
        [self.setCommandKey(k,star,star_key=k.upper(),\
                            alternative=self.standard_inputfile[k])
         for k in add_keys]

        print '** DONE!'
        print '***********************************'
        
        #-- Check the Chemistry database if the model was calculated before
        modelbool = self.checkDatabase()
                
        #-- if no match found in database, calculate new model with new model id 
        #-- if the calculation did not fail, add entry to database for new model
        if not modelbool:
            input_dict = self.command_list.copy()
            input_lines = []
            orig = DataIO.readFile(self.inputfilename)
            for i,s in enumerate(orig):
                split = s.split()
                if s[0] == '!':
                    input_lines.append(s)
                else:
                    input_lines.append(" ".join(split[0:2])+' '+\
                        str(self.command_list[split[0]]))
            output_folder = os.path.join(cc.path.cout,'models',self.model_id)
            DataIO.testFolderExistence(output_folder)
            input_lines.append('OUTPUT_FOLDER = "' +\
                               output_folder+'/"') 
            
            input_filename = os.path.join(cc.path.cout,'models',\
                                          'inputChemistry_%s.txt'%self.model_id)
            
            DataIO.writeFile(filename=input_filename,input_lines=input_lines)
            
            #subprocess.call(' '.join([cc.path.ccode,input_filename]),shell=True)
            subprocess.call(' '.join([os.path.join(cc.path.csource,'csmodel'),input_filename]),shell=True)
            
            
            # files die worden aangemaakt op einde, test of successvol
            testf1 = os.path.join(output_folder,'cscoldens.out')
            testf2 = os.path.join(output_folder,'csfrac.out')
            if os.path.exists(testf1) and os.path.exists(testf2) and \
                    os.path.isfile(testf1) and os.path.isfile(testf2):
                del self.db[self.model_id]['IN_PROGRESS']
                self.db.addChangedKey(self.model_id)
            else:
                print '** Model calculation failed. No entry is added to ' + \
                      'the database.'
                del self.db[self.model_id]
                self.model_id = ''
            if not self.single_session: self.db.sync()

        #- Note that the model manager now adds/changes MUTABLE input keys, 
        #- which MAY be overwritten by the input file inputComboCode.dat
        print '***********************************'
        
