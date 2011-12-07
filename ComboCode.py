# -*- coding: utf-8 -*-

"""
Main code for running GASTRoNOoM (L.Decin) and MCMax (M.Min).

Author: R. Lombaert

"""

import sys
import os
import types
import subprocess
import time



def checkEntryInfo(input_list,number_of_keys,info_type):
    
    '''
    Specific input keywords for ComboCode (currently: MOLECULE, TRANSITION,
    R_POINTS_MASS_LOSS) require multiple arguments separated by a space. This 
    method sorts those arguments and checks if the format is OK.
    
    The method returns a tuple or a list of these arguments depending on 
    multiplicity of input, and the need to create a grid for the parameter.
    
    @param input_list: argument strings for one of the special CC input keys.
    @type input_list: list[string]
    @param number_of_keys: the number of arguments expected
    @type number_of_keys: int
    @param info_type: MOLECULE, TRANSITION or R_POINTS_MASS_LOSS
    @type info_type: string
    @return: the sorted arguments for the requested info_type. The output is 
             given as a list or a tuple, depending on if this gridded parameter
             or not, respectively.
    @rtype: list[tuple(list[string])] or tuple(list[string])
    
    '''
    
    if not info_type in ('MOLECULE','R_POINTS_MASS_LOSS','TRANSITION'):
        raise KeyError('The info_type keyword is unknown. Should be ' + \
                       'MOLECULE, TRANSITION or R_POINTS_MASS_LOSS.')
    input_list = [line.split() for line in input_list]
    if set([len(line) for line in input_list]) != set([number_of_keys]):
        print 'Number of keys should be: %i'%number_of_keys
        print '\n'.join(['%i  for  %s'%(len(line),line) for line in input_list])
        raise IOError('Input for one of the %s lines has wrong number of ' + \
                      'values. Double check, and abort.'%info_type)
    else:
        #- if MOLECULE: only molecule string, 
        #- if R_POINTS_MASS_LOSS: only grid id number, 
        #- if TRANSITION: everything except last entry (n_quad)
        entries = [info_type in ('MOLECULE','R_POINTS_MASS_LOSS') \
                        and line[0] or ' '.join(line[0:-1]) 
                   for line in input_list]  
        unique_entries = set(entries)
        #- ie the defining parameter is never multiply defined ! 
        #- Hence, only a single set of parameters is used here.
        if len(unique_entries) == len(entries):        
            if info_type == 'R_POINTS_MASS_LOSS':
                input_list = ['  '.join(il[1:]) for il in input_list]
            return tuple(input_list)
        else:
            if info_type == 'TRANSITION':
                indices = [i 
                           for i,entry in enumerate(entries) 
                           if entries.count(entry) > 1]
                print 'Identical transition(s) in the transition list: Doubl'+\
                      'es will not be removed even if N_QUAD is the same too!'
                for i in indices:
                    print 'At index %i:  %s' %(i,entries[i])
                raw_input('Abort if identical transitions are not expected.')
            if info_type == 'R_POINTS_MASS_LOSS':
                #- This will be a list of R_POINTS_MASS_LOSS sets, where each 
                #- set is defined as a list of radial grid point parameters
                final = []  
                while input_list:
                    if int(input_list[0][0]) != 1:
                        raise IOError('The grid point ID numbers for the ' + \
                                      'R_POINTS_MASS_LOSS keywords do not ' + \
                                      'follow a correct order. Use ID = 1 ' + \
                                      'as a reset for new set of grid points.')
                    final.append([input_list.pop(0)])
                    while input_list and int(input_list[0][0]) != 1:
                        final[-1].append(input_list.pop(0))
                #- Remove the grid point ID numbers, not relevant anymore
                #- put the rest back into a string
                final = [tuple(['  '.join([num 
                                           for i,num in enumerate(this_point) 
                                           if i != 0]) 
                                for this_point in this_set]) 
                         for this_set in final]      
                return final
            final = [[line 
                      for line,entry in zip(input_list,entries) 
                      if entries.count(entry) == 1]]
            multiples = set([entry 
                             for entry in entries 
                             if entries.count(entry) != 1])
            for entry in multiples:
                this_input = [line 
                              for this_entry,line in zip(entries,input_list) 
                              if this_entry == entry]
                final = [this_list + [extra_line] 
                         for extra_line in this_input 
                         for this_list in final]
            return [tuple(final_list) for final_list in final]



from cc.tools.io import DataIO
from cc.tools.numerical import Gridding
from cc.managers import ModelingManager
from cc.managers import PlottingManager
from cc.managers import Vic
from cc.modeling.objects import Star
from cc.modeling.objects import Transition
from cc.statistics import PeakStats
from cc.statistics import Statistics
from cc.data.instruments import Pacs
from cc.data.instruments import Spire
from cc.data import Sed
from cc.modeling.tools import IceFitter
from cc.modeling.tools import WaterAbundance as wa



class ComboCode(object):
    
    '''
    The interface with which to run the ComboCode package.
    
    '''

    def __init__(self,inputfilename):
        
        '''
        Initializing a ComboCode instance. 
        
        Once this is done, you only need to run startSession(). Then all 
        methods in this class will be called according to your inputfile. Only
        run separate methods of the class if you know what you are doing!
        
        Input is read and parsed, and the parameter objects (Star()) are set 
        
        The instrument and data objects, and the plotting manager are set.
        
        The .spec file is updated here, if requested.
        
        The inputfile can be given on the command line as: 
        python ComboCode.py inputComboCode.dat
        
        In the python or ipython shell you can do:
        >>> import ComboCode
        >>> cc = ComboCode.ComboCode('/home/robinl/ComboCode/CC/inputComboCode.dat')
        
        @param inputfilename: The name of the inputfile. If path is not
                              included, it's set to the default path to CC:
                              ~/ComboCode/. Inputfiles should ideally be in 
                              this folder.
        @type inputfilename: string
        
        '''
        
        self.path_combocode = os.path.join(os.path.expanduser('~'),'ComboCode')
        if not os.path.split(inputfilename)[0]:
            self.inputfilename = os.path.join(self.path_combocode,\
                                              os.path.split(inputfilename)[1])
        else:
            self.inputfilename = inputfilename
        self.readInput()
        self.setGlobalPars()
        self.setLinks()
        self.setPacs()
        self.setSpire()
        self.setSed()
        self.setPlotManager()
        self.setVarPars()
        self.createStarGrid()
        #- Only the extra transition pars will differ across the grid, so grab
        #- the transition list from one of the Star() objects
        if self.update_spec: 
            Transition.updateLineSpec(self.star_grid[0]['GAS_LINES'])
        self.finished = False
        
    
    
    def startSession(self):
         
        '''
        Start a ComboCode session, based on the input read upon initialisation.
        
        The supercomputer and model managers are set and ran.
        
        The plot manager, statistics module, fitter modules are ran if 
        requested.
        
        The session ends by printing some info about the Star() objects.
        
        Once started, the ComboCode object cannot be started again. You will
        have to re-initialize. This will change in the future. 
        
        '''
        
        if not self.finished:    
            self.setVicManager()
            self.setModelManager()
            self.finished = True
            self.runModelManager()
            self.finalizeVic()
            self.runPlotManager()
            self.runStatistics()
            self.doIceFit()
            self.appendResults() 
            if self.write_dust_density:
                [star.writeDensity() for star in self.star_grid]
            self.printStarInfo()
        else:
            print 'This CC session is already finished. Please, create a new one.'
                
                

    def setGlobalPars(self):
         
        '''
        Set the global parameters for this CC session.
        
        '''
        
        default_global = [('mcmax',1),('gastronoom',1),('sphinx',1),('vic',0),\
                          ('iterations',2),('plot_iterative',0),\
                          ('vic_account','vsc30226'),('statistics',0),\
                          ('vic_time_per_sphinx',30),('vic_credits',None),\
                          ('append_results',0),('write_dust_density',0),\
                          ('replace_db_entry',0),('update_spec',0),\
                          ('stat_tolerance',1.),('stat_sigma',None),\
                          ('ln_path_gastronoom',''),('ln_path_mcmax',''),\
                          ('path_gastronoom',''),('path_mcmax',''),\
                          ('print_model_info',1),('stat_mode','chi2'),\
                          ('do_ice_fit',0),('cfg_ice_fit','')]
        global_pars = dict([(k,self.processed_input.pop(k.upper(),v)) 
                            for k,v in default_global])
        self.__dict__.update(global_pars)
        self.star_name = self.processed_input['STAR_NAME']
        if not self.gastronoom or not self.mcmax: self.iterations = 1 
        if (not self.path_mcmax and self.mcmax):
            raise IOError('Please define PATH_MCMAX in your inputfile.')
        if (not self.path_gastronoom and self.gastronoom): 
            raise IOError('Please define PATH_GASTRONOOM in your inputfile.')
                
                
                
    def setPacs(self):
        
        '''
        Collect the PACS relevant parameters from the inputfile and set the PACS
        object.
        
        '''
        
        path_pacs = self.processed_input.pop('PATH_PACS','')
        redo_convolution = self.processed_input.pop('PACS_REDO_CONVOLUTION',0)
        searchstring = self.processed_input.pop('PACS_SEARCHSTRING','')
        oversampling = self.processed_input.pop('PACS_OVERSAMPLING','')
        intrinsic = self.processed_input.pop('PACS_INTRINSIC',1)
        if path_pacs:
            self.pacs = Pacs.Pacs(star_name=self.star_name,\
                                  path_combocode=self.path_combocode,\
                                  path=self.path_gastronoom,\
                                  redo_convolution=redo_convolution,\
                                  oversampling=oversampling,\
                                  path_pacs=path_pacs,\
                                  intrinsic=intrinsic)
            self.pacs.setData(searchstring=searchstring)
        else:
            self.pacs = None
            
            
            
    def setSpire(self):
        
        '''
        Collect the SPIRE relevant parameters from the inputfile and set the SPIRE
        object.
        
        '''
        
        path_spire = self.processed_input.pop('PATH_SPIRE','')
        searchstring = self.processed_input.pop('SPIRE_SEARCHSTRING','')
        resolution = self.processed_input.pop('SPIRE_RESOLUTION',0)
        intrinsic = self.processed_input.pop('SPIRE_INTRINSIC',1)
        if path_spire:
            self.spire = Spire.Spire(star_name=self.star_name,\
                                     path_combocode=self.path_combocode,\
                                     path=self.path_gastronoom,\
                                     resolution=resolution,\
                                     path_spire=path_spire,\
                                     intrinsic=intrinsic)
            self.spire.setData(searchstring=searchstring)
        else:
            self.spire = None



    def setSed(self):
        
        '''
        Collect the SED data and create an Sed() object.
        
        For now done based on a Star() object, but will be changed to the same
        method as used for PACS and SPIRE. Star() eventually will not know 
        anything about data
        
        '''
        
        path_sed = self.processed_input.pop('PATH_SED',None)
        if path_sed:
            cc_path = os.path.join(self.path_combocode,'Data')
            star_index = DataIO.getInputData(cc_path).index(self.star_name)
            ak = DataIO.getInputData(path=cc_path,keyword='A_K')[star_index]
            longitude = DataIO.getInputData(path=cc_path,keyword='LONG')[star_index]
            latitude = DataIO.getInputData(path=cc_path,keyword='LAT')[star_index]
            if abs(longitude) < 5.0 and abs(latitude) < 5.0:
                gal_position = 'GC'
            else:
                gal_position = 'ISM'          
            self.sed = Sed.Sed(star_name=self.star_name,ak=ak,\
                               gal_position=gal_position,\
                               path_combocode=self.path_combocode,\
                               path=path_sed)
        else: 
            self.sed = None
        


    def setVarPars(self):
        
        '''
        Define the list of variable parameters in this CC session.
        
        '''
        
        self.var_pars = [k for k in self.multiplicative_grid.keys() + \
                                    self.additive_grid.keys() 
                           if k[:5] != 'T_MAX']


    def readInput(self):
        
        '''
        Read input for ComboCode and return list.
        
        The MOLECULE, TRANSITION and R_POINTS_MASS_LOSS parameter formats are
        checked for errors in this method. If erroneous, an IOError is raised.
        
        '''
        
        multi_keys = ['MOLECULE','TRANSITION','R_POINTS_MASS_LOSS']
        input_dict = DataIO.readDict(self.inputfilename,convert_floats=1,\
                                     multi_keys=multi_keys)
        #-- keywords in multi_keys require different method
        self.processed_input = dict()
        self.multiplicative_grid = dict()
        self.additive_grid = dict()
        molecules = input_dict.pop('MOLECULE',[])
        transitions = input_dict.pop('TRANSITION',[])
        r_points_mass_loss = input_dict.pop('R_POINTS_MASS_LOSS',[])      
        for k,v in input_dict.items():
            #-- Determine delimiter        
            try:
                if v.find('&') != -1: delimiter = '&'                         
                elif v.find(';') != -1: delimiter = ';'
                elif v.find(',') != -1: delimiter = ','
                elif v.find(':') != -1: delimiter = ':'
                elif v.find('*') != -1: delimiter = '&'
                else: delimiter = ' '                                                        
            except AttributeError: 
                #-- v is already a float, so can't use .find on it => no grids
                #-- no need to check the rest, continue on with the next k/v pair
                self.processed_input[k] = v
                continue
            #-- Expanding '*' entries
            newv = delimiter.join(\
                        [len(value.split('*')) > 1 
                                  and delimiter.join([value.split('*')[0]]*\
                                                      int(value.split('*')[1])) 
                                  or value 
                         for value in v.split(delimiter)])
            #-- Add entries to processed_input, the multiplicative grid or the
            #-- additive grid, depending on the type of delimiter.
            if delimiter == ' ':
                self.processed_input[k] = v 
            elif delimiter == ',':
                newv = [float(value) 
                        for value in newv.split(',')]
                newv = Gridding.makeGrid(*newv)
                self.multiplicative_grid[k] = newv
            else:
                try:
                    if delimiter == '&':
                        newv = tuple([float(value.rstrip()) 
                                      for value in newv.split('&')])
                        self.processed_input[k] = newv
                    elif delimiter == ';':
                        newv = [value.rstrip() == '%' and '%' or float(value.rstrip()) 
                                for value in newv.split(';')]
                        self.multiplicative_grid[k] = newv
                    elif delimiter == ':':
                        newv = [value.rstrip() == '%' and '%' or float(value.rstrip()) 
                                for value in newv.split(':')]
                        self.additive_grid[k] = newv
                except ValueError:
                    if delimiter == '&':
                        newv = tuple([value.rstrip() 
                                      for value in newv.split('&')])
                        self.processed_input[k] = newv
                    elif delimiter == ';':
                        newv = [value.rstrip() 
                                for value in newv.split(';')]
                        self.multiplicative_grid[k] = newv
                    elif delimiter == ':':
                        newv = [value.rstrip() 
                                for value in newv.split(':')]
                        self.additive_grid[k] = newv
        #-- Make sure R_POINTS_MASS_LOSS, Molecule and Transition input makes 
        #-- sense and is correct
        if molecules: 
            molecules = checkEntryInfo(molecules,20,'MOLECULE')
            if type(molecules) is types.ListType:
                self.multiplicative_grid['MOLECULE'] = molecules
            else:
                self.processed_input['MOLECULE'] = molecules
        if r_points_mass_loss: 
            r_points_mass_loss = checkEntryInfo(r_points_mass_loss,4,\
                                                            'R_POINTS_MASS_LOSS')
            if type(r_points_mass_loss) is types.ListType:
                self.additive_grid['R_POINTS_MASS_LOSS'] = r_points_mass_loss
            else:
                self.processed_input['R_POINTS_MASS_LOSS'] = r_points_mass_loss
        if transitions: 
            transitions = checkEntryInfo(transitions,12,'TRANSITION')
            if type(transitions) is types.ListType:
                self.multiplicative_grid['TRANSITION'] = transitions
            else:
                self.processed_input['TRANSITION'] = transitions
            
            
    
    def getStarGrid(self):
         
        '''
        Return the list of Star() objects for this ComboCode session.
        
        @return: The parameter Star() objects are returned.
        @rtype: list[Star()]
                  
        '''
         
        return self.star_grid
            
            
            
    def createStarGrid(self):
         
        '''
        Create a list of Star() objects based on the inputfile that has been
        parsed with cc.readInput().
        
        The list of Star() objects is saved in self.star_grid, and is accessed
        through cc.getStarGrid().
        
        '''
        
        base_star = Star.Star(example_star=self.processed_input,\
                              path_combocode=self.path_combocode,\
                              path_gastronoom=self.path_gastronoom,\
                              path_mcmax=self.path_mcmax)
        if self.additive_grid:
            grid_lengths = [len(v) for v in self.additive_grid.values()]
            if len(set(grid_lengths)) != 1:
                raise IOError('The explicit parameter declaration using <:> '+\
                              'has a variable amount of options (including ' +\
                              'the R_GRID_MASS_LOSS definition). Aborting...')
            else: 
                additive_dicts = [dict([(key,grid[index]) 
                                        for key,grid in self.additive_grid.items()])
                                  for index in xrange(grid_lengths[0])]
                self.star_grid = [Star.Star(example_star=base_star,\
                                         path_combocode=self.path_combocode,\
                                         path_gastronoom=self.path_gastronoom,\
                                         path_mcmax=self.path_mcmax,\
                                         extra_input=d)
                                  for d in additive_dicts]
        else:
            self.star_grid = [base_star]
        for key,grid in self.multiplicative_grid.items():
            self.star_grid = [Star.Star(path_combocode=self.path_combocode,\
                                        path_gastronoom=self.path_gastronoom,\
                                        path_mcmax=self.path_mcmax,\
                                        example_star=star,\
                                        extra_input=dict([(key,value)]))
                              for star in self.star_grid 
                              for value in grid]
        if self.processed_input.has_key('LAST_MCMAX_MODEL'):
            del self.processed_input['LAST_MCMAX_MODEL']
        if self.processed_input.has_key('LAST_GASTRONOOM_MODEL'):
            del self.processed_input['LAST_GASTRONOOM_MODEL']
        
        
    
    def setLinks(self):
        
        '''
        Set the output folder links if requested.
        
        If the folders/links already exist, nothing is done. 
        
        '''
        
        if self.ln_path_gastronoom:
            DataIO.checkLink(path=os.path.join(self.ln_path_gastronoom,\
                                               self.path_gastronoom),\
                             ln_path=os.path.join(os.path.expanduser('~'),\
                                                  'GASTRoNOoM'))
        else:
            DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                       'GASTRoNOoM',self.path_gastronoom))
        if self.ln_path_mcmax:
            DataIO.checkLink(path=os.path.join(self.ln_path_mcmax,\
                                               self.path_mcmax),\
                             ln_path=os.path.join(os.path.expanduser('~'),\
                                                  'MCMax'))
        else:
            DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                       'MCMax',self.path_mcmax))
        
        
    def setVicManager(self):
        
        '''
        Set up the VIC manager.
        
        '''
        
        if self.vic: 
            self.vic_manager = Vic.Vic(path=self.path_gastronoom,\
                                       path_combocode=self.path_combocode,\
                                       account=self.vic_account,\
                                       time_per_sphinx=self.vic_time_per_sphinx,\
                                       credits_acc=self.vic_credits)
            if self.update_spec: self.vic_manager.updateLineSpec()
        else: 
            self.vic_manager = None
        
        
    def setModelManager(self):
        
        '''
        Set up the model manager.
        
        '''    
        
        self.model_manager = ModelingManager.ModelingManager(\
                                       star_name=self.star_name,\
                                       iterations=self.iterations,\
                                       processed_input=self.processed_input,\
                                       var_pars=self.var_pars,\
                                       path_gastronoom=self.path_gastronoom,\
                                       mcmax=self.mcmax,\
                                       gastronoom=self.gastronoom,\
                                       sphinx=self.sphinx,\
                                       iterative=self.plot_iterative,\
                                       num_model_sessions=len(self.star_grid),\
                                       vic_manager=self.vic_manager,\
                                       path_combocode=self.path_combocode,\
                                       replace_db_entry=self.replace_db_entry,\
                                       path_mcmax=self.path_mcmax)
    
    
    def setPlotManager(self):    
        
        '''
        Set up the plot manager.
        
        '''
        
        plot_pars = dict([(k,self.processed_input.pop(k))
                          for k,v in self.processed_input.items()
                          if k[0:5] == 'PLOT_' or k[0:4] == 'CFG_'])
        self.plot_manager = PlottingManager.PlottingManager(\
                                         star_name=self.star_name,\
                                         gastronoom=self.gastronoom,\
                                         mcmax=self.mcmax,\
                                         path_gastronoom=self.path_gastronoom,\
                                         path_mcmax=self.path_mcmax,\
                                         inputfilename=self.inputfilename,\
                                         path_combocode=self.path_combocode,\
                                         pacs=self.pacs,\
                                         spire=self.spire,\
                                         sed=self.sed,\
                                         plot_pars=plot_pars)
        
        
        
    def runModelManager(self):
        
        '''
        Start up the modeling.
        
        '''    
        
        if self.gastronoom or self.mcmax:
            print '***********************************'
            print '** Starting grid calculation.'
            print '***********************************'
            for star_index, star in enumerate(self.star_grid):
                print '***********************************'
                print '** Model #%i out of %i requested models.'\
                      %(star_index+1,len(self.star_grid))
                print '***********************************'
                self.model_manager.startModeling(star,star_index)
                #-- done_mline is True if in previous model an mline calculation 
                #-- was done: Only then do a progress check, because a lot of time 
                #-- has passed, but then a wait time is used to make sure the newly
                #-- queued sphinx models after the mline model are properly queued.
                if self.vic \
                        and self.vic_manager.getQueue() \
                        and self.model_manager.done_mline_list[-1]:    
                    print '***********************************'
                    print '** Current VIC queue:'
                    print self.vic_manager.getQueue()
                    self.vic_manager.checkProgress(wait_qstat=1)



    def finalizeVic(self):
        
        '''
        At the end of a modeling session, allow Vic to be finalized and clean up.
        
        '''
        
        if self.vic:
            #vic_running = vic_manager.checkProgress()
            if self.vic_manager.getQueue():
                vic_running = True
            else:
                vic_running = False
            while vic_running:
                print 'VIC is not yet finished. Waiting 5 minutes before checking again.'
                print self.vic_manager.getQueue()
                try: 
                    time.sleep(300)
                except KeyboardInterrupt: 
                    print 'Ending wait time, continuing with progress check immediately.'
                vic_running = self.vic_manager.checkProgress()
            self.vic_manager.finalizeVic()



    def runPlotManager(self):
        
        '''
        Run the plotting manager.
        
        '''
        
        if self.plot_iterative: 
            print '************************************************'
            print '****** Plotting results for each iterative step of the SED.'
            print '************************************************'
            #- star_grid_old remembers all old models if iterative is on, 
            #- meaning that every list in star_grid_old consists of Star models 
            #- associated with one level of iteration. Following line plots all
            #- iterative steps for one star immutable parameter set
            [self.plot_manager.startPlotting(self.model_manager.star_grid_old[i],\
                                             i+1)
             for i in xrange(len(self.star_grid))]
        print '************************************************'
        print '****** Plotting final results.'
        print '************************************************'
        self.plot_manager.startPlotting(self.star_grid)    
        
        
        
    def appendResults(self):    
        
        '''
        Append results at the end of the inputfile.
        
        '''
        
        print '** Appending results to inputfile and copying to output folders.'
        print '***********************************'
        #-- Check if the transition was intended to be calculated, and if it was 
        #-- successful (ie don't add if it had already been done)
        timestring = '%.4i-%.2i-%.2ih%.2i-%.2i-%.2i'\
                      %(time.gmtime()[0],time.gmtime()[1],time.gmtime()[2],\
                        time.gmtime()[3],time.gmtime()[4],time.gmtime()[5])
        appendage = []
        if self.model_manager.trans_bool_list:
            model_ids_list = [list(set([(trans.molecule.molecule,\
                                         trans.getModelId())
                                        for boolean,trans in zip(trans_bool,\
                                                             star['GAS_LINES']) 
                                        if trans.getModelId() \
                                            and (not trans_bool \
                                            or self.append_results)]))
                              for star,trans_bool in zip(self.star_grid,\
                                           self.model_manager.trans_bool_list)]
            #-- all unique molecules over all stars
            molec_list = list(set([molec 
                                   for model_ids in model_ids_list 
                                   for molec,model_id in model_ids 
                                   if model_ids]))
            #-- all unique modelids for every star separately
            model_id_unique = [list(set([model_id 
                                         for molec,model_id in model_ids])) 
                               for model_ids in model_ids_list]    
            if [modelids for modelids in model_ids_list if modelids] != []:
                appendage += \
                      ['#########################################',\
                       '## Successfully calculated transition model_ids on %s:'\
                       %timestring]
                appendage.extend(['## molecule %s'%molec 
                                  for molec in molec_list])
                for i,(star,model_ids) in enumerate(zip(self.star_grid,\
                                                        model_ids_list)):
                    if model_ids:
                        appendage += ['## For Model %i : cooling id %s'\
                                      %(i+1,star['LAST_GASTRONOOM_MODEL'])] + \
                                     ['#molecule %s #%s' %(molecule,model_id) 
                                      for molecule,model_id in model_ids] + \
                                     ['########']
                for star,model_ids in zip(self.star_grid,model_id_unique):
                    for model_id in model_ids:
                        try:
                            i = 0
                            while True:
                                dummy = DataIO.readFile(\
                                    os.path.join(os.path.expanduser('~'),\
                                        'GASTRoNOoM',self.path_gastronoom,\
                                        'models',model_id,\
                                        os.path.split(self.inputfilename)[1]+\
                                        '_%s_%i'%(model_id,i)))
                                i += 1
                        except IOError:
                            subprocess.call(['cp %s %s'%(self.inputfilename,\
                                        os.path.join(os.path.expanduser('~'),\
                                        'GASTRoNOoM',self.path_gastronoom,\
                                        'models',model_id,\
                                        os.path.split(self.inputfilename)[1]+\
                                        '_%s_%i'%(model_id,i)))],shell=True)
        if self.model_manager.done_mcmax_list:
            model_ids = [star['LAST_MCMAX_MODEL'] 
                         for star,boolean in zip(self.star_grid,\
                                            self.model_manager.done_mcmax_list) 
                         if boolean or self.append_results]
            if model_ids:
                appendage += ['#########################################',\
                        '## MCMax model_ids associated with this grid on %s:'\
                        %timestring]
                appendage += ['#%s'%model_id for model_id in model_ids]
                for model_id in model_ids:
                    try:
                        i = 0
                        while True:
                            dummy = DataIO.readFile(os.path.join(\
                                        os.path.expanduser('~'),'MCMax',\
                                        self.path_mcmax,'models',model_id,\
                                        os.path.split(self.inputfilename)[1]+\
                                        '_%s_%i'%(model_id,i)))
                            i += 1
                    except IOError:
                        subprocess.call(['cp %s %s'%(self.inputfilename,\
                                    os.path.join(os.path.expanduser('~'),\
                                    'MCMax',self.path_mcmax,'models',model_id,\
                                    os.path.split(self.inputfilename)[1]+\
                                    '_%s_%i'%(model_id,i)))],shell=True)
        if appendage: DataIO.writeFile(filename=self.inputfilename,\
                                       input_lines=appendage+['\n'],mode='a')
        
        
        
    def runStatistics(self):
        
        '''
        Run the statistics module.
        
        '''
        
        print '************************************************'
        print '****** Loading Statistics module for %s.'%self.star_name
        print '************************************************'
        if self.statistics and self.pacs <> None:
            self.stats = PeakStats.PeakStats(star_name=self.star_name,\
                                             path_code=self.path_gastronoom,\
                                             path_combocode=self.path_combocode)
            self.stats.setInstrument(instrument='PACS',\
                                     instrument_instance=self.pacs)
            self.stats.setModels(instrument='PACS',star_grid=self.star_grid,\
                                 redo_convolution=self.pacs.redo_convolution)
            self.stats.findRatios(instrument='PACS',\
                                  tolerance=self.stat_tolerance,\
                                  sigma=self.stat_sigma,mode=self.stat_mode)
            self.stats.plotRatioWav(instrument='PACS',\
                                    inputfilename=self.inputfilename,\
                                    tolerance=self.stat_tolerance,\
                                    sigma=self.stat_sigma,mode=self.stat_mode)
        else: 
            self.stats = None    
    
    
    
    def doIceFit(self):
         
         '''
         Run the IceFitter class for this CC session.
         
         '''
         
         if self.do_ice_fit:
              print '** Plotting the 3.1 micron ice feature fit.'
              self.ice = IceFitter.IceFitter(star_grid=self.star_grid,\
                                             spec=self.sed)
              self.ice.prepareModels()
              self.ice.prepareSWS()
              self.ice.show(cfg=self.cfg_ice_fit)
    
    
    
    def printStarInfo(self):
    
        '''
        Print extra Star() info to the shell.
        
        '''
        
        if len(self.star_grid)<20 and self.print_model_info:    
            print '************************************************'
            for star in self.star_grid:
                if star['LAST_MCMAX_MODEL']:
                    print 'Requested MCMax parameters for %s:'\
                          %star['LAST_MCMAX_MODEL']
                    if star.has_key('R_INNER_DUST'): del star['R_INNER_DUST']
                    print '%s = %f R_STAR   (effective inner radius)'%('R_INNER_DUST',star['R_INNER_DUST'])
                    print '%s = %s'%('TEMDUST_FILENAME',star['TEMDUST_FILENAME'])
                    print '%s = %s'%('DUST_TEMPERATURE_FILENAME',star['DUST_TEMPERATURE_FILENAME'])
                    print '%s = %s km/s'%('V_EXP_DUST',star['V_EXP_DUST'])
                    if star.has_key('T_DES_H2O') or star.has_key('T_DES_CH2O')\
                            or star.has_key('T_DES_AH2O'):
                        print ', '.join(['%s = %s K'%(ts,star[ts])
                                         for ts in ['T_DES_H2O','T_DES_CH2O',\
                                                    'T_DES_AH2O']
                                         if star.has_key(ts)])
                        print ', '.join(['%s = %.2f R_STAR = %.2e cm'\
                                         %(rs,star[rs],\
                                           star[rs]*star['R_STAR']*star.r_solar)
                                         for rs in ['R_DES_H2O','R_DES_AH2O',\
                                                    'R_DES_CH2O']
                                         if star.has_key(rs)])
                    print ''
    
                if star['LAST_GASTRONOOM_MODEL']:
                    print 'Requested GASTRoNOoM parameters for %s:'%star['LAST_GASTRONOOM_MODEL']
                    print '%s = %s'%('DENSFILE',star['DENSFILE'])
                    print '%s = %f R_solar = %f AU'%('R_STAR',star['R_STAR'],star['R_STAR']*star.r_solar/star.au)
                    print '%s = %f R_STAR'%('R_OUTER_GAS',star['R_OUTER_EFFECTIVE'])
                    print '%s = %f R_STAR'%('R_INNER_GAS',star['R_INNER_GAS'])
                    print '%s = %f'%('DUST_TO_GAS_INITIAL',star['DUST_TO_GAS_INITIAL'])
                    print '%s = %f'%('DUST_TO_GAS_ITERATED',star['DUST_TO_GAS_ITERATED'])
                    print '%s (semi-empirical) = %f'%('DUST_TO_GAS',star['DUST_TO_GAS'])
                    if star['DUST_TO_GAS_CHANGE_ML_SP']: 
                        print '%s = %s'%('DUST_TO_GAS_CHANGE_ML_SP',star['DUST_TO_GAS_CHANGE_ML_SP'])
                    print '%s = %f'%('[N_H2O/N_H2]',star['F_H2O'])
                    print '%s = %.2f R_STAR = %.2e cm'%('R_OH1612_NETZER',star['R_OH1612_NETZER'],star['R_OH1612_NETZER']*star.r_solar*star['R_STAR'])
                    if star['R_OH1612_AS']: 
                        print '%s = %.2f R_STAR = %.2e cm'%('R_OH1612_OBS',star['R_OH1612'],star['R_OH1612']*star.r_solar*star['R_STAR'])
                        print '-----------------------------------'
                        if star.has_key('R_DES_H2O') or star.has_key('R_DES_CH2O') or star.has_key('R_DES_AH2O'):
                            (nh2o,nh2o_ice,nh2,nh2o_full,nh2_full) = wa.getWaterInfo(star)
                            print 'Ice shell water VAPOUR column density [cm-2]:'
                            print '%.3e'%nh2o
                            print 'Ice shell H2 column density [cm-2]:'
                            print '%.3e'%nh2
                            print 'Total water vapour abundance (ortho + para) wrt H2:'
                            print '%.3e'%(nh2o/nh2)
                            print 'Ice shell water ICE MOLECULAR column density [cm-2]:'
                            print '%.3e'%nh2o_ice
                            print 'Minimum required water vapour abundance N(H2O)/N(H2) for this ice column density:'
                            print '%.3e'%(nh2o_ice/nh2)
                            print 'Ice/vapour fraction in ice shell:'
                            print '%.2f'%(nh2o_ice/nh2o)
                            print 
                            print 'FULL shell water VAPOUR column density [cm-2]:'
                            print '%.3e'%nh2o_full
                            print 'Total water vapour abundance (ortho + para) wrt H2:'
                            print '%.3e'%(nh2o_full/nh2_full)
                            print 
                        else:
                            print 'No water ice present in dust model.'
        print '************************************************'
        


if __name__ == "__main__":
    try:
        inputfilename=sys.argv[1]
    except IndexError:
        raise IOError('Please provide an inputfilename. (syntax in the ' + \
                      'command shell: python ComboCode.py ' + \
                      '/home/robinl/path_combocode/inputComboCode.dat)')
    cc = ComboCode(inputfilename)
    cc.startSession()