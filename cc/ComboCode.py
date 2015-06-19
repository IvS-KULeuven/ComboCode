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

from cc.tools.io import DataIO
from cc.tools.numerical import Gridding
from cc.managers import ModelingManager
from cc.managers import PlottingManager
from cc.managers import Vic
from cc.modeling.objects import Star
from cc.modeling.objects import Transition
from cc.statistics import UnresoStats
from cc.statistics import ResoStats
from cc.statistics import Statistics
from cc.data.instruments import Pacs
from cc.data.instruments import Spire
from cc.data import Sed, Radio
from cc.modeling.tools import ColumnDensity, ContinuumDivision



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
        >>> cc = ComboCode.ComboCode('/home/robinl/ComboCode/aux/inputComboCode.dat')
        
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
        self.setRadio()
        self.setPlotManager()
        self.setVarPars()
        self.createStarGrid()
        self.addRadioData()
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
            self.doContDiv()
            #self.appendResults() 
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
                          ('ln_path_gastronoom',''),('ln_path_mcmax',''),\
                          ('path_gastronoom',''),('path_mcmax',''),\
                          ('print_model_info',1),('stat_chi2','normal'),\
                          ('contdiv_features',[]),('cfg_contdiv',''),\
                          ('show_contdiv',0),('skip_cooling',0),\
                          ('recover_sphinxfiles',0),('stat_print',0),\
                          ('stat_lll_p',None),('stat_method','clipping'),\
                          ('star_name','model'),('opac_path',''),\
                          ('corrflux_path',None)]
        global_pars = dict([(k,self.processed_input.pop(k.upper(),v)) 
                            for k,v in default_global])
        self.__dict__.update(global_pars)
        self.vic = 0
        if not self.gastronoom or not self.mcmax: self.iterations = 1 
        if (not self.path_mcmax and self.mcmax):
            raise IOError('Please define PATH_MCMAX in your inputfile.')
        if (not self.path_gastronoom and self.gastronoom): 
            raise IOError('Please define PATH_GASTRONOOM in your inputfile.')
                
                
                
    def setPacs(self):
        
        '''
        Collect the PACS relevant parameters from the inputfile and set the
        PACS object.
        
        '''
        
        pacs_path = self.processed_input.pop('PACS_PATH','')
        redo_convolution = self.processed_input.pop('PACS_REDO_CONVOLUTION',0)
        searchstring = self.processed_input.pop('PACS_SEARCHSTRING','')
        oversampling = self.processed_input.pop('PACS_OVERSAMPLING','')
        intrinsic = self.processed_input.pop('PACS_INTRINSIC',1)
        linefit = self.processed_input.pop('PACS_LINEFIT','')
        absflux_err = self.processed_input.pop('PACS_ABSFLUX_ERR',0.2)
        if pacs_path:
            self.pacs = Pacs.Pacs(star_name=self.star_name,\
                                  path_combocode=self.path_combocode,\
                                  path=self.path_gastronoom,\
                                  redo_convolution=redo_convolution,\
                                  oversampling=oversampling,\
                                  path_pacs=pacs_path,\
                                  intrinsic=intrinsic,\
                                  path_linefit=linefit,\
                                  absflux_err=absflux_err)
            self.pacs.setData(searchstring=searchstring)
        else:
            self.pacs = None
            
            
            
    def setSpire(self):
        
        '''
        Collect the SPIRE relevant parameters from the inputfile and set the SPIRE
        object.
        
        '''
        
        spire_path = self.processed_input.pop('SPIRE_PATH','')
        searchstring = self.processed_input.pop('SPIRE_SEARCHSTRING','')
        resolution = self.processed_input.pop('SPIRE_RESOLUTION',0)
        intrinsic = self.processed_input.pop('SPIRE_INTRINSIC',1)
        oversampling = self.processed_input.pop('SPIRE_OVERSAMPLING',0)
        linefit = self.processed_input.pop('SPIRE_LINEFIT','')
        absflux_err = self.processed_input.pop('SPIRE_ABSFLUX_ERR',0.2)
        if spire_path:
            self.spire = Spire.Spire(star_name=self.star_name,\
                                     path_combocode=self.path_combocode,\
                                     path=self.path_gastronoom,\
                                     resolution=resolution,\
                                     path_spire=spire_path,\
                                     intrinsic=intrinsic,\
                                     oversampling=oversampling,\
                                     path_linefit=linefit,\
                                     absflux_err=absflux_err)
            self.spire.setData(searchstring=searchstring)
        else:
            self.spire = None



    def setSed(self):
        
        '''
        Collect the SED data and create an Sed() object.
        
        '''
        
        sed_path = self.processed_input.pop('SED_PATH','')
        psuffix = self.processed_input.pop('SED_PHOT_PSUFFIX',\
                                           'Raw/IvS_SEDTool/')
        remove = self.processed_input.pop('SED_PHOT_REMOVE','')
        if not remove: remove = []
        elif type(remove) is types.StringType: remove = [remove]
        else: remove = list(remove)
        
        if sed_path:
            self.sed = Sed.Sed(star_name=self.star_name,\
                               path_combocode=self.path_combocode,\
                               path=sed_path,psuffix=psuffix,remove=remove)
        else: 
            self.sed = None
        
        
    
    def setRadio(self):
        
        '''
        Collect the relevant radio data for the requested star. Only done if 
        the pathname to the data is given.
        
        If a database is not present, it is created. 
        
        If the auto_parse is requested, the data folder will be parsed for new 
        data. 
        
        The data are associated with requested transitions later. Only if also 
        RADIO_AUTOSEARCH is on, these transitions will be automatically added
        to the requested transitions list. 

        '''
        
        self.radio = None
        self.radio_path = self.processed_input.pop('RADIO_PATH','')
        self.radio_autosearch = self.processed_input.pop('RADIO_AUTOSEARCH',0)
        radio_autoparse = self.processed_input.pop('RADIO_AUTOPARSE',0)
        fn = os.path.join(self.radio_path,'radio_data.db')
        if self.radio_path:
            cc_path = os.path.join(self.path_combocode,'Data')
            radio_db = Radio.Radio(path=self.radio_path,cc_path=cc_path,\
                                   auto_parse=radio_autoparse)
            if radio_db.has_key(self.star_name):
                self.radio = radio_db[self.star_name]            
            
            
        
    def addRadioData(self):
        
        '''
        Add radio data to Transition() objects in all Star() objects.
        
        Only done if RADIO_PATH is given and if a file named radio_data.db is 
        present in the given folder.
        
        If the radio_autosearch flag is on, transitions are automatically 
        generated based on the available data. Note that in this case, N_QUAD
        from Star() is taken. 
        
        '''
        
        if self.radio:
            #-- Get the transition definitions (are in the correct format 
            #   automatically, due to the methods in Radio.py) and make sure 
            #   they are all unique.
            radio_trans = sorted(['%s 100'%tr.replace('TRANSITION=','',1) 
                                  for tr in self.radio.keys()])
            radio_trans = DataIO.checkEntryInfo(radio_trans,12,'TRANSITION')
            for star in self.star_grid:
                molecules = [m.molecule for m in star['GAS_LIST']]
                if self.radio_autosearch:
                    n_quad = star['N_QUAD']
                    add_trans = [tr[:-1] + [n_quad] 
                                 for tr in radio_trans
                                 if tr[0] in molecules]
                    if star.has_key('TRANSITION'):
                        star['TRANSITION'] = list(star['TRANSITION'])
                        star['TRANSITION'].extend(add_trans)
                    else:
                        star['TRANSITION'] = add_trans
                for trans in star['GAS_LINES']:
                    if trans:
                        trstr = trans.getInputString(include_nquad=0)
                        if trstr in self.radio.keys():
                            trans.addDatafile(self.radio[trstr],\
                                              path=self.radio_path)

        

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
                                     convert_ints=1,multi_keys=multi_keys)
        #-- keywords in multi_keys require different method
        self.processed_input = dict()
        self.multiplicative_grid = dict()
        self.additive_grid = dict()
        molecules = input_dict.pop('MOLECULE',[])
        transitions = input_dict.pop('TRANSITION',[])
        r_points_mass_loss = input_dict.pop('R_POINTS_MASS_LOSS',[])      
        for k,v in input_dict.items():
            #-- Fortran input is not case sensitive. Note: use the dict
            #   value v, not input_dict[k] because of this transformation.
            k = k.upper()
            #-- Determine delimiter        
            try:
                if v.find('&') != -1: delimiter = '&'                         
                elif v.find(';') != -1: delimiter = ';'
                elif v.find(',') != -1: delimiter = ','
                elif v.find(':') != -1: delimiter = ':'
                #-- * while no ; or , or : means multiple values for ONE model
                #   Only * star for multiplicative grid makes no sense (just 
                #   give the value without delimiter)
                elif v.find('*') != -1: delimiter = '&'
                else: delimiter = ' '                                                        
            except AttributeError: 
                #-- v is already a float, so can't use .find on it => no grids
                #-- no need to check the rest, continue on with the next k/v pair
                self.processed_input[k] = v
                continue
            #-- Expanding '*' entries: Assumes the value is first, the count 
            #   second. Can't be made flexible, because in some cases the value
            #   cannot be discerned from the count (because both are low-value
            #   integers)
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
            molecules = DataIO.checkEntryInfo(molecules,20,'MOLECULE')
            if type(molecules) is types.ListType:
                self.multiplicative_grid['MOLECULE'] = molecules
            else:
                self.processed_input['MOLECULE'] = molecules
        if r_points_mass_loss: 
            r_points_mass_loss = DataIO.checkEntryInfo(r_points_mass_loss,4,\
                                                       'R_POINTS_MASS_LOSS')
            if type(r_points_mass_loss) is types.ListType:
                self.additive_grid['R_POINTS_MASS_LOSS'] = r_points_mass_loss
            else:
                self.processed_input['R_POINTS_MASS_LOSS'] = r_points_mass_loss
        if transitions: 
            transitions = DataIO.checkEntryInfo(transitions,12,'TRANSITION')
            if type(transitions) is types.ListType:
                self.multiplicative_grid['TRANSITION'] = transitions
            else:
                self.processed_input['TRANSITION'] = transitions
            
            
    
    def getStars(self):
         
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
        for star in self.star_grid:
            star.normalizeDustAbundances()
        if self.processed_input.has_key('LAST_MCMAX_MODEL'):
            del self.processed_input['LAST_MCMAX_MODEL']
        if self.processed_input.has_key('LAST_GASTRONOOM_MODEL'):
            del self.processed_input['LAST_GASTRONOOM_MODEL']
        #-- Dust abundance is deleted as it is not changed during the session. 
        #   Hence, it does not need to be re-updated after mutable input is 
        #   removed. The keys need to be deleted to avoid inconsistencies in the
        #   star.normalizeDustAbundances() method. Some A_*** values may not be
        #   variable, while others are. Yet if any of them are variable, and 
        #   have to be rescaled, then in principle they are ALL variable. This
        #   can create a mess, and therefore it is safer to just remove them 
        #   from the input dictionary.
        for akey in [k for k in self.processed_input.keys() if k[0:2] == 'A_']:
            del self.processed_input[akey]
            
        
    
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
        
        if self.vic and self.gastronoom and self.sphinx : 
            self.vic_manager = Vic.Vic(path=self.path_gastronoom,\
                                       path_combocode=self.path_combocode,\
                                       account=self.vic_account,\
                                       time_per_sphinx=self.vic_time_per_sphinx,\
                                       credits_acc=self.vic_credits,\
                                       recover_sphinxfiles=self.recover_sphinxfiles)
            if self.update_spec: 
                self.vic_manager.updateLineSpec()
        else: 
            self.vic_manager = None
        
        
    def setModelManager(self):
        
        '''
        Set up the model manager.
        
        '''    
        
        self.model_manager = ModelingManager.ModelingManager(\
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
                                       path_mcmax=self.path_mcmax,\
                                       skip_cooling=self.skip_cooling,\
                                       recover_sphinxfiles=self.recover_sphinxfiles,\
                                       opac_path=self.opac_path)
    
    
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
                                         corrflux_path = self.corrflux_path,\
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
                #-- mline_done is True if in previous model an mline calculation 
                #-- was done: Only then do a progress check, because a lot of time 
                #-- has passed, but then a wait time is used to make sure the newly
                #-- queued sphinx models after the mline model are properly queued.
                if self.vic_manager \
                        and self.vic_manager.getQueue() \
                        and self.model_manager.mline_done_list[-1]:    
                    print '***********************************'
                    print '** Current VIC queue:'
                    print self.vic_manager.getQueue()
                    self.vic_manager.checkProgress(wait_qstat=1)



    def finalizeVic(self):
        
        '''
        At the end of a modeling session, allow Vic to be finalized and clean up.
        
        '''
        
        if self.vic_manager <> None:
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
        if self.model_manager.mcmax_done_list:
            model_ids = [star['LAST_MCMAX_MODEL'] 
                         for star,boolean in zip(self.star_grid,\
                                            self.model_manager.mcmax_done_list) 
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
        
        self.pacsstats = None
        self.spirestats = None
        self.resostats = None
        if self.statistics and self.pacs <> None:
            print '************************************************'
            print '**** Doing PACS statistics for %s.'%self.star_name
            print '************************************************'
            self.pacsstats = UnresoStats.UnresoStats(star_name=self.star_name,\
                                            path_code=self.path_gastronoom,\
                                            path_combocode=self.path_combocode)
            self.pacsstats.setInstrument(instrument_name='PACS',\
                                         instrument_instance=self.pacs,\
                                         stat_method=self.stat_method)
            self.pacsstats.setModels(star_grid=self.star_grid)
            self.pacsstats.setRatios(chi2_type=self.stat_chi2)
            self.pacsstats.plotRatioWav(inputfilename=self.inputfilename)
        if self.statistics and self.spire <> None:
            print '************************************************'
            print '**** Doing SPIRE statistics for %s.'%self.star_name
            print '************************************************'
            self.spirestats = UnresoStats.UnresoStats(star_name=self.star_name,\
                                            path_code=self.path_gastronoom,\
                                            path_combocode=self.path_combocode)
            self.spirestats.setInstrument(instrument_name='SPIRE',\
                                          instrument_instance=self.spire,\
                                          stat_method=self.stat_method)
            self.spirestats.setModels(star_grid=self.star_grid)
            self.spirestats.setRatios(chi2_type=self.stat_chi2)
            self.spirestats.plotRatioWav(inputfilename=self.inputfilename)
        if self.statistics:
            trans_sel = Transition.extractTransFromStars(self.star_grid,\
                                                         dtype='resolved')
            if not trans_sel:
                return
            print '************************************************'
            print '**** Doing statistics for spectrally resolved lines in %s.'\
                  %self.star_name
            print '**** Use cc_session.resostats for more interactive tools.'
            print '************************************************'
            self.resostats = ResoStats.ResoStats(star_name=self.star_name,\
                                           path_code=self.path_gastronoom,\
                                           path_combocode=self.path_combocode,\
                                           lll_p=self.stat_lll_p)
            self.resostats.setInstrument(trans_sel)
            self.resostats.setModels(star_grid=self.star_grid)
            self.resostats.setIntensities()
            if self.stat_print:
                self.resostats.printStats()
            #bfms = self.resostats.selectBestFitModels(mode='int')
            #self.plot_manager.plotTransitions(star_grid=bfms,fn_suffix='BFM',force=1)
            
    
    
    def doContDiv(self):
         
        '''
        Run the Continuum Division class for this CC session.
         
        '''
        
        if type(self.contdiv_features) is types.StringType:
            self.contdiv_features = [self.contdiv_features]
        for k in self.contdiv_features: 
            features = ['MGS','H2O3.1']
            if k in features:
                print '** Plotting continuum division for the %s feature.'%k
                all_franges = [[20.0,23.,46.,48.5],\
                               [2.6,2.85,3.3,3.7]]
                all_funcs = ['linear',\
                             'power']
                franges = all_franges[features.index(k)]
                func = all_funcs[features.index(k)]
                self.contdiv = ContinuumDivision.ContinuumDivision(\
                                        star_grid=self.star_grid,\
                                        spec=[self.sed],franges=franges,\
                                        plot=self.show_contdiv,func=func,\
                                        path_combocode=self.path_combocode,\
                                        cfg=self.cfg_contdiv)
                self.contdiv.prepareModels()
                self.contdiv.prepareData()
                self.contdiv.show()
                for key,val in self.contdiv.eq_width.items():
                    if key[0:3] == 'sws':
                        key = self.sed.star_name_plots
                        print 'Equivalent width for %s in %s: %f'\
                              %(k,self.sed.star_name_plots,val)
                for star in self.star_grid:
                    model_id = star['LAST_MCMAX_MODEL']
                    print 'Equivalent width for %s in %s: %f'\
                          %(k,model_id,self.contdiv.eq_width[model_id])
    
    
    
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

                    if not int(star['MRN_DUST']):
                        cdcalc = ColumnDensity.ColumnDensity(star)
                        dlist = [sp 
                                 for sp in star.getDustList() 
                                 if 'H2O' not in sp]
                        print 'Dust FULL column densities (if available):'
                        for species in dlist:
                            print 'C_%s = %.3e g/cm^-2'\
                               %(species,cdcalc.dustFullColDens(species))
                        print 'Dust FULL number column densities (if available):'
                        for species in dlist:
                            print 'N_%s = %.3e cm^-2'\
                               %(species,cdcalc.dustFullNumberColDens(species))
                        print 'Dust associated molecular abundances (if available):'
                        print 'Note: Does not use above column densities. See ColumnDensity.py.'
                        for species in dlist:
                            print 'A_%s/A_H2 = %.3e'\
                                %(species,cdcalc.dustMolecAbun(species))
                    if star.has_key('T_DES_H2O') or star.has_key('T_DES_CH2O')\
                            or star.has_key('T_DES_AH2O') \
                            and not int(star['MRN_DUST']):
                        print ', '.join(['%s = %s K'%(ts,star[ts])
                                         for ts in ['T_DES_H2O','T_DES_CH2O',\
                                                    'T_DES_AH2O']
                                         if star.has_key(ts)])
                        print ', '.join(['%s = %.2f R_STAR = %.2e cm'\
                                         %(rs,star[rs],\
                                           star[rs]*star['R_STAR']*star.Rsun)
                                         for rs in ['R_DES_H2O','R_DES_AH2O',\
                                                    'R_DES_CH2O']
                                         if star.has_key(rs)])
                        cdh2o1 = cdcalc.dustFullColDens('CH2O')
                        cdh2o2 = cdcalc.dustFullColDens('AH2O')
                        print 'The FULL water ice column density:'
                        print 'C_H2O_ice = %.3e g cm-2'%(cdh2o1+cdh2o2)
                        ncdh2o1 = cdcalc.dustFullNumberColDens('CH2O')
                        ncdh2o2 = cdcalc.dustFullNumberColDens('AH2O')
                        print 'The FULL water ice column number density:'
                        print 'N_H2O_ice = %.3e cm-2'%(ncdh2o1+ncdh2o2)
                        nh2o1 = cdcalc.dustMolecAbun('CH2O')
                        nh2o2 = cdcalc.dustMolecAbun('AH2O')
                        print 'The associated molecular abundance of water ice:'
                        print 'A_H2O_ice/A_H2 = %.3e'%(nh2o1+nh2o2)
                    print ''
                if star['LAST_GASTRONOOM_MODEL']:
                    print 'Requested GASTRoNOoM parameters for %s:'%star['LAST_GASTRONOOM_MODEL']
                    print '%s = %s'%('DENSFILE',star['DENSFILE'])
                    print '%s = %f Rsun = %f AU'%('R_STAR',star['R_STAR'],star['R_STAR']*star.Rsun/star.au)
                    print '%s = %f R_STAR'%('R_OUTER_GAS',star['R_OUTER_EFFECTIVE'])
                    print '%s = %f R_STAR'%('R_INNER_GAS',star['R_INNER_GAS'])
                    print '%s = %f'%('DUST_TO_GAS_INITIAL',star['DUST_TO_GAS_INITIAL'])
                    print '%s = %f'%('DUST_TO_GAS_ITERATED',star['DUST_TO_GAS_ITERATED'])
                    print '%s (semi-empirical) = %f'%('DUST_TO_GAS',star['DUST_TO_GAS'])
                    if star['DUST_TO_GAS_CHANGE_ML_SP']: 
                        print '%s = %s'%('DUST_TO_GAS_CHANGE_ML_SP',star['DUST_TO_GAS_CHANGE_ML_SP'])
                    print '%s = %f'%('[N_H2O/N_H2]',star['F_H2O'])
                    print '%s = %.2f R_STAR = %.2e cm'%('R_OH1612_NETZER',star['R_OH1612_NETZER'],star['R_OH1612_NETZER']*star.Rsun*star['R_STAR'])
                    if star['R_OH1612_AS']: 
                        print '%s = %.2f R_STAR = %.2e cm'%('R_OH1612_OBS',star['R_OH1612'],star['R_OH1612']*star.Rsun*star['R_STAR'])
                    print '-----------------------------------'
                    #if star.has_key('R_DES_H2O') or star.has_key('R_DES_CH2O') or star.has_key('R_DES_AH2O'):
                        #(nh2o,nh2o_ice,nh2,nh2o_full,nh2_full) = wa.getWaterInfo(star)
                        #print 'Ice shell water VAPOUR column density [cm-2]:'
                        #print '%.3e'%nh2o
                        #print 'Ice shell H2 column density [cm-2]:'
                        #print '%.3e'%nh2
                        #print 'Total water vapour abundance (ortho + para) wrt H2:'
                        #print '%.3e'%(nh2o/nh2)
                        #print 'Ice shell water ICE MOLECULAR column density [cm-2]:'
                        #print '%.3e'%nh2o_ice
                        #print 'Minimum required water vapour abundance N(H2O)/N(H2) for this ice column density:'
                        #print '%.3e'%(nh2o_ice/nh2)
                        #print 'Ice/vapour fraction in ice shell:'
                        #print '%.2f'%(nh2o_ice/nh2o)
                        #print 
                        #print 'FULL shell water VAPOUR column density [cm-2]:'
                        #print '%.3e'%nh2o_full
                        #print 'Total water vapour abundance (ortho + para) wrt H2:'
                        #print '%.3e'%(nh2o_full/nh2_full)
                        #print 
                        #if not int(star['MRN_DUST'):

                    #else:
                        #print 'No water ice present in dust model.'
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