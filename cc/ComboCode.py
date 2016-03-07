# -*- coding: utf-8 -*-

"""
Main code for running GASTRoNOoM (L.Decin) and MCMax (M.Min).

Author: R. Lombaert

"""

import sys
import os
import subprocess
import time

import cc.path
from cc.tools.io import DataIO
from cc.tools.numerical import Gridding
from cc.managers.ModelingManager import ModelingManager as MM
from cc.managers.PlottingManager import PlottingManager as PM
from cc.managers import Vic
from cc.modeling.objects import Star, Transition
from cc.statistics import UnresoStats, ResoStats, SedStats
from cc.data.instruments import Pacs, Spire
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
        >>> cc = ComboCode.ComboCode('/home/robinl/ComboCode/input/inputComboCode.dat')

        @param inputfilename: The name of the inputfile.
        @type inputfilename: string

        '''

        self.inputfilename = inputfilename
        self.readInput()
        self.setGlobalPars()
        self.setOutputFolders()
        self.setPacs()
        self.setSpire()
        self.setSed()
        self.setRadio()
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
                          ('path_gastronoom',''),('path_mcmax',''),\
                          ('print_model_info',1),('stat_chi2','diff'),\
                          ('contdiv_features',[]),('cfg_contdiv',''),\
                          ('show_contdiv',0),('skip_cooling',0),\
                          ('recover_sphinxfiles',0),('stat_print',0),\
                          ('stat_lll_p',None),('stat_method','clipping'),\
                          ('star_name','model'),\
                          ('stat_lll_partial',0),('stat_lll_vcut',0.0)]
        global_pars = dict([(k,self.processed_input.pop(k.upper(),v))
                            for k,v in default_global])
        self.__dict__.update(global_pars)
        self.__setStarName()
        self.vic = 0
        if not self.gastronoom or not self.mcmax: self.iterations = 1
        if (not self.path_mcmax and self.mcmax):
            raise IOError('Please define PATH_MCMAX in your inputfile.')
        if (not self.path_gastronoom and self.gastronoom):
            raise IOError('Please define PATH_GASTRONOOM in your inputfile.')



    def __setStarName(self):

        '''
        Set star_name for the ComboCode object as a tuple.

        Typically this is only one name for a standard modelling session, but
        can be made multiple names as well for a statistical study.

        The ComboCode object keeps track of all the data in dicts.

        '''
        
        if self.multiplicative_grid.has_key('STAR_NAME') \
                or self.additive_grid.has_key('STAR_NAME'):
            raise IOError('STAR_NAME incorrectly defined. Use & for grids.')
        
        if isinstance(self.star_name,str):
            self.star_name = (self.star_name,)
        


    def setPacs(self):

        '''
        Collect the PACS relevant parameters from the inputfile and set the
        PACS object.

        '''
        
        pacs = self.processed_input.pop('PACS',0)
        redo_convolution = self.processed_input.pop('PACS_REDO_CONVOLUTION',0)
        searchstring = self.processed_input.pop('PACS_SEARCHSTRING','')
        oversampling = self.processed_input.pop('PACS_OVERSAMPLING','')
        intrinsic = self.processed_input.pop('PACS_INTRINSIC',1)
        linefit = self.processed_input.pop('PACS_LINEFIT','')
        
        #-- If PACS is not requested, put self.pacs to None. Still popping the
        #   PACS specific keywords to avoid clutter in the Star() objects.
        #   In case of a pure model, no data are available anyway.
        self.pacs = dict()
        for sn in self.star_name: 
            if not pacs or sn == 'model':
                self.pacs[sn] = None
                continue
            self.pacs[sn] = Pacs.Pacs(star_name=sn,\
                                      path=self.path_gastronoom,\
                                      redo_convolution=redo_convolution,\
                                      oversampling=oversampling,\
                                      intrinsic=intrinsic,\
                                      path_linefit=linefit)
            self.pacs[sn].setData(searchstring=searchstring)



    def setSpire(self):

        '''
        Collect the SPIRE relevant parameters from the inputfile and set the SPIRE
        object.

        '''
        
        spire = self.processed_input.pop('SPIRE',0)
        searchstring = self.processed_input.pop('SPIRE_SEARCHSTRING','')
        resolution = self.processed_input.pop('SPIRE_RESOLUTION',0)
        intrinsic = self.processed_input.pop('SPIRE_INTRINSIC',1)
        oversampling = self.processed_input.pop('SPIRE_OVERSAMPLING',0)
        linefit = self.processed_input.pop('SPIRE_LINEFIT','')
        
        #-- If SPIRE is not requested, put self.spire to None. Still popping the
        #   SPIRE specific keywords to avoid clutter in the Star() objects.
        #   In case of a pure model, no data are available anyway.
        self.spire = dict()
        for sn in self.star_name: 
            if not spire or sn == 'model':
                self.spire[sn] = None
                continue
            self.spire[sn] = Spire.Spire(star_name=sn,\
                                         path=self.path_gastronoom,\
                                         resolution=resolution,\
                                         intrinsic=intrinsic,\
                                         oversampling=oversampling,\
                                         path_linefit=linefit)
            self.spire[sn].setData(searchstring=searchstring)



    def setSed(self):

        '''
        Collect the SED data and create an Sed() object.

        '''
        
        sed = self.processed_input.pop('SED',0)
        remove = self.processed_input.pop('SED_PHOT_REMOVE','')
        
        if not remove: remove = []
        elif isinstance(remove,str): remove = [remove]
        else: remove = list(remove)
        
        #-- If SED is not requested, put self.spire to None. Still popping the
        #   SED specific keywords to avoid clutter in the Star() objects.
        #   In case of a pure model, no data are available anyway.
        self.sed = dict()
        for sn in self.star_name: 
            if not sed or sn == 'model':
                self.sed[sn] = None
                continue
            self.sed[sn] = Sed.Sed(star_name=sn,remove=remove)



    def setRadio(self):

        '''
        Collect the relevant radio data for the requested star. Only done if
        the pathname to the data is given.

        If a database is not present, it is created.

        If the radio_autosearch flag is on, transitions are automatically
        generated based on the available data. Note that in this case, N_QUAD
        from Star() is taken.
        
        '''
        
        radio = self.processed_input.pop('RADIO',0)
        radio_autosearch = self.processed_input.pop('RADIO_AUTOSEARCH',0)

        if not radio: 
            self.radio = {sn: None for sn in self.star_name}
            return
            
        #-- If RADIO is not requested, put self.radio to None. Still popping the
        #   RADIO specific keywords to avoid clutter in the Star() objects.
        #   In case of a pure model, the radio db doesn't have an entry anyway.
        radio_db = Radio.Radio()
        
        self.radio = {sn: radio_db.get(sn, None) for sn in self.star_name}
        self.radio_trans = {sn: None for sn in self.star_name}
        
        if not radio_autosearch:
            return

        #-- Get the transition definitions (are in the correct format
        #   automatically, due to the methods in Radio.py). Check entry info
        #   is still ran, eg to get rid of bad 0 offset type sets.
        #-- 11 entries, n_quad is added in the Star() object (to allow for
        #   proper n_quad gridding). N_quad specification through manual
        #   TRANSITION definition is still possible, but then overrides the 
        #   default N_QUAD value.
        radio_trans = sorted(set([tr.replace('TRANSITION=','',1)
                                  for sn in self.star_name
                                  if self.radio[sn]
                                  for tr in self.radio[sn].keys()]))
        radio_trans = DataIO.checkEntryInfo(radio_trans,11,'TRANSITION')
        
        #-- Doubles of transitions not possible from db. radio_trans always
        #   a tuple.
        #-- In case doubles after merge: they will be filtered out by Star()
        if self.additive_grid.has_key('TRANSITION'):
            otrl = self.additive_grid['TRANSITION'] 
            ntrl = [tl + radio_trans for tl in otrl]
            self.additive_grid['TRANSITION'] = ntrl
        elif self.processed_input.has_key('TRANSITION'):
            self.processed_input['TRANSITION'] += radio_trans
        else:
            self.processed_input['TRANSITION'] = radio_trans



    def addRadioData(self,sn):

        '''
        Add radio data to Transition() objects in all Star() objects.

        Only done if RADIO_PATH is given and if a file named radio_data.db is
        present in the given folder.

        This method can be called multiple times per session for different stars
        in which case the data in the transition objects will be replaced.
        
        self.radio_trans maintains a list of "sample transitions" that contain 
        data. These function as data blueprints for transitions in the star_grid
        
        @param sn: The star name for which to add data to the transitions
        @type sn: str
        
        '''
    
        #-- If not data available, reset all transitions so no data linger for a
        #   different star.
        if not self.radio[sn]:
            for star in self.star_grid:
                for tr in star['GAS_LINES']: 
                    tr.resetData()
            return
        
        #-- Extracted transitions (copies) are saved separately, so we only read
        #   data once at most. Hence we can continuously swap between star as we
        #   progress.
        if not self.radio_trans[sn]:        
            #-- First extract all transitions and add/read data for them.
            trans_sel = Transition.extractTransFromStars(self.star_grid,\
                                                         dtype='resolved')

            #-- Note that this is a list of copied transitions. 
            for trans in trans_sel: 
                trstr = trans.getInputString(include_nquad=0)
                if trstr in self.radio[sn].keys():
                    #-- No need to replace. Trans extraction already resets data
                    trans.addDatafile(self.radio[sn][trstr])
                    trans.readData()
                    
            #-- self.radio_trans was defined in setRadio.
            self.radio_trans[sn] = trans_sel

        #-- We can now add these data to all Star() objects. 
        for star in self.star_grid:
            for trans in self.radio_trans[sn]:
                mtrans = star.getTransition(trans)
                if mtrans: mtrans.setData(trans,replace=1)



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
            if isinstance(molecules,list):
                self.additive_grid['MOLECULE'] = molecules
            else:
                self.processed_input['MOLECULE'] = molecules
        if r_points_mass_loss:
            r_points_mass_loss = DataIO.checkEntryInfo(r_points_mass_loss,4,\
                                                       'R_POINTS_MASS_LOSS')
            if isinstance(r_points_mass_loss,list):
                self.additive_grid['R_POINTS_MASS_LOSS'] = r_points_mass_loss
            else:
                self.processed_input['R_POINTS_MASS_LOSS'] = r_points_mass_loss
        if transitions:
            nk = 12
            
            #-- Check if N_QUAD is a grid parameter: remove any manual N_QUAD
            #   definitions to allow gridding to work. N_QUAD is taken from 
            #   Star() by default in this case.
            gkeys = self.multiplicative_grid.keys() + self.additive_grid.keys()
            if 'N_QUAD' in gkeys: 
                transitions = [' '.join(tr.split()[:-1]) for tr in transitions]
                nk = 11
            
            #-- If N_QUAD is a gridding key, transitions will always be a tuple.
            transitions = DataIO.checkEntryInfo(transitions,nk,'TRANSITION')
            if isinstance(transitions,list):
                self.additive_grid['TRANSITION'] = transitions
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
                                         path_gastronoom=self.path_gastronoom,\
                                         path_mcmax=self.path_mcmax,\
                                         extra_input=d)
                                  for d in additive_dicts]
        else:
            self.star_grid = [base_star]
        for key,grid in self.multiplicative_grid.items():
            self.star_grid = [Star.Star(path_gastronoom=self.path_gastronoom,\
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



    def setOutputFolders(self):

        '''
        Set the output folders.

        If the folders do not already exist, they are created.

        The locations are saved in cc.path for later use, but this is generally
        only done inside a ComboCode session. Each module sets these themselves

        '''

        cc.path.gout = os.path.join(cc.path.gastronoom,self.path_gastronoom)
        cc.path.mout = os.path.join(cc.path.mcmax,self.path_mcmax)
        DataIO.testFolderExistence(cc.path.gout)
        DataIO.testFolderExistence(cc.path.mout)


    def setVicManager(self):

        '''
        Set up the VIC manager.

        '''

        if self.vic and self.gastronoom and self.sphinx :
            self.vic_manager = Vic.Vic(path=self.path_gastronoom,\
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

        self.model_manager = MM(iterations=self.iterations,\
                                processed_input=self.processed_input,\
                                var_pars=self.var_pars,\
                                path_gastronoom=self.path_gastronoom,\
                                mcmax=self.mcmax,\
                                gastronoom=self.gastronoom,\
                                sphinx=self.sphinx,\
                                iterative=self.plot_iterative,\
                                num_model_sessions=len(self.star_grid),\
                                vic_manager=self.vic_manager,\
                                replace_db_entry=self.replace_db_entry,\
                                path_mcmax=self.path_mcmax,\
                                skip_cooling=self.skip_cooling,\
                                recover_sphinxfiles=self.recover_sphinxfiles)


    def setPlotManager(self):

        '''
        Set up the plot manager(s) for each star name.

        '''

        plot_pars = dict([(k,self.processed_input.pop(k))
                          for k,v in self.processed_input.items()
                          if k[0:5] == 'PLOT_' or k[0:4] == 'CFG_'])
        fn_add_star = plot_pars.pop('PLOT_FN_ADD_STAR',1)
        self.plot_manager = {sn: PM(star_name=sn,\
                                    gastronoom=self.gastronoom,\
                                    mcmax=self.mcmax,\
                                    path_gastronoom=self.path_gastronoom,\
                                    path_mcmax=self.path_mcmax,\
                                    inputfilename=self.inputfilename,\
                                    pacs=self.pacs[sn],\
                                    spire=self.spire[sn],\
                                    sed=self.sed[sn],\
                                    fn_add_star=fn_add_star,\
                                    plot_pars=plot_pars)
                             for sn in self.star_name}



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

        print '************************************************'
        print '****** Plotting final results.'
        print '************************************************'
        for sn in self.star_name:
            #-- Make sure resolved data files are correctly included. 
            self.addRadioData(sn)
            print '** Plots for %s:'%sn 
            if self.plot_iterative:
                print '** Plotting results for each iterative step.'
                print '************'
                #- star_grid_old remembers all old models if iterative is on,
                #- meaning that every list in star_grid_old consists of Star models
                #- associated with one level of iteration. Following line plots all
                #- iterative steps for one star immutable parameter set
                pm = self.plot_manager[sn]
                for i in range(len(self.star_grid)):
                    pm.startPlotting(self.model_manager.star_grid_old[i],i+1)
            self.plot_manager[sn].startPlotting(self.star_grid)



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
                                    os.path.join(cc.path.gout,'models',model_id,\
                                        os.path.split(self.inputfilename)[1]+\
                                        '_%s_%i'%(model_id,i)))
                                i += 1
                        except IOError:
                            subprocess.call(['cp %s %s'%(self.inputfilename,\
                                        os.path.join(cc.path.gout,\
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
                                        cc.path.mout,'models',model_id,\
                                        os.path.split(self.inputfilename)[1]+\
                                        '_%s_%i'%(model_id,i)))
                            i += 1
                    except IOError:
                        subprocess.call(['cp %s %s'%(self.inputfilename,os.path.join(\
                                    cc.path.mout,'models',model_id,\
                                    os.path.split(self.inputfilename)[1]+\
                                    '_%s_%i'%(model_id,i)))],shell=True)
        if appendage: DataIO.writeFile(filename=self.inputfilename,\
                                       input_lines=appendage+['\n'],mode='a')



    def runStatistics(self):

        '''
        Run the statistics module.

        '''

        self.pacsstats = dict()
        self.spirestats = dict()
        self.unresostats = []
        self.resostats = dict()
        self.sedstats = dict()
        
        #-- Data are handled by these objects themselves for unresolved lines.
        for sn in self.star_name:
            if self.statistics and self.pacs[sn] <> None:
                ss = UnresoStats.UnresoStats(star_name=sn,\
                                             path_code=self.path_gastronoom)
                ss.setInstrument(instrument_name='PACS',\
                                 instrument_instance=self.pacs[sn],\
                                 stat_method=self.stat_method)
                self.pacsstats[sn] = ss
                self.unresostats.append(ss)
            if self.statistics and self.spire[sn] <> None:
                ss = UnresoStats.UnresoStats(star_name=sn,\
                                             path_code=self.path_gastronoom)
                ss.setInstrument(instrument_name='SPIRE',\
                                 instrument_instance=self.spire[sn],\
                                 stat_method=self.stat_method)
                self.spirestats[sn] = ss
                self.unresostats.append(ss)
        
        for ss in self.unresostats:
            instr = ss.instrument.instrument.upper()
            print '************************************************'
            print '** %s statistics for %s.'%(instr,ss.star_name)
            print '************************************************'
            ss.setModels(star_grid=self.star_grid)
            ss.setLineStrengths()
            ss.calcChiSquared(chi2_method=self.stat_chi2)
            ss.plotRatioWav(inputfilename=self.inputfilename)
        
        for sn in self.star_name:
            if self.statistics and self.sed[sn] <> None:
                print '************************************************'
                print '** SED statistics for %s.'%sn
                print '************************************************'
                ss = SedStats.SedStats(star_name=sn,path_code=self.path_mcmax)
                ss.setInstrument(sed=self.sed[sn])
                ss.setModels(star_grid=self.star_grid)
                ss.setModelPhotometry()
                ss.calcChi2(chi2_method=self.stat_chi2)
                self.sedstats[sn] = ss            
            
            if self.statistics and self.radio[sn]:
                #-- Make sure all resolved data properties are set. 
                self.addRadioData(sn)
                if not self.radio_trans[sn]:
                    return
                print '************************************************'
                print '** Statistics for spectrally resolved lines in %s.'%sn
                print '** Use cc_session.resostats[sn] for interactive tools.'
                print '************************************************'
                ss = ResoStats.ResoStats(star_name=sn,\
                                         path_code=self.path_gastronoom,\
                                         lll_p=self.stat_lll_p)
                ss.setInstrument(self.radio_trans[sn])
                ss.setModels(star_grid=self.star_grid)
                ss.setIntensities(partial=self.stat_lll_partial,\
                                  vcut=self.stat_lll_vcut)
                if self.stat_print: ss.printStats()
                self.resostats[sn] = ss
                #bfms = self.resostats.selectBestFitModels(mode='int')
                #self.plot_manager.plotTransitions(star_grid=bfms,fn_suffix='BFM',force=1)



    def doContDiv(self):

        '''
        Run the Continuum Division class for this CC session.

        '''

        if isinstance(self.contdiv_features,str):
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
                      '/home/robinl/inputComboCode.dat)')
    c1m = ComboCode(inputfilename)
    c1m.startSession()