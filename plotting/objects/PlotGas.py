# -*- coding: utf-8 -*-

"""
A plotting environment for gas transitions and all that is associated with that.

Author: R. Lombaert

"""

import os
from scipy import array
import operator
import subprocess
from scipy.interpolate import interp1d

from cc.plotting.PlottingSession import PlottingSession
from cc.tools.io import DataIO
from cc.modeling.objects import Transition
from cc.plotting import Plotting2
from cc.tools.io import LineList
from cc.data.instruments import Pacs
from cc.modeling.objects import Star



class PlotGas(PlottingSession):
    
    """ 
    Class for plotting gas lines and their information.
    
    """    
    
    def __init__(self,star_name,path_gastronoom='codeSep2010',
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),\
                 inputfilename=None,pacs=None,spire=None):
        
        """ 
        Initializing an instance of PlotGas.
        
        @param star_name: name of the star from Star.dat, use default only 
                          when never using any star model specific things 
        @type star_name: string
        
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        @keyword path_gastronoom: Output modeling folder in MCMax home folder
        
                                  (default: 'codeSep2010')
        @type path_gastronoom: string
        @keyword inputfilename: name of inputfile that is also copied to the 
                                output folder of the plots, 
                                if None nothing is copied
                                
                                (default: None)
        @type inputfilename: string
        @keyword pacs: a Pacs object needed for plotting PACS data. None of no 
                       PACS data involved.
                            
                       (default: None)
        @type pacs: Pacs()
        @keyword spire: a Spire object needed for plotting SPIRE data. None of no 
                        SPIRE data involved.
                            
                        (default: None)
        @type spire: Spire()
        
        """
        
        super(PlotGas, self).__init__(star_name=star_name,\
                                      path_combocode=path_combocode,\
                                      path=path_gastronoom,\
                                      code='GASTRoNOoM',\
                                      inputfilename=inputfilename)
        self.star_name_gastronoom = \
                DataIO.getInputData(path=os.path.join(self.path_combocode,\
                                                    'Data'),\
                                  keyword='STAR_NAME_GASTRONOOM')\
                                 [self.star_index]
        self.pacs = pacs
        self.spire = spire
        self.sphinx_flux_list = []

        
        
    def makeExec(self,star,i,sphinx_id):
        
        """
        Making an execGASTRoNOoM file for model calculation.
        
        @param star: The parameter set
        @type star: Star()
        @param i: Index of Star() in the star_grid passed to the Plot## method.
                  Is added to the filename of the plot.
        @type i: int
        @param sphinx_id: The sphinx model id used for this exec file
        @type sphinx_id: string
        
        """
        
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                self.code,self.path,'stars',self.star_name_gastronoom))
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                self.code,self.path,'stars',self.star_name_gastronoom,\
                'figures'))
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                self.code,self.path,'stars',self.star_name_gastronoom,\
                'statistics'))
        exec_file = DataIO.readFile(os.path.join(os.path.expanduser('~'),\
                                    'GASTRoNOoM','execGASTRoNOoM_example'))
        new_exec_file = []
        for index,line in enumerate(exec_file):
            if line.find('DATADIR=') == 0:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                               star['PATH_GAS_DATA']))
            elif line.find('WORKDIR=') == 0:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                               os.path.join(\
                                                    os.path.expanduser('~'),\
                                                    'GASTRoNOoM',\
                                                    self.path + '/')))
            elif line.find('STAR=') == 0:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                               self.star_name_gastronoom))
            elif line.find('MODEL=') == 0:
                new_exec_file.append('%s=%s' %(line.split('=')[0],sphinx_id))
            elif line.find('INDEX=') == 0:
                new_exec_file.append('%s=%s' %(line.split('=')[0],str(i)))
            elif line.find('NY_UP=') == 0 \
                    and star.getMolecule('12C16O') <> None:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                            star.getMolecule('12C16O').ny_up))
            elif line.find('NY_LOW=') == 0 \
                    and star.getMolecule('12C16O') <> None:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                            star.getMolecule('12C16O').ny_low))
            elif line.find('NLINE=') == 0 \
                    and star.getMolecule('12C16O') <> None:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                            star.getMolecule('12C16O').nline))
            elif line.find('N_IMPACT_12CO=') == 0 \
                    and star.getMolecule('12C16O') <> None:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                          star.getMolecule('12C16O').n_impact))
            elif line.find('N_IMPACT_13CO=') == 0 \
                    and star.getMolecule('13C16O') <> None:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                          star.getMolecule('13C16O').n_impact))
            elif line.find('RUNDIR=') == 0:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                               os.path.join(\
                                                    os.path.expanduser('~'),\
                                                    'GASTRoNOoM','src/')))
            elif line.find('SCRIPTDIR=') == 0:
                new_exec_file.append('%s=%s' %(line.split('=')[0],\
                                               os.path.join(\
                                                    os.path.expanduser('~'),\
                                                    'GASTRoNOoM','scripts/')))
            elif (line and line.find('$RUNDIR/exec/') == 0) \
                    or (exec_file[index-1] \
                        and exec_file[index-1].find('$RUNDIR/exec/') == 0) \
                    or (exec_file[index-2] \
                        and exec_file[index-2].find('$RUNDIR/exec/') == 0):
                new_exec_file.append('#' + line)
            else:
                new_exec_file.append(line)    
            DataIO.writeFile(os.path.join(os.path.expanduser('~'),\
                                          'GASTRoNOoM','execGASTRoNOoM'),\
                             new_exec_file)    



    def makeStars(self,models):
        
        '''
        Make a Star list based on either GASTRoNOoM cooling ids or PACS ids.
        
        @param models: model_ids for the MCMax db
        @type models: list(string)
        
        @return: the parameter sets
        @rtype: list[Star()]
        
        '''
        
        star_grid = Star.makeStars(models=models,star_name=self.star_name,\
                                   id_type='pacs' in models[0].lower() \
                                                and 'PACS' \
                                                or 'GASTRoNOoM',\
                                   path_combocode=self.path_combocode,\
                                   code='GASTRoNOoM',path=self.path)
        if 'pacs' in models[0].lower():
            self.pacs.addStarPars(star_grid)
        [star.addCoolingPars() for star in star_grid]
        return star_grid



    def plotVelocity(self,star_grid=[],models=[],cfg=''):
        
        '''
        Plot velocity versus radius for every model in star_grid.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]        
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Plotting Velocity Profiles'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        plot_filenames = []
        for i,star in enumerate(star_grid):
            if star['LAST_GASTRONOOM_MODEL']:    
                radius,vel = star.getGasVelocity()
                radius = radius/star['R_STAR']/star.r_solar
                vel = vel/10.**5
                avgdrift = star.getAverageDrift()/10.**5
                plot_filename = os.path.join(os.path.expanduser('~'),\
                                             'GASTRoNOoM',self.path,'stars',\
                                             self.star_name,self.plot_id,\
                                             'velocity_%s_%i'\
                                             %(star['LAST_GASTRONOOM_MODEL'],\
                                               i))
                plot_title = '%s %s: Velocity Profile for Model %i (%s)'\
                             %(self.plot_id.replace('_','\_'),\
                               self.star_name_plots,i,\
                               star['LAST_GASTRONOOM_MODEL'].replace('_','\_'))
                plot_filenames.append(Plotting2.plotCols(x=radius,\
                                        y=[vel,avgdrift],cfg=cfg,\
                                        filename=plot_filename,\
                                        xaxis='R (R$_*$)',\
                                        yaxis=r'$v$ (km s$^{-1}$)',\
                                        plot_title=plot_title,\
                                        key_location=(0.0,0.0),\
                                        keytags=['Gas Velocity',\
                                                 'Grain-size Weighted Drift'],\
                                        xlogscale=1))
        if plot_filenames and plot_filenames[0][-4:] == '.pdf':    
            new_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                        self.path,'stars',self.star_name,\
                                        self.plot_id,'velocity_profiles.pdf')
            DataIO.joinPdf(old=sorted(plot_filenames),new=new_filename)
            print '** Plots can be found at:'
            print new_filename
            print '***********************************' 
        elif plot_filenames:
            print '** Plots can be found at:'
            print '\n'.join(plot_filenames)
            print '***********************************' 
        else:
            print '** No GASTRoNOoM models were calculated successfully. '+\
                  'No velocity profiles can be plotted.'
            print '***********************************'


    
    def plotTemp(self,star_grid=[],models=[],force_plot=0,cfg=''):
        
        '''
        Plot temperature profiles of all models.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]
        @keyword force_plot: force a plotting if more than models are requested
                             
                             (default: 0)
        @type force_plot: bool
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Plotting Gas Temperature Profiles'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        if len(star_grid) < 20 or force_plot:
            valid_sg = [star 
                        for star in star_grid 
                        if star['LAST_GASTRONOOM_MODEL']]
            radii = [DataIO.getGastronoomOutput(\
                        filename=os.path.join(os.path.expanduser('~'),\
                             'GASTRoNOoM',self.path,'models',\
                             star['LAST_GASTRONOOM_MODEL'],\
                             'coolfgr_all%s.dat'\
                             %star['LAST_GASTRONOOM_MODEL']),\
                        keyword='RADIUS',return_array=1)
                     for star in valid_sg]
            radii_rstar = [rads/star['R_STAR']/star.r_solar 
                                for rads,star in zip(radii,valid_sg)]
            temp = [DataIO.getGastronoomOutput(\
                                filename=os.path.join(os.path.expanduser('~'),\
                                              'GASTRoNOoM',self.path,\
                                              'models',\
                                              star['LAST_GASTRONOOM_MODEL'],\
                                              'coolfgr_all%s.dat'\
                                              %star['LAST_GASTRONOOM_MODEL']),\
                                keyword='TEMP',return_array=1)
                    for star in valid_sg]

            if temp:    
                plot_title = '%s %s: Temperature Profiles'\
                             %(self.plot_id.replace('_','\_'),\
                               self.star_name_plots)
                extension='.eps'
                keytags = star_grid[0].has_key('LAST_PACS_MODEL') \
                            and ['%s,    %s,    Mdot = %.2e'\
                                 %(star['LAST_GASTRONOOM_MODEL']\
                                       .replace('_','\_'),\
                                   star['LAST_PACS_MODEL']\
                                       .replace('_','\_'),\
                                   float(star['MDOT_GAS'])) 
                                 for star in star_grid] \
                            or [star['LAST_GASTRONOOM_MODEL'].replace('_','\_')
                                for star in valid_sg]
                plot_filename = os.path.join(os.path.expanduser('~'),\
                                             'GASTRoNOoM',self.path,\
                                             'stars',self.star_name,\
                                             self.plot_id,\
                                             'temperature_profiles')
                plot_filename_rstar=os.path.join(os.path.expanduser('~'),\
                                                 'GASTRoNOoM',self.path,\
                                                 'stars',self.star_name,\
                                                 self.plot_id,\
                                                 'temperature_profiles_rstar')
                plot_filename_rstar = Plotting2.plotCols(x=radii_rstar,y=temp,\
                            cfg=cfg,xaxis='R (R$_*$)',plot_title='',\
                            filename=plot_filename_rstar,yaxis='T (K)',\
                            key_location=(0.0,0.0),xlogscale=1,ylogscale=1,\
                            keytags=keytags)
                keys_cm = ['Model %i'%(i+1)
                           for i in xrange(len(star_grid))]
                temp = [t[rad < 1e17] for t,rad in zip(temp,radii)]
                radii = [rad[rad < 1e17] for rad in radii]
                plot_filename = Plotting2.plotCols(x=radii,y=temp,cfg=cfg,\
                        filename=plot_filename,xaxis='$r$ (cm)',\
                        yaxis='$T_\mathrm{g}$ (K)',\
                        plot_title='',figsize=(12.5,8),fontsize_ticklabels=26,\
                        key_location=(0.05,0.05),xlogscale=1,ylogscale=1,\
                        keytags=keys_cm,extension=extension,fontsize_axis=26,\
                        linewidth=4,ymin=3,xmax=2e17,ymax=5000,fontsize_key=26)
                print '** Plots can be found at:'
                print plot_filename
                print plot_filename_rstar
                print '***********************************'
            else:
                print '** No GASTRoNOoM models were calculated successfully.'+\
                      'No temperature profiles can be plotted.'
                print '***********************************'



    def plotTransitionsOld(self,star_grid,iterative=0,cfg=''):
        
        """ 
        Creating the classical .ps output file of GASTRoNOoM; model and data!
                
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword iterative: When the iterative mode is on in the CC session
        
                            (default: 0)
        @type iterative: bool
        
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        """
        
        print '***********************************'
        print '** Creating Transition plots.'
        if iterative: pass
        else:
            filenames = []
            for i,star in enumerate(star_grid):
                print '***********************************'
                print '** Plotting model #%i out of %i requested models.'\
                      %(i+1,len(star_grid))
                if not star['LAST_GASTRONOOM_MODEL']: 
                    print '** No cooling model available. No plot is made.'
                elif not [trans 
                          for trans in star['GAS_LINES'] 
                          if trans.getModelId()]:
                    print '** No sphinx lines have been calculated. ' + \
                          'No plot is made.'
                else:
                    these_sphinx_ids = set([trans.getModelId() 
                                            for trans in star['GAS_LINES'] 
                                            if trans.getModelId()])
                    for sphinx_id in these_sphinx_ids:    
                        self.makeExec(star,i=i,sphinx_id=sphinx_id)
                        os.system(os.path.join(os.path.expanduser('~'),\
                                               'GASTRoNOoM',\
                                               './execGASTRoNOoM'))
                        filename = os.path.join(os.path.expanduser('~'),\
                                                'GASTRoNOoM',self.path,\
                                                'stars',\
                                                self.star_name_gastronoom,\
                                                'figures',\
                                                '%s_%i.ps'%(sphinx_id,i))
                        try:
                            plot_ex = open(filename)
                            plot_ex.close()
                            filenames.append(filename)
                            os.system('gv ' + filename + ' &')
                            print '** Transition plots can be found at:'
                            print filename
                        except IOError:
                            print 'Plotting failed for sphinx id %s. '\
                                  %sphinx_id + \
                                  'Double check output from GASTRoNOoM.'
            print '** All transition plots can be found at:'
            print '\n'.join(filenames)
            print '** DONE!'
            print '***********************************'



    def plotTransitions(self,star_grid,cfg='',no_data=0,vg_factor=3,\
                        telescope_label=1,sort_freq=0,sort_molec=0,\
                        no_models=0,limited_axis_labels=0,date_tag=1,\
                        n_max_models=10,fn_suffix='',mfiltered=0,\
                        plot_intrinsic=0):
        
        """ 
        Plotting beam convolved line profiles in Tmb for both model and data if 
        available.
      
        @param star_grid: list of Star objects for which the plotting is done.
        @type star_grid: list[Star]
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        @keyword no_data: Don't include the data
        @type no_data: bool
        @keyword vg_factor: The factor with which the terminal velocity is 
                            multiplied. This determines the xrange of the plots
        @type vg_factor: float
        @keyword telescope_label: Include a label showing the telescope name
        
                                  (default: 1)
        @type telescope_label: bool
        @keyword sort_freq: Sort the lines by frequency rather than wavelength.
                                  
                            (default: 0)
        @type sort_freq: bool
        @keyword sort_molec: Sort the lines by molecule. Can be combined with 
                             sort_freq
                                  
                             (default: 0)
        @type sort_molec: bool
        @keyword no_models: Only show data for the non-PACS lines.
                                  
                            (default: 0)
        @type no_models: bool
        @keyword limited_axis_labels: Remove axis labels not at the left or at 
                                      the bottom of the tiled plot
                                      
                                      (default: 0)
        @type limited_axis_labels: bool
        @keyword date_tag: Add a tag to a plot indicating the date of 
                           observation. Only available for non-intrinsic obs.
                         
                           (default: 1)
        @type date_tag: bool
        @keyword n_max_models: Maximum number of models per tile
        
                               (default: 10)
        @type n_max_models: bool
        @keyword fn_suffix: A suffix that is appended to the filename. For 
                            instance, when running the plot command for a 
                            best fit subgrid of Star() models as to not 
                            overwrite the plot of the full grid.
                            
                            (default: '')
        @type fn_suffix: string
        @keyword mfiltered: Show the models after filtering and scaled with
                            the best_vlsr value, instead of the sphinx output
                            
                            (default: 0)
        @type mfiltered: bool
        @keyword plot_intrinsic: Plot intrinsic line profiles as well, for PACS
                                 data. PACS data are not added to these plots.
                                 By default, this is off as the line profiles
                                 do not give you that much information before
                                 convolution with the wavelength resolution.
                                 
                                 (default: 0)
        @type plot_intrinsic: bool
        
        """
        
        print '***********************************'
        print '** Creating Transition plots.'
        #- Default dimension is (4,3), but can be adapted in cfg
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
        if cfg_dict.has_key('dimensions'):
            x_dim = int(cfg_dict['dimensions'][0])
            y_dim = int(cfg_dict['dimensions'][1])
        else:
            x_dim, y_dim = 4,3
        if cfg_dict.has_key('no_data'):
            no_data = bool(cfg_dict['no_data'])
        if cfg_dict.has_key('vg_factor'):
            vg_factor = float(cfg_dict['vg_factor'])
        if cfg_dict.has_key('telescope_label'):
            telescope_label = int(cfg_dict['telescope_label'])
        if cfg_dict.has_key('sort_freq'):
            sort_freq = int(cfg_dict['sort_freq'])
        if cfg_dict.has_key('sort_molec'):
            sort_molec = int(cfg_dict['sort_molec'])
        if cfg_dict.has_key('no_models'):
            no_models = int(cfg_dict['no_models'])
        if cfg_dict.has_key('limited_axis_labels'):
            limited_axis_labels = cfg_dict['limited_axis_labels']
        if cfg_dict.has_key('date_tag'):
            date_tag = int(cfg_dict['date_tag'])
        if cfg_dict.has_key('n_max_models'):
            n_max_models = int(cfg_dict['n_max_models'])
        if cfg_dict.has_key('plot_intrinsic'):
            plot_intrinsic = int(cfg_dict['plot_intrinsic'])
        if cfg_dict.has_key('mfiltered'):
            mfiltered = int(cfg_dict['mfiltered'])
        if fn_suffix: 
            filename = cfg_dict.get('filename',None)
            if filename <> None: filename = filename + '_%s'%fn_suffix
            cfg_dict['filename'] = filename
        if cfg_dict.has_key('keytags'):
            keytags = cfg_dict['keytags']
            pacs_keytags = keytags
        else:
            keytags = [',\\ '.join(set(\
                        [trans.getModelId() != '' \
                            and str(trans.getModelId())\
                                   .replace('_','\_')
                            or str(None)
                         for trans in star['GAS_LINES']]))
                       for star in star_grid]
            pacs_keytags = list(keytags)
            if not no_data : keytags.append('Data')
             
        #- Check how many non-PACS transitions there are
        trans_list = Transition.extractTransFromStars(star_grid,sort_freq,\
                                                      sort_molec,pacs=0)
        if plot_intrinsic:
            pacs_list  = Transition.extractTransFromStars(star_grid,sort_freq,\
                                                          sort_molec,pacs=1)
        else:
            pacs_list = []
            
        def createTilePlots(trans_list,x_dim,y_dim,no_data,intrinsic,\
                            vg_factor,keytags,telescope_label,no_models,cfg,\
                            star_grid,limited_axis_labels,date_tag,indexi,\
                            indexf,mfiltered):
            
            '''
            Create a tiled plot for a transition list.
            
            A list of transitions is exhausted for every star in star_grid,
            as long as tiles in a single plot are still available. 
            
            @param trans_list: The transition list. Transitions will be removed
                               from this list as tiles are created.
            @type trans_list: list[Transition()]
            @param x_dim: The number of tiles in the horizontal direction
            @type x_dim: int
            @param y_dim: The number of tiles in the vertical direction
            @type y_dim: int
            @param no_data: Include data or not? Will call a function that
                            gathers the data
            @type no_data: bool
            @param intrinsic: Intrinsic line profiles, or convolved with beam
                              profile? Set to True for the former, False for 
                              the latter
            @type intrinsic: bool
            @param star_grid: The grid of models for which the chosen
                              transitions in trans_list are plotted. Can be a 
                              subgrid of the original star_grid. (determined
                              before calling this method)
            @type star_grid: list[Star()]
            @param keytags: list of keys for the different models in star_grid
            @type keytags: list[string]
            @param cfg: The config filename passed to the plotTiles method
            @type cfg: string
            @param vg_factor: The factor with which the terminal velocity is 
                              multiplied. This determines the xrange of the 
                              plots
            @type vg_factor: float
            @param telescope_label: Include a label showing the telescope name
            @type telescope_label: bool
            @param no_models: Only show data for the non-PACS lines.
            @type no_models: bool
            @param limited_axis_labels: Remove axis labels not at the left or 
                                        at the bottom of the tiled plot
                                      
                                        (default: 0)
            @type limited_axis_labels: bool
            @param date_tag: Add a tag to a plot indicating the date of 
                             observation. Only available for non-intrinsic obs.
            @type date_tag: bool
            @param indexi: The start index of the models in the star_grid
            @type indexi: int
            @param indexf: The end index of the models in the star_grid
            @type indexf: int
            @param mfiltered: Show the models after filtering and scaled with
                              the best_vlsr value, instead of the sphinx output
            @type mfiltered: bool
            
            @return: The data list with dictionaries for every tile is returned
            @rtype: list[dict]

            '''
            
            missing_trans = 0
            n_subplots = (x_dim*y_dim) - (keytags and 1 or 0)
            plot_filenames = []
            i = 0
            vexp = max([s['VEL_INFINITY_GAS'] for s in star_grid])
            while trans_list:
                i += 1             
                data = []
                for j in xrange(n_subplots):
                    current_trans = trans_list.pop(0)
                    current_sub = [star.getTransition(current_trans) 
                                   for star in star_grid]
                    if None in current_sub: 
                         missing_trans += 1
                    #-- Just fit the line profile. The data will be read as well
                    if not intrinsic:
                        current_trans.fitLP()
                        vlsr = current_trans.getVlsr()
                        noise = current_trans.getNoise()
                    else:
                        vlsr = 0.0
                        noise = None
                    for trans in current_sub:
                        if trans <> None:
                            trans.readSphinx()
                            #-- Data have been read for current_trans. Don't 
                            #   read again for other objects (same data files),
                            #   but simply set based on the already read data.
                            #   Same with the profile fit results
                            trans.setData(current_trans)
                    ddict = dict()
                    if not no_models:
                        if intrinsic:
                            ddict['x'] = \
                                [(trans <> None and trans.sphinx <> None) \
                                      and list(trans.sphinx.getVelocityIntrinsic())\
                                      or []                          
                                 for trans in current_sub]
                            ddict['y'] = \
                                [(trans <> None and trans.sphinx <> None) \
                                     and list(trans.sphinx.getLPIntrinsic()*10**(23))\
                                     or []
                                 for trans in current_sub]
                        else:
                            ddict['x'] = []
                            ddict['y'] = []
                            for trans in current_sub:
                                bvlsr = trans.getBestVlsr()
                                if mfiltered and trans <> None \
                                        and trans.best_mfilter <> None:
                                    ddict['x'].append(trans.lpdata[0]\
                                                        .getVelocity())
                                    ddict['y'].append(trans.best_mfilter)
                                elif trans <> None and trans.sphinx <> None:
                                    ddict['x'].append(trans.sphinx\
                                                        .getVelocity() + bvlsr)
                                    ddict['y'].append(trans.sphinx.getLPTmb())
                                else:
                                    ddict['x'].append([])
                                    ddict['y'].append([])
                    else:
                        ddict['x'], ddict['y'] = [], []
                    #- Add data, but only if the data filename is known. This 
                    #- will be Tmb, in K. In case of intrinsic==1, you dont 
                    #- even want to check this.
                    
                    if current_trans.lpdata and no_data == 0:
                        ddict['histoplot'] = []
                        n_models = len(ddict['x'])
                        for ilp,lp in enumerate(current_trans.lpdata):
                            ddict['x'].append(lp.getVelocity())
                            ddict['y'].append(lp.getFlux())
                            ddict['histoplot'].append(n_models+ilp)
                        
                    ddict['labels'] = \
                        [('%s'%(current_trans.molecule.molecule_plot),0.05,0.87),\
                         ('%s'%(current_trans.makeLabel()),0.05,0.76)]
                    if telescope_label:
                        if True in [trans.sphinx.nans_present 
                                    for trans in current_sub
                                    if (trans <> None and trans.sphinx <> None)]:
                            telescope_string = '%s*'\
                                %current_trans.telescope.replace('-H2O','')\
                                                        .replace('-CORRB','')
                        else:
                            telescope_string = '%s'\
                                %current_trans.telescope.replace('-H2O','')\
                                                        .replace('-CORRB','')
                        ddict['labels'].append((telescope_string,0.73,0.85))
                    if not intrinsic and date_tag:
                        ddict['labels'].append(\
                            ('; '.join([lp.getDateObs() \
                                        for lp in current_trans.lpdata]),\
                             0.05,0.01))
                    #-- Don't use the fitted vexp for plotting window, keep it 
                    #   the same for all lines, in case higher lines are
                    #   narrower
                    ddict['xmax'] = vlsr + vg_factor * vexp
                    ddict['xmin'] = vlsr - vg_factor * vexp
                    if [yi for yi in ddict['y'] if list(yi)]:
                        ddict['ymax'] = max([max(array(yi)[(array(xi) <= ddict['xmax'])* \
                                                (array(xi) >= ddict['xmin'])]) 
                                             for xi,yi in zip(ddict['x'],\
                                                              ddict['y'])
                                             if list(yi)])*1.3
                        ddict['ymin'] = min([min(array(yi)[(array(xi) <= ddict['xmax'])* \
                                                (array(xi) >= ddict['xmin'])]) 
                                             for xi,yi in zip(ddict['x'],\
                                                              ddict['y'])
                                             if list(yi)])
                        if noise <> None and ddict['ymin'] < -3*noise: 
                            ddict['ymin'] = -3*noise
                    if limited_axis_labels:
                        if j%x_dim == 0:
                            ddict['yaxis'] = intrinsic \
                                                and r'$F_\nu$ (Jy)' \
                                                or '$T_\mathrm{mb}$ (K)'
                        else: 
                            ddict['yaxis'] = ''
                        if j in xrange(n_subplots-x_dim,n_subplots):
                            ddict['xaxis'] = r'$v$ (km s$^{-1}$)'
                        else:
                            ddict['xaxis'] = ''
                    data.append(ddict)
                    if not trans_list:
                        break
                extension = '.pdf'
                filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                        self.path,'stars',self.star_name,\
                                        self.plot_id,\
                                        '%sline_profiles_%i_models_%ito%i'\
                                        %(intrinsic and 'intrinsic_' or '',i,\
                                          indexi,indexf))
                plot_filenames.append(Plotting2.plotTiles(extension=extension,\
                     data=data,filename=filename,keytags=keytags,\
                     xaxis=r'$v$ (km s$^{-1}$)',fontsize_axis=16,cfg=cfg,\
                     yaxis=intrinsic \
                            and r'$F_\nu$ (Jy)' \
                            or '$T_\mathrm{mb}$ (K)',\
                     fontsize_ticklabels=16,dimensions=(x_dim,y_dim),\
                     fontsize_label=20,linewidth=2))
            if missing_trans:
                print 'WARNING! %i requested transitions were not found for a'\
                      %missing_trans+\
                      ' Star(). Within one CC session, this should not be '+\
                      'the case!'
            if plot_filenames and plot_filenames[0][-4:] == '.pdf':    
                new_filename = os.path.join(os.path.expanduser('~'),\
                                            'GASTRoNOoM',\
                                            self.path,'stars',self.star_name,\
                                            self.plot_id,\
                                            '%sline_profiles_models_%ito%i.pdf'\
                                            %(intrinsic and 'intrinsic_' or '',\
                                              indexi,indexf))
                DataIO.joinPdf(old=plot_filenames,new=new_filename)
                print '** %sine profile plots can be found at:'\
                        %(intrinsic and 'Intrinsic l' or 'L')
                print new_filename
                print '***********************************' 
            else:
                print '** %sine profile plots can be found at:'\
                        %(intrinsic and 'Intrinsic l' or 'L')
                print '\n'.join(plot_filenames)
                print '***********************************' 
    
        if trans_list:
            j = 0
            while j < len(star_grid):
                i = 0 
                subgrid = []
                subkeys = []
                while i < n_max_models and i+j < len(star_grid):
                    subgrid.append(star_grid[i+j])
                    if keytags: subkeys.append(keytags[i+j])
                    i += 1
                if keytags: subkeys.append(keytags[-1])
                #- Copying the list so that the destructive loop does not mess
                #- up multiple tile plot runs if n_models > n_max_models
                createTilePlots(trans_list=list(trans_list),\
                                vg_factor=vg_factor,\
                                no_data=no_data,cfg=cfg,star_grid=subgrid,\
                                x_dim=x_dim,y_dim=y_dim,keytags=subkeys,\
                                intrinsic=0,no_models=no_models,\
                                telescope_label=telescope_label,\
                                limited_axis_labels=limited_axis_labels,\
                                date_tag=date_tag,indexi=j,indexf=j+i-1,\
                                mfiltered=mfiltered)
                j += i
        if pacs_list:  
            j = 0
            while j < len(star_grid):
                i = 0 
                subgrid = []
                subkeys = []
                while i < n_max_models and i+j < len(star_grid):
                    subgrid.append(star_grid[i+j])
                    if keytags: subkeys.append(pacs_keytags[i+j])
                    i += 1
                createTilePlots(trans_list=list(pacs_list),\
                                vg_factor=vg_factor,\
                                no_data=1,cfg=cfg,star_grid=subgrid,\
                                intrinsic=1,keytags=subkeys,\
                                x_dim=x_dim,y_dim=y_dim,date_tag=date_tag,\
                                telescope_label=telescope_label,no_models=0,\
                                limited_axis_labels=limited_axis_labels,\
                                indexi=j,indexf=j+i-1,mfiltered=mfiltered)
                j += i            
                


    def createLineLabels(self,star,xmin,xmax):
        
        '''
        Create a list of line labels for all molecules and transitions 
        in the molecular databases requested.
        
        This is used for spectroscopic databases only! Such as JPL, CDMS, LAMDA
        
        @param star: The parameter set
        @type star: Star()
        @param xmin: minimum wavelength
        @type xmin: float
        @param xmax: maximum wavelength
        @type xmax: float
        
        @return: The labels with x location and a molecule index.
        @rtype: list[string, float, index]
        
        '''
        
        cdms = int(star['LL_CDMS'])
        jpl = int(star['LL_JPL'])
        lamda = int(star['LL_LAMDA'])
        min_strength = float(star['LL_MIN_STRENGTH']) \
                            and float(star['LL_MIN_STRENGTH']) or None
        max_exc = float(star['LL_MAX_EXC']) \
                        and float(star['LL_MAX_EXC']) or None
        path = star['LL_PATH']
        xunit = star['LL_UNIT']
        linelisters = []
        if xunit.lower() != 'micron': 
            raise IOError('Units other than micron for the LineLists not ' + \
                          'yet supported!')
        for molecule in star['LL_GAS_LIST']:
            if not 'p1H' in molecule.molecule:
                linelisters.append(\
                    LineList.LineList(molecule=molecule,x_min=xmin,x_max=xmax,\
                                      x_unit=xunit,cdms=cdms,jpl=jpl,\
                                      lamda=lamda,path=path,\
                                      min_strength=min_strength,\
                                      max_exc=max_exc,include_extra=1))
        line_labels = sorted([('%s %s'%(ll.molecule.molecule,\
                                        trans.makeLabel()),\
                               trans.wavelength*10**4,i)  
                              for i,ll in enumerate(linelisters)
                              for trans in ll.makeTransitions()
                              if trans.wavelength*10**4 >= xmin \
                                    and trans.wavelength*10**4 <= xmax],\
                             key=operator.itemgetter(1))
        return line_labels     
    
    
    
    def plotLineLists(self,star_grid,include_sphinx=1,cfg=''):
        
        '''
        Plot linelists along with the indicated data.
        
        @param star_grid: The Parameter sets
        @type star_grid: list[Star()]
        
        @keyword include_sphinx: Include convolved Sphinx models in the plots 
                                 for the star_grid
                                 
                                 (default: 1)
        @type include_sphinx: bool
        @keyword cfg: path to the Plotting2.plotCols config file. If default, the
                          hard-coded default plotting options are used.
                          
                          (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Starting to plot line identifications for %s from databases.'\
              %self.star_name
        if self.pacs is None: 
            print '** No PATH_PACS given. Cannot plot linelists without ' + \
                  'data information. Aborting...'
            return
        xmins = [min(wave_list) for wave_list in self.pacs.data_wave_list]
        xmaxs = [max(wave_list) for wave_list in self.pacs.data_wave_list]
        line_labels = self.createLineLabels(star_grid[0],min(xmins),max(xmaxs))
        plot_filenames = []
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                'GASTRoNOoM',self.path,\
                                                'stars',self.star_name,\
                                                self.plot_id,'LineLists'))
        if include_sphinx:
            if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) \
                            == set([0]) \
                    or set([s['LAST_GASTRONOOM_MODEL'] for s in star_grid]) \
                            == set(['']): 
                include_sphinx = 0
            else: 
                self.setSphinx(star_grid)
        for i_file,(wave,flux,filename,xmin,xmax) in enumerate(\
                    zip(self.pacs.data_wave_list,self.pacs.data_flux_list,\
                        self.pacs.data_filenames,xmins,xmaxs)):
            if include_sphinx:
                sphinx_flux = [sphinx 
                               for sphinx in self.sphinx_flux_list[i_file] 
                               if list(sphinx)]
            else:
                sphinx_flux = []
            these_line_labels = [(label,wl,index) 
                                 for label,wl,index in line_labels 
                                 if wl >= xmin and wl <= xmax]
            plot_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                         self.path,'stars',self.star_name,\
                                         self.plot_id,'LineLists',\
                                         'line_id_'\
                                         +os.path.split(filename)[1]\
                                                    .replace('.dat',''))
            keytags = ['PACS %s'%filename.replace('_','\_')] + \
                      ['Model %i: %s'%(i+1,star['LAST_PACS_MODEL']\
                                            .replace('_','\_')) 
                       for i,star in enumerate(star_grid) 
                       if star['LAST_GASTRONOOM_MODEL'] and include_sphinx]
            plot_filenames.append(Plotting2.plotCols(\
                    x=[wave]*(len(sphinx_flux)+1),y=[flux]+sphinx_flux,\
                    cfg=cfg,filename=plot_filename,keytags=keytags,\
                    plot_title=self.star_name_plots,histoplot=[0],\
                    number_subplots=3,line_labels=these_line_labels,\
                    line_label_color=1,line_label_lines=1,line_label_spectrum=1))
        new_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                    self.path,'stars',self.star_name,\
                                    self.plot_id,'LineLists',\
                                    'line_id_pacs_%s.pdf'%self.star_name)
        DataIO.joinPdf(old=sorted(plot_filenames),new=new_filename)
        print '** Plots can be found at:'
        print new_filename
        print '***********************************'

                                                
   
    def plotAbundanceProfiles(self,star_grid=[],models=[],cfg=''):  
        
        '''
        Plot abundance profiles for all molecules in every model.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Plotting Abundance Profiles'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        plot_filenames = []
        for star in star_grid:
            if star['LAST_GASTRONOOM_MODEL']:    
                radii = [DataIO.getGastronoomOutput(\
                                filename=os.path.join(os.path.expanduser('~'),\
                                      'GASTRoNOoM',self.path,'models',\
                                      molec.getModelId(),'cool1%s_%s.dat'\
                                      %(molec.getModelId(),molec.molecule)),\
                                return_array=1)
                         for molec in star['GAS_LIST'] 
                         if molec.getModelId()]
                h2_ab = [array(DataIO.getGastronoomOutput(\
                                filename=os.path.join(os.path.expanduser('~'),\
                                        'GASTRoNOoM',self.path,'models',\
                                        molec.getModelId(),'cool1%s_%s.dat'\
                                        %(molec.getModelId(),molec.molecule)),\
                                keyword='N(H2)'))
                         for molec in star['GAS_LIST'] 
                         if molec.getModelId()]
                molec_ab = [array(DataIO.getGastronoomOutput(\
                                filename=os.path.join(os.path.expanduser('~'),\
                                        'GASTRoNOoM',self.path,'models',\
                                        molec.getModelId(),'cool1%s_%s.dat'\
                                        %(molec.getModelId(),molec.molecule)),\
                                keyword='N(MOLEC)',key_index=8))
                            for molec in star['GAS_LIST'] 
                            if molec.getModelId()]
                abuns = []
                for rad,molec,amol,ah in zip(radii,star['GAS_LIST'],\
                                             molec_ab,h2_ab):
                     if molec.set_keyword_change_abundance:
                          rfrac,frac = tuple(DataIO.readCols(\
                                            molec.change_fraction_filename))
                          rfrac = rfrac/star.r_solar/star['R_STAR']
                          frac_interpol = interp1d(rfrac,frac)(rad)
                          #- if error happens, catch and print out warning, plus run
                          #- interpolation again with bounds_error=False, fill_value=frac[-1]
                          #- if bounds_error=False and a warning is printed by scipy, 
                          #- then no need to catch error first
                          #- frac_interpol = array(Interpol.doInterpol(x_in=rfrac,\
                          #-            y_in=frac,gridsx=[rad])[0][0])
                     else:
                          frac_interpol = 1
                     abuns.append(amol/ah/molec.abun_factor*frac_interpol) 
                keytags = ['%s'%(molec.molecule_plot) 
                           for molec in star['GAS_LIST'] 
                           if molec.getModelId()]
                linestyles = ['-','-','-','-','-','-',\
                              '--','--','--','--','--','--',\
                              '.-','.-','.-','.-','.-','.-']
                colors = ['r','b','g','m','y','c']
                lt = [ls + col for ls,col in zip(linestyles,3*colors)]
                lt = lt[:len(radii)]
                if star.has_key('R_DES_H2O'):
                    radii.extend([array([star['R_DES_H2O'],star['R_DES_H2O']])])
                    abuns.extend([[1e-2,1e-9]])
                    keytags.append('Condensation radius H$_2$O ice')
                    lt.append('--k')
                if star['R_OH1612']:
                    radii.extend([array([star['R_OH1612'],star['R_OH1612']])])
                    abuns.extend([[1e-2,1e-9]])
                    keytags.append('Max intensity OH maser')
                    lt.append('-k')
                extension = '.eps'
                plot_filename = os.path.join(os.path.expanduser('~'),\
                                    'GASTRoNOoM',self.path,'stars',\
                                    self.star_name,self.plot_id,\
                                    'abundance_profiles_%s'\
                                    %'_'.join(list(set([molec.getModelId() 
                                            for molec in star['GAS_LIST']]))))
                radii = [r*star['R_STAR']*star.r_solar for r in radii]
                plot_filenames.append(\
                     Plotting2.plotCols(x=radii,y=abuns,xaxis='$r$ (cm)',\
                          filename=plot_filename,keytags=keytags,\
                          figsize=(12.5,8.5),cfg=cfg,\
                          yaxis='$n_\mathrm{molec}/n_{\mathrm{H}_2}$',\
                          line_types=lt,transparent=0,\
                          extension=extension,fontsize_axis=28,linewidth=4,\
                          ylogscale=1,key_location=(.05,.1),\
                          fontsize_ticklabels=26,\
                          xlogscale=1,ymin=10**(-5),ymax=10**(-3),\
                          fontsize_key=18,xmax=1e18,\
                          xmin=star['R_INNER_GAS']*star['R_STAR']*star.r_solar))
        if plot_filenames and plot_filenames[0][-4:] == '.pdf':    
            new_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                        self.path,'stars',self.star_name,\
                                        self.plot_id,'abundance_profiles.pdf')
            DataIO.joinPdf(old=plot_filenames,new=new_filename)
            print '** Plots can be found at:'
            print new_filename
            print '***********************************'
        elif plot_filenames:
            print '** Plots can be found at:'
            print '\n'.join(plot_filenames)
            print '***********************************'
        else:
            print '** No GASTRoNOoM models were calculated successfully. ' + \
                    'No temperature profiles can be plotted.'
            print '***********************************'
            
            

    def plotLineContributions(self,star_grid,normalized=1,cfg='',do_sort=1,\
                              include_velocity=1):
        
        '''
        Plot the source function as function of impact parameter for every 
        transition.
        
        @param star_grid: The model parameter sets
        @type star_grid: list[Star()]
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        @keyword normalized: plot the normalized source functions as opposed 
                             to not normalized
                             
                             (default: 1)
        @type normalized: bool
        @keyword do_sort: Sort the transition list according to wavelength. If 
                          off, the original order given in the CC input file is 
                          kept
                          
                          (default: 1)
        @type do_sort: bool
        @keyword include_velocity: Include the velocity profile on the plot
                          
                                   (default: 0)
        @type include_velocity: bool
        
        '''
        
        print '***********************************'
        print '** Plotting Line Contributions'
        if cfg:
             cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
             cfg_dict = dict()
        if cfg_dict.has_key('do_sort'):
             do_sort = int(cfg_dict['do_sort'])
        if cfg_dict.has_key('normalized'):
             normalized = int(cfg_dict['normalized'])
        if cfg_dict.has_key('include_velocity'):
             include_velocity = int(cfg_dict['include_velocity'])
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                'GASTRoNOoM',self.path,'stars',\
                                                self.star_name,self.plot_id,\
                                                'LineContributions'))
        normalized = int(normalized)
        for i,star in enumerate(star_grid):
            extra_pars = dict()
            if do_sort:
                transitions = sorted([trans 
                                      for trans in star['GAS_LINES'] 
                                      if trans.getModelId()],\
                                     key=lambda x:x.wavelength)
            else:
                transitions = sorted([trans 
                                      for trans in star['GAS_LINES'] 
                                      if trans.getModelId()],\
                                     key=lambda x:trans_sel.index(x))
            if include_velocity:
                radius,vel = star.getGasVelocity()
                radius = radius/star['R_STAR']/star.r_solar
                vel = vel/10.**5
                extra_pars['twiny_x'] = [radius]
                extra_pars['twiny_y'] = [vel]
                extra_pars['twiny_keytags'] = [r'$v_\mathrm{g}$']
                extra_pars['twinyaxis'] = r'$v_\mathrm{g}$ (km s$^{-1}$)' 
            [trans.readSphinx() for trans in transitions]
            radii = [trans.sphinx.getImpact() for trans in transitions]
            linecontribs =  [normalized \
                                and list(trans.sphinx.getNormalizedIntensity())\
                                or list(trans.sphinx.getWeightedIntensity())
                             for trans in transitions]
            plot_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                         self.path,'stars',self.star_name,\
                                         self.plot_id,'LineContributions',\
                                         'linecontrib_%s_%i'\
                                         %(star['LAST_GASTRONOOM_MODEL'],i))
            keytags=['%s: %s'%(trans.molecule.molecule,trans.makeLabel())
                     for trans in transitions]
            ymin = normalized and -0.01 or None
            ymax = normalized and 1.02 or None
            plot_title = '%s: Line Contributions for %s'\
                         %(self.star_name_plots,\
                           star['LAST_PACS_MODEL'].replace('_','\_') \
                            or star['LAST_GASTRONOOM_MODEL'].replace('_','\_'))
            plot_filename = Plotting2.plotCols(\
                x=radii,y=linecontribs,filename=plot_filename,cfg=cfg,\
                xaxis='$p$ (R$_*$)',yaxis='$I(p) \\times g(p^2)$',\
                plot_title=plot_title,keytags=keytags,extension='.png',\
                xlogscale=1,figsize=(20,10),fontsize_axis=26,fontsize_title=22,\
                key_location=(0.70,0.1),fontsize_ticklabels=22,\
                fontsize_key=20,linewidth=3,ymin=ymin,ymax=ymax,**extra_pars)
            print '** Plot can be found at:'
            print plot_filename
            print '***********************************'                            
        

        
    def setSphinx(self,star_grid,refresh_sphinx_flux=0):
        
        ''' 
        Prepare Sphinx output in Pacs format (including convolution).
        
        @param star_grid: Parameter sets
        @type star_grid: list[Star()]
        @keyword refresh_sphinx_flux: redo the sphinx flux list by pulling from
                                      db, regardless of whether it's already
                                      been done or not.
                                      
                                      (default: 0)
        @type refresh_sphinx_flux: bool
        
        '''
        
        if not self.sphinx_flux_list or refresh_sphinx_flux:    
            self.pacs.prepareSphinx(star_grid)
            #- The sphinx convolved mosdels always have the same wavelength 
            #- list as their respective data files
            self.sphinx_flux_list = [[os.path.isfile(os.path.join(\
                    os.path.expanduser('~'),'GASTRoNOoM',self.path,'stars',\
                    self.star_name,'PACS_results',star['LAST_PACS_MODEL'],\
                    '_'.join(['sphinx',os.path.split(filename)[1]]))) \
                and DataIO.readCols(os.path.join(os.path.expanduser('~'),\
                            'GASTRoNOoM',self.path,'stars',self.star_name,\
                            'PACS_results',star['LAST_PACS_MODEL'],\
                            '_'.join(['sphinx',os.path.split(filename)[1]])),\
                            make_array=0)[1]
                or []
                    for star in star_grid]
                    for filename in self.pacs.data_filenames]
            self.sphinx_flux_list = [[array(sphinx) 
                                      for sphinx in sphinx_flux] 
                                     for sphinx_flux in self.sphinx_flux_list]



    def createSphinxLineLabels(self,star_grid):
         
         '''
         Create line labels for all transitions in list of Star() models.
         
         @param star_grid: The Star() models
         @type star_grid: list[Star()]
         @return: a sorted list(set) of line labels
         @rtype: list[string]
         
         '''
         
         sphinx_line_labels = list(set(\
                        [('%s %s'\
                          %(trans.molecule.molecule,trans.makeLabel()),\
                          star.c/trans.frequency*10**(4),\
                          trans.molecule.molecule_index)
                         for star in star_grid 
                         for trans in star['GAS_LINES'] 
                         if trans.getModelId()]))
         sphinx_line_labels = sorted(sphinx_line_labels,\
                                     key=operator.itemgetter(1))
         return sphinx_line_labels
         
         

    def plotPacsFull(self,star_grid=[],models=[],no_data=0,cfg=''):
         
        '''
        Plot the full PACS spectrum for a series of models along with the data.
        
        @keyword star_grid: star models for which PACS data will be fetched, 
                            default occurs when model_ids are passed instead, 
                            ie outside a CC modeling session
                               
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: list of pacs_ids or gastronoom model ids, default if 
                         Star models are passed instead
                               
                         (default: [])
        @type models: list[strings]
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                       
                      (default: '')
        @type cfg: string         
        @keyword no_data: if enabled only the sphinx mdels are plotted.
         
                          (default: 0)
        @type no_data: bool
        
        '''
         
        print '***********************************'
        print '** Creating full PACS spectrum plot.'
        if self.pacs is None: 
            print '** No PATH_PACS given. Cannot plot PACS spectra without '+\
                  'data information. Aborting...'
            return
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) == set([0]): 
            return
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
            'GASTRoNOoM',self.path,'stars',self.star_name,self.plot_id,\
            'PACS_results'))
        self.setSphinx(star_grid)
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
        if cfg_dict.has_key('no_data'):
            no_data = bool(cfg_dict['no_data'])
        print '** Plotting now...'
        sphinx_line_labels = self.createSphinxLineLabels(star_grid)
        tiles = []
        
        for wave,flux,sphinx_flux,filename,ordername in \
                        zip(self.pacs.data_wave_list,self.pacs.data_flux_list,\
                            self.sphinx_flux_list,self.pacs.data_filenames,\
                            self.pacs.data_ordernames):
            ddict = dict()
            ddict['x'] = no_data \
                            and [wave]*(len(sphinx_flux)) \
                            or [wave]*(len(sphinx_flux)+1)
            ddict['y'] = no_data and sphinx_flux or [flux]+sphinx_flux
            ddict['labels'] = [(ordername,0.01,0.01)]
            ddict['xmin'] = wave[0]
            ddict['xmax'] = wave[-1]
            ddict['histoplot'] = (not no_data) and [0] or []
            ddict['line_labels'] = [('',x,index) 
                                    for l,x,index in sphinx_line_labels
                                    if x >= wave[0] and x <= wave[-1]]
            tiles.append(ddict)
        filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',self.path,\
                                'stars',self.star_name,self.plot_id,\
                                'PACS_results','PACS_spectrum_full')
        tiles = sorted(tiles,key = lambda x: x['x'][0][0])
        filename = Plotting2.plotTiles(data=tiles,filename=filename,cfg=cfg,\
                                       line_label_color=1,\
                                       line_label_spectrum=1,fontsize_label=30,\
                                       line_label_lines=1,no_line_label_text=1,\
                                       dimensions=\
                                            (1,len(self.pacs.data_filenames)))
        print '** Your plot can be found at:'
        print filename
        print '***********************************'

       

    def plotPacs(self,star_grid=[],models=[],exclude_data=0,fn_plt='',cfg='',\
                 fn_trans_marker='',include_ordername=1):
        
        '''
        Plot PACS data along with Sphinx results, one plot per band.
        
        @keyword star_grid: star models for which PACS data will be fetched, 
                            default occurs when model_ids are passed instead, 
                            ie outside a CC modeling session
                                    
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: list of pacs_ids or gastronoom model ids, default if 
                         Star models are passed instead
                                
                         (default: [])
        @type models: list[strings]            
        @keyword exclude_data: if enabled only the sphinx mdels are plotted.
        
                               (default: 0)
        @type exclude_data: bool
        @keyword fn_plt: A plot filename to which an index is added for each
                         subband.
                         
                         (default: '')
        @type fn_plt: string
        @keyword cfg: path to the Plotting2.plotCols config file. If default, the
                      hard-coded default plotting options are used.
                         
                      (default: '')
        @type cfg: string
        @keyword fn_trans_marker: A file that includes TRANSITION definitions.
                                  These transitions will be marked up in the 
                                  plot. For instance, when indicating a subset 
                                  of transitions for one reason or another.
                                  
                                  (default: '')
        @type fn_trans_marker: string
        @keyword include_ordername: Include a name tag for the band order in 
                                    the plot.
                                    
                                    (default: 1)
        @type include_ordername: bool
        
        '''
        
        print '***********************************'
        print '** Creating PACS + Sphinx plot.'
        if self.pacs is None: 
            print '** No PATH_PACS given. Cannot plot PACS spectra without data'+\
                  ' information. Aborting...'
            return
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) == set([0]): 
            return
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
               'GASTRoNOoM',self.path,'stars',self.star_name,self.plot_id,\
               'PACS_results'))
        self.setSphinx(star_grid)
        print '** Plotting now...'
        sphinx_line_labels = self.createSphinxLineLabels(star_grid)
        plot_filenames = []
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        if cfg_dict.has_key('exclude_data'):
            exclude_data = bool(cfg_dict['exclude_data'])
        if cfg_dict.has_key('fn_trans_marker'):
            fn_trans_marker = cfg_dict['fn_trans_marker']
        if cfg_dict.has_key('include_ordername'):
            include_ordername = bool(cfg_dict['include_ordername'])
        if fn_trans_marker:
            lines = [line.split() 
                     for line in DataIO.readFile(fn_trans_marker) 
                     if line[0] != '#']
            trans_markers = set([Transition.makeTransition(star_grid[0],line) 
                                 for line in lines])
            used_indices = list(set([ll[-1] for ll in sphinx_line_labels]))
            this_index = [ii for ii in range(100) if ii not in used_indices][0]
            extra_line_labels = [('---',t.wavelength*10**4,this_index) 
                                 for t in trans_markers]
            sphinx_line_labels = sphinx_line_labels + extra_line_labels
        for wave,flux,sphinx_flux,filename,ordername in \
                    zip(self.pacs.data_wave_list,\
                        self.pacs.data_flux_list,\
                        self.sphinx_flux_list,\
                        self.pacs.data_filenames,\
                        self.pacs.data_ordernames):
            if fn_plt:
                fn_plt = os.path.splitext(fn_plt)[0]
                this_filename = '%s_%s'%(fn_plt,ordername)
            else:    
                this_filename = os.path.join(os.path.expanduser('~'),\
                                             'GASTRoNOoM',self.path,'stars',\
                                             self.star_name,self.plot_id,\
                                             'PACS_results',\
                                             os.path.split(filename)[1]\
                                                    .replace('.dat',''))
            keytags = ['Model %i: %s'%(i+1,star['LAST_PACS_MODEL']\
                                .replace('_','\_')) 
                       for i,star in enumerate(star_grid)]
            if exclude_data:
                x_list = [wave]*(len(sphinx_flux)) 
                y_list = sphinx_flux
            else:
                x_list = [wave]*(len(sphinx_flux)+1)
                y_list = [flux]+sphinx_flux
                keytags = ['PACS Spectrum'] + keytags
            if include_ordername:
                labels = [(ordername,0.01,0.01)]
            plot_filenames.append(Plotting2.plotCols(x=x_list,y=y_list,\
                    keytags=keytags,number_subplots=3,cfg=cfg_dict,\
                    plot_title='%s: %s - %s'%(self.plot_id.replace('_','\_'),\
                    self.star_name_plots,ordername),\
                    line_labels=sphinx_line_labels,\
                    histoplot=not exclude_data and [0] or [],\
                    filename=this_filename,labels=labels,\
                    line_label_spectrum=1,line_label_color=1))
        if plot_filenames and plot_filenames[0][-4:] == '.pdf':
            if not fn_plt:
                newf = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                    self.path,'stars',self.star_name,\
                                    self.plot_id,'PACS_results',\
                                    'PACS_spectrum.pdf')
            else:
                newf = '%s.pdf'%fn_plt
            DataIO.joinPdf(old=sorted(plot_filenames),new=newf,\
                           del_old=not fn_plt)
            print '** Your plots can be found at:'
            print newf
            print '***********************************'
        else:
            print '** Your plots can be found at:'
            print '\n'.join(plot_filenames)
            print '***********************************'



    def plotPacsSegments(self,star_grid,pacs_segments_path='',mode='sphinx',\
                         include_sphinx=None,no_data=0,cfg=''):
        
        '''
        Plot segments of spectra only.
        
        An inputfile gives the wavelength ranges, given by pacs_path_path.
        
        Can include the sphinx results overplotted with the data, as well as line
        labels generated either for sphinx results (mode == 'sphinx') or from a 
        spectroscopic database (mode == 'll').
                
        @param star_grid: star models for which PACS data will be fetched, 
        @type star_grid: list(Star())
        
        @keyword pacs_segments_path: The path to the file listing pairs of 
                                     wavelength ranges for plotting the 
                                     segments. This par can be passed through 
                                     the cfg file as well. 
        @type pacs_segments_path: string
        @keyword mode: the mode in which this method is used, the string is  
                       added to the outputfilename, can be 'sphinx' or 'll' for
                       now, determines the type of line labels. 'll' gives line
                       labels generated from a spectroscopic database. 'sphinx'
                       gives line labels for all transitions in all Star 
                       objects in star_grid. Can be passed through the cfg file.
                         
                       (default: sphinx)
        @type mode: string
        @keyword include_sphinx: Add sphinx results to the plots
        
                                 (default: None)
        @type include_sphinx: bool
        @keyword no_data: The data are not included when True. 
                                
                          (default: 0) 
        @type no_data: bool
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
        if cfg_dict.has_key('mode'):
            mode = cfg_dict['mode']
        if cfg_dict.has_key('pacs_segments_path'):
            pacs_segments_path = cfg_dict['pacs_segments_path']
        
        if mode == 'll':
            xmins=[min(wave_list) for wave_list in self.pacs.data_wave_list]
            xmaxs=[max(wave_list) for wave_list in self.pacs.data_wave_list]
            line_labels = self.createLineLabels(star_grid[0],min(xmins),\
                                                max(xmaxs))
        elif mode == 'sphinx':
            line_labels = list(set([(trans.molecule.molecule_plot,\
                                        star.c/trans.frequency*10**(4),\
                                        trans.molecule.molecule_index)
                                    for star in star_grid
                                    for trans in star['GAS_LINES']
                                    if trans.getModelId()]))
            line_labels = sorted(line_labels,key=operator.itemgetter(1))
        else:
            print 'Mode for plotting PACS segments not recognized. Aborting...'
            return
        
        if self.pacs is None:
            print 'No PACS data found for plotting PACS segments. Aborting...'
            return
        elif not pacs_segments_path:
            print 'No pacs_segments_path given. Pass in the cfg file or in ' +\
                  'the method call. Aborting...' 
            return
        else:
            self.setSphinx(star_grid)
        
        if include_sphinx is None:
            include_sphinx = self.sphinx_flux_list and 1 or 0
        print '** Plotting spectral segments.'
        for index, (wmin, wmax) in \
                enumerate(zip(*DataIO.readCols(pacs_segments_path))):
            delta = (wmax-wmin)/2.
            for i_file,(wave,flux,filename) in \
                     enumerate(zip(self.pacs.data_wave_list,\
                                   self.pacs.data_flux_list,\
                                   self.pacs.data_filenames)):
                if wmin > wave[0] and wmax < wave[-1]:
                    flux = flux[abs(wave-((wmax+wmin)/2.))<=delta]
                    if not include_sphinx: sphinx_flux = []
                    else: 
                        sphinx_flux = \
                            [f[abs(wave-((wmax+wmin)/2.))<=delta] 
                            for f in self.sphinx_flux_list[i_file] if list(f)]
                    wave = wave[abs(wave-((wmax+wmin)/2.))<=delta]
                    plot_filename = os.path.join(os.path.expanduser('~'),\
                          'GASTRoNOoM',self.path,'stars',self.star_name,\
                          self.plot_id,\
                          mode=='ll' and 'LineLists' or 'PACS_results',\
                          '%s_segment_%.1f-%.1f_'%(mode,wmin,wmax)+ \
                          os.path.split(filename)[1].replace('.dat',''))
                    extra_stats = dict([('line_labels',line_labels),\
                                        ('histoplot',not no_data \
                                                        and [0] \
                                                        or []),\
                                        ('filename',plot_filename)])
                    plot_filename = Plotting2.plotCols(\
                        x=[wave]*(len(sphinx_flux)+(not no_data and 1 or 0)),\
                        y=no_data and sphinx_flux or [flux]+sphinx_flux,\
                        cfg=cfg,**extra_stats)
                    print '** Segment finished and saved at:'
                    print plot_filename
                    
                    

    def plotSpire(self,star_grid=[],models=[],exclude_data=0,cfg=''):
        
        '''
        Plot SPIRE data along with Sphinx results.
        
        @keyword star_grid: star models for which SPIRE data will be fetched, 
                            default occurs when model_ids are passed instead, 
                            ie outside a CC modeling session
                                    
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: list of gastronoom model ids, default if 
                         Star models are passed instead
                                
                         (default: [])
        @type models: list[strings]            
        @keyword exclude_data: if enabled only the sphinx mdels are plotted.
        
                               (default: 0)
        @type exclude_data: bool
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Creating SPIRE + Sphinx plot.'
        if self.spire is None: 
            print '** No PATH_SPIRE given. Cannot plot SPIRE spectra '+\
                  'without data information. Aborting...'
            return
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) == set([0]): 
            return
        print '** Plotting now...'
        sphinx_line_labels = self.createSphinxLineLabels(star_grid)
        plot_filenames = []
        for wave,flux,filename in zip(self.spire.data_wave_list,\
                                      self.spire.data_flux_list,\
                                      self.spire.data_filenames):
            this_filename = os.path.join(os.path.expanduser('~'),\
                'GASTRoNOoM',self.path,'stars',self.star_name,self.plot_id,\
                'spire_spectrum_%s'\
                %os.path.split(filename)[1].replace('.dat',''))
            x_list = exclude_data \
                          and [wave]*(len(sphinx_flux)) \
                          or [wave]*(len(sphinx_flux)+1)
            sphinx_flux = [self.spire.sphinx_convolution\
                                [(star['LAST_GASTRONOOM_MODEL'],i)][filename] 
                           for i,star in enumerate(star_grid)]
            y_list = exclude_data \
                          and sphinx_flux \
                          or sphinx_flux + [flux]
            keytags = ['Model %i: %s'%(i+1,star['LAST_GASTRONOOM_MODEL']\
                                    .replace('_','\_')) 
                       for i,star in enumerate(star_grid)]
            if not exclude_data: 
                keytags = keytags + ['Spire Spectrum']
            plot_filenames.append(Plotting2.plotCols(x=x_list,y=y_list,\
                keytags=keytags,number_subplots=3,line_label_color=1,\
                plot_title='%s: %s' %(self.plot_id.replace('_','\_'),\
                self.star_name_plots),line_labels=sphinx_line_labels,\
                histoplot=[len(sphinx_flux)],filename=this_filename,cfg=cfg,\
                line_label_spectrum=1))
        print '** Your plots can be found at:'
        print '\n'.join(plot_filenames)
        print '***********************************'

