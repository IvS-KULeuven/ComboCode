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
import numpy as np

import cc.path
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
    
    def __init__(self,star_name,path_gastronoom='codeJun2013',\
                 inputfilename=None,pacs=None,spire=None):
        
        """ 
        Initializing an instance of PlotGas.
        
        @param star_name: name of the star from Star.dat, use default only 
                          when never using any star model specific things 
        @type star_name: string
        
        @keyword path_gastronoom: Output modeling folder in MCMax home folder
        
                                  (default: 'codeJun2013')
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
                                      path=path_gastronoom,\
                                      code='GASTRoNOoM',\
                                      inputfilename=inputfilename)
        #-- Convenience path
        cc.path.gout = os.path.join(cc.path.gastronoom,self.path)
        self.pacs = pacs
        if self.pacs:
            fn = os.path.join(cc.path.gout,'stars',self.star_name,\
                              self.plot_id,'PACS_results')
            DataIO.testFolderExistence(fn)
        self.spire = spire
        self.sphinx_flux_list = []

        
        
    def makeStars(self,models):
        
        '''
        Make a Star list based on either GASTRoNOoM cooling ids or PACS ids.
        
        @param models: model_ids for the MCMax db
        @type models: list(string)
        
        @return: the parameter sets
        @rtype: list[Star()]
        
        '''
        
        star_grid = Star.makeStars(models=models,\
                                   id_type='pacs' in models[0].lower() \
                                                and 'PACS' \
                                                or 'GASTRoNOoM',\
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
        pfns = []
        for i,star in enumerate(star_grid):
            if star['LAST_GASTRONOOM_MODEL']:    
                rad = star.getGasRad(unit='rstar')
                vel = star.getGasVelocity()
                vel = vel/10.**5
                avgdrift = star.getAverageDrift()/10.**5
                pfn = os.path.join(self.pplot,'velocity_%s_%i'\
                                   %(star['LAST_GASTRONOOM_MODEL'],i))
                plot_title = '%s %s: Velocity Profile for Model %i (%s)'\
                             %(self.plot_id.replace('_','\_'),\
                               self.star_name_plots,i,\
                               star['LAST_GASTRONOOM_MODEL'].replace('_','\_'))
                pfns.append(Plotting2.plotCols(x=rad,\
                                        y=[vel,avgdrift],cfg=cfg,\
                                        filename=pfn,\
                                        xaxis='R (R$_*$)',\
                                        yaxis=r'$v$ (km s$^{-1}$)',\
                                        plot_title=plot_title,\
                                        key_location=(0.0,0.0),\
                                        keytags=['Gas Velocity',\
                                                 'Grain-size Weighted Drift'],\
                                        xlogscale=1))
        if plot_filenames and plot_filenames[0][-4:] == '.pdf':    
            new_filename = os.path.join(self.pplot,'velocity_profiles.pdf')
            DataIO.joinPdf(old=sorted(pfns),new=new_filename)
            print '** Plots can be found at:'
            print new_filename
            print '***********************************' 
        elif pfns:
            print '** Plots can be found at:'
            print '\n'.join(pfns)
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
            radii = [star.getGasRad(unit='rstar') for star in valid_sg]
            temp = [star.getGasTemperature() for star in valid_sg]

            if temp:    
                keytags = star_grid[0].has_key('LAST_PACS_MODEL') \
                            and ['%s,    %s,    Mdot = %.2e'\
                                 %(star['LAST_GASTRONOOM_MODEL']\
                                       .replace('_','\_'),\
                                   str(star['LAST_PACS_MODEL'])\
                                       .replace('_','\_'),\
                                   float(star['MDOT_GAS'])) 
                                 for star in star_grid] \
                            or [star['LAST_GASTRONOOM_MODEL'].replace('_','\_')
                                for star in valid_sg]
                pfn = os.path.join(self.pplot,'temperature_profiles')
                pfn_rstar=os.path.join(self.pplot,'temperature_profiles_rstar')
                pfn_rstar = Plotting2.plotCols(x=radii_rstar,y=temps,\
                            cfg=cfg,xaxis='R (R$_*$)',\
                            filename=pfn_rstar,yaxis='T (K)',\
                            xlogscale=1,ylogscale=1,keytags=keytags)
                keys_cm = ['Model %i'%(i+1)
                           for i in xrange(len(star_grid))]
                pfn = Plotting2.plotCols(x=radii,y=temps,cfg=cfg,\
                        filename=pfn,xaxis='$r$ (cm)',\
                        yaxis='$T_\mathrm{g}$ (K)',\
                        figsize=(12.5,8),fontsize_ticklabels=26,\
                        key_location=(0.05,0.05),xlogscale=1,ylogscale=1,\
                        keytags=keys_cm,fontsize_axis=26,fontsize_key=26)
                print '** Plots can be found at:'
                print pfn
                print pfn_rstar
                print '***********************************'
            else:
                print '** No GASTRoNOoM models were calculated successfully.'+\
                      'No temperature profiles can be plotted.'
                print '***********************************'



    def plotTransitions(self,star_grid,cfg='',no_data=0,vg_factor=3,\
                        telescope_label=1,sort_freq=0,sort_molec=0,\
                        no_models=0,limited_axis_labels=0,date_tag=1,\
                        n_max_models=10,fn_suffix='',mfiltered=0,\
                        plot_intrinsic=0,plot_unresolved=0,cont_subtract=1):
        
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
        @keyword no_models: Only show data for the resolved lines.
                                  
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
        @keyword plot_intrinsic: Plot the intrinsic profiles instead of the 
                                 beam-convolved profiles (for instance when
                                 comparing GASTRoNOoM models to LIME models)
                                 
                                 (default: 0)
        @type plot_intrinsic: bool
        @keyword plot_unresolved: Plot intrinsic line profiles as well, for 
                                  unresolved data. The data themselves are not 
                                  added. By default, this is off as the line 
                                  profiles do not give you a lot information 
                                  before convolution with the wavelength 
                                  resolution.
        
                                  (default: 0)
        @type plot_unresolved: bool
        
        @keyword cont_subtract: Subtract the continuum from the sphinx line 
                                profile        
        
                                (default: 1)
        @type cont_subtract: bool
        
        """
        
        print '***********************************'
        print '** Creating Transition plots.'
        #- Default dimension is (4,3), but can be adapted in cfg
        cfg_dict = Plotting2.readCfg(cfg)
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
        if cfg_dict.has_key('plot_unresolved'):
            plot_unresolved = int(cfg_dict['plot_unresolved'])
        if cfg_dict.has_key('mfiltered'):
            mfiltered = int(cfg_dict['mfiltered'])
        if cfg_dict.has_key('cont_subtract'):
            cont_subtract = int(cfg_dict['cont_subtract'])
        if fn_suffix: 
            filename = cfg_dict.get('filename',None)
            if filename <> None: filename = filename + '_%s'%fn_suffix
            cfg_dict['filename'] = filename
        if cfg_dict.has_key('keytags'):
            keytags = cfg_dict['keytags']
            unreso_keytags = keytags
        else:
            keytags = [',\\ '.join(set(\
                        [trans.getModelId() != '' \
                            and str(trans.getModelId())\
                                   .replace('_','\_')
                            or str(None)
                         for trans in star['GAS_LINES']]))
                       for star in star_grid]
            unreso_keytags = list(keytags)
            if not no_data : keytags.append('Data')
        #-- Check how many resolved transitions there are (whether they have 
        #   data or not.
        trans_list = Transition.extractTransFromStars(star_grid,sort_freq,\
                                                      sort_molec,\
                                                      dtype='resolved')
                              
        #-- If plot_unresolved is requested, plot all modeled unresolved lines 
        #   (intrinsic)
        if plot_unresolved:
            unreso_list  = Transition.extractTransFromStars(star_grid,\
                                                            sort_freq,\
                                                            sort_molec,\
                                                            dtype='unresolved')
        else:
            unreso_list = []
            
        def createTilePlots(trans_list,x_dim,y_dim,no_data,intrinsic,\
                            vg_factor,keytags,telescope_label,no_models,cfg,\
                            star_grid,limited_axis_labels,date_tag,indexi,\
                            indexf,mfiltered,cont_subtract):
            
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
            @type cfg: string/dict
            @param vg_factor: The factor with which the terminal velocity is 
                              multiplied. This determines the xrange of the 
                              plots
            @type vg_factor: float
            @param telescope_label: Include a label showing the telescope name
            @type telescope_label: bool
            @param no_models: Only show data for the resolved lines.
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
            @param cont_subtract: Subtract the continuum value outside the line
                                  from the whole line profile. 
            @type cont_subtract: bool
            
            @return: The data list with dictionaries for every tile is returned
            @rtype: list[dict]

            '''
            
            if cfg.has_key('filename'):
                filename = cfg.pop('filename')
            else:
                seg = '%slps_%smodels_%ito%i'\
                      %(intrinsic and 'intrinsic_' or '',\
                        no_data and 'nodata_' or '',\
                        indexi,indexf)
                filename = os.path.join(self.pplot,seg)
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
                    #-- Read the data profiles, and set them for all models.
                    if not no_data:
                        current_trans.readData()
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
                                     and list(trans.sphinx.getLPIntrinsic(cont_subtract=cont_subtract)*10**(23))\
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
                                    ddict['y'].append(trans.sphinx.getLPTmb(cont_subtract=cont_subtract))
                                else:
                                    ddict['x'].append([])
                                    ddict['y'].append([])
                    else:
                        ddict['x'], ddict['y'] = [], []
                    
                    #- Add data, but only if the data filename is known. This 
                    #- will be Tmb, in K. In case of intrinsic==1, you dont 
                    #- even want to check this.
                    if current_trans.lpdata and not no_data:
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
                    if not no_data and date_tag:
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
                cfg['filename'] = filename + '_trl%i'%(i)
                plot_filenames.append(Plotting2.plotTiles(extension='pdf',\
                     data=data,filename=filename,keytags=keytags,\
                     xaxis=r'$v$ (km s$^{-1}$)',fontsize_axis=16,cfg=cfg,\
                     yaxis=intrinsic \
                            and r'$F_\nu$ (Jy)' \
                            or '$T_\mathrm{mb}$ (K)',\
                     fontsize_ticklabels=16,dimensions=(x_dim,y_dim),\
                     fontsize_label=20,linewidth=2))
            if missing_trans:
                print 'WARNING! %i requested transitions were '%missing_trans+\
                      'not found for a Star(). Within one CC session, this '+\
                      'should not be the case!'
            print '** %sine profile plots %scan be found at:'\
                    %(intrinsic and 'Intrinsic l' or 'L',\
                        no_data and 'without data ' or '')
            if plot_filenames and plot_filenames[0][-4:] == '.pdf':    
                DataIO.joinPdf(old=plot_filenames,new=filename+'.pdf')
                print filename+'.pdf'
            else:
                print '\n'.join(plot_filenames)
            print '***********************************' 
            
    
        if trans_list:
            j = 0
            while j < len(star_grid):
                if plot_intrinsic == 1: no_data = 1
                i = 0 
                subgrid = []
                subkeys = []
                while i < n_max_models and i+j < len(star_grid):
                    subgrid.append(star_grid[i+j])
                    if keytags: subkeys.append(keytags[i+j])
                    i += 1
                if keytags and not no_data: subkeys.append(keytags[-1])
                #- Copying the list so that the destructive loop does not mess
                #- up multiple tile plot runs if n_models > n_max_models
                createTilePlots(trans_list=list(trans_list),\
                                vg_factor=vg_factor,\
                                no_data=no_data,cfg=cfg_dict,star_grid=subgrid,\
                                x_dim=x_dim,y_dim=y_dim,keytags=subkeys,\
                                intrinsic=plot_intrinsic,no_models=no_models,\
                                telescope_label=telescope_label,\
                                limited_axis_labels=limited_axis_labels,\
                                date_tag=date_tag,indexi=j,indexf=j+i-1,\
                                mfiltered=mfiltered,\
                                cont_subtract=cont_subtract)
                j += i
        if unreso_list: 
            j = 0
            while j < len(star_grid):
                i = 0 
                subgrid = []
                subkeys = []
                while i < n_max_models and i+j < len(star_grid):
                    subgrid.append(star_grid[i+j])
                    if keytags: subkeys.append(unreso_keytags[i+j])
                    i += 1
                createTilePlots(trans_list=list(unreso_list),\
                                vg_factor=vg_factor,\
                                no_data=1,cfg=cfg_dict,star_grid=subgrid,\
                                intrinsic=1,keytags=subkeys,\
                                x_dim=x_dim,y_dim=y_dim,date_tag=date_tag,\
                                telescope_label=telescope_label,no_models=0,\
                                limited_axis_labels=limited_axis_labels,\
                                indexi=j,indexf=j+i-1,mfiltered=mfiltered,\
                                cont_subtract=cont_subtract)
                j += i            
                


    def createLineLabelsFromLineLists(self,star,xmin,xmax,xunit='micron',\
                                      fn_trans_marker='',instrument='PACS'):
        
        '''
        Create a list of line labels for all molecules and transitions 
        in the molecular linelists requested.
        
        This is used for spectroscopic databases only! Such as JPL, CDMS, LAMDA
        
        @param star: The parameter set
        @type star: Star()
        @param xmin: minimum wavelength
        @type xmin: float
        @param xmax: maximum wavelength
        @type xmax: float
        
        @keyword xunit: The unit of the xmax/xmin
                         
                         (default: micron)
        @type xunit: string
        @keyword fn_trans_marker: A file that includes TRANSITION definitions.
                                  These transitions will be marked up in the 
                                  plot. For instance, when indicating a subset 
                                  of transitions for one reason or another.
                                  The line type can be set for this specific 
                                  subset, differently from other lines and 
                                  regardless of the molecule. In combination
                                  with a doubly defined line label (through the
                                  star_grid['GAS_LINES']), lines can be marked
                                  up.
                                  
                                  (default: '')
        @type fn_trans_marker: string
        @keyword instrument: The instrument object for which the line labels 
                             are created. Used to retrieve the v_lsr. Either
                             'PACS' or 'SPIRE'
                             
                             (default: 'PACS')
        @type instrument: str
        
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
        
        linelists = []
        for molecule in star['LL_GAS_LIST']:
            if not 'p1H' in molecule.molecule:
                ll = LineList.LineList(molecule=molecule,x_min=xmin,\
                                       x_unit=xunit,cdms=cdms,jpl=jpl,\
                                       lamda=lamda,x_max=xmax,\
                                       min_strength=min_strength,\
                                       max_exc=max_exc,include_extra=1)
                linelists.append(ll)
        lls = self.createLineLabels(linelists=linelists,\
                                    fn_trans_marker=fn_trans_marker,
                                    instrument=instrument)
        return lls     
    
    
    
    def plotLineLists(self,star_grid,include_sphinx=1,cfg='',\
                      fn_trans_marker='',instrument='PACS'):
        
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
        @keyword fn_trans_marker: A file that includes TRANSITION definitions.
                                  These transitions will be marked up in the 
                                  plot. For instance, when indicating a subset 
                                  of transitions for one reason or another.
                                  The line type can be set for this specific 
                                  subset, differently from other lines and 
                                  regardless of the molecule. In combination
                                  with a doubly defined line label (through the
                                  star_grid['GAS_LINES']), lines can be marked
                                  up.
                                  
                                  (default: '')
        @type fn_trans_marker: string
        @keyword instrument: The unresolved-data instrument for which the line
                             lists are to be plotted. Either 'PACS' or 'SPIRE'
                            
                             (default: 'PACS')
        @type instrument: str
        
        '''
        
        print '***********************************'
        print '** Starting to plot line identifications for %s from databases.'\
              %self.star_name
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('instrument'):
            instrument = cfg_dict['instrument']
        instrument = instrument.upper()
        if instrument == 'SPIRE': 
            instr = self.spire
        else:
            instrument = 'PACS'
            instr = self.pacs
        
        if instr is None: 
            print '** No %s_PATH given. Cannot plot line lists '%instrument + \
                  'without data information. Aborting...'
            return
        if cfg_dict.has_key('fn_trans_marker'):
            fn_trans_marker = cfg_dict['fn_trans_marker']
        if cfg_dict.has_key('include_sphinx'):
            include_sphinx = bool(cfg_dict['include_sphinx'])
        xmins = [min(wave_list) for wave_list in instr.data_wave_list]
        xmaxs = [max(wave_list) for wave_list in instr.data_wave_list]
        lls = self.createLineLabelsFromLineLists(star=star_grid[0],\
                                                 xmin=min(xmins),\
                                                 xmax=max(xmaxs),\
                                                 fn_trans_marker=\
                                                     fn_trans_marker,\
                                                 instrument=instrument)
        plot_filenames = []
        DataIO.testFolderExistence(os.path.join(self.pplot,'LineLists'))
        if include_sphinx:
            if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) \
                            == set([0]) \
                    or set([s['LAST_GASTRONOOM_MODEL'] for s in star_grid]) \
                            == set(['']): 
                include_sphinx = 0
            else: 
                if instrument == 'PACS':
                    self.setSphinxPacs(star_grid)
        for i_file,(wave,flux,filename,xmin,xmax) in enumerate(\
                    zip(instr.data_wave_list,instr.data_flux_list,\
                        instr.data_filenames,xmins,xmaxs)):
            if include_sphinx and instrument == 'PACS':
                sphinx_flux = [sphinx 
                               for sphinx in self.sphinx_flux_list[i_file] 
                               if list(sphinx)]
            elif include_sphinx and instrument == 'SPIRE':
                sphinx_flux = [instr.getSphinxConvolution(star,filename)[1]
                               for star in star_grid]
            else:
                sphinx_flux = []
            seg = 'line_id_'+os.path.split(filename)[1].replace('.dat','')
            plot_filename = os.path.join(self.pplot,'LineLists',seg)
            keytags = ['%s %s'%(instrument,filename.replace('_','\_'))] + \
                      ['Model %i: %s'\
                       %(i+1,instrument=='PACS' \
                             and str(star['LAST_PACS_MODEL']).replace('_','\_')\
                             or star['LAST_GASTRONOOM_MODEL'].replace('_','\_')) 
                       for i,star in enumerate(star_grid) 
                       if star['LAST_GASTRONOOM_MODEL'] and include_sphinx]
            plot_filenames.append(Plotting2.plotCols(\
                    x=[wave]*(len(sphinx_flux)+1),y=[flux]+sphinx_flux,\
                    cfg=cfg_dict,filename=plot_filename,keytags=keytags,\
                    plot_title=self.star_name_plots,histoplot=[0],\
                    number_subplots=3,line_labels=lls,\
                    line_label_color=1,line_label_lines=1,\
                    line_label_spectrum=1))
        seg = 'line_id_%s_%s.pdf'%(instrument.lower(),self.star_name)
        new_filename = os.path.join(self.pplot,'LineLists',seg)
        DataIO.joinPdf(old=sorted(plot_filenames),new=new_filename)
        print '** Plots can be found at:'
        print new_filename
        print '***********************************'

                                                
   
    def plotAbundanceProfiles(self,star_grid=[],models=[],cfg='',\
                              per_molecule=0,unit='cm'):  
        
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
        @keyword per_molecule: Plot one molecule for all models in one figure.
        
                               (default: 0)
        @type per_molecule: bool
        @keyword unit: The radial unit. Can be 'cm', 'au', 'm' or 'rstar'
        
                       (default: cm)
        @type unit: str
        
        '''
        
        print '***********************************'
        print '** Plotting Abundance Profiles'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        pfns = []
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('per_molecule'):
            per_molecule = cfg_dict['per_molecule']
        if cfg_dict.has_key('unit'):
            unit = cfg_dict['unit'].lower()
        
        #-- Some general plot settings
        extra_pars = dict()
        extra_pars['ymin'] = 1e-9
        extra_pars['ymax'] = 1e-3
        extra_pars['ylogscale'] = 1 
        extra_pars['xlogscale'] = 1
        extra_pars['figsize'] = (12.5,8.5)
        
        if unit == 'cm': xaxis = '$r\ \mathrm{(cm)}$'
        elif unit =='au': xaxis = '$r\ \mathrm{(AU)}$'
        elif unit == 'm': xaxis = '$r\ \mathrm{(m)}$'
        else: xaxis = '$r\ \mathrm{(R}_\star\mathrm{)}$'
        extra_pars['xaxis'] = xaxis
        
        #-- Dict to keep track of all data
        ddata = dict()
        for istar,star in enumerate(star_grid):
            if not star['LAST_GASTRONOOM_MODEL']: continue
            ddata[istar] = dict()
            for molec in star['GAS_LIST']: 
                mid = molec.getModelId()
                if not mid: continue
                ddata[istar][molec.molecule] = dict()
                rad = star.getGasRad(unit=unit,ftype='1',mstr=molec.molecule,
                                     modelid=mid)
                nh2 = star.getGasNumberDensity(ftype='1',mstr=molec.molecule,
                                               modelid=mid)
                nmol = star.getGasNumberDensity(molecule=1,ftype='1',\
                                                mstr=molec.molecule,\
                                                modelid=mid)
                ddata[istar][molec.molecule]['rad'] = rad
                ddata[istar][molec.molecule]['nh2'] = nh2
                ddata[istar][molec.molecule]['nmol'] = nmol
                if molec.set_keyword_change_abundance:
                    cff = DataIO.readCols(molec.change_fraction_filename)
                    rfrac,frac = cff[0],cff[1]
                    rfrac = rfrac
                    frac_interpol = interp1d(rfrac,frac)(rad)
                    #- if error happens, catch and print out warning, plus run
                    #- interpolation again with bounds_error=False, fill_value=frac[-1]
                    #- if bounds_error=False and a warning is printed by scipy, 
                    #- then no need to catch error first
                    #- frac_interpol = array(Interpol.doInterpol(x_in=rfrac,\
                    #-            y_in=frac,gridsx=[rad])[0][0])
                else:
                    frac_interpol = 1
                #-- GASTRoNOoM output already takes into account enhance_abundance_factor
                #   abun_factor only takes into account isotope ratios and OPR
                abun = nmol/nh2*frac_interpol*molec.abun_factor     
                ddata[istar][molec.molecule]['abun'] = abun
                ddata[istar][molec.molecule]['key'] = molec.molecule_plot
                ddata[istar][molec.molecule]['id'] = mid
                
            if not per_molecule:
                #-- Collect all data
                radii = [dmol['rad'] for molec,dmol in ddata[istar].items()]
                abuns = [dmol['abun'] for molec,dmol in ddata[istar].items()]
                keytags = [dmol['key'] for molec,dmol in ddata[istar].items()]
                ids = [dmol['id'] for molec,dmol in ddata[istar].items()]
                
                #-- Add additional information if requested
                if star.has_key('R_DES_H2O'):
                    radii.extend([array([star['R_DES_H2O'],star['R_DES_H2O']])])
                    abuns.extend([[1e-2,1e-9]])
                    keytags.append('Condensation radius H$_2$O ice')
                    lt.append('--k')
                if star['R_OH1612']:
                    radii.extend([array([star['R_OH1612'],star['R_OH1612']])])
                    abuns.extend([[1e-2,1e-9]])
                    keytags.append('Location OH maser')
                    lt.append('-k')
                
                #-- Set the yaxis tag
                yaxis = '$n_\mathrm{molec}/n_{\mathrm{H}_2}$'
                
                #-- Make filename
                pfn = os.path.join(self.pplot,'abundance_profiles_%s'\
                                              %'_'.join(list(set(ids))))
                pfns.append(Plotting2.plotCols(x=radii,y=abuns,cfg=cfg,\
                                               filename=pfn,keytags=keytags,\
                                               yaxis=yaxis,**extra_pars))
        
        if per_molecule:
            #-- Collect all data
            molecs = list(set([molec for istar in ddata.keys()
                                     for molec in ddata[istar].keys()]))
            for molec in molecs: 
                #-- Collect data
                radii = [dmol['rad']
                         for istar in ddata.keys()
                         for imolec,dmol in ddata[istar].items()
                         if molec == imolec]
                abuns = [dmol['abun']
                         for istar in ddata.keys()
                         for imolec,dmol in ddata[istar].items()
                         if molec == imolec]
                keytags = [dmol['id'].replace('_','\_') 
                           for istar in ddata.keys()
                           for imolec,dmol in ddata[istar].items()
                           if molec == imolec]

                #-- Set the y axis tag
                strmolec = ddata[0][molec]['key']
                yaxis = '$n_\mathrm{%s}/n_{\mathrm{H}_2}$'%strmolec

                #-- Make filename
                pfn = os.path.join(self.pplot,'abundance_profiles_%s'%molec)
                pfns.append(Plotting2.plotCols(x=radii,y=abuns,yaxis=yaxis,\
                                               filename=pfn,keytags=keytags,\
                                               cfg=cfg,**extra_pars))  
                
        if not per_molecule and pfns and pfns[0][-4:] == '.pdf':    
            newfn = os.path.join(self.pplot,'abundance_profiles.pdf')
            DataIO.joinPdf(old=pfns,new=newfn)
            print '** Plots can be found at:'
            print newfn
            print '***********************************'
        else:
            print '** Plots can be found at:'
            print '\n'.join(pfns)
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
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('do_sort'):
             do_sort = int(cfg_dict['do_sort'])
        if cfg_dict.has_key('normalized'):
             normalized = int(cfg_dict['normalized'])
        if cfg_dict.has_key('include_velocity'):
             include_velocity = int(cfg_dict['include_velocity'])
        DataIO.testFolderExistence(os.path.join(self.pplot,'LCs'))
        normalized = int(normalized)
        for i,star in enumerate(star_grid):
            extra_pars = dict()
            if do_sort:
                transitions = sorted([trans 
                                      for trans in star['GAS_LINES'] 
                                      if trans.getModelId()],\
                                     key=lambda x:x.wavelength)
            else:
                transitions = star['GAS_LINES']
            if include_velocity:
                rad = star.getGasRad(unit='rstar')
                vel = star.getGasVelocity()
                vel = vel/10.**5
                extra_pars['twiny_x'] = [rad]
                extra_pars['twiny_y'] = [vel]
                extra_pars['twiny_keytags'] = [r'$v_\mathrm{g}$']
                extra_pars['twinyaxis'] = r'$v_\mathrm{g}$ (km s$^{-1}$)' 
            
            [trans.readSphinx() for trans in transitions]
            radii = [trans.sphinx.getImpact() for trans in transitions]
            linecontribs =  [normalized \
                                and list(trans.sphinx.getNormalizedIntensity())\
                                or list(trans.sphinx.getWeightedIntensity())
                             for trans in transitions]
            
            pfn = os.path.join(self.pplot,'LCs','linecontrib_%s_%i'\
                               %(star['LAST_GASTRONOOM_MODEL'],i))
            extra_pars['filename'] = pfn
            extra_pars['keytags'] = ['$\mathrm{%s}:$ %s'\
                                      %(trans.molecule.molecule_plot,\
                                        trans.makeLabel())
                                     for trans in transitions]
            extra_pars['key_location'] = 'lower right'
            extra_pars['ymin'] = normalized and -0.01 or None
            extra_pars['ymax'] = normalized and 1.02 or None
            extra_pars['xmin'] = 1
            extra_pars['xaxis'] = '$p\ \mathrm{(R}_\star\mathrm{)}$'
            extra_pars['yaxis'] = '$I(p)\ pdp$'
            extra_pars['linewidth'] = 3
            extra_pars['xlogscale'] = 1
            
            pfn = Plotting2.plotCols(x=radii,y=linecontribs,cfg=cfg_dict,\
                                     **extra_pars)
            print '** Plot can be found at:'
            print pfn
            print '***********************************'                            
        

        
    def setSphinxPacs(self,star_grid,refresh_sphinx_flux=0):
        
        ''' 
        Prepare Sphinx output in Pacs format (including convolution).
        
        @param star_grid: Parameter sets. If empty list, no sphinx models are 
                          set, but an empty list is set for each datafile.
        @type star_grid: list[Star()]
        
        @keyword refresh_sphinx_flux: redo the sphinx flux list by pulling from
                                      db, regardless of whether it's already
                                      been done or not.
                                      
                                      (default: 0)
        @type refresh_sphinx_flux: bool
        
        '''
        
        if not self.sphinx_flux_list or refresh_sphinx_flux:    
            self.pacs.prepareSphinx(star_grid)
            #- The sphinx convolved models always have the same wavelength 
            #- list as their respective data files
            self.sphinx_flux_list = [[self.pacs.getSphinxConvolution(star,fn)[1]
                                      for star in star_grid]
                                     for fn in self.pacs.data_filenames]
            



    def createLineLabels(self,star_grid=[],linelists=[],fn_trans_marker='',\
                         unit='micron',mark_undetected=0,instrument='PACS'):

        '''
        Create line labels for all transitions in Star() objects or in
        LineList() objects or in a TRANSITION definition file. Priority:
        star_grid > linelists. fn_trans_marker is always added in addition.

        @keyword star_grid: The Star() models.

                            (default: [])
        @type star_grid: list[Star()]
        @keyword linelists: The LineList() objects.

                            (default: [])
        @type linelists: list[LineList()]
        @keyword fn_trans_marker: A file that includes TRANSITION definitions.
                                  These transitions will be marked up in the
                                  plot. For instance, when indicating a subset
                                  of transitions for one reason or another.
                                  The line type can be set for this specific
                                  subset, differently from other lines and
                                  regardless of the molecule. In combination
                                  with a doubly defined line label (through
                                  star_grid/linelists), lines can be marked
                                  up.

                                  (default: '')
        @type fn_trans_marker: string
        @keyword mark_undetected: Mark the undetected transitions in the same
                                  way extra marked transitions would be marked
                                  by fn_trans_marker.

                                  (default: 0)
        @type mark_undetected: bool
        @keyword unit: The unit of the location number. Can be 'micron' or
                       'cm-1'.

                       (default: 'micron')
        @type unit: string
        @keyword instrument: The instrument object for which the line labels 
                             are created. Used to retrieve the v_lsr. Either
                             'PACS' or 'SPIRE'
                             
                             (default: 'SPIRE')
        @type instrument: str

        @return: a sorted list(set) of line labels
        @rtype: list[string]

        '''

        instrument = instrument.upper()
        if instrument == 'PACS':
          vlsr = self.pacs.vlsr
        elif instrument == 'SPIRE':
          vlsr = self.spire.vlsr
          
        if star_grid:
            alltrans = Transition.extractTransFromStars(star_grid,\
                                                        dtype=instrument)
        elif linelists:
            alltrans = [t   for ll in linelists
                            for t in ll.makeTransitions()]
        else:
            alltrans = []

        lls = [('%s %s'%(t.molecule.molecule,t.makeLabel()),\
                t.wavelength*10**4*1./(1-vlsr/t.c),\
                t.molecule.molecule_index,\
                t.vup>0)
               for t in alltrans]

        used_indices = list(set([ll[-2] for ll in lls]))
        if fn_trans_marker:
            def_molecs = dict([(ll.molecule.molecule,ll.molecule) 
                               for ll in linelists])
            if star_grid: star = star_grid[0]
            else: star = None
            trl = DataIO.readDict(fn_trans_marker,multi_keys=['TRANSITION'])
            n_entry = len(trl['TRANSITION'][0].split())
            trl_sorted = DataIO.checkEntryInfo(trl['TRANSITION'],n_entry,\
                                               'TRANSITION')
            etrans = [makeTransition(trans=t,def_molecs=def_molecs,star=star) 
                      for t in trl_sorted]
            this_index = max(used_indices)+1
            used_indices = used_indices + [this_index]
            ells = [('%s %s'%(t.molecule.molecule,t.makeLabel()),\
                    t.wavelength*10**4*1./(1-vlsr/t.c),\
                    this_index,\
                    t.vup>0)
                   for t in etrans]
            lls = lls + ells

        if mark_undetected:
            this_index = max(used_indices)+1
            extra_trans = [t 
                           for t in alltrans 
                           if t.getIntIntUnresolved()[0] is None]
            ells = [('%s %s'%(t.molecule.molecule,t.makeLabel()),\
                     t.wavelength*10**4*1./(1-vlsr/t.c),\
                     this_index,\
                     t.vup>0)
                    for t in extra_trans]
            lls = lls + ells

        if unit == 'cm-1':
            lls = [(l,1./w*10**4,i,vib) for l,w,i,vib in lls]
        lls = sorted(lls,key=operator.itemgetter(1))
        return lls 
    
        
        
    def plotPacsLineScans(self,star_grid=[],models=[],exclude_data=0,cfg='',\
                          cont_subtracted=1,fn_trans_marker='',fn_plt='',\
                          dimensions=(5,2),fn_add_star=1,mark_undetected=0,\
                          remove_axis_titles=1,include_band=1):
         
        '''
        Plot PACS line scans.
        
        Data can be in- or excluded, as can models.
        
        Both continuum subtracted data as well as the original spectra can be 
        plotted.
        
        @keyword star_grid: star models for which PACS data will be fetched. 
                               
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: list of pacs_ids or gastronoom model ids. If neither 
                         this or star_grid are defined, only data are plotted.
                         If star_grid is defined, this keyword is ignored.
                         (default: [])
        @type models: list[strings]
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                       
                      (default: '')
        @type cfg: string         
        @keyword fn_add_star: Add the star name to the requested plot filename.
        
                              (default: 1)
        @type fn_add_star: bool
        @keyword fn_plt: A plot filename for the tiled plot.
                         
                         (default: '')
        @type fn_plt: string
        @keyword fn_trans_marker: A file that includes TRANSITION definitions.
                                  These transitions will be marked up in the 
                                  plot. For instance, when indicating a subset 
                                  of transitions for one reason or another.
                                  The line type can be set for this specific 
                                  subset, differently from other lines and 
                                  regardless of the molecule. In combination
                                  with a doubly defined line label (through the
                                  star_grid['GAS_LINES']), lines can be marked
                                  up.
                                  
                                  (default: '')
        @type fn_trans_marker: string
        @keyword mark_undetected: Mark the undetected transitions in the same
                                  way extra marked transitions would be marked
                                  by fn_trans_marker. 
                                  
                                  (default: 0)
        @type mark_undetected: bool
        @keyword remove_axis_titles: Remove axis titles in between different 
                                     tiles in the plot. Only keeps the ones on
                                     the left and the bottom of the full plot.
                                     
                                     (default: 1)
        @type remove_axis_titles: bool
        @keyword include_band: Include a label that names the band. 
        
                               (default: 1)
        @type include_band: bool
        @keyword exclude_data: if enabled only the sphinx mdels are plotted.
         
                               (default: 0)
        @type exclude_data: bool
        @keyword cont_subtracted: Plot the continuum subtracted data.
        
                                  (default: 1)
        @type cont_subtracted: bool
        @keyword dimensions: The number of tiles in the x and y direction is 
                             given: (x-dim,y-dim)
                             
                             (default: (5,2))
        @type dimensions: tuple(int,int)
        
        '''
         
        print '***********************************'
        print '** Plotting line scans.'
        if self.pacs is None: 
            print '** No PATH_PACS given. Cannot plot PACS spectra without '+\
                  'data information. Aborting...'
            return
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
                
        self.setSphinxPacs(star_grid)
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        if cfg_dict.has_key('exclude_data'):
            exclude_data = bool(cfg_dict['exclude_data'])
        if cfg_dict.has_key('fn_trans_marker'):
            fn_trans_marker = cfg_dict['fn_trans_marker']
        if cfg_dict.has_key('cont_subtracted'):
            cont_subtracted = cfg_dict['cont_subtracted']
        if cfg_dict.has_key('fn_add_star'):
            fn_add_star = cfg_dict['fn_add_star']
        if cfg_dict.has_key('mark_undetected'):
            mark_undetected = cfg_dict['mark_undetected']
        if cfg_dict.has_key('remove_axis_titles'):
            remove_axis_titles = cfg_dict['remove_axis_titles']
        if cfg_dict.has_key('dimensions'):
            dimensions = (int(cfg_dict['dimensions'][0]),\
                          int(cfg_dict['dimensions'][1]))
        
        if not star_grid: 
            exclude_data = 0

        lls = self.createLineLabels(star_grid=star_grid,\
                                    fn_trans_marker=fn_trans_marker,\
                                    mark_undetected=mark_undetected,\
                                    instrument='PACS')
        tiles = []            
        print '** Plotting now...'
        for idd,(wave,flux,flux_ori,sphinx_flux,filename,ordername) in \
              enumerate(sorted(zip(self.pacs.data_wave_list,\
                                   self.pacs.data_flux_list,\
                                   self.pacs.data_original_list,\
                                   self.sphinx_flux_list,\
                                   self.pacs.data_filenames,\
                                   self.pacs.data_ordernames),\
                               key=lambda x: x[0][0])):
            ddict = dict()
            ddict['x'] = exclude_data \
                            and [wave]*(len(sphinx_flux)) \
                            or [wave]*(len(sphinx_flux)+1)
            d_yvals = cont_subtracted and [flux] or [flux_ori]
            ddict['y'] = exclude_data and sphinx_flux or d_yvals+sphinx_flux
            if include_band:
                ddict['labels'] = [(ordername,0.08,0.80)]
            ddict['xmin'] = wave[0]
            ddict['xmax'] = wave[-1]
            ddict['ymin'] = 0.95*min([min(ff) for ff in ddict['y'] if ff.size])
            ddict['ymax'] = 1.2*max([max(ff) for ff in ddict['y'] if ff.size])
            ddict['histoplot'] = (not exclude_data) and [0] or []
            ddict['line_labels'] = lls
                                    #[(label,wl,i) 
                                    #for label,wl,i in lls
                                    #if wl <= wave[-1] and wl >= wave[0]]
            if remove_axis_titles: 
                n_tiles = len(self.pacs.data_filenames)
                ddict['xaxis'] = idd in range(n_tiles-dimensions[0],n_tiles)\
                                    and r'$\lambda$ ($\mu$m)' or ''
                ddict['yaxis'] = idd%dimensions[0] == 0 \
                                    and r'$F_\nu$ (Jy)' or ''
            tiles.append(ddict)
        if not fn_plt:
            fn_plt = os.path.join(self.pplot,'PACS_results','PACS_linescans')
        if fn_add_star:
            fn_plt = '_'.join([fn_plt,self.star_name])
        #tiles = sorted(tiles,key = lambda x: x['x'][0][0])
        pfn = Plotting2.plotTiles(filename=fn_plt,data=tiles,cfg=cfg_dict,\
                                  line_label_color=1,fontsize_label=15,\
                                  line_label_lines=1,dimensions=dimensions)
        print '** Your plot can be found at:'
        print pfn
        print '***********************************'

       

    def plotPacs(self,star_grid=[],models=[],exclude_data=0,fn_plt='',cfg='',\
                 fn_trans_marker='',include_band=1,number_subplots=3,\
                 mark_undetected=0,fn_add_star=1):
        
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
                                  The line type can be set for this specific 
                                  subset, differently from other lines and 
                                  regardless of the molecule. In combination
                                  with a doubly defined line label (through the
                                  star_grid['GAS_LINES']), lines can be marked
                                  up.
                                  
                                  (default: '')
        @type fn_trans_marker: string
        @keyword fn_add_star: Add the star name to the requested plot filename.
                              Only relevant if fn_plt is given.
                              
                              (default: 1)
        @type fn_add_star: bool
        @keyword mark_undetected: Mark the undetected transitions in the same
                                  way extra marked transitions would be marked
                                  by fn_trans_marker. 
                                  
                                  (default: 0)
        @type mark_undetected: bool
        @keyword include_band: Include a name tag for the band order in 
                                    the plot.
                                    
                                    (default: 1)
        @type include_band: bool
        @keyword line_label_dashedvib: Use dashed lines for vibrational 
                                       transitions in line label lines. 
                             
                                       (default: 0)
        @type line_label_dashedvib: bool
        
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
        
        self.setSphinxPacs(star_grid)
        print '** Plotting now...'
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        if cfg_dict.has_key('exclude_data'):
            exclude_data = bool(cfg_dict['exclude_data'])
        if cfg_dict.has_key('fn_trans_marker'):
            fn_trans_marker = cfg_dict['fn_trans_marker']
        if cfg_dict.has_key('include_band'):
            include_band = bool(cfg_dict['include_band'])
        if cfg_dict.has_key('mark_undetected'):
            mark_undetected = cfg_dict['mark_undetected']
        if cfg_dict.has_key('fn_add_star'):
            fn_add_star = cfg_dict['fn_add_star']
        if cfg_dict.has_key('labels'):
            labels = bool(cfg_dict['labels'])
        else:
            labels = []
        
        lls = self.createLineLabels(star_grid=star_grid,\
                                    fn_trans_marker=fn_trans_marker,\
                                    mark_undetected=mark_undetected,\
                                    instrument='PACS')
        
        plot_filenames = []
        for wave,flux,sphinx_flux,filename,ordername in \
                    zip(self.pacs.data_wave_list,\
                        self.pacs.data_flux_list,\
                        self.sphinx_flux_list,\
                        self.pacs.data_filenames,\
                        self.pacs.data_ordernames):
            if fn_plt and not fn_add_star:
                fn_plt = os.path.splitext(fn_plt)[0]
                this_filename = '%s_%s'%(fn_plt,ordername)
            elif fn_plt and fn_add_star:
                fn_plt = os.path.splitext(fn_plt)[0]
                this_filename = '%s_%s_%s'%(fn_plt,self.star_name,ordername)
            else:    
                this_filename = os.path.join(self.pplot,'PACS_results',\
                                             os.path.split(filename)[1]\
                                            .replace('.dat',''))
            keytags = ['Model %i: %s'%(i+1,str(star['LAST_PACS_MODEL'])\
                                .replace('_','\_')) 
                       for i,star in enumerate(star_grid)]
            if exclude_data:
                x_list = [wave]*(len(sphinx_flux)) 
                y_list = sphinx_flux
            else:
                x_list = [wave]*(len(sphinx_flux)+1)
                y_list = [flux]+sphinx_flux
                keytags = ['PACS Spectrum'] + keytags
            if include_band:
                elabel = [(ordername,0.05,0.80)]
            else:
                elabel = []
            plot_filenames.append(Plotting2.plotCols(x=x_list,y=y_list,\
                    keytags=keytags,number_subplots=3,\
                    plot_title='%s: %s - %s'%(self.plot_id.replace('_','\_'),\
                    self.star_name_plots,ordername),cfg=cfg_dict,\
                    line_labels=lls,\
                    histoplot=not exclude_data and [0] or [],\
                    filename=this_filename,labels=labels+elabel,\
                    line_label_spectrum=1,line_label_color=1))
        if plot_filenames and plot_filenames[0][-4:] == '.pdf':
            if not fn_plt:
                newf = os.path.join(self.pplot,'PACS_results',\
                                    'PACS_spectrum.pdf')
            elif fn_plt and fn_add_star:
                newf = '%s_%s.pdf'%(fn_plt,self.star_name)
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
                         fn_plt='',fn_trans_marker='',cfg='',\
                         include_sphinx=None,exclude_data=0):
        
        '''
        Plot segments of spectra only.
        
        An inputfile gives the wavelength ranges, given by pacs_segments_path.
        
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
        @keyword exclude_data: The data are not included when True. 
                                
                               (default: 0) 
        @type exclude_data: bool
        @keyword fn_plt: A plot filename to which an index is added for each
                         subband.
                         
                         (default: '')
        @type fn_plt: string
        @keyword fn_trans_marker: A file that includes TRANSITION definitions.
                                  These transitions will be marked up in the 
                                  plot. For instance, when indicating a subset 
                                  of transitions for one reason or another.
                                  The line type can be set for this specific 
                                  subset, differently from other lines and 
                                  regardless of the molecule. In combination
                                  with a doubly defined line label (through the
                                  star_grid['GAS_LINES']), lines can be marked
                                  up.
                                  
                                  (default: '')
        @type fn_trans_marker: string
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('mode'):
            mode = cfg_dict['mode']
        if cfg_dict.has_key('pacs_segments_path'):
            pacs_segments_path = cfg_dict['pacs_segments_path']
        if cfg_dict.has_key('include_sphinx'):
            include_sphinx = cfg_dict['include_sphinx']
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        if cfg_dict.has_key('fn_trans_marker'):
            fn_trans_marker = cfg_dict['fn_trans_marker']
        if cfg_dict.has_key('exclude_data'):
            exclude_data = bool(cfg_dict['exclude_data'])
            
        if self.pacs is None:
            print 'No PACS data found for plotting PACS segments. Aborting...'
            return
        elif not pacs_segments_path:
            print 'No pacs_segments_path given. Pass in the cfg file or in ' +\
                  'the method call. Aborting...' 
            return
        else:
            self.setSphinxPacs(star_grid)

        if mode == 'll':
            xmins=[min(wave_list) for wave_list in self.pacs.data_wave_list]
            xmaxs=[max(wave_list) for wave_list in self.pacs.data_wave_list]
            lls = self.createLineLabelsFromLineLists(star=star_grid[0],\
                                                     xmin=min(xmins),\
                                                     xmax=max(xmaxs),\
                                                     fn_trans_marker=\
                                                            fn_trans_marker,\
                                                     instrument='PACS')
        elif mode == 'sphinx':
            lls = self.createLineLabels(star_grid=star_grid,\
                                        fn_trans_marker=fn_trans_marker,
                                        instrument='PACS')
        else:
            print 'Mode for plotting PACS segments not recognized. Aborting...'
            return
        
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
                    if fn_plt:
                        fn_plt = os.path.splitext(fn_plt)[0]
                        plot_filename = '%s_%s_%.2f_%.2f'%(fn_plt,\
                                                           mode,wmin,wmax)
                    else:    
                        folder = mode=='ll' and 'LineLists' or 'PACS_results'
                        dfn_seg = os.path.split(filename)[1].replace('.dat','')
                        pfn = '%s_segment_%.1f-%.1f_'%(mode,wmin,wmax)+dfn_seg
                        plot_filename = os.path.join(self.pplot,folder,pfn)
                    extra_stats = dict([('line_labels',lls),\
                                        ('histoplot',not exclude_data \
                                                        and [0] or []),\
                                        ('filename',plot_filename)])
                    w = [wave]*(len(sphinx_flux)+(not exclude_data and 1 or 0))
                    f = exclude_data and sphinx_flux or [flux]+sphinx_flux
                    plot_filename = Plotting2.plotCols(x=w,y=f,cfg=cfg_dict,\
                                                       **extra_stats)
                    print '** Segment finished and saved at:'
                    print plot_filename
                    
                    

    def plotSpire(self,star_grid=[],models=[],exclude_data=0,fn_plt='',\
                  fn_trans_marker='',number_subplots=3,cfg='',fn_add_star=1):
        
        '''
        Plot SPIRE data along with Sphinx results. In flux (Jy) vs wavelength
        (micron)
        
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
        @keyword fn_trans_marker: A file that includes TRANSITION definitions.
                                  These transitions will be marked up in the 
                                  plot. For instance, when indicating a subset 
                                  of transitions for one reason or another.
                                  The line type can be set for this specific 
                                  subset, differently from other lines and 
                                  regardless of the molecule. In combination
                                  with a doubly defined line label (through the
                                  star_grid['GAS_LINES']), lines can be marked
                                  up.
                                  
                                  (default: '')
        @type fn_trans_marker: string
        @keyword fn_plt: A plot filename to which an index is added for each
                         subband.
                         
                         (default: '')
        @type fn_plt: string
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        @keyword fn_add_star: Add the star name to the requested plot filename.
                              Only relevant if fn_plt is given.
                              
                              (default: 1)
        @type fn_add_star: bool
        
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

        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        if cfg_dict.has_key('exclude_data'):
            exclude_data = bool(cfg_dict['exclude_data'])
        if cfg_dict.has_key('fn_trans_marker'):
            fn_trans_marker = cfg_dict['fn_trans_marker']
        if cfg_dict.has_key('fn_add_star'):
            fn_add_star = cfg_dict['fn_add_star']
        self.spire.prepareSphinx(star_grid)
        lls = self.createLineLabels(star_grid,instrument='SPIRE')

        if fn_trans_marker:
            used_indices = list(set([ll[-2] for ll in lls]))
            this_index = [ii for ii in range(100) if ii not in used_indices][0]
            ells = self.createLineLabels(fn_trans_marker=fn_trans_marker,\
                                         ilabel=this_index,instrument='SPIRE')
            lls = lls + ells
        plot_filenames = []
        for wave,flux,band,filename in zip(self.spire.data_wave_list,\
                                           self.spire.data_flux_list,\
                                           self.spire.data_ordernames,\
                                           self.spire.data_filenames):
            if fn_plt and not fn_add_star:
                fn_plt = os.path.splitext(fn_plt)[0]
                this_filename = '%s_%s'%(fn_plt,band)
            elif fn_plt and fn_add_star:
                fn_plt = os.path.splitext(fn_plt)[0]
                this_filename = '%s_%s_%s'%(fn_plt,self.star_name,band)
            else:
                pfn = 'spire_spectrum_%s'%band
                this_filename = os.path.join(self.pplot,pfn)
            sphinx_flux = [self.spire.getSphinxConvolution(star,filename)[1]
                           for star in star_grid]
            w = exclude_data \
                          and [wave]*(len(sphinx_flux)) \
                          or [wave]*(len(sphinx_flux)+1)
            f = exclude_data and sphinx_flux or [flux] + sphinx_flux 
            keytags = ['Model %i: %s'%(i+1,star['LAST_GASTRONOOM_MODEL']\
                                    .replace('_','\_')) 
                       for i,star in enumerate(star_grid)]
            if not exclude_data: 
                keytags = ['Spire Spectrum'] + keytags
            plot_filenames.append(Plotting2.plotCols(x=w,y=f,\
                keytags=keytags,number_subplots=3,\
                line_label_color=1,line_labels=lls,\
                plot_title='%s: %s' %(self.plot_id.replace('_','\_'),\
                                      self.star_name_plots),\
                histoplot= not exclude_data and [0] or [],\
                filename=this_filename,cfg=cfg_dict,\
                line_label_spectrum=1))
        if plot_filenames and plot_filenames[0][-4:] == '.pdf':
            if not fn_plt:
                newf = os.path.join(self.pplot,'SPIRE_spectrum.pdf')
            elif fn_plt and fn_add_star:
                newf = '%s_%s.pdf'%(fn_plt,self.star_name)
            else:
                newf = '%s.pdf'%fn_plt
            DataIO.joinPdf(old=sorted(plot_filenames,reverse=True),new=newf,\
                           del_old=not fn_plt)
            print '** Your plots can be found at:'
            print newf
            print '***********************************'
        else:
            print '** Your plots can be found at:'
            print '\n'.join(plot_filenames)
            print '***********************************'

