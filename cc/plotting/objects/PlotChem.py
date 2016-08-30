# -*- coding: utf-8 -*-

"""
A plotting environment for Chemistry output.

Author: M. Van de Sande

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


class PlotChem(PlottingSession):
    
    """ 
    Class for plotting gas lines and their information.
    
    """    
    
    def __init__(self,star_name,path_chemistry='OutputClumpy',\
                 inputfilename=None,fn_add_star=1):
        
        """ 
        Initializing an instance of PlotGas.
        
        @param star_name: name of the star from Star.dat, use default only 
                          when never using any star model specific things 
        @type star_name: string
        
        @keyword path_chemistry: Output modeling folder in MCMax home folder
        
                                  (default: 'OutputClumpy')
        @type path_chemistry: string
        @keyword inputfilename: name of inputfile that is also copied to the 
                                output folder of the plots, 
                                if None nothing is copied
                                
                                (default: None)
        @type inputfilename: string

        @keyword fn_add_star: Add the star name to the requested plot filename.
                              Only relevant if fn_plt is given in a sub method.
                              
                              (default: 1)
        @type fn_add_star: bool
                            
        
        """
        
        super(PlotChem, self).__init__(star_name=star_name,\
                                      path=path_chemistry,\
                                      code='Chemistry',\
                                      inputfilename=inputfilename,\
                                      fn_add_star=fn_add_star)
        #-- Convenience path
        cc.path.cout = os.path.join(cc.path.chemistry,self.path)


        
    def makeStars(self,models):
        
        '''
        Make a Star list based on the Chemistry model ids.
        
        @param models: model_ids for the Chemistry db
        @type models: list(string)
        
        @return: the parameter sets
        @rtype: list[Star()]
        
        '''
        
        star_grid = Star.makeStars(models=models,\
                                   id_type='Chemistry',\
                                   code='Chemistry',path=self.path)
        [star.addCoolingPars() for star in star_grid]
        return star_grid



    def plotVelocity(self,star_grid=[],models=[],fn_plt='',force_plot=0,cfg=''):
        
        '''
        Plot velocity versus radius for every model in star_grid.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]    
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string        
        @keyword force_plot: force a plotting if more than models are requested
                             
                             (default: 0)
        @type force_plot: bool        
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Plotting Gas Velocity Profiles'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
        
        if len(star_grid) < 20 or force_plot:
            
            valid_sg = [star for star in star_grid 
                        if star['LAST_CHEMISTRY_MODEL']]
            
            folders = [os.path.join(cc.path.cout,'models',\
                    star['LAST_CHEMISTRY_MODEL'])+'/' for star in valid_sg]
            
            radii = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'RADIUS') \
                for folder in folders]
            temps = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'VELOCITY') \
                for folder in folders]

            if temps:    
                keytags = [star['LAST_CHEMISTRY_MODEL'].replace('_','\_')
                                for star in valid_sg]
                
                #-- Set filenames
                pfn = fn_plt if fn_plt else 'velocity_profiles'
                pfn = self.setFnPlt(pfn)
                
                #-- Run the two plots
                keys_cm = ['Model %i'%(i+1)
                           for i in xrange(len(star_grid))]
                pfn = Plotting2.plotCols(x=radii,y=temps,cfg=cfg_dict,\
                        filename=pfn,xaxis='$r$ (cm)',\
                        yaxis=r'$v$ (km s$^{-1}$)',\
                        figsize=(12.5,8),fontsize_ticklabels=26,\
                        key_location=(0.05,0.05),xlogscale=1,ylogscale=1,\
                        keytags=keys_cm,fontsize_axis=26,fontsize_key=26)
                print '** Plots can be found at:'
                print pfn
                print '***********************************'
        

    def plotTemp(self,star_grid=[],models=[],fn_plt='',force_plot=0,cfg=''):
        
        '''
        Plot temperature profiles of all models.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string    
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
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
            
        if len(star_grid) < 20 or force_plot:
            
            valid_sg = [star for star in star_grid 
                        if star['LAST_CHEMISTRY_MODEL']]
            
            folders = [os.path.join(cc.path.cout,'models',\
                    star['LAST_CHEMISTRY_MODEL'])+'/' for star in valid_sg]
            
            radii = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'RADIUS') \
                for folder in folders]
            temps = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'TEMP') \
                for folder in folders]

            if temps:    
                keytags = [star['LAST_CHEMISTRY_MODEL'].replace('_','\_')
                                for star in valid_sg]
                
                #-- Set filenames
                pfn = fn_plt if fn_plt else 'temperature_profiles'
                pfn = self.setFnPlt(pfn)
                
                #-- Run the two plots
                keys_cm = ['Model %i'%(i+1)
                           for i in xrange(len(star_grid))]
                pfn = Plotting2.plotCols(x=radii,y=temps,cfg=cfg_dict,\
                        filename=pfn,xaxis='$r$ (cm)',\
                        yaxis='$T_\mathrm{g}$ (K)',\
                        figsize=(12.5,8),fontsize_ticklabels=26,\
                        key_location=(0.05,0.05),xlogscale=1,ylogscale=1,\
                        keytags=keys_cm,fontsize_axis=26,fontsize_key=26)
                print '** Plots can be found at:'
                print pfn
                print '***********************************'
            else:
                print '** No Chemistry models were calculated successfully.'+\
                      'No temperature profiles can be plotted.'
                print '***********************************'



    def plotVisualExtinction(self,star_grid=[],models=[],fn_plt='',force_plot=0,cfg=''):
        
        '''
        Plot temperature profiles of all models.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string    
        @keyword force_plot: force a plotting if more than models are requested
                             
                             (default: 0)
        @type force_plot: bool
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Plotting Visual Extinction  Profiles'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
            
        if len(star_grid) < 20 or force_plot:
            
            valid_sg = [star for star in star_grid 
                        if star['LAST_CHEMISTRY_MODEL']]
            
            folders = [os.path.join(cc.path.cout,'models',\
                    star['LAST_CHEMISTRY_MODEL'])+'/' for star in valid_sg]
            
            radii = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'RADIUS') \
                for folder in folders]
            avs = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'A_V') \
                for folder in folders]

            if avs:    
                keytags = [star['LAST_CHEMISTRY_MODEL'].replace('_','\_')
                                for star in valid_sg]
                
                #-- Set filenames
                pfn = fn_plt if fn_plt else 'visual_extinction'
                pfn = self.setFnPlt(pfn)
                
                #-- Run the two plots
                keys_cm = ['Model %i'%(i+1)
                           for i in xrange(len(star_grid))]
                pfn = Plotting2.plotCols(x=radii,y=avs,cfg=cfg_dict,\
                        filename=pfn,xaxis='$r$ (cm)',\
                        yaxis='$A_V$ (mag)',\
                        figsize=(12.5,8),fontsize_ticklabels=26,\
                        key_location=(0.05,0.05),xlogscale=1,ylogscale=1,\
                        keytags=keys_cm,fontsize_axis=26,fontsize_key=26)
                print '** Plots can be found at:'
                print pfn
                print '***********************************'
            else:
                print '** No Chemistry models were calculated successfully.'+\
                      'No temperature profiles can be plotted.'
                print '***********************************'


    def plotPhotodissociationrateCo(self,star_grid=[],models=[],fn_plt='',force_plot=0,cfg=''):
        
        '''
        Plot temperature profiles of all models.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string    
        @keyword force_plot: force a plotting if more than models are requested
                             
                             (default: 0)
        @type force_plot: bool
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Plotting Photodissociation Rate CO'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
            
        if len(star_grid) < 20 or force_plot:
            
            valid_sg = [star for star in star_grid 
                        if star['LAST_CHEMISTRY_MODEL']]
            
            folders = [os.path.join(cc.path.cout,'models',\
                    star['LAST_CHEMISTRY_MODEL'])+'/' for star in valid_sg]
            
            radii = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'RADIUS') \
                for folder in folders]
            cos = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'COK(PHOT)') \
                for folder in folders]

            if cos:    
                keytags = [star['LAST_CHEMISTRY_MODEL'].replace('_','\_')
                                for star in valid_sg]
                
                #-- Set filenames
                pfn = fn_plt if fn_plt else 'CO_photodissociation_rate'
                pfn = self.setFnPlt(pfn)
                
                #-- Run the two plots
                keys_cm = ['Model %i'%(i+1)
                           for i in xrange(len(star_grid))]
                pfn = Plotting2.plotCols(x=radii,y=cos,cfg=cfg_dict,\
                        filename=pfn,xaxis='$r$ (cm)',\
                        yaxis='$k_{CO,photo}$ (s$^{-1}$)',\
                        figsize=(12.5,8),fontsize_ticklabels=26,\
                        key_location=(0.05,0.05),xlogscale=1,ylogscale=1,\
                        keytags=keys_cm,fontsize_axis=26,fontsize_key=26)
                print '** Plots can be found at:'
                print pfn
                print '***********************************'
            else:
                print '** No Chemistry models were calculated successfully.'+\
                      'No temperature profiles can be plotted.'
                print '***********************************'



    def plotRadiationField(self,star_grid=[],models=[],fn_plt='',force_plot=0,cfg=''):
        
        '''
        Plot temperature profiles of all models.
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string    
        @keyword force_plot: force a plotting if more than models are requested
                             
                             (default: 0)
        @type force_plot: bool
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        '''
        
        print '***********************************'
        print '** Plotting Radiation Field'
        if not star_grid and models:
            star_grid = self.makeStars(models=models)
        elif (not models and not star_grid) or (models and star_grid):
            print '** Input is undefined or doubly defined. Aborting.'
            return
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
            
        if len(star_grid) < 20 or force_plot:
            
            valid_sg = [star for star in star_grid 
                        if star['LAST_CHEMISTRY_MODEL']]
            
            folders = [os.path.join(cc.path.cout,'models',\
                    star['LAST_CHEMISTRY_MODEL'])+'/' for star in valid_sg]
            
            radii = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'RADIUS') \
                for folder in folders]
            cos = [DataIO.getChemistryPhysPar(folder+'csphyspar.out', 'COK(PHOT)') \
                for folder in folders]

            if cos:    
                keytags = [star['LAST_CHEMISTRY_MODEL'].replace('_','\_')
                                for star in valid_sg]
                
                #-- Set filenames
                pfn = fn_plt if fn_plt else 'radiation_field'
                pfn = self.setFnPlt(pfn)
                
                #-- Run the two plots
                keys_cm = ['Model %i'%(i+1)
                           for i in xrange(len(star_grid))]
                pfn = Plotting2.plotCols(x=radii,y=cos,cfg=cfg_dict,\
                        filename=pfn,xaxis='$r$ (cm)',\
                        yaxis='$I_D$ (photons /s /cm$^2$ /Hz / sr) ',\
                        figsize=(12.5,8),fontsize_ticklabels=26,\
                        key_location=(0.05,0.05),xlogscale=1,ylogscale=1,\
                        keytags=keys_cm,fontsize_axis=26,fontsize_key=26)
                print '** Plots can be found at:'
                print pfn
                print '***********************************'
            else:
                print '** No Chemistry models were calculated successfully.'+\
                      'No temperature profiles can be plotted.'
                print '***********************************'


    
    def plotAbundanceProfiles(self,star_grid=[],models=[],force_plot=0,cfg='',\
                              fn_plt='',molecules=[],per_molecule=0,frac=1):  
        
        '''
        Plot abundance profiles for all molecules in every model.
        
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
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string
        @keyword molecules: Molecules to be plotted.
        
                               (default: [])
        @type molecules: bool
        @keyword per_molecule: Plot one molecule for all models in one figure.
        
                               (default: 0)
        @type per_molecule: bool
        @keyword per_model: Plot all molecules for one model in one figure.
        
                               (default: 0)
        @type per_model: bool
        @keyword frac: Plot the fractional abundances. If not frac, plot number
                       densities.
                       
                                (default: 1)

        
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
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
        if cfg_dict.has_key('molecules'):
            molecules = cfg_dict.pop('molecules')
        if cfg_dict.has_key('per_molecule'):
            per_molecule = cfg_dict['per_molecule']
        if cfg_dict.has_key('per_model'):
            per_model = cfg_dict['per_model']            
        
        #-- Some general plot settings
        extra_pars = dict()
        extra_pars['ymin'] = 1e-9
        extra_pars['ymax'] = 1e-3
        extra_pars['ylogscale'] = 1 
        extra_pars['xlogscale'] = 1
        extra_pars['figsize'] = (12.5,8.5)
        extra_pars['xaxis'] = 'cm'
        
        #-- Dict to keep track of all data
        ddata = dict()
        for istar,star in enumerate(star_grid):
            if not star['LAST_CHEMISTRY_MODEL']: continue
            ddata[istar] = dict()
            
            folder = os.path.join(cc.path.cout,'models',\
                        star['LAST_CHEMISTRY_MODEL'])+'/'
            ddata[istar]['rad'] = DataIO.getChemistryPhysPar(folder+\
                'csphyspar.out', 'RADIUS')
            ddata[istar]['id'] = star['LAST_CHEMISTRY_MODEL']
            if frac:
                species = DataIO.getChemistryAbundances(folder+'csfrac.out')
            else:
                species = DataIO.getChemistryAbundances(folder+'csnum.out')
            
            for molec in molecules: 
                ddata[istar][molec] = species[molec]
                
            if not per_molecule:
                #-- Collect all data
                radii = [ddata[istar]['rad']]*len(molecules)
                abuns = [ddata[istar][molec] for molec in molecules]
                keytags = molecules
                #ids = star['LAST_CHEMISTRY_MODEL']
                ids = ddata[istar]['id']
                
                #-- Set the yaxis tag
                yaxis = '$n_\mathrm{molec}/n_{\mathrm{H}_2}$'
                
                #-- Set filename
                pfn = fn_plt if fn_plt else 'abundance_profiles'
                suff = '_'.join(list(set(ids)))
                pfn = self.setFnPlt(pfn,fn_suffix=suff)

                pfns.append(Plotting2.plotCols(x=radii,y=abuns,cfg=cfg_dict,\
                                               filename=pfn,keytags=keytags,\
                                               plot_title=ids.replace('_','\_'),\
                                               yaxis=yaxis,**extra_pars))
        
        if per_molecule:
            #-- Collect all data
            #molecs = list(set([molec for istar in ddata.keys()
                                     #for molec in ddata[istar].keys()]))
            for molec in molecules: 
                #-- Collect data
                radii = [dstar['rad']
                         for istar,dstar in ddata.items()]
                abuns = [dstar[molec]
                         for istar,dstar in ddata.items()]
                keytags = [dstar['id'].replace('_','\_') 
                           for istar,dstar in ddata.items()]

                #-- Set the y axis tag
                #strmolec = ddata[0][molec]['key']
                yaxis = '$n_\mathrm{%s}/n_{\mathrm{H}_2}$'%str(molec)

                #-- Make filename
                pfn = fn_plt if fn_plt else 'abundance_profiles'
                pfn = self.setFnPlt(pfn,fn_suffix=molec)

                pfns.append(Plotting2.plotCols(x=radii,y=abuns,yaxis=yaxis,\
                                               filename=pfn,keytags=keytags,\
                                               cfg=cfg_dict,**extra_pars))  
        
        if not per_molecule and pfns and pfns[0][-4:] == '.pdf':    
            pfn = fn_plt if fn_plt else 'abundance_profiles'
            pfn = self.setFnPlt(pfn) + '.pdf'
            DataIO.joinPdf(old=pfns,new=pfn)
            print '** Plots can be found at:'
            print pfn
            print '***********************************'
        else:
            print '** Plots can be found at:'
            print '\n'.join(pfns)
            print '***********************************'
            
            



        
