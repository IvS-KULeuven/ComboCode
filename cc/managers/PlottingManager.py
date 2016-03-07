# -*- coding: utf-8 -*-

"""
Interface for plotting.

Author: R. Lombaert

"""

import os

import cc.path
from cc.plotting.objects import PlotGas
from cc.plotting.objects import PlotDust


class PlottingManager():
    
    """ 
    An interface for managing requested plots for a CC session.
    
    """
    
    def __init__(self,star_name,mcmax=False,gastronoom=False,pacs=None,\
                 path_gastronoom='codeJun2010',path_mcmax='codeJun2010',\
                 inputfilename='inputComboCode.dat',spire=None,fn_add_star=1,\
                 plot_pars=dict(),sed=None):
                
        """ 
        Initializing a PlottingManager instance.
        
        @keyword star_name: name of the star from Star.dat, use default only 
                            when never using any star model specific things 
                                  
                            (default: "model")
        @type star_name: string
        @keyword inputfilename: name of inputfile that is also copied to the 
                                output folder of the plots, 
                                if None nothing is copied
                                
                                (default: None)
        @type inputfilename: string
        @keyword mcmax: Running MCMax?
        
                        (default: 0)
        @type mcmax: bool
        @keyword gastronoom: Running GASTRoNOoM?
        
                             (default: 0)
        @type gastronoom: bool
        @keyword path_mcmax: modeling folder in MCMax home
        
                             (default: 'runTest')
        @type path_mcmax: string
        @keyword path_gastronoom: modeling folder in GASTRoNOoM home
        
                                  (default: 'runTest')
        @type path_gastronoom: string
        @keyword pacs: A Pacs() object for managing data and model handling for
                       PACS spectra. None if not applicable
                            
                       (default: None)
        @type pacs: Pacs()
        @keyword spire: A Spire() object for managing data and model handling
                        for SPIRE spectra. None if not applicable
                            
                        (default: None)
        @type spire: Spire()
        @keyword sed: The SED for managing data and model handling
                      for SED spectra/photometry. None if not applicable
                          
                      (default: None)
        @type sed: Sed()
        @keyword fn_add_star: Add the star name to the requested plot filename.
                              Only relevant if fn_plt is given in a sub method.
                              
                              (default: 1)
        @type fn_add_star: bool
        @keyword plot_pars: dictionary with all the plotting parameters that
                            turn on or off plotting modules. By default they
                            are all turned off.
                                  
                            (default: dict())
        @type plot_pars: dict
        
        """
        
        self.dust_pars = dict()
        self.dust_cfg = dict()
        self.gas_pars = dict()
        self.gas_cfg = dict()
        for k,v in plot_pars.items():
            if k[0:9] == 'PLOT_GAS_' and v:
                self.gas_pars[k.replace('_GAS','',1)] = v
            elif k[0:8] == 'CFG_GAS_' and v:
                self.gas_cfg[k.replace('_GAS','',1)] = v
            elif k[0:10] == 'PLOT_DUST_' and v:
                self.dust_pars[k.replace('_DUST','',1)] = v
            elif k[0:9] == 'CFG_DUST_' and v:
                self.dust_cfg[k.replace('_DUST','',1)] = v
        self.mcmax = mcmax
        self.gastronoom = gastronoom
        if self.mcmax: 
            self.plotter_dust = PlotDust.PlotDust(star_name=star_name,\
                                                  path_mcmax=path_mcmax,\
                                                  inputfilename=inputfilename,
                                                  sed=sed,\
                                                  fn_add_star=fn_add_star)
        else: 
            self.plotter_dust = None
        if self.gastronoom or self.gas_pars.has_key('PLOT_LINE_LISTS')\
                or 'PLOT_TRANSITIONS' in self.gas_pars:
            self.plotter_gas = PlotGas.PlotGas(star_name=star_name,pacs=pacs,\
                                               path_gastronoom=path_gastronoom,\
                                               inputfilename=inputfilename,\
                                               spire=spire,\
                                               fn_add_star=fn_add_star)                                     
        else:
            self.plotter_gas = None
        

    
    def startPlotting(self,star_grid,iterative=0):
        
        """ 
        Start plotting PLOT_INPUT requests for those models that are available 
        (i.e. gastronoom and/or mcmax).
        
        @param star_grid: list of stars to be plotted
        @type star_grid: list[Star()]
        
        @keyword iterative: if true the old grids are plotted on a 
                            model_iteration per model_iteration basis, only
                            works for MCMax models for now. 
                                  
                            (default: 0)
        @type iterative: int
        
        """
        
        if iterative:
            if self.mcmax and self.dust_pars.has_key('PLOT_SED'):
                self.plotter_dust.plotSed(star_grid=star_grid,\
                                          iterative=iterative,\
                                          cfg=self.dust_cfg.has_key('CFG_SED')\
                                                and self.dust_cfg['CFG_SED'] \
                                                or '')
            return
        if self.mcmax:
            for k in self.dust_pars:
                method_name = 'plot' + \
                              ''.join([w.capitalize() 
                                       for w in k.replace('PLOT_','')\
                                                 .split('_')])
                thisMethod = getattr(self.plotter_dust,method_name)
                thisMethod(star_grid=star_grid,\
                           cfg=self.dust_cfg.get(k.replace('PLOT_','CFG_'),''))
        if self.gastronoom or self.gas_pars.has_key('PLOT_LINE_LISTS') \
                 or 'PLOT_TRANSITIONS' in self.gas_pars:
            for k in self.gas_pars:
                method_name = 'plot' + \
                              ''.join([w.capitalize() 
                                       for w in k.replace('PLOT_','')\
                                                 .split('_')])
                thisMethod = getattr(self.plotter_gas,method_name)
                thisMethod(star_grid=star_grid,\
                           cfg=self.gas_cfg.get(k.replace('PLOT_','CFG_'),''))
    
    
    
    def plotTransitions(self,star_grid,force=0,cfg='',fn_suffix=''):
        
        '''
        Run the plotTransitions method in the PlotGas object of this manager. 
        
        If available, the cfg is taken from the plot input.
        
        Only done if PLOT_GAS_TRANSITIONS is indeed 1 in the CC inputfile. This
        can be forced through the keyword force.
        
        This method is mainly used if you want to plot a sub selection of 
        models from the main input grid, for instance after statistical 
        selection. 
        
        @param star_grid: The collection of models to be plotted
        @type star_grid: list[Star()]
        
        @keyword force: Plot the transitions regardless of the 
                        PLOT_GAS_TRANSITIONS keyword in cc input, e.g. when the
                        manager is ran outside a cc session.
                        
                        (default: 0)
        @type force: bool
        @keyword cfg: The config file. If not given, the manager checks if a
                      cfg file is available in the CC input. 
                      
                      (default: None)
        @type cfg: string
        @keyword fn_suffix: A suffix that is appended to the filename. For 
                            instance, when running the plot command for a 
                            best fit subgrid of Star() models as to not 
                            overwrite the plot of the full grid.
                            
                            (default: '')
        @type fn_suffix: string
        
        
        '''
        
        if 'PLOT_TRANSITIONS' in self.gas_pars or force: 
            if not cfg: cfg = self.gas_cfg.get('CFG_TRANSITIONS','')
            self.plotter_gas.plotTransitions(star_grid,cfg=cfg,\
                                             fn_suffix=fn_suffix)
