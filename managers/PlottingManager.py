# -*- coding: utf-8 -*-

"""
Interface for plotting.

Author: R. Lombaert

"""

import os

from cc.plotting.objects import PlotGas
from cc.plotting.objects import PlotDust
from cc.tools.io import DataIO



class PlottingManager():
    
    """ 
    An interface for managing requested plots for a CC session.
    
    """
    
    def __init__(self,star_name,mcmax=False,gastronoom=False,pacs=None,\
                 path_gastronoom='codeJun2010',path_mcmax='codeJun2010',\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),\
                 inputfilename='inputComboCode.dat',spire=None,\
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
        @keyword path_combocode: CC home folder
        
                                 (default: '/home/robinl/ComboCode')
        @type path_combocode: string
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
        @keyword plot_pars: dictionary with all the plotting parameters that
                            turn on or off plotting modules. By default they
                            are all turned off.
                                  
                            (default: dict())
        @type plot_pars: dict
        
        """
        
        self.path_combocode = path_combocode
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
        self.sed=sed
        if self.mcmax: 
            self.plotter_dust = PlotDust.PlotDust(star_name=star_name,\
                                                path_combocode=path_combocode,\
                                                path_mcmax=path_mcmax,\
                                                inputfilename=inputfilename)
        else: 
            self.plotter_dust = None
        if self.gastronoom or self.gas_pars.has_key('PLOT_LINE_LISTS'):
            self.plotter_gas = PlotGas.PlotGas(star_name=star_name,pacs=pacs,\
                                               path_combocode=path_combocode,\
                                               path_gastronoom=path_gastronoom,\
                                               inputfilename=inputfilename,\
                                               spire=spire)                                     
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
                if k == 'PLOT_SED':
                     self.plotter_dust.plotSed(star_grid=star_grid,\
                                        spec=self.sed,\
                                        cfg=self.dust_cfg.has_key('CFG_SED') \
                                                and self.dust_cfg['CFG_SED'] \
                                                or '')
                else:
                     method_name = 'plot' + \
                                   ''.join([w.capitalize() 
                                            for w in k.replace('PLOT_','')\
                                                      .split('_')])
                     thisMethod = getattr(self.plotter_dust,method_name)
                     thisMethod(star_grid=star_grid,\
                                cfg=self.dust_cfg\
                                        .get(k.replace('PLOT_','CFG_'),''))
        if self.gastronoom:
            if self.gas_pars.pop('PLOT_PACS_SEGMENTS',0):
                path_segments = self.gas_pars.pop('PLOT_PACS_SEGMENTS_PATH','')
            else:
                path_segments = ''
            for k in self.gas_pars:
                method_name = 'plot' + \
                              ''.join([w.capitalize() 
                                       for w in k.replace('PLOT_','')\
                                                 .split('_')])
                thisMethod = getattr(self.plotter_gas,method_name)
                thisMethod(star_grid=star_grid,\
                           cfg=self.gas_cfg.get(k.replace('PLOT_','CFG_'),''))
                if method_name == 'plotPacs' and path_segments:
                    self.plotter_gas.plotPacsSegments(star_grid=star_grid,\
                            pacs_segments_path=path_segments,mode='sphinx',\
                            cfg=self.gas_cfg.get('CFG_PACS_SEGMENTS',''))
                elif method_name == 'plotLineLists' and path_segments:    
                    self.plotter_gas.plotPacsSegments(star_grid=star_grid,\
                            pacs_segments_path=path_segments,mode='ll',\
                            cfg=self.gas_cfg.get('CFG_PACS_SEGMENTS',''))
                
