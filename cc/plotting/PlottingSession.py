# -*- coding: utf-8 -*-

"""
Interface for creating plotting environments.

Author: R. Lombaert

"""

import os
from time import gmtime
import cPickle
import subprocess

import cc.path
from cc.tools.io import DataIO



class PlottingSession(object):
    
    """
    Plotting environment for data and models of circumstellar shells.
    
    """
        
    def __init__(self,star_name='model',inputfilename=None,\
                 path='',code='GASTRoNOoM'):
        
        """ 
        Initializing an instance of PlottingSession.
        
        @keyword star_name: name of the star from Star.dat, use default only 
                            when never using any star model specific things 
                                  
                            (default: "model")
        @type star_name: string
        @keyword path: Output modeling folder in code home folder
        
                       (default: '')
        @type path: string
        @keyword inputfilename: name of inputfile that is also copied to the 
                                output folder of the plots, 
                                if None nothing is copied
                                
                                (default: None)
        @type inputfilename: string
        @keyword code: the modeling code
        
                       (default: GASTRoNOoM)
        @type code: string
        
        """
        
        self.inputfilename = inputfilename
        self.star_name = star_name
        self.star_index = DataIO.getInputData(path=cc.path.usr).index(star_name)
        self.star_name_plots = DataIO.getInputData(path=cc.path.usr,
                                                   keyword='STAR_NAME_PLOTS',\
                                                   remove_underscore=1,\
                                                   rindex=self.star_index)
        
        #-- Can't use convenience paths here through cc.path, because the 
        #   module is not code specific. Within a single python session, there
        #   may be multiple instances of PlottingSession
        if not path:
            print('Warning! %s model output folder not set.'%code)
        self.path = path
        fn_mcm = os.path.join(cc.path.aux,'Mutable_Parameters_MCMax.dat')
        self.mutable_mcmax = [line[0] 
                              for line in DataIO.readFile(fn_mcm,delimiter=' ')
                              if ''.join(line).strip()]
        self.mutable_mcmax = [line 
                              for line in self.mutable_mcmax 
                              if line[0] != '#']
        fn_gas = os.path.join(cc.path.aux,'Mutable_Parameters_GASTRoNOoM.dat')
        self.mutable_gastronoom = [line[0] 
                                   for line in DataIO.readFile(fn_gas,\
                                                               delimiter=' ')
                                   if ''.join(line).strip()]
        self.mutable_gastronoom = [line 
                                   for line in self.mutable_gastronoom 
                                   if line[0] != '#']
        self.mutable = self.mutable_mcmax + self.mutable_gastronoom
        self.plot_id = 'plot_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i' \
                       %(gmtime()[0],gmtime()[1],gmtime()[2],\
                         gmtime()[3],gmtime()[4],gmtime()[5])
        self.code = code
        
        #-- Folder management and copying inputfile to plot output folder
        pout = os.path.join(getattr(cc.path,self.code.lower()),self.path,\
                            'stars')
        pstar = os.path.join(pout,self.star_name)
        self.pplot = os.path.join(pstar,self.plot_id)
        for pp in [pout,pstar,self.pplot]:
            DataIO.testFolderExistence(pp)
        
        if self.inputfilename <> None:
            ipfn = os.path.split(self.inputfilename)[1]
            newf = os.path.join(self.pplot,ipfn)
            subprocess.call(['cp %s %s'%(self.inputfilename,newf)],shell=True)
            
            
            
    def checkChangedPars(self,star_grid):
        
        """ 
        Check which parameters change throughout the modeling in the mutable 
        list.
        
        @param star_grid: The parameter sets
        @type star_grid: list[Star()]
        
        @return: parameters in the mutable list that change
        @rtype: list[string]
        
        """
        
        all_values = [[par.find('FILE')==-1 \
                            and (star[par],) \
                            or (bool(star[par]),) 
                       for star in star_grid]
                      for par in self.mutable ]
        return [self.mutable[par_index] 
                for par_index,value_list in enumerate(all_values) 
                if len(set(value_list)) > 1]
        #- set(list) returns a 'set' which removes duplicate items from list. 
        #- Hence, there should be only one entry if all
        #- values in the list are the same
        
        
          
    def makeModelList(self,star_grid,id_type):
        
        '''
        Return a list of model id's in the star_grid.
                
        @param star_grid: The parameter sets
        @type star_grid: list[Star()]
        @param id_type: the type of model id (MCMAX or GASTRONOOM or PACS)
        @type id_type: string
        
        @return: the model_ids
        @rtype: list[string]
        
        '''
        
        return [star['LAST_'+id_type.upper()+'_MODEL'] 
                for star in star_grid 
                if star['LAST_'+id_type.upper()+'_MODEL']]
        