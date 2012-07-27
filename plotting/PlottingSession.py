# -*- coding: utf-8 -*-

"""
Interface for creating plotting environments.

Author: R. Lombaert

"""

import os
from time import gmtime
import cPickle
import subprocess

from cc.tools.io import DataIO



class PlottingSession(object):
    
    """
    Plotting environment for data and models of circumstellar shells.
    
    """
        
    def __init__(self,star_name='model',inputfilename=None,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),\
                 path='codeJun2010',code='GASTRoNOoM'):
        
        """ 
        Initializing an instance of PlottingSession.
        
        @keyword star_name: name of the star from Star.dat, use default only 
                            when never using any star model specific things 
                                  
                            (default: "model")
        @type star_name: string
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        @keyword path: Output modeling folder in code home folder
        
                       (default: 'codeJun2010')
        @type path: string
        @keyword inputfilename: name of inputfile that is also copied to the 
                                output folder of the plots, 
                                if None nothing is copied
                                
                                (default: None)
        @type inputfilename: string
        @keyword code: the code itself with which we're working
        
                       (default: GASTRoNOoM)
        @type code: string
        
        """
        
        self.inputfilename = inputfilename
        self.star_name = star_name
        self.path_combocode = path_combocode
        self.star_index = DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'))\
                              .index(self.star_name)
        self.star_name_plots = DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='STAR_NAME_PLOTS',remove_underscore=1)\
                               [self.star_index]
        self.vlsr = DataIO.getInputData(\
                        path=os.path.join(self.path_combocode,'Data'),\
                        keyword='V_LSR',remove_underscore=1)\
                        [self.star_index]
        self.path = path
        self.mutable_mcmax = [line[0] 
                              for line in DataIO.readFile(\
                                   os.path.join(self.path_combocode,'CC',\
                                                'Mutable_Parameters_MCMax.dat'),\
                                   delimiter=' ')
                              if ''.join(line).strip()]
        self.mutable_mcmax = [line 
                              for line in self.mutable_mcmax 
                              if line[0] != '#']
        self.mutable_gastronoom = [line[0] 
                                   for line in DataIO.readFile(\
                                        os.path.join(self.path_combocode,'CC',\
                                        'Mutable_Parameters_GASTRoNOoM.dat'),\
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
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                self.code,self.path,'stars'))
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                self.code,self.path,'stars',\
                                    self.star_name))
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                self.code,self.path,'stars',\
                                                self.star_name,self.plot_id))
        if self.inputfilename <> None:
            newf = os.path.join(os.path.expanduser('~'),self.code,self.path,\
                                'stars',self.star_name,self.plot_id,\
                                os.path.split(self.inputfilename)[1])
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
        
        
          
    def makeModelList(self,star_grid,code):
        
        '''
        Return a list of all model id's in the star_grid.
                
        @param star_grid: The parameter sets
        @type star_grid: list[Star()]
        @param code: the code itself with which we're working
        
                     (default: GASTRoNOoM)
        @type code: string
        
        @return: the model_ids
        @rtype: list[string]
        
        '''
        
        return [star['LAST_'+code.upper()+'_MODEL'] 
                for star in star_grid 
                if star['LAST_'+code.upper()+'_MODEL']]
        