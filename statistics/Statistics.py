# -*- coding: utf-8 -*-

"""
Interface for statistics.

Author: R. Lombaert

"""

import os
import scipy
import scipy.stats

from cc.tools.io import DataIO
from cc.data.instruments import Pacs



class Statistics(object):
    
    """
    An evironment for various statistical calculations for transitions.
    
    """
        
    def __init__(self,star_name,path_code='codeSep2010',\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),\
                 code='GASTRoNOoM'):
        
        """ 
        Initializing an instance of Statistics.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2010')
        @type path_code: string
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        @keyword code: the code used for producing your output 
        @type code: string
        
        """

        self.star_name = star_name
        self.path_combocode = path_combocode
        self.instruments = dict()
        self.star_grid = []
        self.path_code = path_code
        self.code = code
        self.data_stats = dict()
        


    def setInstrument(self,instrument,data_path='',searchstring='',\
                      data_filenames=[],instrument_instance=None,\
                      pacs_oversampling=None):
       
        '''
        Set an instrument as a source of data.
        
        @param instrument: The instrument (such as 'PACS')
        @type instrument: string
        
        @keyword data_path: full path to the data folder, excluding star_name. 
                            Not required if instrument_instance is given.
                            
                            (default: '')
        @type data_path: string
        @keyword data_filenames: The data filenames. If empty, auto-search is 
                                 done, with searchstring.
                                 
                                 (default: [])
        @type data_filenames: list[string]
        @keyword searchstring: the searchstring conditional for the auto-search
        
                               (default: '')
        @type searchstring: string
        @keyword instrument_instance: If None a new instance is made, otherwise
                                      the instance is given here for the 
                                      instrument
                                      
                                      (default: None)
        @type instrument_instance: Instrument()
        @keyword pacs_oversampling: The PACS instrumental oversampling, for 
                                    correct convolution of the Sphinx output.
                                    Not required if instrument_instance is 
                                    given
                                    
                                    (default: None)
        @type pacs_oversampling: int
        
        '''
        
        if instrument.upper() == 'PACS':
            if instrument_instance <> None:
                self.instruments['PACS'] = instrument_instance
                self.instruments['PACS'].setData(data_filenames=data_filenames,\
                                                 searchstring=searchstring)
            else:
                self.instruments['PACS'] = Pacs.Pacs(star_name=self.star_name,\
                                                     pacs_oversampling\
                                                        =pacs_oversampling,\
                                                     path=self.path_code,\
                                                     path_pacs=data_path,\
                                                     path_combocode\
                                                        =self.path_combocode)
                self.instruments['PACS'].setData(data_filenames=data_filenames,\
                                                 searchstring=searchstring)
            self.doDataStats('PACS')


            
    def setModels(self,instrument,star_grid=[],models=[],**extra_keywords):
        
        '''
        Load the models and remember them.
        
        @param instrument: The instrument (such as 'PACS')
        @type instrument: string
        
        @keyword star_grid: The parameter sets, if not given: model ids needed
        
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: the PACS ids, only relevant if star_grid == []
        
                         (default: [])
        @type models: list[string]
        @keyword extra_keywords: any instrument specific keywords that you need
        
                                 (default: dict())
        @type extra_keywords: dict
        
        '''
        if instrument.upper() == 'PACS':
            print '** Checking Sphinx models for comparison with PACS.'
            if not star_grid and not models:
                raise IOError('Statistics.setModels requires either ' + \
                              'star_grid or models to be defined.')
            elif not star_grid:
                star_grid = self.instruments['PACS'].makeStars(models=models)
            if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) == set([0]): 
                return
            self.instruments['PACS'].prepareSphinx(star_grid)
        else:
            raise IOError('Instruments other than PACS not yet implemented.')
        self.star_grid = [star 
                          for star in star_grid if star['LAST_PACS_MODEL']]
        print '***********************************'


    
    def doDataStats(self,instrument):
        
        '''
        Calculate means, tmeans, stds and tstds for the data files.
        
        They are saved in self.data_stats[instrument].
        
        @param instrument: The instrument (such as 'PACS')
        @type instrument: string
        
        '''
        
        stats = dict()
        for filename,data_flux in zip(self.instruments[instrument]\
                                          .data_filenames,\
                                      self.instruments[instrument]\
                                          .data_flux_list):
            filename = os.path.split(filename)[1]
            these_stats = dict()
            these_stats['mean'] \
                = scipy.mean(data_flux[scipy.isfinite(data_flux)])             
            these_stats['std'] \
                = scipy.std(data_flux[scipy.isfinite(data_flux)])
            these_stats['tmean'] \
                = scipy.stats.tmean(data_flux[scipy.isfinite(data_flux)],\
                                    limits=(these_stats['mean']\
                                                -these_stats['std'],\
                                            these_stats['mean']\
                                                +these_stats['std']))
            these_stats['tstd'] \
                = scipy.stats.tstd(data_flux[scipy.isfinite(data_flux)],\
                                    limits=(these_stats['mean']\
                                                -these_stats['std'],\
                                            these_stats['mean']\
                                                +these_stats['std']))
            print '* Data statistics for %s:'%filename
            print 'mean = ',these_stats['mean'],' Jy'
            print 'std = ',these_stats['std'],' Jy'      
            print 'tmean = ',these_stats['tmean'],' Jy'      
            print 'tstd = ',these_stats['tstd'],' Jy'    
            stats[filename] = these_stats
        self.data_stats[instrument] = stats
        print '***********************************'
        