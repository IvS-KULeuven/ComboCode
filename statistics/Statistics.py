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
from cc.data import Data



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
        self.data_info = dict()
        try:
            self.star_index = DataIO.getInputData(path=os.path.join(\
                                                        self.path_combocode,\
                                                        'Data'))\
                                                 .index(self.star_name)
            self.v_lsr = DataIO.getInputData(path=os.path.join(\
                                                        self.path_combocode,\
                                                        'Data'),\
                                             keyword='V_LSR')[self.star_index]
        except KeyError,ValueError: 
            print 'No (correct) star name has been supplied. ' + \
                  'No data stats will be calculated for freq-resolved lines.'
            self.v_lsr = None
            


    def setInstrument(self,instrument,data_path='',searchstring='',\
                      data_filenames=[],instrument_instance=None,\
                      pacs_oversampling=None,sample_transitions=[]):
       
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
        @keyword sample_transitions: Sample transitions used as reference for 
                                     the data files. Only relevant if 
                                     instrument=='FREQ_RESO'.
                                     
                                     (default: [])
        @type sample_transitions: list
        
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
        
        elif instrument.upper() == 'FREQ_RESO' and self.v_lsr <> None:
            #- Make copy so that any changes do not translate to whatever the 
            #- original list of transitions might be. 
            self.instruments['FREQ_RESO'] = list(sample_transitions)
            [t.readData(self.v_lsr) for t in self.instruments['FREQ_RESO']]
            if not sample_transitions:
                print 'WARNING! No sample transitions given for Statistics ' +\
                      'module with instrument == FREQ_RESO. No stats will ' + \
                      'be calculated.'
            

            
    def setModels(self,instrument,star_grid=[],models=[],**extra_keywords):
        
        '''
        Load the models and remember them.
        
        @param instrument: The instrument (such as 'PACS')
        @type instrument: string
        
        @keyword star_grid: The parameter sets, if not given: model ids needed
        
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: the PACS ids, only relevant if star_grid == [] and if 
                         instrument == PACS. In all other cases a star_grid is
                         required.
        
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
        elif instrument.upper() == 'FREQ_RESO':
            print '** Checking Sphinx models for comparison with freq-resolved data.'
            if not star_grid:
                raise IOError('Statistics.setModels requires a ' + \
                              'star_grid to be defined for freq-resolved data.')
            if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) == set([0]): 
                return
            [[t.readSphinx() for t in s['GAS_LINES']] for s in star_grid]
        else:
            raise IOError('Instruments other than PACS or freq-resolved data not yet implemented.')
        self.star_grid = star_grid
        print '***********************************'


    
    def doDataStats(self,instrument):
        
        '''
        Calculate means, tmeans, stds and tstds for the data files.
        
        They are saved in self.data_stats[instrument].
        
        @param instrument: The instrument (such as 'PACS')
        @type instrument: string
        
        '''
        
        stats = dict()
        instrument = instrument.upper()
        if instrument == 'PACS':
            self.getDataInfo(instrument)
            for filename,data_wave,data_flux,band in \
                        zip(self.instruments[instrument].data_filenames,\
                            self.instruments[instrument].data_wave_list,\
                            self.instruments[instrument].data_flux_list,\
                            self.instruments[instrument].data_ordernames):
                this_index = self.data_info[instrument]['bands'].index(band)
                w_std_min = self.data_info[instrument]['w_std_min'][this_index]
                w_std_max = self.data_info[instrument]['w_std_max'][this_index]
                sigma = self.data_info[instrument]['sigmas'][this_index]
                filename = os.path.split(filename)[1]
                these_stats = dict()
                these_stats['mean'] = Data.getMean(wave=data_wave,\
                                                   flux=data_flux,\
                                                   wmin=w_std_min,\
                                                   wmax=w_std_max)
                these_stats['std'] = Data.getStd(wave=data_wave,\
                                                 flux=data_flux,\
                                                 wmin=w_std_min,\
                                                 wmax=w_std_max)
                these_stats['rms'] = Data.getRMS(wave=data_wave,\
                                                 flux=data_flux-these_stats['mean'],\
                                                 wmin=w_std_min,\
                                                 wmax=w_std_max)
                these_stats['sigma'] = sigma
                #these_stats['tmean'] \
                    #= scipy.stats.tmean(data_flux[scipy.isfinite(data_flux)],\
                                        #limits=(these_stats['mean']\
                                                    #-these_stats['std'],\
                                                #these_stats['mean']\
                                                    #+these_stats['std']))
                #these_stats['tstd'] \
                    #= scipy.stats.tstd(data_flux[scipy.isfinite(data_flux)],\
                                        #limits=(these_stats['mean']\
                                                    #-these_stats['std'],\
                                                #these_stats['mean']\
                                                    #+these_stats['std']))
                print '* Data statistics for %s:'%filename
                print '* Taken between %.2f and %.2f micron.'%(w_std_min,w_std_max)
                print 'mean = ',these_stats['mean'],' Jy'
                print 'std = ',these_stats['std'],' Jy'      
                print 'RMS = ',these_stats['rms'],' Jy'      
                #print 'tmean = ',these_stats['tmean'],' Jy'      
                #print 'tstd = ',these_stats['tstd'],' Jy'    
                stats[filename] = these_stats
        self.data_stats[instrument] = stats
        print '***********************************'
        
        
        
    def getDataInfo(self,instrument):
        
        '''
        Read data info from the ComboCode/Data/Data.dat file. Will return
        information such as sigma levels and wavelength range for determining
        the std value. 
        
        @param instrument: The instrument for which you want this information.
        @type instrument: string
        
        '''
        
        instrument = instrument.upper()
        if self.data_info.has_key(instrument): return
        self.data_info[instrument] = dict()
        indices = [i
                   for i,band in enumerate(DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='INSTRUMENT',filename='Data.dat'))
                   if band.find(instrument) == 0]
        bands = [band.replace('%s_'%instrument,'')
                 for band in DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='INSTRUMENT',filename='Data.dat')
                 if band.find(instrument) == 0]
        w_std_min = [wmin
                     for i,wmin in enumerate(DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='W_STD_MIN',filename='Data.dat'))
                     if i in indices]
        w_std_max = [wmax
                     for i,wmax in enumerate(DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='W_STD_MAX',filename='Data.dat'))
                     if i in indices]
        sigmas = [sigma
                  for i,sigma in enumerate(DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='SIGMA',filename='Data.dat'))
                  if i in indices]
        self.data_info[instrument]['sigmas'] = sigmas
        self.data_info[instrument]['bands'] = bands
        self.data_info[instrument]['w_std_min'] = w_std_min
        self.data_info[instrument]['w_std_max'] = w_std_max
        