# -*- coding: utf-8 -*-

"""
Interface for statistics.

Author: R. Lombaert

"""

import os
import scipy.stats
import numpy as np

import cc.path
from cc.tools.io import DataIO
from cc.data.instruments import Pacs
from cc.data import Data
from cc.modeling.objects import Star



class Statistics(object):
    
    """
    An evironment for various statistical calculations for transitions.
    
    """
        
    def __init__(self,star_name,path_code='codeSep2010',code='GASTRoNOoM'):
        
        """ 
        Initializing an instance of Statistics.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2010')
        @type path_code: string
        @keyword code: the code used for producing your output 
        
                       (default: 'GASTRoNOoM')
        @type code: string
        
        """

        self.star_name = star_name
        self.star_grid = []
        self.instrument = None
        self.path_code = path_code
        self.code = code
        self.data_stats = dict()
        self.data_info = dict()



    def setInstrument(self,instrument_name,searchstring='',\
                      data_filenames=[],instrument_instance=None,\
                      sample_transitions=[],\
                      stat_method='clipping',**kwargs):
       
        '''
        Set an instrument as a source of data. Define the instrument_name on 
        calling the method: PACS, SPIRE, SED, FREQ_RESO.
        
        FREQ_RESO requires no additional keywords nor an instrument_instance. 
        
        For the other options, the instrument instance that contains the data 
        references and other info can be passed as an object, or can be set in 
        situ. If the latter, include required keywords for creating the 
        instrument object. 
        
        For PACS: resolution
        For SPIRE: resolution, oversampling
        
        See the respective modules __init__ methods for more information on 
        these keywords as well as other optional keys.
        
        Spire and Pacs objects have to set data manually. Provide the 
        data_filenames or the searchstring for this. Not needed for FREQ_RESO or
        SED. 
        
        @param instrument_name: The instrument (such as 'PACS', 'SPIRE', 'SED',
                                or 'FREQ_RESO')
        @type instrument_name: string
        
        @keyword data_filenames: The data filenames. If empty, auto-search is 
                                 done, with searchstring. Only for Pacs and 
                                 Spire
                                 
                                 (default: [])
        @type data_filenames: list[string]
        @keyword searchstring: the searchstring conditional for the auto-search,
                               only for Pacs and Spire.
        
                               (default: '')
        @type searchstring: string
        @keyword instrument_instance: If None a new instance is made, otherwise
                                      the instance is given here for the 
                                      instrument. Remember to pass the relevant
                                      keywords to the method for the creation.
                                      
                                      (default: None)
        @type instrument_instance: Instrument()
        @keyword stat_method: Std/Mean/RMS determination method for unresolved
                              data (Pacs and Spire) Options:
                                - 'clipping': 1-sigma clipping based on std/... 
                                              of full spectrum, then takes 
                                              std/... of clipped spectrum. 
                                - 'preset': wavelength ranges for std/... 
                                            determination taken from Data.dat.
                         
                              (default: 'clipping')
        @type stat_method: string
        
        '''
        instrument_name = instrument_name.upper()
        if instrument_name == 'PACS':
            if instrument_instance <> None:
                self.instrument = instrument_instance
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
            else:
                self.instrument = Pacs.Pacs(star_name=self.star_name,\
                                            path=self.path_code,**kwargs)
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
        
        elif instrument_name == 'SPIRE':
            if instrument_instance <> None:
                self.instrument = instrument_instance
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
            else:
                self.instrument = Spire.Spire(star_name=self.star_name,\
                                              path=self.path_code,**kwargs)
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
        
        elif instrument_name == 'SED':
            if instrument_instance <> None:
                self.instrument = instrument_instance
            else:
                self.instrument = Sed.Sed(star_name=self.star_name,**kwargs)
        
        elif instrument_name == 'FREQ_RESO':
            self.instrument = None
            
        #-- Do data stats for unresolved data. Resolved or SED data will never 
        #   do this
        self.doDataStats(method=stat_method)
        
        
        
    def setModels(self,star_grid=[],models=[]):
        
        '''
        Load the models and remember them.
        
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

        #-- The SED case
        if self.instrument and self.instrument.instrument == 'SED':
            if not star_grid: 
                raise IOError('Statistics.setModels requires a ' + \
                              'star_grid to be defined for SED data.')
            if set([s['LAST_MCMAX_MODEL'] for s in star_grid]) == set(['']):
                return
            #-- Get rid of models that were not calculated successfully
            self.star_grid = [s for s in star_grid if s['LAST_MCMAX_MODEL']]
            self.star_grid = np.array(self.star_grid)
            
        #-- The unresolved-data case
        elif self.instrument:
            instr = self.instrument.instrument.upper()
            print '***********************************'
            print '* Checking Sphinx models for comparison with %s.'%instr
            if instr == 'SPIRE': models = []
            if not star_grid and not models:
                raise IOError('Statistics.setModels requires either ' + \
                              'star_grid or models to be defined.')
            #-- star_grid can be reconstructed for PACS through its database
            elif not star_grid:
                star_grid = Star.makeStars(models=models,id_type='PACS',\
                                           path=path_code,code=self.code)
                self.instrument.addStarPars(star_grid)
            if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) == set([0]): 
                return
            self.instrument.prepareSphinx(star_grid)
            self.star_grid = [star 
                              for star in star_grid 
                              if star['LAST_%s_MODEL'%instr] <> None]
        #-- The FREQ_RESO case
        else:
            if not star_grid:
                raise IOError('Statistics.setModels requires a ' + \
                              'star_grid to be defined for freq-resolved data.')
            if set([s['MOLECULE'] and 1 or 0 for s in star_grid]) == set([0]): 
                return
            self.star_grid = star_grid
            
            
    
    def doDataStats(self,method='clipping'):
        
        '''
        Calculate mean, std, rms for the data files.
        
        They are saved in self.data_stats.
        
        @keyword method: Std/Mean/RMS determination method. Can be 'clipping'
                         or 'preset'. The former 1-sigma clips based on 
                         mean/std/... of full spectrum, then takes mean/std/...
                         of clipped spectrum. Latter takes preset values from
                         Data.dat.
                         
                         (default: 'clipping')
        @type method: string
        
        '''
        
        method = method.lower()
        if method not in ['preset','clipping']: method = 'clipping'
        self.data_stats = dict()
        if self.instrument \
                and self.instrument.instrument.upper() in ['PACS','SPIRE']:
            if not self.data_info: self.setDataInfo()
            for filename,data_wave,data_flux,band in \
                        zip(self.instrument.data_filenames,\
                            self.instrument.data_wave_list,\
                            self.instrument.data_flux_list,\
                            self.instrument.data_ordernames):
                these_stats = dict()
                iband = self.data_info['bands'].index(band)
                these_stats['sigma'] = self.data_info['sigmas'][iband]
                if method == 'preset':
                    w_std_min = self.data_info['w_std_min'][iband]
                    w_std_max = self.data_info['w_std_max'][iband]
                    if w_std_min < data_wave[0] or w_std_max > data_wave[-1]:
                        method = 'clipping'
                    test_mean = Data.getMean(wave=data_wave,flux=data_flux,\
                                             wmin=w_std_min,wmax=w_std_max)
                    if test_mean is None:
                        method = 'clipping'
                if method == 'clipping':
                    totmean = Data.getMean(wave=data_wave,flux=data_flux)
                    totstd = Data.getStd(wave=data_wave,flux=data_flux)
                    limits = (-totstd,totstd)
                    these_stats['mean'] = Data.getMean(flux=data_flux,\
                                                       limits=limits)
                    these_stats['std'] = Data.getStd(flux=data_flux,\
                                                     limits=limits)
                    these_stats['rms'] = Data.getRMS(flux=data_flux-\
                                                          these_stats['mean'],\
                                                     limits=limits)                    
                else:
                    these_stats['mean'] = Data.getMean(wave=data_wave,\
                                                       flux=data_flux,\
                                                       wmin=w_std_min,\
                                                       wmax=w_std_max) 
                    these_stats['std'] = Data.getStd(wave=data_wave,\
                                                     flux=data_flux,\
                                                     wmin=w_std_min,\
                                                     wmax=w_std_max)
                    these_stats['rms'] = Data.getRMS(wave=data_wave,\
                                                     flux=data_flux-\
                                                          these_stats['mean'],\
                                                     wmin=w_std_min,\
                                                     wmax=w_std_max)
                print '* Data statistics for %s using "%s" method:'\
                      %(filename,method)
                if method == 'preset':
                    print '* Taken between %.2f and %.2f micron.'\
                          %(w_std_min,w_std_max)
                print 'mean = ',these_stats['mean'],' Jy'
                print 'std = ',these_stats['std'],' Jy'      
                print 'RMS = ',these_stats['rms'],' Jy'      
                self.data_stats[filename] = these_stats
            print '***********************************'
        
        
        
    def setDataInfo(self):
        
        '''
        Read data info from the ComboCode/usr/Data.dat file. Will return
        information such as sigma levels and wavelength range for determining
        the std value. 
        
        '''
        
        if not self.instrument: return
        self.data_info = dict()
        instbands = DataIO.getInputData(keyword='INSTRUMENT',\
                                        filename='Data.dat')
        instbands = [band.upper() for band in instbands]
        indices = [i
                   for i,band in enumerate(instbands)
                   if band.find(self.instrument.instrument.upper()) == 0]
        bands = [band.replace('%s_'%self.instrument.instrument.upper(),'')
                 for i,band in enumerate(instbands)
                 if i in indices]
        w_std_min = [wmin
                     for i,wmin in enumerate(DataIO.getInputData(\
                                      keyword='W_STD_MIN',filename='Data.dat'))
                     if i in indices]
        w_std_max = [wmax
                     for i,wmax in enumerate(DataIO.getInputData(\
                                      keyword='W_STD_MAX',filename='Data.dat'))
                     if i in indices]
        sigmas = [sigma
                  for i,sigma in enumerate(DataIO.getInputData(\
                                          keyword='SIGMA',filename='Data.dat'))
                  if i in indices]
        self.data_info['sigmas'] = sigmas
        self.data_info['bands'] = bands
        self.data_info['w_std_min'] = w_std_min
        self.data_info['w_std_max'] = w_std_max
        