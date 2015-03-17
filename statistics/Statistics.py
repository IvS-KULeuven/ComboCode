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
from cc.modeling.objects import Star



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
        self.star_grid = []
        self.instrument = None
        self.path_code = path_code
        self.code = code
        self.data_stats = dict()
        self.data_info = dict()



    def setInstrument(self,instrument_name,data_path='',searchstring='',\
                      data_filenames=[],instrument_instance=None,\
                      pacs_oversampling=None,sample_transitions=[],\
                      redo_convolution=0,resolution=None,absflux_err=0.2,\
                      stat_method='clipping'):
       
        '''
        Set an instrument as a source of data.
        
        @param instrument_name: The instrument (such as 'PACS', 'SPIRE', 
                                or 'FREQ_RESO')
        @type instrument_name: string
        
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
        @keyword oversampling: The instrumental oversampling, for 
                               correct convolution of the Sphinx output.
                               Not required if instrument_instance is given
                                   
                               (default: None)
        @type oversampling: int
        @keyword redo_convolution: if you want to do the convolution of the 
                                   sphinx models for PACS regardless of what's 
                                   already in the database. The pacs id will 
                                   not change, nor the entry in the db, and the
                                   old convolution will be copied to a backup
                                   Only relevant if instrument==PACS. Not 
                                   required if instrument_instance is given
                                 
                                   (default: 0)
        @type redo_convolution: bool
        @keyword resolution: The spectral resolution of the SPIRE apodized 
                             spectra. Only required when instrument_name is 
                             SPIRE.
                             
                             (default: None)
        @type resolution: float
        
        @keyword stat_method: Std/Mean/RMS determination method for unresolved
                              data. Options:
                              - 'clipping': 1-sigma clipping based on std/... 
                                            of full spectrum, then takes 
                                            std/... of clipped spectrum. 
                              - 'preset': wavelength ranges for std/... 
                                          determination taken from Data.dat.
                         
                              (default: 'clipping')
        @type stat_method: string
        @keyword absflux_err: The absolute flux calibration uncertainty of the
                              instrument. Only relevant for PACS or SPIRE and 
                              if no instrument_instance is given.
                               
                              (default: 0.2)
        @type absflux_err: float
        
        '''
        
        instrument_name = instrument_name.upper()
        if instrument_name == 'PACS':
            if instrument_instance <> None:
                self.instrument = instrument_instance
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
            else:
                self.instrument = Pacs.Pacs(star_name=self.star_name,\
                                            oversampling=oversampling,\
                                            path=self.path_code,\
                                            path_pacs=data_path,\
                                            path_combocode=self.path_combocode,\
                                            redo_convolution=redo_convolution,\
                                            absflux_err=absflux_err)
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
        
        elif instrument_name == 'SPIRE':
            if instrument_instance <> None:
                self.instrument = instrument_instance
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
            else:
                self.instrument = Spire.Spire(star_name=self.star_name,\
                                              oversampling=oversampling,\
                                              path=self.path_code,\
                                              path_pacs=data_path,\
                                              resolution=resolution,\
                                              path_combocode=self.path_combocode,\
                                              absflux_err=absflux_err)
                self.instrument.setData(data_filenames=data_filenames,\
                                        searchstring=searchstring)
        
        elif instrument_name == 'FREQ_RESO':
            self.instrument = None
            
        #-- Do data stats for unresolved data. Resolved data will never do this
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

        #-- The unresolved-data case
        if self.instrument:
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
                                           path=path_code,code='GASTRoNOoM',\
                                           path_combocode=self.path_combocode)
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
        if self.instrument:
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
        Read data info from the ComboCode/Data/Data.dat file. Will return
        information such as sigma levels and wavelength range for determining
        the std value. 
        
        @param instrument: The instrument for which you want this information.
        @type instrument: string
        
        '''
        
        if not self.instrument: return
        self.data_info = dict()
        pathcc = os.path.join(self.path_combocode,'Data')
        instbands = DataIO.getInputData(path=pathcc,keyword='INSTRUMENT',\
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
        self.data_info['sigmas'] = sigmas
        self.data_info['bands'] = bands
        self.data_info['w_std_min'] = w_std_min
        self.data_info['w_std_max'] = w_std_max
        