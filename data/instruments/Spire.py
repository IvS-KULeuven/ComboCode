# -*- coding: utf-8 -*-

"""
Producing SPIRE-related output

Author: R. Lombaert

"""

import os
import numpy as np
from scipy import array,sqrt,trapz,log, argmin

from cc.data import Data
from cc.tools.io import Database
from cc.tools.io import DataIO
from cc.data.instruments.Instrument import Instrument



class Spire(Instrument):
    
    """
    Environment with several tools for Spire molecular spectroscopy managing.
    
    """
    
    def __init__(self,star_name,resolution,path_spire,oversampling,\
                 path='codeSep2010',intrinsic=1,path_linefit='',\
                 absflux_err=0.1,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):
        
        '''
        Initializing an instance of Spire().
        
        @param star_name: Name of the star from Star.dat
        @type star_name: string
        @param resolution: The spectral resolution of the SPIRE apodized 
                           spectra
        @type resolution: float
        @param path_spire: full path to SPIRE data folder, excluding star_name
        @type path_spire: string
        @param oversampling: The SPIRE instrumental oversampling, for correct
                             convolution of the Sphinx output. This is the 
                             oversampling factor for the apodized spectra!
        @type oversampling: int
        
        @keyword path: Output folder in the code's home folder
                       
                       (default: 'codeSep2010')
        @type path: string
        @keyword intrinsic: Use the intrinsic Sphinx line profiles for 
                            convolving with the spectral resolution? Otherwise
                            the beam convolved line profiles are used.
                            
                            (default: 1)
        @type intrinsic: bool
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/'
        @type path_combocode: string
        @keyword path_linefit: The folder name for linefit results from Hipe
                               (created by Maries script, assuming her syntax).
                               The folder is located in
                               $path_pacs$/$star_name$/. If no folder is given,
                               no linefits are processed.
                               
                               (default: '')
        @type path_linefit: string
        @keyword absflux_err: The absolute flux calibration uncertainty of the
                               instrument. 
                               
                               (default: 0.1)
        @type absflux_err: float
        
        '''
        
        super(Spire,self).__init__(star_name=star_name,code='GASTRoNOoM',\
                                   path=path,path_combocode=path_combocode,\
                                   path_instrument=path_spire,\
                                   absflux_err=absflux_err,\
                                   oversampling=oversampling,\
                                   path_linefit=path_linefit,\
                                   instrument_name='SPIRE',intrinsic=intrinsic)
        #- resolution is given in cm^-1
        self.resolution = float(resolution)
        self.sigma = self.resolution/(2.*sqrt(2.*log(2.)))
        self.sphinx_convolution = dict()
        if not self.resolution:
            print 'WARNING! SPIRE resolution is undefined!'
        self.readLineFit()

    
    def prepareSphinx(self,star_grid,redo_sphinx_prep=0):
        
        '''
        
        Prepare Sphinx SPIRE models by checking if the convolved spectrum is 
        already present in the Spire object and if not calling convolveSphinx.
        
        @param star_grid: list of Star() instances
        @type star_grid: list[Star()]
        @keyword redo_sphinx_prep: redo the sphinx prep, regardless of whether 
                                   this Spire instance already did it once, 
                                   for instance in case the star_grid changed
                                            
                                   (default: 0)
        @type redo_sphinx_prep: bool

        '''

        if not self.sphinx_convolution or redo_sphinx_prep:
            print 'Convolving with Gaussian and rebinning to data ' + \
                  'wavelength grid.'            
            for i,star in enumerate(star_grid):
                print '* Sphinx model %i out of %i.' %(i+1,len(star_grid))
                if not star['LAST_GASTRONOOM_MODEL']: 
                    print '* No cooling model found.'
                else:
                    self.sphinx_convolution[i] = dict()
                    star['LAST_SPIRE_MODEL'] = i
                    self.__convolveSphinx(star=star)
                    if self.sphinx_convolution[i]:
                        print '* Model %i with cooling id %s is done!'\
                              %(i,star['LAST_GASTRONOOM_MODEL'])
                    else:
                        star['LAST_SPIRE_MODEL'] = None
                      
                      
                      
    def getSphinxConvolution(self,star,fn):
        
        '''
        Return the sphinx convolution and return if it has already been done. 
        
        Returns None if the convolution is not available. 
        
        @param star: The Star() object
        @type star: Star()
        @param fn: The filename of the dataset (band) for which the convolution
                   is to be returned.
        @type fn: str
        
        @return: The sphinx convolution result. (wavelength,flux)
        @rtype: array
        
        '''
        
        if star['LAST_SPIRE_MODEL'] is None:
            return ([],[])
        ifn = self.data_filenames.index(fn)
        return (self.data_wave_list[ifn],\
                self.sphinx_convolution[star['LAST_SPIRE_MODEL']][fn])




    def __convolveSphinx(self,star):
        
        '''
        Convolve the Sphinx output with the SPIRE resolution. The convolution
        is done in wave number (cm^-1).
        
        @param star: The Star() object for which Sphinx profiles are loaded
        @type star: Star()
        
        '''     

        #- Get sphinx model output and merge, for all star models in star_grid
        if not self.resolution: 
            print '* Resolution is undefined. Cannot convolve Sphinx.'
            return
        print '* Reading Sphinx model and merging.'
        sphinx_wav,sphinx_flux = star['LAST_GASTRONOOM_MODEL'] \
                                        and self.mergeSphinx(star) \
                                        or [[],[]]
        if not sphinx_wav: 
            print '* No Sphinx data found.'
            return
        sphinx_wav = 1./array(sphinx_wav)*10**(4)
        sphinx_flux = array(sphinx_flux)
        sphinx_wav = sphinx_wav[::-1]
        sphinx_flux = sphinx_flux[::-1]
        
        #-- eliminate some of the zeroes in the grid to reduce calculation time
        #   (can reduce the array by a factor up to 100!!)
        s = self.sigma
        lcs = array(sorted([1./line.wavelength 
                            for line in star['GAS_LINES']]))
        new_wav, new_flux = [sphinx_wav[0]],[sphinx_flux[0]]
        for w,f in zip(sphinx_wav[1:],sphinx_flux[1:]):
            if f != 0 or (w < 5*s+lcs[argmin(abs(lcs-w))] \
                                and w > lcs[argmin(abs(lcs-w))]-5*s):
                new_wav.append(w)
                new_flux.append(f)
        new_wav, new_flux = array(new_wav), array(new_flux)
        
        #-- convolve the model fluxes with a gaussian and constant sigma(spire)
        print '* Convolving Sphinx model for SPIRE.'
        convolution = Data.convolveArray(new_wav,new_flux,s)
        
        for data_wav,fn in zip(self.data_wave_list,self.data_filenames):
            rebinned = []
            #-- Convert wavelengths to wave number for integration, and reverse
            data_cm = data_wav[::-1]
            data_cm = 1./data_cm*10**4
            rebinned = [trapz(y=convolution[abs(new_wav-wavi)<=self.resolution/self.oversampling],\
                              x=new_wav[abs(new_wav-wavi)<=self.resolution/self.oversampling])\
                            /(self.resolution/self.oversampling)
                        for wavi in data_cm]
            #-- Reverse the rebinned fluxes so they match up with the 
            #   wavelength grid.
            rebinned = array(rebinned)[::-1]
            self.sphinx_convolution[star['LAST_SPIRE_MODEL']][fn] = rebinned



    def readLineFit(self):
        
        '''
        Read the data from the line fit procedure done with Maries Hipe 
        script.
        
        Assumes structure and syntax as given in the example file
        /home/robinl/Data/SPIRE/rdor/lineFit/lineFitResults
        
        The line fit results are saved in self.linefit as a np.recarray.
        
        The relative FWHM is given with respect to the intrinsic SPIRE 
        resolution for unapodized spectra, i.e. 0.04 cm^-1.
        
        The columns include (with unit indicated): 
        band
        freq_in (GHz),
        freq_fit (GHz),
        wave_fit (micron), 
        line_flux (W/m2),
        line_flux_err (W/m2), 
        line_flux_rel, 
        line_peak (W/m2 m), 
        fwhm_fit_freq (GHz),
        fwhm_fit (micron),
        fwhm_rel
        
        '''
        
        dd = super(Spire,self).readLineFit(**dict([('start_row',2)]))
        if dd is None: return
                  
        #-- Remove +/-
        del dd[5] 
        names = ('band','freq_in','freq_fit','wave_fit','line_flux',\
                 'line_flux_err','line_flux_rel','line_peak','fwhm_fit_ghz',\
                 'fwhm_fit', 'fwhm_rel')
        self.linefit = np.recarray(shape=(len(dd[-1]),),\
                                   dtype=zip(names,['|S3']+[float]*10))
        for n,d in zip(names,dd):
            self.linefit[n] = d
