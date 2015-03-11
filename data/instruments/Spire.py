# -*- coding: utf-8 -*-

"""
Producing SPIRE-related output

Author: R. Lombaert

"""

import os
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
                 path='codeSep2010',\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),intrinsic=1):
        
        '''
        Initializing an instance of Spire().
        
        @param star_name: Name of the star from Star.dat
        @type star_name: string
        @param resolution: The spectral resolution of the SPIRE data
        @type resolution: float
        @param path_spire: full path to PACS data folder, excluding star_name
        @type path_spire: string
        
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
        
        '''
        
        super(Spire,self).__init__(star_name=star_name,code='GASTRoNOoM',\
                                   path=path,path_combocode=path_combocode,\
                                   path_instrument=path_spire,\
                                   instrument_name='SPIRE',intrinsic=intrinsic)
        #- resolution is given in cm^-1
        self.resolution = float(resolution)
        self.sigma = self.resolution/(2.*sqrt(2.*log(2.)))
        self.sphinx_convolution = dict()
        self.oversampling = int(oversampling)
        if not self.resolution:
            print 'WARNING! SPIRE resolution is undefined!'


    
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
                    self.sphinx_convolution[(star['LAST_GASTRONOOM_MODEL'],i)]\
                          = dict()
                    self.__convolveSphinx(star=star,istar=i)
                    if self.sphinx_convolution:
                        print '* Model %i with id %s is done!'\
                              %(i,star['LAST_GASTRONOOM_MODEL'])
                      


    def __convolveSphinx(self,star,istar):
        
        '''
        Convolve the Sphinx output with the SPIRE resolution. The convolution
        is done in wave number (cm^-1).
        
        @param star: The Star() object for which Sphinx profiles are loaded
        @type star: Star()
        @param istar: The index of the Star() object
        @type istar: int
        
        '''     

        #- Get sphinx model output and merge, for all star models in star_grid
        print '* Reading Sphinx model and merging.'
        sphinx_wav,sphinx_flux = star['LAST_GASTRONOOM_MODEL'] \
                                        and self.mergeSphinx(star) \
                                        or [[],[]]
        if not sphinx_wav: 
            print '* No Sphinx data found.'
            return
        if not self.resolution: 
            print '* Resolution is undefined. Cannot convolve Sphinx.'
            return
        sphinx_wav = 1./array(sphinx_wav)*10**(4)
        sphinx_flux = array(sphinx_flux)
        sphinx_wav = sphinx_wav[::-1]
        sphinx_flux = sphinx_flux[::-1]
        DataIO.writeCols('/home/mariev/test1.dat',[sphinx_wav,sphinx_flux])
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
        DataIO.writeCols('/home/mariev/test2.dat',[new_wav,new_flux])
        #-- convolve the model fluxes with a gaussian and constant sigma(spire)
        convolution = Data.convolveArray(new_wav,new_flux,s)
        DataIO.writeCols('/home/mariev/test3.dat',[new_wav,convolution])
        
        for data_wav,fn in zip(self.data_wave_list,self.data_filenames):
            print '* Convolving Sphinx model between %.2f and %.2f cm^-1.'\
                  %(data_wav[0],data_wav[-1])
            rebinned = []
            [rebinned.append(\
                trapz(y=convolution[abs(new_wav-wavi)<=self.resolution/self.oversampling],\
                      x=new_wav[abs(new_wav-wavi)<=self.resolution/self.oversampling])\
                /(self.resolution/self.oversampling))
             for wavi in data_wav]
            self.sphinx_convolution[(star['LAST_GASTRONOOM_MODEL'],istar)][fn]\
                        = rebinned
            DataIO.writeCols('/home/mariev/test4_%s.dat'%(fn),[data_wav,rebinned])
    
