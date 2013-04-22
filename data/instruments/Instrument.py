# -*- coding: utf-8 -*-

"""
Interface for Instrument-specific methods.

Author: R. Lombaert

"""

import os
from glob import glob
from scipy import argmax,array
import scipy

from cc.tools.io import DataIO
from cc.modeling.objects import Star



class Instrument(object):
    
    """
    Instrument-specific calculations.
    
    """
        
    def __init__(self,star_name,path_instrument,instrument_name,\
                 code='GASTRoNOoM',path='codeSep2010',intrinsic=1,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):        
        
        """ 
        Initializing an instance of Instrument.
        
        @param star_name: The name of the star from Star.py
        @type star_name: string
        @param path_instrument: The data folder for this instrument. Data are
                                  located in ~/<path_instrument>/<star_name>/
                                  cont_subtracted/
                                  If None
        @type path_instrument: string                          
        @param instrument_name: The name of the instrument (SPIRE or PACS)
        @type instrument_name: string
                                  
        @keyword code: The radiative transfer code used to model the data
        
                       (default: 'GASTRoNOoM')
        @type code: string
        @keyword path: Output folder in the code's home folder
        
                       (default: 'codeSep2010')
        @type path: string                                
        @keyword intrinsic: Use the intrinsic Sphinx line profiles for 
                            convolving with the spectral resolution? Otherwise
                            the beam convolved line profiles are used.
                            
                            (default: 1)
        @type intrinsic: bool
        @keyword path_combocode: CC home folder
        
                                 (default: '/home/robinl/ComboCode/')
        @type path_combocode: string        
        
        """
        
        self.path = path
        self.path_combocode = path_combocode
        self.code = code
        self.star_name = star_name
        self.path_instrument = path_instrument
        self.instrument = instrument_name.upper()
        self.intrinsic = intrinsic
        self.data_filenames = []
        ccd = os.path.join(self.path_combocode,'Data')
        istar = DataIO.getInputData(path=ccd,keyword='STAR_NAME').index(star_name)
        #-- Set relevant velocities in cm/s
        self.c = 2.99792458e10 
        self.vlsr = DataIO.getInputData(path=ccd,keyword='V_LSR')[istar]*10**5
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                self.code,self.path,'stars'))
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                self.code,self.path,'stars',\
                                                self.star_name))
                                            


    def makeModelList(self,star_grid,id_type):
        
        '''
        Return a list of all model id's in the star_grid.
                
        @param star_grid: The grid of parameter sets
        @type star_grid: list[Star()]
        @param id_type: the type of model id (MCMAX or GASTRONOOM or PACS)
        @type id_type: string
        @return: the model_ids
        @rtype: list[string]
        
        '''
        
        return [star['LAST_'+id_type.upper()+'_MODEL'] for star in star_grid]



    def getDataFilenames(self,searchstring=''):
        
        '''
        Retrieve the data filenames that you wish to be plotted.
        
        @keyword searchstring: the searchstring conditional for selecting data 
                               files.
                                      
                               (default: '')
        @type searchstring: string
        @return: the data filenames found in the search
        @rtype: list[string]
        
        '''
        
        return [f 
                for f in glob(os.path.join(self.path_instrument,\
                                           self.star_name,'cont_subtracted',\
                                           '*'+searchstring+'*'))
                if f[-1] != '~' and f[-7:] != '.tar.gz']
                     


    def setData(self,data_filenames=[],searchstring=''):
        
        '''
        Read data from given filenames and remember them.
        
        The data are saved as wavelength and flux lists in the object.
        
        @keyword data_filenames: The data filenames. If empty, auto-search is 
                                 done
                                 
                                 (default: [])
        @type data_filenames: list[string]
        @keyword searchstring: the searchstring conditional for the auto-search
        
                               (default: '')
        @type searchstring: string
        
        '''
        if not self.data_filenames:
            if not data_filenames:
                data_filenames = self.getDataFilenames(searchstring=\
                                                       searchstring)
            if not data_filenames:
                print 'No data filenames or star object given: ' + \
                      'Cannot retrieve data. No data will be set.'
                return
            self.data_filenames = data_filenames
            print '** Reading %s data.'%self.instrument
            self.readData(data_filenames)
            
            

    def readData(self,data_filenames):
        
        '''
        Read in data, taking special care of NaNs. 
        
        Two colums are taken as input!
        
        @param data_filenames: list of filenames to read
        @type data_filenames: list[string]
        
        '''
        
        self.data_wave_list = [DataIO.readCols(filename=filename,nans=1)[0] 
                               for filename in data_filenames]
        self.data_flux_list = [DataIO.readCols(filename=filename,nans=1)[1] 
                               for filename in data_filenames]

            
            
    def mergeSphinx(self,star):
        
        '''
        Merge Sphinx output line profiles on a zero-continuum.
        
        For now only done in wavelength units of micron for PACS/SPIRE spectra.
        
        @param star: The Star object for which all lines are collected + merged
        @type star: Star()
        @return: wave list in micron and flux list in Jy
        @rtype: (list,list)
        
        '''
                
        #read all sphinx output data
        sphinx_transitions = [trans 
                              for trans in star['GAS_LINES'] 
                              if trans.getModelId()]
        #- If no sphinx output found, this list will be empty and no convolution
        #- should be done. 
        if not sphinx_transitions: 
            return [[],[]]
        
        [trans.readSphinx() for trans in sphinx_transitions]
        sphinx_input = [self.intrinsic \
                            and (trans.sphinx.getVelocityIntrinsic(),\
                                 trans.sphinx.getLPIntrinsic())
                            or (trans.sphinx.getVelocity,\
                                trans.sphinx.getLPConvolved())
                        for trans in sphinx_transitions]
        
        #- convert km/s to cm/s to micron and flux to Jy 
        #- doppler shift (1-(v_source - v_observer=delta_v)/c)*f_zero converted 
        #- to wavelength in micron
        sphinx_input = [(1/(1.-(vel*10**5/star.c))*trans.wavelength*10**(4),\
                         flux*10**(23)) 
                        for (vel,flux),trans in zip(sphinx_input,\
                                                    sphinx_transitions)]
        
        #-- In case the wavelength/freq scale is counting down, reverse arrays
        sphinx_input = [wav[0] > wav[-1] \
                            and (wav[::-1],flux[::-1]) or (wav,flux) 
                        for wav,flux in sphinx_input]
        sphinx_input = [(list(wav),list(flux)) for wav,flux in sphinx_input]    

        #- Make sure all sphinx segments are increasing in wavelength/frequency
        sphinx_input.sort()
        
        #- Check if overlap between lines is present: 
        #- False if THIS line is overlapping with the one before it. 
        #- Hence the first line will always be true.
        overlap_bools = [1] + \
                        [sphinx_input[i-1][0][-1] <= sphinx_input[i][0][0] 
                         for i in xrange(1,len(sphinx_input))]
        
        if False in overlap_bools:
            print 'WARNING! There is overlap between emission lines in ' + \
                  'Sphinx output. Overlap is included by simple addition only!'
        
        #- Add zeroes on a grid before the first line to make sure the 
        #- convolution goes right
        sphinx_final = [(sphinx_input[0][0][0] - 0.1 + d/1000.,0.0) 
                        for d in xrange(1,100)]  
        
        #- For now stitch up the segments with zeroes. Later on: 
        #- Make grid with zeroes, then add segments if no overlap, 
        #- if overlap, make sure the addition works well
        for i in xrange(len(sphinx_input)-1):
            #- if True, no overlap with the sphinx file before i and it has to 
            #- be added
            #- if False, overlap with the sphinx file before i and it has been 
            #- added already
            if overlap_bools[i]:  
                j = 1
                this_wave = sphinx_input[i][0]
                this_flux = sphinx_input[i][1]
                blend_wave = []
                blend_flux = []
                #- If False, overlap: interpolation is needed, select all the 
                #- sphinx results for which there is a subsequent blend
                while not i+j == len(sphinx_input)-1 \
                        and not overlap_bools[i+j]:
                    blend_wave.append(sphinx_input[i+j][0])
                    blend_flux.append(sphinx_input[i+j][1])
                    j += 1
                
                #- note that at i+j there is no more blend, and it is a 
                #- discrete sphinx result: This should not be added yet
                #- if there are blends here, blend_wave is not []
                if blend_wave:
                    #- remember where the largest wavelength is to be found
                    max_index = argmax(array([wav[-1] 
                                              for wav in blend_wave]))+i+1
                    if sphinx_input[max_index][0][-1] < sphinx_input[i][0][-1]: 
                        max_index = i
                    #- Making a new wavelength grid beyond the current line at 
                    #- i; there should be no need to expand the grid
                    #- to lower wavelengths, since sphinx_final was sorted. 
                    delta_lambda = (this_wave[-1]-this_wave[0])/len(this_wave)
                    #- expand until all wavelenghts in all blends are covered 
                    while this_wave[-1] < max([wav[-1] for wav in blend_wave]):
                        this_wave.append(this_wave[-1]+delta_lambda)
                        this_flux.append(0.0)
                    interpolations = []
                    #- interpolate to get values for the new grid, where the 
                    #- original wavelength grid for every blend is expanded
                    #- such that all wavelengths in the new grid are covered, 
                    #- adding zeroes in this way does not change the result
                    for x,y in zip(blend_wave,blend_flux):
                        this_x = array([this_wave[0],x[0]-(x[1]-x[0])] + x + \
                                       [x[-1]+(x[-1]-x[-2]),this_wave[-1]])
                        this_y = array([0,0] + y + [0,0])
                        interpolations.append(array(scipy.interpolate.interp1d\
                                                   (this_x,this_y)(this_wave)))
                    #- add together 
                    this_flux=array(this_flux)
                    for blend in interpolations: this_flux += blend
                else:
                    max_index = i
                
                #- extend result with the new wavelengths and fluxes, 
                #- including blends if needed
                sphinx_final.extend([(x,y) 
                                     for x,y in zip(this_wave,this_flux)])
                #- extend result further with zeroes up until the next sphinx 
                #- line, after the last line in the blend 
                #- Take delta equal to average delta in last Sphinx segment, 
                #- this may be different from delta in next segment, 
                #- but since they're zeroes you don't really care 
                sphinx_final.append((2*sphinx_final[-1][0]\
                                        -sphinx_final[-2][0],\
                                     0.0))
                while sphinx_final[-1][0] + 0.001 < sphinx_input[i+j][0][0]:
                    sphinx_final.append((sphinx_final[-1][0]+0.001,0.0))
                sphinx_final.append((2*sphinx_input[i+j][0][0]-\
                                        sphinx_input[i+j][0][1],\
                                     0.0))
        sphinx_final.extend([(sphinx_input[-1][0][-1] + d/1000.,0.0) 
                             for d in xrange(1,100)])
        if not overlap_bools[-1]:
            sphinx_final.extend([(w,f) 
                                 for w,f in zip(sphinx_input[-1][0],\
                                                sphinx_input[-1][1])])
        #- No final sort for now, it should be aranged already, even if line 
        #- blends are present
        return [s[0] 
                for s in sphinx_final],\
               [s[1] 
                for s in sphinx_final]
                
        