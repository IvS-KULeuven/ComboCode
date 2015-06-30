# -*- coding: utf-8 -*-

"""
Interface for Instrument-specific methods.

Author: R. Lombaert

"""

import os
from glob import glob
from scipy import argmin,argmax,array,sqrt
import scipy

from cc.tools.io import DataIO
from cc.modeling.objects import Star



class Instrument(object):
    
    """
    Instrument-specific calculations.
    
    """
        
    def __init__(self,star_name,path_instrument,instrument_name,oversampling,\
                 absflux_err,\
                 code='GASTRoNOoM',path=None,intrinsic=1,path_linefit='',\
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
        @param oversampling: The instrumental oversampling, for correct
                             convolution of the Sphinx output.
        @type oversampling: int                         
        @param absflux_err: The absolute flux calibration uncertainty of the
                            instrument. 
        @type absflux_err: float
        
        @keyword code: The radiative transfer code used to model the data
        
                       (default: 'GASTRoNOoM')
        @type code: string
        @keyword path: Output folder in the code's home folder. Used to locate 
                       eg PACS database. If None, no model info required (eg 
                       for line measurement matching/identification)
                       
                       (default: None)
        @type path: string                        
        @keyword intrinsic: Use the intrinsic Sphinx line profiles for 
                            convolving with the spectral resolution? Otherwise
                            the beam convolved line profiles are used.
                            
                            (default: 1)
        @type intrinsic: bool
        @keyword path_combocode: CC home folder
        
                                 (default: '/home/robinl/ComboCode/')
        @type path_combocode: string        
        @keyword path_linefit: The folder name for linefit results from Hipe
                               (created by Pierre, assuming his syntax). The 
                               folder is located in $path_pacs$/$star_name$/.
                               If no folder is given, no linefits are 
                               processed.
                               
                               (default: '')
        @type path_linefit: string

        """
        
        self.path = path
        self.path_combocode = path_combocode
        self.code = code
        self.star_name = star_name
        self.path_instrument = path_instrument
        self.path_linefit = path_linefit
        self.instrument = instrument_name.lower()
        self.intrinsic = intrinsic
        self.absflux_err = absflux_err
        self.oversampling = int(oversampling)
        self.data_filenames = []
        ccd = os.path.join(self.path_combocode,'Data')
        istar = DataIO.getInputData(path=ccd,keyword='STAR_NAME').index(star_name)
        #-- Set relevant velocities in cm/s
        self.c = 2.99792458e10 
        self.vlsr = DataIO.getInputData(path=ccd,keyword='V_LSR')[istar]*10**5
        if self.path <> None:
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
                print '** No data filenames or star object given: ' + \
                      'Cannot retrieve data. No %s data will be read.'\
                      %self.instrument.upper()
            else:
                print '** Reading %s data.'%self.instrument.upper()
            self.data_filenames = data_filenames
            self.readData()
            if self.instrument == 'pacs':
                bands = ['B2A','B3A','B2B','R1A','R1B']
            elif self.instrument == 'spire':
                bands = ['SSW','SLW']
            self.data_ordernames = [[word 
                                     for word in os.path.split(f)[1].split('_') 
                                     if word.upper() in bands][0] 
                                    for f in self.data_filenames]
            if len(self.data_ordernames) != len(self.data_filenames): 
                raise IOError('Could not match number of ordernames to ' + \
                              'number of filenames when selecting PACS ' + \
                              'datafiles. Check filenames for missing or ' + \
                              'extra order indications between "_".')

            

    def readData(self):
        
        '''
        Read in data, taking special care of NaNs. 
        
        Four colums are taken as input! wave - contsub - original - continuum
        
        Two columns still works, but may result in errors in other places in 
        the code. 
        
        Data are always read in Jy versus micron, for both SPIRE and PACS.
        
        '''
        
        self.data_wave_list = [] 
        self.data_flux_list = []
        self.data_original_list = [] 
        self.data_continuum_list = [] 
        for filename in self.data_filenames:
            data = DataIO.readCols(filename=filename,nans=1)
            self.data_wave_list.append(data[0])
            self.data_flux_list.append(data[1])
            if len(data) == 2: 
                continue
            self.data_original_list.append(data[2])
            self.data_continuum_list.append(data[3])



    def readLineFit(self,**kwargs):
        
        '''
        Read the data from the line fit procedure.
        
        @keyword kwargs: Extra keywords for the readCols method.
                        
                         (default: dict())
        @type kwargs: dict
        
        @return: The line fit columns are returned.
        @rtype: list[array]
        
        '''
        
        fn = os.path.join(self.path_instrument,self.star_name,\
                          self.path_linefit,'lineFitResults')
        if not self.path_linefit or not os.path.isfile(fn):
            self.linefit = None
            return
        dd = DataIO.readCols(fn,make_array=0,**kwargs)
        return dd
        
    
    
    def intIntMatch(self,trans_list,ifn):
        
        '''
        Match the wavelengths of integrated intensities with transitions.
        
        Checks if a blend might be present based on the data, as well as the 
        available transitions. 
        
        Note that if a line is observed multiple times in the same band (eg in
        a line scan), it cannot at this moment be discerned. In other words,
        there is no point including a line fitted twice in two line scans, as 
        it will not be taken into account for the int matching algorithm.
        
        @param trans_list: The list of transitions for which the check is done.
        @type trans_list: list[Transition()]
        @param ifn: The index of the filename in self.data_filenames. This is 
                    done per file! 
        @type ifn: int
                
        '''
        
        if self.linefit is None:
            return
        dwav = self.data_wave_list[ifn]
        ordername = self.data_ordernames[ifn]
        fn = self.data_filenames[ifn]
        #   1) Prep: Create a list of sample transitions and get central
        #            wavs corrected for vlsr, in micron. Select fit results
        #            for this particular band.
        lf = self.linefit[self.linefit['band'] == ordername]
        #      No info available for band, so don't set anything.
        if not list(lf.wave_fit):
            return
        
        #   2) Collect all transitions with doppler shifted central wavelength
        #      within the wavelength region of the datafile selected here.
        strans = [(t,t.wavelength*10**4*1./(1-self.vlsr/t.c))
                  for t in trans_list 
                  if t.wavelength*10**4*1./(1-self.vlsr/t.c) >= dwav[0]\
                    and t.wavelength*10**4*1./(1-self.vlsr/t.c) <= dwav[-2]]
                
        #   3) Check if the wav of a trans matches a wav in the fitted
        #      intensities list, within the fitted_fwhm of the line with 
        #      respect to the fitted central wavelength, on BOTH sides.
        #      Note that there cannot be 2 fitted lines with central wav
        #      closer than fwhm/2. Those lines would be inseparable!
        #      With the exception of lines superimposed, which should typically
        #      be avoided. Line matching will not be very accurate in this 
        #      case.
        imatches = [argmin(abs(lf.wave_fit-mwav))
                    for (st,mwav) in strans]
        matches = [(mwav <= lf.wave_fit[ii] + lf.fwhm_fit[ii] \
                        and mwav >= lf.wave_fit[ii] - lf.fwhm_fit[ii]) \
                    and (lf.wave_fit[ii],ii) or (None,None)
                   for (st,mwav),ii in zip(strans,imatches)]
        #   4) If match found, check if multiple mtrans fall within   
        #      fitted_FWHM/2 from the fitted central wavelength of the
        #      line. These are blended IN MODEL and/or IN DATA.
        #      Check for model blend is done by looking for ALL transitions 
        #      that have been matched with a single fitted wavelength.
        matches_wv = array([mm[0] for mm in matches])
        wf_blends = [list(array([st for st,mwav in strans])[matches_wv==wv]) 
                     for wv in lf.wave_fit]
        #      Use the wave_fit array indices from matches[:][1] to check  
        #      if indeed multiple transitions were found for the same wav.
        #      If positive, include True if the particular transition is  
        #      the first among the blended ones, else include False. The 
        #      False ones are not taken into account as blended because: 
        #      1) no match was found at all, 
        #      2) only one match was found, 
        #      3) if multiple matches have been found it was not the first. 
        #      
        blended = [ii <> None and len(wf_blends[ii]) > 1. \
                                and wf_blends[ii].index(st) != 0
                   for (st,mwav),(match,ii) in zip(strans,matches)]
        for (st,mwav),blend,(match,ii) in zip(strans,blended,matches):
            #   5) No match found in linefit for this band: no 
            #      integrated intensity is set for this filename in this  
            #      transition. Simply move on to the next transtion. (not 
            #      setting gives None when asking Trans for integrated line 
            #      of unresolved lines)
            if match is None:
                continue
            #   6) Line is blended with other line that is already added. Just
            #      for bookkeeping purposes, all blended lines involved are 
            #      added here as well.
            elif blend:
                st.setIntIntUnresolved(fn,'inblend',None,self.vlsr,wf_blends[ii])
            #   7) Match found with a wave_fit value once. Check for line 
            #      blend IN DATA: Check the ratio fitted FWHM/PACS FWHM. If
            #      larger by 30% or more, put the int int negative. 
            elif len(wf_blends[ii]) == 1:
                err = sqrt((lf.line_flux_rel[ii]/100)**2+self.absflux_err**2)
                factor = lf.fwhm_rel[ii] >= 1.2 and -1 or 1
                st.setIntIntUnresolved(fn,factor*lf.line_flux[ii],err,self.vlsr)
            #   8) If multiple matches, give a selection of strans included
            #      in the blend (so they can be used to select model 
            #      specific transitions later for addition of the 
            #      integrated intensities of the mtrans). Make the integrated 
            #      line strength negative to indicate a line blend. 
            #      The *OTHER* transitions found this way are not compared 
            #      with any data and get None. (see point 4) )
            else: 
                err = sqrt((lf.line_flux_rel[ii]/100)**2+self.absflux_err**2)
                st.setIntIntUnresolved(fn,-1.*lf.line_flux[ii],err,self.vlsr,
                                       wf_blends[ii])
                
                
                
    def mergeSphinx(self,star):
        
        '''
        Merge Sphinx output line profiles on a zero-continuum.
        
        For now only done in wavelength units of micron for PACS/SPIRE spectra.
        
        @param star: The Star object for which all lines are collected + merged
        @type star: Star()
        @return: wave list in micron and flux list in Jy
        @rtype: (list,list)
        
        '''
                
        #read all sphinx output data for the relevant instrument
        sphinx_transitions = [trans 
                              for trans in star['GAS_LINES'] 
                              if trans.getModelId() \
                                and self.instrument.upper() in trans.telescope]
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
        sphinx_final = [(sphinx_input[0][0][0] - 1 + d/1000.,0.0) 
                        for d in xrange(1,1000)]  

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
        #-- Remember overlap_bools contains FALSE if there IS an overlap. 
        #   Counterintuitive, I know. 
        if overlap_bools[-1]:
            sphinx_final.extend([(w,f) 
                                 for w,f in zip(sphinx_input[-1][0],\
                                                sphinx_input[-1][1])])
        sphinx_final.extend([(sphinx_input[-1][0][-1] + d/100.,0.0) 
                             for d in xrange(1,1000)])

        #- No final sort for now, it should be aranged already, even if line 
        #- blends are present
        return [s[0] 
                for s in sphinx_final],\
               [s[1] 
                for s in sphinx_final]
                
        