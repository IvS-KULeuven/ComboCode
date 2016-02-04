# -*- coding: utf-8 -*-

"""
Performing statistics on the peak and integrated intensity of resolved 
molecular lines.

Author: R. Lombaert

"""

import os
import scipy
from scipy import argmin,ones
from scipy import array
from scipy import sum
import operator
import types
import numpy as np

import cc.path
from cc.tools.io import DataIO
from cc.statistics.Statistics import Statistics
import matplotlib.pyplot as plt
from cc.modeling.objects import Transition
from cc.statistics import BasicStats    
from time import gmtime

class ResoStats(Statistics):
    
    """
    Environment with several tools to perform statistics on the peak and 
    integrated intensities of resolved molecular lines.
    
    """
        
    def __init__(self,star_name,code='GASTRoNOoM',path_code='codeJun2013',\
                 lll_p=None):        
        
        """ 
        Initializing an instance of IntIntStats.
        
        Then run setInstrument, setModels and setIntensities. The rest of the 
        class works interactively.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
        @keyword code: the code used for producing your output 
        
                       (default: 'GASTRoNOoM')
        @type code: string
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2010')
        @type path_code: string
        @keyword lll_p: The number of variable parameters in the model grid. If 
                        None, the automatic loglikelihood determination is not 
                        done. You can still set the threshold values with a 
                        class bound method.
        
                        (default: None)
        @type lll_p: int
        
        """
        
        super(ResoStats,self).__init__(star_name=star_name,\
                                       code=code,path_code=path_code)
        
        #-- List of template transitions
        self.translist = []
        #-- Dict of lists of models (value) for each template transition (key)
        self.trans_models = dict()
        #-- Similarly, keep track of each Star() model that is valid.
        #   This way one can keep track of those models that have a successful
        #   cooling result
        self.star_selection = dict()
        
        #-- Dicts keeping integrated and peak intensities of datafiles. 
        #   key: the sample Transition(), 
        #   value: list of values for first datafile included
        #   Note that for every model the data value is recalculated, based on 
        #   the best vlsr guess which depends on the model.
        #   As a result, every element in the value list corresponds to the 
        #   Star() object with the same index in self.star_grid.
        self.dinttmb = dict()
        self.dpeaktmb = dict()
        #-- Dicts keeping integrated and peak intensities of sphinx models, as 
        #   well as the loglikelihood measure of fitness between model and data
        #   The peak and integrated intensity ratio dicts follow same pattern
        #   key: the sample Transition(), 
        #   value: list of values for every model included
        #   Every element in the value list corresponds to the Star() object 
        #   with the same index in self.star_grid.
        self.minttmb = dict()
        self.mpeaktmb = dict()
        self.loglikelihood = dict()
        self.ratiopeak = dict()
        self.ratioint = dict()
        self.ratiocombo = dict()
        #-- Only set to True if something failed somewhere. Likely not yet 
        #   implemented/resolved issues.
        self.no_stats = False
        self.stats = dict([('peak',self.ratiopeak),('int',self.ratioint),\
                           ('combo',self.ratiocombo)])
        #-- Default Flux calibration uncertainties from Telescope.dat
        #   I suggest to change these values based on what you want to use 
        #   yourself, depending on the transition under consideration.
        telescopes = DataIO.getInputData(keyword='TELESCOPE',start_index=5,\
                                         filename='Telescope.dat')
        abs_errs = DataIO.getInputData(keyword='ABS_ERR',start_index=5,\
                                       filename='Telescope.dat')
        self.tele_uncertainties = dict(zip(telescopes,abs_errs))
        
        #-- The uncertainties and loglikelihoods for each template transition 
        #   is kept here, but the key is the INDEX of the template transition.
        self.trans_uncertainties = dict()
        self.lll_threshold = dict()
        #-- A dict keeping track of noisy data. It is possible to set this to
        #   False anyway, if one want to use the dataset regardless of it being
        #   too noisy. Structure as eg lll_threshold
        self.noisy = dict()
        
        #-- Set the number of variables in the model grid, needed for automatic
        #   loglikelihood threshold calculation
        if lll_p is None or lll_p == '' or int(lll_p) < 1 or int(lll_p) > 20:
            self.lll_p = None
            print 'No STAT_LLL_P given. Not calculating the loglikelihood '+\
                  'threshold automatically.'
        else: 
            self.lll_p = int(lll_p)
        
        
    def setInstrument(self,sample_transitions):
       
        '''
        Set and read the data objects for this statistics module. 
        
        In this case, a list of sample transitions in which the data will be 
        read is the only input required.
        
        Make sure sample_transitions are copies of the originals, such as 
        returned by Transition.extractTransFromStars().
    
        @param sample_transitions: Sample transitions used as reference for 
                                   the data files. 
        @type sample_transitions: list[Transition()]
        
        '''
        
        super(ResoStats,self).setInstrument(instrument_name='FREQ_RESO')
        #-- Make copy so that any changes do not translate to whatever the 
        #   original list of transitions might be. If transitions are 
        #   selected through Transition.extractTransFromStars, this is 
        #   taken care of.
        if not sample_transitions:
            print 'WARNING! No sample transitions given for Statistics ' +\
                    'module with instrument == FREQ_RESO. No stats will ' + \
                    'be calculated.'
        [t.readData() for t in sample_transitions]
        self.sample_trans = sample_transitions    

        
    def setIntensities(self,use_bestvlsr=1, partial = 0, vcut = 0):
        """
        The data intensities are stored in the dictionary 
        self.dintint/self.dpeakint, the model intensities in 
        self.mintint/self.mpeakint. 
        
        They keys in the dictionaries are the Transition to which the sphinx 
        or datafiles belong. 
        
        The model values are lists in accordance with self.star_grid, the data
        values are single floats for the first dataset for the transition.

        In addition, the loglikelihoods are calculated.

        @keyword use_bestvlsr: Use the fitted best-guess for the v_lsr when 
                               determining the velocity grid for the model. If 
                               not, the vlsr from the Star.dat file or the fits
                               file is used. 
                               
                               (default: 1)
        @type use_bestvlsr: bool
        
        """
        self.vcut = vcut
        self.partial = partial
        
        if partial != 0:
            print 'Calculating loglikelihood statistics using a partial line profile.'
            print 'Cutoff velocity = '+str(vcut)
        
        
        tnodata = [t 
                   for t in self.sample_trans
                   if not t.lpdata]
        self.translist = [t 
                          for t in self.sample_trans
                          if t.lpdata]
        self.includedtrans = [i for i in range(len(self.translist))]

        #- For every sample transition (st), collect the equivalent transitions
        #- in the model grid. Then retrieve all integrated and peak tmb values,
        #- for both data and model. 
        for ist,st in enumerate(self.translist):
            #-- make sure the noise value is set in the data object.
            noise = st.getNoise()
            
            #-- Check which default uncertainty is needed for this transition
            for k,v in self.tele_uncertainties.items(): 
                if k in st.telescope: 
                    self.trans_uncertainties[ist] = v
                    continue
            self.lll_threshold[ist] = None
            
            #-- Collect all models for this transition that have a successful
            #   cooling subcode result. If not successful, they will be
            #   excluded everywhere.
            self.trans_models[st]  = [star.getTransition(st)
                                      for star in self.star_grid
                                      if star['LAST_GASTRONOOM_MODEL']]
            
            #-- Keep track of which Star() models are valid for each sample 
            #   transition
            self.star_selection[st] = [star
                                       for star in self.star_grid
                                       if star['LAST_GASTRONOOM_MODEL']]
            
            #-- If None's are still in the list of results, it means either 
            #   mline or sphinx failed, while cooling didn't. I simply did not
            #   yet take this into account as a possibility. TBI.
            all_ids = [bool(t.getModelId()) for t in self.trans_models[st]]
            if False in all_ids: 
                self.no_stats = True
                print 'One of the sphinx models was not calculated properly,'+\
                      ' even though cooling succeeded. Cannot do statistics.'
                return
                
            #-- Set all data according to the template transition, to avoid too
            #   much overhead reading, setting and fitting data.
            for mt in self.trans_models[st]: 
                mt.setData(st)
                
            #-- Collect the model integrated and peak Tmbs
            self.minttmb[st] = array([mt.getIntTmbSphinx() 
                                      for mt in self.trans_models[st]])
            self.mpeaktmb[st] = array([mt.getPeakTmbSphinx() 
                                       for mt in self.trans_models[st]])
            
            #-- Set the data integrated and peak Tmb for this dataset
            self.dinttmb[st],abs_err = st.getIntTmbData()
            self.dpeaktmb[st] = st.getPeakTmbData() 
            
            if self.dpeaktmb[st] <= 3*noise:
                self.noisy[ist] = True
            else: 
                self.noisy[ist] = False
            
            #-- Collect the loglikelihoods for all models
            self.loglikelihood[st] = array([mt.getLoglikelihood(use_bestvlsr, \
                                                partial = partial, vcut = vcut) \
                                                for mt in self.trans_models[st]])
            
            #-- Calculate the ratios for integrated and peak Tmbs (model/data)
            self.ratioint[st] = self.minttmb[st]/self.dinttmb[st]
            self.ratiopeak[st] = self.mpeaktmb[st]/self.dpeaktmb[st]
            self.ratiocombo[st] = zip(self.ratiopeak[st],self.ratioint[st])

        self.calcLoglikelihoodThreshold()
    
    
    def calcLoglikelihoodThreshold(self,bfms=[],ist=None):
        
        '''
        Determine the loglikelihood thresholds from the model grid, assuming 
        the amount of free parameters in the model grid is known. 
        
        The latter is given by STAT_LLL_P in the CC input.
        
        If a selection of models is given, and the index of a sample transition
        the max lll of the selection of models and threshold value based on 
        that --- for this sample trans --- are returned. No values are 
        remembered in this case.
        
        @keyword bfms: Subselection of models for which the LLL threshold is 
                       calculated. If None, ist is ignored, and the LLL 
                       threshold is calculated for the whole collection of 
                       models, the value for which is saved in 
                       self.lll_threshold
                       
                       (default: [])
        @type bfms: list[Star]
        @keyword ist: If bfms is not an empty list, this gives the index of the 
                      sample trans for which the lll threshold is calculated 
                      based on subgrid of Star() models, given by bfms. If 
                      bfms == [], ist is ignored.
        
                      (default: None)
        @type ist: int
        
        @return: The max lll and lll threshold for this sample transistion and 
                 subset of models. None if no subset of models was given.
        @rtype: (float,float)
        
        '''
        
        #-- Determine the loglikelihood 95% quantile of the Chi^2_p 
        #   distribution, with p degrees of freedom. 
        #   See Decin et al 2007 Table 2.
        quantiles = dict([(1,3.8414588),(2,5.9914645),(3,7.8147279),\
                          (4,9.4877290),(5,11.070498),(6,12.591587),\
                          (7,14.067140),(8,15.507313),(9,16.918978),\
                          (10,18.307038),(11,19.675138),(12,21.026070),\
                          (13,22.362032),(14,23.684791),(15,24.995790),\
                          (16,26.296228),(17,27.587112),(18,28.869299),\
                          (19,30.143527),(20,31.410433)])
        if self.lll_p is None:
            return 
        quant = quantiles[self.lll_p]
        self.lll_quant = quant
        if not bfms:
            for ist,st in enumerate(self.translist):
                #-- See Decin et al. 2007 section 3.1.3.
                maxlll = max(self.loglikelihood[st])
                self.lll_threshold[ist] = -quant/2.+maxlll
            return
        else: 
            st = self.translist[ist]
            maxlll = max([lll
                          for lll,smodel in zip(self.loglikelihood[st],\
                                                self.star_selection[st]) 
                          if smodel in bfms])
            lll_threshold = -quant/2.+maxlll
            return (maxlll,lll_threshold)
            
  
    def setNoisy(self,ist,isnoisy):
        
        '''
        Set a transition given by its index to be included or excluded in the
        proper best fit determination, through either integrated or peak Tmbs.
        The other option is to merely look at the 3 sigma level compared to the
        model peak value. 
        
        By default, the script checks if dtmb is above 3 sigma in the dataset.
        
        Sometimes if dtmb is between 2 sigma and 3 sigma, one may want to 
        include it properly anyway. Also if the peak Tmb determination is not 
        very accurate just because of the noisy data, it may be necessary to 
        include it regardless. This is possible here.
        
        Note that the state of the 'noisiness' has to be included when calling
        this method. False or True works, of course.
        
        @param ist: The index/indices of the transitions to be included.
        @type ist: int/list
        @param isnoisy: Is the dataset noisy or not?
        @type isnoisy: bool
        
        '''
        
        isnoisy = int(isnoisy)
        self.noisy[ist] = isnoisy



    def resetNoisy(self, factor, use_fit = True):
        '''
        Change criterion for noisy lines.  By default, the script checks if
        the peak of the line profile is above 3 sigma in the dataset. 
        With this method, the factor can be changed. 
        
        With use_fit, it is possible to use the gaussian or parabolic fit to 
        the data to derive the integrated intensity instead of the line itself 
        for noisy lines. 
        
        @keyword factor: Sigma level. If the peak intensity of a line is below
                         this level, it is marked as noisy.
        @type factor: int
        
        @keyword use_fit: Use gaussian or parabolic fit instead of the line
                          to derive the integrated intensity of noisy lines
        @type use_fit: bool
        '''
        
        for ist,st in enumerate(self.translist):
            noise = st.getNoise()
            if self.dpeaktmb[st] <= factor*noise:
                self.noisy[ist] = True
            else:
                self.noisy[ist] = False

        
            if use_fit == True:
                if self.noisy[ist] == False:
                    self.dinttmb[st],abs_err = st.getIntTmbData(use_fit = 1)
                else:
                    self.dinttmb[st],abs_err = st.getIntTmbData(use_fit = 0)
                self.ratioint[st] = self.minttmb[st]/self.dinttmb[st]
            else:
                self.dinttmb[st],abs_err = st.getIntTmbData(use_fit = 0)
                

        
    def includeTrans(self,ist):
        
        '''
        Include a transition in selecting the best fit models. 
        
        By default, all transitions that have data and sphinx models associated
        with them are selected.
        
        Use the index of the transition to include them.
        
        @param ist: The index/indices of the transitions to be included.
        @type ist: int/list
        
        '''
        
        if type(ist) is types.IntType:
            ist = [ist]
        self.includedtrans.extend(ist)
        self.includedtrans = sorted(list(set(self.includedtrans)))
        
        
        
    def excludeTrans(self,ist):
        
        '''
        Exclude a transition in selecting the best fit models. 
        
        By default, all transitions that have data and sphinx models associated
        with them are selected.
        
        Use the index of the transition to exclude them.
        
        @param ist: The index/indices of the transitions to be included.
        @type ist: int/list
        
        '''
        
        if type(ist) is types.IntType:
            ist = [ist]
        self.includedtrans = [i for i in self.includedtrans if i not in ist]
    
    
    
    def listTrans(self):
        
        '''
        List all the transitions that have both data and sphinx models 
        associated with them.
        
        Can be used to figure the index of each transition so you can in- or 
        exclude them.
        
        Also indicates whether it is a noisy line or not. 
        '''
        
        strings = ['%i [%scluded]: %s \t at %.2f uncertainty. The line is %snoisy!%s'\
                    %(ist,ist in self.includedtrans and 'IN' or 'EX',\
                      DataIO.fillOutSpaces(str(st),60),\
                      self.trans_uncertainties[ist],\
                      (not self.noisy[ist] and 'NOT ' or ''),\
                      self.lll_threshold[ist] <> None and \
                        ' The lll threshold is %.2e'%self.lll_threshold[ist] or\
                        '')
                   for ist,st in enumerate(self.translist)]
        print '\n'.join(strings)
        
        
        
    def setUncertainty(self,ist,value):
        
        '''
        Change the default uncertainty for a given transition.
        
        Makes use of the template transition index! (see listTrans method)
        
        The value is given as a decimal number, i.e. x% / 100.
        
        @param ist: The index of transition 
        @type ist: int
        @param value: The new uncertainty requested 
        @type value: float
        
        '''
        
        self.trans_uncertainties[ist] = float(value)
        
    
    
    def setLoglikelihoodThreshold(self,ist,value):
        
        '''
        Set the loglikelihood below which models are considered best fit models
        
        Makes use of the template transition index!
        
        An automatically calculated value is determined as well from the lll 
        values of all models. This assumes the amount of variable parameters
        in the model grid is known, and set through the CC input! (STAT_LLL_P)
        
        @param ist: The index of transition 
        @type ist: int
        @param value: The new loglikelihood threshold requested 
        @type value: float
        
        '''
        
        self.lll_threshold[ist] = float(value)
        
        
    
    def setLLL(self,ist,value):
        
        '''
        Convenience method that runs setLogelikelihoodThreshold().
        
        @param ist: The index of transition 
        @type ist: int
        @param value: The new loglikelihood threshold requested 
        @type value: float
        
        '''
        
        self.setLoglikelihoodThreshold(ist,value)
            
            
    def printStats(self,bfms=[]):
          

        '''
        Print the statistics for all transitions and models. 
        
        @keyword bfms: A list of a subselection of Star() models. If specified
                       only the models in this list will be printed.
                       
                       (default: [])
        @type bfms: list[Star]
        
        '''
        
        bfms = list(bfms)
        if self.no_stats: return
        for ist,st in enumerate(self.translist):
            if ist in self.includedtrans:
            
                #-- No need for the vexp value here, we know the noise is already
                #   set.
                noise = st.getNoise()
                print '*******************************************************'
                print 'Statistics for %i: %s:'%(ist,str(st))
                print '-------------------------------------'
                print 'Data intensities [integrated --- peak --- STD (noise)]:'
                print '%f K km/s \t---\t %f K \t---\t %f K'\
                        %(self.dinttmb[st],self.dpeaktmb[st],noise)
                print '-------------------------------------'
                print 'Model/Data intensities [integrated --- peak --- lll]:'
                lines = ['- %s: \t %.3f \t---\t %.3f \t---\t %.2e\
                            '%(str(tr.getModelId()),ri,rp,lll)
                        for tr,ri,rp,lll,s in zip(self.trans_models[st],\
                                                self.ratioint[st],\
                                                self.ratiopeak[st],\
                                                self.loglikelihood[st],\
                                                self.star_selection[st])
                        if (not bfms or s in bfms)]
                print '\n'.join(lines)
                if not bfms:
                    print 'The max LLL for this model is : %.2e.'\
                        %max(self.loglikelihood[st])
                    if self.lll_p <> None: 
                        print 'This leads to a threshold of %.2e with %i free parameters.'\
                          %(self.lll_threshold[ist],self.lll_p)            
                else:
                    if self.lll_p <> None:
                        maxlll, lll_thr = self.calcLoglikelihoodThreshold(bfms,ist)
                        print 'The max LLL for this subselection of models is: ' +\
                            '%.2e.' %maxlll
                        print 'This leads to a threshold of %.2e with %i free ' \
                            %(lll_thr,self.lll_p) + 'parameters.'
                    else:
                        maxlll = max([lll
                                    for lll,smodel in zip(self.loglikelihood[st],\
                                                        self.star_selection[st]) 
                                    if smodel in bfms])
                        print 'The max LLL for this subselection of models is: ' +\
                            '%.2e.' %maxlll
                
            
        
    def selectBestFitModels(self,mode='int',use_lll=1,output=1):
            
        '''
        Returns a list of Star() models for which the uncertainty criteria
        are satisfied for the selected transitions.
        
        Transitions can be selected or deselected by the includeTrans and
        excludeTrans methods. 
        
        The uncertainty criteria are preset but can be changed by the 
        setUncertainty method.
        
        Note that for lines in which the peak Tmb is less than 3*sigma, the 
        goodness of fit is determined based on the comparison between peak of
        sphinx model and the 3*sigma level. Less: included, more: excluded.
        Note that the sphinx peak values are scaled down by the uncertainty on 
        the data, to take into account data uncertainties as well. 
        
        @keyword mode: The mode for model selection. Include: int, peak, combo
                       int: based on the integrated intensity ratios, 
                       peak: based on the peak intensity ratios
                       combo: Based on both
                       
                       (default: 'int')
        @type mode: string
        @keyword use_lll: Also use the loglikelihood statistics in selecting 
                          best fit models. For this it is required to set a 
                          threshold loglikelihood value: Anything smaller than 
                          this is included in the best fit model list. Use
                          the setLoglikelihoodThreshold() or setLLL()
                          
                          (default: 1)
        @type use_lll: bool
        
        @return: The 'best fitting' Star() models according to the 
                 uncertainty criteria. 
        @rtype: list[Star()]
        
        '''
        
        if self.no_stats: return
        if output:
            print 'Selecting best fit models in <%s> mode.'%mode

        self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId() \
            for ii in range(len(self.star_grid))]
        stars = array(self.modellist)

        bfbools = ones(len(stars),dtype='bool')
        for ist,st in enumerate(self.translist):
            if ist in self.includedtrans:
                dpeak = self.dpeaktmb[st]
                noise = st.getNoise()
                err = self.trans_uncertainties[ist]
                lll_thresh = self.lll_threshold[ist]
                for i,(rat,lll,mt) in enumerate(zip(self.stats[mode][st],\
                                                    self.loglikelihood[st],\
                                                    self.trans_models[st])):
                    if self.noisy[ist]:
                        this_bool = mt.getPeakTmbSphinx()*(1.-err) <= 3.*noise
                        if not this_bool: 
                            bfbools[i] = False
                    elif mode == 'combo':
                        this_bool = (rat[0] < 1.+err and rat[0] > 1.-err) and\
                                    (rat[1] < 1.+err and rat[1] > 1.-err)
                        if not this_bool: bfbools[i] = False
                    else:
                        this_bool = rat < 1.+err and rat > 1.-err
                        if not this_bool: 
                            bfbools[i] = False
                            
                    ##-- Loglikelihood is maximized by best fitting model
                    if lll_thresh <> None and not self.noisy[ist] \
                            and use_lll and lll < lll_thresh:
                        bfbools[i] = False
                        
        self.bfm = stars[bfbools]
        self.bfm = list(self.bfm)
        return self.bfm
        
    
    
    def selectBestFitModelsLLL(self, output = 1):
        '''
        Same function as selectBestFitModels, but uses only the 
        loglikelihood statistic.
        
        @return: The 'best fitting' Star() models according to the 
                 loglikelihood statistic. 
        @rtype: list[Star()]        
        '''
        if self.no_stats: return
        if output:
            print 'Selecting best fit models, only based on loglikelihood statistic.'    
        
        self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId() \
            for ii in range(len(self.star_grid))]
        stars = array(self.modellist)

        bfbools = ones(len(stars),dtype='bool')
        for ist,st in enumerate(self.translist):
            if ist in self.includedtrans:
                lll_thresh = self.lll_threshold[ist]
                for i,lll in enumerate(self.loglikelihood[st]):
                    ##-- Loglikelihood is maximized by best fitting model
                    #if lll_thresh <> None and not self.noisy[ist] \
                            #and lll < lll_thresh:
                        #bfbools[i] = False
                    if lll_thresh <> None  \
                            and lll < lll_thresh:
                        bfbools[i] = False
        self.bfmlll = stars[bfbools]      
        self.bfmlll = list(self.bfmlll)
        
        return self.bfmlll    


            
    def findLLLRange(self, includedTrans = 1):
        '''
        Determines whether a model lies within a 95% confidence interval around 
        the best fitting model (per line). As the best fitting model is included
        in the arrays, at least one model within the range is to be expected.
        
        Based on Decin et al. 2007: l_m \leq l \leq l_m - quant/2
                                    l_m \leq l \leq lll_threshold
        
        
        @return confidenceLLL: dictionary containing the difference between the 
        calculate loglikelihood and the threshold value
        @rtype: dict(list[])
        
        @return confidenceLLL_verdict: dictionary containing the result (1 or 0) 
        for every model per transition
        @rtype: dict(list[])
        
        @return confidenceLLL_models: models that fit all included transitions
        @rtype: array([])
        '''
        
        self.confidenceLLL = dict()
        self.confidenceLLL_verdict = dict()
        
        for i,st in enumerate(self.translist):
            self.confidenceLLL[st] = []
            self.confidenceLLL_verdict[st] = []
            
            maxlll = max(self.loglikelihood[st])
            
            for kk in range(len(self.loglikelihood[st])):
                if maxlll >= self.loglikelihood[st][kk] \
                   and self.loglikelihood[st][kk] >= self.lll_threshold[i]:
                    self.confidenceLLL[st].append(self.loglikelihood[st][kk]-self.lll_threshold[i])
                    self.confidenceLLL_verdict[st].append(1)
                else:
                    self.confidenceLLL[st].append(self.loglikelihood[st][kk]-self.lll_threshold[i])
                    self.confidenceLLL_verdict[st].append(0)
        
        total = []
        if includedTrans == 1:
            includedTrans = self.includedtrans
        
        for itr in includedTrans:
            for i in range(len(self.star_grid)):
                if self.confidenceLLL_verdict[self.translist[itr]][i] == 1:
                    total.append(i)
        
        counted = [total.count(k) for k in set(total)]
        self.confidenceLLL_models = np.where(array(counted) == max(counted))[0]




    def selectBestFitperLine(self, use_lll = 1, excludevib = 1):
        '''
        Checks which lines are fitted by which models, using all four selection
        criteria. For every criterion, a list 'ocurrences_bfmxxx' is made. This 
        list contains a sublist for each model, len(list) = # models. Every sublist
        contains a one or zero, indicating wether a line is fitted by the model 
        according to the criterion, len(sublist) = # lines.
        
        The lists fittedLines_bfmxxx contain how often a line was modelled according
        to the criterion, len(list) = # lines.
        
        The lists fittedModels_bfmxxx contain the number of lines each model fitted
        according to the criterion, len(list) = # models.
        
        '''
        if self.no_stats: return
        print "Checking which model occurs most as a best fit per line..."
        
        trans = self.translist
        T = len(trans)
        orig_included = self.includedtrans
              
        #-- Excluded vibrational transitions, if necessary      
        if excludevib:
            self.excludeTrans([i for i in range(T) if str(trans[i]).split()[1] != '0'])
        
        #-- Perform selection routine for every line separately
        all_bfmint = []
        all_bfmpeak = []
        all_bfmcombo = []
        all_bfmlll = []    

        
        for ii in orig_included:
            self.excludeTrans([i for i in range(T)])
            self.includeTrans(ii)
            all_bfmint.append(self.selectBestFitModels('int', use_lll, output = 0))
            all_bfmpeak.append(self.selectBestFitModels('peak', use_lll, output = 0))
            all_bfmcombo.append(self.selectBestFitModels('combo', use_lll, output = 0))
            all_bfmlll.append(self.selectBestFitModelsLLL(output = 0))

        #-- Make the main list
        self.occurences_bfmint = [[all_bfmint[x].count(self.modellist[i]) \
            for x in range(len(all_bfmint))] \
            for i in range(len(self.modellist))]
        
        self.occurences_bfmpeak = [[all_bfmpeak[x].count(self.modellist[i]) \
            for x in range(len(all_bfmpeak))] \
            for i in range(len(self.modellist))]    
        
        self.occurences_bfmcombo = [[all_bfmcombo[x].count(self.modellist[i]) \
            for x in range(len(all_bfmcombo))] \
            for i in range(len(self.modellist))] 
        
        self.occurences_bfmlll = [[all_bfmlll[x].count(self.modellist[i]) \
            for x in range(len(all_bfmlll))] \
            for i in range(len(self.modellist))]

        #-- Set included_trans back to original state
        self.includedtrans = orig_included

        #-- How often is a line fitted by a model?
        self.fittedLines_bfmint = sum(self.occurences_bfmint, axis = 0)
        self.fittedLines_bfmpeak = sum(self.occurences_bfmpeak, axis = 0)
        self.fittedLines_bfmcombo = sum(self.occurences_bfmcombo, axis = 0)
        self.fittedLines_bfmlll = sum(self.occurences_bfmlll, axis = 0)
        
        #-- How many lines did a model fit?
        self.fittedModels_bfmint = sum(self.occurences_bfmint, axis = 1)
        self.fittedModels_bfmpeak = sum(self.occurences_bfmpeak, axis = 1)
        self.fittedModels_bfmcombo = sum(self.occurences_bfmcombo, axis = 1)
        self.fittedModels_bfmlll = sum(self.occurences_bfmlll, axis = 1)        
        

    def selectBestFitperLineLLL(self, excludevib = 1, plot = 0):
        '''
        Checks which lines are fitted by which models according to the loglikelihood
        statistic. A  list 'ocurrences_bfmlll' is made. This 
        list contains a sublist for each model, len(list) = # models. Every sublist
        contains a one or zero, indicating wether a line is fitted by the model 
        according to the criterion, len(sublist) = # lines.
        
        The list fittedLines_bfmlll contains how often a line was modelled according
        to the criterion, len(list) = # lines.
        
        The list fittedModels_bfmlll contains the number of lines each model fitted
        according to the criterion, len(list) = # models.
        
        '''
        if self.no_stats: return
        print "Checking which model occurs most as a best fit per line..."
        
        trans = self.translist
        T = len(trans)
        orig_included = self.includedtrans
              
        #-- Excluded vibrational transitions, if necessary      
        if excludevib:
            self.excludeTrans([i for i in range(T) if str(trans[i]).split()[1] != '0'])
        
        #-- Perform selection routine for every line separately
        all_bfmlll = []    

        
        for ii in orig_included:
            self.excludeTrans([i for i in range(T)])
            self.includeTrans(ii)
            all_bfmlll.append(self.selectBestFitModelsLLL(output = 0))

        #-- Make the main list        
        self.occurences_bfmlll = [[all_bfmlll[x].count(self.modellist[i]) \
            for x in range(len(all_bfmlll))] \
            for i in range(len(self.modellist))]

        #-- Set included_trans back to original state
        self.includedtrans = orig_included

        #-- How often is a line fitted by a model?
        self.fittedLines_bfmlll = sum(self.occurences_bfmlll, axis = 0)
        
        #-- How many lines did a model fit?
        self.model_bfmlll = sum(self.occurences_bfmlll, axis = 1)        
    
        if plot:
            plot_id = 'plot_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i' \
                %(gmtime()[0],gmtime()[1],gmtime()[2],\
                    gmtime()[3],gmtime()[4],gmtime()[5])
            self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId().replace('_','-') \
                for ii in range(len(self.star_grid))]
            
            plt.clf()
            fig = plt.figure(2, figsize = (15, 10))
            ax = fig.add_subplot(111)
            ax.set_xticks(np.arange(len(self.translist))-0.5)
            ax.set_yticks(np.arange(len(self.modellist)))
            ax.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope if self.noisy[ist]== False\
                else '\\textbf{'+str(st.jup)+'-'+str(st.jlow)+' '+st.telescope+'}' for ist,st in enumerate(self.translist)], rotation = 60)
            ax.yaxis.set_ticklabels([i for i in self.modellist])
            ax.imshow(self.occurences_bfmlll, interpolation='nearest', origin='upper', cmap = 'Blues')
            plt.tight_layout()
            fig.suptitle('Loglikehood criterium', size = 18)
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename_combo = os.path.join(path, 'resostats','lll-%s_len_%s-%s'%(self.modellist[0],(len(self.modellist)),plot_id))
            fig.savefig(filename_combo+'.pdf')   
            print '*** Plot of stats can be found at:'
            print filename_combo+'.pdf'
  
  
    
    def calcLLL(self, plot = 0):
        
        self.line_lll = dict()
        self.line_lll_range = dict()
        
        self.model_lll = []
        self.model_lll_range = []
        
        translist = [self.translist[i] for i in self.includedtrans]
        
        
        for ist,st in enumerate(translist):
            self.line_lll[st] = []
            self.line_lll_range[st] = []
            
            maxlll = max(self.loglikelihood[st])
            
            for kk in range(len(self.loglikelihood[st])):
                if self.loglikelihood[st][kk] >= self.lll_threshold[ist]:
                    self.line_lll[st].append(1)
                else:
                    self.line_lll[st].append(0)
                
                if maxlll >= self.loglikelihood[st][kk] \
                   and self.loglikelihood[st][kk] >= self.lll_threshold[ist]:
                    self.line_lll_range[st].append(1)
                else:
                    self.line_lll_range[st].append(0)
            
        for ii in range(len(self.star_grid)):
            self.model_lll.append([self.line_lll[tr][ii] for tr in translist])
            self.model_lll_range.append([self.line_lll_range[tr][ii] for tr in translist])
        
        self.verdict_model_lll = sum(self.model_lll, axis = 1)        
        self.verdict_model_lll_range = sum(self.model_lll_range, axis = 1)
        
        if plot:
            plot_id = 'plot_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i' \
                %(gmtime()[0],gmtime()[1],gmtime()[2],\
                    gmtime()[3],gmtime()[4],gmtime()[5])
            self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId().replace('_','-') \
                for ii in range(len(self.star_grid))]
            
            plt.clf()
            fig = plt.figure(1, figsize = (15, 10))
            ax1 = fig.add_subplot(121)
            ax1.set_xticks(np.arange(len(translist))-0.5)
            ax1.set_yticks(np.arange(len(self.modellist)))
            ax1.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope if self.noisy[ist]== False\
                else '\\textbf{'+str(st.jup)+'-'+str(st.jlow)+' '+st.telescope+'}' for ist,st in enumerate(translist)], rotation = 60)
            ax1.yaxis.set_ticklabels([i for i in self.modellist])
            ax1.imshow(self.model_lll, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax1.set_title('LLL criterion')

            
            
            ax2 = fig.add_subplot(122)
            ax2.set_xticks(np.arange(len(translist))-0.5)
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope if self.noisy[ist]== False\
                else '\\textbf{'+str(st.jup)+'-'+str(st.jlow)+' '+st.telescope+'}' for ist,st in enumerate(translist)], rotation = 60)
            ax2.yaxis.set_ticklabels([i for i in self.modellist])
            ax2.imshow(self.model_lll_range, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_title('LLL 95\% confidence interval')
            
            plt.tight_layout()
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename = os.path.join(path, 'resostats','LLL+range-%s_len_%s-%s'%(self.modellist[0],(len(self.modellist)),plot_id))
            fig.savefig(filename+'.pdf')   
            print '*** Plot of stats can be found at:'
            print filename+'.pdf'
  
        
        
    
    def calcRatioIntTmb(self, useNoisy = 1, useRms = 0, err = 0.2, err_noisy = 0.3, plot = 0):
        '''
        Calculate the ratio of the integrated main beam intensity per transition.
        
        @keyword useNoisy: assign a larger error to noisy lines
        @type useNoisy: bool
        
        @keyword useRms: use statistical noise in addition to instrumental error
        @type useRms: bool
        
        @keyword err: error on data
        @type err: float
        
        @keyword err_noisy: error on noisy data (only needed when useNoisy = 1)
        @type err_noisy: float
        
        @keyword plot: Plot the output
        @type plot: bool
        
        @return self.verdict_ratioint: dictionary containing whether a line satisfies 
                                       the condition or not. One list per transition,
                                       each list containing the verdict per model.
        @rtype verdictpermodel_int: dict(list[])
        
        @return self.model_ratioint: list containing the verdict per transition, per model.
                                     Number of lists = number of models. Length of each 
                                     list = number of transitions.                                    
        @rtype self.model_ratioint: list[list[]]
        
        @return self.model_ratioint_verdict: list containing the total number of transitions
                                             that satisfy the condition, per model. Number 
                                             of lists = number of models. Each list contains
                                             a single number.
        '''
        
        print 'Calculating ratios of integrated main beam intensities...'
        print 'Error on data = '+str(err)
        if useNoisy:
            print 'Error on noisy data = '+str(err_noisy)
        
        self.line_ratioint = dict()
        
        translist = [self.translist[i] for i in self.includedtrans]
        
        for ist, st in enumerate(translist):
            if useNoisy == 1 and useRms == 1:
                if self.noisy[ist] ==  True:
                    self.line_ratioint[st] = [1 if x <= (1. + err_noisy + st.getNoise()) \
                        and x >= (1. - err_noisy - st.getNoise()) else 0 for x in self.ratioint[st]]
                else:
                    self.line_ratioint[st] = [1 if x <= (1. + err + st.getNoise()) \
                        and x >= (1. - err - st.getNoise()) else 0 for x in self.ratioint[st]]
        
            elif useNoisy == 1 and useRms == 0:
                if self.noisy[ist] ==  True:
                    self.line_ratioint[st] = [1 if x <= (1. + err_noisy) \
                        and x >= (1. - err_noisy) else 0 for x in self.ratioint[st]]
                else:
                    self.line_ratioint[st] = [1 if x <= 1. + err \
                        and x >= 1. - err else 0 for x in self.ratioint[st]]
        
            elif useNoisy == 0 and useRms == 1:    
                self.line_ratioint[st] = [1 if x <= (1. + err + st.getNoise()) \
                    and x >= (1. - err - st.getNoise()) else 0 for x in self.ratioint[st]]
            
            else:
                self.line_ratioint[st] = [1 if x <= 1. + err \
                    and x >= 1. - err else 0 for x in self.ratioint[st]]
        
        
        self.model_ratioint = []
        
        for ii in range(len(self.star_grid)):
            self.model_ratioint.append([self.line_ratioint[tr][ii] for tr in translist])
        
        self.verdict_model_ratioint = sum(array(self.model_ratioint), axis = 1)

        
        if plot:
            plot_id = 'plot_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i' \
                %(gmtime()[0],gmtime()[1],gmtime()[2],\
                    gmtime()[3],gmtime()[4],gmtime()[5])
            self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId().replace('_','-') \
                for ii in range(len(self.star_grid))]
            
            plt.clf()
            fig = plt.figure(2, figsize = (15, 10))
            ax = fig.add_subplot(111)
            ax.set_xticks(np.arange(len(translist))-0.5)
            ax.set_yticks(np.arange(len(self.modellist)))
            ax.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope if self.noisy[ist]== False\
                else '\\textbf{'+str(st.jup)+'-'+str(st.jlow)+' '+st.telescope+'}' for ist,st in enumerate(translist)], rotation = 60)
            ax.yaxis.set_ticklabels([i for i in self.modellist])
            ax.imshow(self.model_ratioint, interpolation='nearest', origin='upper', cmap = 'Blues')
            plt.tight_layout()
            fig.suptitle('Ratio of integrated intensities \\ Error = '+str(err*100.)+'\%, error noisy lines = '+str(err_noisy*100.)+'\%', size = 18)
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename_combo = os.path.join(path, 'resostats','int-%s_len_%s-%s'%(self.modellist[0],(len(self.modellist)),plot_id))
            fig.savefig(filename_combo+'.pdf')   
            print '*** Plot of stats can be found at:'
            print filename_combo+'.pdf'
            
    
    def combineRatioIntLLL(self, useNoisy = 1, useRms = 0, err = 0.2, err_noisy = 0.3, plot = 0):
        
        self.calcRatioIntTmb(useNoisy=useNoisy,useRms=useRms,err=err,err_noisy=err_noisy)
        self.combRatioIntLLL = []
        self.calcLLL()
        
        translist = [self.translist[i] for i in self.includedtrans]
        
        for ii in range(len(self.star_grid)):
            self.combRatioIntLLL.append([self.model_lll[ii][x] \
                if self.model_lll[ii][x] == self.model_ratioint[ii][x] else 0 \
                for x in range(len(translist))])
        
        self.verdict_combRatioIntLLL = sum(array(self.combRatioIntLLL), axis = 1)
        
        if plot:
            plt.clf()
            plot_id = 'plot_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i' \
                %(gmtime()[0],gmtime()[1],gmtime()[2],\
                    gmtime()[3],gmtime()[4],gmtime()[5])
            self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId().replace('_','-') \
                for ii in range(len(self.star_grid))]
            
            fig = plt.figure(2, figsize = (15, 10))
            ax1 = fig.add_subplot(131)
            ax1.set_xticks(np.arange(len(translist))-0.5)
            ax1.set_yticks(np.arange(len(self.modellist)))
            ax1.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope if self.noisy[ist]== False\
                else '\\textbf{'+str(st.jup)+'-'+str(st.jlow)+' '+st.telescope+'}' for ist,st in enumerate(translist)], rotation = 60)
            ax1.yaxis.set_ticklabels([i for i in self.modellist])
            ax1.imshow(self.combRatioIntLLL, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax1.set_title('Combination')
            
            ax2 = fig.add_subplot(132)
            ax2.set_xticks(np.arange(len(translist))-0.5)
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope if self.noisy[ist]== False\
                else '\\textbf{'+str(st.jup)+'-'+str(st.jlow)+' '+st.telescope+'}' for ist,st in enumerate(translist)], rotation = 60)
            ax2.yaxis.set_ticklabels([i for i in self.modellist])
            ax2.imshow(self.model_ratioint, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_title('Ratio integrated intensities')
            
            ax3 = fig.add_subplot(133)
            ax3.set_xticks(np.arange(len(translist))-0.5)
            ax3.set_yticks(np.arange(len(self.modellist)))
            ax3.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope if self.noisy[ist]== False\
                else '\\textbf{'+str(st.jup)+'-'+str(st.jlow)+' '+st.telescope+'}' for ist,st in enumerate(translist)], rotation = 60)
            ax3.yaxis.set_ticklabels([i for i in self.modellist])
            ax3.imshow(self.model_lll, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax3.set_title('Loglikelihood')
            
            plt.tight_layout()
            fig.suptitle('Models complying to integrated intensity and loglikelihood criteria \\ Error = '+str(err*100.)+'\%, error noisy lines = '+str(err_noisy*100.)+'\%', size = 18)
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename = os.path.join(path, 'resostats','intLLL-%s_len_%s_partial_%s_vcut_%s-%s'\
                %(self.modellist[0],(len(self.modellist)),self.partial,self.vcut,plot_id))
            fig.savefig(filename+'.pdf')   
            print '*** Plot of stats can be found at:'
            print filename+'.pdf'
        
        
        
    def makeInput(self, model_array):
        
        for m in model_array:
            print 'MOLECULE='+(' ').join(self.star_grid[m]['MOLECULE'][0])
            
        
        
    
    
    def calcRatioPeakTmb(self, useNoisy = 1, useRms = 0, err = 0.2, err_noisy = 0.3, plot = 0):
        '''
        Calculate the ratio of the peak main beam intensity per transition.
        
        @keyword useNoisy: assign a larger error to noisy lines
        @type useNoisy: bool
        
        @keyword useRms: use statistical noise in addition to instrumental error
        @type useRms: bool
        
        @keyword err: error on data
        @type err: float
        
        @keyword err_noisy: error on noisy data (only needed when useNoisy = 1)
        @type err_noisy: float
        
        @return self.verdict_ratiopeak: dictionary containing whether a line satisfies 
                                       the condition or not. One list per transition,
                                       each list containing the verdict per model.
        @rtype verdictpermodel_peak: dict(list[])
        
        @return self.model_ratiopeak: list containing the verdict per transition, per model.
                                     Number of lists = number of models. Length of each 
                                     list = number of transitions.                                    
        @rtype self.model_ratiopeak: list[list[]]
        
        @return self.model_ratiopeak_verdict: list containing the total number of transitions
                                             that satisfy the condition, per model. Number 
                                             of lists = number of models. Each list contains
                                             a single number.
        '''
        
        print 'Calculating ratios of peak main beam intensities...'
        print 'Error on data = '+str(err)
        if useNoisy:
            print 'Error on noisy data = '+str(err_noisy)
            
        self.verdict_ratiopeak = dict()
        
        for ist, st in enumerate(self.translist):
            if useNoisy == 1 and useRms == 1:
                if self.noisy[ist] ==  True:
                    self.verdict_ratiopeak[st] = [1 if x <= (1. + err_noisy + st.getNoise()) \
                        and x >= (1. - err_noisy - st.getNoise()) else 0 for x in self.ratiopeak[st]]
                else:
                    self.verdict_ratiopeak[st] = [1 if x <= (1. + err + st.getNoise()) \
                        and x >= (1. - err - st.getNoise()) else 0 for x in self.ratiopeak[st]]
        
            elif useNoisy == 1 and useRms == 0:
                if self.noisy[ist] ==  True:
                    self.verdict_ratiopeak[st] = [1 if x <= (1. + err_noisy) \
                        and x >= (1. - err_noisy) else 0 for x in self.ratiopeak[st]]
                else:
                    self.verdict_ratiopeak[st] = [1 if x <= 1. + err \
                        and x >= 1. - err else 0 for x in self.ratiopeak[st]]
        
            elif useNoisy == 0 and useRms == 1:    
                self.verdict_ratiopeak[st] = [1 if x <= (1. + err + st.getNoise()) \
                    and x >= (1. - err - st.getNoise()) else 0 for x in self.ratiopeak[st]]
            
            else:
                self.verdict_ratiopeak[st] = [1 if x <= 1. + err \
                    and x >= 1. - err else 0 for x in self.ratiopeak[st]]
        
        
        self.model_ratiopeak = []
        
        for ii in range(len(self.star_grid)):
            self.model_ratiopeak.append([self.verdict_ratiopeak[tr][ii] for tr in self.translist])
        
        self.model_ratiopeak_verdict = sum(array(self.model_ratiopeak), axis = 1)    
        
        
        if plot:
            self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId() \
                for ii in range(len(self.star_grid))]
            tr = [self.verdict_ratiopeak[st] for st in self.translist]
            
            plt.clf()
            fig = plt.figure(1, figsize = (15, 10))
            ax = fig.add_subplot(111)
            ax.set_xticks(np.arange(len(tr))-0.5)
            ax.set_yticks(np.arange(len(self.modellist)))
            ax.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope for st in self.translist], rotation = 60)
            ax.yaxis.set_ticklabels([i for i in self.modellist])
            ax.imshow(self.model_ratiopeak, interpolation='nearest', origin='upper', cmap = 'Blues')
            plt.tight_layout()
            fig.suptitle('Ratio of peak intensities', size = 18)
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename_combo = os.path.join(path, 'resostats','int-%s_len_%s-%s'%(self.modellist[0],(len(self.modellist)),plot_id))
            fig.savefig(filename_combo+'.pdf')   
            print '*** Plot of stats can be found at:'
            print filename_combo+'.pdf'





    def calcRatioComboTmb(self, useNoisy = 1, useRms = 0, err = 0.2, err_noisy = 0.3, plot = 0):
        '''
        Combine the ratio of the integrated and peak main beam intensity per transition.
        
        @keyword useNoisy: assign a larger error to noisy lines
        @type useNoisy: bool
        
        @keyword useRms: use statistical noise in addition to instrumental error
        @type useRms: bool
        
        @keyword err: error on data
        @type err: float
        
        @keyword err_noisy: error on noisy data (only needed when useNoisy = 1)
        @type err_noisy: float
        
        @return self.verdict_ratiocombo: dictionary containing whether a line satisfies 
                                       the condition or not. One list per transition,
                                       each list containing the verdict per model.
        @rtype verdictpermodel_combo: dict(list[])
        
        @return self.model_ratiocombo: list containing the verdict per transition, per model.
                                     Number of lists = number of models. Length of each 
                                     list = number of transitions.                                    
        @rtype self.model_ratiocombo: list[list[]]
        
        @return self.model_ratiocombo_verdict: list containing the total number of transitions
                                             that satisfy the condition, per model. Number 
                                             of lists = number of models. Each list contains
                                             a single number.
        '''
        
        print 'Calculating ratios of main beam intensities, combining int and peak... '
        print 'Error on data = '+str(err)
        if useNoisy:
            print 'Error on noisy data = '+str(err_noisy)
            
        self.verdict_ratiocombo = dict()
        self.y = dict()
        
        for ist, st in enumerate(self.translist):
            if useNoisy == 1 and useRms == 1:
                if self.noisy[ist] ==  True:
                    self.y[st] = [[1 if x <= (1. + err_noisy + st.getNoise()) \
                        and x >= (1. - err_noisy - st.getNoise()) else 0 for x in self.ratiocombo[st][i]] \
                        for i in range(len(self.ratiocombo[st]))]
                else:
                    self.y[st] = [[1 if x <= (1. + err + st.getNoise()) \
                        and x >= (1. - err - st.getNoise()) else 0 for x in self.ratiocombo[st][i]] \
                        for i in range(len(self.ratiocombo[st]))]
        
            elif useNoisy == 1 and useRms == 0:
                if self.noisy[ist] ==  True:
                    self.y[st] = [[1 if x <= (1. + err_noisy) \
                        and x >= (1. - err_noisy) else 0 for x in self.ratiocombo[st][i]] \
                        for i in range(len(self.ratiocombo[st]))]
                else:
                    self.y[st] = [[1 if x <= 1. + err \
                        and x >= 1. - err else 0 for x in self.ratiocombo[st][i]] \
                        for i in range(len(self.ratiocombo[st]))]
        
            elif useNoisy == 0 and useRms == 1:    
                self.y[st] = [[1 if x <= (1. + err + st.getNoise()) \
                    and x >= (1. - err - st.getNoise()) else 0 for x in self.ratiocombo[st][i]] \
                    for i in range(len(self.ratiocombo[st]))]
            
            else:
                self.y[st] = [[1 if x <= 1. + err and x >= 1. - err \
                    else 0 for x in self.ratiocombo[st][i]] \
                    for i in range(len(self.ratiocombo[st]))]
                        
        for st in self.translist:
            self.verdict_ratiocombo[st] = [0 if 0 in self.y[st][i] else 1 for i in range(len(self.y[st]))]
        
        
        
        self.model_ratiocombo = []
        
        for ii in range(len(self.star_grid)):
            self.model_ratiocombo.append([self.verdict_ratiocombo[tr][ii] for tr in self.translist])
        
        self.model_ratiocombo_verdict = sum(array(self.model_ratiocombo), axis = 1)    
        
        
        if plot:
            self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId() \
                for ii in range(len(self.star_grid))]
            tr = [self.verdict_ratiocombo[st] for st in self.translist]
            
            plt.clf()
            fig = plt.figure(3, figsize = (15, 10))
            ax = fig.add_subplot(111)
            ax.set_xticks(np.arange(len(tr))-0.5)
            ax.set_yticks(np.arange(len(self.modellist)))
            ax.xaxis.set_ticklabels([str(st.jup)+'-'+str(st.jlow)+' '+st.telescope for st in self.translist], rotation = 60)
            ax.yaxis.set_ticklabels([i for i in self.modellist])
            ax.imshow(self.model_ratiocombo, interpolation='nearest', origin='upper', cmap = 'Blues')
            plt.tight_layout()
            fig.suptitle('Combo - ratio of integrated intensities and peak intensities', size = 18)
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename_combo = os.path.join(path, 'resostats','int-%s_len_%s-%s'%(self.modellist[0],(len(self.modellist)),plot_id))
            fig.savefig(filename_combo+'.pdf')   
            print '*** Plot of stats can be found at:'
            print filename_combo+'.pdf'

    
    
    def calcChiSquared(self, P = 2, useTeleUncertainties = 1, useNoisy = 1, err = 0.2, err_noisy = 0.3):
        '''
        Calculate the (reduced) chi squared of the integrated main beam intensities.
        
        @keyword P: Degrees of freedom = len(data) - P
        @type P: int
        
        @keyword useTeleUncertainties: Use telescope uncertainties instead of a fixed error.
        @type useTeleUncertainties: bool
        
        @keyword useNoisy: Use a higher uncertainty for noisy lines
        @type useNoisy: bool
        
        @keyword err: Fixed uncertainty on lines (if not telescope uncertainties)
        @type err: float
        
        @keyword err_noisy: Uncertainty on noisy lines
        @type err_noisy: float 
        
        @return chiSquared: Chi squared values of each model
        @type chiSquared: list[]
        
        @return redChiSquared: Reduced chi squared values of each model
        @type redChiSquared: list[]
        
        @return errRedChiSquared: Error on reduced chi squared 
        @type errRedChiSquared: float
        
        @return redChiSquaredWithinThreeSigma: Models that have a reduced chi squared within
                                               three sigma of the best model (ie the model with
                                               the lowest reduced chi squared)
        @type redChiSquaredWithinThreeSigma: list[]
        '''
        self.modellist = [self.star_grid[ii]['GAS_LINES'][0].getModelId() \
            for ii in range(len(self.star_grid))]
       
        self.chiSquared = []
        self.redChiSquared = []
        data = array([self.dinttmb[st] for st in self.translist])
        dof = P - 1

        for mm in range(len(self.modellist)):
            model = array([self.minttmb[st][mm] for st in self.translist])
           
            if useTeleUncertainties == 1 and useNoisy == 1:
                noise = array([self.tele_uncertainties[t.telescope]*self.dinttmb[t] if not self.noisy[i] else err_noisy*self.dinttmb[t] \
                        for i,t in enumerate(self.translist)])
                cs = BasicStats.calcChiSquared(data,model,noise,ndf=dof)
                self.redChiSquared.append(cs)
                self.chiSquared.append(cs*(len(data) - dof - 1))
                
            if useTeleUncertainties == 0 and useNoisy == 1: 
                noise = array([err*self.dinttmb[t] if not self.noisy[i]*self.dinttmb[t] else err_noisy*self.dinttmb[t] \
                        for i,t in enumerate(self.translist)])
                cs = BasicStats.calcChiSquared(data,model,noise,ndf=dof)
                self.redChiSquared.append(BasicStats.calcChiSquared(data,model,noise,ndf=dof))
                self.chiSquared.append(cs*(len(data) - dof - 1))
                                           
            if useTeleUncertainties == 1 and useNoisy == 0:
                noise = [self.tele_uncertainties[t.telescope] for t in self.translist]*data
                cs = BasicStats.calcChiSquared(data,model,noise,ndf=dof)
                self.redChiSquared.append(BasicStats.calcChiSquared(data,model,noise,ndf=dof))
                self.chiSquared.append(cs*(len(data) - dof - 1))
                                           
            if useTeleUncertainties == 0 and useNoisy == 0:
                noise = err*data
                cs = BasicStats.calcChiSquared(data,model,noise,ndf=dof)
                self.redChiSquared.append(BasicStats.calcChiSquared(data,model,noise,ndf=dof))
                self.chiSquared.append(cs*(len(data) - dof - 1))


        self.errRedChiSquared = (2.0/len(self.translist))**0.5
       
        best = np.where(np.array(self.redChiSquared) == min(self.redChiSquared))[0][0]
       
        self.redChiSquaredWithinThreeSigma = [ii for ii in range(len(self.modellist)) \
                                            if self.redChiSquared[ii] > (min(self.redChiSquared) - 3*self.errRedChiSquared) \
                                            and self.redChiSquared[ii] < (min(self.redChiSquared) + 3*self.errRedChiSquared)]
      
           
       
       
        print '*******************************************************'
        print '****************** Chi Squared ************************'
        print '********  Chi squared  -  Reduced chi squared  ********'
        print '*******************************************************'
        for mm in range(len(self.modellist)):
            print str(self.modellist[mm]) + '\t' + str('%.3f' %self.chiSquared[mm]) + '\t' + str('%.3f' %self.redChiSquared[mm]) + '+-' + str('%.3f' %self.errRedChiSquared)
        print '\n'
        print 'Best model = ' + str(best) +', '+ str(self.modellist[best]) + ', with ' + str('%.3f' %self.redChiSquared[best])
        print '*******************************************************'

    
    

    
    
    def plotBestFit(self, plot_int = 1, plot_peak = 0, plot_combo = 0):
        '''
        Can only be perfomed after self.selectBestFitperLine().
        
        Visualization of self.occurences_bfm(mode), self,verdict_(mode)_soft,
        and self.occurences_bfmlll.        
        '''
        
        jup = [t.jup for t in self.translist]
        
        if plot_int:
            plt.close()
            fig1 = plt.figure(1, figsize = (20,11))
            ax1 = fig1.add_subplot(131)
            ax1.set_xticks(np.arange(len(jup)))
            ax1.xaxis.set_ticklabels(jup)
            ax1.set_yticks(np.arange(len(self.modellist)))
            ax1.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.fittedModels_bfmint[i])]) for i in range(len(self.modellist))])
            ax1.imshow(self.occurences_bfmint, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax1.set_xlabel('$J_{up}$')
            ax1.set_title('Int')

            ax2 = fig1.add_subplot(132)
            ax2.set_xticks(np.arange(len(jup)))
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels(jup)
            ax2.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.verdictpermodel_int[i])]) \
                for i in range(len(self.modellist))])
            ax2.imshow(self.verdict_int_soft, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_xlabel('$J_{up}$')
            ax2.set_title('$\int T_\mathrm{mb}$')

            ax3 = fig1.add_subplot(133)
            ax3.set_xticks(np.arange(len(jup)))
            ax3.set_yticks(np.arange(len(self.modellist)))
            ax3.xaxis.set_ticklabels(jup)
            ax3.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.fittedModels_bfmlll[i])]) for i in range(len(self.modellist))])
            ax3.imshow(self.occurences_bfmlll, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax3.set_xlabel('$J_{up}$')
            ax3.set_title('LLL')
            
            plt.tight_layout()
            
            fig1.suptitle('ResoStats mode = int',size = 18)
            
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename_int = os.path.join(path, 'resostats','int-%s_len_%i'%(self.modellist[0],(len(self.modellist))))
            fig1.savefig(filename_int+'.pdf')    
            
            print '** Plot of ResoStats Int can be found at:'
            print filename_int+'.pdf'
            print '***********************************'                            
        
        if plot_peak:
            plt.close()
            fig2 = plt.figure(2, figsize = (20,11))
            ax1 = fig2.add_subplot(131)
            ax1.set_xticks(np.arange(len(jup)))
            ax1.xaxis.set_ticklabels(jup)
            ax1.set_yticks(np.arange(len(self.modellist)))
            ax1.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.fittedModels_bfmpeak[i])]) for i in range(len(self.modellist))])
            ax1.imshow(self.occurences_bfmpeak, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax1.set_xlabel('$J_{up}$')
            ax1.set_title('Peak')

            ax2 = fig2.add_subplot(132)
            ax2.set_xticks(np.arange(len(jup)))
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels(jup)
            ax2.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.verdictpermodel_peak[i])]) \
                for i in range(len(self.modellist))])
            ax2.imshow(self.verdict_peak_hard, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_xlabel('$J_{up}$')
            ax2.set_title('$\int T_\mathrm{mb}$')

            ax3 = fig2.add_subplot(133)
            ax3.set_xticks(np.arange(len(jup)))
            ax3.set_yticks(np.arange(len(self.modellist)))
            ax3.xaxis.set_ticklabels(jup)
            ax3.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.fittedModels_bfmlll[i])]) for i in range(len(self.modellist))])
            ax3.imshow(self.occurences_bfmlll, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax3.set_xlabel('$J_{up}$')
            ax3.set_title('LLL')
            
            plt.tight_layout()
            
            fig2.suptitle('ResoStats mode = peak',size = 18)
            
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename_peak = os.path.join(path, 'resostats','peak-%s_len_%i'%(self.modellist[0],(len(self.modellist))))
            fig2.savefig(filename_peak+'.pdf')    
            
            print '** Plot of ResoStats Peak can be found at:'
            print filename_peak+'.pdf'
            print '***********************************'              
            
    
        if plot_combo:
            plt.close()
            fig3 = plt.figure(3, figsize = (20,11))
            ax1 = fig3.add_subplot(131)
            ax1.set_xticks(np.arange(len(jup)))
            ax1.xaxis.set_ticklabels(jup)
            ax1.set_yticks(np.arange(len(self.modellist)))
            ax1.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.fittedModels_bfmcombo[i])]) for i in range(len(self.modellist))])
            ax1.imshow(self.occurences_bfmcombo, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax1.set_xlabel('$J_{up}$')
            ax1.set_title('Combo')

            ax2 = fig3.add_subplot(132)
            ax2.set_xticks(np.arange(len(jup)))
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels(jup)
            ax2.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.verdictpermodel_combo[i])]) \
                for i in range(len(self.modellist))])
            ax2.imshow(self.verdict_combo_hard, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_xlabel('$J_{up}$')
            ax2.set_title('$\int T_\mathrm{mb}$')

            ax3 = fig3.add_subplot(133)
            ax3.set_xticks(np.arange(len(jup)))
            ax3.set_yticks(np.arange(len(self.modellist)))
            ax3.xaxis.set_ticklabels(jup)
            ax3.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.fittedModels_bfmlll[i])]) for i in range(len(self.modellist))])
            ax3.imshow(self.occurences_bfmlll, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax3.set_xlabel('$J_{up}$')
            ax3.set_title('LLL')
            
            plt.tight_layout()
            
            fig3.suptitle('ResoStats mode = combo',size = 18)
            
            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'resostats'))
            filename_combo = os.path.join(path, 'resostats','combo-%s_len_%i'%(self.modellist[0],(len(self.modellist))))
            fig3.savefig(filename_combo+'.pdf')    
            
            print '** Plot of ResoStats Combo can be found at:'
            print filename_combo+'.pdf'
            print '***********************************'              
    
        
    def printBestFitperLine(self, bfm_int = 1, bfm_peak = 0, bfm_combo = 0, bfm_lll = 1):
        '''
        Print output of selectBestFitperLine(), using fittedLines_bfmxxx.
        '''
        calcTrans = [str(self.translist[x]) for x in self.includedtrans]

        if bfm_int == 1:
            print '-------------------------------------------------------------------'
            print 'How often was a line modeled according to selectBestFitModels?'
            print 'Mode  = int'
            print '-------------------------------------------------------------------'
            print 'Number of models calculated = ' + str(len(self.modellist))
            print ''
            print 'Transition                  -                 Fitted by # models'
            for ii in range(len(calcTrans)):
                print calcTrans[ii] + "\t" + str(self.fittedLines_bfmint[ii])
        
        if bfm_peak == 1:
            print '-------------------------------------------------------------------'
            print 'How often was a line modeled according to selectBestFitModels?'
            print 'Mode  = peak'
            print '-------------------------------------------------------------------'
            print 'Number of models calculated = ' + str(len(self.modellist))
            print ''
            print 'Transition                  -                 Fitted by # models'
            for ii in range(len(calcTrans)):
                print calcTrans[ii] + "\t" + str(self.fittedLines_bfmpeak[ii])

        if bfm_combo == 1:
            print '-------------------------------------------------------------------'
            print 'How often was a line modeled according to selectBestFitModels?'
            print 'Mode  = combo'
            print '-------------------------------------------------------------------'
            print 'Number of models calculated = ' + str(len(self.modellist))
            print ''
            print 'Transition                  -                 Fitted by # models'
            for ii in range(len(calcTrans)):
                print calcTrans[ii] + "\t" + str(self.fittedLines_bfmcombo[ii])
            
        if bfm_lll == 1:
            print '----------------------------------------------------------------------'
            print 'How often was a line modeled according to selectBestFitModelsLLL?'
            print '----------------------------------------------------------------------'
            print 'Number of models calculated = ' + str(len(self.modellist))
            print ''
            print 'Transition                  -                 Fitted by # models'
            for ii in range(len(calcTrans)):
                print calcTrans[ii] + "\t" + str(self.fittedLines_bfmlll[ii])            
            
            
    def printBestFitperModel(self, bfm_int = 1, bfm_peak = 0, bfm_combo = 0, bfm_lll = 1):
        '''
        Print output of selectBestFitperLine(), using fittedModels_bfmxxx.
        '''
        if bfm_int == 1:
            print '---------------------------------------------------------------------'
            print 'Which model modelled the most lines according to selectBestFitModels?'
            print 'Mode = int'
            print '---------------------------------------------------------------------'
            print 'Number of lines included = ' + str(len(self.includedtrans))
            print ''
            print 'Model          -         # lines fitted '
            for ii in range(len(self.modellist)):
                print self.modellist[ii] + "\t" + str(self.fittedModels_bfmint[ii])
            print ''

        if bfm_peak == 1:
            print '---------------------------------------------------------------------'
            print 'Which model modelled the most lines according to selectBestFitModels?'
            print 'Mode = peak'
            print '---------------------------------------------------------------------'
            print 'Number of lines included = ' + str(len(self.includedtrans))
            print ''
            print 'Model          -         # lines fitted '
            for ii in range(len(self.modellist)):
                print self.modellist[ii] + "\t" + str(self.fittedModels_bfmpeak[ii])
            print ''
            
        if bfm_combo == 1:
            print '---------------------------------------------------------------------'
            print 'Which model modelled the most lines according to selectBestFitModels?'
            print 'Mode = combo'
            print '---------------------------------------------------------------------'
            print 'Number of lines included = ' + str(len(self.includedtrans))
            print ''
            print 'Model          -         # lines fitted '
            for ii in range(len(self.modellist)):
                print self.modellist[ii] + "\t" + str(self.fittedModels_bfmcombo[ii])
            print ''
            
        if bfm_lll == 1:
            print '-----------------------------------------------------------------------'
            print 'Which model modelled the most lines according to selectBestFitModelsLLL?'
            print '-----------------------------------------------------------------------'
            print 'Number of lines included = ' + str(len(self.includedtrans))
            print ''
            print 'Model          -         # lines fitted '
            for ii in range(len(self.modellist)):
                print self.modellist[ii] + "\t" + str(self.fittedModels_bfmlll[ii])
            print ''
        