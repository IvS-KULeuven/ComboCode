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
        #-- The default uncertainties were taken from the manuals of the 
        #   respective telescopes, and give a safe estimate based on the given
        #   values. 
        #   APEX - http://www.apex-telescope.org/~mdumke/publications/2010/spie2010/node4.html
        #   Herschel - See respective manuals (times 2, for other sources of err)
        #   JCMT - See Kemper et al. 2003 p611
        #   Other instrument flux calibration uncertainties are arbitrary! 
        #   I suggest to change these values based on what you want to use 
        #   yourself, depending on the transition under consideration.
        self.tele_uncertainties = dict([('APEX',.15),\
                                        ('CSO',.20),\
                                        ('FCRAO',.20),\
                                        ('HIFI',.20),\
                                        ('IRAM',.20),\
                                        ('JCMT',.30),\
                                        ('MOPRA',.20),\
                                        ('NRAO',.20),\
                                        ('OSO',.20),\
                                        ('SEST',.20)])
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

        
    def setIntensities(self,use_bestvlsr=1):
        
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
            self.dinttmb[st] = st.getIntTmbData()
            self.dpeaktmb[st] = st.getPeakTmbData() 
            
            if self.dpeaktmb[st] <= 3*noise:
                self.noisy[ist] = True
            else: 
                self.noisy[ist] = False
            
            #-- Collect the loglikelihoods for all models
            self.loglikelihood[st] = array([mt.getLoglikelihood(use_bestvlsr) 
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
        
        #self.modellist = list(set([t.getModelId() for x in range(len(self.translist)) \
              #for t in self.trans_models[self.translist[x]]]))
        #stars = array(self.modellist)
        
        
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
                    if lll_thresh <> None and not self.noisy[ist] \
                            and lll < lll_thresh:
                        bfbools[i] = False
        self.bfmlll = stars[bfbools]      
        self.bfmlll = list(self.bfmlll)
        
        return self.bfmlll    

    
    def findLLLMaximizer(self):
        
        self.difflll = [[] for _ in range(len(self.modellist))]
        
        #self.dlll = dict()
        
        for ist,st in enumerate(self.translist):
            lll_thresh = self.lll_threshold[ist]
            for i,lll in enumerate(self.loglikelihood[st]):
                if lll_thresh <> None and not self.noisy[ist]:
                    self.difflll[i].append(lll-lll_thresh) 
                    #self.dlll[st].append(lll-lll_thresh)
                else:
                    self.difflll[i].append(0)
                    #self.dlll[st].append(None)
        
        
    
    

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
        self.all_bfmint = []
        self.all_bfmpeak = []
        self.all_bfmcombo = []
        self.all_bfmlll = []
        for ii in orig_included:
            self.excludeTrans([i for i in range(T)])
            self.includeTrans(ii)
            self.all_bfmint.append(self.selectBestFitModels('int', use_lll, output = 0))
            self.all_bfmpeak.append(self.selectBestFitModels('peak', use_lll, output = 0))
            self.all_bfmcombo.append(self.selectBestFitModels('combo', use_lll, output = 0))
            self.all_bfmlll.append(self.selectBestFitModelsLLL(output = 0))
        
        #-- Make the main list
        self.occurences_bfmint = [[self.all_bfmint[x].count(self.modellist[i]) \
            for x in range(len(self.all_bfmint))] \
            for i in range(len(self.modellist))]
        
        self.occurences_bfmpeak = [[self.all_bfmpeak[x].count(self.modellist[i]) \
            for x in range(len(self.all_bfmpeak))] \
            for i in range(len(self.modellist))]    
        
        self.occurences_bfmcombo = [[self.all_bfmcombo[x].count(self.modellist[i]) \
            for x in range(len(self.all_bfmcombo))] \
            for i in range(len(self.modellist))] 
        
        self.occurences_bfmlll = [[self.all_bfmlll[x].count(self.modellist[i]) \
            for x in range(len(self.all_bfmlll))] \
            for i in range(len(self.modellist))]

        #-- Set included_trans back to original state
        self.includedtrans = orig_included
        
        #-- Initialise arrays
        self.occurences_bfmint = array(self.occurences_bfmint)
        self.occurences_bfmpeak = array(self.occurences_bfmpeak)
        self.occurences_bfmcombo = array(self.occurences_bfmcombo)
        self.occurences_bfmlll = array(self.occurences_bfmlll)
        
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
    
        self.maxint = np.where(self.fittedModels_bfmint == max(self.fittedModels_bfmint))[0]
        self.maxpeak = np.where(self.fittedModels_bfmpeak == max(self.fittedModels_bfmpeak))[0]
        self.maxcombo = np.where(self.fittedModels_bfmcombo == max(self.fittedModels_bfmcombo))[0]
        self.maxlll = np.where(self.fittedModels_bfmlll == max(self.fittedModels_bfmlll))[0]
    
    
    def findBestModelIntersect(self):
        '''
        Find which models fit the largest number of lines using mode = 'int' and
        mode = 'peak'.
        
        @return: Indices of the models in self.modellist
        @rtype: list[int]
        '''
        self.bestModel_index = np.intersect1d(self.maxint,self.maxpeak)
        
        return self.bestModel_index
    
    
    #def findBestModel(self):
        
        #self.maxint = np.where(self.fittedModels_bfmint == max(self.fittedModels_bfmint))[0]
        #self.maxpeak = np.where(self.fittedModels_bfmpeak == max(self.fittedModels_bfmpeak))[0]
        #self.maxcombo = np.where(self.fittedModels_bfmcombo == max(self.fittedModels_bfmcombo))[0]
        #self.maxlll = np.where(self.fittedModels_bfmlll == max(self.fittedModels_bfmlll))[0]
        
        
        
        
        
    def listCorrespondingAbundance(self, to_print):
        '''
        List the input abundances for model(s) of self.modellist.
        
        @param to_print: Indices of the models in self.modellist
        @type to_print: list[int]
        '''
        for ii in to_print:
            print self.modellist[ii], self.star_grid[ii]['MOLECULE'][0][6:10], \
                self.star_grid[ii]['VELOCITY_BETA']
        
        
   
    def checkTmb(self):
        
        S = len(self.star_grid)
        T = len(self.translist)
        
        ###- Integrated main beam temperature
        verdict_int_hard = [[] for _ in range(S)]
        verdict_int_soft = [[] for _ in range(S)]

        ratios = []
        for st in self.translist:
            ratios.append(self.ratioint[st])
        ratios_int = np.array(ratios).T
   
        for ii in range(S):
            for jj in range(T):
                if ratios_int[ii][jj] <= 1.20 and ratios_int[ii][jj] >= 0.80:
                    verdict_int_hard[ii].append(1)
                else:
                    verdict_int_hard[ii].append(0)
         
        for ii in range(S):
            for jj in range(T):
                if ratios_int[ii][jj] < 1.15 and ratios_int[ii][jj] > 0.85:
                    verdict_int_soft[ii].append(1)
                elif ratios_int[ii][jj] >= 1.15 and ratios_int[ii][jj] <= 1.20\
                    or ratios_int[ii][jj] >= 0.80 and ratios_int[ii][jj] <= 0.85:
                    verdict_int_soft[ii].append(0.5)
                else:
                    verdict_int_soft[ii].append(0)
                    
        self.ratios_int = ratios_int
        self.verdict_int_hard = verdict_int_hard
        self.verdict_int_soft = verdict_int_soft
        self.verdictpermodel_int = sum(array(self.verdict_int_hard), axis = 1)
    
        ###- Peak main beam temperature
        verdict_peak_hard = [[] for _ in range(S)]
        verdict_peak_soft = [[] for _ in range(S)]
        
        ratios = []
        for st in self.translist:
            ratios.append(self.ratiopeak[st])
        ratios_peak = np.array(ratios).T
   
        for ii in range(S):
            for jj in range(T):
                if ratios_peak[ii][jj] <= 1.20 and ratios_peak[ii][jj] >= 0.80:
                    verdict_peak_hard[ii].append(1)
                else:
                    verdict_peak_hard[ii].append(0)
                            
        for ii in range(S):
            for jj in range(T):
                if ratios_peak[ii][jj] < 1.15 and ratios_peak[ii][jj] > 0.85:
                    verdict_peak_soft[ii].append(1)
                elif ratios_peak[ii][jj] >= 1.15 and ratios_peak[ii][jj] <= 1.20\
                    or ratios_peak[ii][jj] >= 0.80 and ratios_peak[ii][jj] <= 0.85:
                    verdict_peak_soft[ii].append(0.5)
                else:
                    verdict_peak_soft[ii].append(0)
                    
        self.ratios_peak = ratios_peak
        self.verdict_peak_hard = verdict_peak_hard
        self.verdict_peak_soft = verdict_peak_soft
        self.verdictpermodel_peak = sum(array(self.verdict_peak_hard), axis = 1)
        
        
        ###- Combination of int and peak
        
        verdict_combo_hard = []
        
        for ii in range(S):
            verdict_combo_hard.append([self.verdict_int_hard[ii][x] \
                if self.verdict_int_hard[ii][x] == self.verdict_peak_hard[ii][x] else 0 \
                for x in range(len(self.translist))])
        
        self.verdict_combo_hard = verdict_combo_hard
        self.verdictpermodel_combo = sum(array(self.verdict_combo_hard), axis = 1)
        
        
    
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
    
        
    
    
    def plotCheckTmb(self, plot_int = 1, plot_peak = 0, plot_combo = 0):
        
        '''
        Can only be perfomed after self.checkTmb().
        
        Visualization of self.verdict_(mode)_soft. 
        Dark blue: model \pm 15%, light blue: model \pm 15-20% (so more at the edge).
        '''
        
        
        jup = [t.jup for t in self.translist]
        
        plt.close()
        
        if plot_int:
            fig = plt.figure(2, figsize = (20,10))

            ax2 = fig.add_subplot(111)
            ax2.set_xticks(np.arange(len(jup)))
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels(jup)
            ax2.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.verdictpermodel_int[i])]) \
                for i in range(len(self.modellist))])
            ax2.imshow(self.verdict_int_soft, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_xlabel('$J_{up}$')
            ax2.set_title('Dark blue: model = data $\pm$ 15% \n \
                Light blue: model = data $\pm$ 15-20%')
            
            fig.suptitle('Ratio of integrated main beam temperature: model/data', size = 18)


            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'StatsTmb'))
            filename = os.path.join(path, 'StatsTmb','%s_len_%i'%(self.modellist[0],(len(self.modellist))))
        
            
            print '***********************************'   
            print '** Plot of Stats Tmb can be found at:'
            print filename+'.pdf'
            print '***********************************'                            
            
            fig.savefig(filename+'.pdf')    

        
        if plot_peak:
            fig = plt.figure(2, figsize = (20,10))

            ax2 = fig.add_subplot(111)
            ax2.set_xticks(np.arange(len(jup)))
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels(jup)
            ax2.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.verdictpermodel_peak[i])]) \
                for i in range(len(self.modellist))])
            ax2.imshow(self.verdict_peak_soft, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_xlabel('$J_{up}$')
            ax2.set_title('Dark blue: model = data $\pm$ 15% \n \
                Light blue: model = data $\pm$ 15-20%')
            
            fig.suptitle('Ratio of integrated main beam temperature: model/data', size = 18)


            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'StatsTmb'))
            filename = os.path.join(path, 'StatsTmb','%s_len_%i'%(self.modellist[0],(len(self.modellist))))
        
            
            print '***********************************'   
            print '** Plot of Stats Tmb can be found at:'
            print filename+'.pdf'
            print '***********************************'                            
            
            fig.savefig(filename+'.pdf')            

        
        if plot_combo:
            fig = plt.figure(2, figsize = (20,10))

            ax2 = fig.add_subplot(111)
            ax2.set_xticks(np.arange(len(jup)))
            ax2.set_yticks(np.arange(len(self.modellist)))
            ax2.xaxis.set_ticklabels(jup)
            ax2.yaxis.set_ticklabels([" -- ".join([self.modellist[i][6:], str(self.verdictpermodel_combo[i])]) \
                for i in range(len(self.modellist))])
            ax2.imshow(self.verdict_combo_soft, interpolation='nearest', origin='upper', cmap = 'Blues')
            ax2.set_xlabel('$J_{up}$')
            ax2.set_title('Dark blue: model = data $\pm$ 15% \n \
                Light blue: model = data $\pm$ 15-20%')
            
            fig.suptitle('Ratio of integrated main beam temperature: model/data', size = 18)


            path = os.path.join(getattr(cc.path,self.code.lower()), self.path_code,'stars', self.star_name)
            DataIO.testFolderExistence(os.path.join(path,'StatsTmb'))
            filename = os.path.join(path, 'StatsTmb','%s_len_%i'%(self.modellist[0],(len(self.modellist))))
        
            
            print '***********************************'   
            print '** Plot of Stats Tmb can be found at:'
            print filename+'.pdf'
            print '***********************************'                            
            
            fig.savefig(filename+'.pdf')          
        
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
        