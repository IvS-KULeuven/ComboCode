# -*- coding: utf-8 -*-

"""
Performing simple peak-to-peak statistics.

Author: R. Lombaert

"""

import os
import scipy
from scipy import argmin,array,sqrt
import operator

from cc.tools.io import DataIO
from cc.modeling.objects import Transition
from cc.statistics.Statistics import Statistics
from cc.statistics.BasicStats import calcChiSquared
from cc.plotting import Plotting2



class PeakStats(Statistics):
    
    """
    Environment with several tools to perform statistics on peak-to-peak ratios
    for unresolved lines.
    
    """
        
    def __init__(self,star_name,instrument,path_code='codeSep2010',\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):        
        
        """ 
        Initializing an instance of PeakStats.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        @param instrument: The instrument with which the data were taken 
                           (PACS or SPIRE)
        @type instrument: string
        
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2010')
        @type path_code: string
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        
        """
        
        super(PeakStats,self).__init__(star_name=star_name,\
                                       path_combocode=path_combocode,\
                                       code='GASTRoNOoM',path_code=path_code)
        
        #-- Set the type of instrument used in this instance
        self.instrument = instrument.upper()
        if self.instrument == 'FREQ_RESO':
            raise IOError('Peak2Peak ratios cannot yet be calculated for FREQ_RESO.')
        
        #-- Remember sample transitions per band and their central wavelengths. 
        #   key: filename, value: list[Transition()]
        self.sample_trans = dict()
        self.central_mwav = dict()
        
        #-- Remember the peak and integrated intensity ratios, and chi^2. 
        ##   key: filename
        #   value: dict([instrument based id, list[float] or float])
        self.peak_ratios = dict()     
        self.int_ratios = dict()
        self.int_ratios_err = dict()
        self.chi2 = dict()
        
        #-- Set the upper and lower data noise limits
        upper_limits = dict()
        upper_limits['PACS'] = 1.3
        lower_limits = dict()
        lower_limits['PACS'] = 0.7
        self.upper_limit = upper_limits[self.instrument]
        self.lower_limit = lower_limits[self.instrument]
        
        
    
    def setInstrument(self,*args,**kwargs):
        
        '''
        Set an instrument, see Statistics.py.
        
        '''
        
        super(PeakStats,self).setInstrument(instrument=self.instrument,*args,\
                                            **kwargs)
        
        
        
    def setModels(self,*args,**kwargs):
        
        '''
        Set the model grid, see Statistics.py.
        
        '''
        
        super(PeakStats,self).setModels(instrument=self.instrument,*args,\
                                        **kwargs)
        
        

    def setRatios(self):
        
        ''' 
        Find the peak to peak ratios of data versus model.
        
        The result are saved in the self.peak_ratios dictionary, see 
        __init__.__doc__()
        
        '''
        
        print '***********************************'
        print '** Calculating integrated/peak intensity ratios for %s.'\
              %self.instrument
        
        #-- Make a sample selection of Transition()s. 
        sample_trans = Transition.extractTransFromStars(self.star_grid,pacs=1)
        inst = self.instruments[self.instrument]
        
        #-- Set the tolerance: a factor that is multiplied with half the 
        #   wavelength resolution of data when looking for the peak value 
        #   around a central wavelength in the data. As a default, this is 
        #   equal to the PACS_OVERSAMPLING. No point in changing this.
        self.tolerance = inst.oversampling
        
        for ifn,(fn,dwav) in enumerate(zip(inst.data_filenames,\
                                           inst.data_wave_list)):
            fn = os.path.split(fn)[1]
            #-- Create a list of sample transitions
            self.sample_trans[fn] = [trans
                                     for trans in sample_trans
                                     if trans.wavelength*10**4 >= dwav[0]\
                                    and trans.wavelength*10**4 <= dwav[-2]]
            #-- Get the central wavelength of the lines, corrected for 
            #   Doppler shift due to vlsr of the central source. In micron.
            self.central_mwav[fn] = [t.wavelength*10**4*1./(1-inst.vlsr/t.c)
                                     for t in self.sample_trans[fn]]
            
            self.__setPeakRatios(ifn,fn)
            if self.instrument == 'PACS' and inst.linefit <> None:
                self.__setIntRatios(ifn,fn)
        print '***********************************'
                
                
                
    def __setIntRatios(self,ifn,fn):
        
        '''
        Calculate ratios of integrated intensities, if requested. 
        
        Assumes there is a linefit filename available in the Instrument() 
        object. Only done for those line present in this file, based on 
        Doppler shifted wavelength.
        
        @param ifn: Index of the data band in self.instruments[self.instrument]
                    lists.
        @type ifn: int
        @param fn: The filename of the data set. Needed for book keeping.
        @type fn: string
        
        '''
        
        #-- Get some data properties, and extract data wavelength and flux
        inst = self.instruments[self.instrument]
        self.int_ratios[fn] = dict()
        self.int_ratios_err[fn] = dict() 
        
        #-- Comparing integrated intensities between PACS and models.
        #   Comparisons only made for the same band! 
        #   1) Select the info belonging to this particular band    
        ordername = inst.data_ordernames[ifn]
        lf = inst.linefit[inst.linefit['band'] == ordername]
        
        #   2) Check if the wav of an mtrans matches a wav in the fitted
        #      intensities list, within the fitted_fwhm/2 of the line with 
        #      respect to the fitted central wavelength.
        #      Note that there cannot be 2 fitted lines with central wavelength
        #      closer than fwhm/2. Those lines would be inseparable!
        imatches = [argmin(abs(lf.wave_fit-mwav))
                    for mwav in self.central_mwav[fn]]
        matches = [(mwav <= lf.wave_fit[ii] + lf.fwhm_fit[ii]/2. \
                        and mwav >= lf.wave_fit[ii] - lf.fwhm_fit[ii]/2.) \
                    and (lf.wave_fit[ii],ii) or (None,None)
                   for mwav,ii in zip(self.central_mwav[fn],imatches)]
        
        #   3) If match found, check if multiple mtrans fall within   
        #      fitted_FWHM/2 from the fitted central wavelength of the
        #      line. These are blended IN MODEL and/or IN DATA.
        matches_wv = array([mm[0] for mm in matches])
        wf_blends = [list(array(self.sample_trans[fn])[matches_wv==wv])
                     for wv in lf.wave_fit]
        blended = [ii <> None and len(wf_blends[ii]) > 1. \
                        and wf_blends[ii].index(st) != 0
                   for st,(match,ii) in zip(self.sample_trans[fn],matches)]
        
        for star in self.star_grid:
            #--  From here on, we start extracting the model specific int ints.
            this_id = star['LAST_%s_MODEL'%self.instrument]
            mtrans = array([star.getTransition(t) 
                            for t in self.sample_trans[fn]])
            these_ratios = []
            these_errs = []
            for mt,blend,(match,ii) in zip(mtrans,blended,matches):
                #   4) No trans == sample_trans found for this model, or no
                #      match found in linefit for this band, or line is blended
                #      with other line that is already added: ratio is None.
                if mt is None or match is None or blend:
                    these_ratios.append(None)
                    these_errs.append(None)
               
                #   5) Match found with a wave_fit value once. Get int ratio 
                #      m/d + check for line blend IN DATA:
                #      Check the ratio fitted FWHM/PACS FWHM. If larger by 20% 
                #      or more, put the int ratio negative. 
                elif len(wf_blends[ii]) == 1:
                    this_ratio = mt.getIntIntIntSphinx()/lf.line_flux[ii]
                    factor = lf.fwhm_rel[ii] >= 1.2 and -1 or 1
                    these_ratios.append(factor*this_ratio)
                    err = sqrt(lf.line_flux_rel[ii]/100**2+0.2**2)
                    these_errs.append(err*this_ratio)
                
                #   4) If multiple matches, add up the integrated intensities 
                #      of the mtrans and make m/d based on that sum. Make 
                #      negative to indicate a line blend. 
                #      The *OTHER* transitions found this way are not compared 
                #      with any data and get None. 
                else: 
                    blended_mt = [star.getTransition(st) 
                                  for st in wf_blends[ii]
                                  if star.getTransition(st) <> None]
                    total_mflux = sum([bmt.getIntIntIntSphinx() 
                                       for bmt in blended_mt])
                    this_ratio = total_mflux/lf.line_flux[ii]
                    these_ratios.append(-1.*this_ratio)
                    err = sqrt(lf.line_flux_rel[ii]/100**2+0.2**2)
                    these_errs.append(err*this_ratio)
            self.int_ratios[fn][this_id] = these_ratios
            self.int_ratios_err[fn][this_id] = these_errs


    def __setPeakRatios(self,ifn,fn):
        
        '''
        Calculate Peak ratios for all models included in the star grid. 
        
        Done per filename. 
        
        @param ifn: Index of the data band in self.instruments[self.instrument]
                    lists.
        @type ifn: int
        @param fn: The filename of the data set. Needed for book keeping.
        @type fn: string
        
        '''
        
        #-- Get some data properties, and extract data wavelength and flux
        inst = self.instruments[self.instrument]
        dwav = inst.data_wave_list[ifn]
        dflux = inst.data_flux_list[ifn]
        
        #-- Get some data statistics
        d_mean = self.data_stats[self.instrument][fn]['mean']
        d_std = self.data_stats[self.instrument][fn]['std']
        #-- this sigma is used for determining whether a peak flux is
        #   significant compared to the std value of the band.
        d_sigma = self.data_stats[self.instrument][fn]['sigma']
        
        self.peak_ratios[fn] = dict()
        self.chi2[fn] = dict() 
        
        
        for star in self.star_grid:
            #-- Read the convolved sphinx model
            this_id = star['LAST_%s_MODEL'%self.instrument]
            sphinx_file = os.path.join(os.path.expanduser('~'),\
                                       'GASTRoNOoM',self.path_code,\
                                       'stars',self.star_name,\
                                       '%s_results'%self.instrument,this_id,\
                                       '%s_%s'%('sphinx',fn))
            mwav, mflux = DataIO.readCols(sphinx_file)
            if list(mflux[mflux < 0]) != []: 
                print 'There are negative sphinx flux values! They will '+\
                      'not be taken into account.'
            
            #-- Calculate the chi-squared value for this model
            self.chi2[fn][this_id] = calcChiSquared(dflux[mflux>0],\
                                                    mflux[mflux>0],d_std)
            
            #-- Calculate the peak-to-peak ratios. 
            #   1) Central wavelengths of mtrans are set in previous method
            #   2) Get the central model flux, which should coincide exactly 
            #      with the Doppler shifted rest wavelength of the line. If the 
            #      model flux is negative, the value is not used.
            central_mflux = [mflux[argmin(abs(mwav-wav))] > 0 \
                                and mflux[argmin(abs(mwav-wav))] \
                                or None
                             for wav in self.central_mwav[fn]]
            #   3) Get the central data flux, at the Doppler shifted central 
            #      wavelength expected from the model. The maximum flux is 
            #      taken in the wavelength bin tolerance*wav_resolution/2.
            #      allowing for small wave shifts due to instrumental effects.
            central_dflux = [max(dflux[abs(dwav-wav)<= \
                             self.tolerance/2.*( dwav[argmin(abs(dwav-wav))+1]\
                                                -dwav[argmin(abs(dwav-wav))])])
                             for wav in self.central_mwav[fn]]
            #   4) Check if the data flux point is actually significant 
            #      compared to the noise in the spectrum. Compare with dstd, 
            #      given d_sigma from path_combocode/Data/Data.dat .
            #      Insignificant values are multiplied by -1, to indicate they
            #      are upper limits in the data at that wavelength.
            central_dflux = [d >= d_mean+(d_std*d_sigma) \
                                and d \
                                or -1*abs(d_mean+(d_std*d_sigma)) 
                             for d in central_dflux]
            #   5) Calculate the ratios, only if the model flux is not None 
            #      (was a negative model flux value: We don't want that)
            #      Negative ratios are possible, in case of ratio lower limits 
            self.peak_ratios[fn][this_id] = [m <> None and m/d or None
                                             for m,d in zip(central_mflux,\
                                                            central_dflux)]
                                                            
    
                                                                                
    def getRatios(self,this_id,sel_type='peak_ratios',\
                  data_type='peak_ratios',return_negative=0,filename=None):
        
        '''
        Return all ratios for all filenames, including only the values that are
        not None.
        
        @param this_id: The requested instrument id
        @type this_id: string 
        
        @keyword sel_type: The type of data used to select values, one of 
                           'peak_ratios' or 'int_ratios'
                              
                           (default: 'peak_ratios')
        @type sel_type: string
        @keyword data_type: the type of data returned, one of 'peak_ratios' or
                            'central_wav'
                            
                            (default: 'peak_ratios')
        @type data_type: string
        @keyword return_negative: only return the negative peak ratios or 
                                  equivalent (lower limits). 
                                  
                                  (default: 0)
        @type return_negative: bool
        @keyword filename: the filename for which you want to return the list, 
                           if None, all filenames are used and the lists are 
                           merged into one
                           
                           (default: None)
        @type filename: string
        
        @return: The ratios requested
        @rtype: array
        
        '''

        inst = self.instruments[self.instrument]
        filenames = filename is None and inst.data_filenames or [filename]
        filenames = [os.path.split(fn)[1] for fn in filenames]
        if data_type == 'central_wav':
            return array([v2     
                          for fn in filenames
                          for v,v2 in zip(getattr(self,sel_type)[fn][this_id],\
                                          self.central_mwav[fn])
                          if v <> None \
                            and ((return_negative and v < 0) \
                                 or (not return_negative and v > 0))])
        else:
            return array([v2  
                          for fn in filenames
                          for v,v2 in zip(getattr(self,sel_type)[fn][this_id],\
                                          getattr(self,data_type)[fn][this_id])
                          if v <> None \
                            and ((return_negative and v < 0) \
                                 or (not return_negative and v > 0))])



    def plotRatioWav(self,inputfilename,no_peak=False):
        
        '''
        Plot peak ratios as a function of their central wavelength.
        
        @param inputfilename: the input filename for the grid, which will be 
                              attached to the final plot filename
        @type inputfilename: string
        
        @keyword no_peak: Hide the peak ratios. 
                          
                          (default: False)
        @type no_peak: bool
        
        '''

        this_grid = self.sortStarGrid()
        plot_filenames = []
        inst = self.instruments[self.instrument]
        for star in this_grid:
            this_id = star['LAST_%s_MODEL'%self.instrument]
            chi2 = sum([vals[this_id] for vals in self.chi2.values()])
            lp = [] 
            waves = []
            ratios = []
            ratios_err = []
            
            #-- the ratios included in the statistics
            if not no_peak:
                this_wav_inc = self.getRatios(data_type='central_wav',\
                                            this_id=this_id)
                this_ratio_inc = self.getRatios(this_id=this_id)
                if list(this_wav_inc):    
                    waves.append(this_wav_inc)
                    ratios.append(this_ratio_inc)
                    ratios_err.append(None)
                    lp.append('ob')
                
                #-- ratios replaced by lower limits if data point is in the noise
                #   ie data point can only be smaller than or equal to used value
                this_wav_lower = self.getRatios(data_type='central_wav',\
                                                this_id=this_id,return_negative=1)
                this_ratio_lower = self.getRatios(return_negative=1,\
                                                this_id=this_id)
                this_ratio_lower = [abs(r) for r in this_ratio_lower]
                if list(this_wav_lower):
                    waves.append(this_wav_lower)
                    ratios.append(this_ratio_lower)
                    ratios_err.append(None)
                    lp.append('dg')
            
            if self.instrument == 'PACS' and inst.linefit <> None:
                #-- If integrated intensities are available for the instrument, get
                #   the integrated intensity ratios
                this_wav_int = self.getRatios(sel_type='int_ratios',\
                                              data_type='central_wav',\
                                              this_id=this_id)
                this_ratio_int = self.getRatios(sel_type='int_ratios',\
                                                data_type='int_ratios',\
                                                this_id=this_id)
                this_ratio_int_err = self.getRatios(sel_type='int_ratios',\
                                                    data_type='int_ratios_err',\
                                                    this_id=this_id)
                if list(this_wav_int):
                    waves.append(this_wav_int)
                    ratios.append(this_ratio_int)
                    ratios_err.append(this_ratio_int_err)
                    lp.append('or')
                
                #-- Get the ratios that are lower limits due to line blends.
                #   Line blends detected due to fitted FHWM/PACS FHWM > 120%
                #   ie model int can only be larger than or equal to used value
                this_wav_lowerint = self.getRatios(sel_type='int_ratios',\
                                                   data_type='central_wav',\
                                                   this_id=this_id,\
                                                   return_negative=1)
                this_ratio_lowerint = self.getRatios(sel_type='int_ratios',\
                                                     data_type='int_ratios',\
                                                     this_id=this_id,\
                                                     return_negative=1)
                this_ratio_lowerint_err = self.getRatios(sel_type='int_ratios',\
                                                    data_type='int_ratios_err',\
                                                    this_id=this_id)
                this_ratio_lowerint = [abs(r) for r in this_ratio_lowerint]
                if list(this_wav_lowerint):
                    waves.append(this_wav_lowerint)
                    ratios.append(this_ratio_lowerint)
                    ratios_err.append(this_ratio_lowerint_err)
                    lp.append('dm')
                    
            #- prepping input for the plot command
            xmin = min([min(x) for x in waves])
            xmax = max([max(x) for x in waves])
            waves.extend([[0.5*xmin,1.5*xmax]]*3)
            ratios.extend([[1,1],\
                           [self.lower_limit,\
                            self.lower_limit],\
                           [self.upper_limit,\
                            self.upper_limit]])
            ratios_err.extend([None,None,None])
            lp.extend(['-k','--k','--k'])
            plot_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                         self.path_code,'stars',\
                                         self.star_name,\
                                         '%s_results'%self.instrument,this_id,\
                                         'ratio_wav_%s'%this_id)
            labels = [('Mdot = %.2e Msolar/yr'%star['MDOT_GAS'],0.05,0.05),\
                      ('Teff = %.1f K'%star['T_STAR'],0.05,0.1),\
                      ('$\psi$ = %0.2e'%star['DUST_TO_GAS_CHANGE_ML_SP'],0.05,\
                       0.15),\
                      ('A$_{H_2O}$/A$_{H_2}$ = %0.2e'%star['F_H2O'],0.05,0.2),\
                      ('R$_(o,H_2O)$ = %i'%int(star['R_OUTER_H2O']),0.05,0.25)]
            if star.getMolecule('1H1H16O') \
                  and star.getMolecule('1H1H16O').set_keyword_change_abundance:
                labels.append(('$H_2O$ profile = %s'\
                               %os.path.split(star.getMolecule('1H1H16O')\
                                                   .change_fraction_filename\
                                                   .replace('_','\_'))[1],\
                               0.05,0.30))
            plot_title = '%s: $\chi^2$ %.4f'%(this_id.replace('_','\_'),chi2)
            plot_filenames.append(Plotting2.plotCols(\
                    filename=plot_filename,x=waves,y=ratios,yerr=ratios_err,\
                    yaxis=r'$F_{\nu,p,m}/F_{\nu,p,d}$',\
                    plot_title=plot_title,labels=labels,extension='pdf',\
                    xlogscale=0,ylogscale=1,line_types=lp,xmin=xmin*0.9,\
                    xmax=xmax*1.03,figsize=(10.*scipy.sqrt(2.), 10.),\
                    linewidth=2,fontsize_title=20,fontsize_label=16))
        inputf_short = os.path.splitext(os.path.split(inputfilename)[1])[0]
        new_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                    self.path_code,'stars',self.star_name,\
                                    '%s_results'%self.instrument,\
                                    'ratio_wav_%s.pdf'%inputf_short)
        DataIO.joinPdf(old=plot_filenames,new=new_filename)
        print '** Stat plots can be found at:'
        print new_filename
        print '***********************************'



    def sortStarGrid(self):
        
        '''
        Return a sorted list of the star grid in this instance according to 
        the chi^2 value determined from convolved model vs data.
        
        @return: A sorted list of models according to chi^2
        @rtype: list[Star]
        
        '''
        
        chi2s = dict()
        for s in self.star_grid:
            this_id = s['LAST_%s_MODEL'%self.instrument]
            chi2s[this_id] = sum([vals[this_id] 
                                  for vals in self.chi2.values()])
        
        return sorted(self.star_grid,\
                      key=lambda x: chi2s[x['LAST_%s_MODEL'%self.instrument]])