# -*- coding: utf-8 -*-

"""
Performing simple peak-to-peak statistics.

Author: R. Lombaert

"""

import os
import scipy
from scipy import argmin,array
import operator

from cc.tools.io import DataIO
from cc.statistics.Statistics import Statistics
from cc.plotting import Plotting2



class PeakStats(Statistics):
    
    """
    Environment with several tools to perform statistics on peak-to-peak ratios
    for unresolved lines.
    
    """
        
    def __init__(self,star_name,path_code='codeSep2010',\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):        
        
        """ 
        Initializing an instance of PeakStats.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
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
        #- key: (instrument_id,filename,tolerance), value: peak ratios
        self.peak_ratios = dict()     
        #- key: (instrument_id), value: central wavelengths in mic
        self.central_wav = dict()          
        #- key: (instrument_id,filename,tolerance), value: central fluxes in Jy
        self.central_modelflux = dict() 
        #- key: (instrument_id,filename,tolerance), value: central fluxes in Jy
        self.central_dataflux = dict()    
        #- key: (instrument_id,tolerance,sigma), value: (mean,std)
        self.star_stats = dict()
        #- keeps track of all the data included in this class
        self.depos = dict()              
        self.depos['peak_ratios'] = self.peak_ratios
        self.depos['central_wav'] = self.central_wav
        self.depos['central_modelflux'] = self.central_modelflux
        self.depos['central_dataflux'] = self.central_dataflux
        self.depos['star_stats'] = self.star_stats
        #- key: instrument, value: the upper data noise limits
        self.upper_limit = dict()      
        self.upper_limit['PACS'] = 1.3
        #- key: instrument, value: the lower data noise limits
        self.lower_limit = dict()      
        self.lower_limit['PACS'] = 0.7
        


    def findRatios(self,instrument,tolerance=1,sigma=None,mode='chi2'):
        
        ''' 
        Find the peak to peak ratios of data versus model for a 
        specific instrument and with a tolerance.
        
        The result are saved in the self.peak_ratios dictionary, see 
        __init__.__doc__()
        
        @param instrument: The instrument (such as 'PACS')
        @type instrument: string
        
        @keyword tolerance: a factor that is multiplied with half the 
                            wavelength resolution of data/model when looking 
                            for the peak value around a central wavelength in 
                            the data. If default, no tolerance is added, cannot
                            be smaller than one
                            
                            (default: 1)
        @type tolerance: float
        @keyword sigma: apply sigma clipping to the determined ratios if not 
                        None, in order to deal with outliers
                        
                        (default: none)
        @type sigma: int
        @keyword mode: The statistical method used to calculate "chi squared" 
                       vals can be: chi2, log10_severecut, chi_abs
                       
                       (default: chi2)
        @type mode: string
        
        '''
        
        instrument = instrument.upper()
        print '** Calculating peak (data) to peak (model) ratios... for %s.'\
              %instrument
        if tolerance < 1: tolerance = 1
        for star in self.star_grid:
            this_id = star['LAST_%s_MODEL'%instrument]
            #-reset for these parameters
            self.star_stats[this_id,tolerance,sigma] = []                
            self.central_wav[this_id] = [trans.wavelength*10**(4) 
                                         for trans in star['GAS_LINES']]
            for filename,data_wav,data_flux in \
                    zip(self.instruments[instrument].data_filenames,\
                        self.instruments[instrument].data_wave_list,\
                        self.instruments[instrument].data_flux_list):  
                filename = os.path.split(filename)[1]
                sphinx_file = os.path.join(os.path.expanduser('~'),\
                                           'GASTRoNOoM',self.path_code,\
                                           'stars',self.star_name,\
                                           '%s_results'%instrument,this_id,\
                                           '%s_%s'%('sphinx',filename))
                model_wav = DataIO.readCols(sphinx_file)[0]
                model_flux = DataIO.readCols(sphinx_file)[1]
                if list(model_flux[model_flux < 0]) != []: 
                    print 'There are negative sphinx flux values! They will '+\
                          'not be taken into account.'
                self.central_modelflux[this_id,filename,tolerance] = \
                    [(wav>model_wav[0] \
                            and wav<model_wav[-2] \
                            and model_flux[argmin(abs(model_wav-wav))] > 0) \
                        and model_flux[argmin(abs(model_wav-wav))] \
                        or None
                     for wav in self.central_wav[this_id]]
                d_tmean = self.data_stats[instrument][filename]['tmean']
                d_tstd = self.data_stats[instrument][filename]['tstd']
                
                self.central_dataflux[this_id,filename,tolerance] = \
                    [(wav>=data_wav[0] and wav<=data_wav[-2]) 
                        and max(data_flux\
                                [abs(data_wav-wav)<tolerance\
                                    *(data_wav[argmin(abs(data_wav-wav))+1]\
                                      -data_wav[argmin(abs(data_wav-wav))])\
                                    /2.])
                        or None
                     for wav in self.central_wav[this_id]]
                
                #- If the data point is smaller than the noise, then u don't 
                #- want to include it in statistics: replace by mean + std, but
                #- make it negative to see the difference
                self.central_dataflux[this_id,filename,tolerance] = \
                    [(d >= d_tmean+d_tstd or d is None) \
                        and d \
                        or -1*abs(d_tmean+d_tstd) 
                     for d in self.central_dataflux[this_id,filename,\
                                                    tolerance]]
                
                #- negative peak ratios are possible: these are upper limits, 
                #- where data point was replaced by mean + std
                self.peak_ratios[(this_id,filename,tolerance)] = \
                    [not (m is None or d is None) and m/d or None
                     for m,d in zip(self.central_modelflux[this_id,filename,\
                                                           tolerance],\
                                    self.central_dataflux[this_id,filename,\
                                                          tolerance])]
            
            #- Calculate the chi^2 values according to the requested method
            if mode == 'chi2':
                self.star_stats[this_id,tolerance,sigma,mode] = \
                    (scipy.sqrt(sum(\
                         [(val-1)**2/self.data_stats[instrument]\
                                                    [os.path.split(f)[1]]\
                                                    ['tstd']**2
                           for f in self.instruments[instrument].data_filenames
                           for val in self.getRatios(instrument,this_id,\
                                                     tolerance,sigma=sigma,\
                                                     filename=f)])),\
                     scipy.std(scipy.log10(self.getRatios(instrument,this_id,\
                                                          tolerance,\
                                                          sigma=sigma))))
            
            elif mode == 'chi_abs':
                self.star_stats[this_id,tolerance,sigma,mode] = \
                    (scipy.sqrt(sum(\
                         [abs(val-1)/self.data_stats[instrument]\
                                                    [os.path.split(f)[1]]\
                                                    ['tstd']
                           for f in self.instruments[instrument].data_filenames
                           for val in self.getRatios(instrument,this_id,\
                                                     tolerance,sigma=sigma,\
                                                     filename=f)])),\
                     scipy.std(scipy.log10(self.getRatios(instrument,this_id,\
                                                          tolerance,\
                                                          sigma=sigma))))
            
            elif mode == 'log10_severecut':
                self.star_stats[this_id,tolerance,sigma,mode] = \
                    (scipy.sqrt(sum(\
                         [((val< 1.5 \
                               and abs(scipy.log10(val)) \
                               or scipy.log10(val*10))/self.data_stats\
                                                        [instrument]\
                                                        [os.path.split(f)[1]]\
                                                        ['tstd'])
                           for f in self.instruments[instrument].data_filenames 
                           for val in self.getRatios(instrument,this_id,\
                                                     tolerance,sigma=sigma,\
                                                     filename=f)])),\
                     scipy.std(scipy.log10(self.getRatios(instrument,this_id,\
                                                          tolerance,\
                                                          sigma=sigma))))
            
            #Pearson chi square
            #self.star_stats[this_id,tolerance,sigma] = (sum(\
                #[(d-m)**2/m
                    #for d,m in zip(self.getRatios(instrument,this_id,tolerance,\
                                #sigma=sigma,data_type='central_dataflux'),\
                            #self.getRatios(instrument,this_id,tolerance,\
                                #sigma=sigma,data_type='central_modelflux'))\
                                                                          #]),\
                #scipy.std(scipy.log10(self.getRatios(instrument,this_id,\
                                                    #tolerance,sigma=sigma))))
            ##print '** The logarithmic chi-squared value for %s is:'%this_id
            #print self.star_stats[this_id,tolerance,sigma]
            #print '** The mean and standard deviation of the log10 peak ratios for %s are:'%this_id
            #print self.star_stats[this_id,tolerance,sigma]
            print '** The "chi-squared" and std values for %s are:'%this_id
            print '"chi-squared": %.4f,  std: %.4f'\
                  %self.star_stats[this_id,tolerance,sigma,mode]  
        id_closest_match \
            = self.star_grid[argmin(array([self.star_stats[star['LAST_%s_MODEL'\
                                                                %instrument],\
                                                           tolerance,sigma,\
                                                           mode][0] 
                             for star in self.star_grid 
                             if self.star_stats[star['LAST_%s_MODEL'\
                                                     %instrument],\
                                                tolerance,sigma,mode][0]]))]\
                            ['LAST_%s_MODEL'%instrument]
        id_lowest_spread \
            = self.star_grid[argmin(array([self.star_stats[star['LAST_%s_MODEL'\
                                                                %instrument],\
                                                           tolerance,sigma,\
                                                           mode][1]
                             for star in self.star_grid 
                             if self.star_stats[star['LAST_%s_MODEL'\
                                                     %instrument],\
                                                tolerance,sigma,mode][0]]))]\
                            ['LAST_%s_MODEL'%instrument]
        if not id_closest_match:
            return
        print '** The minimum statistical value is found for:'
        print id_closest_match
        print '"chi-squared": %.4f, std: %.4f'\
              %self.star_stats[id_closest_match,tolerance, sigma,mode]
        print '** The lowest spread is found for:'
        print id_lowest_spread
        print '"chi-squared": %.4f, std: %.4f'\
              %self.star_stats[id_lowest_spread,tolerance, sigma,mode]
        #id_closest_match = self.star_grid[argmin(array([abs(self.star_stats[star['LAST_%s_MODEL'%instrument],tolerance,sigma][0]) for star in self.star_grid]))]['LAST_%s_MODEL'%instrument]
        #id_lowest_spread = self.star_grid[argmin(array([self.star_stats[star['LAST_%s_MODEL'%instrument],tolerance,sigma][1] for star in self.star_grid]))]['LAST_%s_MODEL'%instrument]
        #print '** The best match (ie mean log10 closest to 0) is found for:'
        #print id_closest_match
        #print 'with mean: %s, and std: %s'%(self.star_stats[id_closest_match,tolerance, sigma])
        #print '** The lowest spread is found for:'
        #print id_lowest_spread
        #print 'with mean: %s, and std: %s'%(self.star_stats[id_lowest_spread,tolerance, sigma])
        print '***********************************'

##############################
######### GET RATIOS #########
##############################
                                                                                
    def getRatios(self,instrument,this_id,tolerance,data_type='peak_ratios',\
                  return_negative=0,sigma=None,filename=None,return_clipped=0):
        
        '''
        Return all ratios for all filenames for one instrument and one 
        tolerance, including only the values that are not None.
        
        @param instrument: The name of the instrument (such as 'PACS')
        @type instrument: string
        @param this_id: The requested model id
        @type this_id: string 
        @param tolerance: The tolerance from self.findRatios()
        @type tolerance: float
        
        @keyword data_type: the type of data returned, one of 'peak_ratios' or
                            'central_wav'
                            
                            (default:'peak_ratios')
        @type data_type: string
        @keyword return_negative: only return the negative peak ratios or 
                                  equivalent (upper limits)
                                  
                                  (default: 0)
        @type return_negative: bool
        @keyword sigma: if None, this parameter is ignored, 
                        if an int only the values within this sigma interval 
                        are included
                        
                        (default: None)
        @type sigma: int
        @keyword return_clipped: return a list of the values that were sigma 
                                 clipped, only relevant if sigma is not None
        
                                 (default: 0)
        @type return_clipped: bool
        @keyword filename: the filename for which you want to return the list, 
                           if None, all filenames are used and the lists are 
                           merged into one
                           
                           (default: None)
        @type filename: string
        
        @return: The ratios requested
        @rtype: array
        
        '''
        
        filenames = filename is None \
                        and self.instruments[instrument].data_filenames \
                        or [os.path.split(filename)[1]]
        if sigma is None or return_negative:
            if data_type == 'central_wav':
                return array([val2     
                              for filename in filenames
                              for val,val2 in \
                                    zip(self.peak_ratios\
                                            [(this_id,\
                                              os.path.split(filename)[1],\
                                              tolerance)],\
                                        self.depos[data_type][this_id])
                              if val <> None \
                                    and ((return_negative and val < 0) \
                                    or (not return_negative and val > 0))])
            else:
                return array([val2  
                              for filename in filenames
                              for val,val2 in \
                                    zip(self.peak_ratios\
                                            [(this_id,\
                                              os.path.split(filename)[1],
                                              tolerance)],\
                                        self.depos[data_type]\
                                            [(this_id,\
                                              os.path.split(filename)[1],\
                                              tolerance)])
                              if val <> None \
                                    and ((return_negative and val < 0) \
                                    or (not return_negative and val > 0))])
        else:
            if data_type == 'central_wav':
                returned_vals = array([val2
                                       for filename in filenames
                                       for val,val2 in \
                                            zip(self.peak_ratios\
                                                 [(this_id,\
                                                   os.path.split(filename)[1],\
                                                   tolerance)],\
                                                self.depos[data_type][this_id])
                                       if val <> None and val > 0])
            else:
                returned_vals = array([val2     
                                       for filename in filenames
                                       for val,val2 in \
                                            zip(self.peak_ratios\
                                                 [(this_id,\
                                                   os.path.split(filename)[1],\
                                                   tolerance)],\
                                                self.depos[data_type]\
                                                 [(this_id,\
                                                   os.path.split(filename)[1],\
                                                   tolerance)])
                                       if val <> None and val > 0])
            ratios = array([val     
                            for filename in filenames
                            for val in self.peak_ratios\
                                            [(this_id,\
                                              os.path.split(filename)[1],\
                                              tolerance)]
                            if val <> None and val > 0])
            this_mean = scipy.mean(scipy.log10(ratios))
            this_std = scipy.std(scipy.log10(ratios))
            if return_clipped:
                returned_vals = array([val2 
                                       for val,val2 in zip(ratios,\
                                                           returned_vals) 
                                       if scipy.log10(val) > this_mean \
                                                + this_std*sigma \
                                            or scipy.log10(val) < this_mean \
                                                - this_std*sigma])
            else:
                returned_vals = array([val2 
                                       for val,val2 in zip(ratios,\
                                                           returned_vals) 
                                       if scipy.log10(val) <= this_mean \
                                                + this_std*sigma \
                                            and scipy.log10(val) >= this_mean \
                                                - this_std*sigma])
            return returned_vals



    def plotRatioWav(self,instrument,inputfilename,tolerance=1,sigma=None,\
                     mode='chi2'):
        
        '''
        Plot peak ratios as a function of their central wavelength.
        
        @param instrument: The name of the instrument (such as 'PACS')
        @type instrument: string
        @param inputfilename: the input filename for the grid, which will be 
                              attached to the final plot
        @type inputfilename: string

        @keyword tolerance: The tolerance from self.findRatios()
                        
                          (default: 1)
        @type tolerance: float
        @keyword sigma: if None, this parameter is ignored, 
                        if an int only the values within this sigma interval 
                        are included
                        
                        (default: None)
        @type sigma: int
        @keyword mode: The statistical method used to calculate "chi squared" 
                       vals. Can be: chi2, log10_severecut, chi_abs
                       
                       (default: chi2)
                       
        @type mode: string
        
        '''
        
        instrument = instrument.upper()
        plot_filenames = []
        this_grid = self.sortStarGrid(instrument,tolerance,sigma,mode)
        plot_filenames = []
        for star in this_grid:
            this_id = star['LAST_%s_MODEL'%instrument]
            statvals = self.star_stats[this_id,tolerance,sigma,mode]
            line_types_inc = ['or','ob','oc']
            #- upper limits in data
            line_types_upper = ['dr','db','dc']        
            lp = [] 
            waves = []
            ratios = []
            
            #- the ratios included in the statistics
            this_wav_inc = self.getRatios(instrument=instrument,\
                                          tolerance=tolerance,sigma=sigma,\
                                          data_type='central_wav',\
                                          this_id=this_id)
            this_ratio_inc = self.getRatios(instrument=instrument,\
                                            tolerance=tolerance,sigma=sigma,\
                                            this_id=this_id)
            if list(this_wav_inc):    
                waves.append(this_wav_inc)
                ratios.append(this_ratio_inc)
                lp.append('ob')
            
            #- the ratios replaced by upper limit if data point is in the noise
            this_wav_upper = self.getRatios(instrument=instrument,\
                                            tolerance=tolerance,sigma=sigma,\
                                            data_type='central_wav',\
                                            this_id=this_id,return_negative=1)
            this_ratio_upper = self.getRatios(instrument=instrument,\
                                              tolerance=tolerance,sigma=sigma,\
                                              return_negative=1,\
                                              this_id=this_id)
            this_ratio_upper = [abs(r) for r in this_ratio_upper]
            if list(this_wav_upper):
                waves.append(this_wav_upper)
                ratios.append(this_ratio_upper)
                lp.append('dr')
            
            #- the ratios excluded through sigma clipping
            if sigma <> None:    
                this_wav_exc = self.getRatios(instrument=instrument,\
                                              tolerance=tolerance,sigma=sigma,\
                                              data_type='central_wav',\
                                              this_id=this_id,return_clipped=1)
                this_ratio_exc = self.getRatios(instrument=instrument,\
                                                tolerance=tolerance,\
                                                sigma=sigma,return_clipped=1,\
                                                this_id=this_id)
                if list(this_wav_exc):
                    waves.append(this_wav_exc)
                    ratios.append(this_ratio_exc)
                    lp.append('xk')
            
            #- prepping input for the plot command
            xmin = min([min(x) for x in waves])
            xmax = max([max(x) for x in waves])
            waves.extend([[0.5*xmin,1.5*xmax]]*3)
            ratios.extend([[1,1],\
                           [self.lower_limit[instrument],\
                            self.lower_limit[instrument]],\
                           [self.upper_limit[instrument],\
                            self.upper_limit[instrument]]])
            lp.extend(['-k','--k','--k'])
            plot_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                         self.path_code,'stars',\
                                         self.star_name,\
                                         '%s_results'%instrument,this_id,\
                                         'ratio_wav_%s_sig%s_tol%s'\
                                         %(this_id,str(sigma),str(tolerance)))
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
            plot_title = '%s: sigma '%this_id.replace('_','\_') + \
                         '%s, '%(sigma <> None and str(sigma) or str(None)) + \
                         'tol %i, '%int(tolerance) + \
                         '"chi-squared" %.4f, '%statvals[0] + \
                         'std: %.4f'%statvals[1]
            plot_filenames.append(Plotting2.plotCols(\
                    filename=plot_filename,x=waves,y=ratios,\
                    yaxis=r'$F_{\nu,p,m}/F_{\nu,p,d}$',\
                    plot_title=plot_title,labels=labels,\
                    xlogscale=0,ylogscale=1,line_types=lp,xmin=xmin*0.9,\
                    xmax=xmax*1.03,figsize=(10.*scipy.sqrt(2.), 10.),\
                    linewidth=2,fontsize_title=20,fontsize_label=16))
        new_filename = os.path.join(\
                os.path.expanduser('~'),'GASTRoNOoM',self.path_code,'stars',\
                self.star_name,'%s_results'%instrument,\
                'ratio_wav_sig%s_tol%s_%s.pdf'\
                %(str(sigma),str(tolerance),\
                  os.path.splitext(os.path.split(inputfilename)[1])[0]))
        DataIO.joinPdf(old=plot_filenames,new=new_filename)
        print '** Stat plots can be found at:'
        print new_filename
        print '***********************************'



    def sortStarGrid(self,instrument,tolerance,sigma,mode):
        
        '''
        Return a sorted list of the star grid in this instance according to 
        best fit parameters (chi-squared).
        
        @param instrument: The name of the instrument (such as 'PACS')
        @type instrument: string
        @param tolerance: The tolerance from self.findRatios()
        @type tolerance: float
        @param sigma: if None, this parameter is ignored, 
                      if an int only the values within this sigma interval 
                      are included
        @type sigma: int
        @param mode: The statistical method used to calculate "chi squared" 
                     vals. Can be: chi2, log10_severecut, chi_abs
        @type mode: string
        
        '''
        
        return sorted(self.star_grid,\
                      key=lambda x: self.star_stats\
                                        [x['LAST_%s_MODEL'%instrument],\
                                         tolerance,sigma,mode])
                                         