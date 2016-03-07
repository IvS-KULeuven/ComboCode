# -*- coding: utf-8 -*-

"""
A plotting environment for dust information such as SEDs and all that is 
associated with that.

Author: R. Lombaert

"""

import os
import subprocess
import pyfits
import glob
from scipy import array
import numpy as np

import cc.path
from cc.data import Sed
from cc.plotting.PlottingSession import PlottingSession
from cc.plotting import Plotting2
from cc.tools.io import DataIO, KappaReader
from cc.modeling.objects import Star
from cc.modeling.codes import MCMax
from cc.modeling.tools import Profiler,Reddening 

from ivs.sed.model import synthetic_flux
from ivs.units import conversions 

class PlotDust(PlottingSession):
    
    """ 
    Plotting environment for SEDs and all dust parameters.
    
    """
    
    def __init__(self,star_name='model',sed=None,path_mcmax='',\
                 inputfilename=None,fn_add_star=0):
        
        '''
        Initializing PlotDust session.
        
        @keyword star_name: name of the star from Star.dat, use default only 
                            when never using any star model specific things 
                                  
                            (default: "model")
        @type star_name: string
        @keyword path_mcmax: Output modeling folder in MCMax home folder
        
                             (default: '')
        @type path_mcmax: string
        @keyword inputfilename: name of inputfile that is also copied to the 
                                output folder of the plots, 
                                if None nothing is copied
                                
                                (default: None)
        @type inputfilename: string
        @keyword sed: an Sed object needed for plotting the SED. None if not 
                      plotting an Sed, but just dust parameters such as \
                      temperature
                            
                      (default: None)
        @type sed: Sed()
        @keyword fn_add_star: Add the star name to the requested plot filename.
                              Only relevant if fn_plt is given in a sub method.
                              
                              (default: 0)
        @type fn_add_star: bool
        
        '''

        super(PlotDust, self).__init__(star_name=star_name,\
                                       path=path_mcmax,\
                                       code='MCMax',\
                                       inputfilename=inputfilename,\
                                       fn_add_star=fn_add_star)
        #-- Convenience path
        cc.path.mout = os.path.join(cc.path.mcmax,self.path)
        self.sed = sed



    def plotSed(self,star_grid=[],fn_plt='',cfg='',iterative=0,no_models=0,\
                show_phot_filter=0,**kwargs):
        
        """ 
        Creating an SED with 0, 1 or more models and data. 
        
        Includes data preparation on the spot.
        
        Additional plotCols keywords can be passed through kwargs. They 
        overwrite cfg input.
        
        @keyword star_grid: list of Star() models to plot. If star_grid is [], 
                            only data are plotted.
                            
                            (default: [])
        @type star_grid: list[Star()]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                        
                      (default: '')
        @type cfg: string
        @keyword iterative: add an extra suffix to the filename for each 
                            iteratively calculated model, with this number 
                            giving the model muber (index in star_grid), 
                            0 if not used.
                                  
                            (default: 0)
        @type iterative: int
        @keyword no_models: Only show data.
                                  
                            (default: 0)
        @type no_models: bool
        @keyword show_phot_filter: Show the wavelength band of the photometric
                                   filters as an x error bar on the model phot
                                   
                                   (default: 0)
        @type show_phot_filter: bool
        
        """
        
        if self.sed is None:
            print 'No dsed given in Path.dat. Cannot plot SED. Aborting...'
            return
        print '***********************************'
        print '** Creating SED plot.'
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('no_models'):
            no_models = cfg_dict['no_models']
        if cfg_dict.has_key('show_phot_filter'):
            show_phot_filter = bool(cfg_dict['show_phot_filter'])
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
        
        data_labels = dict([(dt,(n,ls))
                            for n,dt,ls in zip(DataIO.getInputData(path=cc.path.usr,\
                                                        keyword='PLOT_NAMES',\
                                                        filename='Sed.dat',\
                                                        remove_underscore=1),\
                                               DataIO.getInputData(path=cc.path.usr,\
                                                        keyword='DATA_TYPES',\
                                                        filename='Sed.dat'),\
                                               DataIO.getInputData(path=cc.path.usr,\
                                                        keyword='LINE_TYPES',\
                                                        filename='Sed.dat'))])

        plot_title='SED %s'%self.star_name_plots
       
        #- prepare and collect data, keytags and line types
        keytags = []
        data_x = []
        data_y = []
        data_xerr = []
        data_yerr = []
        line_types = []
        for (dt,fn),tdata in sorted([dset
                                     for dset in self.sed.data.items()
                                     if 'PHOT' not in dset[0][0].upper()]):
             keytags.append(data_labels[dt][0])
             data_x.append(tdata[0])
             data_y.append(tdata[1])
             #data_err.append(tdata[2])
             #-- For now, no error bars for spectra. 
             data_xerr.append(None)
             data_yerr.append(None)
             line_types.append(data_labels[dt][1])
        
        for (dt,fn),(w,f,err) in sorted([dset
                                         for dset in self.sed.data.items()
                                         if 'PHOT' in dset[0][0].upper()]):
             keytags.append(data_labels[dt][0])
             data_x.append(w)
             data_xerr.append(None)
             data_y.append(f)
             data_yerr.append(err)
             line_types.append(data_labels[dt][1])
        
        #-- Set line types manually, so the photometry can be plotted with the 
        #   right colors.
        elp = Plotting2.getLineTypes()
        elp = [lp for lp in elp if lp not in line_types]
        
        #- Collect model data as well as keytags and set line types
        model_ids_mcm = [s['LAST_MCMAX_MODEL'] 
                         for s in star_grid
                         if s['LAST_MCMAX_MODEL']]
        
        #- Only if the model_ids list is not empty, MCMax models are available
        #- Otherwise the ray tracing keyword is unnecessary.
        if no_models:
            model_ids_mcm = []
        
        wmodels = []
        fmodels = []
        for model_id,s in zip(model_ids_mcm,star_grid):
            dpath = os.path.join(cc.path.mout,'models',model_id)
            fn_spec = 'spectrum{:04.1f}.dat'.format(s['RT_INCLINATION'])
            w,f = MCMax.readModelSpectrum(dpath,s['RT_SPEC'],fn_spec)
            if s['REDDENING']:
                print 'Reddening models to correct for interstellar extinction.'
                ak = self.sed.getAk(s['DISTANCE'],s['REDDENING_MAP'],\
                                    s['REDDENING_LAW'])
                f = Reddening.redden(w,f,ak,law=s['REDDENING_LAW'])
            wmodels.append(w)
            fmodels.append(f)
        
        for i,(w,f,model_id) in enumerate(zip(wmodels,fmodels,model_ids_mcm)):
            data_x.append(w)
            data_y.append(f)
            data_xerr.append(None)
            data_yerr.append(None)
            keytags.append(model_id.replace('_','\_'))
            line_types.append(elp[i])
        
        if self.sed.photbands.size:
            filts = self.sed.filter_info  
            for i,(w,f) in enumerate(zip(wmodels,fmodels)):
                mphot = Sed.calcPhotometry(w,f,self.sed.photbands)
                data_x.append(self.sed.photwave)
                data_y.append(mphot)
                if show_phot_filter:
                    data_xerr.append([self.sed.photwave-filts.wlower,\
                                      filts.wupper-self.sed.photwave])
                else:
                    data_xerr.append(None)
                data_yerr.append(None)
                line_types.append('o%s'%Plotting2.splitLineType(elp[i])[1])

        keytags = [tag.replace('#','') for tag in keytags]
        extra_pars = dict()
        extra_pars['line_types'] = line_types
        extra_pars['keytags'] = keytags
        extra_pars['ymax'] = 1.3*max([max(dy[np.isfinite(dy)]) 
                                      for dy in data_y])
        extra_pars['ymin'] = 0.5*min([min(dy[np.isfinite(dy)]) 
                                      for dy in data_y])
        extra_pars['xmin'] = 2
        extra_pars['xmax'] = 200
        extra_pars['fontsize_key'] = 16
        extra_pars['xlogscale'] = 1
        extra_pars['ylogscale'] = 0
        extra_pars['xerr_capsize'] = 0
        extra_pars['xerr_linewidth'] = 4
        extra_pars.update(cfg_dict)
        extra_pars.update(kwargs)
        
        #-- Set filename plot
        pfn = fn_plt if fn_plt else 'SED'
        suff = 'iterative_{}'.format(iterative) if iterative else ''
        pfn = self.setFnPlt(pfn,fn_suffix=suff)
        
        pfn = Plotting2.plotCols(x=data_x,y=data_y,xerr=data_xerr,\
                                 yerr=data_yerr,filename=pfn,**extra_pars)
        print '** Your SED plots can be found at:'
        print pfn
        print '***********************************'
         
         
         
    def plotCorrflux(self,star_grid=[],fn_plt='',cfg='',no_models=0):
        
        """ 
        Plot correlated fluxes with 0, 1 or more models and data. 
        
        Includes data preparation on the spot.
        
        @keyword star_grid: list of Star() models to plot. If star_grid is [], 
                            only data are plotted.
                            
                            (default: [])
        @type star_grid: list[Star()]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                        
                      (default: '')
        @type cfg: string
        @keyword no_models: Only show data.
                                  
                            (default: 0)
        @type no_models: bool
        
        """
        
        print '***********************************'
        print '** Creating Correlated Fluxes plot.'
        if not cc.path.dcflux:
            print 'No dcflux given in Path.dat. No data are plotted.'
            
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('no_models'):
            no_models = cfg_dict['no_models']
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
        
        #-- Select MIDI data. Assumes baseline at the end of the filename.
        ssd = os.path.join(cc.path.dcflux,self.star_name,\
                           '_'.join([self.star_name,'MIDI','*.fits']))
        files = [os.path.splitext(gi)[0] for gi in glob.glob(ssd)]
        ggd = dict([(float(gi.split('_')[-1].strip('m')),gi+'.fits') 
                    for gi in files])

        #-- Read the models
        models = []
        if not no_models:
            for s in star_grid:
                model_id = s['LAST_MCMAX_MODEL']
                dpath = os.path.join(cc.path.mout,'models',model_id)
                fn_vis = 'visibility{:04.1f}.dat'.format(s['RT_INCLINATION'])
                model = MCMax.readVisibilities(dpath=dpath,fn_vis=fn_vis)
                models.append(model)
            real_models = [model for model in models if model]
            if not real_models: 
                no_models = 1 
                baselines = sorted(ggd.keys())
            else:
                baselines = sorted(real_models[0]['baseline'].keys())
        
        #-- prepare and collect data, keytags and line types
        data = []
        for bl in baselines:
            ddict = dict()
            data.append(ddict)

            #-- Extract data from the fits file
            if ggd.has_key(bl): 
                dfits = pyfits.open(ggd[bl])
                x = 1e6*dfits['OI_WAVELENGTH'].data['EFF_WAVE'][::-1]
                ddict['x'] = [x]
                ddict['y'] = [dfits['OI_VIS'].data['VISAMP'][0]][::-1]
                ddict['yerr'] = [dfits['OI_VIS'].data['VISAMPERR'][0]][::-1]
                dfits.close()
            else:
                ddict['x'] = []
                ddict['y'] = []
                ddict['yerr'] = []
                                
            if no_models:
                continue
           
            #-- Extract models from the model folders
            for model in models:
                ddict['yerr'].append(None)
                if not model: 
                    ddict['x'].append(np.empty(0))
                    ddict['y'].append(np.empty(0))
                    continue
                ddict['x'].append(model['wavelength'])
                ddict['y'].append(model['flux']*model['baseline'][bl])
            
            #-- Set some plot limits
            ddict['xmin'] = 8
            ddict['xmax'] = 13
            ddict['ymin'] = -0.1
            ddict['labels'] = [('MIDI %.1f m'%bl,0.05,0.9)]
            #-- Wavelength limits between 8 and 13 micron, limits of the N band
            #   atmospheric transmission. Outside these ranges, the flux is not
            #   relevant
            ddict['ymax'] = 1.1*max([max(iy[(ix<=13.)*(ix>=8.)]) 
                                     for ix,iy in zip(ddict['x'],ddict['y'])
                                     if iy.size])
        kwargs = dict()
        kwargs['keytags'] = ['MIDI'] if ggd else []
        if not no_models:
            kwargs['keytags'].extend([s['LAST_MCMAX_MODEL'].replace('_','\_') 
                                      for s in star_grid])
        kwargs['xaxis'] = '$\lambda$ ($\mu$m)'
        kwargs['yaxis'] = 'Corr.~FLux (Jy)'
        kwargs['dimensions'] = (1,len(data)+1)
        kwargs['figsize'] = (10,15)
        kwargs['fontsize_axis'] = 20
        kwargs['fontsize_ticklabels'] = 20
        kwargs['fontsize_key'] = 18
        kwargs['fontsize_label'] = 14
        kwargs['linewidth'] = 3
        kwargs['cfg'] = cfg_dict
        kwargs['extension'] = '.pdf'
        kwargs['hspace'] = 0.3
        kwargs['ws_bot'] = 0.01
        kwargs['ws_top'] = 0.99
        kwargs['ws_left'] = 0.10
        kwargs['ws_right'] = 0.98
        
        #-- Set filename plot
        pfn = fn_plt if fn_plt else 'CorrFlux'
        pfn = self.setFnPlt(pfn)
        
        pfn = Plotting2.plotTiles(data=data,filename=pfn,**kwargs)
        print '** Your Correlated Flux plots can be found at:'
        print pfn
        print '***********************************'
                 
                 
                 
    def plotVisibilities(self,star_grid=[],fn_plt='',cfg='',no_models=0):
        
        """ 
        Plot visibilities as a function of baseline.
        
        Wavelengths plotted are what is requested in the ray tracing

        Includes data preparation on the spot.
        
        Data location is that of correlated flux (for the visibilities), but 
        also requires an sed object to retrieve the MIDI spectrum. If one of 
        them is not available, models are be plotted without data.

        @keyword star_grid: list of Star() models to plot. If star_grid is [], 
                            only data are plotted.
                            
                            (default: [])
        @type star_grid: list[Star()]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                        
                      (default: '')
        @type cfg: string
        @keyword no_models: Only show data.
                                  
                            (default: 0)
        @type no_models: bool
        @keyword fn_add_star: Add the star name to the requested plot filename.
        
                              (default: 1)
        @type fn_add_star: bool
        
        """
        
        if not cc.path.dcflux:
            print 'No dcflux given in Path.dat. Aborting...'
            return
        
        print '***********************************'
        print '** Creating Visibilities plot.'
        if not self.sed or 'MIDI' not in self.sed.data_types:
            print 'No dsed given in Path.dat or no MIDI spectral data found.'

        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('no_models'):
            no_models = cfg_dict['no_models']
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')

        #-- Read the models. Wavelengths are taken from the ray-tracing output
        models = []
        if not no_models:
            for s in star_grid:
                model_id = s['LAST_MCMAX_MODEL']
                dpath = os.path.join(cc.path.mout,'models',model_id)
                fn_vis = 'basevis{:04.1f}.dat'.format(s['RT_INCLINATION'])
                model = MCMax.readVisibilities(dpath=dpath,fn_vis=fn_vis)
                models.append(model)
            real_models = [model for model in models if model]
            if not real_models: 
                no_models = 1 
            else:
                wavelengths = sorted(real_models[0]['wavelength'].keys())
        
        if no_models: wavelengths = (8.,10.,13.)
        
        #-- Grab the MIDI spectrum
        if self.sed and 'MIDI' in self.sed.data_types:
            fn = self.sed.data_filenames[self.sed.data_types.index('MIDI')]
            midi_flux = self.sed.data[('MIDI',fn)][1]
            midi_err = self.sed.data[('MIDI',fn)][2]
            midi_relerr = (midi_err/midi_flux)**2
        
            #-- Select MIDI data. Assumes baseline at the end of the filename.
            ssd = os.path.join(cc.path.dcflux,self.star_name,\
                               '_'.join([self.star_name,'MIDI','*.fits']))
            files = [os.path.splitext(gi)[0] for gi in glob.glob(ssd)]
            ggd = dict([(float(gi.split('_')[-1].strip('m')),gi+'.fits') 
                        for gi in files])
        else: 
            ggd = dict()
                    
        #-- Collect MIDI data from the fits file and calculate visibilities
        ddf = dict()
        for k,v in sorted(ggd.items()):
            ddf[k] = dict()
            dfits = pyfits.open(v)
            
            #-- Read the wavelength
            cwave = 1e6*dfits['OI_WAVELENGTH'].data['EFF_WAVE'][::-1]
            
            #-- Read flux + err and select the right range
            cflux = dfits['OI_VIS'].data['VISAMP'][0][::-1]
            cflux = cflux[(cwave<=13.)*(cwave>=8.)]
            cflux_err = dfits['OI_VIS'].data['VISAMPERR'][0][::-1]
            cflux_err = cflux_err[(cwave<=13.)*(cwave>=8.)]
            
            #-- The visibilities are correlated flux divided by real flux
            ddf[k]['y'] = cflux/midi_flux
            
            #-- Error propagation
            cflux_relerr = (cflux_err/cflux)**2
            yerr = np.sqrt(midi_relerr + cflux_relerr)*cflux/midi_flux
            ddf[k]['yerr'] = yerr
            
            #-- Wavelength grid
            ddf[k]['x'] = cwave[(cwave<=13.)*(cwave>=8.)]
            dfits.close()
        
        #-- prepare and collect plot data, keytags and line types
        data = []
        for w in wavelengths:
            ddict = dict()
            data.append(ddict)
            
            #-- Set the plot x and y if data are available 
            if ddf:
                bls = [k for k in sorted(ddf.keys())]
                ddict['x'] = [[bl for bl in bls]]
                ddict['y'] = [[ddf[bl]['y'][np.argmin(abs(ddf[bl]['x']-w))]
                               for bl in bls]]
                ddict['yerr'] = [[ddf[bl]['yerr'][np.argmin(abs(ddf[bl]['x']-w))]
                                  for bl in bls]]
            else:
                ddict['x'] = []
                ddict['y'] = []
                ddict['yerr'] = []
            
            #-- Set labels
            ddict['labels'] = [('MIDI %s $\\mu$m'%w,0.85,0.9)]
            
            if no_models:
                continue
            
            #-- Extract models from the model folders
            for model in models:
                ddict['yerr'].append(None)
                if not model: 
                    ddict['x'].append(np.empty(0))
                    ddict['y'].append(np.empty(0))
                    continue
                ddict['x'].append(model['baseline'])
                ddict['y'].append(model['wavelength'][w])
                            
        kwargs = dict()
        kwargs['keytags'] = ['MIDI'] if ddf else []
        if not no_models:
            kwargs['keytags'].extend([s['LAST_MCMAX_MODEL'].replace('_','\_') 
                                      for s in star_grid])
        kwargs['xaxis'] = 'Baseline (m)'
        kwargs['yaxis'] = 'Visibility'
        kwargs['dimensions'] = (1,len(data)+1)
        kwargs['figsize'] = (10,15)
        kwargs['fontsize_axis'] = 20
        kwargs['fontsize_ticklabels'] = 20
        kwargs['fontsize_key'] = 18
        kwargs['fontsize_label'] = 14
        kwargs['linewidth'] = 3
        kwargs['cfg'] = cfg_dict
        kwargs['extension'] = '.pdf'
        kwargs['hspace'] = 0.3
        kwargs['ws_bot'] = 0.01
        kwargs['ws_top'] = 0.99
        kwargs['ws_left'] = 0.10
        kwargs['ws_right'] = 0.98

        #-- Set filename plot
        pfn = fn_plt if fn_plt else 'Visibilities'
        pfn = self.setFnPlt(pfn)

        pfn = Plotting2.plotTiles(data=data,filename=pfn,**kwargs)
        print '** Your Correlated Flux plots can be found at:'
        print pfn
        print '***********************************'

                                   

    def plotDens(self,star_grid=[],models=[],fn_plt='',unit='cm',cfg=''):
        
        """ 
        Plotting the temperature stratification of the dust.
        
        All models are shown in one plot.
        
        @keyword star_grid: parameter sets, if [], the parameter
                            sets are determined from the model ids
        
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model_ids, if [], the parameter sets are expected
                         in star_grid
                         
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword unit: The unit of the plotted radial grid. Can be 'cm','rstar',
                       'au', 'm'
        
                       (default: 'cm')
        @type unit: str      
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        """
        
        print '***********************************'
        print '** Starting to plot the dust density profile.'
        if not star_grid and not models:
            print 'Input is undefined. Aborting.'
            return        
        elif not star_grid and models:
            star_grid = self.makeMCMaxStars(models=models)
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('unit'):
            unit = cfg_dict['unit']
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')

        rads = [s.getDustRad(unit=unit) for s in star_grid]
        denss = [s.getDustDensity() for s in star_grid]
        keys = [s['LAST_MCMAX_MODEL'].replace('_','\_') for s in star_grid]
        
        ppars = dict()
        ppars['yaxis'] = r'$\rho_\mathrm{d}\ \mathrm{(g cm}^{-3}\mathrm{)}$'
        if unit == 'rstar': ppars['xaxis'] = '$r\ \mathrm{(R}_\star\mathrm{)}$'
        elif unit =='au': ppars['xaxis'] = '$r\ \mathrm{(AU)}$'
        elif unit == 'm': ppars['xaxis'] = '$r\ \mathrm{(m)}$'
        else: ppars['xaxis'] = '$r\ \mathrm{(cm)}$'
        ppars['xlogscale'] = 1
        ppars['ylogscale'] = 1
            
        #-- Set plot filename
        pfn = fn_plt if fn_plt else 'dens'
        pfn = self.setFnPlt(pfn)
        
        pfn = Plotting2.plotCols(x=rads,y=denss,filename=pfn,\
                                 keytags=keys,cfg=cfg_dict,**ppars)
        print '** Your plot can be found at:'
        print pfn
        print '***********************************'
            

                                    
    def plotTemp(self,star_grid=[],models=[],fn_plt='',power=[],unit='cm',\
                 cfg=''):
        
        """ 
        Plotting the temperature stratification of the dust.
        
        All models are shown in one plot.
        
        @keyword star_grid: parameter sets, if [], the parameter
                            sets are determined from the model ids
        
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model_ids, if [], the parameter sets are expected
                         in star_grid
                         
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword power: A list of values for s in below formula. If [] no power
                        law is included. Power law parameters  are taken from 
                        star_grid[0].
                                
                        See Thesis p32, where power is s in 
                        T(r) = T_eff*(2*r/R_STAR)**(-2/(4+s)). This value is 
                        typically 1. 
                
                        (default: [])
        @type power: list
        @keyword unit: The unit of the plotted radial grid. Can be 'cm','rstar',
                       'au', 'm'
        
                       (default: 'cm')
        @type unit: str        
        @keyword cfg: path to the Plotting2.plotCols config file. If default,
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        """
        
        print '***********************************'
        print '** Starting to plot dust temperature stratification.'
        if not star_grid and not models:
            print 'Input is undefined. Aborting.'
            return        
        elif not star_grid and models:
            star_grid = self.makeMCMaxStars(models=models)
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('power'):
            power = cfg_dict['power']
        if cfg_dict.has_key('unit'):
            unit = cfg_dict['unit']
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')

        rads = []
        temps = []
        keytags = []
        for star in star_grid:
            rad = star.getDustRad(unit=unit)
            temp,key = star.getDustTemperature(add_key=1) 
            rads.append(rad)
            temps.append(temp)
            keytags.append(key)
        
        #-- Add power laws if requested
        for s in power:
            rad_rstar = star_grid[0].getDustRad(unit='rstar')
            rad = star_grid[0].getDustRad(unit=unit)
            tstar = star_grid[0]['T_STAR']
            temp,key = Profiler.dustTemperaturePowerLaw(rad=rad_rstar,\
                                                        add_key=1,\
                                                        tstar=tstar,s=s)
            rads.append(rad)
            temps.append(temp)
            keytags.append(key)
        
        ppars = dict()
        ppars['yaxis'] = '$T_\mathrm{d}$ (K)'
        if unit == 'rstar': ppars['xaxis'] = '$r\ \mathrm{(R}_\star\mathrm{)}$'
        elif unit =='au': ppars['xaxis'] = '$r\ \mathrm{(AU)}$'
        elif unit == 'm': ppars['xaxis'] = '$r\ \mathrm{(m)}$'
        else: ppars['xaxis'] = '$r\ \mathrm{(cm)}$'
        ppars['plot_title'] = 'Average Dust Temperature Stratification for %s'\
                              %(self.star_name_plots)
        
        #-- Set plot filename
        pfn = fn_plt if fn_plt else 'Td_avg'
        pfn = self.setFnPlt(pfn)
        
        pfn = Plotting2.plotCols(x=rads,y=temps,filename=pfn,\
                                 key_location=(0.05,0.05),cfg=cfg_dict,\
                                 xlogscale=1,ylogscale=1,\
                                 keytags=keytags,**ppars)
        print '** Your plots can be found at:'
        print pfn
        print '***********************************'
            


    def plotTempSpecies(self,star_grid=[],models=[],fn_plt='',include_total=1,\
                        unit='cm',power=[],cfg=''):
        
        """ 
        Plotting the temperature stratification of the dust for the species 
        separately, per model.
        
        @keyword star_grid: parameter sets, if [], the parameter
                            sets are determined from the model ids
        
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model_ids, if [], the parameter sets are expected
                         in star_grid
                         
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword include_total: Include the sum of all temperature profiles as 
                                well for comparison. 
                                        
                                (default: 0)
        @type include_total: bool
        @keyword unit: The unit of the plotted radial grid. Can be 'cm','rstar',
                       'au', 'm'
        
                       (default: 'cm')
        @type unit: str      
        @keyword power: A list of values for s in below formula. If [] no power
                        law is included. Power law parameters are taken from 
                        star_grid[0].
                                
                        See Thesis p32, where power is s in 
                        T(r) = T_eff*(2*r/R_STAR)**(-2/(4+s)). This value is 
                        typically 1. 
                
                        (default: [])
        @type power: list        
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
                
        """            
        
        print '***********************************'
        print '** Starting to plot dust temperature for separate species.'
        if not star_grid and not models:
            print 'Input is undefined. Aborting.'
            return        
        elif not star_grid and models:
            star_grid = self.makeMCMaxStars(models=models,id_type='MCMax')
            raise IOError('Reading dust species temperatures from a model id'+\
                          ' list only, not yet implemented.')
            #- Requires star.dust_list and T_CONTACT to be taken from the log file. 
            #- It's possible, but needs some programming

        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('power'):
            power = cfg_dict['power']
        if cfg_dict.has_key('unit'):
            unit = cfg_dict['unit']
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')

        plot_filenames = []
        for star in star_grid:
            if not int(star['T_CONTACT']):
                rads = [star.getDustRad(species=species,unit=unit)
                        for species in star.getDustList()]
                temps = [star.getDustTemperature(species=species)
                         for species in star.getDustList()]
                rads = [r[t<=star['T_DES_%s'%sp]] 
                        for r,t,sp in zip(rads,temps,star.getDustList())]
                temps = [t[t<=star['T_DES_%s'%sp]] 
                        for t,sp in zip(temps,star.getDustList())]
                keytags = list(star.getDustList())
            else:
                include_total = 1
                print 'Thermal contact is on. All dust species share the ' + \
                      'same temperature profile. Vertical lines indicate ' + \
                      'inner radii of dust species.'
                rads, temps, keytags = [], [], []
            
            if include_total:
                rad = star.getDustRad(unit=unit)
                temp, key = star.getDustTemperature(add_key=1)
                rads.append(rad[rad>star['R_INNER_DUST']\
                                 *star.Rsun*star['R_STAR']])
                temps.append(temp[rad>star['R_INNER_DUST']\
                                  *star.Rsun*star['R_STAR']])
                keytags.append(key)
        
            #-- Add power laws if requested
            for s in power:
                rad_rstar = star_grid[0].getDustRad(unit='rstar')
                rad  = star_grid[0].getDustRad(unit=unit)
                tstar = star_grid[0]['T_STAR']
                temp,key = Profiler.dustTemperaturePowerLaw(rad=rad_rstar,\
                                                            add_key=1,\
                                                            tstar=tstar,s=s)
                rads.append(rad)
                temps.append(temp)
                keytags.append(key)
            
            ppars = dict()
            ppars['yaxis'] = '$T_\mathrm{d}$ (K)'
            if unit == 'rstar': ppars['xaxis'] = '$r\ \mathrm{(R}_\star\mathrm{)}$'
            elif unit =='au': ppars['xaxis'] = '$r\ \mathrm{(AU)}$'
            elif unit == 'm': ppars['xaxis'] = '$r\ \mathrm{(m)}$'
            else: ppars['xaxis'] = '$r\ \mathrm{(cm)}$'
            
            #-- Set plot filename
            pfn = fn_plt if fn_plt else 'Td_species'
            suff = star['LAST_MCMAX_MODEL']
            pfn = self.setFnPlt(pfn,fn_suffix=suff)

            plot_filenames.append(Plotting2.plotCols(x=rads,y=temps,\
                                  cfg=cfg_dict,filename=pfn,\
                                  keytags=keytags,xlogscale=1,ylogscale=1,\
                                  fontsize_key=16,**ppars))
        
        if len(plot_filenames) != len(star_grid):
            print 'At least one of the models does not yet have a MCMax model.'        
        if plot_filenames[0][-4:] == '.pdf':
            pfn = fn_plt if fn_plt else 'Td_species'
            pfn = self.setFnPlt(pfn) + '.pdf'
            DataIO.joinPdf(old=plot_filenames,new=pfn)
            print '** Your plots can be found at:'
            print pfn
            print '***********************************'
        else:
            print '** Plots can be found at:'
            print '\n'.join(plot_filenames)
            print '***********************************'
            


    def plotOpacities(self,star_grid=[],fn_plt='',scaling=0,species=['AMC'],\
                      cfg='',index=0,*args,**kwargs):
        
        """ 
        Plotting wavelength dependent mass extinction coefficients 
        (ie opacities).
        
        If based on star_grid or modelslist, they are scaled with abundances if 
        wanted.
        
        If no model info is given, the input opacities are plotted.
        
        Args and kwargs can be given straight to the plot command.
    
        @keyword star_grid: The input Star() models. If default, the MCMax 
                            input opacities are plotted.
                                  
                            (default: [])
        @type star_grid: list(Star())
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword scaling: allow species abundance scaling of opacities
                                
                          (default: 0)
        @type scaling: bool
        @keyword species: If no star_grid or model list are given, this gives 
                          the species requested to be plotted from Dust.dat
                            
                          (default: ['AMC'])
        @type species: list(string)
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        @keyword index: The index of the kappas in the .opacity/.particle file. 
                        0: extinction, 1: absorption, 2: scattering
                        
                        (default: 0)
        @type index: int
                
        """
        
        print '***********************************'
        print '** Starting to plot dust opacities.'
        
        #-- Set the filename
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
        
        if not star_grid and not fn_plt:
            fn_plt = os.path.join(cc.path.mopac,\
                                  'dust_opacities_%s'%'_'.join(species))

        #-- Set some plot parameters
        ppars = dict()
        ppars['xaxis'] = '$\lambda$ ($\mu \mathrm{m}$)'
        ppars['yaxis'] = '$\kappa_\lambda$ ($\mathrm{cm}^2\mathrm{/g}$)'
        ppars['xlogscale'] = 1
        ppars['ylogscale'] = 1
        ppars['key_location'] = (0.05,0.05)
        ppars.update(kwargs)
        ppars.update(cfg_dict)
        
        #-- Check if raw opacities or modeling results are requested
        if not star_grid:
            kr = KappaReader.KappaReader()
            wl_list = [kr.getKappas(sp)[0] for sp in species]
            q_list = [kr.getKappas(sp)[1] for sp in species]
            if not ppars.has_key('keytags'):
                ppars['keytags'] = species
            
            #-- Set plot filename
            pfn = fn_plt if fn_plt else 'opacities_species'
            pfn = self.setFnPlt(pfn)
            pfn = Plotting2.plotCols(x=wl_list,y=q_list,filename=pfn,\
                                     *args,**ppars)
            print '** Your plot can be found at:'
            print pfn
        else:    
            fns = []
            for star in star_grid:        
                try:    
                    wave,opacities = star.readKappas()
                except IOError:
                    continue
                opacities = [(opacities[i]+opacities[i+len(star.getDustList())]) 
                             for i in range(len(star.getDustList()))]
                if scaling:
                    opacities = [op*star['A_%s'%sp]
                                 for op,sp in zip(opacities,star.getDustList())]
                
                #-- Set plot filename
                pfn = fn_plt if fn_plt else 'opacities_species' 
                suff = star['LAST_MCMAX_MODEL']
                pfn = self.setFnPlt(pfn,fn_suffix=suff)
                
                keys = ['%s with $A$ = %s and $T_{des} = %i$ K'\
                         %(sp,str(star['A_%s'%sp]),int(star['T_DES_%s'%sp])) 
                        for sp in star.getDustList()]
                fns.append(Plotting2.plotCols(x=wave,y=opacities,keytags=keys,\
                                              filename=pfn,*args,**ppars))
            if len(fns) != len(star_grid):
                print 'At least one of the models requested does not yet ' + \
                      'have a MCMax model.'
            print '** Your plots can be found at:'
            if fns and fns[-1][-4] == '.pdf':
                pfn = fn_plt if fn_plt else 'opacities_species'
                pfn = self.setFnPlt(pfn) + '.pdf'
                DataIO.joinPdf(old=fns,new=pfn)
                print pfn
            else:
                print '\n'.join(fns)
        print '***********************************'
        


    def plotExtinction(self,star_grid=[],models=[],fn_plt='',plot_default=1,\
                       cfg=''):
        
        """ 
        Plotting wavelength dependent extinction efficiencies wrt grain size.
        
        This always depends on a star_grid or one created from a list of MCMax 
        model ids.
        
        Plotted are the total efficiencies, including relative weights between 
        the included dust species. This is the input for GASTRoNOoM!
        
        @keyword star_grid: List of Star() instances. If default, model ids 
                            have to be given.
                                  
                            (default: [])
        @type star_grid: list[Star()]
        @keyword models: The model ids, only required if star_grid is []
        
                         (default: [])
        @type models: list[string]
        @keyword fn_plt: A base plot filename. Includes folder. If not, a 
                         default is added
                         
                         (default: '')
        @type fn_plt: string  
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
        
        """
        
        print '***********************************'
        print '** Plotting Q_ext/a.'
        if not star_grid and not models:
            print 'Input is undefined. Aborting.'
            return      
        elif not star_grid and models:
            star_grid = self.makeMCMaxStars(models=models)
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict.pop('filename')
            
        x = []
        y = []
        keys = []
        for star in star_grid:        
            try:
                inputfile = os.path.join(cc.path.gdata,star['TEMDUST_FILENAME'])
                opacities = DataIO.readCols(filename=inputfile)
                x.append(opacities[0])
                y.append(opacities[1])
                keys.append('$Q_\mathrm{ext}/a$ for MCMax %s'\
                            %star['LAST_MCMAX_MODEL'].replace('_','\_'))
            except IOError: 
                pass
        
        #-- Set plot filename
        pfn = fn_plt if fn_plt else 'gastronoom_opacities'
        pfn = self.setFnPlt(pfn)

        title = 'GASTRoNOoM Extinction Efficiencies in %s'\
                 %(self.star_name_plots)
        pfn = Plotting2.plotCols(x=x,y=y,cfg=cfg_dict,filename=pfn,\
                                 xaxis='$\lambda$ ($\mu$m)',keytags=keys,\
                                 yaxis='$Q_{ext}/a$ (cm$^{-1}$)',\
                                 plot_title=title,key_location=(0.7,0.6),\
                                 xlogscale=1,ylogscale=1,fontsize_key=20)
        print '** The extinction efficiency plot can be found at:'
        print pfn
        print '***********************************'  
            
            
            
    def makeMCMaxStars(self,models):
        
        '''
        Set parameters for star_list taken from the MCMax database.
        
        Based on the model id of MCMax.
        
        @param models: model_ids for the MCMax db
        @type models: list(string)
        @return: The model instances 
        @rtype: list(Star())
        
        '''
        
        star_grid = Star.makeStars(models=models,code='MCMax',id_type='MCMax',\
                                   path=self.path)
        for star,model in zip(star_grid,models):    
            filepath = os.path.join(cc.path.mout,'models',\
                                    star['LAST_MCMAX_MODEL'])
            denstemp = os.path.join(filepath,'denstemp.dat')
            logfile = os.path.join(filepath,'log.dat')
            grid_shape = DataIO.getMCMaxOutput(filename=denstemp,incr=1,\
                                               keyword='NGRAINS',single=0)[0]
            star.update({'NTHETA':int(grid_shape[1]),\
                         'NRAD':int(grid_shape[0]),\
                         'T_STAR':float(DataIO.getMCMaxOutput(filename=logfile,\
                                                incr=0,\
                                                keyword='STELLAR TEMPERATURE',\
                                                single=0)[0][2]),\
                         'R_STAR':float(DataIO.getMCMaxOutput(filename=logfile,\
                                                incr=0,\
                                                keyword='STELLAR RADIUS',\
                                                single=0)[0][2])})            
        return star_grid  
        
'''
#Some titles, units etc strings for title and label representation
        
pars_units = dict([('T_STAR',('T_{*}','K','%i')),\
                                ('L_STAR',('L_{*}','L_{@{&{i}_{/=12 \307}}O}','%i')),\
                                ('R_STAR',('R_*','R_{@{&{i}_{/=12 \307}}O}','%.1f')),\
                                ('M_STAR',('M_{*}','M_{@{&{i}_{/=12 \307}}O}','%.1f')),\
                                ('DISTANCE',('d','pc','%.1f')),\
                                ('MDOT_DUST',('@^{&{/=8 o}{/=13 \264}}M_{d}','M_{@{&{i}_{/=11 \307}}O}/yr','%.2e')),\
                                ('MDOT_GAS',('@^{&{/=8 o}{/=13 \264}}M_{g}','M_{@{&{i}_{/=11 \307}}O}/yr','%.2e')),\
                                ('R_INNER_DUST',('R_{inner,d}','R_{*}','%.2f')),\
                                ('R_INNER_GAS',('R_{inner,g}','R_{*}','%.2f')),\
                                ('R_OUTER_DUST',('R_{outer,d}','R_{*}','%i')),\
                                ('R_OUTER_GAS',('R_{outer,g}','R_{*}','%i')),\
                                ('T_INNER_DUST',('T_{inner,d}','K','%i')),\
                                ('DUST_TO_GAS',('{/Symbol Y}_{semi-emp}','','%.4f')),\
                                ('V_EXP_DUST',('v_{{/Symbol \245},d}','km/s','%.2f')),\
                                ('MDUST',('M_{d}','M_{@{&{i}_{/=12 \307}}O}','%.2e')),\
                                ('VEL_INFINITY_GAS',('v_{{/Symbol \245},g}','km/s','%.2f')),\
                                ('DRIFT',('w_{/Symbol \245}','km/s','%.2f')),\
                                ('SED',('Ray-traced SED','','%i')),\
                                ('T_CONTACT',('Thermal contact','','%i')),\
                                ('PHOTON_COUNT',('Photon count','','%i')),\
                                ('LAM1',('{/Symbol l} grid (min)','{/Symbol m}m','%.1f')),\
                                ('LAM2',('{/Symbol l} grid (max)','{/Symbol m}m','%.1f')),\
                                ('NLAM',('{/Symbol l} grid','Steps','%i')),\
                                ('ZLAM1',('{/Symbol l} subgrid (min)','{/Symbol m}m','%.2f')),\
                                ('ZLAM2',('{/Symbol l} subgrid (max)','{/Symbol m}m','%.2f')),\
                                ('NZLAM',('{/Symbol l} subgrid','Steps','%i')),\
                                ('USE_MARCS',('MARCS','','%i')),\
                                ('MARCS_TYPE',('MARCS type','','%s')),\
                                ('MARCS_KERNEL',('MARCS kernel','','%i')),\
                                ('DENSTYPE',('Density type','','%s')),\
                                ('NRAD',('R grid','Steps','%s')),\
                                ('NTHETA',('{/Symbol q} grid','Steps','%s')),\
                                ('DUST_TEMPERATURE_FILENAME',('T_{d}(r)','','%s')),\
                                ('OPA_FILE',('MCMax qpr','','%s')),\
                                ('TEMDUST_FILE',('MCMax temdust.k','','%s')),\
                                ('DENSFILE',('{/Symbol r}_{d}(r)','','%s')),\
                                ('TDESITER',('T_{des,iter}','','%i')),\
                                ('FLD',('FLD','','%i')),\
                                ('N_QUAD',('N_{quad}','','%s')),\
                                ('RATIO_12C_TO_13C',('^{12}C/^{13}C','','%i')),\
                                ('OPR',('o-H_2O/p-H_2O','','%i')),\
                                ('KEYWORD_DUST_TEMPERATURE_TABLE',('Consistent T_{d}','','%i')),\
                                ('NUMBER_INPUT_DUST_TEMP_VALUES',('len(T_d)','','%i')),\
                                ('MOLECULE',('Molecule','','%s'))])
        dust_species = DataIO.getInputData(keyword='SPECIES_SHORT',filename='Dust.dat')
        pars_units.update(dict([('A_' + species,('A_{' + species + '}','','%.2f')) for species in dust_species]))
        pars_units.update(dict([('T_DESA_' + species,('T_{desA,' + species + '}','','%.3f')) for species in dust_species]))
        pars_units.update(dict([('T_DESB_' + species,('T_{desB,' + species + '}','','%.3f')) for species in dust_species]))
        pars_units.update(dict([('T_DES_' + species,('T_{des,' + species + '}','','%.3f')) for species in dust_species]))
        pars_units.update(dict([('T_MIN_' + species,('T_{min,' + species + '}','K','%i')) for species in dust_species]))
        pars_units.update(dict([('T_MAX_' + species,('T_{max,' + species + '}','K','%i')) for species in dust_species]))
        pars_units.update(dict([('R_MIN_' + species,('R_{min,' + species + '}','R_*','%.2f')) for species in dust_species]))
        pars_units.update(dict([('R_MAX_' + species,('R_{max,' + species + '}','R_*','%.2f')) for species in dust_species]))
        '''