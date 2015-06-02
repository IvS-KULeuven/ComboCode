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

from cc.plotting.PlottingSession import PlottingSession
from cc.plotting import Plotting2
from cc.tools.io import DataIO
from cc.modeling.objects import Star
from cc.modeling.codes import MCMax
from cc.modeling.tools import Profiler



class PlotDust(PlottingSession):
    
    """ 
    Plotting environment for SEDs and all dust parameters.
    
    """
    
    def __init__(self,star_name='model',sed=None,corrflux_path=None,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),\
                 path_mcmax='runTestDec09',inputfilename=None):
        
        '''
        Initializing PlotDust session.
        
        @keyword star_name: name of the star from Star.dat, use default only 
                            when never using any star model specific things 
                                  
                            (default: "model")
        @type star_name: string
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        @keyword path_mcmax: Output modeling folder in MCMax home folder
        
                             (default: 'runTestDec09')
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
        @keyword corrflux_path: The full path to the folder containing 
                                correlated fluxes, such as for MIDI. 
                                
                                (default: None)
        @type corrflux_path: str
        
        '''

        super(PlotDust, self).__init__(star_name=star_name,\
                                       path_combocode=path_combocode,\
                                       path=path_mcmax,\
                                       code='MCMax',\
                                       inputfilename=inputfilename)
        self.sed = sed
        self.corrflux_path = corrflux_path



    def plotSed(self,star_grid=[],cfg='',iterative=0,no_models=0,\
                fn_add_star=0):
        
        """ 
        Creating an SED with 0, 1 or more models and data. 
        
        Includes data preparation on the spot.
        
        @keyword star_grid: list of Star() models to plot. If star_grid is [], 
                            only data are plotted.
                            
                            (default: [])
        @type star_grid: list[Star()]
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
        @keyword fn_add_star: Add the star name to the requested plot filename.
        
                              (default: 1)
        @type fn_add_star: bool
        
        """
        
        if self.sed is None:
            print 'No PATH_SED given. Cannot plot SED. Aborting...'
            return
        print '***********************************'
        print '** Creating SED plot.'
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('no_models'):
            no_models = cfg_dict['no_models']
        if cfg_dict.has_key('fn_add_star'):
            fn_add_star = bool(cfg_dict['fn_add_star'])
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        else:
            fn_plt = ''
        
        ccpath = os.path.join(self.path_combocode,'Data')
        data_labels = dict([(dt,(n,ls))
                            for n,dt,ls in zip(DataIO.getInputData(path=ccpath,\
                                                        keyword='PLOT_NAMES',\
                                                        filename='Sed.dat',\
                                                        remove_underscore=1),\
                                               DataIO.getInputData(path=ccpath,\
                                                        keyword='DATA_TYPES',\
                                                        filename='Sed.dat'),\
                                               DataIO.getInputData(path=ccpath,\
                                                        keyword='LINE_TYPES',\
                                                        filename='Sed.dat'))])

        #- filename settings and copying inputfiles to plot output folder
        if not fn_plt:
            fn_plt = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                  'stars',self.star_name,self.plot_id,\
                                  'SED_%s'%self.star_name)
        if fn_add_star:
            fn_plt = '_'.join([fn_plt,self.star_name])
        if iterative:
            fn_plt = fn_plt + '_iterative_%i'%iterative
        
        if self.inputfilename <> None:
            subprocess.call(['cp ' + self.inputfilename + ' ' + \
                             os.path.join(os.path.expanduser('~'),'MCMax',\
                                self.path,'stars',self.star_name,self.plot_id,\
                                os.path.split(self.inputfilename)[1])],\
                            shell=True)

        plot_title='SED %s'%self.star_name_plots
       
        #- prepare and collect data, keytags and line types
        keytags = []
        data_x = []
        data_y = []
        line_types = []
        for (dt,fn),(w,f) in sorted([dset
                                     for dset in self.sed.data.items()
                                     if 'PHOT' not in dset[0][0].upper()]):
             keytags.append(data_labels[dt][0])
             data_x.append(self.sed.data[(dt,fn)][0])
             data_y.append(self.sed.data[(dt,fn)][1])
             line_types.append(data_labels[dt][1])
        
        for (dt,fn),(w,f) in sorted([dset 
                                     for dset in self.sed.data.items()
                                     if 'PHOT' in dset[0][0].upper()]):
             keytags.append(data_labels[dt][0])
             data_x.append(self.sed.data[(dt,fn)][0])
             data_y.append(self.sed.data[(dt,fn)][1])
             line_types.append(data_labels[dt][1])
        
        #- Collect model data as well as keytags and set line types
        model_ids_mcm = [s['LAST_MCMAX_MODEL'] 
                         for s in star_grid
                         if s['LAST_MCMAX_MODEL']]
        #- Only if the model_ids list is not empty, MCMax models are available
        #- Otherwise the ray tracing keyword is unnecessary.
        if no_models:
            model_ids_mcm = []
        if model_ids_mcm: 
            rt_sed = star_grid[0]['RT_SED']
        for model_id in model_ids_mcm:
            dpath = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                 'models',model_id)
            w,f = MCMax.readModelSpectrum(dpath,rt_sed)
            data_x.append(w)
            data_y.append(f)
            keytags.append(model_id.replace('_','\_'))

        line_types += [0]*len(star_grid)
        keytags = [tag.replace('#','') for tag in keytags]
        extra_pars = dict()
        try:
            extra_pars['ymax'] = 1.3*max([max(dy) for dy in data_y])
        except ValueError:
            pass
        try:    
            extra_pars['ymin'] = 0.5*min([min(dy) for dy in data_y])
        except ValueError:
            pass        
        filename = Plotting2.plotCols(x=data_x,y=data_y,filename=fn_plt,\
                                      figsize=(20,10),number_subplots=1,\
                                      plot_title=plot_title,fontsize_axis=20,\
                                      keytags=keytags,fontsize_title=24,\
                                      linewidth=3,key_location=(0.0,0.75),\
                                      xlogscale=1,transparent=0,cfg=cfg_dict,\
                                      line_types=line_types,ylogscale=0,\
                                      fontsize_ticklabels=20,fontsize_key=18,\
                                      xmin=2,xmax=200,extension='.pdf',\
                                      **extra_pars)
        print '** Your SED plots can be found at:'
        print filename
        print '***********************************'
         
         
         
    def plotCorrflux(self,star_grid=[],cfg='',no_models=0,fn_add_star=0):
        
        """ 
        Plot correlated fluxes with 0, 1 or more models and data. 
        
        Includes data preparation on the spot.
        
        @keyword star_grid: list of Star() models to plot. If star_grid is [], 
                            only data are plotted.
                            
                            (default: [])
        @type star_grid: list[Star()]
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
        
        if not self.corrflux_path:
            print 'No CORRFLUX_PATH given. Cannot plot Correlated Fluxes. Aborting...'
            return
        print '***********************************'
        print '** Creating Correlated Fluxes plot.'
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('no_models'):
            no_models = cfg_dict['no_models']
        if cfg_dict.has_key('fn_add_star'):
            fn_add_star = bool(cfg_dict['fn_add_star'])
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        else:
            fn_plt = ''
        
        #- filename settings and copying inputfiles to plot output folder
        if not fn_plt:
            fn_plt = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                  'stars',self.star_name,self.plot_id,\
                                  'CorrFlux')
        if fn_add_star:
            fn_plt = '_'.join([fn_plt,self.star_name])
        
        if self.inputfilename <> None:
            subprocess.call(['cp ' + self.inputfilename + ' ' + \
                             os.path.join(os.path.expanduser('~'),'MCMax',\
                                self.path,'stars',self.star_name,self.plot_id,\
                                os.path.split(self.inputfilename)[1])],\
                            shell=True)

        #-- Currently only MIDI is implemented. 
        #   Assumption: inc=01 is baseline 46.5m   
        #               inc=02 is baseline 51.4m
        #               inc=03 is baseline 60.6m
        baseline = dict([('46.5m','visibility01.0.dat'),\
                         ('54.4m','visibility02.0.dat'),\
                         ('60.6m','visibility03.0.dat')])

        ssd = os.path.join(self.corrflux_path,self.star_name,\
                           '_'.join([self.star_name,'MIDI','*.fits']))
        ggd = dict([(gi[-10:-5],gi) for gi in glob.glob(ssd)
                    if gi[-10:-5] in ('46.5m','54.4m','60.6m')])
        
        #-- prepare and collect data, keytags and line types
        data = []
        for k,v in sorted(ggd.items()):
            ddict = dict()
            data.append(ddict)
            
            #-- Extract data from the fits file
            dfits = pyfits.open(v)
            x = 1e6*dfits['OI_WAVELENGTH'].data['EFF_WAVE']
            ddict['x'] = [x]
            ddict['y'] = [dfits['OI_VIS'].data['VISAMP'][0]]
            ddict['yerr'] = [dfits['OI_VIS'].data['VISAMPERR'][0]]
            dfits.close()
            #-- Wavelength limits between 8 and 13 micron, limits of the N band
            #   atmospheric transmission. Outside these ranges, the flux is not
            #   relevant
            ddict['xmin'] = 8
            ddict['xmax'] = 13
            ddict['ymin'] = -0.1
            ddict['ymax'] = 1.3*max(ddict['y'][0][(x<13)*(x>8)])
            ddict['labels'] = [('MIDI %s'%k,0.05,0.9)]
            
            if no_models:
                continue
            #-- Extract models from the model folders
            for s in star_grid:
                model_id = s['LAST_MCMAX_MODEL']
                dpath = os.path.join(os.path.expanduser('~'),'MCMax',\
                                     self.path,'models',model_id)
                model = MCMax.readVisibilities(dpath=dpath,fn_vis=baseline[k])
                ddict['x'].append(model[0])
                if model[1] != []:
                    ddict['y'].append(model[1]*model[2])
                else:
                    ddict['y'].append([])
                ddict['yerr'].append(None)
            
        kwargs = dict()
        kwargs['keytags'] = ['MIDI']
        if not no_models:
            kwargs['keytags'].extend([s['LAST_MCMAX_MODEL'].replace('_','\_') 
                                      for s in star_grid])
        kwargs['xaxis'] = '$\lambda$ ($\mu$m)'
        kwargs['yaxis'] = 'Corr.~FLux (Jy)'
        kwargs['dimensions'] = (1,4)
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
        filename = Plotting2.plotTiles(data=data,filename=fn_plt,**kwargs)
        print '** Your Correlated Flux plots can be found at:'
        print filename
        print '***********************************'
                 
                 
                 
    def plotCorrfluxBaseline(self,wav,star_grid=[],cfg='',no_models=0,\
                             fn_add_star=0):
        
        """ 
        Plot correlated fluxes as a function of baseline.
        
        Requires wavelengths to be requested.
        
        Includes data preparation on the spot.
        
        Data location is that of correlated flux, but not required. Models can
        be plotted without data.

        @param wav: The wavelengths (in micron) at which the visibilities
                    are plotted. 
        @type wav: list[float]
        
        @keyword star_grid: list of Star() models to plot. If star_grid is [], 
                            only data are plotted.
                            
                            (default: [])
        @type star_grid: list[Star()]
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
        
        print '***********************************'
        print '** Creating Visibilities plot.'
        
        cfg_dict = Plotting2.readCfg(cfg)
        if cfg_dict.has_key('no_models'):
            no_models = cfg_dict['no_models']
        if cfg_dict.has_key('fn_add_star'):
            fn_add_star = bool(cfg_dict['fn_add_star'])
        if cfg_dict.has_key('wav'):
            wav = cfg_dict['wav']        
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']
        else:
            fn_plt = ''
        
        #- filename settings and copying inputfiles to plot output folder
        if not fn_plt:
            fn_plt = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                  'stars',self.star_name,self.plot_id,\
                                  'CorrFluxBaseline')
        if fn_add_star:
            fn_plt = '_'.join([fn_plt,self.star_name])
        
        if self.inputfilename <> None:
            subprocess.call(['cp ' + self.inputfilename + ' ' + \
                             os.path.join(os.path.expanduser('~'),'MCMax',\
                                self.path,'stars',self.star_name,self.plot_id,\
                                os.path.split(self.inputfilename)[1])],\
                            shell=True)

        #-- Currently only MIDI is implemented. 
        #   Assumption: inc=01 is baseline 46.5m   
        #               inc=02 is baseline 51.4m
        #               inc=03 is baseline 60.6m
        baseline = dict([('46.5m','visibility01.0.dat'),\
                         ('54.4m','visibility02.0.dat'),\
                         ('60.6m','visibility03.0.dat')])

        ssd = os.path.join(self.corrflux_path,self.star_name,\
                           '_'.join([self.star_name,'MIDI','*.fits']))
        ggd = dict([(gi[-10:-5],gi) for gi in glob.glob(ssd)
                    if gi[-10:-5] in ('46.5m','54.4m','60.6m')])
        
        #-- Collect MIDI data from the fits file
        ddf = dict()
        for k,v in sorted(ggd.items()):
            dfits = pyfits.open(v)
            ddf[k] = dict()
            ddf['x'] = 1e6*dfits['OI_WAVELENGTH'].data['EFF_WAVE']
            ddf['y'] = dfits['OI_VIS'].data['VISAMP'][0]
            ddf['yerr'] = dfits['OI_VIS'].data['VISAMPERR'][0]
            dfits.close()
        
        #-- prepare and collect plot data, keytags and line types
        data = []
        for w in wav:
            ddict = dict()
            data.append(ddict)
            
            #-- Set the plot x and y
            bls = [k for k in sorted(ddf.keys())]
            ddict['x'] = [float(bl.strip('m')) for bl in bls]
            ddict['y'] = [ddf[bl]['y'][argmin(abs(ddf[bl]['x']-w))] 
                          for bl in bls]
            ddict['yerr'] = [ddf[bl]['yerr'][argmin(abs(ddf[bl]['x']-w))] 
                             for bl in bls]
            
            ddict['ymin'] = -0.1
            ddict['ymax'] = 1.3*max(ddict['y'][0][(x<13)*(x>8)])
            ddict['labels'] = [('MIDI %s'%k,0.05,0.9)]
            
            if no_models:
                continue
            
            #-- Extract models from the model folders
            for s in star_grid:
                model_id = s['LAST_MCMAX_MODEL']
                dpath = os.path.join(os.path.expanduser('~'),'MCMax',\
                                     self.path,'models',model_id)
                model = MCMax.readVisibilities(dpath=dpath,fn_vis=baseline[k])
                ddict['x'].append(model[0])
                if model[1] != []:
                    ddict['y'].append(model[1]*model[2])
                else:
                    ddict['y'].append([])
                ddict['yerr'].append(None)
            
        kwargs = dict()
        kwargs['keytags'] = ['MIDI']
        if not no_models:
            kwargs['keytags'].extend([s['LAST_MCMAX_MODEL'].replace('_','\_') 
                                      for s in star_grid])
        kwargs['xaxis'] = '$\lambda$ ($\mu$m)'
        kwargs['yaxis'] = 'Corr.~FLux (Jy)'
        kwargs['dimensions'] = (1,4)
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
        filename = Plotting2.plotTiles(data=data,filename=fn_plt,**kwargs)
        print '** Your Correlated Flux plots can be found at:'
        print filename
        print '***********************************'

                                   
                                    
    def plotTemp(self,star_grid=[],models=[],power=[1],fn_plt='',cfg=''):
        
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
        @keyword power: A list of values for s in below formula. If [] no power
                        law is included. Power law parameters  are taken from 
                        star_grid[0].
                                
                        See Thesis p32, where power is s in 
                        T(r) = T_eff*(2*r/R_STAR)**(-2/(4+s)).
                
                        (default: [1])
        @type power: list        
        @keyword fn_plt: A plot filename for the tiled plot.
                         
                         (default: '')
        @type fn_plt: string
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
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']    
        else:
            fn_plt = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                  'stars',self.star_name,self.plot_id,\
                                  'Td_avg')
        radii = []
        temps = []
        keytags = []
        for star in star_grid:
            rad = star.getDustRad()
            temp,key = star.getDustTemperature(add_key=1) 
            radii.append(rad)
            temps.append(temp)
            keytags.append(key)
        
        #-- Add power laws if requested
        for s in power:
            rad = star_grid[0].getDustRad(unit='rstar')
            tstar = star_grid[0]['T_STAR']
            temp,key = Profiler.dustTemperaturePowerLaw(rad=rad,add_key=1,\
                                                        tstar=tstar,s=s)
            radii.append(rad)
            temps.append(temp)
            keytags.append(key)
            
        title = 'Average Dust Temperature Stratification for %s'\
                %(self.star_name_plots)
        filename = Plotting2.plotCols(x=radii,y=temps,filename=fn_plt,\
                                      yaxis='$T_\mathrm{d}$ (K)',\
                                      plot_title=title,xaxis='$R$ (cm)',\
                                      key_location=(0.05,0.05),cfg=cfg_dict,\
                                      xlogscale=1,ylogscale=1,fontsize_key=20,\
                                      keytags=keytags,fontsize_axis=26,\
                                      figsize=(12.5,8),linewidth=3,\
                                      fontsize_ticklabels=26,)
        print '** Your plots can be found at:'
        print filename
        print '***********************************'
            


    def plotTempSpecies(self,star_grid=[],models=[],include_total=1,\
                        power=[1],fn_plt='',cfg=''):
        
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
        @keyword include_total: Include the sum of all temperature profiles as 
                                well for comparison. 
                                        
                                (default: 0)
        @type include_total: bool
        @keyword power: A list of values for s in below formula. If [] no power
                        law is included. Power law parameters  are taken from 
                        star_grid[0].
                                
                        See Thesis p32, where power is s in 
                        T(r) = T_eff*(2*r/R_STAR)**(-2/(4+s)).
                
                        (default: [1])
        @type power: list        
        @keyword fn_plt: A plot filename for the tiled plot.
                         
                         (default: '')
        @type fn_plt: string
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
        if cfg_dict.has_key('filename'):
            fn_plt = cfg_dict['filename']
            del cfg_dict['filename']    
        else:
            fn_plt = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                  'stars',self.star_name,self.plot_id,\
                                  'Td_species')
        plot_filenames = []
        for star in star_grid:
            if not int(star['T_CONTACT']):
                rads = [star.getDustRad(species=species)
                        for species in star.getDustList()]
                temps = [star.getDustTemperature(species=species)
                         for species in star.getDustList()]
                rads = [r[t<=star['T_DES_%s'%sp]] 
                        for r,t,sp in zip(rads,temps,star.getDustList())]
                temps = [t[t<=star['T_DES_%s'%sp]] 
                        for t,sp in zip(temps,star.getDustList())]
                keytags = star.getDustList()
                vert_lines = []
            else:
                include_total = 1
                print 'Thermal contact is on. All dust species share the ' + \
                      'same temperature profile. Vertical lines indicate ' + \
                      'inner radii of dust species.'
                radii, temps, keytags = [], [], []
                vert_lines = [star['R_DES_%s'%d]*star.Rsun*star['R_STAR'] 
                              for d in star.getDustList()]
            rad = star.getDustRad()
            if include_total:
                temp, key = star.getDustTemperature(add_key=1)
                radii.append(rad[rad>star['R_INNER_GAS']\
                                 *star.Rsun*star['R_STAR']])
                temps.append(temp[rad>star['R_INNER_GAS']\
                                  *star.Rsun*star['R_STAR']])
                keytags.append(key)
        
            #-- Add power laws if requested
            for s in power:
                rad = star_grid[0].getDustRad(unit='rstar')
                tstar = star_grid[0]['T_STAR']
                temp,key = Profiler.dustTemperaturePowerLaw(rad=rad,add_key=1,\
                                                            tstar=tstar,s=s)
                radii.append(rad)
                temps.append(temp)
                keytags.append(key)
                
            filename = '_'.join([fn_plt,star['LAST_MCMAX_MODEL']])
            plot_filenames.append(Plotting2.plotCols(x=radii,y=temps,\
                        cfg=cfg_dict,filename=filename,xaxis='$r$ (cm)',\
                        yaxis='$T_\mathrm{d}$ (K)',keytags=keytags,\
                        xmax=star['R_OUTER_DUST']*star.Rsun*star['R_STAR'],\
                        xmin=star['R_STAR']*star.Rsun,fontsize_axis=26,\
                        xlogscale=1,ylogscale=1,fontsize_key=16,\
                        figsize=(12.5,8),transparent=0,linewidth=3,\
                        fontsize_ticklabels=26,\
                        vert_lines=[star['R_INNER_DUST']\
                                        *star.Rsun*star['R_STAR']]))
        if len(plot_filenames) != len(star_grid):
            print 'At least one of the models does not yet have a MCMax model.'        
        if plot_filenames[0][-4:] == '.pdf':
            new_filename = fn_plt + '.pdf'
            DataIO.joinPdf(old=plot_filenames,new=new_filename)
            print '** Your plots can be found at:'
            print new_filename
            print '***********************************'
        else:
            print '** Plots can be found at:'
            print '\n'.join(plot_filenames)
            print '***********************************'
            


    def plotOpacities(self,star_grid=[],models=[],scaling=1,species=['AMC'],\
                      cfg=''):
        
        """ 
        Plotting wavelength dependent mass extinction coefficients 
        (ie opacities).
        
        If based on star_grid or modelslist, they are scaled with abundances if 
        wanted.
        
        If no model info is given, the input opacities are plotted.
    
        @keyword star_grid: The input Star() models. If default, the MCMax 
                            input opacities are plotted.
                                  
                            (default: [])
        @type star_grid: list(Star())
        @keyword models: MCMax model_ids. Can be given instead of star_grid. 
                         Not yet implemented!
                              
                         (default: [])
        @type models: list(string)
        @keyword scaling: allow species abundance scaling of opacities
                                
                          (default: 1)
        @type scaling: bool
        @keyword species: If no star_grid or model list are given, this gives 
                          the species requested to be plotted from Dust.dat
                            
                          (default: ['AMC'])
        @type species: list(string)
        @keyword cfg: path to the Plotting2.plotCols config file. If default, 
                      the hard-coded default plotting options are used.
                          
                      (default: '')
        @type cfg: string
                
        """
        
        print '***********************************'
        print '** Starting to plot dust opacities.'
        if not star_grid and models:
            #star_grid = self.makeMCMaxStars(models=models,id_type='MCMax')
            raise IOError('Reading dust opacities from a model id list ' + \
                          'only, not yet implemented.')
        filenames = []
        if not star_grid and not models:
            species_index = [DataIO.getInputData(keyword='SPECIES_SHORT',\
                                    filename='Dust.dat',\
                                    path=os.path.join(self.path_combocode,\
                                                      'Data')).index(sp)
                             for sp in species]
            filename_species = [DataIO.getInputData(keyword='PART_FILE',\
                                    filename='Dust.dat',\
                                    path=os.path.join(self.path_combocode,\
                                                      'Data'))[i] 
                                for i in species_index]
            part_file = [DataIO.readFile(filename=\
                                os.path.join(os.path.expanduser('~'),'MCMax',\
                                             'src',fn),\
                                         delimiter=' ') 
                         for fn in filename_species]
            wl_list = [array([float(wl[0]) 
                              for wl in pf if len(wl) == 4]) 
                       for pf in part_file]
            q_list = [array([float(q[1]) 
                             for q in pf if len(q) == 4]) 
                      for pf in part_file]
            filename = os.path.join(os.path.expanduser('~'),'MCMax',\
                                    'dust_opacities_%s'%'_'.join(species))
            filename = Plotting2.plotCols(x=wl_list,y=q_list,\
                                          filename=filename,\
                                          xaxis='$\lambda$ ($\mu \mathrm{m}$)',cfg=cfg,\
                                          yaxis='$\kappa_\lambda$ ($\mathrm{cm}^2\mathrm{/g}$)',\
                                          keytags=species,fontsize_key=20,\
                                          plot_title='Dust Opacities',\
                                          key_location=(0.05,0.05),\
                                          number_subplots=1,xlogscale=1,\
                                          ylogscale=1)
            print '** Your plot can be found at:'
            print filename
        else:    
            for star in star_grid:        
                try:    
                    wave,opacities = star.readKappas()
                except IOError:
                    continue
                opacities = [(opacities[i]+opacities[i+len(star.getDustList())]) 
                             for i,species in enumerate(star.getDustList())]
                if scaling:
                    opacities = [opa*float(star['A_%s'%species]) 
                                 for opa in opacities]
                filename = os.path.join(os.path.expanduser('~'),'MCMax',\
                                        self.path,'stars',self.star_name,\
                                        self.plot_id,\
                                        'opacities_species_%s'\
                                        %star['LAST_MCMAX_MODEL'])
                title = 'Dust Opacities in %s (%s)' \
                        %(self.star_name_plots,\
                            star['LAST_MCMAX_MODEL'].replace('_','\_'))
                keytags = ['%s with $A$ = %s and $T_{des} = %i$ K'\
                           %(sp,str(star['A_%s'%sp]),int(star['T_DES_%s'%sp])) 
                           for sp in star.getDustList()]
                filenames.append(Plotting2.plotCols(x=wave,y=opacities,\
                                 xaxis='$\lambda$ ($\mu \mathrm{m}$)',\
                                 yaxis='$\kappa_\lambda$ ($\mathrm{cm}^2\mathrm{/g}$)',\
                                 keytags=keytags,plot_title=title,\
                                 key_location=(0.05,0.05),filename=filename,\
                                 cfg=cfg,number_subplots=1,xlogscale=1,\
                                 ylogscale=1,fontsize_key=20))
            if len(filenames) != len(star_grid):
                print 'At least one of the models requested does not yet ' + \
                      'have a MCMax model.'
            print '** Your plots can be found at:'
            if filenames[-1][-4] == '.pdf':
                new_file = os.path.join(os.path.expanduser('~'),'MCMax',\
                                        self.path,'stars',self.star_name,\
                                        self.plot_id, 'dust_opacities.pdf')
                DataIO.joinPdf(old=filenames,new=new_file)
                print new_file
            else:
                print '\n'.join(filenames)
        print '***********************************'
        


    def plotExtinction(self,star_grid=[],models=[],plot_default=1,cfg=''):
        
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
        @keyword plot_default: Include the default extinction efficiencies for 
                               amorphous silicates (temdust.kappa)
                               [NYI]
                                      
                               (default: 1)
        @type plot_default: bool
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
        if cfg_dict.has_key('plot_default'):
            plot_default = int(cfg['plot_default'])
        x = []
        y = []
        keys = []
        for star in star_grid:        
            try:
                inputfile = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                         'src','data',star['TEMDUST_FILENAME'])
                opacities = DataIO.readCols(filename=inputfile)
                x.append(opacities[0])
                y.append(opacities[1])
                keys.append('$Q_\mathrm{ext}/a$ for MCMax %s'\
                            %star['LAST_MCMAX_MODEL'].replace('_','\_'))
            except IOError: 
                pass
        if plot_default:
             print 'Including default amorphous silicate extinction ' + \
                   'efficiencies not yet implemented.'
        filename = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                'stars',self.star_name,self.plot_id,\
                                'gastronoom_opacities_%s'\
                                %star['LAST_MCMAX_MODEL'])
        title = 'GASTRoNOoM Extinction Efficiencies in %s'\
                 %(self.star_name_plots)
        filename = Plotting2.plotCols(x=x,y=y,cfg=cfg,filename=filename,\
                                      xaxis='$\lambda$ ($\mu$m)',keytags=keys,\
                                      yaxis='$Q_{ext}/a$ (cm$^{-1}$)',\
                                      plot_title=title,key_location=(0.7,0.6),\
                                      xlogscale=1,ylogscale=1,fontsize_key=20)
        print '** The extinction efficiency plot can be found at:'
        print filename
        print '***********************************'  
            
            
            
    def makeMCMaxStars(self,models,\
                      data_path=os.path.join(os.path.expanduser('~'),'MCMax',\
                                             'Data')):
        
        '''
        Set parameters for star_list taken from the MCMax database.
        
        Based on the model id of MCMax.
        
        @param models: model_ids for the MCMax db
        @type models: list(string)
        @keyword data_path: path to the data used here
        
                            (default: ~/MCMax/Data)
        @type data_path: string
        @return: The model instances 
        @rtype: list(Star())
        
        '''
        
        star_grid = Star.makeStars(models=models,\
                                   code='MCMax',id_type='MCMax',path=self.path)
        for star,model in zip(star_grid,models):    
            filepath = os.path.join(os.path.expanduser('~'),'MCMax',\
                                    self.path,'models',\
                                    star['LAST_MCMAX_MODEL'])
            denstemp = os.path.join(filepath,'denstemp.dat')
            logfile = os.path.join(filepath,'log.dat')
            grid_shape = DataIO.getMCMaxOutput(filename=denstemp,incr=1,\
                                               keyword='NGRAINS',single=0)[0]
            star.update({'PATH_DUST_DATA':data_path,\
                         'NTHETA':int(grid_shape[1]),\
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
        dust_species = DataIO.getInputData(path=os.path.join(self.path_combocode,'Data'),keyword='SPECIES_SHORT',\
                                                    filename='Dust.dat')
        pars_units.update(dict([('A_' + species,('A_{' + species + '}','','%.2f')) for species in dust_species]))
        pars_units.update(dict([('T_DESA_' + species,('T_{desA,' + species + '}','','%.3f')) for species in dust_species]))
        pars_units.update(dict([('T_DESB_' + species,('T_{desB,' + species + '}','','%.3f')) for species in dust_species]))
        pars_units.update(dict([('T_DES_' + species,('T_{des,' + species + '}','','%.3f')) for species in dust_species]))
        pars_units.update(dict([('T_MIN_' + species,('T_{min,' + species + '}','K','%i')) for species in dust_species]))
        pars_units.update(dict([('T_MAX_' + species,('T_{max,' + species + '}','K','%i')) for species in dust_species]))
        pars_units.update(dict([('R_MIN_' + species,('R_{min,' + species + '}','R_*','%.2f')) for species in dust_species]))
        pars_units.update(dict([('R_MAX_' + species,('R_{max,' + species + '}','R_*','%.2f')) for species in dust_species]))
        '''