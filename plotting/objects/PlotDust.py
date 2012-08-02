# -*- coding: utf-8 -*-

"""
A plotting environment for dust information such as SEDs and all that is 
associated with that.

Author: R. Lombaert

"""

import os
import subprocess
from scipy import array

from cc.plotting.PlottingSession import PlottingSession
from cc.plotting import Plotting2
from cc.tools.io import DataIO
from cc.modeling.objects import Star
from cc.modeling.codes import MCMax



class PlotDust(PlottingSession):
    
    """ 
    Plotting environment for SEDs and all dust parameters.
    
    """
    
    def __init__(self,star_name='model',sed=None,\
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
        
        '''

        super(PlotDust, self).__init__(star_name=star_name,\
                                       path_combocode=path_combocode,\
                                       path=path_mcmax,\
                                       code='MCMax',\
                                       inputfilename=inputfilename)
        self.sed = sed
                


    def plotSed(self,star_grid,cfg='',iterative=0,no_models=0):
        
        """ 
        Creating an SED with 1 or more models and data. 
        
        Includes data preparation on the spot.
        
        @param star_grid: list of Star() models to plot. If star_grid is [], 
                          only data are plotted.
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
                
        """
        
        if self.sed is None:
            print 'No PATH_SED given. Cannot plot SED. Aborting...'
            return
        print '***********************************'
        print '** Creating SED plot.'
        
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
        if cfg_dict.has_key('no_models'):
            no_models = cfg_dict['no_models']
        ccpath = os.path.join(self.path_combocode,'Data')
        data_labels = dict([(dt,(n,ls))
                            for n,dt,ls in zip(DataIO.getInputData(path=ccpath,\
                                                        keyword='PLOT_NAMES',\
                                                        filename='Dust.dat',\
                                                        remove_underscore=1),\
                                               DataIO.getInputData(path=ccpath,\
                                                        keyword='DATA_TYPES',\
                                                        filename='Dust.dat'),\
                                               DataIO.getInputData(path=ccpath,\
                                                        keyword='LINE_TYPES',\
                                                        filename='Dust.dat'))])
        
        #- filename settings and copying inputfiles to plot output folder
        filename = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                'stars',self.star_name,self.plot_id,\
                                'SED_%s'%self.star_name)
        if iterative:
            filename = filename + '_iterative_%i'%iterative
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
        for k,(w,f) in sorted([datatype 
                               for datatype in self.sed.data.items()
                               if 'PHOT' not in datatype[0].upper()]):
             keytags.append(data_labels[k][0])
             data_x.append(self.sed.data[k][0])
             data_y.append(self.sed.data[k][1])
             line_types.append(data_labels[k][1])
        
        for k,(w,f) in sorted([datatype 
                               for datatype in self.sed.data.items()
                               if 'PHOT' in datatype[0].upper()]):
             keytags.append(data_labels[k][0])
             data_x.append(self.sed.data[k][0])
             data_y.append(self.sed.data[k][1])
             line_types.append(data_labels[k][1])
        
        #- Collect model data as well as keytags and set line types
        model_ids = [s['LAST_MCMAX_MODEL'] 
                     for s in star_grid
                     if s['LAST_MCMAX_MODEL']]
        #- Only if the model_ids list is not empty, MCMax models are available
        #- Otherwise the ray tracing keyword is unnecessary.
        if no_models:
            model_ids = []
        if model_ids: 
            rt_sed = star_grid[0]['RT_SED']
        for model_id in model_ids:
            w,f = MCMax.readModelSpectrum(self.path,model_id,rt_sed)
            data_x.append(w)
            data_y.append(f)
            keytags.append(model_id.replace('_','\_'))
        line_types += [0]*len(star_grid)
        keytags = [tag.replace('#','') for tag in keytags]
        try:
            ymax = 1.3*max([max(dy) for dy in data_y])
        except ValueError:
            pass
        try:    
            ymin = 0.5*min([min(dy) for dy in data_y])
        except ValueError:
            pass        
        filename = Plotting2.plotCols(x=data_x,y=data_y,filename=filename,\
                                      figsize=(20,10),number_subplots=1,\
                                      plot_title=plot_title,fontsize_axis=20,\
                                      keytags=keytags,fontsize_title=24,\
                                      linewidth=3,key_location=(0.0,0.75),\
                                      xlogscale=1,transparent=0,cfg=cfg,\
                                      line_types=line_types,ylogscale=0,\
                                      fontsize_ticklabels=20,fontsize_key=18,\
                                      xmin=2,xmax=200,ymin=ymin,ymax=ymax,\
                                      extension='.pdf')
        print '** Your SED plots can be found at:'
        print filename
        print '***********************************'
         

                                    
    def plotTemp(self,star_grid=[],models=[],powerlaw=[0.4],cfg=''):
        
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
        @keyword powerlaw: A list of power laws to include on the plot. If [], 
                           no power law is included. Power laws are taken from 
                           star_grid[0].
                                
                           See Thesis p32, where power is p in 
                           T(r) = T_eff*(2*r/R_STAR)**(-p).
                
                           (default: [0.4])
        @type powerlaw: list        
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
        radii = []
        temps = []
        keytags = []
        for star in star_grid:
            rad,temp,key = star.getDustTemperature() 
            radii.append(rad)
            temps.append(temp)
            keytags.append(key)
        if powerlaw:      #take T_STAR from the logfile of model in models
            for power in powerlaw:
                rad,temp,key = star_grid[0].getDustTemperaturePowerLaw(power)
                radii.append(rad)
                temps.append(temp)
                keytags.append(key)
        filename = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                'stars',self.star_name,self.plot_id,\
                                'Td_gastronoom')
        title = 'Average Dust Temperature Stratification for %s'\
                %(self.star_name_plots)
        filename = Plotting2.plotCols(x=radii,y=temps,filename=filename,\
                                      yaxis='$T_\mathrm{d}$ (K)',\
                                      plot_title=title,xaxis='$R$ (cm)',\
                                      key_location=(0.05,0.05),cfg=cfg,\
                                      xlogscale=1,ylogscale=0,fontsize_key=20,\
                                      keytags=keytags)
        print '** Your plots can be found at:'
        print filename
        print '***********************************'
            


    def plotTempSpecies(self,star_grid=[],models=[],include_total=1,\
                        powerlaw=[0.4],cfg=''):
        
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
        @keyword powerlaw: A list of power laws to include on the plot. If [], 
                           no power law is included. Power laws are taken from 
                           star_grid[0].
                                
                           See Thesis p32, where power is p in 
                           T(r) = T_eff*(2*r/R_STAR)**(-p).
                
                           (default: [0.4])
        @type powerlaw: list(float)
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
            #- Requires DUST_LIST and T_CONTACT to be taken from the log file. 
            #- It's possible, but needs some programming
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
        if cfg_dict.has_key('powerlaw'):
            powerlaw = cfg_dict['powerlaw']
        plot_filenames = []
        for star in star_grid:
            if not int(star['T_CONTACT']):
                radii,temps,keytags = star.getDustTemperatureSpecies()
                vert_lines = []
            else:
                include_total = 1
                print 'Thermal contact is on. All dust species share the ' + \
                      'same temperature profile. Vertical lines indicate ' + \
                      'inner radii of dust species.'
                radii, temps, keytags = [], [], []
                vert_lines = [star['R_DES_%s'%d]*star.r_solar*star['R_STAR'] 
                              for d in star['DUST_LIST']]
            if include_total:
                rad, temp, key = star.getDustTemperature()
                radii.append(rad[rad>star['R_INNER_GAS']\
                                *star.r_solar*star['R_STAR']])
                temps.append(temp[rad>star['R_INNER_GAS']\
                                *star.r_solar*star['R_STAR']])
                keytags.append(key)
            if powerlaw:
                for power in powerlaw:
                     rad,temp,key = star.getDustTemperaturePowerLaw(power)
                     radii.append(rad)
                     temps.append(temp)
                     keytags.append(key)
            filename = os.path.join(os.path.expanduser('~'),'MCMax',self.path,\
                                    'stars',self.star_name,self.plot_id,\
                                    'Td_species_%s'%star['LAST_MCMAX_MODEL'])
            extension = '.eps'
            title = 'Dust Temperature in %s'%(self.star_name_plots)
            plot_filenames.append(Plotting2.plotCols(x=radii,y=temps,cfg=cfg,\
                        filename=filename,xaxis='$r$ (cm)',\
                        yaxis='$T_\mathrm{d}$ (K)',keytags=keytags,\
                        key_location=(.65,.45),ymin=20,ymax=3000,\
                        xmax=star['R_OUTER_DUST']*star.r_solar*star['R_STAR'],\
                        xmin=star['R_STAR']*star.r_solar,\
                        xlogscale=1,ylogscale=1,fontsize_key=18,\
                        figsize=(12.5,8),transparent=0,linewidth=4,\
                        plot_title='',fontsize_title=22,\
                        extension=extension,fontsize_axis=26,\
                        fontsize_ticklabels=26,\
                        vert_lines=[star['R_INNER_DUST']\
                                        *star.r_solar*star['R_STAR']]))
        if len(plot_filenames) != len(star_grid):
            print 'At least one of the models does not yet have a MCMax model.'        
        if plot_filenames[0][-4:] == '.pdf':
            new_filename = os.path.join(os.path.expanduser('~'),'MCMax',\
                                        self.path,'stars',self.star_name,\
                                        self.plot_id,'Td_species.pdf')
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
                                          xaxis='$\lambda$ ($\mu$m)',cfg=cfg,\
                                          yaxis='$\kappa_\lambda$ (cm$^2$/g)',\
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
                opacities = [(opacities[i]+opacities[i+len(star['DUST_LIST'])]) 
                             for i,species in enumerate(star['DUST_LIST'])]
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
                           for sp in star['DUST_LIST']]
                filenames.append(Plotting2.plotCols(x=wave,y=opacities,\
                                 xaxis='$\lambda$ ($\mu$m)',\
                                 yaxis='$\kappa_\lambda$ (cm$^2$/g)',\
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
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
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
        
        star_grid = Star.makeStars(models=models,star_name=self.star_name,\
                                   code='MCMax',id_type='MCMax',path=self.path)
        for star,model in zip(star_grid,models):    
            grid_shape = star.getMCMaxOutput(incr=1,keyword='NGRAINS',\
                                             single=0)[0]
            star.update({'PATH_DUST_DATA':data_path,\
                         'NTHETA':int(grid_shape[1]),\
                         'NRAD':int(grid_shape[0]),\
                         'T_STAR':float(star.getMCMaxOutput(incr=0,\
                                                keyword='STELLAR TEMPERATURE',\
                                                filename='log.dat',\
                                                single=0)[0][2]),\
                         'R_STAR':float(star.getMCMaxOutput(incr=0,\
                                                keyword='STELLAR RADIUS',
                                                filename='log.dat',\
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