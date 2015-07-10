# -*- coding: utf-8 -*-

"""
A tool for comparing continuum divided data and models in a given wavelength 
range. Useful for fitting dust emission and absorption features.

Author: R. Lombaert

"""

import os,types
from scipy.integrate import trapz

import cc.path
from cc.tools.numerical import Interpol
from cc.data import Data
from cc.modeling.codes import MCMax
from cc.plotting import Plotting2
from cc.tools.io import DataIO

class ContinuumDivision(object):
    
    '''
    A class for dividing a spectrum by the continuum and comparing models and 
    data.
    
    '''
    
    def __init__(self,star_grid=[],spec=[],franges=[2.6,2.85,3.3,3.7],plot=0,\
                 func='power',cfg=''):
        
        '''
        Initializing a ContinuumDivision instance.
        
        @keyword star_grid: The parameter sets for which the fitting is done
                            
                            (default: [])
        @type star_grid: list[Star()]
        @keyword spec: information concerning dust data available for star
        
                       (default: [])
        @type spec: list[Sed()]
        @keyword franges: The fitting ranges for this instance. 4 values in mic,
                          delimiting 2 parts of the spectrum blueward and 
                          redward of the dust feature. Default is for 3.1 
                          micron ice feature.
                          Can also give a list of lists, in which case the 
                          franges are defined for every Star() and Sed() 
                          included, such that 
                          len(franges) == len(spec)+len(star_grid). The model
                          franges are listed first, then the sed franges.
                          Can also be given through the cfg file, and takes 
                          priority in that case. 
                        
                          (default: [2.6,2.85,3.3,3.7])
        @type franges: list[float] or list[list[float]]
        @keyword plot: Show the continuum division and fitting plots.
        
                       (default: 0)
        @type plot: bool
        @keyword func: The function used for fitting the continuum. Can be a 
                        list of functions as well, of len == len(star_grid) + 
                        len(spec)
                        Can also be given through the cfg file, and takes 
                        priority in that case. 
                        
                        (default: power)
        @type func: string or list[string]
        @keyword cfg: a configuration file for Plotting2.py. 
        @type cfg: string
        
        '''
    
        self.star_grid = star_grid
        self.spec = spec
        self.cont_division = dict()
        self.eq_width = dict()
        self.plot = plot
        self.cfg = cfg
        if cfg:
            cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
        else:
            cfg_dict = dict()
        if cfg_dict.has_key('franges'):
            franges = cfg_dict['franges']
            if type(franges[0]) is types.TupleType:
                franges = [list(f) for f in franges]
        if cfg_dict.has_key('func'):
            func = cfg_dict['func']
        if type(franges[0]) is types.ListType:
            if not len(franges) == (len(spec) + len(star_grid)):
                print 'Not enough sets of franges defined for all Star() and'+\
                      ' Sed() objects. Taking first set for all.'
                self.franges = [sorted(franges[0]) 
                                for i in range(len(spec)+len(star_grid))]
            else:
                self.franges = [sorted(fr) for fr in franges]
        else:
            self.franges = [sorted(franges) 
                            for i in range(len(spec)+len(star_grid))]
        
        self.frmin = min([min(fr) for fr in self.franges])
        self.frmax = max([max(fr) for fr in self.franges])
        if type(func) is types.ListType:
            if not len(func) == (len(spec) + len(star_grid)):
                print 'Not enough functions defined for all Star() and'+\
                      ' Sed() objects. Taking first func for all.'
                self.func = [func[0] for i in range(len(spec)+len(star_grid))]
            else:
                self.func = func
        else:
            self.func = [func 
                         for i in range(len(spec)+len(star_grid))]



    def show(self):
        
        '''
        Show the model calculations in a plot.
                
        '''
        
        x = []
        y = []
        keys = []
        for i,sed in enumerate(self.spec):
            x.append(self.cont_division['sws%i'%i]['w_feat'])
            y.append(self.cont_division['sws%i'%i]['f_division'])
            keys.append(sed.star_name_plots)        
        for i,star in enumerate(self.star_grid):
            if star['LAST_MCMAX_MODEL']:
                x.append(self.cont_division[star['LAST_MCMAX_MODEL']]\
                                           ['w_feat'])
                y.append(self.cont_division[star['LAST_MCMAX_MODEL']]\
                                           ['f_division'])
                #star = [s 
                    #for s in self.star_grid 
                    #if k == s['LAST_MCMAX_MODEL']][0]
                keys.append(star['LAST_MCMAX_MODEL'].replace('_','\_'))
                #keys.append(', '.join(['%s = %.2f'%(par.replace('_','\_'),\
                                                #star[par]) 
                                    #for par in star.keys() 
                                    #if par[0:2] == 'A_' and \
                                        #self.species in par]))
        fn = Plotting2.plotCols(x=x,y=y,xmin=self.frmin*0.9,\
                                xmax=self.frmax*1.1,keytags=keys,\
                                key_location=(0.6,0.02),cfg=self.cfg)
        if fn <> None:
            print 'Your Continuum Division plot can be found at '
            print fn



    def calcEqWidth(self,dtype,frindex):
        
        '''
        Calculate the equivalent width after continuum division.
        
        @param dtype: data type (only 'model','sws' for now)
        @type dtype: string
        @param frindex: The index in the franges list for this entry.
        @type frindex: int
        
        '''
        
        fr1 = self.franges[frindex][0]
        fr4 = self.franges[frindex][3]
        w = self.cont_division[dtype]['w_feat']
        w_feat = w[(w>fr1)*(w<fr4)]
        f_feat = self.cont_division[dtype]['f_division'][(w>fr1)*(w<fr4)]
        self.eq_width[dtype] = trapz(x=w_feat,y=(1-f_feat))

        
    def prepareModels(self):
        
        '''
        Prepare models for dust feature continuum division.
        
        The equivalent width of the feature is also calculated here.
        
        '''
        
        model_ids = [s['LAST_MCMAX_MODEL'] 
                     for s in self.star_grid
                     if s['LAST_MCMAX_MODEL']]
        if model_ids: 
            rt_sed = self.star_grid[0]['RT_SED']
            path_mcmax = self.star_grid[0].path_mcmax
        for i,model_id in enumerate(model_ids):
            w,f = MCMax.readModelSpectrum(path_mcmax,model_id,rt_sed)
            self.divideContinuum(w,f,dtype=model_id,frindex=i)
            self.calcEqWidth(dtype=model_id,frindex=i)
                
        
    def prepareData(self):
        
        '''
        Prepare data for dust feature continuum division.
        
        For now only SWS is available.
        
        The equivalent width of the feature is also calculated here.
        
        '''
        
        for i,sp in enumerate(self.spec):
            sws_type = [k for k in sp.data.keys() if 'SWS' in k][0]
            w = sp.data[sws_type][0]
            f = sp.data[sws_type][1]
            self.divideContinuum(w,f,dtype='sws%i'%i,\
                                 frindex=len(self.star_grid)+i)
            self.calcEqWidth(dtype='sws%i'%i,frindex=len(self.star_grid)+i)        
                
        
    def divideContinuum(self,w,f,dtype,frindex):
        
        '''
        Divide flux by the continuum flux in a dust feature. 
        
        @param w: The wavelength grid
        @type w: list/array
        @param f: The flux grid
        @type f: list/array
        @param dtype: data type (only 'model','sws' for now)
        @type dtype: string
        @param frindex: The index in the franges list for this entry.
        @type frindex: int
        
        @keyword plot: Show a plot of the continuum division 
        @type plot: bool
        
        '''
        
        fr1,fr2 = self.franges[frindex][0],self.franges[frindex][1]
        fr3,fr4 = self.franges[frindex][2],self.franges[frindex][3]
        w_in = list(w[(w>fr1) * (w<fr2)]) + list(w[(w>fr3) * (w<fr4)])
        f_in = list(f[(w>fr1) * (w<fr2)]) + list(f[(w>fr3) * (w<fr4)])
        w_cont = w[(w>self.frmin*0.9)*(w<self.frmax*1.1)]
        f_ori = f[(w>self.frmin*0.9)*(w<self.frmax*1.1)]
        f_cont = Interpol.fitFunction(x_in=w_in,y_in=f_in,x_out=w_cont,\
                                      func=self.func[frindex])
        
        self.cont_division[dtype] = dict()
        self.cont_division[dtype]['w_feat'] = w_cont
        self.cont_division[dtype]['f_feat'] = f_ori
        self.cont_division[dtype]['w_fitsel_feat'] = w_in
        self.cont_division[dtype]['f_fitsel_feat'] = f_in
        self.cont_division[dtype]['f_interp'] = f_cont
        self.cont_division[dtype]['f_division'] = f_ori/f_cont
        
        if self.plot:
            x = [w_cont,w_cont]
            y = [f_ori,f_cont]
            Plotting2.plotCols(x=x,y=y,xmin=self.frmin*0.9,xmax=self.frmax*1.1,\
                               ylogscale=0,xlogscale=0)
        