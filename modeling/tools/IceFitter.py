# -*- coding: utf-8 -*-

"""
Comparing the continuum subtracted data and models in the 3.1 micron H2O ice feature.

Author: R. Lombaert

"""

from scipy.optimize import leastsq

from cc.tools.numerical import Interpol
from cc.data import Data
from cc.modeling.codes import MCMax
from cc.plotting import Plotting2

class IceFitter(object):
    
    '''
    A class for fitting the ice feature in SWS data. 
    
    '''
    
    def __init__(self,star_grid,spec,plot=0):
        
        '''
        Initializing an IceFitter instance.
        
        @param star_grid: The parameter sets for which the fitting is done
        @type star_grid: list[Star()]
        @param spec: information concerning dust data available for star
        @type spec: Sed()
        
        @keyword plot: Show the continuum division and fitting plots.
        @type plot: bool
        
        '''
        
        self.path_mcmax = star_grid[0].path_mcmax
        self.path_gastronoom = star_grid[0].path_mcmax 
        self.path_combocode = star_grid[0].path_combocode
        self.star_grid = star_grid
        self.spec = spec
        self.cont_division = dict()
        self.plot = plot


    def show(self,cfg=''):
        
        '''
        Show the model calculations in a plot.
        
        @keyword cfg: a configuration file for Plotting2.py. 
        @type cfg: string
        
        '''
        
        x = []
        y = []
        keys = []
        for k in sorted(self.cont_division.keys()):
            x.append(self.cont_division[k]['w_ice'])
            y.append(self.cont_division[k]['f_division'])
            if 'model' in k:
                star = [s for s in self.star_grid if k == s['LAST_MCMAX_MODEL']][0]
                keys.append(', '.join(['%s = %.2f'%(par.replace('_','\_'),star[par]) 
                                       for par in star.keys() 
                                       if par[0:2] == 'A_' and 'H2O' in par]))
            else:
                keys.append(k)
        Plotting2.plotCols(x=x,y=y,xmin=2.5,xmax=4,vert_lines=[2.6,2.85,3.3,3.7],\
                           keytags=keys,key_location=(0.6,0.02),cfg=cfg)
        
        
    def prepareModels(self):
        
        '''
        Prepare model information for ice feature continuum subtraction.
        
        '''
        
        model_ids = [s['LAST_MCMAX_MODEL'] 
                     for s in self.star_grid
                     if s['LAST_MCMAX_MODEL']]
        rt_sed = self.star_grid[0]['RT_SED']
        for model_id in model_ids:
            w,f = MCMax.readModelSpectrum(self.path_mcmax,model_id,rt_sed)
            self.divideContinuum(w,f,dtype=model_id)
        
        
        
        
    def prepareSWS(self):
        
        '''
        Prepare SWS data for ice feature continuum subtraction.
        
        '''
        
        sws_type = [k for k in self.spec.data.keys() if 'SWS' in k][0]
        w = self.spec.data[sws_type][0]
        f = self.spec.data[sws_type][1]
        self.divideContinuum(w,f,dtype='sws')
                
        
    def divideContinuum(self,w,f,dtype):
        
        '''
        Divide flux by the continuum flux in the 3.1 micron ice feature. 
        
        @param w: The wavelength grid
        @type w: list/array
        @param f: The flux grid
        @type f: list/array
        @param dtype: data type ('model','sws')
        @type dtype: string
        
        @keyword plot: Show a plot of the continuum division 
        @type plot: bool
        
        '''
        
        w_in = list(w[(w>2.6) * (w<2.85)]) + list(w[(w>3.3) * (w<3.7)])
        f_in = list(f[(w>2.6) * (w<2.85)]) + list(f[(w>3.3) * (w<3.7)])
        w_cont = w[(w>2.6)*(w<3.7)]
        f_ori = f[(w>2.6)*(w<3.7)]
        f_cont = Interpol.fitFunction(x_in=w_in,y_in=f_in,x_out=w_cont,\
                                      func='power')
        
        self.cont_division[dtype] = dict()
        self.cont_division[dtype]['w_ice'] = w_cont
        self.cont_division[dtype]['f_ice'] = f_ori
        self.cont_division[dtype]['w_fitsel_ice'] = w_in
        self.cont_division[dtype]['f_fitsel_ice'] = f_in
        self.cont_division[dtype]['f_interp'] = f_cont
        self.cont_division[dtype]['f_division'] = f_ori/f_cont
        
        if self.plot:
            x = [w_cont,w_cont,w_cont]
            y = [f_ori,f_cont,f_ori/f_cont]
            Plotting2.plotCols(x=x,y=y,xmin=2.5,xmax=4.0,ymin=0.5,ylogscale=1,\
                               xlogscale=1,vert_lines=[2.6,2.85,3.3,3.7])
        