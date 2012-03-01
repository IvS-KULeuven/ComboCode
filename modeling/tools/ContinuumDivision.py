# -*- coding: utf-8 -*-

"""
A tool for comparing continuum divided data and models in a given wavelength 
range. Useful for fitting dust emission and absorption features.

Author: R. Lombaert

"""

from cc.tools.numerical import Interpol
from cc.data import Data
from cc.modeling.codes import MCMax
from cc.plotting import Plotting2

class ContinuumDivision(object):
    
    '''
    A class for dividing a spectrum by the continuum and comparing models and 
    data.
    
    '''
    
    def __init__(self,star_grid,spec,franges=[2.6,2.85,3.3,3.7],plot=0,func='power'):
        
        '''
        Initializing a ContinuumDivision instance.
        
        @param star_grid: The parameter sets for which the fitting is done
        @type star_grid: list[Star()]
        @param spec: information concerning dust data available for star
        @type spec: Sed()
        
        @keyword franges: The fitting ranges for this instance. 4 values in mic, 
                        delimiting 2 parts of the spectrum blueward and redward
                        of the dust feature. Default is for 3.1 micron ice 
                        feature.
                        
                        (default: [2.6,2.85,3.3,3.7])
        @type franges: list[float]
        @keyword plot: Show the continuum division and fitting plots.
        
                       (default: 0)
        @type plot: bool
        @keyword func: The function used for fitting the continuum
        
                       (default: power)
        @type func: string
        
        '''
        
        self.path_mcmax = star_grid[0].path_mcmax
        self.path_gastronoom = star_grid[0].path_mcmax 
        self.path_combocode = star_grid[0].path_combocode
        self.star_grid = star_grid
        self.spec = spec
        self.cont_division = dict()
        self.plot = plot
        self.franges = sorted(franges)
        self.func = func
        if len(franges) != 4:
            print 'The franges keyword is not defined correctly. Aborting ' + \
                  'ContinuumDivision call.'
            return



    def show(self,cfg=''):
        
        '''
        Show the model calculations in a plot.
        
        @keyword cfg: a configuration file for Plotting2.py. 
        @type cfg: string
        
        '''
        
        fr1,fr2 = self.franges[0],self.franges[1]
        fr3,fr4 = self.franges[2],self.franges[3]
        x = []
        y = []
        keys = []
        for k in sorted([ik for ik in self.cont_division.keys()
                            if 'model' not in ik]):
            x.append(self.cont_division[k]['w_feat'])
            y.append(self.cont_division[k]['f_division'])
            keys.append(k)        
        for star in self.star_grid:
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
        fn = Plotting2.plotCols(x=x,y=y,xmin=fr1*0.9,xmax=fr4*1.1,\
                                keytags=keys,vert_lines=self.franges,\
                                key_location=(0.6,0.02),cfg=cfg)
        if fn <> None:
            print 'Your Continuum Division plot can be found at '
            print fn


        
    def prepareModels(self):
        
        '''
        Prepare models for dust feature continuum division.
        
        '''
        
        model_ids = [s['LAST_MCMAX_MODEL'] 
                     for s in self.star_grid
                     if s['LAST_MCMAX_MODEL']]
        rt_sed = self.star_grid[0]['RT_SED']
        for model_id in model_ids:
            w,f = MCMax.readModelSpectrum(self.path_mcmax,model_id,rt_sed)
            self.divideContinuum(w,f,dtype=model_id)
        
                
        
    def prepareData(self):
        
        '''
        Prepare data for dust feature continuum division.
        
        For now only SWS is available.
        
        '''
        
        sws_type = [k for k in self.spec.data.keys() if 'SWS' in k][0]
        w = self.spec.data[sws_type][0]
        f = self.spec.data[sws_type][1]
        self.divideContinuum(w,f,dtype='sws')
                
                
        
    def divideContinuum(self,w,f,dtype):
        
        '''
        Divide flux by the continuum flux in a dust feature. 
        
        @param w: The wavelength grid
        @type w: list/array
        @param f: The flux grid
        @type f: list/array
        @param dtype: data type (only 'model','sws' for now)
        @type dtype: string
        
        @keyword plot: Show a plot of the continuum division 
        @type plot: bool
        
        '''
        
        fr1,fr2 = self.franges[0],self.franges[1]
        fr3,fr4 = self.franges[2],self.franges[3]
        w_in = list(w[(w>fr1) * (w<fr2)]) + list(w[(w>fr3) * (w<fr4)])
        f_in = list(f[(w>fr1) * (w<fr2)]) + list(f[(w>fr3) * (w<fr4)])
        w_cont = w[(w>fr1)*(w<fr4)]
        f_ori = f[(w>fr1)*(w<fr4)]
        f_cont = Interpol.fitFunction(x_in=w_in,y_in=f_in,x_out=w_cont,\
                                      func=self.func)
        
        self.cont_division[dtype] = dict()
        self.cont_division[dtype]['w_feat'] = w_cont
        self.cont_division[dtype]['f_feat'] = f_ori
        self.cont_division[dtype]['w_fitsel_feat'] = w_in
        self.cont_division[dtype]['f_fitsel_feat'] = f_in
        self.cont_division[dtype]['f_interp'] = f_cont
        self.cont_division[dtype]['f_division'] = f_ori/f_cont
        
        if self.plot:
            x = [w_cont,w_cont,w_cont]
            y = [f_ori,f_cont,f_ori/f_cont]
            Plotting2.plotCols(x=x,y=y,xmin=fr1*0.9,xmax=fr4*1.1,ylogscale=1,\
                               xlogscale=1,vert_lines=self.franges)
        