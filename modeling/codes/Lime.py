# -*- coding: utf-8 -*-

"""
Preparing input for LIME.

Author: R. Lombaert

"""

import os

from cc.tools.io import DataIO

def prepInput(star,path,repl_str=''):
    
    '''
    Prepare inputfiles for LIME from a GASTRoNOoM and MCMax model. 
    
    Input is taken from the model ouput and written into files. The output
    folder can be chosen. The model_id tag can be replaced with an arbitrary
    string.
    
    Input is converted to SI units, except opacities, which are in cm2/g.
    
    @param star: The model object including the GASTRoNOoM and MCMax model ids.
    @type star: Star()
    @param path: The target folder for the files.
    @type path: str
    
    @keyword repl_str: Replacement string for the model_id tag.
    
                       (default: '')
    @type repl_str: str
    
    '''
    
    #-- First opacities and t_dust, which are not part of the GASTRoNOoM output
    mcmid = star['LAST_MCMAX_MODEL']
    fnopac = os.path.join(path,'opac_%s.dat'%(repl_str and repl_str or mcmid))
    DataIO.writeCols(fnopac,star.readWeightedKappas())
    
    #-- Write dust temperature to a file.
    fntd = os.path.join(path,'td_%s.dat'%(repl_str and repl_str or mcmid))
    rad,td = star.getDustTemperature()
    rad = rad*star['R_STAR']*star.Rsun*10**-2
    DataIO.writeCols(fntd,[rad,td])
    
    #-- Finally the gas properties
    props = dict([('nh2','N(H2)'),('vg','VEL'),('rad','RADIUS')])
    if not repl_str: repl_str = star['LAST_GASTRONOOM_MODEL']
    fngas = star.getCoolFn('fgr_all')
    fncogas = star.getCoolFn('1',mstr='12C16O')
    
    #-- Radius for all files
    rad = DataIO.getGastronoomOutput(fngas,keyword='RADIUS',return_array=1)
    rad = rad*10**-2
    
    #-- h2 number density
    nh2 = DataIO.getGastronoomOutput(fngas,keyword='N(H2)',return_array=1)
    nh2 = nh2*10**6
    fnnh2 = os.path.join(path,'nh2_%s.dat'%repl_str)
    DataIO.writeCols(fnnh2,[rad,nh2])
    
    #-- Velocity profile
    vel = DataIO.getGastronoomOutput(fngas,keyword='VEL',return_array=1)
    vel = vel*10**-2
    fnvel = os.path.join(path,'vg_%s.dat'%repl_str)
    DataIO.writeCols(fnvel,[rad,vel])
    
    #-- CO abundance
    rad = DataIO.getGastronoomOutput(fncogas,keyword='RADIUS',return_array=1)
    rad = rad*star['R_STAR']*star.Rsun*10**-2
    ah2 = DataIO.getGastronoomOutput(fncogas,keyword='N(H2)',return_array=1)
    aco = DataIO.getGastronoomOutput(fncogas,keyword='N(MOLEC)',key_index=8,\
                                     return_array=1)
    aco = aco/ah2
    fnaco = os.path.join(path,'aco_%s.dat'%repl_str)
    DataIO.writeCols(fnaco,[rad,aco])
    
    
    