# -*- coding: utf-8 -*-

"""
Preparing input for LIME.

Author: R. Lombaert

"""

import os
import numpy as np

from cc.tools.io import DataIO
from cc.data import Data

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
    rad = star.getDustRad(unit='m')
    td = star.getDustTemperature()
    DataIO.writeCols(fntd,[rad,td])
    
    #-- Finally the gas properties
    if not repl_str: repl_str = star['LAST_GASTRONOOM_MODEL']
    
    #-- Radius for nh2 and vel
    rad = star.getGasRad(unit='m',ftype='fgr_all')
    if rad.size > 10000:
        icutoff = np.argmin(abs(rad-1e14))
        rad = Data.reduceArray(rad,20,1e14,'remove')
    else: 
        icutoff = None
    
    #-- h2 number density
    nh2 = star.getGasNumberDensity(ftype='fgr_all')
    nh2 = nh2*10**6
    if icutoff <> None: 
        nh2 = Data.reduceArray(nh2,20,nh2[icutoff],'remove')
    fnnh2 = os.path.join(path,'nh2_%s.dat'%repl_str)
    DataIO.writeCols(fnnh2,[rad,nh2])
    
    #-- Velocity profile
    vel = star.getGasVelocity(ftype='fgr_all')
    vel = vel*10**-2
    if icutoff <> None: 
        vel = Data.reduceArray(vel,20,vel[icutoff],'remove')
    fnvel = os.path.join(path,'vg_%s.dat'%repl_str)
    DataIO.writeCols(fnvel,[rad,vel])
    
    #-- CO abundance
    rad = star.getGasRad(ftype='1',mstr='12C16O',unit='m')
    nh2 = star.getGasNumberDensity(ftype='1',mstr='12C16O')
    nco = star.getGasNumberDensity(ftype='1',mstr='12C16O',molecule=1)
    aco = nco/nh2
    fnaco = os.path.join(path,'aco_%s.dat'%repl_str)
    DataIO.writeCols(fnaco,[rad,aco])

    