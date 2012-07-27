# -*- coding: utf-8 -*-

"""
Some tools for measuring line profiles of emission lines.

Author: R. Lombaert

"""

from scipy import mean
from scipy import argmin
from scipy.integrate import trapz

from cc.tools.numerical import Interpol
from cc.tools.io import FitsReader, TxtReader


def integrateLPData(filename,vexp,vlsr=None,window=1.2):
    
    """
    Integrate a line profile read from a fits file or a txt file. 
    
    Requires a terminal expansion velocity to determine the window in which the
    integration is done. 
    
    If a fits file is read, the source velocity is taken from it (if available)
    If a text file is read, a source velocity needs to be given as input.
    
    @param filename: The filename to the data file of the line profile
    @type filename: string
    @param vexp: The terminal gas velocity
    @type vexp: float
    
    @keyword vlsr: the source velocity, only needed if fits file or text file
                   do not contain vlsr. 
                   
                   (default: None)
    @type vlsr: float
    @keyword window: The factor with which vexp is multiplied when selecting
                     the integration window in the velocity grid. 
                     
                     (default: 1.2)
    @type window: float
    
    @return: The integrated intensity of the profile
    @rtype: float
    
    """
    
    if filename[-5:] == '.fits':
        lprof = FitsReader.FitsReader(filename,vlsr)
    else:
        lprof = TxtReader.TxtReader(filename,vlsr)
    
    vel = lprof.getVelocity()
    flux = lprof.getFlux()
    vlsr = lprof.getVlsr()
    vexp = float(vexp)
    window = float(window)
    
    vel_sel = vel[(vel > vlsr - window*vexp)*(vel < vlsr + window*vexp)]
    flux_sel = flux[(vel > vlsr - window*vexp)*(vel < vlsr + window*vexp)]
    return trapz(x=vel_sel,y=flux_sel)
    


def getPeakLPData(filename,vlsr=None):
    
    """
    Calculate the peak value of a line profile read from a fits file or a txt 
    file. 
    
    If a fits file is read, the source velocity is taken from it (if available)
    If a text file is read, a source velocity needs to be given as input.
    
    The peak value is defined as the mean of central 5 flux points around the
    source velocity. 
    
    @param filename: The filename to the data file of the line profile
    @type filename: string
    
    @keyword vlsr: the source velocity, only needed if fits file or text file
                   do not contain vlsr. 
                   
                   (default: None)
    @type vlsr: float
    
    @return: The peak intensity of the profile
    @rtype: float
    
    """
    
    if filename[-5:] == '.fits':
        lprof = FitsReader.FitsReader(filename,vlsr)
    else:
        lprof = TxtReader.TxtReader(filename,vlsr)
    
    vel = lprof.getVelocity()
    flux = lprof.getFlux()
    vlsr = lprof.getVlsr()
    
    i_mid = argmin(abs(vel-vlsr))
    return mean(flux[i_mid-2:i_mid+3])
    