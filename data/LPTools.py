# -*- coding: utf-8 -*-

"""
Some tools for measuring line profiles of emission lines.

Author: R. Lombaert

"""

import os
import numpy as np
from scipy import mean,sqrt,log, std
from scipy import argmin,argmax
from scipy.integrate import trapz
from matplotlib import pyplot as plt

from cc.tools.io import FitsReader, TxtReader

from ivs.sigproc import fit, funclib


def readLineProfile(filename,info_path=os.path.join(os.path.expanduser('~'),\
                                                      'ComboCode','Data')):
        
    '''
    Read a line profile, independent of the type of extension. 
    
    @keyword filename: The filename to the data file of the line profile. If 
                       None a line profile object is expected. 
                       
                       (default: None)
    @type filename: string
    
    @keyword info_path: The path to the folder containing the info file on
                        stars, called Star.dat. 
                        
                        (default: ~/ComboCode/Data)
    @type info_path: string
    
    @return: The line profile
    @rtype: LPDataReader()
    
    '''
    
    if filename[-5:] == '.fits':
        lprof = FitsReader.FitsReader(filename,info_path)
    else:
        lprof = TxtReader.TxtReader(filename,info_path)
    return lprof
    
    

def integrateLPData(vexp,filename=None,lprof=None,window=1.2,\
                    info_path=os.path.join(os.path.expanduser('~'),\
                                           'ComboCode','Data')):
    
    """
    Integrate a line profile read from a fits file or a txt file. 
    
    Requires a terminal expansion velocity to determine the window in which the
    integration is done. 
    
    If a fits file is read, the source velocity is taken from it (if available)
    If a text file is read, a source velocity needs to be given as input.
    
    @param vexp: The terminal gas velocity
    @type vexp: float
    
    @keyword filename: The filename to the data file of the line profile. If 
                       None a line profile object is expected. 
                       
                       (default: None)
    @type filename: string
    @keyword lprof: A line profile object (LPDataReader or inheriting classes)
                    If None, a filename is expected!
                    
                    (default: None)
    @type lprof: LPDataReader()
    @keyword window: The factor with which vexp is multiplied when selecting
                     the integration window in the velocity grid. 
                     
                     (default: 1.2)
    @type window: float
    @keyword info_path: The path to the folder containing the info file on
                        stars, called Star.dat. 
                        
                        (default: ~/ComboCode/Data)
    @type info_path: string
    
    @return: The integrated intensity of the profile
    @rtype: float
    
    """
    
    if lprof is None:
        lprof = readLineProfile(filename,file_path)
    vel = lprof.getVelocity()
    flux = lprof.getFlux()
    vlsr = lprof.getVlsr()
    vexp = float(vexp)
    window = float(window)
    
    vel_sel = vel[(vel > vlsr - window*vexp)*(vel < vlsr + window*vexp)]
    flux_sel = flux[(vel > vlsr - window*vexp)*(vel < vlsr + window*vexp)]
    return trapz(x=vel_sel,y=flux_sel)
    


def getPeakLPData(filename=None,lprof=None,\
                  info_path=os.path.join(os.path.expanduser('~'),\
                                         'ComboCode','Data')):
    
    """
    Calculate the peak value of a line profile read from a fits file or a txt 
    file. 
    
    If a fits file is read, the source velocity is taken from it (if available)
    If a text file is read, a source velocity needs to be given as input.
    
    The peak value is defined as the mean of central 5 flux points around the
    source velocity. 
    
    @keyword filename: The filename to the data file of the line profile. If 
                       None a line profile object is expected. 
                       
                       (default: None)
    @type filename: string
    @keyword lprof: A line profile object (LPDataReader or inheriting classes)
                    If None, a filename is expected!
                    
                    (default: None)
    @type lprof: LPDataReader()
    @keyword info_path: The path to the folder containing the info file on
                        stars, called Star.dat. 
                        
                        (default: ~/ComboCode/Data)
    @type info_path: string
    
    @return: The peak intensity of the profile
    @rtype: float
    
    """
    
    if lprof is None:
        lprof = readLineProfile(filename,file_path)
    vel = lprof.getVelocity()
    flux = lprof.getFlux()
    vlsr = lprof.getVlsr()
    
    i_mid = argmin(np.abs(vel-vlsr))
    return mean(flux[i_mid-2:i_mid+3])
    
    

def fitSoftParabola(vel,flux,initial):
    
    """
    Fit a soft parabola to a line profile.
    
    @param vel: The velocity grid
    @type vel: array
    @param flux: The intensity grid
    @type flux: array
    @param initial: initial parameters [int,vlsr,vexp,gamma]
    @type initial: list
    
    @return: The model after minimization
    @type: funclib.soft_parabola
    """
    
    #-- fit only soft parabola
    #   1. setup model    
    mymodel = funclib.soft_parabola()
    #   2. Initial values: [int,vlsr,vexp,gamma] 
    mymodel.setup_parameters(values=initial)
    #   3. minimize and evaluate fit
    result = fit.minimize(vel,flux,mymodel)
    return mymodel
    
    
    
def checkLPShape(vel,flux,vlsr,vexp,show=0):
    
    """
    Check the shape of the line profile, to see if any irregular patterns are
    present. 
    
    Based on the diff() method, and the std on the result of it. If sharp 
    patterns are present, they should be picked up by this method, and an extra
    component can be included in fitting the line.
    
    Detects absorption irregularities, not emission! At least, it tries to 
    filter the emission effects out. If an emission effect is stronger than 
    an absorption effect, it should also not detect any irregularities.
    
    Still being tested!
    
    @param vel: The velocity grid
    @type vel: array
    @param flux: The flux grid
    @type flux: array
    @param vlsr: the central source velocity
    @type vlsr: float
    @param vexp: The expected gas terminal velocity
    @type vexp: float
    
    @keyword show: Show the results of the diff and std methods.
                    
                   (default: 0)
    @type show: bool
    
    @return: The details of the absorption irregularity. Assuming a Gaussian,
             a depth, mid point velocity, a width and a continuum value are 
             returned. If no irregularity is found, None is returned
    @rtype: tuple(float)
    
    """
    
    fdf = np.diff(flux)
    velfdf = vel[1:]
    fluxfdf = flux[1:]
    noisedf = std(fdf[np.abs(velfdf-(vlsr-2.*vexp))<vexp])
    fdf_large = np.abs(fdf) > 3*noisedf
    if show:
        plt.clf()
        plt.step(vel,flux,'r-',lw=3,label='Observed profile')
        plt.plot(velfdf,fdf,'b-',lw=3,label='diff(profile)')
        plt.plot([vel[0],vel[-1]],[noisedf*3,noisedf*3],'--k',label='Upper df limit')
        plt.plot([vel[0],vel[-1]],[-noisedf*3,-noisedf*3],'--k',label='Lower df limit')
        leg = plt.legend(loc='best',fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.show()
    if fdf_large.any():
        imin = argmin(fdf)
        imax = argmax(fdf)
        if imin < imax: 
            #- The FWHM of the irregularity is roughly vel_imax - vel_imin
            #  as the max and min in df give the strongest decrease and increase of 
            #  the profile
            width = velfdf[imax] - velfdf[imin] 
            #- This leads to a gaussian sigma of
            sigma = width/(2.*sqrt(2.*log(2.)))
            #- velmid is where the irregularity reaches its most extreme point
            velmid = (velfdf[imax] + velfdf[imin])/2.
            #- The depth of the irregularity is roughly equal to the difference
            #  in max and min flux in an interval between imin-1/2*width and
            #  imax+1/2*width, ie between velmid-width and velmid+width 
            interval = np.abs(velfdf-velmid) <= width
            fmax = max(fluxfdf[interval])
            fmin = min(fluxfdf[interval])
            depth = -abs(fmax - fmin)
            return (depth,velmid,sigma,0)
        else:
            #- the found change in slope is likely due to the line itself
            return None
    else:
        return None
    return fdf_large.any()



def fitLP(filename=None,lprof=None,show=0,\
          info_path=os.path.join(os.path.expanduser('~'),'ComboCode','Data')):
    
    '''
    Fit a line profile with a soft parabola, and a Gaussian component if 
    required. 
    
    The method automatically checks if a second component is needed (eg an 
    extra absorption component). An estimate of the expansion velocity (width 
    of the profile) and an improved guess of the vlsr are given. 
    
    A guess for the gas terminal velocity is returned. 
    
    @keyword filename: The filename to the data file of the line profile. If 
                       None a line profile object is expected. 
                       
                       (default: None)
    @type filename: string
    @keyword lprof: A line profile object (LPDataReader or inheriting classes)
                    If None, a filename is expected!
                    
                    (default: None)
    @type lprof: LPDataReader()
    @keyword info_path: The path to the folder containing the info file on
                        stars, called Star.dat. 
                        
                        (default: ~/ComboCode/Data)
    @type info_path: string
    @keyword show: Show the results of the fit
                    
                   (default: 0)
    @type show: bool
    
    @return: the gas terminal velocity in km/s estimated from the fit
    @rtype: float
    
    '''
    
    print '***********************'
    if filename <> None: print '** Fitting line in %s.'%filename
    if lprof is None:
        lprof = readLineProfile(filename,info_path)
    vel = lprof.getVelocity()
    flux = lprof.getFlux()
    vlsr = lprof.getVlsr()
    #-- First guess, only relevant for initial window selection. 
    #-- vexp is narrowed down later.
    vexp = 60
    #-- select line
    keep = np.abs(vel-vlsr)<=(1.5*vexp)
    velsel,fluxsel = vel[keep],flux[keep]
    #-- Initial values: [peak tmb,vlsr,vexp,gamma] 
    #   vexp ~ 20, order should be OK. Large enough, in case an absorption
    #   seems to make a line narrower
    #   For the central peak value, get a first guess from the data
    peak = getPeakLPData(lprof=lprof)
    initial = [peak,vlsr,20.,1.]
    firstguess = fitSoftParabola(velsel,fluxsel,initial)
    #-- Improve the vexp guess
    print firstguess.param2str(accuracy=5)
    pars = firstguess.get_parameters()
    vexp = pars[0][2]
    keep = np.abs(vel-vlsr)<=(3.*vexp)
    velsel,fluxsel = vel[keep],flux[keep]
    flux_firstguess = firstguess.evaluate(velsel)
    
    #-- Check whether irregularities are present in the profile. 
    #   Initial parameters for a gaussian are returned if true. 
    include_gauss = checkLPShape(velsel,fluxsel,vlsr,vexp,show=show)
    #-- Do the fit of the line again, including a gaussian if irregularities
    #   are present. 
    initial = [peak,vlsr,20.,1.]
    if include_gauss <> None:
        #-- fit soft parabola + gaussian
        #   1. setup model
        soft_parabola = funclib.soft_parabola()
        soft_parabola.setup_parameters(values=initial)
        gaussian = funclib.gauss()
        #   initial guesses assuming an interstellar absorption line from the
        #   checkLPShape method
        gaussian.setup_parameters(values=include_gauss)
        #   2. minimize and evaluate fit
        mymodel = fit.Model(functions=[soft_parabola,gaussian])
        result = fit.minimize(vel,flux,mymodel)
    else: 
        #-- fit soft parabola
        #   1. setup model
        #   2. evaluate
        soft_parabola = fitSoftParabola(velsel,fluxsel,initial)
        mymodel = soft_parabola
    pars = soft_parabola.get_parameters()
    vexp = pars[0][2]
    print mymodel.param2str(accuracy=10)
    
    #-- Collect some data for a plot.
    flux_sf = soft_parabola.evaluate(velsel)
    vel_highres = np.linspace(velsel[0],velsel[-1],10000)
    fitted_flux = mymodel.evaluate(velsel)
    fitted_flux_highres = mymodel.evaluate(vel_highres)
    
    #-- compute numerical integrations
    print('I_mb (emission line data): %f'%trapz(y=fluxsel,x=velsel))
    print('I_mb (SF -- initial guess): %f'\
          %trapz(y=flux_firstguess,x=velsel))
    print('I_mb (SF -- improved guess): %f'\
          %trapz(y=flux_sf,x=velsel))
    print('I_mb (emission line fit): %f'\
          %trapz(y=fitted_flux,x=velsel))
    
    #-- plot
    if show: 
        plt.clf()
        plt.step(velsel,fluxsel,'-r',lw=3,label='Observed profile')
        plt.plot(velsel,flux_firstguess,'b--',lw=3,label='First guess fit')
        plt.plot(vel_highres,fitted_flux_highres,'g-',lw=3,\
                 label='SP+Gauss to whole profile')
        plt.plot(velsel,flux_sf,'g--',lw=2,\
                 label='SP+Gauss to whole profile (only SP)')
        leg = plt.legend(loc='best',fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.show()
    return vexp