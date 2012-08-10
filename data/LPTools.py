# -*- coding: utf-8 -*-

"""
Some tools for measuring line profiles of emission lines.

Author: R. Lombaert

"""

import os
import numpy as np
from scipy import mean,sqrt,log, std
from scipy import argmin,argmax,array
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
    
    
    
def varyInitialFit(vel,flux,initial,index,values,\
                   function=funclib.soft_parabola,vary_window=0): 

    """
    Fit a function to a line profile for different initial guesses of a single
    parameter.
    
    The best fit function is returned. Best fit is determined by checking for 
    the smallest relative error for the parameter in question.
    
    @param vel: The velocity grid
    @type vel: array
    @param flux: The intensity grid
    @type flux: array
    @param initial: initial parameters (eg [int,vlsr,vexp,gamma] for sp)
    @type initial: list
    @param index: Index of initial parameter to be varied
    @type index: int
    @param values: The values used for the variable initial parameter
    @type vary_init: tuple
    
    @keyword function: The function to be fitted
    
                       (default: funclib.soft_parabola)
    @type function: funclib.function
    @keyword vary_window: Vary the velocity window based on the third initial
                          parameter. 1.5 for soft parabola, 3 for other functions
                          centered on the second initial parameter. 
                          (mu/vlsr and sigma/vexp respectively)
                          
                          (default: 0) 
    @type window: bool
                          
    @return: The model after minimization
    @type: funclib.soft_parabola
    
    """
    
    values, initial, vary_window = list(values), list(initial), int(vary_window)
    all_init = [[p]*len(values) for i,p in enumerate(initial) if i != index]
    all_init.insert(index,values)
    if vary_window:
        #-- The window for soft parabola can be quite narrow. For Gaussians 
        #   it should be wider, taking into account broader wings, and that 
        #   sigma/2 < vexp = fwhm/2
        window = function == funclib.soft_parabola and 1.5 or 3.
        results = [fitFunction(vel[np.abs(vel-initi[1])<=(initi[2]*window)],\
                               flux[np.abs(vel-initi[1])<=(initi[2]*window)],\
                               initi,function)
                   for initi in zip(*all_init)]
    else: 
        results = [fitFunction(vel,flux,initi,function)
                   for initi in zip(*all_init)]
    results = [fg 
               for fg in results
               if fg.get_parameters()[1][index] != 0.]
    fitvalues = array([fg.get_parameters()[0][index] 
                       for fg in results])
    fiterrors = array([fg.get_parameters()[1][index] 
                       for fg in results])
    bestresult = results[argmin(abs(fiterrors/fitvalues))]
    return bestresult
    
 

def fitFunction(x,y,initial,function):
    
    """
    Fit a function to a set of x and y values.
    
    @param vel: The x grid
    @type vel: array
    @param flux: The y grid
    @type flux: array
    @param initial: initial parameters
    @type initial: list
    @param function: The function to be fitted
    @type function: funclib.function (e.g. funclib.soft_parabola, funclib.gauss )
    
    @return: The model after minimization
    @type: funclib.######## (some function)
    """
    
    #-- fit only soft parabola
    #   1. setup model    
    #mymodel = funclib.soft_parabola()
    mymodel = function()
    #   2. Initial values: e.g. for SP [int,vlsr,vexp,gamma] 
    mymodel.setup_parameters(values=initial)
    #   3. minimize and evaluate fit
    result = fit.minimize(x,y,mymodel)
    return mymodel
    
    
    
def checkLPShape(vel,flux,vlsr,vexp,window=2.,show=0):
    
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
    
    @keyword window: The window for line selection. Value is multiplied with 
                     0.6. For the usual window of 2., this is 1.2 (SP). For the
                     Gaussian window of 3., this is 1.8.
                     
                     (default: 2.)
    @type window: float
    @keyword show: Show the results of the diff and std methods.
                    
                   (default: 0)
    @type show: bool
        
    @return: The details of the absorption irregularity. Assuming a Gaussian,
             a depth, mid point velocity, a width and a continuum value are 
             returned. If no irregularity is found, None is returned
    @rtype: tuple(float)
    
    """
    
    #-- Get diff of flux to have a measure of change in slope throughout lprof
    fdf = np.diff(flux)
    velfdf = vel[1:]
    fluxfdf = flux[1:]
    #-- Select emission line. Use 0.8 as factor instead of 0.6 (results in 1.4 
    #   for window == 2), to take into account the possible underestimation of
    #   vexp
    fdfwhereline = np.abs(velfdf-vlsr)<=window*0.8*vexp
    #-- create window for noise calculation:
    #       - close to the emission line
    #       - contains enough points to be statistically significant
    i = 2
    fdfwindow = np.abs(velfdf-vlsr)<=i*vexp
    while len(velfdf[-fdfwhereline*fdfwindow]) < 60 and i != 6:
        i += 1
        fdfwindow = np.abs(velfdf-vlsr)<=i*vexp
    #-- Check if flux outliers (outside the line) are present: 
    #-- ignore for diff noise calculations
    noiseflux = std(fluxfdf[-fdfwhereline*fdfwindow])
    fdfcleanflux = fluxfdf<3.*noiseflux
    #-- Calc noise of diff in region where 'not emission line', but also not 
    #   too far out, and where no flux outliers are present
    noisedf = std(fdf[-fdfwhereline*fdfwindow*fdfcleanflux])

    #-- Show the result
    if show:
        plt.clf()
        velfdfplot = velfdf[fdfwindow]
        fluxfdfplot = fluxfdf[fdfwindow]
        plt.step(velfdfplot,fluxfdfplot,'-r',lw=3,where='mid',\
                 label='Observed profile')
        cleanflux = fluxfdf.copy() 
        cleanflux[-(-fdfwhereline*fdfwindow*fdfcleanflux)] = np.nan
        plt.step(velfdf,cleanflux,'-k',lw=2,where='mid',label='Clean noise')
        plt.plot(velfdf[fdfwindow],fdf[fdfwindow],'b-',lw=3,\
                 label='diff(profile)')
        plt.plot([velfdfplot[0],velfdfplot[-1]],[noisedf*3,noisedf*3],'--k',\
                 label='Upper df limit')
        plt.plot([velfdfplot[0],velfdfplot[-1]],[-noisedf*3,-noisedf*3],'--k',\
                 label='Lower df limit')
        ax = plt.gca()
        ax.axvline(x=vlsr-window*0.8*vexp,linewidth=2, ls='--', color='k')
        ax.axvline(x=vlsr+window*0.8*vexp,linewidth=2, ls='--', color='k')    
        leg = plt.legend(loc='best',fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.show()
    
    #-- See where the change in profile is significantly large, but only do 
    #   this where the line is located
    fdf_large = np.abs(fdf[fdfwhereline]) >= 3.*noisedf       
    #-- If any large values are detected, check the min and max values
    if fdf_large.any():
        imin = argmin(fdf[fdfwhereline])
        imax = argmax(fdf[fdfwhereline])
        #-- Check if both the min and max are actually larger than 3*noisedf
        if not (fdf_large[imin] and fdf_large[imax]):
            return None
        #-- If imin > imax, it's likely due to a very sharp profile: You do not 
        #   want to change these.
        #   Onle select changes in profile when they are reversed!
        if imin < imax: 
            #- The FWHM of the irregularity is roughly vel_imax - vel_imin
            #  as the max and min in df give the strongest decrease and increase of 
            #  the profile
            velline = velfdf[fdfwhereline]
            width = velline[imax] - velline[imin] 
            #- This leads to a gaussian sigma of
            sigma = width/(2.*sqrt(2.*log(2.)))
            #- velmid is where the irregularity reaches its most extreme point
            velmid = (velline[imax] + velline[imin])/2.
            #- The depth of the irregularity is roughly equal to the difference
            #  in max and min flux in an interval between imin-1/2*width and
            #  imax+1/2*width, ie between velmid-width and velmid+width 
            interval = np.abs(velfdf-velmid) <= width
            fmax = max(fluxfdf[interval])
            fmin = min(fluxfdf[interval])
            depth = -abs(fmax - fmin)
            #-- Return half the sigma, to make sure the gauss fit doesnt try 
            #   too wide gaussians (eg whya co10 IRAM)
            return (depth,velmid,sigma/2.,0)
        else:
            #- the found change in slope is likely due to the line itself
            return None
    else:
        return None
    


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
    
    @return: (vexp,evexp,function,bool) The soft parabola or gaussian fit to the 
             line, excluding the extra gaussian component if it has been added.
             The vexp and its error are returned as well, to avoid confusion 
             if the function is a gaussian. Lastly, a bool is returned indicating
             if an extra gauss was needed in addition to the soft parabola.
    @rtype: tuple(float,float,funclib.Function(),bool)
    
    '''
    
    print '***********************'
    if filename <> None: print '** Fitting line in %s.'%filename
    if lprof is None:
        lprof = readLineProfile(filename,info_path)
    vel = lprof.getVelocity()
    flux = lprof.getFlux()
    vel = vel[-np.isnan(flux)]
    flux = flux[-np.isnan(flux)]
    vlsr = lprof.getVlsr()
    
    #-- First guess, only relevant for initial window selection. 
    #   vexp is narrowed down later.
    #vexp = 80
    #-- select line
    #keep = np.abs(vel-vlsr)<=(vexp)
    #velsel,fluxsel = vel[keep],flux[keep]
    #-- Initial values: [peak tmb,vlsr,vexp,gamma] 
    #   For the central peak value, get a first guess from the data
    #   Attempt multiple vexp values and return the best fitting case. 
    #   The initial values are given, with an arbitrary value for the vexp key
    peak = getPeakLPData(lprof=lprof)
    ivexp = array([50.,40.,30.,25.,20.,15.,10.])
    #-- varyInitialFit adapts the velocity window itself. No more assumptions 
    #   needed for the expansion velocity
    firstguess = varyInitialFit(vel,flux,[peak,vlsr,0.,1.],index=2,\
                                values=ivexp,function=funclib.soft_parabola,\
                                vary_window=1)
    vexp = abs(firstguess.get_parameters()[0][2])
    window = 2.
    print 'First guess fit, using a soft parabola:'
    print firstguess.param2str(accuracy=5)
    
    #-- If vexp > 100, replace with 50. This value is unrealistic, and might be
    #   improved with an extra Gaussian. If not, it will be replaced with a 
    #   full Gaussian fit anyway
    if vexp > 100: vexp = 50.
    #-- Check whether irregularities are present in the profile. 
    #   Initial parameters for a gaussian are returned if true. 
    #   For this, a broad selection is made of the profile, to allow for a 
    #   decent noise determination outside the line profile
    keep = np.abs(vel-vlsr)<=(2*window*vexp)
    velsel,fluxsel = vel[keep],flux[keep]
    include_gauss = checkLPShape(velsel,fluxsel,vlsr,vexp,window,show=show)
    
    #-- Do the fit of the line again, including an extra gaussian if 
    #   irregularities are present. 
    if include_gauss <> None:
        #-- fit soft para model + gaussian
        #   1. Set up new soft parabola for several guesses of vexp
        ivexp = list(ivexp)
        initial = [peak,vlsr,0.,1.]
        all_init = [[p]*len(ivexp) for i,p in enumerate(initial) if i != 2]
        all_init.insert(2,ivexp)
        functions = [funclib.soft_parabola() for i in ivexp]
        [ff.setup_parameters(values=initi) 
         for ff,initi in zip(functions,zip(*all_init))]
        #   2. setup gaussian
        gaussians = [funclib.gauss() for ff in functions]
        #   initial guesses assuming an interstellar absorption line from the
        #   checkLPShape method
        [gg.setup_parameters(values=include_gauss) for gg in gaussians]
        #   3. combine soft para + gaussian, and minimize fit
        mymodels = [fit.Model(functions=[ff,gg]) 
                    for ff,gg in zip(functions,gaussians)]
        [fit.minimize(vel[np.abs(vel-vlsr)<=(init[2]*1.5)],\
                      flux[np.abs(vel-vlsr)<=(init[2]*1.5)],\
                      mymodel) 
         for mymodel,init in zip(mymodels,zip(*all_init))]
        #   4. Select the best fitting result based on the error on vexp
        mymodels = [fg 
                   for fg in mymodels
                   if fg.get_parameters()[1][2] != 0.]
        functions = [ff 
                     for fg,ff in zip(mymodels,functions)
                     if fg.get_parameters()[1][2] != 0.]
        gaussians = [gg
                     for fg,gg in zip(mymodels,gaussians)
                     if fg.get_parameters()[1][2] != 0.]
        fitvalues = array([fg.get_parameters()[0][2] 
                           for fg in mymodels])
        fiterrors = array([fg.get_parameters()[1][2] 
                           for fg in mymodels])
        mymodel = mymodels[argmin(abs(fiterrors/fitvalues))]
        finalfit = functions[argmin(abs(fiterrors/fitvalues))]
        gaussian = gaussians[argmin(abs(fiterrors/fitvalues))]
        print 'Improved fit, including extra Gaussian:'
        print mymodel.param2str(accuracy=5)
    else: 
        #-- firstguess is best we can do at the moment
        finalfit = firstguess
        
    #-- If the relative error on vexp is larger than 30%, usually something 
    #   funky is going on in the emission line. Try a Gaussian instead.
    vexp = abs(finalfit.get_parameters()[0][2])
    evexp = abs(finalfit.get_parameters()[1][2])
    gamma = finalfit.get_parameters()[0][3]
    #-- Gamma has to be positive. If it isnt, dont bother with Gaussian
    #   (double peaked line profile will not be fitted well with a Gaussian!)
    if (evexp/vexp > 0.40 and gamma > 0) or (evexp/vexp > 0.20 and vexp> 30.):
        #-- Go back to default window to try a Gaussian fit
        #keep = np.abs(vel-vlsr)<=(80)
        #velselg,fluxselg = vel[keep],flux[keep]
        include_gauss = None
        #-- FWHM is twice vexp!
        sigmas = 2*ivexp/(2.*sqrt(2.*log(2.)))
        finalfit = varyInitialFit(vel,flux,[peak,vlsr,0.,0.],index=2,\
                                  values=sigmas,function=funclib.gauss,\
                                  vary_window=1)
        vexp = abs(finalfit.get_parameters()[0][2])*(2.*sqrt(2.*log(2.)))/2.
        evexp = abs(finalfit.get_parameters()[1][2])*(2.*sqrt(2.*log(2.)))/2.
        window = 3.
        print 'Improved fit, using a gaussian instead of soft parabola:'
        print finalfit.param2str(accuracy=5)
        
    #-- Compute numerical integrations.
    #   After fitting, window for integration should be 0.6*window. vexp is
    #   not expected to be too small anymore as in checkLPShape
    keep = np.abs(vel-vlsr)<=(0.6*window*vexp)
    velsel = vel[keep]
    flux_first = firstguess.evaluate(velsel)
    flux_fg = finalfit.evaluate(velsel)
    print('I_mb (emission line data): %f'\
          %trapz(y=flux[keep],x=velsel))
    print('I_mb (SP -- initial guess): %f'\
          %trapz(y=flux_first,x=velsel))
    print('I_mb (SP -- improved guess): %f'\
          %trapz(y=flux_fg,x=velsel))
    if include_gauss <> None:
        fitted_flux = mymodel.evaluate(velsel)
        print('I_mb (SP + Gauss fit): %f'\
              %trapz(y=fitted_flux,x=velsel))
    print('Final v_exp guess: %.2f +/- %.2f km/s'%(vexp,evexp))
    
    #-- plot
    if show:
        plt.clf()    
        #-- improve velocity window for plotting
        keep = np.abs(vel-vlsr)<=(1.5*window*vexp)
        velsel,fluxsel = vel[keep],flux[keep]
        vel_highres = np.linspace(velsel[0],velsel[-1],10000)
        flux_fg_highres = finalfit.evaluate(vel_highres)
        flux_first_highres = firstguess.evaluate(vel_highres)    
        plt.step(velsel,fluxsel,'-r',where='mid',lw=3,\
                 label='Observed profile')
        plt.plot(vel_highres,flux_first_highres,'b-',lw=3,\
                 label='First guess')
        plt.plot(vel_highres,flux_fg_highres,'g--',lw=3,\
                 label='Improved guess')
        if include_gauss <> None:
            flux_full_highres = mymodel.evaluate(vel_highres)
            plt.plot(vel_highres,flux_full_highres,'g-',lw=2,\
                     label='Full fit (including Gaussian)')
        leg = plt.legend(loc='best',fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.show()
    return vexp,evexp,finalfit,not include_gauss is None