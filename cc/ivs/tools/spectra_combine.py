# -*- coding: utf-8 -*-

"""
The combine method from the ivs.spectra.tools module.

"""

import logging
import numpy as np
import scipy.stats
from cc.ivs.units import conversions

logger = logging.getLogger("SPEC.TOOLS")

def combine(list_of_spectra,R=200.,lambda0=(950.,'AA'),lambdan=(3350.,'AA')):
    """
    Combine and weight-average spectra on a common wavelength grid.
    
    C{list_of_spectra} should be a list of lists/arrays. Each element in the
    main list should be (wavelength,flux,error).
    
    If you have FUSE fits files, use L{cc.ivs.fits.read_fuse}.
    If you have IUE FITS files, use L{cc.ivs.fits.read_iue}.
    
    After Peter Woitke.
    
    @param R: resolution
    @type R: float
    @param lambda0: start wavelength, unit
    @type lambda0: tuple (float,str)
    @param lambdan: end wavelength, unit
    @type lambdan: tuple (float,str)
    @return: binned spectrum (wavelengths,flux, error)
    @rtype: array, array, array
    """
    l0 = conversions.convert(lambda0[1],'AA',lambda0[0])
    ln = conversions.convert(lambdan[1],'AA',lambdan[0])
    #-- STEP 1: define wavelength bins
    Delta = np.log10(1.+1./R)
    x = np.arange(np.log10(l0),np.log10(ln)+Delta,Delta)
    x = 10**x
    lamc_j = 0.5*(np.roll(x,1)+x)

    #-- STEP 2: rebinning of data onto newly defined wavelength bins
    Ns = len(list_of_spectra)
    Nw = len(lamc_j)-1
    binned_fluxes = np.zeros((Ns,Nw))
    binned_errors = np.inf*np.ones((Ns,Nw))

    for snr,(wave,flux,err) in enumerate(list_of_spectra):
        wave0 = np.roll(wave,1)
        wave1 = np.roll(wave,-1)
        lam_i0_dc = 0.5*(wave0+wave)
        lam_i1_dc = 0.5*(wave1+wave)
        dlam_i = lam_i1_dc-lam_i0_dc
        
        for j in range(Nw):
            A = np.min(np.vstack([lamc_j[j+1]*np.ones(len(wave)),lam_i1_dc]),axis=0)
            B = np.max(np.vstack([lamc_j[j]*np.ones(len(wave)),lam_i0_dc]),axis=0)
            overlaps = scipy.stats.threshold(A-B,threshmin=0)
            norm = np.sum(overlaps)
            binned_fluxes[snr,j] = np.sum(flux*overlaps)/norm
            binned_errors[snr,j] = np.sqrt(np.sum((err*overlaps)**2))/norm
    
    #-- STEP 3: all available spectra sets are co-added, using the inverse
    #   square of the bin uncertainty as weight
    binned_fluxes[np.isnan(binned_fluxes)] = 0
    binned_errors[np.isnan(binned_errors)] = 1e300
    weights = 1./binned_errors**2
    
    totalflux = np.sum(weights*binned_fluxes,axis=0)/np.sum(weights,axis=0)
    totalerr = np.sqrt(np.sum((weights*binned_errors)**2,axis=0))/np.sum(weights,axis=0)
    totalspec = np.sum(binned_fluxes>0,axis=0)
    
    #-- that's it!
    return x[:-1],totalflux,totalerr,totalspec

