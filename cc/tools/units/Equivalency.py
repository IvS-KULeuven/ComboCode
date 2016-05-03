# -*- coding: utf-8 -*-

"""
Equivalencies for the astropy.units module for converting several quantities in
both directions.

Author: R. Lombaert

"""

from astropy import units as u
from astropy import constants as cst
import numpy as np


def Tmb(diameter):

    '''
    Converting main-beam temperature in K to Fnu in cgs and vice versa.
    
    Requires the diameter to be given in m.
     
    >>> equiv = Tmb(12.)
    >>> Fnu = (1*u.K).to(u.W/u.m/u.m/u.Hz,equivalencies=equiv)
    gives the length in astronomical units.
    
    Note that you don't have to give arcsec or au. You can also ask for u.cm, 
    u.arcminute, etc. The astropy unit conversion module takes care of these 
    conversions as long as no equivalency is needed.
    
    @param distance: The distance to the object in parsec
    @type distance: float
    
    @return: list of (unit in, unit out, forward conversion, reverse conversion)
    @rtype: list(tuple)
    
    '''
    
    
    d = float(diameter)*u.m
    factor = np.pi*d**2/8./cst.k_B
    return [(u.K,u.W/u.m/u.m/u.Hz, lambda x: x/factor, lambda x: x*factor)]



def angularSize(distance): 

    '''
    Equivalency for the astropy.units module for converting angular size to 
    physical length and vice versa. 
    
    Requires the distance to be given in pc.
    
    Can be used to convert angular sizes and lengths using the units module. 
    >>> equiv = angularSize(1000.)
    >>> rad = (0.001*u.arcsec).to(u.au,equivalencies=equiv)
    gives the length in astronomical units.
    
    Note that you don't have to give arcsec or au. You can also ask for u.cm, 
    u.arcminute, etc. The astropy unit conversion module takes care of these 
    conversions as long as no equivalency is needed.
    
    @param distance: The distance to the object in parsec
    @type distance: float
    
    @return: list of (unit in, unit out, forward conversion, reverse conversion)
    @rtype: list(tuple)
    
    '''
    
    #- 1 AU at 1 pc is 1 as on the sky, hence 1 AU at x pc is 1/x as on the sky
    #- hence [real_rad/1 AU] = [angrad/1 as] * [distance/1 pc]
    distance = float(distance)
    return [(u.arcsec,u.au, lambda x: x*distance, lambda x: x/distance)]