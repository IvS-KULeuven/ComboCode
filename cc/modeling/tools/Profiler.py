# -*- coding: utf-8 -*-

"""
Tools for making profiles of any kind.

Author: R. Lombaert

"""

from scipy import array 
import scipy
import os 

import cc.path
from cc.tools.io import DataIO


def dustTemperaturePowerLaw(rad,tstar,s=1,add_key=0):
    
    '''
    Return a dust temperature power law of the form as suggested by 
    observational evidence. Does not depend on MCMax output.
    
    The radial grid is given in Rstar!
            
    See Thesis p32, where power is -2/(4+s) in 
    T(r) = T_eff*(2*r/R_STAR)**(-2/(4+s))
    and s is typically 1 in the Rayleigh limit and optically thin case.
    
    @param rad: The radial grid for the power law in Rstar
    @type rad: array
    @param tstar: The stellar effective temperature
    @type tstar: float
    
    @keyword s: The s parameter in the power law T(r) given above.
    
                (default: 1)
    @type s: int
    @keyword add_key: Add a key for a legend to the ouput as third tuple
                      element.
                        
                      (default: 0)
    @type add_key: bool

    @return: The temperature profile (K) as well as a key, if requested.
    @rtype: (array,string)
    
    '''
    
    s, tstar, rad = int(s), float(tstar), array(rad)
    temp = tstar*(2.*rad)**(-2./(4+s))
        
    if add_key:
        key = '$T_\mathrm{d}(r) = %i \\left(\\frac{2r}'%int(tstar) + \
              '{\mathrm{R}_\star}\\right)^{\\frac{2}{4+%i}}$'%int(s)
        return temp, key
    else: 
        return temp
    


def waterFraction1StepProfiler(model_id,path_gastronoom,fraction,rfrac):

    '''
    Create a 1-step fractional profile for water.
    
    The original water abundance profile is taken from the output of the 
    original model without fractional abundances. 
    
    These fraction profiles can be used for CHANGE_ABUNDANCE_FRACTION in mline
    
    @param model_id: The model id of the original cooling model
    @type model_id: string
    @param path_gastronoom: The model subfolder in ~/GASTRoNOoM/
    @type path_gastronoom: string
    @param fraction: the fraction used
    @type fraction: float
    @param rfrac: the radius at the step to the fractional abundance [cm]
    @type rfrac: float
    
    '''
    
    rfrac = float(rfrac)
    fraction = float(fraction)
    filename = os.path.join(cc.path.gastronoom,path_gastronoom,'models',\
                            model_id,'coolfgr_all%s.dat'%model_id)
    rad = Gastronoom.getGastronoomOutput(filename=filename,keyword='RADIUS',\
                                         return_array=1)
    fraction_profile = scipy.ones(len(rad))
    step_index = scipy.argmin(abs(rad-rfrac))
    fraction_profile[step_index:] = fraction
    output_filename = os.path.join(cc.path.gastronoom,path_gastronoom,\
                                   'profiles',\
                                   'water_fractions_%s_%.2f_r%.3e.dat'\
                                   %(model_id,fraction,rfrac))
    DataIO.writeCols(output_filename,[rad,fraction_profile])
    