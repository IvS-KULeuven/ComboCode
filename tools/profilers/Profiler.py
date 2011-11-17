# -*- coding: utf-8 -*-

"""
Tools for making profiles of any kind.

Author: R. Lombaert

"""

import scipy
import os 

from cc.modeling.objects import Star
from cc.tools.io import DataIO



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
    filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                            path_gastronoom,'models',model_id,\
                            'coolfgr_all%s.dat'%model_id)
    rad = Gastronoom.getGastronoomOutput(filename=filename,keyword='RADIUS',\
                                         return_array=1)
    fraction_profile = scipy.ones(len(rad))
    step_index = scipy.argmin(abs(rad-rfrac))
    fraction_profile[step_index:] = fraction
    output_filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                   path_gastronoom,'profiles',\
                                   'water_fractions_%s_%.2f_r%.3e.dat'\
                                   %(model_id,fraction,rfrac))
    DataIO.writeCols(output_filename,[rad,fraction_profile])
    