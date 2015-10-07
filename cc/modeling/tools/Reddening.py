# -*- coding: utf-8 -*-

"""
Toolbox for interstellar Reddening of models or dereddening of data, and 
interstellar extinction laws.

Author: R. Lombaert

How does reddening in ComboCode work? 

First and foremost, we define the interstellar extinction in the Johnson K band.
The magnitude of extinction, Ak, is derived from extinction models. Either from
Marshall et al. 2006, or from Drimmel et al. 2003 if not available from the 
former. Ak requires the distance to be given, before it can be applied to a 
model. All of our modeling requires a distance as input, so this is not an issue
in practice.

The interstellar extinction law of preference is that of Chiar and Tielens 2006.
This law is given normalized per Ak and can be directly combined with the 
interstellar extinction given from Marshall or Drimmel. We use the curve for the
local ISM. The alternative is a curve for the galactic center, but even in the 
direction of the galactic center the local ISM doesn't change much in terms of 
dust extinction, except at very large distances on the order of 5kpc or more. We
don't work with sources at those distances for now, so we can safely ignore it.
For completeness, the GC curve is made as well and provided as an option in the
reddening module of IvS repo.

However, while Marshall gives Ak and presents no issue, Drimmel gives Av. To
convert Av to Ak, we have to convert the V-band normalization of Drimmel to 
K-band normalization. Chiar and Tielens, however, derived a law only in the IR
hence no V-band normalization can be defined. We need a different interstellar 
reddening law in V-band to be compared with the infrared law of the
former. The most recent V-band reddening law is given by Fitzpatrick 2004.

We therefore created our own interstellar reddening law from the combination of 
Fitzpatrick 2004 up to the wavelength where Chiar and Tielens 2006 begins. They 
match almost identically in the overlapping region, following a power law of the 
form lambda**-1.8. From there, the combined law follows Chiar and Tielens, and 
is extrapolated to further wavelengths with the same power law with power -1.8 
as mentioned by Chiar & Tielens 2006 between 2 and 5 micron. At long 
wavelengths, the extinction becomes negligible, so the extrapolation is barely 
noticeable, but maintains consistency. 

To convert Fitzpatrick 2004 to Ak, we do need to assume a Av to Ak conversion
that does not take into account Chiar and Tielens. The latter suggest their law
can be converted back to Av with a factor of ak/av = 0.09, which is in very
good agreement with the factor derived from fitzpatrick 2004 itself: 0.088.

Using this self-constructed reddening law, we can now convert Av to Ak from 
Drimmel, and then apply that Ak together with the self-constructed reddening law
to redden our models. We use the IvS repository for the reddening.

"""

import os
from scipy import hstack, array
import numpy as np

import cc.path
from cc.tools.io import DataIO

import ivs.sed.reddening as red
import ivs.sed.extinctionmodels as em


def getAk(ll,bb,distance=None,law='Fitz2004Chiar2006',lawtype='ism'):

    '''
    Find the Johnson K-band interstellar extinction at given longitude and 
    latitude.
    
    @param ll: The galactic longitude of the star
    @type ll: float
    @param bb: The galactic latitude of the star
    @type bb: float
    
    @keyword distance: Distance to the star. Default is None, in which case the
                       full extinction to infinity in a given direction is 
                       returned
                       
                       (default: None)
    @type distance: float    
    @keyword law: The reddening law
                
                  (default: 'Fitz2004Chiar2006')
    @type law: str
    @keyword lawtype: The type of Chiar & Tielens reddening law (either ism or 
                      gc)
                      
                      (default: 'ism')
    @type lawtype: str
    
    @return: The interstellar extinction magnitude in K-band 
    @rtype: float 
    
    '''
        
    ak = em.findext_marshall(ll=ll,bb=bb,distance=distance,redlaw=law,\
                             curve=lawtype,norm='Ak')
    if not ak:
        ak = em.findext_drimmel(lng=ll,lat=bb,distance=distance,redlaw=law,\
                                curve=lawtype,norm='Ak')[0]
    return ak



def redden(wave,flux,ak,law='Fitz2004Chiar2006',lawtype='ism'):
    
    '''
    Redden model fluxes, correcting for interstellar extinction. 
    
    Flux is assumed to be flux, and not magnitudes!
    
    For dereddening, pass -ak instead.
    
    The reddening law can be chosen, but should probably be Fitz2004Chiar2006 as
    it was tailored to infrared reddening of AGB sources in the Solar 
    neighbourhood.
    
    @param wave: The wavelength grid
    @type wave: array
    @param flux: The flux from the models
    @type flux: array
    @param ak: The interstellar reddening magnitude in Johnson K-band
    @type ak: float 
    
    @keyword law: The reddening law
                
                  (default: 'Fitz2004Chiar2006')
    @type law: str
    @keyword lawtype: The type of Chiar & Tielens reddening law (either ism or 
                      gc)
                      
                      (default: 'ism')
    @type lawtype: str
    
    @return: The reddened fluxes
    @rtype: array
     
    '''
    
    wave,a_ak = red.get_law(name=law,wave=wave,curve=lawtype,\
                            norm='Ak',wave_units='micron')
    return flux / 10**(a_ak*ak/2.5)
    
    
    
def combineRedLaw(ofn,chiar_curve='ism',power=-1.8):

    '''
    A method to combine the Fitzpatrick 2004 and Chiar & Tielens 2006 reddening
    laws as well as to extrapolate Chiar and Tielens 2006 to longer wavelengths.
    
    The result is saved in a file and used by the IvS repository as a valid 
    reddening law. 
    
    @param ofn: The output filename with path
    @type ofn: str
    
    @keyword chiar_curve: The curve type for Chiar & Tielens 2004. Either 'gc' 
                          or 'ism'.
                          
                          (default: 'ism')
    @type chiar_curve: str
    @keyword power: The power for the power law extrapolation. Default is taken
                    from Chiar and Tielens 2006, as a typical value for local
                    ISM between 2 and 5 micron. gc may require different value
                    but not very important.
                    
                    (default: -1.8)
    @type power: float

    '''
    
    chiar_curve = chiar_curve.lower()
    
    #-- Extract the two relevant extinction laws.
    xchiar, a_ak_chiar = red.get_law('chiar2006',norm='Ak',wave_units='micron',\
                                     curve=chiar_curve)
    xfitz, a_ak_fitz = red.get_law('fitzpatrick2004',norm='Ak',\
                                   wave_units='micron')
    
    #-- Define a power law for the extrapolation
    def power_law(x,scale,power): return scale*(x)**power
        
    #-- Determine the scaling factor from specific chiar/tielens law
    scale = a_ak_chiar[-1]/(xchiar[-1]**power)

    #-- Create an x grid for longer wavelengths.
    xlong = np.linspace(xchiar[-1]+0.1,1000,1000)
    a_ak_long = power_law(xlong,scale,power)

    #-- Combine the three sections
    xcom = hstack([xfitz[xfitz<xchiar[0]],xchiar,xlong])
    a_ak_com = hstack([a_ak_fitz[xfitz<xchiar[0]],a_ak_chiar,a_ak_long])

    #-- Write the result to a file
    comments = '#-- wavelength (micron)   A_lambda/A_k\n'
    DataIO.writeCols(filename=ofn,cols=[[comments]])
    DataIO.writeCols(filename=ofn,cols=[xcom,a_ak_com],mode='a')