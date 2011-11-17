# -*- coding: utf-8 -*-

"""
Toolbox for merging opacities.

Author: R. Lombaert

"""

import os
from glob import glob
from scipy import argmin

import DataIO



def mergeOpacity(species,lowres='nom_res',highres='high_res'):
    
    '''
    Merge high-res opacities into a grid of low-res opacities.
    
    The wavelength range of the inserted high res opacities is taken from the 
    given high res grid. 
        
    @param species: The dust species for which this is done. This is also the 
                    name of the folder in ~/MCMax/DustOpacities/ that contains
                    the data files.
    @type species: string
    @keyword lowres: The subfolder in ~/MCMax/DustOpacities/species containing
                     the low resolution datafiles.
                   
                     (default: low_res)
    @type lowres: string
    @keyword highres: The subfolder in ~/MCMax/DustOpacities/species containing
                      the high resolution datafiles.
                      
                      (default: high_res)
    @type highres: string
    
    '''
    
    path = os.path.join(os.path.expanduser('~'),'MCMax','DustOpacities',species)
    lowres_files = [f 
                    for f in glob(os.path.join(path,lowres,'*')) 
                    if f[-5:] == '.opac']
    highres_files = [f 
                     for f in glob(os.path.join(path,highres,'*')) 
                     if f[-5:] == '.opac']
    files = set([os.path.split(f)[1] for f in lowres_files] + \
                [os.path.split(f)[1] for f in highres_files])
    
    for f in files:
        hdfile = os.path.join(path,highres,f)
        ldfile = os.path.join(path,lowres,f)
        if os.path.isfile(ldfile) and os.path.isfile(hdfile):
            hd = DataIO.readCols(hdfile)
            ld = DataIO.readCols(ldfile)
            hdw = hd[0]
            ldw = ld[0]
            wmin = hdw[0]
            wmax = hdw[-1]
            ld_low = [list(col[ldw<wmin]) for col in ld]
            ld_high = [list(col[ldw>wmax]) for col in ld]
            hd = [list(col) for col in hd]
            merged = [ld_low[i] + hd[i] + ld_high[i] 
                      for i in range(len(hd))]
            DataIO.writeCols(filename=os.path.join(path,f),cols=merged)
            