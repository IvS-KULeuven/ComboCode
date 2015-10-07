# -*- coding: utf-8 -*-

"""
Toolbox for reading spectra and spectrum manipulation.

Author: R. Lombaert

"""

import os
from scipy import argsort
from scipy.integrate import trapz
from numpy.core.defchararray import rfind,ljust
import numpy.lib.recfunctions as recfunc
import numpy as np
from glob import glob
import operator

from ivs.sed import builder, filters
from ivs.sed.model import synthetic_flux
from ivs.units import conversions

import cc.path
from cc.tools.io import DataIO
from cc.modeling.tools import Reddening


def calcPhotometry(w,f,photbands):

    ''' 
    Calculate the (model) photometry for a given set of photometric bands. 
    
    Reddening is assumed to have been done before this.
    
    @param w: the wavelengths in micron
    @type w: array()
    @param f: Flux grid in Jy
    @type f: array()
    @param photbands: the photometric bands
    @type photbands: array(str)
    
    @return: The photometry in Jy
    @rtype: array
    
    '''

    mflam = conversions.convert('Jy','erg/s/cm2/AA',f,wave=(w,'micron'))
    mlam = conversions.convert('micron','AA',w)
    mphot = synthetic_flux(mlam,mflam,photbands,units=['Fnu']*len(photbands))
    mphot = conversions.convert('erg/s/cm2/Hz','Jy',mphot)
    return mphot



def buildPhotometry(star_name,fn='Photometric_IvS',remove=[]):
    '''
    Retrieve the photometry of a star through the IvS repo's SED builder. 
    
    Save the photometry in a target location. 
    
    @param star_name: Name of the star (cc name from Star.dat)
    @type star_name: str
    
    @keyword fn: Output filename of the photometry file. Always appends 
                 '_STAR.dat' where STAR is star_name. The file is saved in dp.
    
                 (default: 'Photometric_IvS')
    @type fn: str
    @keyword remove: Photometry to be removed from the output file. 
                     e.g. ['WISE','DENIS'] 
    
                     (default: [])
    @type remove: list[str]
    
    '''
    
    #-- Get the SIMBAD name of the star
    si = DataIO.getInputData(filename='Star.dat').index(star_name)
    sn_sim = DataIO.getInputData(keyword='STAR_NAME_PLOTS')[si]
    sn_sim = sn_sim.replace('$','').replace('\\','')
    
    #-- Define the outputfolder of the Raw data
    ofn_raw = os.path.join(cc.path.dphot,'_'.join([star_name,'phot.txt']))
    
    ## Retrieve p = builder.SED()
    star = builder.SED(sn_sim,photfile=ofn_raw)
    star.get_photometry()
    
    ## Read in photometric data
    columns = ['wave', 'meas', 'emeas', 'photband', 'bibcode', 'comments']
    data = np.genfromtxt(star.photfile,usecols = (9,10,11,4,15,16),\
                         names=columns,dtype = None)
    
    #-- Find and remove WISE and DENIS data
    for photband in remove:
        selection = rfind(data['photband'],photband)
        data = data[np.where(selection==-1)]
    
    #-- Remove colour measurements (NaN in wavelength)
    data = data[np.where(np.isfinite(data['wave']))]
    
    #-- Fill photbands with whitespace for formatting
    data['photband'] = ljust(data['photband'],15)
    
    #-- Flux: erg/s/cm**2/AA -> Jy (10**-23 erg/s/cm**2/Hz)
    #   From http://astro.wku.edu/strolger/UNITS.txt (non-linear!)
    c = 2.99792458e18       # in AA/s 
    data['meas'] = data['meas']*data['wave']**2/c*1e23
    data['emeas'] = data['emeas']*data['wave']**2/c*1e23
    
    #-- Wavelength: angstrom -> micron
    data['wave'] = data['wave']*10**-4
    
    #-- Sort the data on wavelength
    data.sort()
    
    ofn_final = os.path.join(cc.path.dsed,'_'.join([fn,star_name+'.dat']))
    hdr = 'Photometry extracted with ivs.sed.builder. \nWave (micron) Flux '+\
          '(Jy) Error Flux (Jy) Photband Bibcode Comments'
    np.savetxt(ofn_final,data,fmt=['%.8e']*3+['%s']*3,header=hdr)
    


def normalizeSed(sed):
    
    '''
    Normalize an SED by division by integrated energy.
    
    Input is given by F_nu, integration of the SED is done for nu*F_nu.
    
    @param sed: The input SED wavelenth and F_nu.
    @type sed: (list,list)
    
    '''
    
    c = 2.99792458e14          #in micron/s

    return [sed[0],sed[1]/trapz(x=sed[0],y=sed[1]*(c/sed[0]))]
    
    
    
class Sed(object):
    
    '''
    Tools for:
        - reading data
        - downloading photometry
    
    '''
    
    def __init__(self,star_name,remove=[]):
        
        ''' 
        Initializing an Sed instance. 
        
        Setting starting parameters from a star object.
        
        @param star_name: The star name of the object
        @type star_name: string
        
        @keyword remove: Photometry to be removed from the output file. 
                         e.g. ['WISE','DENIS'] 
    
                         (default: [])
        @type remove: list[str]
        
        '''
        
        self.instrument = 'SED'
        self.photbands = np.empty(0)
        self.star_name = star_name
        self.data = dict()
        self.setStarPars()
        self.setData(remove=remove)
        self.readData()
        self.readPhotInfo()
        self.ak = dict()


    def setStarPars(self):
        
        """
        Set some standard stellar parameters such as Ak and galactic position.
        
        """
        
        self.star_index = DataIO.getInputData().index(self.star_name)
        self.ll = DataIO.getInputData(keyword='LONG',rindex=self.star_index)
        self.bb = DataIO.getInputData(keyword='LAT',rindex=self.star_index)
        snp = DataIO.getInputData(keyword='STAR_NAME_PLOTS',\
                                  remove_underscore=1,rindex=self.star_index)
        self.star_name_plots = snp    
       
       
       
    def setData(self,**kwargs):
        
        '''
        Select available data.
        
        Based on the data file types in Sed.dat and the available data files.
        
        Also calls the buildPhotometry method to create a photometry file from
        the IvS Sed builder tool
        
        Any keywords required for buildPhotometry can be passed here.
        
        '''
        
        data_types = DataIO.getInputData(keyword='DATA_TYPES',\
                                         filename='Sed.dat')
                    
        if 'Photometric_IvS' in data_types:
            buildPhotometry(self.star_name,**kwargs)
        
        self.data_types = []
        self.data_filenames = []
        for dt in data_types:
            searchpath = os.path.join(cc.path.dsed,'%s_*%s*.dat'\
                                                   %(dt,self.star_name))
            add_files = glob(searchpath)
            for ff in add_files: 
                if ff not in self.data_filenames:
                    self.data_filenames.append(ff)
                    self.data_types.append(dt)



    def readData(self):
        
        '''
        Read the raw SED data. 
        
        '''
        
        for dt,fn in zip(self.data_types,self.data_filenames):
            data = DataIO.readCols(fn,nans=0)
            #-- Currently, error bars only available for these types of data.
            if 'Photometric' in dt or 'MIDI' in dt or 'Sacha' in dt: 
                #-- Sort MIDI data
                if 'MIDI' in dt: 
                    cdat = [dd[(data[0]<=13.)*(data[0]>=8.)] for dd in data]
                    i = argsort(cdat[0])
                    self.data[(dt,fn)] = (cdat[0][i],cdat[1][i],cdat[2][i])
                else:
                    self.data[(dt,fn)] = (data[0],data[1],data[2])
                if 'Photometric_IvS' in dt:
                    self.photbands = data[3]
            else:
                #-- Still sorting for PACS. Obsolete when separate bands for 
                #   PACS are available. 
                i = argsort(data[0])
                self.data[(dt,fn)] = (data[0][i],data[1][i])



    def readPhotInfo(self,level=.5):
    
        '''
        Read the photometry band information associated with photometry of this
        SED.
        
        @keyword level: The level at which the cut off for significant 
                        transmission of the photometric bands is placed.
                        
                        (default: 0.5)        
        @type level: float
        
        '''
        
        #-- Get photometry bands info from IvS repo. recarray structure same as 
        #   self.photbands_ivs
        filter_info = filters.get_info()
        keep = np.searchsorted(filter_info['photband'],self.photbands)
        self.filter_info = filter_info[keep]
        self.filter_info.eff_wave = self.filter_info.eff_wave/1e4
        
        response = [filters.get_response(photband) 
                    for photband in self.photbands]
        selection = [waver[transr/max(transr)>level]/1e4
                     for waver,transr in response]
        wlower = [effw-sel[0] 
                  for sel,effw in zip(selection,self.filter_info.eff_wave)]
        wupper = [sel[-1]-effw 
                  for sel,effw in zip(selection,self.filter_info.eff_wave)]
        self.filter_info = recfunc.append_fields(self.filter_info,\
                                                 ['wlower','wupper'],\
                                                 [wlower,wupper],usemask=0,\
                                                 asrecarray=1)
    
    def getAk(self,distance):
    
        '''
        Retrieve the Ak extinction magnitude for a given distance.
        
        The Ak values are saved in the sed object, to cut down the overhead in 
        subsequent calls.
        
        @param distance: The distance to the source. If the default, the total
                         extinction is calculated in given direction. 
        @type distance: float
        
        @return: The extinction magnitude in K-band for the requested distance.
        @rtype: float
        
        '''
        
        distance = float(distance)
        if not self.ak.has_key(distance):
            self.ak[distance] = Reddening.getAk(self.ll,self.bb,distance)
        
        return self.ak[distance]