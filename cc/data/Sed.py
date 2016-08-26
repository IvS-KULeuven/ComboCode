# -*- coding: utf-8 -*-

"""
Toolbox for reading spectra and spectrum manipulation.

Author: R. Lombaert

"""

import os
from numpy import argsort,array
from scipy.integrate import trapz
from scipy.interpolate import interp1d
from numpy.core.defchararray import rfind,ljust
import numpy.lib.recfunctions as recfunc
import numpy as np
from glob import glob
import operator

from astropy import units as u

from cc.ivs.sed import builder, filters
import cc.ivs.sed.reddening as ivs_red
from cc.ivs.sed.model import synthetic_flux
#from ivs.units import conversions

import cc.path
from cc.tools.io import DataIO
from cc.modeling.tools import Reddening
from cc.modeling.codes import MCMax


def getCFlux(wav,seds=[],star_grid=[],nans=1,deredden=[],\
             law='Fitz2004Chiar2006',lawtype='ism',map='marshall'):
    
    ''' 
    Retrieve the continuum flux at a given wavelength from either a model
    spectrum or an observation. If both seds and star_grid are given, 
    values from the seds are returned first in the array, then the models.
    
    For now assumes the observation is either an ISO SWS spectrum OR that the 
    continuum point is given in a dictionary that is property of the Sed() 
    object (sed.cflux) with wavelengths as keys. 
    
    star_grid are all models! Works also when seds or star_grid are empty.
    
    Reddening is taken into account when requested in the models and parameters
    are taken from the model objects. However, this is only allowed if only one
    data object is given (otherwise model reddening doesn't make sense)
    
    Dereddening of data is also possible (and extra arguments can be passed to 
    the reddening law), in which case distances have to be given for the seds. 
    If any model reddening is requested and only sed is given, sed dereddening 
    is always turned off.
    
    @param wav: The continuum wavelength point
    @type wav: float
    
    @keyword seds: The SEDs of the data objects. Number of SEDs sets the amount
                   of Star() objects represent data versus number of models.
                 
                   (default: [])
    @type seds: list(Sed())
    @keyword star_grid: The data + model objects
                        
                        (default: [])
    @type star_grid: list(Star())
    @keyword nans: Set undefined line strengths as nans. Errors are set as a 
                   nan if it concerns mode==dint. Otherwise, they are not set.
                   
                   (default: 1)
    @type nans: bool
    @keyword deredden: Deredden the SEDs with distances given here. This option 
                       is turned off automatically if any reddening is requested
                       in the models and only one sed is given to avoid double 
                       correction. Number of distances given must be equal to 
                       number of SEDs given. If not, it is also turned off. 
                       
                       (default: []) 
    @type deredden: list
    @keyword law: The reddening law for DEREDDENING
                
                  (default: 'Fitz2004Chiar2006')
    @type law: str
    @keyword lawtype: The type of Chiar & Tielens reddening law (either ism or 
                      gc) for DEREDDENING
                      
                      (default: 'ism')
    @type lawtype: str
    @keyword map: The galactic 3d extinction model for DEREDDENING. 
    
                      (default: 'marshall')
    @type map: str
    
    @return: The continuum fluxes in W/m2/Hz with length that of star_grid, as 
             well as errors if applicable. If a combo mode is requested, errors 
             are given when available, and listed as None/nan if not available 
             (Plotting2 module knows how to deal with this).
    @rtype: (array[float],array[float])
    
    '''
    
    all_cflux = []
    all_eflux = []
    rlaw = []
    
    dtype = ''
    if seds:
        ddict = dict([('SWS',(2.4,45.0)),('PACS',(55.1,189.))])
        for k,v in ddict.items():
            if v[0] <= wav and wav <= v[1]:
                dtype = k
                break
    
    #-- Only allow model reddening if one sed is given. Multiple: makes no sense
    #   to redden models. None: No galactic coordinates given to redden.
    redden = [s['REDDENING'] for s in star_grid] if len(seds) == 1 else []
    
    #-- Only allow dereddening if enough distances are given, and if model 
    #   reddening for one sed is not requested.
    if len(deredden) != len(seds) or np.any(redden):
        deredden = []
    
    #-- Interpolate the reddening law once.
    if deredden or redden: 
        wave_arr,rlaw = ivs_red.get_law(name=law,wave=wav,curve=lawtype,\
                                        norm='Ak',wave_units='micron')
        
    #-- First all data objects
    for ised,sed in enumerate(seds):
        dtypes = []
        if dtype:
            dtypes = [dt for dt in sed.data.keys() if dtype in dt[0].upper()]
        #-- If no SWS spectrum find, check if the flux is available in sed.flux
        if not dtypes:
            if sed.cflux.has_key(wav): 
                tflux = sed.cflux[wav]
                if deredden:
                    ak = sed.getAk(deredden[ised],map=map,law=law)
                    #-- deredden so increase flux, as opposed to redden
                    tflux = tflux * 10**(rlaw*ak/2.5)
                all_cflux.append(tflux)
                all_eflux.append(sed.eflux[wav])
            else:
                all_cflux.append(nans and float('nan') or None)
                all_eflux.append(nans and float('nan') or None)
            continue
        #-- At least one dtype spectrum found, take the first one.
        dt = dtypes[0]
        abs_err = sed.abs_err[dt[0]]
        dwave = sed.data[dt][0]
        dflux = sed.data[dt][1]
        
        #-- Check if the data object gives the standard deviation
        if len(sed.data[dt]) > 2:
            i = np.argmin(abs(dwave-wav))
            ilow = i if wav>dwave[i] else i-1
            iup = i if wav<dwave[i] else i+1
            errs = sed.data[dt][2]
            deflux = np.sqrt((errs/dflux)[ilow]**2+(errs/dflux)[iup]**2)
        else:
            deflux = 0.0
        
        #-- Interpolate for flux, and set the error taking into account abs flux
        #   calib uncert.
        interp = interp1d(dwave,dflux)
        tflux = interp(wav)
        if deredden:
            ak = sed.getAk(deredden[ised],map=map,law=law)
            #-- deredden so increase flux
            tflux = tflux * 10**(rlaw*ak/2.5)
        all_cflux.append(tflux)
        all_eflux.append(np.sqrt(deflux**2+abs_err**2))
        
    #-- Then all model objects
    all_eflux.extend([nans and float('nan') or None]*len(star_grid))
    for s in star_grid:
        if not s['LAST_MCMAX_MODEL']:
            all_cflux.append(nans and float('nan') or None)
            continue
        cc.path.mout = os.path.join(cc.path.mcmax,s.path_mcmax)
        dpath = os.path.join(cc.path.mcmax,s.path_mcmax,'models',\
                             s['LAST_MCMAX_MODEL'])
        w,f = MCMax.readModelSpectrum(dpath,rt_sed=1)
        interp = interp1d(w,f)
        tflux = interp(wav)
        if s['REDDENING'] and seds:
            #-- Only one SED is supposed to be given.
            ak = seds[0].getAk(s['DISTANCE'],map=s['REDDENING_MAP'],\
                               law=s['REDDENING_LAW'])
            #-- redden so decrease flux
            tflux = tflux / 10**(rlaw*ak/2.5)
        all_cflux.append(tflux)
    
    #-- All fluxes for SED type data or models are given in Jy. Convert to 
    #   W/m2/Hz. Errors are given in relative numbers, so no conversion needed.
    all_cflux = array(all_cflux)*1e-26
    all_eflux = array(all_eflux)
    
    return (all_cflux,all_eflux)



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
    
    #-- Convert wavelength from micron to angstrom. astropy conversion module 
    #   returns a "Quantity" that has properties "unit" and "value"
    mlam = (w*u.micron).to(u.AA).value
    
    #-- Convert Jy flux density to erg/s/cm2/aa for windowed integration
    #   Requires equivalency because conversion depends on the wavelength where
    #   the flux is measured. Uses spectral_density equivalency.
    mflam = (f*u.Jy).to(u.erg/u.s/u.cm**2/u.AA,\
                        equivalencies=u.spectral_density(w*u.micron)).value

#    mflam = conversions.convert('Jy','erg/s/cm2/AA',f,wave=(w,'micron'))
#    mlam = conversions.convert('micron','AA',w)
    mphot = synthetic_flux(mlam,mflam,photbands,units=['Fnu']*len(photbands))
    
    mphot = (mphot*u.erg/u.s/u.Hz/u.cm**2).to(u.Jy).value
#    mphot = conversions.convert('erg/s/cm2/Hz','Jy',mphot)
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
        self.photbands = np.empty(0,dtype=str)
        self.photwave = np.empty(0)
        self.star_name = star_name
        self.data = dict()
        self.setStarPars()
        self.setData(remove=remove)
        self.readData()
        self.readPhotInfo()
        self.ak = dict()
        
        #-- Extra property defining continuum flux points that can be set
        #   externally. Assumed to give (wavelength,flux) as (key,value) pair
        #   Used by, e.g., the Sed.getCFlux() method. An error must also be 
        #   given (in relative number).
        self.cflux = dict()
        self.eflux = dict()        
        
        
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
        abs_errs = DataIO.getInputData(keyword='ABS_ERR',filename='Sed.dat')
        
        if 'Photometric_IvS' in data_types:
            buildPhotometry(self.star_name,**kwargs)
        
        self.data_types = []
        self.data_filenames = []
        self.abs_err = dict()
        for dt,ierr in zip(data_types,abs_errs):
            searchpath = os.path.join(cc.path.dsed,'%s_*%s*.dat'\
                                                   %(dt,self.star_name))
            add_files = glob(searchpath)
            for ff in add_files: 
                if ff not in self.data_filenames:
                    self.data_filenames.append(ff)
                    self.data_types.append(dt)
                    self.abs_err[dt] = ierr



    def readData(self):
        
        '''
        Read the raw SED data. 
        
        '''
        
        for dt,fn in zip(self.data_types,self.data_filenames):
            print dt, fn
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
                if  dt == 'Photometric_IvS':
                    self.photbands = data[3]
                    self.photwave = data[0]
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
        wlower = [sel[0] for sel in selection]
        wupper = [sel[-1] for sel in selection]
        self.filter_info = recfunc.append_fields(self.filter_info,\
                                                 ['wlower','wupper'],\
                                                 [wlower,wupper],usemask=0,\
                                                 asrecarray=1)
    
    def getAk(self,distance=None,map='marshall',law='Fitz2004Chiar2006'):
    
        '''
        Helper method to retrieve the Ak extinction magnitude for a given 
        distance.
        
        The Ak values are saved in the sed object, to cut down the overhead in 
        subsequent calls. 
        
        The law and map have to be passed as well. Defaults are marshall and 
        Fitz2004Chiar2006 respectively.
        
        @param distance: The distance to the source. If the default, the total
                         extinction is calculated in given direction. 
                         
                         (default: None)
        @type distance: float
        @keyword map: The galactic 3d extinction model. 
    
                      (default: 'marshall')
        @type map: str
        @keyword law: The reddening law
                
                      (default: 'Fitz2004Chiar2006')
        @type law: str
    
        @return: The extinction magnitude in K-band for the requested distance.
        @rtype: float
        
        '''
        
        distance = float(distance)
        if not self.ak.has_key((distance,map,law)):
            aki = Reddening.getAk(self.ll,self.bb,distance,map,law)
            self.ak[(distance,map,law)] = aki
        
        return self.ak[(distance,map,law)]