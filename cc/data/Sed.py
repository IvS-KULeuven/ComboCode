# -*- coding: utf-8 -*-

"""
Toolbox for reading spectra and spectrum manipulation.

Author: R. Lombaert

"""

import os
from scipy import array, hstack, argsort
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from numpy.core.defchararray import rfind,ljust
import numpy as np
from glob import glob
import operator

from ivs.sed import extinctionmodels as em 
from ivs.sed import builder

import cc.path
from cc.tools.numerical import Interpol
from cc.tools.io import DataIO
from cc.plotting import Plotting2


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
        - dereddening data
    
    '''
    
    def __init__(self,star_name,distance=None,plot_extrapol_extinction=0,\
                 remove=[]):
        
        ''' 
        Initializing an Sed instance. 
        
        Setting starting parameters from a star object.
        
        @param star_name: The star name of the object
        @type star_name: string
        
        @keyword distance: The distance to the star. If not None, the K-band
                           absorption coefficient is taken from the Marshall
                           or Drimmel catalogs.
                           
                           (default: None)
        @type distance: float
        @keyword remove: Photometry to be removed from the output file. 
                         e.g. ['WISE','DENIS'] 
    
                         (default: [])
        @type remove: list[str]
        @keyword plot_extrapol_extinction: Plot and show the result of the 
                                           extrapolated interstellar extinction
                                           law by chiar and tielens (2006)
                                           
                                           (default: 0)
        @type plot_extrapol_extinction: bool
        
        '''
        
        self.star_name = star_name
        self.distance = distance 
        self.ak = None
        self.plot_extrapol_extinction = plot_extrapol_extinction
        self.data = dict()
        self.data_raw = dict()
        self.setStarPars()
        self.setData(remove=remove)
        self.readData()
        self.dereddenData()



    def setStarPars(self):
        
        """
        Set some standard stellar parameters such as Ak and galactic position.
        
        """
        
        self.star_index = DataIO.getInputData().index(self.star_name)
        ll = DataIO.getInputData(keyword='LONG',rindex=self.star_index)
        bb = DataIO.getInputData(keyword='LAT',rindex=self.star_index)
        if self.distance <> None:
            self.ak = em.findext_marshall(ll=ll,bb=bb,distance=self.distance,\
                                          norm='Ak')
            if self.ak is None:
                self.ak = em.findext_drimmel(lng=ll,lat=bb,norm='Ak',\
                                             distance=self.distance)
        if self.ak is None:
            self.ak = DataIO.getInputData(keyword='A_K',rindex=self.star_index)
        snp = DataIO.getInputData(keyword='STAR_NAME_PLOTS',\
                                  remove_underscore=1,rindex=self.star_index)
        self.star_name_plots = snp    
        if (abs(ll) < 5.0 or ll > 355.0) and abs(bb) < 5.0:
            self.gal_position = 'GC'
        else:
            self.gal_position = 'ISM'    
        
       
       
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
            data = DataIO.readCols(fn,nans=1)
            #-- Currently, error bars only available for these types of data.
            if 'Photometric' in dt or 'MIDI' in dt or 'Sacha' in dt: 
                #-- Sort MIDI data
                if 'MIDI' in dt: 
                    cdat = [dd[(data[0]<=13.)*(data[0]>=8.)] for dd in data]
                    i = argsort(cdat[0])
                    self.data_raw[(dt,fn)] = (cdat[0][i],cdat[1][i],cdat[2][i])
                else:
                    self.data_raw[(dt,fn)] = (data[0],data[1],data[2])
            else:
                #-- Still sorting for PACS. Obsolete when separate bands for 
                #   PACS are available. 
                i = argsort(data[0])
                self.data_raw[(dt,fn)] = (data[0][i],data[1][i])
            
        #print '** LWS is not being aligned automatically with SWS for now.'
        #- Check if SWS and LWS are present (for now, LWS is not aligned automatically)
        # included_sws = [(data_type,i) 
        #                for i,(data_type,boolean) in enumerate(zip(data_types,\
        #                                                           which_data)) 
        #                if boolean and data_type.lower().find('sws') > -1]
        #included_lws = [] #[(data_type,i) for i,(data_type,boolean) in enumerate(zip(data_types,which_data))
        #                #   if boolean and data_type.lower().find('sws') == -1 and data_type.lower().find('lws') > -1]
        #align if necessary (SWS/LWS for now); can become iteration for multiple alignments. alignY is already suited for that
        #
        #if included_sws and included_lws:
            #align_bound_min = 43.0
            #align_bound_max = 45.15
        
            #aligned_y,shifts = Data.alignY([sorted(datagrids[included_sws[0][1]],key=operator.itemgetter(0))]+\
                                #[sorted(datagrids[lwslist[1]],key=operator.itemgetter(0)) for lwslist in included_lws],\
                                #[align_bound_min for lwslist in included_lws],\
                                #[align_bound_max for lwslist in included_lws])
            #print '\n'.join(['The ' + datatype[0] + ' data have been scaled with a scaling factor of ' + str(shift) 
                                #for datatype,shift in zip([included_sws[0]]+included_lws,shifts)])
            #new_datagrids = [i in [datatype[1] for datatype in [included_sws[0]] + included_lws] and aligned_y.pop(0) or data 
                                    #for i,data in enumerate(datagrids)]
        #else:
            #new_datagrids = datagrids
        
       
       
    def dereddenData(self):
        
        '''
        Deredden the data. 
        
        The interstellar extinction curve by Chiar and Tielens (2006) is used 
        with an Av to Ak conversion factor of 0.112.
        
        The correction is done by interpolating the curve, and extrapolating to
        longer wavelengths.
        
        The extrapolation is done with a power law fitted to the long
        wavelength region (lambda > 22 micron).
        
        For short wavelength points below 1.24 micron, the cardelli extinction
        curve is used.
        
        '''
    
        #- Read the extinction curve, then inter/extrapolate it
        ext_x,ext_y = getExtinctionCurve(self.gal_position,'chiar_tielens',\
                                         0.112)
        ext_car_x, ext_car_y = getExtinctionCurve(curve_type='cardelli',\
                                                  av_to_ak_conv=0.112)
        #- Initial param guess for extrapolation of long wavelength extinction curve
        p0 = [ -2,  0.01 ,  1, -1]
        #- Assuming a power law for the extrapolation, and fitting from 22 mic  
        deredfunc = 'power'
        extrapol_xmin = 22.0
        chiar_min = ext_x[0]
        #- Fit the power law to the > 22 micron wavelength range
        plsq = leastsq(Interpol.getResiduals,p0,\
                       args=(ext_x[ext_x>=extrapol_xmin],\
                             ext_y[ext_x>=extrapol_xmin],\
                             deredfunc),maxfev=20000)[0]
        #- Calculate the extrapolation and interpolation for the datagrids
        #- Then combine and apply the correction to the data
        for (dt,fn),data in self.data_raw.items():
            if len(data) == 3: 
                data_x, data_y, data_ey = data[0], data[1], data[2]
            else:
                data_x, data_y = data[0], data[1]
            
            extra = Interpol.pEval(data_x[data_x>=extrapol_xmin],plsq,deredfunc)
            inter = interp1d(ext_x,ext_y)(data_x[(data_x<extrapol_xmin)\
                                                    *(data_x>=chiar_min)])
            short = interp1d(ext_car_x,ext_car_y)(data_x[data_x<chiar_min])
            corr = hstack([short,inter,extra])
            if self.plot_extrapol_extinction: 
                Plotting2.plotCols(x=[ext_x,data_x],y=[ext_y,corr],\
                                   xlogscale=1,ylogscale=1)
            if len(data) == 3: 
                self.data[(dt,fn)] = (data_x,\
                                      data_y*10**(corr*self.ak*0.4),\
                                      data_ey*10**(corr*self.ak*0.4))
            else:
                self.data[(dt,fn)] = (data_x,data_y*10**(corr*self.ak*0.4))
                    
        

def getExtinctionCurve(gal_position='ism',curve_type='chiar_tielens',\
                       av_to_ak_conv=0.112):
    
    """
    Read extinction curve.
    
    @keyword gal_position: galactic position star: 'ism' or 'gc'. Difference of 
                           extinction curve for ISM or galactic center 
                           (as for chiar_tielens), if no difference 
                           (as for cardelli) this parameter is ignored. 
                            
                           (default: 'ism')
    @type gal_position: string
    @keyword curve_type: type of extinction curve 
                         ('chiar_tielens' or 'cardelli')
                            
                         (default: 'chiar_tielens')
    @type curve_type: string
    @keyword av_to_ak_conv: V-band to K-band conversion factor for extinction 
                            correction
                            
                            (default: 0.112)
    @type av_to_ak_conv: float
    
    @return: The wavelength grid and extinction curve values
    @rtype: (array,array)
    
    """
    
    extcurve_input = DataIO.readFile(os.path.join(cc.path.aux,\
                                                  curve_type + '.dat'),\
                                     delimiter=' ')
    extcurve_x = [float(row[0]) for row in extcurve_input if row[0][0] != '#']
    if gal_position.lower() == 'ism' and curve_type.lower() == 'chiar_tielens':
        extcurve_y = [float(row[2]) 
                      for row in extcurve_input 
                      if row[0][0] != '#' and len(row) == 3]
        #- Cut off excess wavelengths
        extcurve_x = extcurve_x[:len(extcurve_y)]                
    elif curve_type.lower() == 'cardelli':
        #- putting wavelength in micron
        extcurve_x = array(extcurve_x)*10**(-4)                  
        extcurve_y = [float(row[1]) 
                      for row in extcurve_input 
                      if row[0][0] != '#']
        #- Put values to Alambda/Ak
        extcurve_y = array(extcurve_y)/(3.1*av_to_ak_conv)       
    else:
        extcurve_y = [float(row[1]) 
                      for row in extcurve_input 
                      if row[0][0] != '#']
    
    return (array(extcurve_x),array(extcurve_y))
