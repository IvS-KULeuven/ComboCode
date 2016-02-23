# -*- coding: utf-8 -*-

"""
Performing statistics on SED data. Currently includes photometry only.

Author: R. Lombaert

"""

import os
from scipy.interpolate import interp1d
import operator
import types
import numpy as np

import cc.path
from cc.tools.io import DataIO
from cc.data import Sed
from cc.modeling.codes import MCMax
from cc.modeling.tools import Reddening
from cc.statistics import BasicStats
from cc.statistics.Statistics import Statistics



class SedStats(Statistics):
    
    """
    Environment with several tools to perform statistics on SED photometry.
    
    Then run several methods to set data, and models. e.g. (requires a 
    star_grid to be defined):
    
    >>> sed = Sed.Sed('rscl')
    >>> sedstats = SedStats.SedStats(star_name='rscl',path_code='codeSep2015')
    >>> sedstats.setInstrument(sed=sed)
    >>> sedstats.setModels(star_grid=star_grid)
    >>> sedstats.setModelPhotometry()
    
    Alternatively, you can take the SedStats object from a ComboCode object
    given that STATISTICS=1 in the inputfile for CC. Then the above is 
    already done for you. 
    
    >>> model = ComboCode.ComboCode('inputfile.dat')
    >>> model.startSession()
    >>> sedstats = model.sedstats
    
    Now you can go ahead and calculate chi^2 values and write them to a 
    file:
    
    >>> chi2 = sedstats.calcChi2()
    >>> sedstats.writeChi2('myfile.dat')
    
    The writeChi2 method writes the latest chi2 calculation to a file with
    given filename. Several options are available for selecting and sorting
    photometry during the calculation, and sorting and adding extra columns
    of model parameters when writing the results to a file. Check the 
    docstrings of each method for more information.
    
    """
        
    def __init__(self,star_name,code='MCMax',path_code='codeSep2015'):        
        
        """ 
        Initializing an instance of SedStats.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
        @keyword code: the code used for producing your output 
        
                       (default: 'MCMax')
        @type code: string
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2015')
        @type path_code: string
        
        """
        
        super(SedStats,self).__init__(star_name=star_name,\
                                       code=code,path_code=path_code)
        
        #-- Initiating a series of model related variables
        #   *_ivs: IvS repo related photometry for which photbands are available
        #   *_other: Other photometry, for which interpolation is needed
        self.mphot_ivs = []
        self.mphot_other = dict()
        
        #-- Initiating a series of data related variables
        #   Only one possible IvS phot file, multiple for other phot
        self.dphot_ivs = dict([('wave',np.empty(0)),\
                               ('phot',np.empty(0)),\
                               ('ephot',np.empty(0))])
        self.dphot_other = dict()
        self.photbands = np.empty(0)
        self.sed = None
        
        self.chi2 = np.empty(0)
        
        
        
    def setInstrument(self,sed=None,**kwargs):
       
        '''
        Set and read the data objects for this statistics module. 
        
        The Sed object is made on the spot if sed is None. Extra arguments for
        making the Sed can be passed through the method in kwargs.
        
        @keyword sed: The Sed object. In case of default, it is made on the spot
                      for given star with extra arguments kwargs if needed
                     
                      (default: None)
        @type sed: Sed()
        
        '''
        
        super(SedStats,self).setInstrument(instrument_name='SED',\
                                           instrument_instance=sed,**kwargs)
        self.sed = self.instrument
        self.photbands = sed.photbands
        
        for (dt,fn),tdata in sorted([dset
                                     for dset in self.sed.data.items()
                                     if 'PHOT' in dset[0][0].upper()]):
            if dt.lower() == 'photometric_ivs':
                self.dphot_ivs['wave'] = tdata[0]
                self.dphot_ivs['phot'] = tdata[1]
                self.dphot_ivs['ephot'] = tdata[2]
            else:
                self.dphot_other[fn] = dict()
                self.dphot_other[fn]['wave'] = tdata[0]
                self.dphot_other[fn]['phot'] = tdata[1]
                self.dphot_other[fn]['ephot'] = tdata[2]
            
            
    def setModelPhotometry(self):
    
        '''
        Prepare the model photometry to be compared with the data. 
        
        Two kinds: The IvS photometry with proper photometric bands, and other
        photometry to be compared with the interpolated model spectrum. 
        
        '''
        
        
        #- Collect model data and set ak if needed
        mids = [s['LAST_MCMAX_MODEL'] for s in self.star_grid]
        if not mids:
            print "No successfully calculated MCMax models found."
            return
        
        self.mwave = []
        self.mflux = []
        for model_id,s in zip(mids,self.star_grid):
            dpath = os.path.join(cc.path.mout,'models',model_id)
            fn_spec = 'spectrum{:04.1f}.dat'.format(s['RT_INCLINATION'])
            w,f = MCMax.readModelSpectrum(dpath,s['RT_SPEC'],fn_spec)
            if s['REDDENING']:
                print 'Reddening models to correct for interstellar extinction.'
                ak = self.sed.getAk(s['DISTANCE'],s['REDDENING_MAP'],\
                                    s['REDDENING_LAW'])
                f = Reddening.redden(w,f,ak,law=s['REDDENING_LAW'])
            self.mwave.append(w)
            self.mflux.append(f)
        
        if self.photbands.size:
            self.mphot_ivs = [Sed.calcPhotometry(w,f,self.photbands)
                              for w,f in zip(self.mwave,self.mflux)]
            
        for fn in self.dphot_other.keys():
            self.mphot_other[fn] = []
        if self.dphot_other.keys():
            for w,f in zip(self.mwave,self.mflux):
                interp = interp1d(w,f)
                for fn in self.dphot_other.keys():
                    finter = interp(self.dphot_other[fn]['wave'])
                    self.mphot_other[fn].append(finter)
        
        
        
    def calcChi2(self,ndf=0,fns=None,cwave=0.0,phot_ivs=1,sort=1,\
                 chi2_method='diff'):
    
        '''
        Calculate, save and return the chi^2 values for a given set of models.
        
        For now, only photometry (both with and without photometric bands 
        available, the latter being associated with the Photometric_IvS file) 
        is taken into account. 
        
        You can specify the number of degrees of freedom for the reduced chi^2 
        in case you are comparing models probing a grid across ndf different 
        parameters. 
        
        The method overrides previously calculated chi2 in the object. Write the
        chi^2 to a file before running the method if you want to save the values
        with self.writeChi2()
        
        
        @keyword ndf: Number of degrees of freedom. Default in case of 
                      calculating for one single model. Typically the number of
                      variable grid parameters in a grid calculation.
                  
                      (default: 0) 
        @type ndf: int
        @keyword fns: The photometric filenames to be included in the chi^2 calc
                      Default if all are to be included. Empty list if none are
                      to be included.
                      
                      (default: None)
        @type fns: list(str)
        @keyword cwave: A cut-off wavelength in micron used to calculate the 
                        chi^2 stat for all photometry above or equal to the 
                        wavelength in micron. Use a negative value to grab all
                        wavelengths lower than cwave. Default if no cut-off.
                        
                        (default: 0.0)
        @type cwave: float
        @keyword phot_ivs: Include the Photometric_IvS photometry as well.
        
                           (default: 1)
        @type phot_ivs: bool
        @keyword sort: Sort the chi^2 values from lowest to highest. 
        
                       (default: 1)
        @type sort: bool
        @keyword chi2_method: Method for calculating chi^2. Can be diff or 
                              log (see BasicStats for more information)
        
                              (default: 'diff')
        @type chi2_method: str
        
        @return: The list of chi^2 values with the same length as the model grid
                 of valid MCMax models. 
        @rtype: list
        
        '''
        
        #-- Don't use [[]]*len(self.star_grid). This only multiplies the pointer
        #   to the same list, it doesn't do copy([]).
        mphot = [[] for s in self.star_grid]
        dphot = []
        ephot = []
        
        if fns is None:
            fns = self.dphot_other.keys()
        else:
            fns = [fn if os.path.split(fn)[0] else os.path.join(cc.path.dsed,fn)
                   for fn in fns]
            fns = [fn for fn in fns if fn in self.dphot_other.keys()]
        
        #-- When no valid filenames are given, and IvS phot is not requested, 
        #   just forget the filename selection. 
        if not fns and not phot_ivs: 
            fns = self.dphot_other.keys()
            
        for fn in fns:
            #-- Includes cwave == 0.0, in which case everything is selected
            if cwave >= 0.0:
                keep = self.dphot_other[fn]['wave'] >= cwave
            else:
                keep = self.dphot_other[fn]['wave'] < -cwave
            
            #-- Check if anything is kept at all, otherwise move on
            if not keep.any():
                continue
                
            #-- Add the model and data points to the list
            for i,iphot in enumerate(self.mphot_other[fn]):
                mphot[i].extend(iphot[keep])
            dphot.extend(self.dphot_other[fn]['phot'][keep])
            ephot.extend(self.dphot_other[fn]['ephot'][keep])
        
        if phot_ivs and self.photbands.size:
            if cwave >= 0.0:
                keep = self.dphot_ivs['wave'] >= cwave
            else:
                keep = self.dphot_ivs['wave'] < -cwave
            
            #-- Check if anything is kept at all, otherwise move on
            if keep.any():
                for i,iphot in enumerate(self.mphot_ivs):
                    mphot[i].extend(iphot[keep])
                dphot.extend(self.dphot_ivs['phot'][keep])
                ephot.extend(self.dphot_ivs['ephot'][keep])
        
        #-- If no data are selected at all, return an empty array
        if not dphot:
            self.chi2 = np.empty(0)
            return self.chi2
            
        self.chi2 = []
        for iphot in mphot:
            self.chi2.append(BasicStats.calcChiSquared(dphot,iphot,ephot,ndf,\
                                                       chi2_method))
        self.chi2 = np.array(self.chi2)
        
        return np.sort(self.chi2) if sort else self.chi2
    
    
    
    def getStarGrid(self,sort=1):
        
        '''
        Return the star grid for which the chi^2 is calculated. 
        
        This grid may differ from the original ComboCode grid because the local
        star_grid is filtered for unavailable MCMax models (either not 
        calculated or not loaded). 
        
        Each object in the star_grid is itself a dictionary that contains all
        the parameters of the model.
        
        @keyword sort: Sort the star_grid according to the chi^2 values from 
                       lowest to highest. Requires calcChi2 to be ran first.
        
                       (default: 1)
        @type sort: bool
        
        @return: The star grid
        @rtype: list(Star())
        
        '''
        
        if not self.chi2.size:
            sort = 0
        if sort:
            isort = np.argsort(self.chi2)
        
        return self.star_grid[isort] if sort else self.star_grid



    def writeChi2(self,fn,sort=1,parameters=[]):
        
        '''
        Write the Chi^2 values to a file. Lists the model id in the first column
        with the chi^2 value in the second. 
        
        The chi^2 values can be requested to be sorted.
        
        Parameters from the Star() objects can be added as additional columns.
        Given parameters must be valid.
        
        @param fn: The output filename
        @type fn: str
        
        @keyword sort: Sort the star_grid according to the chi^2 values from 
                       lowest to highest. Requires calcChi2 to be ran first.
        
                       (default: 1)
        @type sort: bool
        @keyword parameters: The additional model parameters to be added as 
                             columns in the file. 
                             
                             (default: [])
        @type parameters: list(str)
        
        '''
        
        #-- If no chi^2 was calculated, do nothing
        if not self.chi2.size:
            return
            
        #-- Write the header
        comments = ['# '] + ['ID','RedChi^2'] + parameters + ['\n']
        DataIO.writeFile(filename=fn,input_lines=comments,delimiter='\t')
        
        #-- Define the columns
        cols = [[s['LAST_MCMAX_MODEL'] for s in self.getStarGrid(sort=sort)]]
        if sort: 
            isort = np.argsort(self.chi2) 
            cols.append(self.chi2[isort])
        else: 
            cols.append(self.chi2)
        
        #-- Add additional model parameters if requested
        for par in parameters:
            cols.append([s[par] for s in self.getStarGrid(sort=sort)])

        #-- Append the columns to the file after the header
        DataIO.writeCols(filename=fn,cols=cols,mode='a')