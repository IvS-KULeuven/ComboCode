# -*- coding: utf-8 -*-

"""
A class for reading and managing Sphinx output.

Author: R. Lombaert

"""

import os
from scipy import array 
from scipy import isnan
from scipy import isfinite

from cc.tools.io.Reader import Reader
from cc.tools.io import DataIO



class SphinxReader(Reader):
    
    '''
    A Reader for Sphinx output files. 
    
    '''
    
    def __init__(self,filename):
        
        '''
        Creating a Sphinx object ready for reading Sphinx output.
        
        Reading sphinx output through full filename, including filepath.
        
        The sphinxfile number is given as *, which is then replaced as needed.
        
        @param filename: The sphinx filename, including filepath. The sphinx
                         file number is given by *.
        @type filename: string
        
        '''
        
        super(SphinxReader, self).__init__()
        self.filename = filename.replace('sph1','sph*').replace('sph2','sph*')
        self.nans_present=False
        self.parseImpact()
        self.parseProfile()
        
    
    def parseImpact(self):
        
        ''' 
        Parse sphinx file 1, which includes all the impact parameter info. 
        
        The output is stored in dict self.sph1.
        
        '''
        
        
        self.sph1 = dict()
        self.contents['sph1'] = self.sph1
        data = DataIO.readCols(self.filename.replace('*','1'),start_row=1)
        self.sph1['impact'] = data[0]
        self.sph1['norm_intens'] = data[1] 
        self.sph1['weighted_intens'] = data[2] 
        self.sph1['sum_intens_p'] = data[3]
        self.sph1['sum_intens'] = data[4]
        
        
        
    def parseProfile(self):
        
        '''
        Parse the sphinx file 2, which includes all line profile info.
        
        The output is stored in dict self.sph2.
        
        '''
        
        self.sph2 = dict()
        self.contents['sph2'] = self.sph2
        self.sph2['nobeam'] = dict()
        self.sph2['beam'] = dict()
        self.sph2['nobeam_cont'] = dict()
        self.sph2['beam_cont'] = dict()
        data = self.getFile(self.filename.replace('*','2'))
        data_col_1 = [d[0] for d in data]
        data_i = 6
        data_j = DataIO.findString(data_i,data_col_1)
        self.sph2['nobeam']['velocity'] = array([float(line[0]) 
                                              for line in data[data_i:data_j]])
        #-- Reverse this flux grid. Sphinx output files give the mirrored
        #   flux grid for the associated velocity grid.
        self.sph2['nobeam']['flux'] =     array([DataIO.convertFloat(line[-1],\
                                                                     nans=1) 
                                              for line in data[data_i:data_j]])
        self.sph2['nobeam']['flux'] = self.sph2['nobeam']['flux'][::-1]
        data_k = data_j + 4
        data_l = DataIO.findString(data_k,data_col_1)
        self.sph2['beam']['velocity'] =  array([float(line[0]) 
                                              for line in data[data_k:data_l]])
        self.sph2['beam']['flux'] =      array([float(line[-1]) 
                                              for line in data[data_k:data_l]])  
        self.sph2['beam']['norm_flux'] = array([float(line[1]) 
                                              for line in data[data_k:data_l]])
        self.sph2['beam']['tmb'] =       array([float(line[2]) 
                                              for line in data[data_k:data_l]])
        
        self.setContinuum('nobeam','flux')
        for lp in ['flux','norm_flux','tmb']:
            self.setContinuum('beam',lp)
        if self.sph2['beam']['velocity'][0] > self.sph2['beam']['velocity'][-1]:
            self.sph2['beam']['velocity'] = self.sph2['beam']['velocity'][::-1]
            self.sph2['beam']['flux'] = self.sph2['beam']['flux'][::-1]
            self.sph2['beam']['norm_flux'] = self.sph2['beam']['norm_flux'][::-1]
            self.sph2['beam']['tmb'] = self.sph2['beam']['tmb'][::-1]
        if self.sph2['nobeam']['velocity'][0] > self.sph2['nobeam']['velocity'][-1]:
            self.sph2['nobeam']['velocity'] = self.sph2['nobeam']['velocity'][::-1]
            self.sph2['nobeam']['flux'] = self.sph2['nobeam']['flux'][::-1]
        if True in list(isnan(self.sph2['nobeam']['flux'])):
            self.nans_present = True
            print "WARNING! There are NaN's in the intrinsic line profile " + \
                  "with model id %s:"\
                  %(os.path.split(os.path.split(self.filename)[0])[1])
            print os.path.split(self.filename.replace('sph*','sph2'))[1]
    
    
    def setContinuum(self,beam,lp):
        
        '''
        Set the continuum value for this line profile. 
        
        @param beam: Either 'beam' or 'nobeam'.
        @type beam: str
        @param lp: The type of line profile
        @type lp: str
        
        '''
        
        flux = self.sph2[beam][lp] 
        flux = flux[isfinite(flux)]
        continuum = (flux[0] + flux[-1])/2.
        self.sph2[beam+'_cont'][lp] = continuum
                
        
    def getLPIntrinsic(self,cont_subtract=1):
        
        '''
        Return the intrinsic flux line profile.
        
        Continuum subtraction can be requested. 
        
        @keyword cont_subtract: Subtract the continuum value outside the line
                                from the whole line profile. 
        @type cont_subtract: bool
        
        @return: The intrinsic flux
        @rtype: list
        
        '''
        
        flux = self.sph2['nobeam']['flux']
        if cont_subtract:
            flux = flux - self.sph2['nobeam_cont']['flux']
        
        return flux
        
    
    def getVelocityIntrinsic(self):
        
        '''
        Return the velocity grid for the intrinsic line profile. 
        
        @return: The velocity grid
        @rtype: list
        
        '''
        
        return self.sph2['nobeam']['velocity']
    

    def getVelocity(self):
        
        '''
        Return the velocity grid for the convolved line profile. 
        
        @return: The velocity grid
        @rtype: list
        
        '''
        
        return self.sph2['beam']['velocity']
        
        
    def getLPTmb(self,cont_subtract=1):
        
        '''
        Return the main beam temperature line profile.
        
        Continuum subtraction can be requested. 
        
        @keyword cont_subtract: Subtract the continuum value outside the line
                                from the whole line profile. 
                                
                                (default: 1)
        @type cont_subtract: bool
        
        @return: The main beam temperature
        @rtype: list
        
        '''
        
        tmb = self.sph2['beam']['tmb']
        if cont_subtract:
            tmb = tmb - self.sph2['beam_cont']['tmb']

        return tmb
        
        
        
    def getLPConvolved(self,cont_subtract=1):
        
        '''
        Return the line profile after convolution with the beam profile. 
        
        Continuum subtraction can be requested. 
        
        @keyword cont_subtract: Subtract the continuum value outside the line
                                from the whole line profile. 
        @type cont_subtract: bool
        
        @return: The convolved flux
        @rtype: list
        
        '''
        
        flux = self.sph2['beam']['flux']
        if cont_subtract:
            flux = flux - self.sph2['beam_cont']['flux']
        
        return flux
        
        
        
    def getLPNormalized(self,cont_subtract=1):
        
        '''
        Return the normalized line profile after convolution with the beam 
        profile.
        
        Continuum subtraction can be requested. 
        
        @keyword cont_subtract: Subtract the continuum value outside the line
                                from the whole line profile. 
        @type cont_subtract: bool
        
        @return: The normalized + convolved flux
        @rtype: list
        
        '''
        
        flux = self.sph2['beam']['norm_flux']
        if cont_subtract:
            flux = flux - self.sph2['beam_cont']['norm_flux']
        
        return flux        
        
    
    def getWeightedIntensity(self):
        
        '''
        Return the weighted intensity with respect to impact parameter squared.
        
        @return: The weighted intensity list.
        @rtype: list
        
        '''
        
        return self.sph1['weighted_intens']
        
        
        
    def getNormalizedIntensity(self):
        
        '''
        Return the normalized weighted intensity with respect to 
        impact parameter squared.
        
        @return: The normalized weighted intensity lists.
        @rtype: list
        
        '''
        
        return self.sph1['norm_intens']
    
    
    
    def getImpact(self):
        
        '''
        Return the impact parameter grid.
        
        @return: the grid of impact parameters in the sphinx outputfile 1.
        @rtype: list
        
        '''
        
        return self.sph1['impact']