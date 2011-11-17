# -*- coding: utf-8 -*-

"""
A class for reading and managing Sphinx output.

Author: R. Lombaert

"""

import os
from scipy import array 
from scipy import isnan

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
        data = self.getFile(self.filename.replace('*','2'))
        data_col_1 = [d[0] for d in data]
        data_i = 6
        data_j = DataIO.findString(data_i,data_col_1)
        self.sph2['nobeam']['velocity'] = array([float(line[0]) 
                                              for line in data[data_i:data_j]])
        self.sph2['nobeam']['flux'] =     array([DataIO.convertFloat(line[-1],\
                                                                     nans=1) 
                                              for line in data[data_i:data_j]])
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
        if True in list(isnan(self.sph2['nobeam']['flux'])):
            self.nans_present = True
            print "WARNING! There are NaN's in the intrinsic line profile " + \
                  "with model id %s:"\
                  %(os.path.split(os.path.split(self.filename)[0])[1])
            print os.path.split(self.filename.replace('sph*','sph2'))[1]
    
    
    def getLPIntrinsic(self):
        
        '''
        Return the intrinsic flux line profile.
        
        @return: The intrinsic flux
        @rtype: list
        
        '''
        
        return self.sph2['nobeam']['flux']
        
    
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
        
        
    def getLPTmb(self):
        
        '''
        Return the main beam temperature line profile.
        
        If the continuum is non-zero, it is subtracted from the line profile. 
        
        @return: The main beam temperature
        @rtype: list
        
        '''
        
        tmb = self.sph2['beam']['tmb']
        continuum = (tmb[0] + tmb[-1])/2.
        return tmb - continuum
        
        
        
    def getLPConvolved(self):
        
        '''
        Return the line profile after convolution with the beam profile. 
        
        @return: The convolved flux
        @rtype: list
        
        '''
        
        return self.sph2['beam']['flux']
        
        
        
    def getLPNormalized(self):
        
        '''
        Return the normalized line profile after convolution with the beam profile.
        
        @return: The normalized + convolved flux
        @rtype: list
        
        '''
        
        return self.sph2['beam']['norm_flux']
        
        
    
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