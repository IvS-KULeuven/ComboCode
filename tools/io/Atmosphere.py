# -*- coding: utf-8 -*-

"""
Reading stellar model atmospheres for use with MCMax and GASTRoNOoM.

Author: R. Lombaert

"""

from glob import glob
import pyfits, os
from scipy import rec,array,argmin



class Atmosphere(object):
    
    """
    Reading stellar model atmospheres.
    
    """
    
    def __init__(self,modeltype,filename=None,\
                 folder='/STER/pieterd/IVSDATA/sedtables/modelgrids/'):
        
        """
        Initializing an Atmosphere() object. 
        
        The stellar model atmospheres are read from a given folder and the 
        preferred model is chosen after initialisation.
        
        @param modeltype: Type of model atmosphere (any from IVSDATA folder)
        @type modeltype: string
        
        @keyword filename: The filename of the model atmosphere in the 
                           modelgrids folder. If None, the filename can be set
                           later.
                           
                           (default: None)
        @type filename: string
        @keyword folder: The folder where the model grids are stored
        
                         (default: '/STER/100/pieterd/IVSDATA/sedtables/modelgrids/')
        @type folder: string
                
        """
        
        self.modeltype = modeltype
        self.folder = folder
        self.filename = filename
        if self.filename <> None:
            self.filepath = os.path.join(self.folder,self.filename)
        self.modellist = glob(os.path.join(self.folder,modeltype+'*'))
        self.modelgrid = None
        self.teff_actual = None
        self.logg_actual = None



    def getModelList(self):
        
        """
        Retrieve a list of models available in the modelgrids folder
        for the requested model type.
        
        When setting a model, the model is identified by its index in this list.
        
        @return: The filepaths available for given modeltype. 
        @rtype: list[string]
        
        """
        
        return self.modellist
        
    
    
    def setModelGrid(self,index):
        
        """
        Set the model atmosphere to be read, identified by the index in the 
        model list.
        
        Once set here, or if given upon initialisation of the Atmosphere() 
        object, the model cannot be changed for this instance.
        
        """
        
        if self.filename is None:
            self.filename = os.path.split(self.modellist[index])[1]
            self.filepath = os.path.join(self.folder,self.filename)


    
    def readModelGrid(self):
        
        """
        Read the model atmosphere fits file.
        
        """
        
        self.ff = pyfits.open(self.filepath)
        self.header = self.ff[0].header
        teffs = array([self.ff[i].header['TEFF'] 
                       for i in xrange(1,len(self.ff))])
        loggs = array([self.ff[i].header['LOGG'] 
                       for i in xrange(1,len(self.ff))])
        self.modelgrid = rec.fromarrays([array(range(1,len(self.ff))),teffs,\
                                         loggs],\
                                        names=['INDEX','TEFF','LOGG'])
                                        
    
    
    def getHeader(self):
        
        """
        Return the main header of the fits file.
        
        @return: The header information
        @rtype: pyfits.NP_pyfits.Header()
        
        """
        
        return self.header
        


    def getModelGrid(self):
        
        """
        Return the TEFF and LOGG grid for this model. 
        
        @return: the grid values of the teff and logg parameters
        @rtype: recarray
        
        """
        
        return self.modelgrid
        
        
        
    def getModel(self,teff,logg):
        
        """
        Return the model atmosphere for given effective temperature and log g.
        
        Not yet scaled to the distance!
        
        Units returned are (micron,Jy)
        
        @param teff: the stellar effective temperature
        @type teff: float
        @param logg: the log g value
        @type logg: float
        
        @return: The model spectrum in (micron,Jy)
        @rtype: recarray
        
        """
        
        c = 2.99792458e18          #in angstrom/s
        if self.modelgrid is None:
            self.readModelGrid()
        mg = self.modelgrid
        #- Find the closest temperature in the grid
        teff_prox = mg['TEFF'][argmin(abs(mg['TEFF']-teff))]
        #- Select all models with that temperature
        mgsel = mg[mg['TEFF']==teff_prox]
        #- Select the closest log g in the selection
        logg_prox = mgsel['LOGG'][argmin(abs(mgsel['LOGG']-logg))]
        #- Get the index of the model closest to teff and logg
        imodel = mgsel[mgsel['LOGG']==logg_prox]['INDEX'][0]
        
        self.teff_actual = teff_prox
        self.logg_actual = logg_prox        
        
        wave = self.ff[imodel].data.field('wavelength')
        flux = self.ff[imodel].data.field('flux')
        if self.header['FLXUNIT'] == 'erg/s/cm2/A':
            #- Go to erg/s/cm2/Hz, lFl = nFn, then to Jy (factor 10**(23))
            flux = flux * wave**2 / c * 10**(23)
        else:
            raise Error('Flux unit unknown in atmosphere model fits file.')
        if self.header['WAVUNIT'] == 'angstrom':
            wave = wave * 10**(-4)
        else:
            raise Error('Wavelength unit unknown in atmosphere model fits file.')
        
        model = rec.fromarrays([wave,flux],names=['wave','flux'])        
        return model 
        
        