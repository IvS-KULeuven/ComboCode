# -*- coding: utf-8 -*-

"""
Tools for making and managing opacity files.

Author: R. Lombaert

"""

from scipy import array,argmin
from glob import glob
from math import pi
import os
import types

import cc.path
from cc.tools.io import DataIO
from cc.tools.numerical import Interpol


def massFractionGSD(acut,amin=0.01,amax=100.,slope=-3.5):
    
    '''
    Calculate the mass fraction of a subset of grain sizes from a grain size
    distribution.
    
    A minimum of two subsets is created, so two masses are returned. Can be 
    anything higher than that, depending on how many cuts are requested.
    
    @param acut: The grain size(s) at which the cuts are made. Can be a float 
                 or a list. 
    @type acut: list(float)
    
    @keyword amin: The minimum grain size in the distribution
                    
                   (default: 0.01)
    @type amin: float
    @keyword amax: The maximum grain size in the distribution
                    
                   (default: 100)
    @type amax: float
    @keyword slope: The slope of the a-dependence of the distribution. Default
                    is the MRN value of -3.5
                    
                    (default: -3.5)
    @type slope: float
    
    @return: The mass fraction of the first slice, the second slice, etc, 
             depending on how many cuts are requested. Minimum of two.
    @rtype: [float,float,...]         
    
    '''
    
    if not type(acut) is types.ListType:
        acut = [float(acut)]
    acut = sorted(acut)
    #-- Slope is negative typically. If slope == -4, this doesn't work
    if slope == -4.: return np.empty(0)
    #-- n(a) = Cst * Int(a**3 * a**-3.5 da) ==> Cst * [a**0.5]^amax_amin
    p = (3.+slope)+1.
    subsets = []
    #-- first subset:
    subsets.append(acut[0]**p-amin**p)    
    #-- additional subsets if len(acut) > 1
    for i in xrange(len(acut[1:])):
        subsets.append(acut[i+1]**p-acut[i]**p)
    #-- final subset
    subsets.append(amax**p-acut[-1]**p)
    #-- Determine the fractions with respect to the total number density.
    #   Constant does not matter because of it being ratios.
    fractions = array(subsets)/sum(subsets)
    return fractions
    


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
    
    path = os.path.join(cc.path.mopac,species)
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


class CustomOpacity():
    """
    An interface for creating custom opacity files by taking the original and 
    tinkering with it.
    
    """
    
    def __init__(self,species):
        
        """
        Initializing an instance of the custom opacities interface.
    
        @param species: the species short name
        @type species: string
        
        """
        
        self.setSpecies(species)
        self.opacity_file = False

    
    def setSpecies(self,species):
        
        """
        Change the species of the current CustomOpacity instance.
        
        @param species: the species short name
        @type species: string
        
        """
        
        self.species = species
        self.index = DataIO.getInputData(keyword='SPECIES_SHORT',\
                                         filename='Dust.dat')\
                                        .index(self.species)
        self.filename =  DataIO.getInputData(keyword='PART_FILE',\
                                             filename='Dust.dat',\
                                             rindex=self.index)
        fn = os.path.join(cc.path.mopac,self.filename)
        self.input_data = DataIO.readFile(filename=fn,delimiter=' ') 
        
        
        
    def makeOpa(self,mode='ZERO',**args):
        
        """
        Making custom .particle files.
        
        Every method called here will put the results in self.output_data.
        
        @keyword mode: type of extrapolation (ZERO,FUNCTION,HONY,ZHANG)
                        
                       (default: 'ZERO')
        @type mode: string
        @keyword args: Optional keywords required for the other methods of the 
                       class
        @type args: dict
        
        """
        
        self.output_data = []
        mode = mode.upper()
        if hasattr(self,'do' + mode):
            getattr(self,'do' + mode)(**args)
            self.output_data = [' '.join(str(line)) 
                                for line in self.output_data]
            output_filename = '_'.join(['customOpacity',mode] + \
                                       sorted(args.values()) + \
                                       [self.filename])
            if self.opacity_file:
                output_filename.replace('.particle','.opacity')
            DataIO.writeFile(filename=os.path.join(cc.path.mopac,\
                                                   output_filename),\
                             input_lines=self.output_data)
            new_short = self.species + mode
            #- filename is already present: the new file should have the same
            #- parameters and the short name can be kept, 
            #- nothing is changed in the Dust.dat file
            try:    
                DataIO.getInputData(keyword='PART_FILE',filename='Dust.dat')\
                                   .index(output_filename)
            #- filename is not present: do the normal procedure, ie check if 
            #- short name is already present
            except ValueError:        
                i=0
                while ' '.join(DataIO.getInputData(keyword='SPECIES_SHORT',\
                                                   filename='Dust.dat'))\
                                                  .find(new_short) != -1:
                    i+=1    
                    new_short = new_short + str(i)
                adding_line = [new_short] + \
                              [str(DataIO.getInputData(keyword=key,\
                                                       filename='Dust.dat',\
                                                       rindex=self.index))
                               for key in ['SPEC_DENS','T_DES','T_DESA','T_DESB']]
                adding_line.insert(2,output_filename)
                adding_line = '\t\t'.join(adding_line)
                DataIO.writeFile(os.path.join(cc.path.usr,'Dust.dat'),\
                                 [adding_line+'\n'],mode='a')
        else:
            print 'Mode "' + mode + '" not available. Aborting.'
        
        
        
    def doZERO(self,wl=10.0):
        
        """
        Replace all opacities by zero below a cut-off wavelength.
    
        @keyword wl: the wavelength of the discontinuity to zero
        
                     (default: 10.0) 
        @type wl: float
        
        """
        
        for line in self.input_data:
            if len(line) == 4 and float(line[0]) < wl:
                self.output_data.append([line[0],str(2e-60),str(1e-60),str(1e-60)])
            else:
                self.output_data.append(line)



    def doFUNCTION(self,wl_min=10.0,wl_max=20.0,function='linear'):
        
        """
        Replace all opacities by an extrapolation below a cut-off wavelength.
    
        @keyword wl_min: lower boundary interpolation range and the cutoff 
                         wavelength for the extrapolation
                         
                         (default: 10.0)
        @type wl_min: float
        @keyword wl_max: upper boundary interpolation range
        
                         (default: 20.0)
        @type wl_max: float
        @keyword function: type of function used for the inter/extrapolation
                           See Interpol.pEval for function types    
                        
                           (default: 'linear')
        @type function: string 
        
        """
        
        #raise TypeError("WARNING! The CustumOpa().doFUNCTION method is obsolete! New version NYI.")
        
        #- Select relevant inputlines (not saving the scattering matrices)
        self.opacity_file = True
        inputsel = [line for line in self.input_data if len(line) == 4]
        inputsel = array([[float(f) for f in line] for line in inputsel])
        wl = inputsel[:,0]
        function = function.lower()
        
        #- Select the extrapolation and interpolation regions.
        wl_low = wl[wl < wl_min]
        wlinter = wl[(wl >= wl_min) * (wl <= wl_max)]
        abs_inter = inputsel[:,2][(wl >= wl_min) * (wl <= wl_max)]
        sca_inter = inputsel[:,3][(wl >= wl_min) * (wl <= wl_max)]
        
        #- Fit the interpolation region and extrapolate this to lower wavelength
        abs_p0 = [abs_inter[0],wl[0],1,1]
        sca_p0 = [sca_inter[0],wl[0],1,1]
        abs_low = Interpol.fitFunction(x_in=wlinter,y_in=abs_inter,\
                                       x_out=wl_low,func=function,\
                                       initial=abs_p0,show=1)        
        abs_low = array([a > 0 and a or 0 for a in abs_low])
        sca_low = Interpol.fitFunction(x_in=wlinter,y_in=sca_inter,\
                                       x_out=wl_low,func=function,\
                                       initial=sca_p0,show=1)
        sca_low = array([s > 0 and s or 0 for s in sca_low])
        ext_low = abs_low + sca_low
        
        #- Set the output data
        self.output_data = [array([w,el,al,sl])
                            for w,el,al,sl in zip(wl_low,ext_low,\
                                                  abs_low,sca_low)]
        self.output_data.extend(inputsel[wl >= wl_min])
                            
                            
                            
    def doHONY(self,wl1 = 1.0,wl2=2.0,wl3=10.0,a_mod=0.05,q_cst=1.0):
        
        """
        Replace all opacities below a cut-off wavelength with values as assumed 
        by Hony et al (2003, 2004)
    
        @keyword wl1: first discontinuity in micron
        
                      (default: 1.0)
        @type wl1: float
        @keyword wl2: second discontinuity in micron
                      
                      (default: 2.0)
        @type wl2: float
        @keyword wl3: third discontinuity in micron
        
                      (default: 10)
        @type wl3: float
        @keyword a_mod: model grain size, 0.01 according to Hony 2004 in micron
        
                        (default: 0.05)
        @type a_mod: float
        @keyword q_cst: constant extinction efficiency at wavelengths < wl1
        
                        (default: 1.0)
        @type q_cst: float
        
        """

        spec_dens = DataIO.getInputData(keyword='SPEC_DENS',rindex=self.index,\
                                        filename='Dust.dat')
        opa_cst = q_cst/4.0*3.0/spec_dens/(a_mod*10**(-4))      
        for line in self.input_data:
            if len(line) == 4 and float(line[0]) < wl1:
                self.output_data.append([line[0],str(opa_cst+1e-60),\
                                         str(opa_cst),str(1e-60)])
            elif len(line)==4 and float(line[0])>=wl1 and float(line[0]) < wl2:
                opa = Interpol.linInterpol([wl1,wl2],[opa_cst,1e-60],\
                                            float(line[0]))
                self.output_data.append([line[0],str(opa+1e-60),str(opa),\
                                         str(1e-60)])
            elif len(line)==4 and float(line[0])>=wl2 and float(line[0]) < wl3:
                self.output_data.append([line[0],str(2e-60),str(1e-60),\
                                         str(1e-60)])
            else:
                self.output_data.append(line)
         
         
         
    def doZHANG(self,wl=10.0,a_mod=0.05,q_cst=1.6,metallic=True):
        
        """
        Replace all opacities below a cut-off wavelength with values as assumed
        by Zhang et al (2009)
    
        @keyword wl: wavelength of the discontinuity
        
                     (default: 10.0)
        @type wl: float
        @keyword a_mod: model grain size, following Zhang et al 2009
        
                        (default: 0.05)
        @type a_mod: float
        @keyword q_cst: the constant extinction efficieny at short wavelengths
        
                        (default: 1.6)
        @type q_cst: float
        @keyword metallic: if metallic dust, cut off wavelength is calculated 
                           differently than if dielectric dust (ie metallic=0)
                           
                           (default: True)
        @type metallic: bool
        
        """
        
        spec_dens = DataIO.getInputData(keyword='SPEC_DENS',rindex=self.index,\
                                        filename='Dust.dat')
        wl1 = metallic and pi*a_mod or 2*pi*a_mod
        opa_cst = q_cst/4.0*3.0/spec_dens/(a_mod*10**(-4))
        i = 0
        while float(self.input_data[i][0]) <= wl:
            i += 1
        opa_10 = float(self.input_data[i][1])
        for line in self.input_data:
            if len(line) == 4 and float(line[0]) <= wl1:
                self.output_data.append([line[0],str(opa_cst+1e-60),\
                                         str(opa_cst),str(1e-60)])
            elif len(line)==4 and float(line[0]) > wl1 and float(line[0]) < wl:
                eps = min(1,float(line[0])/wl)
                opa = (1-eps) * opa_cst*wl1/float(line[0]) + eps*opa_10
                
                self.output_data.append([line[0],str(opa+1e-60),str(opa),\
                                         str(1e-60)])
            else:
                self.output_data.append(line)
         
            