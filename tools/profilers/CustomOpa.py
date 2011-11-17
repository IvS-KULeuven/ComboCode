# -*- coding: utf-8 -*-

"""
Tools for making opacity files.

Author: R. Lombaert

"""

from scipy import array
import os
from math import pi

from cc.tools.io import DataIO
from cc.tools.numerical import Interpol
from cc.modeling.objects import Star



class CustomOpa():
    """
    An interface for creating custom opacity files by taking the original and 
    tinkering with it.
    
    """
    
    def __init__(self,species,\
                 path_cc=os.path.join(os.path.expanduser('~'),'ComboCode',\
                                      'Data'),\
                 path_mcmax=os.path.join(os.path.expanduser('~'),'MCMax')):
        
        """
        Initializing an instance of the custom opacities interface.
    
        @param species: the species short name
        @type species: string
        @keyword path_cc: The CC home folder
        
                          (default: ~/ComboCode/)
        @type path_cc: string
        @keyword path_mcmax: The MCMax home folder
                    
                     (default: ~/MCMax/)
        @type path_mcmax: string
        
        """
        
        self.path_cc = path_cc
        self.path_mcmax = path_mcmax
        self.setSpecies(species)
        self.opacity_file = False

    
    def setSpecies(self,species):
        
        """
        Change the species of the current CustomOpa instance.
        
        @param species: the species short name
        @type species: string
        
        """
        
        self.species = species
        self.index = Star.getInputData(keyword='SPECIES_SHORT',\
                                       filename='Dust.dat',\
                                       path=self.path_cc).index(self.species)
        self.filename =  Star.getInputData(keyword='PART_FILE',\
                                           filename='Dust.dat',\
                                           path=self.path_cc)[self.index]
        self.input_data = DataIO.readFile(\
                                filename=os.path.join(self.path_mcmax,'src',\
                                                      self.filename),\
                                delimiter=' ') 
        
        
        
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
            output_filename = '_'.join(['customOpa',mode] + \
                                       sorted(args.values()) + \
                                       [self.filename])
            if self.opacity_file:
                output_filename.replace('.particle','.opacity')
            DataIO.writeFile(filename=os.path.join(self.path_mcmax,'src',\
                                                   output_filename),\
                             input_lines=self.output_data)
            new_short = self.species + mode
            #- filename is already present: the new file should have the same
            #- parameters and the short name can be kept, 
            #- nothing is changed in the Dust.dat file
            try:    
                Star.getInputData(keyword='PART_FILE', filename='Dust.dat',\
                                  path=self.path_cc).index(output_filename)
            #- filename is not present: do the normal procedure, ie check if 
            #- short name is already present
            except ValueError:        
                i=0
                while ' '.join(Star.getInputData(keyword='SPECIES_SHORT',\
                                                 filename='Dust.dat',\
                                                 path=self.path_cc))\
                                                 .find(new_short) != -1:
                    i+=1    
                    new_short = new_short + str(i)
                adding_line = [new_short] + \
                              [str(Star.getInputData(keyword=key,\
                                                     filename='Dust.dat',\
                                                     path=self.path_cc)\
                                    [self.index])
                               for key in ['SPEC_DENS','T_DES','T_DESA','T_DESB']]
                adding_line.insert(2,output_filename)
                adding_line = '\t\t'.join(adding_line)
                DataIO.writeFile(os.path.join(self.path_cc,'Dust.dat'),\
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

        spec_dens = Star.getInputData(keyword='SPEC_DENS',filename='Dust.dat',\
                                      path=self.path_cc)[self.index]
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
        
        spec_dens = Star.getInputData(keyword='SPEC_DENS',filename='Dust.dat',\
                                      path=self.path_cc)[self.index]
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
         
            