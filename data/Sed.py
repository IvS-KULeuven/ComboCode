# -*- coding: utf-8 -*-

"""
Toolbox for reading spectra and spectrum manipulation.

Author: R. Lombaert

"""

import os
from scipy import array, hstack
from scipy.optimize import leastsq
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from glob import glob
import operator

from cc.tools.numerical import Interpol
from cc.tools.io import DataIO
from cc.plotting import Plotting2



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
    
    def __init__(self,star_name,path,plot_extrapol_extinction=0,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):
        
        ''' 
        Initializing an Sed instance. 
        
        Setting starting parameters from a star object.
        
        @param star_name: The star name of the object
        @type star_name: string
        @param path: The path to the folder containing the SED data
        @type path: string
        
        @keyword path_combocode: path to the combocode folder
        
                               (default: ~/ComboCode/)
        @type path_combocode: string
        @keyword plot_extrapol_extinction: Plot and show the result of the 
                                           extrapolated interstellar extinction
                                           law by chiar and tielens (2006)
                                           
                                           (default: 0)
        @type plot_extrapol_extinction: bool
        
        '''
        
        self.star_name = star_name
        self.path = path
        self.path_combocode = path_combocode
        self.plot_extrapol_extinction = plot_extrapol_extinction
        self.data = dict()
        self.data_raw = dict()
        self.setStarPars()
        self.setData()
        self.readData()
        self.dereddenData()



    def setStarPars(self):
        
        """
        Set some standard stellar parameters such as Ak and galactic position.
        
        """
        
        cc_path = os.path.join(self.path_combocode,'Data')
        star_index = DataIO.getInputData(cc_path).index(self.star_name)
        self.ak = DataIO.getInputData(path=cc_path,keyword='A_K')[star_index]
        longitude = DataIO.getInputData(path=cc_path,keyword='LONG')[star_index]
        latitude = DataIO.getInputData(path=cc_path,keyword='LAT')[star_index]
        if (abs(longitude) < 5.0 or longitude > 355.0) and abs(latitude) < 5.0:
            self.gal_position = 'GC'
        else:
            self.gal_position = 'ISM'    
        self.star_index = DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'))\
                              .index(self.star_name)
        self.star_name_plots = DataIO.getInputData(\
                                path=os.path.join(self.path_combocode,'Data'),\
                                keyword='STAR_NAME_PLOTS',remove_underscore=1)\
                               [self.star_index]
    
       
       
    def setData(self):
        
        '''
        Select available data.
        
        Based on the data file types in Dust.dat and the available data files.
        
        '''
        
        cc_path = os.path.join(self.path_combocode,'Data')
        all_data_types = DataIO.getInputData(path=cc_path,keyword='DATA_TYPES',\
                                             filename='Dust.dat')
        searchpath = os.path.join(self.path,'*_%s.dat'%self.star_name)
        data_list = [os.path.split(f)[1].replace('_%s.dat'%self.star_name,'') 
                     for f in glob(searchpath)] 
        self.data_types = [dt for dt in data_list if dt in all_data_types]
            
        
    
    def readData(self):
        
        '''
        Read the raw SED data. 
        
        '''
        
        for dt in self.data_types:
            inputfile = os.path.join(self.path,'%s_%s.dat'%(dt,self.star_name))
            data = DataIO.readCols(inputfile,nans=1)
            data_sorted = array([(x,y) 
                                 for x,y in sorted(zip(data[0],data[1]),\
                                                   key=operator.itemgetter(0)) 
                                 if y])
            self.data_raw[dt] = (data_sorted[:,0],data_sorted[:,1])
            
        print '** LWS is not being aligned automatically with SWS for now.'
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
        
        '''
    
        #- Read the extinction curve, then inter/extrapolate it
        ext_x,ext_y = getExtinctionCurve(self.gal_position,'chiar_tielens',\
                                         0.112)
        #- Initial param guess for extrapolation of long wavelength extinction curve
        p0 = [ -2,  0.01 ,  1, -1]
        #- Assuming a power law for the extrapolation, and fitting from 22 mic  
        deredfunc = 'power'
        extrapol_xmin = 22.0
        #- Fit the power law to the > 22 micron wavelength range
        plsq = leastsq(Interpol.getResiduals,p0,\
                       args=(ext_x[ext_x>=extrapol_xmin],\
                             ext_y[ext_x>=extrapol_xmin],\
                             deredfunc),maxfev=20000)[0]
        #- Calculate the extrapolation and interpolation for the datagrids
        #- Then combine and apply the correction to the data
        for dt,(data_x, data_y) in self.data_raw.items():
            extra = Interpol.pEval(data_x[data_x>=extrapol_xmin],plsq,deredfunc)
            inter = interp1d(ext_x,ext_y)(data_x[data_x<extrapol_xmin])
            corr = hstack([inter,extra])
            if self.plot_extrapol_extinction: 
                Plotting2.plotCols(x=[ext_x,data_x],y=[ext_y,corr],\
                                   xlogscale=1,ylogscale=1)
            self.data[dt] = (data_x,data_y*10**(corr*self.ak*0.4))
                    
        

def getExtinctionCurve(gal_position='ism',curve_type='chiar_tielens',\
                       av_to_ak_conv=0.112,\
                       path=os.path.join(os.path.expanduser('~'),'MCMax',\
                                         'Data','Extinction_Curves')):
    
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
    @keyword path: path for extinction curve data
    
                   (default: '~/MCMax/extinction_curves/')
    @type path: string
    
    @return: The wavelength grid and extinction curve values
    @rtype: (array,array)
    
    """
    
    extcurve_input = DataIO.readFile(os.path.join(path,curve_type + '.dat'),\
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
        #- Put values to Alambda/Ak
        extcurve_y = array(extcurve_y)/(3.1*av_to_ak_conv)       
    else:
        extcurve_y = [float(row[1]) 
                      for row in extcurve_input 
                      if row[0][0] != '#']
    
    extcurve = [[x,y] for x,y in zip(extcurve_x,extcurve_y)]

    return (array(extcurve_x),array(extcurve_y))
