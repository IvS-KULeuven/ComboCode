# -*- coding: utf-8 -*-

"""
Tools for trend analysis.

Author: R. Lombaert

"""

from cc.modeling.objects import Transition
from cc.plotting import Plotting2
import os
import operator
from numpy import array
import numpy as np
import types
from numpy.random import normal
from matplotlib import pyplot as plt


def makeDiagnosticPlot(sg,molec,scaling=[],escaling=[],combine_water=0,\
                       edists=[]):
    
    '''
    Make a diagnostic plot for a series of stars, where Imb for transitions of 
    a given molecule are compared with their upper level energies.
    
    Line strengths are always scaled with distance squared. Additional scaling
    can be requested.
    
    Three plots are made: One with scaling versus distance and two with scaling
    versus distance and the CO line strength of the J=15-14 and J=30-29 lines
    respectively. The comparison with CO line strengths is not scaled except
    with distance.
    
    @param sg: The stellar models, in which the transitions have been matched
               with integrated line strengths. 
    @type sg: list[Star()]
    @param molec: The molecule for which this is done, shorthand notation.
    @type molec: str
    
    @keyword scaling: Scale the line strengths also with a keyword given here.
                      'MDOT_GAS', etc. Assumes a Star() object knows the 
                      keyword.
                      
                      (default:[])
    @type scaling: list[string]
    @keyword combine_water: Combine both ortho and para water in a single plot
    
                            (default: 0)
    @type combine_water: bool
    @keyword edists: Include errors for distance estimates of stars here. 
    
                    (default: [])
    @type edists: list[float]
    @keyword escaling: Include relative errors of extra scaling parameters 
                       as lists. len is len of scaling, len of an element
                       is len of sg
    
                       (default: [])
    @type escaling: list[list]
    
    '''
    
    if '1H1H16O' not in molec: combine_water = 0
    if len(scaling) != len(escaling):
        print 'No errors on the scaling factors taken into account.'
        escaling = []
    if len(edists) != len(sg):
        print 'No errors on the distance taken into account.'
        edists = [0]*len(sg)
    allints = []
    allenergies = []
    allerrors = []
    allwaves = []
    allints_co15 = []
    allints_co30 = []
    allerrs_co15 = []
    allerrs_co30 = []
    
    #-- Select all CO line strengths for two typical transitions: 15-14, 30-29
    co = sg[0].getMolecule('12C16O')
    co1514 = Transition.Transition(molecule=co,telescope='PACS',jup=15,\
                                   jlow=14,path_gastronoom='codeJun2013')
    co3029 = Transition.Transition(molecule=co,telescope='PACS',jup=30,\
                                   jlow=29,path_gastronoom='codeJun2013')
    trl_co15 = Transition.getTransFromStarGrid(sg,co1514,'sample')
    trl_co30 = Transition.getTransFromStarGrid(sg,co3029,'sample')
    ls_co15,errs_co15 = Transition.getLineStrengths(trl_co15,'dint')
    ls_co30,errs_co30 = Transition.getLineStrengths(trl_co30,'dint')
    
    for istar,(s,intco15,intco30,errco15,errco30) in \
            enumerate(zip(sg,ls_co15,ls_co30,errs_co15,errs_co30)):
        sn = s['STAR_NAME']
        allints.append([])
        allenergies.append([])
        allerrors.append([])
        allwaves.append([])
        allints_co15.append([])
        allints_co30.append([])
        allerrs_co15.append([])
        allerrs_co30.append([])
        if combine_water: \
            trans = s.getTransitions('1H1H16O') + s.getTransitions('p1H1H16O')
        else:
            trans = s.getTransitions(molec)
        
        ls_trans,errls_trans = Transition.getLineStrengths(trans,'dint')
        for t,lst,elst in zip(trans,ls_trans,errls_trans):
            if not np.isnan(lst):
                allints_co15[-1].append(abs(lst/intco15))
                allints_co30[-1].append(abs(lst/intco30))
                allerrs_co15[-1].append(np.sqrt(sum([elst**2]+[errco15**2]))\
                                        *abs(lst/intco15))
                allerrs_co30[-1].append(np.sqrt(sum([elst**2]+[errco30**2]))\
                                        *abs(lst/intco30))
                
                tint = abs(lst*s['DISTANCE']**2/100.**2 )
                for par in scaling:
                    tint *= 1./s[par]
                totalerr = np.sqrt(sum([elst**2]+\
                                    [2*edists[istar]**2]+\
                                    [esca[istar]**2 for esca in escaling]))
                allints[-1].append(tint)
                allenergies[-1].append(t.getEnergyUpper())
                allerrors[-1].append(totalerr*tint)
                allwaves[-1].append(t.wavelength*10**4)
            else:
                for tlist in [allints_co15,allints_co30,allerrs_co15,\
                              allerrs_co30,allints,allerrors]:
                    tlist[-1].append(float('nan'))
                allenergies[-1].append(float(t.getEnergyUpper()))
                allwaves[-1].append(t.wavelength*10**4)
        isort = np.argsort(allenergies[-1])
        for tlist in [allints_co15,allints_co30,allerrs_co15,allwaves,\
                      allerrs_co30,allints,allerrors,allenergies]:
            tlist[-1] = array(tlist[-1])[isort]
    
    isort = np.argsort([s['MDOT_GAS'] for s in sg])
    for tlist in [allints_co15,allints_co30,allerrs_co15,allwaves,\
                  allerrs_co30,allints,allerrors,allenergies,sg]:
        tlist = array(tlist)[isort]
    
    pfn_path = os.path.join(os.path.expanduser('~'),'Projects','Cstars_h2o',\
                            'DiagnosticPlots')
    linestyles = ['o-','o-','o-','o-','o-','o-','o-',\
                  '--x','--x','--x','--x','--x','--x','--x',\
                  '-.s','-.s','-.s','-.s','-.s','-.s','-.s']
    colors = ['r','b','k','g','m','y','c']
    line_types = [ls + col for ls,col in zip(linestyles,3*colors)]
    line_types = line_types[:len(sg)]

    xmin = 0
    xmax = len(allenergies[0])+1
    plot_title = ', '.join(['%i: %.1f cm$^{-1}$ - %.1f $\mu$m'\
                             %(i+1,t.getEnergyUpper(),t.wavelength*10**4) 
                            for i,t in enumerate(sg[0].getTransitions(molec))])
    if scaling:
        plot_title += 'Extra scaling: ' + ', '.join([sca.replace('_','\_') 
                                                     for sca in scaling])
    x = []
    y = []
    yerr = []
    y_co15 = []
    ye15 = []
    y_co30 = []
    ye30 = []
    for istar in range(len(sg)):
        x.append([])
        y.append([])
        yerr.append([])
        y_co15.append([])
        y_co30.append([])
        ye15.append([])
        ye30.append([])
        for i,(iint,eint,iint_co15,iint_co30,eint_co15,eint_co30) in \
                enumerate(zip(allints[istar],allerrors[istar],\
                              allints_co15[istar],allints_co30[istar],\
                              allerrs_co15[istar],allerrs_co30[istar])):
            if not np.isnan(iint):
                x[-1].append(i+1)
                y[-1].append(iint)
                yerr[-1].append(eint)
            if not np.isnan(iint_co15):
                y_co15[-1].append(iint_co15)
                ye15[-1].append(eint_co15)
            if not np.isnan(iint_co30):
                y_co30[-1].append(iint_co30)
                ye30[-1].append(eint_co30)
    extension = 'pdf'
    
    path = os.path.join(pfn_path,'ints_vs_eul_%s'%molec)
    for par in scaling:
        if par == scaling[0]:
            path += '_extrascaling'
        path += '_%s'%par
    print Plotting2.plotCols(filename=path,\
                             x=x,y=y,extension=extension,\
                             yerr=yerr,\
                             xmin=xmin,xmax=xmax,plot_title=plot_title,\
                             yaxis='$I_\mathrm{int}$ (W/m$^2$)',\
                             xaxis='Index Energy Upper Level',\
                             line_types=line_types,\
                             keytags=['%.1e -- %s'%(s['MDOT_GAS'],s['STAR_NAME']) for s in sg],\
                             key_location=(0.87,0.01),ylogscale=1,\
                             linewidth=3,fontsize_key=26,fontsize_title=20)
    
    if molec != '12C16O' and not scaling:
        path = os.path.join(pfn_path,'ints_co15_vs_eul_%s'%molec)
        print Plotting2.plotCols(filename=path,\
                                 x=x,y=y_co15,extension=extension,\
                                 yerr=ye15,\
                                 xmin=xmin,xmax=xmax,plot_title=plot_title,\
                                 yaxis='$I_\mathrm{int}/I_\mathrm{CO 15-14}$ (W/m$^2$)',\
                                 xaxis='Index Energy Upper Level',\
                                 line_types=line_types,\
                                 keytags=['%.1e -- %s'%(s['MDOT_GAS'],s['STAR_NAME']) for s in sg],\
                                 key_location=(0.87,0.01),ylogscale=1,\
                                 linewidth=3,fontsize_key=26,fontsize_title=20)
                        
        path = os.path.join(pfn_path,'ints_co30_vs_eul_%s'%molec)
        print Plotting2.plotCols(filename=path,\
                                 x=x,y=y_co30,extension=extension,\
                                 yerr=ye30,\
                                 xmin=xmin,xmax=xmax,plot_title=plot_title,\
                                 yaxis='$I_\mathrm{int}/I_\mathrm{CO 30-29}$ (W/m$^2$)',\
                                 xaxis='Index Energy Upper Level',\
                                 line_types=line_types,\
                                 keytags=['%.1e -- %s'%(s['MDOT_GAS'],s['STAR_NAME']) for s in sg],\
                                 key_location=(0.87,0.01),ylogscale=1,\
                                 linewidth=3,fontsize_key=26,fontsize_title=20)


    
def makeParamPlot(sg,xpar,ypar,expar=[],eypar=[],xratios=[],yratios=[],\
                  emdot=[],exparlog=0,eyparlog=0,edists=[],mode='dint',\
                  n_data=0,extra_mcon=[],extra_dcon=[],cfg='',\
                  add_linear_fit=0,alf_xmin=None,alf_xmax=None,efcont63=[],\
                  **kwargs):
    
    '''
    Make a diagnostic plot of either measured line strengths or intrinsic 
    parameters versus measured line strengths or intrinsic parameters. 
    
    Ratios are possible for line strengths. Not for intrinsic parameters.
    
    Requires preparatory work done for the Pacs() and the Star() objects.
    
    @param sg: The stellar models, in which the transitions have been matched
               with integrated line strengths. If both models and data are 
               combined, the data Star() objects are assumed to be listed 
               first.
    @type sg: list[Star()]
    @param xpar: The parameter on the x-axis. Can be either a string (Star() 
                 keyword), or an index (of the transition in the first object 
                 in the sg list) for line strengths. In a combo mode (cint or 
                 ctmb) this means it is the index in the transition list of the
                 data objects rather than the model objects. Transitions in 
                 objects other than the first 1 can have different indices.
    @type xpar: string/int
    @param ypar: The parameter on the y-axis. Can be either a string (Star() 
                 keyword), or an index (of the transition in the first object 
                 in the sg list) for line strengths. In a combo mode (cint or 
                 ctmb) this means it is the index in the transition list of the
                 data objects rather than the model objects. Transitions in 
                 objects other than the first 1 can have different indices.
    @type ypar: string/int
    
    @keyword xratios: If xpar is a line strength, multiple ratios can be 
                      requested to be plotted in succession. 
                      Therefore, this gives the indices (if int, refers to the 
                      1st Star() object in sg) or 'mdot' (if ratio wrt Mdot)  
                      or 'fcont63' (if ratio wrt 6.3 mic continuumflux) for the
                      x-axis ratio.
                      
                      (default: [])
    @type xratios: list[int/str]
    @keyword yratios: If ypar is a line strength, multiple ratios can be 
                      requested to be plotted in succession.
                      Therefore, this gives the indices (if int, refers to the 
                      1st Star() object in sg) or 'mdot' (if ratio wrt Mdot)  
                      or 'fcont63' (if ratio wrt 6.3 mic continuumflux) for the
                      x-axis ratio.
                                            
                      (default: [])
    @type yratios: list[int/str]
    @keyword emdot: Include errors for the x/yratio quantity if it is mdot. Not 
                    used for error estimation on mdot as a parameter! The mdot 
                    errors are given in log scale.
    
                    (default: [])
    @type emdot: list[float]
    @keyword efcont63: Include errors for the x/yratio quantity if it is 
                       fcont63. Not used for error estimation on fcont63 as a 
                       parameter! The fcont63 errors are given in linear scale
                       
                       (default: [])
    @type efcont63: list[float]
    @keyword expar: The error on the x-parameter if it is a Star() key and if 
                    mode is cint or dint. Number of entries in array is equal 
                    to the number of data Star() objects.
    
                    (default: [])
    @type expar: array
    @keyword eypar: The error on the y-parameter if it is a Star() key and if 
                    mode is cint or dint. Number of entries in array is equal 
                    to the number of data Star() objects.
    
                    (default: [])
    @type eypar: array
    @keyword exparlog: The xpar error is given in logscale. Only relevant for
                       the d and c modes.
                                  
                       (default: 0)
    @type exparlog: bool
    @keyword eyparlog: The ypar error is given in logscale. Only relevant for
                       the d and c modes.
                                  
                       (default: 0)
    @type eyparlog: bool
    @keyword edists: Include errors for distance estimates of stars here. These
                     distances are only used to rescale line strengths if they 
                     are not in a ratio.
    
                     (default: [])
    @type edists: list[float]
    @keyword mode: The mode in which line strengths are selected, ie either 
                   from data or from models. Either 'dint', 'mint', 'mtmb' or 
                   'dtmb' values. A combination of both is possible by setting
                   this key to 'cint' or 'ctmb'. Then the extra keyword 
                   'n_data' is required, which indicates how many Star() 
                   objects are associated with data. The method assumes they 
                   are the first objects in the list of Star() objects.
                   
                   (default: 'dint')
    @type mode: str
    @keyword n_data: The number of data Star() objects, assuming they are the 
                     first in the star_grid. Only required if mode == 'combo'.
                     
                     (default: 0)
    @type n_data: int
    @keyword extra_mcon: If extra conditionals are requested for models, the 
                         plot is colour coded based on them. For instance, 
                         Star() object keywords can serve as conditionals.
                         Note that these are only applied when mode == mtmb, 
                         mint, cint or ctmb.
                        
                         (default: [])
    @type extra_mcon: list[string]
    @keyword extra_dcon: If extra conditionals are requested for data, the plot  
                         is colour coded based on them. For instance, Star() 
                         object keywords can serve as conditionals. Note that 
                         these are only applied when mode == dtmb, dint, cint 
                         or ctmb.
                                
                         (default: [])
    @type extra_dcon: list[string]
    @keyword cfg: config filename read as a dictionary, can replace any keyword 
                  given to plotCols. Can also be a dictionary itself, in which
                  case no file is read and the kwargs are updated with the 
                  content of cfg
                  
                  (default: '')
    @type cfg: string/dict
    @keyword add_linear_fit: Add a linear fit to the figures. The fit is done 
                             through corrSG method, of which extra arguments 
                             can be given in kwargs. (xpar_co, ypar_co)
                             Only works in dint or cint mode if xratios or 
                             yratios has len less than 2.
                             
                             (default: 0)
    @type add_linear_fit: bool
    @keyword alf_xmin: The minimum x value for the linear fit plot, if 
                       requested. (This is not the cut off value for the 
                       fitting routine itself!) Has to be given if a linear fit
                       is requested.
                       
                       (default: None)
    @type alf_xmin: float
    @keyword alf_xmax: The maximum x value for the linear fit plot, if 
                       requested. (This is not the cut off value for the 
                       fitting routine itself!) Has to be given if a linear fit
                       is requested.
                       
                       (default: None)
    @type alf_xmax: float
     
    @keyword **kwargs: extra keywords needed for the linear fit, if requested.
    @type **kwargs: dict
    
    @return: The filename of the produced plot is returned. 
    @rtype: str
    
    '''
    
    x_titles = dict([('MDOT_GAS',r'$\log$ $\left[\dot{M}_\mathrm{g}\ (\mathrm{M}_\odot/\mathrm{yr})\right]$'),\
                     ('MDOT_DUST',r'$\log$ $\left[\dot{M}_\mathrm{d}\ (\mathrm{M}_\odot/\mathrm{yr})\right]$'),\
                     ('VEL_INFINITY_GAS','$v_{\infty\mathrm{,g}}$ ($\mathrm{km} \mathrm{s}^{-1}$)'),\
                     ('SHELLMASS',r'$\log$ $\left[\bar{M_\mathrm{s}}\ (\mathrm{g}\ \mathrm{cm}^{-1})\right]$'),\
                     ('SHELLDENS',r'$\log$ $\left[\bar{\rho}\ (\mathrm{g}\ \mathrm{cm}^{-3})\right]$'),\
                     ('SHELLCOLDENS',r'$\log$ $\left[\bar{m}\ (\mathrm{g}\ \mathrm{cm}^{-2})\right]$'),\
                     ('SHELLDENS2',r'$\sqrt{\bar{\rho}^2 R_\star}$ (g/cm$^{5/2}$)'),\
                     ('L_STAR','$L_\star$ (L$_\odot$)'),\
                     ('P_STAR',r'$\log$ $\left[P\ (\mathrm{days})\right]$'),\
                     ('T_STAR','$T_\star$ (K)'),\
                     ('Q_STAR','$Q_\star$ (days)'),\
                     ('R_INNER_GAS','$R_\mathrm{i,g}$ (R$_\star$)'),\
                     ('F_H2O',r'$\log$ $\left[A_{\mathrm{H}_2\mathrm{O}}/A_{\mathrm{H}_2}\right]$'),\
                     ('F_CONT_63',r'$\log$ $\left[F_\mathrm{6.3\ \mu m}\ (\mathrm{Jy})\right]$'),\
                     ('FD2_CONT_63',r'$\log$ $\left[F_\mathrm{6.3\ \mu m}\times D^2\ (\mathrm{Jy}\ \mathrm{pc}^2)\right]$'),\
                     ('FD2M_CONT_63',r'$\log$ $\left[F_\mathrm{6.3\ \mu m}\times D^2 / \dot{M}_\mathrm{g}\ (\mathrm{Jy}\ \mathrm{pc}^2\ yr/\mathrm{M}_\odot \right]$)'),\
                     ])
    pfn_parts = dict([('MDOT_GAS','mg'),\
                      ('MDOT_DUST',r'md'),\
                      ('STARTYPE','startype'),\
                      ('A_SICB','asicb'),\
                      ('A_AMCSPH','aamcsph'),\
                      ('VEL_INFINITY_GAS','vg'),\
                      ('SHELLMASS','shellmass'),\
                      ('SHELLDENS','dens'),\
                      ('SHELLCOLDENS','coldens'),\
                      ('SHELLDENS2','dens3-2'),\
                      ('L_STAR','lstar'),\
                      ('P_STAR','period'),\
                      ('Q_STAR','qstar'),\
                      ('T_STAR','tstar'),\
                      ('P_TYPE','ptype'),\
                      ('F_H2O','ah2o'),\
                      ('DUST_TO_GAS_CHANGE_ML_SP','d2g'),\
                      ('TEMPERATURE_EPSILON_GAS','eps1'),\
                      ('TEMPERATURE_EPSILON2_GAS','eps2'),\
                      ('TEMPERATURE_EPSILON3_GAS','eps3'),\
                      ('RADIUS_EPSILON2_GAS','rt12'),\
                      ('RADIUS_EPSILON3_GAS','rt23'),\
                      ('R_INNER_GAS','rig'),\
                      ('MDOT_CLASS','mdotgrad'),\
                      ('SCD_CLASS','scdgrad'),\
                      ('ABUN_O','abuno'),\
                      ('L_CLASS','lclass'),\
                      ('T_CLASS','tclass'),\
                      ('VG_CLASS','vgclass'),\
                      ('F_CONT_63','fcont63'),\
                      ('FD2_CONT_63','fd2cont63'),\
                      ('F_CONT_63_TYPE','fcont63type'),\
                      ('FD2M_CONT_63','fd2mcont63'),\
                      ('ENHANCE_ABUNDANCE_FACTOR_H2O','h2oabunfac'),\
                      ('ABUNDANCE_FILENAME_H2O','h2oabunfile')])
    keynames = dict([('MDOT_GAS','$\dot{M}_\mathrm{g}$'),\
                     ('MDOT_DUST',r'$\dot{M}_\mathrm{d}$'),\
                     ('A_SICB','A(SICB)'),\
                     ('STARTYPE','StarType'),\
                     ('A_AMCSPH','A(AMCSPH)'),\
                     ('VEL_INFINITY_GAS','$v_{\infty\mathrm{,g}}$'),\
                     ('SHELLMASS','$\bar{M_\mathrm{s}}$'),\
                     ('SHELLDENS',r'$\bar{\rho}$'),\
                     ('SHELLCOLDENS',r'$\bar{m}$'),\
                     ('SHELLDENS2',r'$\sqrt{\bar{\rho}^2 R_\star}}$'),\
                     ('L_STAR','$L_\star$'),\
                     ('P_STAR','$P$'),\
                     ('Q_STAR','$Q_\star$'),\
                     ('T_STAR','$T_\star$'),\
                     ('P_TYPE','Var.~Type'),\
                     ('F_H2O','$n_{\mathrm{H}_2\mathrm{O}}/n_{\mathrm{H}_2}$'),\
                     ('DUST_TO_GAS_CHANGE_ML_SP','$\psi$'),\
                     ('TEMPERATURE_EPSILON_GAS',r'$\epsilon$'),\
                     ('TEMPERATURE_EPSILON2_GAS',r'$\epsilon_2$'),\
                     ('TEMPERATURE_EPSILON3_GAS',r'$\epsilon_3$'),\
                     ('RADIUS_EPSILON2_GAS','$R_\mathrm{T, 12}$'),\
                     ('RADIUS_EPSILON3_GAS','$R_\mathrm{T, 23}$'),\
                     ('R_INNER_GAS','$R_\mathrm{i,g}$'),\
                     ('MDOT_CLASS',''),\
                     ('SCD_CLASS',''),\
                     ('ABUN_O','$n_{\mathrm{O}}/n_{\mathrm{H}_\mathrm{tot}}$'),\
                     ('L_CLASS',''),\
                     ('T_CLASS',''),\
                     ('VG_CLASS',''),\
                     ('F_CONT_63',r'$F_\mathrm{6.3\ \mu m}$'),\
                     ('FD2_CONT_63',r'$F_\mathrm{6.3\ \mu m}\times D^2$'),\
                     ('F_CONT_63_TYPE','$Type F_\mathrm{6.3\ \mu m}$'),\
                     ('FD2M_CONT_63',r'$F_\mathrm{6.3\ \mu m}\times D^2 / \dot{M}_\mathrm{g}$'),\
                     ('ENHANCE_ABUNDANCE_FACTOR_H2O','h2oAbunFac'),\
                     ('ABUNDANCE_FILENAME_H2O','h2oAbunFile')])
    keyunits = dict([('MDOT_GAS','$\mathrm{M}_\odot\ \mathrm{yr}^{-1}$'),\
                     ('MDOT_DUST','$\mathrm{M}_\odot\ \mathrm{yr}^{-1}$'),\
                     ('STARTYPE',''),\
                     ('A_SICB',''),\
                     ('A_AMCSPH',''),\
                     ('VEL_INFINITY_GAS','$\mathrm{km\;s}^{-1}$'),\
                     ('SHELLMASS','$\mathrm{g\;cm}^{-1}$'),\
                     ('SHELLDENS','$\mathrm{g\;cm}^{-3}$'),\
                     ('SHELLCOLDENS','$\mathrm{g\;cm}^{-2}$'),\
                     ('SHELLDENS2','$\mathrm{g\;cm}^{5/2}$'),\
                     ('L_STAR','$\mathrm{L}_\odot$'),\
                     ('P_STAR','$\mathrm{days}$'),\
                     ('Q_STAR','$\mathrm{days}$'),\
                     ('T_STAR','$\mathrm{K}$'),\
                     ('P_TYPE',''),\
                     ('F_H2O',''),\
                     ('DUST_TO_GAS_CHANGE_ML_SP',''),\
                     ('TEMPERATURE_EPSILON_GAS',''),\
                     ('TEMPERATURE_EPSILON2_GAS',''),\
                     ('TEMPERATURE_EPSILON3_GAS',''),\
                     ('RADIUS_EPSILON2_GAS','$\mathrm{R}_\star$'),\
                     ('RADIUS_EPSILON3_GAS','$\mathrm{R}_\star$'),\
                     ('R_INNER_GAS','$\mathrm{R}_\star$'),\
                     ('MDOT_CLASS',''),\
                     ('SCD_CLASS',''),\
                     ('ABUN_O',''),\
                     ('T_CLASS',''),\
                     ('VG_CLASS',''),\
                     ('L_CLASS',''),\
                     ('F_CONT_63','$\mathrm{Jy}$'),\
                     ('FD2_CONT_63','$\mathrm{Jy}$ $\mathrm{pc}^2$'),\
                     ('F_CONT_63_TYPE',''),\
                     ('FD2M_CONT_63',r'$\mathrm{Jy}$ $\mathrm{pc}^2$ $\mathrm{yr}\ \mathrm{M}_{\odot}^{-1}$'),\
                     ('ENHANCE_ABUNDANCE_FACTOR_H2O',''),\
                     ('ABUNDANCE_FILENAME_H2O','')])
    makeints = dict([('MDOT_GAS',0),\
                     ('MDOT_DUST',0),\
                     ('STARTYPE',0),\
                     ('A_SICB',0),\
                     ('A_AMCSPH',0),\
                     ('VEL_INFINITY_GAS',1),\
                     ('SHELLMASS',0),\
                     ('SHELLDENS',0),\
                     ('SHELLCOLDENS',0),\
                     ('SHELLDENS2',0),\
                     ('L_STAR',1),\
                     ('P_STAR',1),\
                     ('Q_STAR',0),\
                     ('T_STAR',1),\
                     ('P_TYPE',0),\
                     ('F_H2O',0),\
                     ('DUST_TO_GAS_CHANGE_ML_SP',0),\
                     ('TEMPERATURE_EPSILON_GAS',0),\
                     ('TEMPERATURE_EPSILON2_GAS',0),\
                     ('TEMPERATURE_EPSILON3_GAS',0),\
                     ('RADIUS_EPSILON2_GAS',0),\
                     ('RADIUS_EPSILON3_GAS',0),\
                     ('R_INNER_GAS',0),\
                     ('MDOT_CLASS',0),\
                     ('SCD_CLASS',0),\
                     ('ABUN_O',0),\
                     ('L_CLASS',0),\
                     ('T_CLASS',0),\
                     ('VG_CLASS',0),\
                     ('F_CONT_63',0),\
                     ('FD2_CONT_63',0),\
                     ('F_CONT_63_TYPE',0),\
                     ('FD2M_CONT_63',0),\
                     ('ENHANCE_ABUNDANCE_FACTOR_H2O',0),\
                     ('ABUNDANCE_FILENAME_H2O',0)])
    
    edists,emdot,efcont63 = array(edists), array(emdot), array(efcont63)    
    n_data = int(n_data)
    expar, eypar = array(expar), array(eypar)
    if type(extra_mcon) is types.StringType:
        extra_mcon = [extra_mcon]
    if type(extra_dcon) is types.StringType:
        extra_dcon = [extra_dcon]
    if type(xpar) is types.StringType:
        xratios = []
    if type(ypar) is types.StringType:
        yratios = []
    ratios = [xratios,yratios]
    pfn_path = os.path.join(os.path.expanduser('~'),'Projects','Cstars_h2o',\
                            'DiagnosticPlots')
    sg_dists = array([s['DISTANCE'] for s in sg])
    sg_mdot = array([s['MDOT_GAS'] for s in sg])
        
    #-- Collect x and y information to simplify coding later on, as all data
    #   collection and error handling is the same for x and y.
    pars = [xpar,ypar]
    epars = [expar,eypar]
    eparlogs = [exparlog,eyparlog]
    
    if mode[0] == 'm':
        #-- In model mode, no errors can be given for any of the parameters.
        add_linear_fit = 0
        expar = array([])
        eypar = array([])
        n_data = 0
        extra_dcon = []
        extra_con = extra_mcon
        current_con = 'm'
    elif mode[0] == 'd':
        n_data = len(sg)
        extra_mcon = []
    if mode[0] != 'm':
        extra_con = extra_dcon
        current_con = 'd'
    for istar,s in enumerate(sg):
        if n_data and istar == n_data:
            extra_con = extra_mcon
            current_con = 'm'
        s['EC'] = (current_con,\
                    tuple([s[con]
                            for con in extra_con 
                            if s[con] <> None]))
    ecl = sorted(list(set([s['EC'] for s in sg])))
    ecl_num = []
    for ec in ecl:
        isg = array([s['EC'] == ec for s in sg])
        isgd = isg[:n_data]
        isgm = isg[n_data:]
        nsgd = len(isg[:n_data][isgd])
        nsgm = len(isg[n_data:][isgm])
        ecl_num.append((ec,isg,isgd,isgm,nsgd,nsgm))
    
    if len(xratios) > 1 or len(yratios) > 1:
        #-- Linear fits can only be added if only one ratio is requested.
        add_linear_fit = 0
    if add_linear_fit:
        add_linear_fit = dict()
        xratio, yratio = None, None
        if xratios: 
            xratio = xratios[0]
        if yratios:
            yratio = yratios[0]
        if yratio == 'fcont63':
            eyratio = efcont63
        else:
            eyratio = []
        if xratio == 'fcont63':
            exratio = efcont63
        else:
            exratio = []
        results = corrSG(sg=sg[:n_data],xpar=xpar,ypar=ypar,expar=expar,\
                         eypar=eypar,xratio=xratio,yratio=yratio,edist=edists,\
                         show=0,eyratio=eyratio,exratio=exratio,**kwargs)
        fitcoef = results['results']
        xgrid = []
        ygrid = []
        this_x = array([alf_xmin,alf_xmax])
        add_linear_fit['xmean'] = this_x
        add_linear_fit['ymean'] = results['intercept']+results['slope']*this_x
        for n in range(0, fitcoef.shape[0],4):
            this_y = fitcoef[n,1] + fitcoef[n,0] * this_x
            xgrid.append(this_x)
            ygrid.append(this_y)
        add_linear_fit['xgrid'] = xgrid
        add_linear_fit['ygrid'] = ygrid

    #-- Check if 'mdot' is requested. Split these up in y and x mdot, because
    #   the error estimate on the ratio depends on the x or y line strength.
    #   This is not the case for fcont63 because the error estimate is simpler
    #   there, hence either y or x fcont63 ratios can remain the same.
    if 'mdot' in xratios: 
        xratios[xratios.index('mdot')] = 'xmdot'
    if 'mdot' in yratios: 
        yratios[yratios.index('mdot')] = 'ymdot'
        
    #-- Select all line strengths and errors for the ratios. fcont63 if 
    #   requested is added as well. These dicts hold the info for both x and y.
    ls_ratios = dict()
    els_ratios = dict()
    for i in set(yratios+xratios): 
        #-- mdot must be done separately, due to the cumbersome error estimate
        if i == 'xmdot' or i == 'ymdot': continue
        elif i == 'fcont63':
            #-- Convert to W/m2/Hz for unit consistency later on 
            #   (LS/fcont63 is in Hz)
            sg_fcont63 = array([s['F_CONT_63']*1e-26 for s in sg])
            ls_ratios[i] = sg_fcont63
            els_ratios[i] = efcont63
        else:
            i = int(i)
            ratsample = sg[0]['GAS_LINES'][i]
            rattrans = Transition.getTransFromStarGrid(sg,ratsample,'sample')
            ls,els = Transition.getLineStrengths(rattrans,mode,n_data=n_data)
            ls_ratios[i] = ls
            els_ratios[i] = els
    
    #-- Set the dictionaries for x and y that will hold all to be plotted data
    #   No extra keys are added if x/yratios is empty. (xratios empty by 
    #   default if xpar/ypar is a string).
    x, xerr, xblend = dict([('x',[])]), dict([('x',[])]), dict([('x',[])])
    for k in xratios:
        x[k], xerr[k], xblend[k] = [], [], []
    y, yerr, yblend = dict([('y',[])]), dict([('y',[])]), dict([('y',[])])
    for k in yratios:
        y[k], yerr[k], yblend[k] = [], [], []
    
    blends = [xblend,yblend]
    xy = [x,y]
    errs = [xerr,yerr]
    axes = ['x','y']
    #-- Select the x/y parameters, making a difference between a Star() key or 
    #   line strengths. Star keys are never used in ratios, except MDOT_GAS and
    #    F_CONT_63, but those are handled separately anyway.
    sample, seltrans, allint, allerr = [] , [], [], []
    for par,epar,eparlog,blend,xyi,err,axisstr in zip(pars,epars,eparlogs,\
                                                      blends,xy,errs,axes):
        if type(par) is types.StringType:
            sg_par = array([s[par] for s in sg])
            sample.append(None)
            seltrans.append(None)
            allint.append(None)
            allerr.append(None)
            for (ec,isg,isgd,isgm,nsgd,nsgm) in ecl_num:
                #-- No blends possible for parameters, so add False for all. 
                blend[axisstr].append(np.zeros(nsgd) != 0)
                xyi[axisstr].append(np.log10(sg_par[isg]))
                #-- Check if errors are given (in d or c mode) and check for log
                if epar.size and eparlog:
                    err[axisstr].append(np.concatenate([epar[isgd],np.zeros(nsgm)]))
                elif epar.size:
                    ll = np.concatenate([-np.log10(1-epar[isgd]),np.zeros(nsgm)])
                    ul = np.concatenate([np.log10(1+epar[isgd]),np.zeros(nsgm)])
                    err[axisstr].append([ll,ul])
        else:
            #-- Select the line strengths and errors of the main transition for 
            #   both axes. This information is only used when par is not a str.
            #   Later, when ratios are set, this also is only done if par is 
            #   not a str. So no mix-ups can happen. The sample and selection 
            #   transitions are remembered for later.
            isample = sg[0]['GAS_LINES'][par]
            iseltrans = Transition.getTransFromStarGrid(sg,isample,'sample')
            iallint,iallerr = Transition.getLineStrengths(iseltrans,mode,\
                                                          n_data=n_data)
            sample.append(isample)
            seltrans.append(iseltrans)
            allint.append(iallint)
            allerr.append(iallerr)
                
    #-- If an Mdot ratio is requested, set the second component of the ratio
    #   and the errors here. (Takes a bit of calc time to estimate errors for 
    #   these ratios)
    for iratios,iallint,iallerr,axisstr in zip(ratios,allint,allerr,axes):
        if axisstr+'mdot' in iratios:
            line1 = abs(iallint[:n_data])*sg_dists[:n_data]**2/100.**2
            line1_err = line1*np.sqrt(iallerr[:n_data]**2+4*edists**2)
            line2 = np.log10(sg_mdot[:n_data])
            line2_err = emdot
            mratios = guessRatio(line1,line1_err,line2,line2_err,\
                                line2_log=1,positive=1,n_fit=10000)
            emrat = array([np.std(mratios[:,istar])/np.mean(mratios[:,istar])
                        for istar in range(len(line1))])
            ls_ratios[axisstr+'mdot'] = sg_mdot
            els_ratios[axisstr+'mdot'] = emrat
    
    #-- Complete the x and y data dicts.
    #   And make sure to remember which stars (in the datagrid) go where.
    for (ec,isg,isgd,isgm,nsgd,nsgm) in ecl_num:
        #-- Set the main line strength for the x/yaxis if par is not a string
        for par,blend,xyi,err,axisstr,iallint,iallerr,iratios \
                in zip(pars,blends,xy,errs,axes,allint,allerr,ratios):
            if type(par) is not types.StringType:
                blend[axisstr].append(iallint[isgd] < 0)
                #-- This irat1 is used later as well but only in this for loop
                irat1 = abs(iallint[isg])
                xyi[axisstr].append(np.log10(irat1*sg_dists[isg]**2/100.**2))
                #-- Set the errors in case data are involved.
                #   Note that the error bars in log scale can be calculated without
                #   knowing the real value, since it is a constant added to the real 
                #   value by matplotlib. Unfortunately, going to log space, means the 
                #   upper and lower limits will differ so two lists are needed. Add a 
                #   minus sign to the lower limit, as matplotlib will subtract the ll.
                if mode[0] != 'm':
                    etot = np.sqrt(iallerr[isgd]**2+4*edists[isgd]**2)
                    ll = np.concatenate([-np.log10(1-etot),np.zeros(nsgm)])
                    ul = np.concatenate([np.log10(1+etot),np.zeros(nsgm)])
                    err[axisstr].append([ll,ul])

            for k in iratios:
                if k == axisstr+'mdot':
                    #-- Just append the bool array for the 1st component LS
                    blend[k].append(blend[axisstr][-1])
                    xyi[k].append(xyi[axisstr][-1]-np.log10(ls_ratios[k][isg]))
                    if mode[0] != 'm': 
                        etot = els_ratios[k][isgd]
                elif k =='fcont63':
                    #-- Just append the bool array for the 1st component LS
                    blend[k].append(blend[axisstr][-1])
                    xyi[k].append(np.log10(irat1/ls_ratios[k][isg]))
                    if mode[0] != 'm': 
                        etot = np.sqrt(iallerr[isgd]**2+els_ratios[k][isgd]**2)
                else:        
                    blend[k].append((blend[axisstr][-1])+(ls_ratios[k][isgd]<0))
                    xyi[k].append(np.log10(irat1/abs(ls_ratios[k][isg])))
                    if mode[0] != 'm': 
                        etot = np.sqrt(iallerr[isgd]**2+els_ratios[k][isgd]**2)
                if mode[0] != 'm':
                    ll = np.concatenate([-np.log10(1-etot),np.zeros(nsgm)])
                    ul = np.concatenate([np.log10(1+etot),np.zeros(nsgm)])
                    err[k].append([ll,ul])
    
    #-- Set a number of parameters for the figures.
    cfg_dict = Plotting2.readCfg(cfg)
    extra_pars = dict()
    #-- Set the title, depending on if LS are requested vs pars.
    pt = ''
    if not type(xpar) is types.StringType:
        pt += 'VS %s: E$_\mathrm{ul,x}$ = %.1f - %.2f'\
              %(str(sample[0]),sample[0].getEnergyUpper(),\
                sample[0].wavelength*10**4)
    if not type(ypar) is types.StringType:
        pt += '%s: E$_\mathrm{ul,y}$ = %.1f - %.2f'\
              %(str(sample[1]),sample[1].getEnergyUpper(),\
                sample[1].wavelength*10**4)
    extra_pars['plot_title'] = pt
    extra_pars['fontsize_title'] = 20
    extra_pars['figsize'] = (8*np.sqrt(2),8)
    extra_pars['extension'] = '.pdf'
    extra_pars['fontsize_key'] = 14
    extra_pars['linewidth'] = 2
    
    #-- Set the keytags and linestyles based on if data or models are
    #   plotted.
    if not ecl in [[('m',()),('d',())],[('m',())],[('d',())]]:
        keytags = []
        for currcon,ec in ecl:
            this_con = currcon == 'm' and extra_mcon or extra_dcon
            k = []
            for con,v in zip(this_con,ec):
                if 'CLASS' in con:
                    kstr = v[1]
                elif con == 'P_TYPE':
                    kstr = '$\mathrm{%s}$'%v 
                elif con == 'SHELLCOLDENS':
                    kstr = '%s = $%.2f$ %s'%(keynames[con],v,keyunits[con])
                else:
                    kstr = '%s = $%s$ %s'\
                           %(keynames[con],\
                             makeints[con] and str(int(v)) or str(v),\
                             keyunits[con])
                k.append(kstr)
            keytags.append(', '.join(k))
        extra_pars['key_location'] = 'best'
    mlinestyles = ['-x','-x','-x','-x','-x','-x','-x',\
                   '--s','--s','--s','--s','--s','--s','--s',\
                   '-.+','-.+','-.+','-.+','-.+','-.+','-.+',\
                   '--p','--p','--p','--p','--p','--p','--p',\
                   'o-','o-','o-','o-','o-','o-','o-']
    dlinestyles = ['o','o','o','o','o','o','o',\
                   'x','x','x','x','x','x','x',\
                   's','s','s','s','s','s','s']
    colors = ['r','b','g','k','m','y','c']
    dline_types = [ls + col for ls,col in zip(dlinestyles,3*colors)]
    colors.reverse()
    mline_types = [ls + col for ls,col in zip(mlinestyles,5*colors)]
    if mode[0] == 'm':
        line_types = mline_types[:len(ecl)]
        zorder = range(len(ecl))
    elif mode[0] == 'd':
        line_types = dline_types[:len(ecl)]
        zorder = range(len(ecl))
    else:
        d_ecl = len([ec for ec in ecl if ec[0] == 'd'])
        m_ecl = len([ec for ec in ecl if ec[0] == 'm'])
        line_types = dline_types[:d_ecl] \
                        + mline_types[:m_ecl]
        zorder = range(10,10+d_ecl) + range(-m_ecl,0)
    markersize = [6]*len(keytags)
    pfn_ecl = '_'.join([pfn_parts[ec] for ec in extra_dcon+extra_mcon])
    
    #-- Avoid overhead: If ratios are requested, you generally dont want the 
    #   separate line strengths outside a ratio.
    if xratios:
        del x['x']
    if yratios: 
        del y['y']
        
    #-- Loop over the X-AXIS KEYS
    for xk in x.keys():
        #-- extract the ratio transition if applicable
        if xk not in ['x','xmdot','fcont63']:
            xratsample = sg[0]['GAS_LINES'][xk]
        
        #-- Set the x-axis title, xmin, xmax and pfn_xtag.
        if type(xpar) is types.StringType:
            extra_pars['xaxis'] = x_titles[xpar]
            extra_pars['xmin'] = min([min(xi) for xi in x[xk]])-0.2
            extra_pars['xmax'] = max([max(xi) for xi in x[xk]])+0.2
            pfn_xtag = pfn_parts[xpar]
            pfn_xrat = ''
        else: 
            #-- Adapt the xaxis title based on the xratios.
            s1 = sample[0].makeAxisLabel()
            if xk == 'xmdot':
                extra_pars['xaxis'] = r'$\log$ $\left[%s/\dot{M}_\mathrm{'%s1+\
                                      r'g}\ (\mathrm{W}/\mathrm{m}^2\ \mathrm{yr}/\mathrm{M}_\odot)\right]$'
            if xk == 'fcont63':
                extra_pars['xaxis'] = r'$\log$ $\left[%s/F_\mathrm{6.3\ '%s1+\
                                      r'\mu m}\ (\mathrm{Hz})\right]$'
            elif xk == 'x':
                extra_pars['xaxis'] = r'$\log$ $\left[%s\ (\mathrm{W}/\mathrm{m}^2)\right]$'%s1
            else:
                iml = sample[0].molecule.molecule != xratsample.molecule.molecule
                s1 = sample[0].makeAxisLabel(iml)
                s2 = xratsample.makeAxisLabel(iml)
                extra_pars['xaxis'] = r'$\log$ $\left[%s/%s\right]$'%(s1,s2)
            #-- Change the xaxis name/min/max based on each plot
            #   ie if no errors are given, just take min and max and scale.
            if not xerr[xk]:
                extra_pars['xmin'] = min([min(yi) for yi in y[xk]])-0.2
                extra_pars['xmax'] = max([max(yi) for yi in y[xk]])+0.2
            #-- If errors are given, full plot, unless they are line-strength 
            #   ratios.
            elif xk not in ['x','xmdot','fcont63'] \
                 and xratsample.molecule.molecule != sample[0].molecule.molecule:
                extra_pars['xmin'] = -2
                extra_pars['xmax'] = 1
            else:
                if extra_pars.has_key('xmin'): del extra_pars['xmin']
                if extra_pars.has_key('xmax'): del extra_pars['xmax']
            
            pfn_xtag = '%s_eul_%i_wl_%.1f'\
                       %(sample[0].molecule.molecule,\
                         int(sample[0].getEnergyUpper()),\
                         float(sample[0].wavelength*10**4))
            if xk =='xmdot':
                pfn_xrat = '_mdot'
            elif xk == 'fcont63':
                pfn_xrat = '_fcont63'
            elif xk == 'x':
                pfn_xrat = ''
            else:
                ms = xratsample.molecule.isWater() and 'h2o' or 'co'
                ts = not xratsample.molecule.isWater() and xratsample.jup or xk
                pfn_xrat = '_%s%i'%(ms,ts)
            
        #-- Loop over the Y-AXIS KEYS
        for yk in y.keys():
            if yk not in ['y','ymdot','fcont63']:
                yratsample = sg[0]['GAS_LINES'][yk]
            if type(ypar) is types.StringType:
                extra_pars['yaxis'] = x_titles[ypar]
                extra_pars['ymin'] = min([min(yi) for yi in y[yk]])-0.2
                extra_pars['ymax'] = max([max(yi) for yi in y[yk]])+0.2
                pfn_ytag = pfn_parts[ypar]
                pfn_yrat = ''
            else: 
                s1 = sample[1].makeAxisLabel()
                if yk == 'ymdot':
                    extra_pars['yaxis'] = r'$\log$ $\left[%s/\dot{M}_\mathrm{'%s1+\
                                        r'g}\ (\mathrm{W}/\mathrm{m}^2\ \mathrm{yr}/\mathrm{M}_\odot)\right]$'
                if yk == 'fcont63':
                    extra_pars['yaxis'] = r'$\log$ $\left[%s/F_\mathrm{6.3\ '%s1+\
                                        r'\mu m}\ (\mathrm{Hz})\right]$'
                elif yk == 'y':
                    extra_pars['yaxis'] = r'$\log$ $\left[%s\ (\mathrm{W}/\mathrm{m}^2)\right]$'%s1
                else:
                    iml = sample[1].molecule.molecule != yratsample.molecule.molecule
                    s1 = sample[1].makeAxisLabel(iml)
                    s2 = yratsample.makeAxisLabel(iml)
                    extra_pars['yaxis'] = r'$\log$ $\left[%s/%s\right]$'%(s1,s2)
                
                #-- Change the yaxis name/min/max based on each plot
                #   ie if no errors are given, just take min and max and scale.
                if not yerr[yk]:
                    extra_pars['ymin'] = min([min(yi) for yi in y[yk]])-0.2
                    extra_pars['ymax'] = max([max(yi) for yi in y[yk]])+0.2
                #-- If errors are given, full plot, unless they are line-strength 
                #   ratios.
                elif yk not in ['y','ymdot','fcont63'] \
                    and yratsample.molecule.molecule != sample[1].molecule.molecule:
                    extra_pars['ymin'] = -2
                    extra_pars['ymax'] = 1
                else:
                    if extra_pars.has_key('ymin'): del extra_pars['ymin']
                    if extra_pars.has_key('ymax'): del extra_pars['ymax']
                
                pfn_ytag = '%s_eul_%i_wl_%.1f'\
                        %(sample[1].molecule.molecule,\
                            int(sample[1].getEnergyUpper()),\
                            float(sample[1].wavelength*10**4))
                if yk =='ymdot':
                    pfn_yrat = '_mdot'
                elif yk == 'fcont63':
                    pfn_yrat = '_fcont63'
                elif yk == 'y':
                    pfn_yrat = ''
                else:
                    ms = yratsample.molecule.isWater() and 'h2o' or 'co'
                    ts = not yratsample.molecule.isWater() and yratsample.jup or yk
                    pfn_yrat = '_%s%i'%(ms,ts)
            
            #-- Make an extra list of blends.
            xb, yb = [[]], [[]]            
            for blend1,blend2,xi,yi in zip(xblend[xk],yblend[yk],x[xk],y[yk]):
                blended = blend1 + blend2
                xb[-1].extend(xi[blended])
                yb[-1].extend(yi[blended])
            
            if xb[-1]: 
                extra_pars['keytags'] = keytags + ['$\mathrm{Blended}$']
                extra_pars['line_types'] = line_types + ['xk']
                extra_pars['markersize'] = markersize + [14]
                extra_pars['zorder'] = zorder + [max(zorder)+1]
            else: 
                xb, yb, = [], []
                extra_pars['keytags'] = keytags
                extra_pars['line_types'] = line_types
                extra_pars['markersize'] = markersize
                extra_pars['zorder'] = zorder 
            
            if add_linear_fit: 
                extra_pars['keytags'] = extra_pars['keytags'] + ['Mean Linear fit']
                extra_pars['line_types'] = extra_pars['line_types'] + ['-g'] + ['-k']*len(add_linear_fit['xgrid'])
                extra_pars['markersize'] = extra_pars['markersize'] + [4] + [4]*len(add_linear_fit['xgrid'])
                extra_pars['zorder'] = extra_pars['zorder'] + [min(zorder)-1] + [min(zorder)-2]*len(add_linear_fit['xgrid'])
                extra_pars['alpha'] = [1]*len(x[xk])+[1]*len(xb)+[1]+[0.002]*len(add_linear_fit['xgrid'])
                xb.append(add_linear_fit['xmean'])
                yb.append(add_linear_fit['ymean'])
                xb.extend(add_linear_fit['xgrid'])
                yb.extend(add_linear_fit['ygrid'])
                
            
            #-- Update from cfg file, in case any setting (except plotfile) has
            #   to be overridden. 
            extra_pars.update(cfg_dict)
            pfn = os.path.join(pfn_path,'%s_%s%s_vs_%s%s_%s'\
                                        %(mode,pfn_ytag,pfn_yrat,pfn_xtag,\
                                          pfn_xrat,pfn_ecl))
            extra_pars['filename'] = pfn
            ff = Plotting2.plotCols(x=xb and x[xk]+xb or x[xk],\
                                    y=yb and y[yk]+yb or y[yk],\
                                    yerr=yb and yerr[yk]+[None]*len(yb) or yerr[yk],\
                                    xerr=xb and xerr[xk]+[None]*len(xb) or xerr[xk],\
                                    **extra_pars)
            print ff
            
            if n_data > 0:                    
                print 'Stars plotted (in order of x):'
                for xi,yi,(ec,isg,isgd,isgm,nsgd,nsgm) \
                        in zip(x[xk],y[yk],ecl_num):
                    if ec[0] == 'm': continue
                    ifin = np.isfinite(xi) * np.isfinite(yi)
                    isort = np.argsort(xi[ifin])
                    sgsort = array(sg)[isgd][ifin][isort]
                    k = ', '.join(['%s = %s%s%s'
                                    %(keynames[con],\
                                      makeints[con] and str(int(v)) or str(v),\
                                      keyunits[con] and ' ' or '',\
                                      keyunits[con])
                                   for con,v in zip(extra_dcon,ec[1])])
                    print k, ': %s'%', '.join([s['STAR_NAME'] for s in sgsort])
            
            return ff
            

def guessRatio(line1,line1_err,line2,line2_err,line1_log=0,line2_log=0,\
               n_fit=10000,positive=0):
    
    '''
    Guess a ratio of given values with error bars a given number of times. 
    
    For both components of the ratio, values are drawn from a Gaussian 
    distribution around the given value, with the error as sigma. 
    
    Can be used for error analysis, or for estimating errors on ratios in case
    the errors on the separate components are difficult to propagate properly. 
    
    The option to guess a value within error bars in log space is possible. The 
    resulting value is then converted back to linear space, after which the 
    ratio is taken.
    
    Negative values can occur, due to the Gaussian nature of the guesses. Keep
    this in mind when taking the log of output values. If you do not want
    negative values, this van be requested via the positive keyword.
    
    A guess of the ratio, and a standard deviation, can be calculated by taking
    the mean and std of the columns in the ouput array.

    @param line1: Values of the first parameter on the y-axis
    @type line1: array
    @param line1_err: Uncertainties on line1, assuming they are in a normal 
                      distribution (1-sigma)
    @type line1_err: array
    @param line2: Values of the second parameter on the y-axis. Can be an empty
                  array in case you want to fit a correlation between par and 
                  line1 without any ratio involved. Pass an empty array if you
                  simply want to randomize a single array, instead of a ratio.
    @type line2: array
    @param line2_err: Uncertainties on line2, assuming they are in a normal 
                      distribution (1-sigma)
    @type line2_err: array
    @keyword line1_log: If line 1 is in log scale. In the ratio, 10**line1 is 
                        then taken.
                        
                        (default: 0)
    @type line1_log: bool
    @keyword line2_log: if line 2 is in log scale. In the ratio, 10**line1 is 
                        then taken.
                        
                        (default: 0)
    @type line2_log: bool
    @keyword n_fit: The number of times the correlation is fitted.
    
                    (default: 10000)
    @type n_fit: int
    @keyword positive: In some cases, you may want to disallow negative values,
                       eg when you take the log of the results. This switch 
                       allows you to exclude negative values from the output. 
                       Use this with caution! In some case, this will severely
                       affect the Gaussian distribution.
                       
                       (default: 0)
    @type positive: bool
    
    @return: The n_fit guesses of the requested ratio. 
    @rtype: array((n_fit,len(line1)))
    
    '''
   
    n_fit = int(n_fit)
    line1_log, line2_log = bool(line1_log), bool(line2_log)
    line1, line1_err = array(line1), array(line1_err)
    line2, line2_err = array(line2), array(line2_err)
    yarr = np.empty((n_fit,len(line1)))
    if positive:
        for n in range(n_fit):
            while True:
                guess1 = normal(line1, line1_err)
                if line1_log:
                    guess1 = 10**guess1
                    break
                elif False not in (guess1[np.isfinite(guess1)] > 0):
                    break
                    
            while True:
                if line2.size != 0:
                    guess2 = normal(line2, line2_err)
                else:
                    guess2 = 1
                    break
                if line2_log:
                    guess2 = 10**guess2
                    break
                elif False not in (guess2[np.isfinite(guess2)] > 0):
                    break
            yarr[n] = guess1/guess2
    else:
        for n in range(n_fit):
            guess1 = normal(line1, line1_err)
            if line1_log:
                guess1 = 10**guess1
            if line2.size != 0:
                guess2 = normal(line2, line2_err)
            else:
                guess2 = 1
            if line2_log and line2.size != 0:
                guess2 = 10**guess2
            yarr[n] = guess1/guess2
    return yarr

        
def fitCorrPolyLog(par1,par1_err,par2,par2_err,line1,line1_err,line2,line2_err,\
                   par1_log=0,par2_log=0,line1_log=0,line2_log=0,n_fit=10000,\
                   poly_degree=1,show=0,fn_plt='',x_for_yratio=0,\
                   y_for_xratio=0):

    '''
    Fit a polynomial to a data set.
    
    The data set can consist of straight up values or of ratios on both the x 
    and y axis (in log space). 
    
    Takes into account errors in both dimensions.
    
    Can be used for e.g. error estimation on a correlation. 
    
    @param par1: Values of the first parameter on x-axis.
    @type par1: array
    @param par1_err: Uncertainties on par, assuming they are in a normal 
                    distribution (1-sigma)
    @type par1_err: array
    @param par2: Values of the second parameter on the x-axis. Can be an empty
                 array in case you don't want a ratio on the x-axis.
    @type par2: array
    @param par2_err: Uncertainties on par2, assuming they are in a normal 
                     distribution (1-sigma)
    @type par2_err: array
    @param line1: Values of the first parameter on the y-axis
    @type line1: array
    @param line1_err: Uncertainties on line1, assuming they are in a normal 
                      distribution (1-sigma)
    @type line1_err: array
    @param line2: Values of the second parameter on the y-axis. Can be an empty
                  array in case you don't want a ratio on the x-axis.
    @type line2: array
    @param line2_err: Uncertainties on line2, assuming they are in a normal 
                      distribution (1-sigma)
    @type line2_err: array
    @keyword par1_log: If par is in log scale. If not, the log will be taken of 
                      par, since this method fits log log correlations.
                      
                      (default: 0)
    @type par1_log: bool
    @keyword par2_log: If par2 is in log scale. If not, the log will be taken  
                       of par2, since this method fits log log correlations.
                      
                       (default: 0)
    @type par2_log: bool
    @keyword line1_log: If line 1 is in log scale. In the ratio, 10**line1 is 
                        then taken.
                        
                        (default: 0)
    @type line1_log: bool
    @keyword line2_log: if line 2 is in log scale. In the ratio, 10**line1 is 
                        then taken.
                        
                        (default: 0)
    @type line2_log: bool
    @keyword n_fit: The number of times the correlation is fitted.
    
                    (default: 10000)
    @type n_fit: int
    @keyword poly_degree: The degree of the polynomial that is fitted. 
                          
                          (default: 1)
    @type poly_degree: int                          
    @keyword show: Show a plot with the results. If cfg is given, the plot is 
                   adapted, including the filename. 
                   
                   (default: 0)
    @type show: bool
    @keyword fn_plt: The filename of the plot, in case show is True, and a 
                     saved plot is requested. 
                     
                     (default: '')
    @type fn_plt: str
    @keyword x_for_yratio: Use the par grid as the second component in the y
                           ratio. This can be useful for instance if the ratio
                           has Mdot as numerator, while Mdot is also on the x
                           axis. In this case, you want to use the same random
                           value for the same point on both x and y.
                           
                           (default: 0)
    @type x_for_yratio: bool
    @keyword y_for_xratio: Use the line1 grid as the second component in the x
                           ratio. This can be useful for instance if the ratio
                           has Mdot as numerator, while Mdot is also on the y
                           axis. In this case, you want to use the same random
                           value for the same point on both x and y.
                           
                           (default: 0)
    @type y_for_xratio: bool
    
    @return: The fit results are returned for all n_fit fitted functions. The 
             parameters are the output of np.polyfit and the amount depends on 
             the polynomial degree.
    @rtype: array
    
    '''
    
    poly_degree = int(poly_degree)

    fitcoef = np.empty((n_fit, poly_degree+1))
    if y_for_xratio:
        xarr = guessRatio(par1,par1_err,[],[],line1_log=par1_log,n_fit=n_fit,\
                          positive=1)
        y1 = guessRatio(line1,line1_err,[],[],line1_log=line1_log,n_fit=n_fit,\
                        positive=1)
    else:
        xarr = guessRatio(par1,par1_err,par2,par2_err,line1_log=par1_log,\
                          line2_log=par2_log,n_fit=n_fit,positive=1)

    if x_for_yratio:
        yarr = guessRatio(line1,line1_err,[],[],line1_log,\
                          n_fit=n_fit,positive=1)
        x1 = guessRatio(par1,par1_err,[],[],line1_log=par1_log,\
                        n_fit=n_fit,positive=1)
    else:
        yarr = guessRatio(line1,line1_err,line2,line2_err,line1_log,line2_log,\
                          n_fit=n_fit,positive=1)
    
    for n,x,y in zip(range(n_fit),xarr,yarr):
        #-- Set up the dataset of x and y values.
        #   The x-values are drawn using gaussian distributed par values.
        #   The y-values are drawn using gaussian distributed line1/line2 ratio
        #   values. 
        #   For both x and y, checks are done for negative values, since the 
        #   log10 is taken of both of them. 
        xl = np.log10(x)
        yl = np.log10(y)
        if x_for_yratio:
            yl = yl - np.log10(x1[n])
        if y_for_xratio:
            xl = xl - np.log10(y1[n])
        fitcoef[n] = np.polyfit(xl, yl, poly_degree)
            
    if show and poly_degree == 1:
        #-- Plot a bunch of stuff        
        plt.figure(1)
        plt.clf()
        
        #-- Create a scatter plot of all fitted polynomials.
        plt.subplot(221)
                   
        x1 = par1
        if par2.size != 0:
            x2 = par2
        else:
            x2 = 1
        if par1_log:
            x1 = 10**x1
        if par2_log and par2.size != 0:
            x2 = 10**x2
        x = np.log10(x1/x2)
        
        y1 = line1
        if line2.size != 0:
            y2 = line2
        else:
            y2 = 1
        if line1_log: 
            y1 = 10**y1
        if line2_log and line2.size != 0:
            y2 = 10**y2
        y = np.log10(y1/y2)
        
        plt.scatter(x, y, color='blue', marker='o')
        
        x_grid = np.linspace(1.05*x.min(),0.95*x.max(),100)
        for n in range(0, fitcoef.shape[0], 2):
            y_grid = fitcoef[n,1] + fitcoef[n,0] * x_grid
            plt.plot(x_grid, y_grid, color="red", alpha = 0.006)
        
        
        plt.xlabel("log(X)")
        plt.ylabel("log(Y)")

        plt.subplot(222)
        plt.hexbin(fitcoef[:,0], fitcoef[:,1], bins=40)
        plt.xlabel("Slope")
        plt.ylabel("Intercept")

        plt.subplot(223)
        plt.hist(fitcoef[:,0], bins=40)
        plt.xlabel("Slope")
        plt.ylabel("N")

        plt.subplot(224)
        plt.hist(fitcoef[:,1], bins=40)
        plt.xlabel("Intercept")
        plt.ylabel("N")

        if not fn_plt:
            plt.show()
        else:
            plt.savefig(fn_plt)
            
    return fitcoef

    
def selectDataSG(sg,par,epar,par_co=None,edist=[]):
    
    '''
    Tool for selecting data from a star_grid. Mainly used by corrStarGrid() to
    define input arrays for the correlation study. 
    
    @param sg: The grid of Star() objects.
    @type sg: list[Star()]
    @param par: The requested parameter. If a string, a Star() keyword is used
                and if MDOT_GAS then the error bars are set as log10(3)/3. If 
                an integer, it is the index of a Transition() in the Star() 
                and the line strength and -error on it- of the line is taken. 
    @type par: str/int
    @param epar: The 1-sigma error bar of the parameter. Only relevant if par 
                 is a string Star() keyword that is not MDOT_GAS. 
    @type epar: array
    @keyword par_co: Define cutoff values here. Always given as an array of 
                     size 2. If a lower and/or upper boundary is not needed, 
                     it is set as None. In case of MDOT_GAS==xpar, these values 
                     are converted to log scale. Can be set as None if no 
                     cutoff is needed. Only relevant if par is a string.
                          
                     (default: None)
    @type par_co: array
    
    @keyword edist: Give the relative error of the distance here. Used to 
                     estimate an uncertainty on the rescaled line strengths 
                     according to distance (down to 100 pc). Not relevant when 
                     par is a string. An empty array implies no scaling. 
                    
                     (default: [])
    @type edist: array

    @return: The values for the requested parameter, as well as the 
             uncertainty, the cutoff values and the log request, ie all 
             keywords needed for the fitCorrPolyLog method.
    @rtype: (arr,arr,arr,bool)
    
    '''
    
    edist = array(edist)
    if type(par) is types.StringType:
        vals = array([s[par] for s in sg])
        if par_co <> None:
            if par_co[0] is None:
                par_co[0] = min(vals)
            if par_co[1] is None:
                par_co[1] = max(vals)    
            par_co = array(par_co,dtype=np.float)
        if par == 'MDOT_GAS':
            evals = np.ones(len(sg))*np.log10(3.)/3.
            vals = np.log10(vals)
            if par_co <> None: par_co = np.log10(par_co)
            vals_log = 1
        else:
            evals = epar
            vals_log = 0
    else:
        sample = sg[0]['GAS_LINES'][par]
        trans = Transition.getTransFromStarGrid(sg,sample,'sample')
        vals,evals = Transition.getLineStrengths(trans,mode='dint')
        #-- For now blends are ignored! Negative values cannot be included for
        #   the trend analysis and have to be handled more properly at a later
        #   stage if blends are somehow to be incorporated.
        vals = abs(vals)
        vals_log = 0
        if edist.size: 
            dists = array([s['DISTANCE'] for s in sg])
            vals = vals*dists**2/100.**2
            evals = vals*np.sqrt(4*edist**2+evals**2)
        else:
            evals = vals*evals    
        #-- Set a default for par_co since a cutoff is never needed for LS.
        finvals = vals[np.isfinite(vals)]
        par_co = array([min(finvals),max(finvals)])
    return (vals,evals,par_co,vals_log)
    
    
    
def corrSG(sg,xpar,ypar,expar=[],eypar=[],xratio=None,yratio=None,\
           eyratio=[],exratio=[],edist=[],xpar_co=(None,None),\
           ypar_co=(None,None),**kwargs):
    
    '''
    A method focused on finding correlations between parameters and/or data
    of multiple Star() objects. 
    
    @param sg: The grid of Star() objects.
    @type sg: list[Star()]
    @param xpar: The parameter on x-axis. If a string, a Star() keyword is used
                 and if MDOT_GAS then the error bars are set as log10(3)/3. If 
                 an integer, it is the index of a Transition() in the Star() 
                 and the line strength and -error on it- of the line is taken. 
                 If the same as yratio, the same guesses are used for both.
    @type xpar: str/int
    @param ypar: The parameter on y-axis. If a string, a Star() keyword is used
                 and if MDOT_GAS then the error bars are set as log10(3)/3. If 
                 an integer, it is the index of a Transition() in the Star()
                 and the line strength and -error on it- of the line is taken. 
                 If the same as xratio, the same guesses are used for both.
    @type ypar: int/str

    @keyword expar: The 1-sigma error bar of the parameter on the xaxis. Only 
                    relevant if xpar is a string Star() keyword that is not 
                    MDOT_GAS. 
    
                    (default: [])
    @type expar: array
    @keyword eypar: The 1-sigma error bar of the parameter on the yaxis. Only 
                    relevant if ypar is a string Star() keyword that is not 
                    MDOT_GAS. 
    
                    (default: [])
    @type eypar: array
    @keyword xratio: If a ratio on the x-axis is requested, this is the second 
                     component. Input syntax same as xpar.
                   
                     (default: None)
    @type xratio: int/str
    @keyword exratio: The relative error bars on the xratio parameter are given
                      here. Only relevant if xratio is a string Star() keyword
                      that is not MDOT_GAS.
                      
                      (default: [])
    @type exratio: array
    @keyword yratio: If a ratio on the y-axis is requested, this is the second 
                     component. Input syntax same as ypar.
                     
                     (default: None)
    @type yratio: int/str
    @keyword eyratio: The relative error bars on the yratio parameter are given
                      here. Only relevant if yratio is a string Star() keyword
                      that is not MDOT_GAS.
                      
                      (default: [])
    @type eyratio: array
    @keyword edist: Give the relative error of the distance here. Used to 
                    estimate an uncertainty on the rescaled line strengths 
                    according to distance (down to 100 pc). Not relevant when a
                    line ratio is requested.
                    
                    (default: [])
    @type edist: array                    
    @keyword xpar_co: Define cutoff values here. Always given as an array of 
                      size 2. If a lower and/or upper boundary is not needed, 
                      it set as None. In case of MDOT_GAS==xpar, these values 
                      are converted to log scale.
                          
                      (default: array(None,None))
    @type xpar_co: array
    @keyword ypar_co: Define cutoff values here. Always given as an array of 
                      size 2. If a lower and/or upper boundary is not needed, 
                      it set as None. In case of MDOT_GAS==ypar, these values 
                      are converted to log scale.
                          
                      (default: array(None,None))
    @type ypar_co: array
    @keyword kwargs: Extra keywords relevant for fitLinearCorr or 
                     fitCorrPolyLog
    @type kwargs: dict
    
    @return: The resulting fit parameters are returned. NYI if poly_degree!=1.
    @rtype: dict()
    
    '''
    
    #-- Get the line strengths of the correlated transition where not in every 
    #   case a scaling with distance is required. The absolute uncertainty is
    #   set here as well.
    expar, eypar, edist = array(expar), array(eypar), array(edist)
    exratio, eyratio = array(exratio), array(eyratio)
    xpar_co, ypar_co = array(xpar_co), array(ypar_co)
    ep = dict() 
    
    #-- The x-axis parameter is set here. 
    if xratio is None or type(xratio) is types.StringType:
        xv, exv, xpar_co, ep['par1_log'] = selectDataSG(sg,xpar,expar,xpar_co,\
                                                        edist)
    else:
        xv, exv, xpar_co, ep['par1_log'] = selectDataSG(sg,xpar,expar,xpar_co)
    
    #-- The x-ratio is set here. 
    if xratio is None or xratio == ypar:
        if xratio <> None: ep['y_for_xratio'] = 1
        xrat, exrat, ep['par2_log'] = array([]),array([]),0
    else:
        xrat, exrat, dummy, ep['par2_log'] = selectDataSG(sg,xratio,exratio)
    
    #-- The y-axis parameter is set here
    if yratio is None or type(yratio) is types.StringType:
        yv, eyv, ypar_co, ep['line1_log'] = selectDataSG(sg,ypar,eypar,\
                                                         ypar_co,edist)
    else:
        yv, eyv, ypar_co, ep['line1_log'] = selectDataSG(sg,ypar,eypar,ypar_co)
    
    #-- The y-ratio is set here.
    if yratio is None or yratio == xpar:
        if yratio <> None: ep['x_for_yratio'] = 1
        yrat, eyrat, ep['line2_log'] = array([]),array([]),0
    else:
        yrat, eyrat, dummy, ep['line2_log'] = selectDataSG(sg,yratio,eyratio)
    
    #-- Sort the grids according to xpar.
    isort = np.argsort(xv)
    xv, exv = xv[isort], exv[isort]
    yv, eyv = yv[isort], eyv[isort]
    if xrat.size:
        xrat, exrat = xrat[isort], exrat[isort]
    if yrat.size:
        yrat, eyrat = yrat[isort], eyrat[isort]
    
    #-- Make sure there are no NaNs in the grids.
    if xrat.size:
        xselfinite = np.isfinite(xv/xrat)
    else:
        xselfinite = np.isfinite(xv)
    if yrat.size:
        yselfinite = np.isfinite(yv/yrat)
    else:
        yselfinite = np.isfinite(yv)
    selfinite = yselfinite * xselfinite
    if xrat.size:
        xrat, exrat = xrat[selfinite], exrat[selfinite]
    if yrat.size:
        yrat, eyrat = yrat[selfinite], eyrat[selfinite]
    xv, exv = xv[selfinite], exv[selfinite]
    yv, eyv = yv[selfinite], eyv[selfinite]
    
    #-- Select subset based on cutoff vals for x and/or y (1st component only).
    xbools = (xv>=xpar_co[0]) * (xv<=xpar_co[1])
    ybools = (yv>=ypar_co[0]) * (yv<=ypar_co[1])
    bools = xbools*ybools
    xv, exv = xv[bools], exv[bools]
    yv, eyv = yv[bools], eyv[bools]
    if xrat.size:
        xrat, exrat = xrat[bools], exrat[bools]
    if yrat.size:
        yrat, eyrat = yrat[bools], eyrat[bools]
    kwargs.update(ep)

    allcoef = fitCorrPolyLog(par1=xv,par1_err=exv,par2=xrat,par2_err=exrat,\
                             line1=yv,line1_err=eyv,line2=yrat,\
                             line2_err=eyrat,**kwargs)
    results = dict()
    if kwargs.get('poly_degree',1) == 1:
        results['n_points'] = len(xv)
        results['slope'] = allcoef[:,0].mean()
        results['eslope'] = allcoef[:,0].std()
        results['intercept'] = allcoef[:,1].mean()
        results['eintercept'] = allcoef[:,1].std()
        
        for x, x_err, name in [(results['slope'],results['eslope'],"Slope"),\
                               (results['intercept'],results['eintercept'],\
                                "Intercept")]:
            print("{0} = {1} +/- {2}".format(name, x, x_err))
        
        corrcoef = np.corrcoef(allcoef[:,0], allcoef[:,1])
        results['corrcoef'] = corrcoef[1,0]
        results['covariance'] = results['corrcoef']*results['eslope']\
                                    *results['eintercept']
        results['results'] = allcoef
        print("The correlation coefficient between slope & intercept is {0}."\
              .format(results['corrcoef']))
        print("This leads to a covariance of {0} for slope & intercept."\
              .format(results['covariance']))
        print("Finally, {0} data points were available to produce this fit."\
              .format(results['n_points']))
    else:
        results['poly_degree'] = kwargs.get('poly_degree')
        results['results'] = allcoef
    return results