# -*- coding: utf-8 -*-

"""
Methods for water abundance determination.

Author: R. Lombaert

"""

import os
from scipy.integrate import trapz
from scipy import average
import math

from cc.modeling.codes import Gastronoom



def getWaterInfo(star):
    
    '''
    Return abundance information for water for a Star().
    
    @param star: The Star() object for which to collect all water info.
    @type star: Star()
    
    @return: Column densities for: h2o vapor in ice shell, h2o ice, h2 in ice 
             shell, h2o vapor in full shell and h2 in full shell
    @rtype: (float,float,float,float,float)
    
    '''
    
    rinds = [star[rs]*star['R_STAR']*star.r_solar
            for rs in ['R_DES_H2O','R_DES_AH2O','R_DES_CH2O']
            if star.has_key(rs)]
    ah2os = [star[ah]
            for ah in ['A_H2O','A_AH2O','A_CH2O']
            if star.has_key(ah) and star[ah] != 0]
    print 'The average condensation radius for water ice is %.2f.'\
          %(average(rinds)/star['R_STAR']/star.r_solar)
    rin = star['R_INNER_GAS']*star['R_STAR']*star.r_solar
    rout = star['R_OUTER_DUST']*star['R_STAR']*star.r_solar
    rstar = star['R_STAR']*star.r_solar
    g_id = star['LAST_GASTRONOOM_MODEL']
    mdot_dust = star['MDOT_DUST']
    mdot_gas = star['MDOT_GAS']
    vexp_dust = star['V_EXP_DUST']
    vexp_gas = star['VEL_INFINITY_GAS']
    nh2 = getH2ColDens(rin=average(rinds),rout=rout,mdot_gas=mdot_gas,\
                       vexp_gas=vexp_gas)
    nh2o = getWaterColDens(g_id,star.path_gastronoom,average(rinds),rout,\
                           star['OPR'])
    nh2o_full = getWaterColDens(g_id,star.path_gastronoom,rin,rout,star['OPR'])
    nh2_full = getH2ColDens(rin,rout=rout,mdot_gas=mdot_gas,vexp_gas=vexp_gas)
    nh2o_ice = 0
    for ah,rind in zip(ah2os,rinds):
        nh2o_ice += getWaterIceColDens(rind,rout,mdot_dust,ah,vexp_dust)
    return (nh2o,nh2o_ice,nh2,nh2o_full,nh2_full)



def getWaterColDens(model_id,path_gastronoom,rin,rout,opr):
    
    '''
    Get the water number density from GASTRoNOoM output.
    
    @param model_id: the model_id of the GASTRoNOoM model
    @type model_id: string
    @param path_gastronoom: the subfolder in GASTRoNOoM for model output
    @type path_gastronoom: string
    @param rin: the inner radius from which to integrate over the whole shell 
                (cm)
    @type rin: float
    @param rout: the outer radius from which to integrate over the whole shell 
                 (cm)
    @type rout: float
    @param opr: The ortho-to-para water ratio
    @type opr: float
    
    @return: the water vapout column density in the shell between rin and rout 
             (cm-2)
    @rtype: float
    
    '''
    
    r_solar = 6.955e10
    opr = float(opr)
    rin = float(rin)
    rout = float(rout)
    rmed = (rout+rin)/2.
    rdelta = (rout-rin)/2.
    filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
               path_gastronoom,'models',model_id,'coolfgr_all%s.dat'%model_id)
    rad = Gastronoomn.getGastronoomOutput(filename=filename,keyword='RADIUS',\
                                   return_array=1)
    nh2o = Gastronoom.getGastronoomOutput(filename=filename,keyword='n(h2o)',\
                                    return_array=1)
    nh2o_sel = nh2o[abs(rad - rmed ) <= rdelta]
    rad_sel = rad[abs(rad - rmed ) <= rdelta]
    nh2o_tot = trapz(x=rad_sel,y=nh2o_sel)
    nh2o_op = nh2o_tot + nh2o_tot/opr
    return nh2o_op
    
    
    
def getWaterIceColDens(rin,rout,mdot_dust,ah2o,vexp_dust):
    
    '''
    Calculate the water molecular content of an ice shell.
    
    @param rin: the inner radius from which to integrate over the whole shell 
                (cm)
    @type rin: float
    @param rout: the outer radius from which to integrate over the whole shell 
                 (cm)
    @type rout: float
    @param mdot_dust: The total dust mass-loss  (Msolar/yr)
    @type mdot_dust: float
    @param ah2o: The dust mass fraction in water ice
    @type ah2o: float
    @param vexp_dust: The expansion dust velocity
    @type vexp_dust: float
    
    @return: the molecular ice column density in the shell between rin and rout
             (cm-2)
    @rtype: float
    
    '''
    
    rin = float(rin)
    rout = float(rout)
    mdot_dust = float(mdot_dust)
    ah2o = float(ah2o)
    vexp_dust = float(vexp_dust)
    r_solar = 6.955e10       #in cm
    m_solar = 1.98892e33     #in g
    year = 31557600. 
    mh2o = 18.
    avogadro = 6.022e23
    sigma = (1./rin-1./rout)*ah2o*mdot_dust \
            /vexp_dust/4./math.pi*m_solar/100000./year
    nh2o_ice = sigma * avogadro / mh2o
    return nh2o_ice
    
def getH2ColDens(rin,rout,mdot_gas,vexp_gas):
    
    '''
    Calculate the H2 content of a shell.
    
    @param rin: the inner radius from which to integrate over the whole shell 
                (cm)
    @type rin: float
    @param rout: the outer radius from which to integrate over the whole shell 
                 (cm)
    @type rout: float
    @param mdot_gas: The total gas mass-loss  (Msolar/yr)
    @type mdot_gas: float
    @param vexp_gas: The expansion dust velocity
    @type vexp_gas: float
    
    @return: the H2 column density in the shell between rin and rout (cm-2)
    @rtype: float
    
    '''
    
    rin = float(rin)
    rout = float(rout)
    mdot_gas = float(mdot_gas)
    vexp_gas = float(vexp_gas)
    r_solar = 6.955e10       #in cm
    m_solar = 1.98892e33     #in g
    year = 31557600. 
    mh2 = 2.
    avogadro = 6.022e23
    sigma = (1./rin-1./rout)*mdot_gas \
            /vexp_gas/4./math.pi*m_solar/100000./year
    nh2 = sigma * avogadro / mh2
    return nh2