# -*- coding: utf-8 -*-

"""
Performing basic statistics.

Author: R. Lombaert

"""

import types
import numpy as np


def calcChiSquared(data,model,noise,ndf=0,mode='diff'):
    
    """
    Calculate the reduced chi-squared value of a data array minus a model array,
    taking into account the noise in the data array.
    
    @param data: The data set. Must have same dimensions as model!
    @type data: array
    @param model: The model array. Must have same dimensions as data!
    @type model: array
    @param noise: the noise in the data array. Give one value for overall noise
                  or individual values for every entry in data/model. 
    @type noise: float/array

    @keyword ndf: Number of degrees of freedom. Default in case of calculating
                  for one single model. Typically the number of variable grid
                  parameters in a grid calculation.
                  
                  (default: 0) 
    @type ndf: int
    @keyword mode: The method used for the chi^2 calculation. 'diff' is the 
                   standard differentiation of the chi^2. 'log' redistributes 
                   the ratio of data and model points on a logarithmic scale 
                   such that lower than 1 or larger than 1 are essentially 
                   equivalent. This removes bias in either direction of 1. Other
                   than the input array distribution the chi^2 'log' method is 
                   mathematically equivalent to the differentiation.
                   
                   (default: 'diff')
    @type mode: str
    
    @return: The chi squared value
    @rtype: float
    
    """
    
    mode = str(mode).lower()
    if type(data) not in [types.ListType,scipy.ndarray]:
        data = [data]
    data, model, noise = np.array(data), np.array(model), np.array(noise) 
    if mode == 'diff':
        chi2 = ((data - model)**2./noise**2.).sum()/(len(data)-ndf-1)
    elif mode == 'log':
        chi2 = ((10**abs(np.log10(data/model))-1)**2./(noise/model)**2.).sum()
        chi2 /= (len(data)-ndf-1)
    else:
        print 'Chi^2 mode not recognized.'
        chi2 = None
    return chi2
    
    

def calcLoglikelihood(data,model,noise):
    
    """
    Calculate the loglikelihood value of a data array minus a model array,  
    taking into account the noise in the data array.
    
    @param data: The data set. Must have same dimensions as model!
    @type data: array
    @param model: The model array. Must have same dimensions as data!
    @type model: array
    @param noise: the noise in the data array. 
    @type noise: float/array

    @return: The loglikelihood value
    @rtype: float
    
    """
    
    data, model, noise = np.array(data), np.array(model), np.array(noise) 
    lll = (-np.log(np.sqrt(2.*np.pi)) - np.log(noise) - 1./2.*((data-model)/noise)**2.).sum()
    return lll
    
    

                