# -*- coding: utf-8 -*-

"""
Tools for inter- and extrapolation.

Author: R. Lombaert

"""

from scipy import array, hstack
from scipy import exp
from scipy.optimize import leastsq
from scipy import isnan

from cc.plotting import Plotting2


def getResiduals(p,x,y,func='power'):
    
    """
    Calculate residuals between function evaluation and given y-values.
    
    @param p: parameters for function evaluation
    @type p: list
    @param x: the x-grid evaluated by function
    @type x: array
    @param y: The given y-values
    @type y: array
    
    @keyword func: type of function f(x), can be ['power','exp','power10',\
                   'square','linear'] for now.
                      
                   (default: 'power')
    @type func: string
    @return: residuals as a function of the x-grid
    @rtype: array
    
    """    
    
    x,y,p = array(x),array(y),array(p)
    err = y - pEval(x,p,func)                        
    return err
    


def pEval(x,p,func='power',x_unit=None):
    
    """
    Evaluate a function prescription.
    
    @param x: The x-grid for function evaluation
    @type x: array
    @param p: parameters for function evaluation
    @type p: list

    @keyword func: type of function f(x), can be ['power','exp','power10',\
                   'square','linear'] for now.
                       
                   (default: 'power')
    @type func: string
    @keyword x_unit: the unit of the x array, only relevant if func=='planck'. 
                     Units other than micron not yet available.
    
                     (default: None)
    @type x_unit: string
    
    @return: f(x)
    @rtype: array
    
    """
    
    x,p = array(x),array(p)
    if func.lower() == 'power':
        return p[0] + p[2]*(x+p[1])**(p[3])
    if func.lower() == 'exp':
        return p[0] + p[2]*(exp(p[3]*(x+p[1])))
    elif func.lower() == 'power10':
        return p[0] + p[2]*(10**(p[3]*(x+p[1])))
    elif func.lower() == 'square':
        return p[0] + p[2]*(x+p[1]) + p[3]*(x+p[1])**2
    elif func.lower() == 'linear':
        return p[0] + p[2]*(x + p[1])
    elif func.lower() == 'planck':
        pass
        if x_unit != 'micron':
            raise IOError('Keyword x_unit must be micron. Other options not yet available.')
        c = 2.99792458e10          #in cm/s
        h = 6.62606957e-27         #in erg*s Planck constant
        k = 1.3806488e-16          #in erg/K Boltzmann constant 
        f = c/x
        #planck = 2.*h/c**2.*f**3.*1./(exp(h*f/k/T)-1.)
    
        
        
def fitFunction(x_in,y_in,x_out,show=0,func='linear',initial=[1,0,0.5,-1.5],\
                maxfev=20000):
     
    '''
    Fit a template function to input data.
    
    @param x_in: The input x values
    @type x_in: list/array
    @param y_in: The input y values
    @type y_in: list/array
    @param x_out: The output x grid
    @type x_out: list/array

    @keyword show: show the fit result
    
                   (default: 0)
    @type show: bool
    @keyword func: The function template (see pEval)
    
                   (default: 'linear')
    @type func: string
    @keyword initial: initialisation parameters for your function, number of
                      parameters depends on type of function (see pEval)
    
                      (default: [1,0,0.5,-1.5])
    @type initial: list(float)
    @keyword maxfev: Maximum number of calls to the function for the leastsq 
                     method. Zero means (1+number of elements in x_in)*100
                     
                     (default: 20000)
    @return: The output y values
    @rtype: list/array
    
    '''
    
    p_lsq = leastsq(getResiduals,initial,args=(x_in,y_in,func),\
                    maxfev=20000)[0]
    y_out = pEval(x_out,p_lsq,func=func)
    if show:
        print p_lsq
        Plotting2.plotCols(x=[x_in,hstack([x_out,x_in])],\
                           y=[y_in,hstack([y_out,pEval(x_in,p_lsq,func=func)])],\
                           xlogscale=1,ylogscale=1)
    return y_out
          
          
              
def linInterpol(x_in,y_in,x):
    
    """
    Linear interpolation between two points.
    
    @param x_in: 2 x-points for interpolation
    @type x_in: list
    @param y_in: 2 y-points for interpolation corresponding to x_in
    @type y_in: list
    @param x: x value for interpolation
    @type x: float/array
    
    @return: f(x), where f is linear interpolation from x_in and y_in
    @rtype: float/array
    
    """
    
    x = array(x)
    try:
        return y_in[0] + (y_in[0]-y_in[1])/(x_in[0]-x_in[1]) * (x-x_in[0])
    except ZeroDivisionError:
        print 'Identical x-coordinates were submitted: Division by zero. ' + \
              'Aborting.'
        return





#def doInterpol(x_in,y_in,gridsx,xmin=0,xmax=0,func='power',p0=[1,0,0.5,-1.5],stitch=False,x_cutoff=10.0,y_cutoff=60.0):
    #"""Inter/extrapolation of a dataset.
    
    #Input:    x_in=list/array, input x coordinates to be interpolated
        #y_in=list/array, input y coordinates to be interpolated
        #gridsx=list, includes multiple lists with x coordinates for which new y coordinates must be calculated
            #from the inter/extrapolation of (x_in,y_in), if list=[0] it is ignored and a corresponding
            #y-list [0] is added to the final list of lists giving the results
        #OPTIONAL xmin=int/float/..., default=0, lower boundary xrange for least squares fit in case of extrapolation, 
            #if x_in[0] > xmin then xmin is x_in[0] by default, the default option is 0 in case no
            #extrapolation will be needed and the parameter can be ignored
        #OPTIONAL xmax=int/float/..., default=0, upper boundary xrange for least squares fit in case of extrapolation, 
            #if x_in[len(x_in)-1] < xmax then xmax is x_in[len(x_in)-1] by default, the default option is 0 in case no
            #extrapolation will be needed and the parameter can be ignored
        #OPTIONAL func=string, default='power', function used for the Levenberg-Marquardt
            #least squares minimisation procedure (scipy.optimization.leastsq) in case of extrapolation
        #OPTIONAL p0=list/array, default=[1,0,0.5,-1.5], parameters for function definition in case of extrapolation:
            #[a,b,c,d] in (for instance) power law a+c*(x+b)**d for least squares fit
        #OPTIONAL stitch=boolean, default=0, if func=='linear' then stitch==1 will shift the resulting y-values, such that
            #they connect smoothly with given data in gridsx at the connection point
        #OPTIONAL x_cutoff=float, default=10.0, in stitch mode, this x value is where the extrapolated data are connected
            #to the real input data
        #OPTIONAL y_cutoff=float, default=10.0, in stitch mode, this y value is where the extrapolated data are connected
            #to the real input data at x_cutoff and is equal to the real data value at this x value
    #Output:    list of lists with resulting y-coordinates corresponding to x-coordinates
            #after the interpolation.
        #list of best least squares fit parameters for the input (x_in,y_in) coordinates in indicated xrange
    #"""
    ##calculate least squares best fit in case extrapolation is needed     
    
    #x_in,y_in,p0 = array(x_in),array(y_in),array(p0)
    #if not (xmin == 0 and xmax == 0):
        #try:
            #p_lsq = leastsq(getResiduals,p0,args=([x for x in x_in if x >= xmin and x <= xmax],[y for x,y in zip(x_in,y_in) if x >= xmin and x <= xmax],func),maxfev=20000)[0]
        #except TypeError:
            #print 'Only 2 or less points have been submitted for least squares fit: This is not possible.'
            #print 'Please submit a different inter/extrapolation dataset or change boundaries for fitting. Aborting now.'
            #return
    #elif (xmin == 0 and xmax == 0):
        #try:
            #p_lsq = leastsq(getResiduals,p0,args=(x_in,y_in,func),maxfev=20000)[0]
        #except TypeError:
            #print 'Only 2 or less points have been submitted for least squares fit: This is not possible.'
            #print 'Please submit a different inter/extrapolation dataset or change boundaries for fitting. Aborting now.'
            #return
    #else:
        #p_lsq = [0,0,0,0]
    #gridsy = []
    #for index,x in enumerate(gridsx):
        #y = []    
        #if [xi for xi in x if isnan(xi)]:
             #print "Warning! There are NaNs in the x-grid. The following may fail!"
        #for i,xi in enumerate(x):
            #if xi == 0 and len(x) == 1:
                #y = [0]
            #elif (xi < x_in[0]) or (xi > x_in[-1]):        #extrapolation ONLY if linear interpolation is not possible    
                #if stitch:
                    #y.append(float(pEval(xi,[y_cutoff,-x_cutoff,p_lsq[2],p_lsq[3]],func)))
                #else:
                    #y.append(float(pEval(xi,p_lsq,func)))
            #elif xi == x_in[len(x_in)-1]:      #equality is checked for last element, as that case is not covered by next step
                #y.append(y_in[len(x_in)-1])
            #else:
                ##Find boundaries in between which xi can be found, and remember index!
                ##lower boundary may be equal to xi, upper boundary may not! 
                ##Hence the boundaries can never be equal, such as in the case that one or more items in x_in are duplicate
                #boundaries = [[x_in[j],j] for j in xrange(len(x_in)) if (x_in[j] <= xi and x_in[j+1] > xi) or (x_in[j-1] <= xi and x_in[j] > xi)] 
                #y.append(y_in[boundaries[0][1]] + (xi-x_in[boundaries[0][1]]) * ((y_in[boundaries[1][1]]-y_in[boundaries[0][1]])/float(x_in[boundaries[1][1]]-x_in[boundaries[0][1]])))
        #gridsy.append(y)
    #return gridsy,p_lsq
    
    