import numpy as np
from math import sin,cos,exp
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pylab as pl

from cc.plotting import Plotting2



def expDensFunc(r,rSw,exponent):

    return exp(-(r/rSw)**exponent)



@np.vectorize
def densityFunc(r,theta,A,B,C,D,E,F,rmin,rSw):

    exponent = -B*(1+C*(sin(theta)**F)*(expDensFunc(r,rSw,D)/\
                                        expDensFunc(rmin,rSw,D)))
    rho = (r/rmin)**exponent
    rho = rho * (1+A*((1-cos(theta))**F)*(expDensFunc(r,rSw,E)/\
                                          expDensFunc(rmin,rSw,E)))

    return rho



def plotDens(A,B,C,D,E,F,rPlot,rMin,rSw,nRad,nTheta,filename=None,\
             extension=None,landscape=0,show=1):
    

    R = np.logspace(np.log10(rMin), np.log10(rPlot), num=nRad, endpoint=True)
    thetaArray = np.arange(0, np.pi/2.0+1.0/float(nTheta), (np.pi/2.0)/(nTheta))
    R,thetaArray = np.meshgrid(R,thetaArray)
    rho = densityFunc(R,thetaArray,A,B,C,D,E,F,rMin,rSw)

    rho = rho/rho.max()

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    pl.pcolor(R*np.sin(thetaArray),R*np.cos(thetaArray),np.log10(rho))
    pl.pcolor(-R*np.sin(thetaArray),R*np.cos(thetaArray),np.log10(rho))
    pl.pcolor(-R*np.sin(thetaArray),-R*np.cos(thetaArray),np.log10(rho))
    pl.pcolor(R*np.sin(thetaArray),-R*np.cos(thetaArray),np.log10(rho))
    pl.colorbar()
    
    if filename: filename = Plotting2.saveFig(filename,extension,landscape)
    if show: pl.show()
    
    return filename



if __name__ == '__main__':

    A = 100.0
    B = 2.0
    C = 1.0
    D = 0.5
    E = 0.0
    F = 2.0
    rPlot = 1E16
    rMin = 1E13
    rSw = 5E13
    nTheta = 50
    nRad = 50

    plotDens(A,B,C,D,E,F,rPlot,rMin,rSw,nRad,nTheta,filename='test')
