import numpy as np
from math import sin,cos,exp
import matplotlib.pyplot as plt
import pylab as pl

from cc.plotting import Plotting2



def expDensFunc(r,rSw,exponent):

    return exp(-(r/rSw)**exponent)



@np.vectorize
def densityFunc(r,theta,A,B,C,D,E,F,rMin,rSw):

    exponent = -B*(1+C*(sin(theta)**F)*(expDensFunc(r,rSw,D)/\
                                        expDensFunc(rMin,rSw,D)))
    rho = (r/rMin)**exponent
    rho = rho * (1+A*((1-cos(theta))**F)*(expDensFunc(r,rSw,E)/\
                                          expDensFunc(rMin,rSw,E)))

    return rho



def plotDens(A,B,C,D,E,F,rPlot,rMin,rSw,nRad,nTheta,filename=None,\
             extension='pdf',landscape=0,show=0,figsize=(5.5, 10)):
    
    #-- Define the radial and angular grids
    R = np.logspace(np.log10(rMin), np.log10(rPlot), num=nRad, endpoint=True)
    Theta = np.arange(0, np.pi/2.0+1.0/float(nTheta), (np.pi/2.0)/(nTheta))
    R,Theta = np.meshgrid(R,Theta)
    
    #-- Calculate the Meixner density distribution
    rho = densityFunc(R,Theta,A,B,C,D,E,F,rMin,rSw)

    #-- Figure consists of two subplots, one with the density map, and one with
    #   profiles at 0 and 90 degrees.
    fig = plt.figure(figsize=figsize)
    
    #-- Make the color scale density plot
    ax = fig.add_subplot(211)
    pl.pcolor(R*np.sin(Theta),R*np.cos(Theta),np.log10(rho/rho.max()))
    pl.pcolor(-R*np.sin(Theta),R*np.cos(Theta),np.log10(rho/rho.max()))
    pl.pcolor(-R*np.sin(Theta),-R*np.cos(Theta),np.log10(rho/rho.max()))
    pl.pcolor(R*np.sin(Theta),-R*np.cos(Theta),np.log10(rho/rho.max()))
    fig.subplots_adjust(right=0.8)
    fig.subplots_adjust(bottom=0.1)
    cbar_ax = fig.add_axes([0.85, 0.54, 0.05, 0.36])
    pl.colorbar(cax=cbar_ax)
    
    #-- Add density distribution at 0 and 90 degrees, and a r^-2 distribution
    ax2 = fig.add_subplot(212)
    ax2.plot(np.log10(R[0]),np.log10(rho[0]/rho[0].max()),'black',\
             lw=2,marker='x',label=r'$\theta = 0^\circ$')
    ax2.plot(np.log10(R[-1]),np.log10(rho[-1]/rho[-1].max()),'magenta',\
             lw=2,marker='|', label=r'$\theta = 90^\circ$')
    ax2.plot(np.log10(R[0]),np.log10(rMin*rMin/(R[0]*R[0])),'green',lw=2,\
             lw=2,label=r'$r^{-2}$')
    ax2.legend()
    label = r'$\rho_0[r_{\rm min}] / \rho_{90}[r_{\rm min}] = $'+\
             '{:10.3e}'.format(rho[0][0]/rho[-1][0])
    ax2.text(0.25,1.02,label,transform=ax2.transAxes)

    
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

    plotDens(A,B,C,D,E,F,rPlot,rMin,rSw,nRad,nTheta,show=1)
