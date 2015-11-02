import numpy as np
import scipy.integrate as intg
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

def LEeqn(n):
    '''
    Lane-Emden equation.
    Given an n this function returns a new function which is the L-E equation for that n. That function rReturns an array of Theta and Theta' as a function of Theta and z
    '''
    def LEeq(th,z):
        #use abs so dont get a complex number when th<0 and n is not an integer
        #dont care about any values when th is negative anyway
        return np.array([th[1],(-2./z)*th[1]-(np.abs(th[0])**n)])
    return LEeq

def thPN(th,z):
    #3rd order taylor expansion of the L-E equation
    return np.array([(-1./3.)*z+(1./30.)*z**3, (-1./3.)+(1./10.)*z**2])

def poly(n,zi,dz,filename):
    '''
    Integrates the L-E equation and returns and writes to a file the important values.
    n is the index, zi is the interval to use the polynomial solution, dz is the step size, filename is the prefix for the output file (filename+.fits)
    '''
    #set up equation and variables
    LEeq = LEeqn(n)
    zs = np.arange(0.,zi,dz)
    nz = np.size(zs)

    #integrate over initial interval using polynomial
    y = intg.odeint(thPN,np.array([1.,0.]),zs)

    #integrate until th<0
    while y[-1,0]>0.:
        zs2 = np.arange(zs[-1],zs[-1]+200.*dz,dz)
        zs = np.append(zs,zs2)
        y = np.append(y,intg.odeint(LEeq,y[-1],zs2),axis=0)
        #write progress
        sys.stdout.write("\rz={0},th={1}".format(zs[-1],y[-1,0]))
        sys.stdout.flush()
    sys.stdout.write("\n")

    #only keep values where th>0
    ngood= np.size(np.where(y[:,0]>0.))
    zs = zs[:ngood]
    y = y[:ngood]

    #prepare output array
    data = np.array([zs,y[:,0],-(zs**2)*y[:,1],(-3./zs)*y[:,1]])

    #write and return data
    fits.writeto(filename+'.fits',data)
    return data
