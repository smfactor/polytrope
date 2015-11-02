import numpy as np
import scipy.integrate as intg
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

def LEeqn(n):
    def LEeq(th,z):
        return np.array([th[1],(-2./z)*th[1]-(np.abs(th[0])**n)])
    return LEeq

def thPN(th,z):
    return np.array([(-1./3.)*z+(1./30.)*z**3, (-1./3.)+(1./10.)*z**2])

def poly(n,zi,dz,filename):
    LEeq = LEeqn(n)
    zs = np.arange(0.,zi,dz)
    nz = np.size(zs)

    y = intg.odeint(thPN,np.array([1.,0.]),zs)

    while y[-1,0]>0.:
        zs2 = np.arange(zs[-1],zs[-1]+200.*dz,dz)
        zs = np.append(zs,zs2)
        y = np.append(y,intg.odeint(LEeq,y[-1],zs2),axis=0)
        sys.stdout.write("\rz={0},th={1}".format(zs[-1],y[-1,0]))
        sys.stdout.flush()
    sys.stdout.write("\n")

    ngood= np.size(np.where(y[:,0]>0.))
    zs = zs[:ngood]
    y = y[:ngood]

    data = np.array([zs,y[:,0],-(zs**2)*y[:,1],(-3./zs)*y[:,1]])

    fits.writeto(filename+'.fits',data)
    return data
