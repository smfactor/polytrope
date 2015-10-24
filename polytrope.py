import numpy as np
import scipy.integrate as intg
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

def LEeqn(th,z,n):
    return np.array([th[1],(-2./z)*th[1]-(th[0]**n)])

def thPN(th,z):
    return np.array([(-1./3.)*z+(1./30.)*z**3, (-1./3.)+(1./10.)*z**2])

def poly(n,zi,dz,filename):
    zs = np.arange(0.,zi,dz)
    nz = np.size(zs)

    y = intg.odeint(thPN,np.array([1.,0.]),zs)

    while y[-1,0]>0.:
        nz = np.size(zs)
        zs2 = np.arange(zs[-1],zs[-1]+100.*dz,dz)
        zs = np.append(zs,zs2)
        y = np.append(y,intg.odeint(LEeqn,y[nz-1],zs2,args=tuple([n])),axis=0)
        sys.stdout.write("\rz={0},th={1}".format(zs[-1],y[-1,0]))
        sys.stdout.flush()
    sys.stdout.write("\n")

    nover = np.size(np.where(y[:,0]<0.)-1
    print nover
    zs = zs[:-nover]
    y = y[:-nover]

    data = np.array([zs,y[:,0],-(zs**2)*y[:,1],(-3./zs)*y[:,1]])

    fits.writeto(filename+'.fits',data)
    return data
