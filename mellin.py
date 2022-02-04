import sys,os
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import gamma
from mpmath import fp
import params as par
from numba import jit

npts=500
x,w=np.polynomial.legendre.leggauss(npts)

a,b=0,1
u=0.5*(b-a)*x+0.5*(a+b)
ju=1/(1-u)**2    

Z=u/(1-u)
W=w
JAC=0.5*(b-a)*ju

#--gen mellin contour
c=1.9 
phi=3/4*np.pi

N=c+Z*np.exp(complex(0,phi)) 
phase= np.exp(complex(0,phi))

@jit(nopython=True)
def invert(x,F):
    return np.sum(np.imag(phase * x**(-N) * F)/np.pi * W * JAC)
