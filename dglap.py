import sys,os
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import gamma
from mpmath import fp
import params as par
from numba import jit

import params as par
import alfa
import mellin
N     = mellin.N

Nf    = 1
beta0 = 4*Nf/3
beta1 = 4*Nf

def get_norm(N):
    return np.sqrt(N.real**2+N.imag**2)

def _psi1(z):
    if get_norm(z)<1e3:
        return fp.psi(1,z)
    else:
        return 1/z+1/2/z**2+1/6/z**3-1/30/z**5+1/42/z**7-1/30/z**9+5/66/z**11-691/2730/z**13+7/6/z**15

psi0 = lambda N: np.array([fp.psi(0,complex(n.real,n.imag)) for n in N],dtype=complex)
psi1 = lambda N: np.array([_psi1(complex(n.real,n.imag)) for n in N],dtype=complex)

S1   = fp.euler + psi0(N+1)

aini = alfa.get_alfa(par.me2)
f1   = - psi0(N)**2 - 2*fp.euler*psi0(N) \
       - psi0(N + 2)**2 - 2*fp.euler*psi0(N + 2)\
       + psi1(N) + psi1(N + 2) - np.pi**2/3 - 2*fp.euler**2 + 7/2

d1   = f1 -2*psi1(N) - 2*psi1(N + 2) - 5/2 + 2*np.pi**2/3

f1*= aini/(2*np.pi)     
d1*= aini/(2*np.pi)     

P  = (3.0+2.0/N/(N+1)-4.0*S1) #--PQQ
R  = P/beta0
aini=aini

@jit(nopython=True)
def evolve(mu2,order):
    afin = alfa.get_alfa(mu2)
    L    = np.log(afin/aini)
    U=np.exp(R*L)
    if order==0:
        F=U 
        D=U 
    elif order==1:
        F=U*(1 + f1) 
        D=U*(1 + d1)
    return F,D

@jit(nopython=True)
def get_LDF(x,F,mod=False): 
    #return (1-x)**(-0.9)    
    if mod==False:
        return mellin.invert(x,F)
    else:
        return mellin.invert(x,F/(N-1))

@jit(nopython=True)
def get_LFF(z,D,mod=False): 
    #return (1-z)**(-0.9) 
    if mod==False:
        return mellin.invert(z,D)
    else:
        return mellin.invert(z,D/(N-1))
