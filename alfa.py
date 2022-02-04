import sys,os
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.special import gamma
from mpmath import fp
import params as par
from numba import jit

@jit(nopython=True)
def beta_func(a,Nf,order):
    beta0 = 4*Nf/3
    beta1 = 4*Nf
    betaf = beta0
    if order>=1: betaf+=a*beta1
    return betaf*a**2

@jit(nopython=True)
def evolve_a(Q20,a,Q2,Nf,order):
    # Runge-Kutta implemented in pegasus  
    LR = np.log(Q2/Q20)/20.0
    for k in range(20):
        XK0 = LR * beta_func(a,Nf,order)
        XK1 = LR * beta_func(a + 0.5 * XK0,Nf,order)
        XK2 = LR * beta_func(a + 0.5 * XK1,Nf,order)
        XK3 = LR * beta_func(a + XK2,Nf,order)
        a+= (XK0 + 2.* XK1 + 2.* XK2 + XK3) * 0.166666666666666
    return a

order=1
ae=par.alfa0/(4*np.pi)
#am=evolve_a(par.me2,ae,par.mm2,1,order)
#at=evolve_a(par.mm2,am,par.mt2,2,order)

@jit(nopython=True)
def get_a(Q2):
    Q20,a0,Nf=par.me2,ae,1
    return evolve_a(Q20,a0,Q2,Nf,order)

@jit(nopython=True)
def get_alfa(Q2):
    return get_a(Q2)*4*np.pi
