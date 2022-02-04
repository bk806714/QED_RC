import numpy as np
from numba import jit
import params as par
from stflib.special import beta

@jit(nopython=True)
def get_s(Q2):
    lam2 = par.lam2
    Q20  = par.Q20
    return np.log(np.log(Q2/lam2)/np.log(Q20/lam2))

@jit(nopython=True)
def __get_shape(x,p):
    return p[0]*x**p[1]*(1-x)**p[2]*(1+p[3]*x+p[4]*x**2)

@jit(nopython=True)
def _get_shape(x,p,s):
    N=p[0] + p[5] * s
    a=p[1] + p[6] * s
    b=p[2] + p[7] * s
    c=p[3] + p[8] * s
    d=p[4] + p[9] * s
    n= 1
    return __get_shape(x,np.array([N/n,a,b,c,d]))

@jit(nopython=True)
def get_shape(x,Q2,p1,p2):
    s=get_s(Q2)
    shape=_get_shape(x,p1,s)
    shape2=_get_shape(x,p2,s)
    if p2[0]!=0: shape=shape*(1+shape2)
    n1 = beta(p1[1]+2,p1[2]+1) + p1[3]*beta(p1[1]+3,p1[2]+1) + p1[4]*beta(p1[1]+4,p1[2]+1)
    n2 = p2[0] * (beta(p1[1]+p2[1]+2,p1[2]+p2[2]+1) + p1[4]*p2[4]*beta(p1[1]+p2[1]+6,p1[2]+p2[2]+1))
    n3 = p2[0] * (p1[3]+p2[3])*beta(p1[1]+p2[1]+3,p1[2]+p2[2]+1)
    n4 = p2[0] * (p1[3]*p2[4]+p2[3]*p1[4])*beta(p1[1]+p2[1]+5,p1[2]+p2[2]+1)
    n5 = p2[0] * (p1[3]*p2[3]+p1[4]+p2[4])*beta(p1[1]+p2[1]+4,p1[2]+p2[2]+1)
    n=n1+n2+n3+n4+n5
    shape=shape/n
    return shape

@jit(nopython=True)
def get_collinear(x,Q2,shape1,shape2):
    N = np.zeros(11)
    for i in range(11):
        N[i] = get_shape(x,Q2,shape1[i],shape2[i])
    return N

@jit(nopython=True)
def charge_conj(p):
    out=np.zeros(len(p))
    cnt=0
    for i in [0,2,1,4,3,6,5,8,7,10,9]:
        out[cnt]=p[i]
        cnt+=1
    return out

@jit(nopython=True)
def p2n(p):
    out=np.zeros(len(p))
    cnt=0
    for i in [0,3,4,1,2,5,6,7,8,9,10]:
        out[cnt]=p[i]
        cnt+=1
    return out

