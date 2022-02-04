import sys,os
import numpy as np
from numba import jit
from stflib import interp
import params as par

#--Collinear parts
L =open(os.environ['stflib']+'/tables/CJ15lo_0000.dat').readlines()
xx=np.array([float(val) for val in L[3].split()])
QQ=np.array([float(val) for val in L[4].split()])
iflav=[int(val) for val in L[5].split()]
L=[l.split() for l in L[6:] if '---' not in l]
data=[[float(val) for val in l] for l in L]
data=np.transpose(data)

nx=len(xx)
nQ=len(QQ)
nf=len(iflav)

table=np.zeros((len(iflav),nx,nQ))
for i in range(nf): 
    cnt=-1
    for ix in range(nx):
        for iQ in range(nQ):
            cnt+=1
            table[i,ix,iQ]=data[i][cnt]

@jit(nopython=True)
def get_f(x,Q2):
    #iflav=[-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21]
    #        0   1   2   3   4  5  6  7  8, 9, 10
    #self.IFLAV=[21,2,-2,1,-1,3,-3,4,-4,5,-5]
    iflav=[10,6,3,5,4,7,2,8,1,9,0]
    return np.array([interp.intepolate2d(xx,QQ, table[_,:,:],x,np.sqrt(Q2),5,5)/x for _ in iflav])

@jit(nopython=True)
def get_C(x, Q2):
    return get_f(x,Q2)
    
#--TMD parts
data=np.load(os.environ['stflib']+'/jam20_params.npy',allow_pickle=True).item(0)        
widths1 = data['pdf']['widths1']
widths2 = data['pdf']['widths2']

@jit(nopython=True)
def get_widths(Q2):
    s=np.log(Q2/par.Q20 )
    return np.abs(widths1+s*widths2)
