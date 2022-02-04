import os
import numpy as np
from numba import jit
import stflib.params as par
from stflib import core

# free parameters
data=np.load(os.environ['stflib']+'/jam20_params.npy',allow_pickle=True).item(0)        
widths1 = data['transversity']['widths1']
widths2 = data['transversity']['widths2']
shape1  = data['transversity']['shape1']
shape2  = data['transversity']['shape2']

@jit(nopython=True)
def get_C(x, Q2):
    return core.get_collinear(x, Q2,shape1,shape2)

@jit(nopython=True)
def get_widths(Q2):
    s=np.log(Q2/par.Q20 )
    return np.abs(widths1+s*widths2)
