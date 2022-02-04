import numpy as np
from numba import jit
from stflib import tmd_collins,tmd_transversity
from stflib.core import charge_conj, p2n
from stflib.params import e2
import stflib.params as par
    
@jit(nopython=True)
def _get_FUT(x,z,Q2,pT,tar,had,F,D,w_tar,w_had):

    Mh=par.Mpi
 
    if had.endswith('+'):
        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        K = 2 * x * z * pT * Mh / wq
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2*K*F*D*gauss)

    elif had.endswith('-'):
        D=charge_conj(D)
        w_had=charge_conj(w_had)
        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        K = 2 * x * z * pT * Mh / wq
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2*K*F*D*gauss)

    elif had.endswith('0'):

        Dp=D
        Dm=charge_conj(D)
        w_hadp=w_had
        w_hadm=charge_conj(w_had)

        wqp = z**2 * np.abs(w_tar) + np.abs(w_hadp)
        wqm = z**2 * np.abs(w_tar) + np.abs(w_hadm)
        gaussp = np.exp(-pT**2 / wqp) / (np.pi * wqp)
        gaussm = np.exp(-pT**2 / wqm) / (np.pi * wqm)
        Kp = 2 * x * z * pT * Mh / wqp
        Km = 2 * x * z * pT * Mh / wqm

        FUTp=np.sum(e2*Kp*F*Dp*gaussp)
        FUTm=np.sum(e2*Km*F*Dm*gaussm)
        return 0.5*(FUTp+FUTm)

@jit(nopython=True)
def get_FUT(x,z,Q2,pT,tar,had):

    # get collinear parts (proton and positive hadrons)
    F = tmd_transversity.get_C(x, Q2)
    D = tmd_collins.get_C(z, Q2)
    F[0],D[0]=0,0  # set glue to zero

    # get widths (proton and positive hadrons)
    w_tar=tmd_transversity.get_widths(Q2)
    w_had=np.abs(tmd_collins.get_widths(Q2))
        
    # build structure function
    K = x
    if tar=='p':

        return _get_FUT(x,z,Q2,pT,tar,had,F,D,w_tar,w_had)

    elif tar=='n':

        F=p2n(F)
        w_tar=p2n(w_tar)
        return  _get_FUT(x,z,Q2,pT,tar,had,F,D,w_tar,w_had)
  
    elif tar=='d':

        return 0.5*(get_FUT(x,z,Q2,pT,'p',had)+get_FUT(x,z,Q2,pT,'n',had))
