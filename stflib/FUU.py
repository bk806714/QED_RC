import numpy as np
from numba import jit
from stflib import tmd_ff,tmd_pdf
from stflib.core import charge_conj, p2n
from stflib.params import e2


@jit(nopython=True)
def _get_FUU(x,z,Q2,pT,tar,had,F,D,w_tar,w_had):

    K = x

    if had.endswith('+'):
        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
        gauss = np.exp(-pT**2 / wq) / (np.pi * wq)
        return np.sum(e2*K*F*D*gauss)

    elif had.endswith('-'):
        D=charge_conj(D)
        w_had=charge_conj(w_had)
        wq = z**2 * np.abs(w_tar) + np.abs(w_had)
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

        FUUp=np.sum(e2*K*F*Dp*gaussp)
        FUUm=np.sum(e2*K*F*Dm*gaussm)
        return 0.5*(FUUp+FUUm)

@jit(nopython=True)
def get_FUU(x,z,Q2,pT,tar,had):

    # get collinear parts (proton and positive hadrons)
    F = tmd_pdf.get_C(x, Q2)
    D = tmd_ff.get_C(z, Q2)
    F[0],D[0]=0,0  # set glue to zero

    # get widths (proton and positive hadrons)
    w_tar=tmd_pdf.get_widths(Q2)
    w_had=np.abs(tmd_ff.get_widths(Q2))


    # build structure function

    if tar=='p':

        return _get_FUU(x,z,Q2,pT,tar,had,F,D,w_tar,w_had)

    elif tar=='n':

        F=p2n(F)
        w_tar=p2n(w_tar)
        return  _get_FUU(x,z,Q2,pT,tar,had,F,D,w_tar,w_had)

    elif tar=='d':

        return 0.5*(get_FUU(x,z,Q2,pT,'p',had)+get_FUU(x,z,Q2,pT,'n',had))

    elif tar=='He':

        return 1.5*(2*get_FUU(x,z,Q2,pT,'p',had)+get_FUU(x,z,Q2,pT,'n',had))
