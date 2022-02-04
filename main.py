import sys,os,time
import numpy as np
from numba import jit
import vegas

import params as par
import alfa
import dglap
from tools import lprint,checkdir

os.environ["stflib"] = "./stflib"
from stflib import FUU,FUTsiv,FUTcol


@jit(nopython=True)
def get_int0(xi,zeta,kinematics,iupol,isiv,icol,rotation=1):
    """
    _..._ stands for non RC quantities 
    """
    M  =par.M 
    M_h=par.Mpi
    
    #--retrieve & build non RC quantities
    _s_        =kinematics[0]
    _Q_        =kinematics[1]
    _x_        =kinematics[2]
    _y_        =kinematics[3]
    _z_        =kinematics[4]
    _PhT_      =kinematics[5]
    _cos_phih_ =kinematics[6]
    _cos_phis_ =kinematics[7]
    _sin_phih_ =kinematics[8]
    _sin_phis_ =kinematics[9]
    _ST_       =kinematics[10]

    _l1T2_   = _Q_**2*(-M**2*_x_**2*_y_**2 - _Q_**2*_y_ + _Q_**2)/(_y_**2*(4*M**2*_x_**2 + _Q_**2))
    _l1T_     = np.sqrt(_l1T2_)

    _qdotPh_ = _Q_*(_Q_**3*_z_ - np.sqrt(-(4*M**2*_x_**2 + _Q_**2)*(4*M**2*M_h**2*_x_**2 + 4*M**2*_PhT_**2*_x_**2 - _Q_**4*_z_**2)))/(4*M**2*_x_**2)
    _l1dotPh_= (-4*M**2*_PhT_*_l1T_*_x_**2*_y_*_cos_phih_ + 2*M**2*_x_**2*_y_*_qdotPh_ - \
               _PhT_*_Q_**2*_l1T_*_y_*_cos_phih_ - _Q_**4*_y_*_z_/2 + _Q_**4*_z_ + _Q_**2*_qdotPh_)/(_y_*(4*M**2*_x_**2 + _Q_**2)) 
    _l2dotPh_=_l1dotPh_-_qdotPh_ 

    _qdotS_  = -_Q_*np.sqrt(-(_ST_ - 1)*(_ST_ + 1)*(4*M**2*_x_**2 + _Q_**2))/(2*M*_x_)    
    _l1dotS_ = (-4*M**2*_ST_*_l1T_*_x_**2*_y_*_cos_phis_ + 2*M**2*_x_**2*_y_*_qdotS_\
               - _Q_**2*_ST_*_l1T_*_y_*_cos_phis_ + _Q_**2*_qdotS_)/(_y_*(4*M**2*_x_**2 + _Q_**2))  

    _l2dotS_=_l1dotS_-_qdotS_ 

    _eps_Pl1l2Ph_=-_PhT_*_Q_**2*_l1T_*np.sqrt(4*M**2*_x_**2/_Q_**2 + 1)*_sin_phih_/(2*_x_)
    _eps_Pl1l2S_ =-_ST_ *_Q_**2*_l1T_*np.sqrt(4*M**2*_x_**2/_Q_**2 + 1)*_sin_phis_/(2*_x_)


    #--build RC kinematic quantities
    s        = xi*_s_
    x        = _Q_**2*_x_*xi*_y_/(_Q_**2*(xi*zeta + _y_ - 1))
    y        = (xi*zeta + _y_ - 1)/(xi*zeta)
    z        = _y_*_z_*zeta/(xi*zeta + _y_ - 1)
    Q2       = _Q_**2*xi/zeta
    Q        = np.sqrt(Q2)        
    gam2     = (2*M*x)**2/Q2
    gam      = np.sqrt(gam2)

    l1T2     = Q**2*(-M**2*x**2*y**2 - Q**2*y + Q**2)/(y**2*(4*M**2*x**2 + Q**2))
    l1T      = np.sqrt(l1T2)

    l1dotPh  = xi*_l1dotPh_
    l2dotPh  = _l2dotPh_/zeta
    l1dotS   = xi*_l1dotS_
    l2dotS   = _l2dotS_/zeta

    eps_Pl1l2Ph = (xi/zeta)*_eps_Pl1l2Ph_
    eps_Pl1l2S  = (xi/zeta)*_eps_Pl1l2S_

    qdotPh   = l1dotPh-l2dotPh 
    qdotS    = l1dotS-l2dotS 

    PhT2     = (-4*M**2*M_h**2*Q**2*x**2 - 4*M**2*x**2*qdotPh**2 - M_h**2*Q**4 + Q**6*z**2 + 2*Q**4*z*qdotPh)/(Q**2*(4*M**2*x**2 + Q**2))
    PhT      = np.sqrt(np.abs(PhT2))
    if  rotation==0: PhT=_PhT_

    ST2      = (4*M**2*Q**2*x**2 - 4*M**2*x**2*qdotS**2 + Q**4)/(Q**2*(4*M**2*x**2 + Q**2))
    ST       = np.sqrt(np.abs(ST2))
    Spar     = 2*M*x*qdotS/(Q**2*np.sqrt(4*M**2*x**2/Q**2 + 1))

    cos_phih = (4*M**2*x**2*y*qdotPh - Q**4*y*z + 2*Q**2*(Q**2*z + qdotPh)- 2*y*(4*M**2*x**2 + Q**2)*l1dotPh)/(2*PhT*l1T*y*(4*M**2*x**2 + Q**2))                    
    cos_phis = (2*M**2*x**2*y*qdotS + Q**2*qdotS-y*(4*M**2*x**2 + Q**2)*l1dotS)/(ST*l1T*y*(4*M**2*x**2 + Q**2))

    sin_phih = -2*x*eps_Pl1l2Ph/(PhT*Q**2*l1T*np.sqrt(4*M**2*x**2/Q**2 + 1))
    sin_phis = -2*x*eps_Pl1l2S /( ST*Q**2*l1T*np.sqrt(4*M**2*x**2/Q**2 + 1))

    if iupol==1: FUU_ = FUU.get_FUU(x,z,Q2,PhT,'p','pi+')
    else:        FUU_ = 0
        
    if isiv==1: FUTsiv_ = FUTsiv.get_FUT(x,z,Q2,PhT,'p','pi+')
    else:       FUTsiv_ = 0
    
    if icol==1: FUTcol_ = FUTcol.get_FUT(x,z,Q2,PhT,'p','pi+')    
    else:       FUTcol_ = 0
    
    eps  = (1-y-0.25*gam2*y**2)/(1-y+0.5*y**2+0.25*gam2*y**2)
    jac  = x/_x_/xi/zeta 
    norm = alfa.get_alfa(Q2)**2/(x*y*Q2) * y**2/2/(1-eps) * (1+gam2/2/x) 

    phase_sivers  = sin_phih*cos_phis - sin_phis*cos_phih
    phase_collins = sin_phih*cos_phis + sin_phis*cos_phih


    out=jac*norm*(FUU_ + ST*phase_sivers*FUTsiv_ + ST*eps*phase_collins*FUTcol_)
    #print(out)
    return out
    
@jit(nopython=True)
def get_xi_min(zeta,kinematics):
    y,x,z=kinematics[3],kinematics[2],kinematics[4]
    return max([(1-y+z*zeta*y)/zeta,(1-y)/(zeta-x*y)])

@jit(nopython=True)
def get_zeta_min(kinematics):
    y,x,z=kinematics[3],kinematics[2],kinematics[4]
    return max([(1-y)/(1-z*y),1-y+x*y])   

@jit(nopython=True)
def get_box(xi,zeta,kinematics,F,jac_xi,mu2,iupol,isiv,icol,rotation=1):
    xi_min =get_xi_min(zeta,kinematics)
    LDF    =dglap.get_LDF(xi,F,mod=False)
    LDF2   =dglap.get_LDF(xi_min,F,mod=True)
    int0_xi_zeta =get_int0(xi,zeta,kinematics,iupol,isiv,icol,rotation=1)
    int0_one_zeta=get_int0(1,zeta,kinematics,iupol,isiv,icol,rotation=1)
    return LDF*(int0_xi_zeta-int0_one_zeta)*jac_xi +int0_one_zeta*xi_min*LDF2

@jit(nopython=True)
def get_integrand0(X,kinematics,F,D,iupol,isiv,icol,rotation=1):
    mu2=kinematics[1]**2
    zeta_min=get_zeta_min(kinematics)
    jac_zeta=1-zeta_min
    zeta=zeta_min + X[0]*jac_zeta

    xi_min=get_xi_min(zeta,kinematics)
    jac_xi=1-xi_min
    xi=xi_min + X[1]*jac_xi
    

    LFF =dglap.get_LFF(zeta,D,mod=False)
    LFF2=dglap.get_LFF(zeta_min,D,mod=True)
        
    box_xi_zeta=get_box(xi,zeta,kinematics,F,jac_xi,mu2,iupol,isiv,icol,rotation)
    box_xi_one =get_box(xi,1,kinematics,F,jac_xi,mu2,iupol,isiv,icol,rotation)

    theta=1
    if xi<get_xi_min(1,kinematics): theta=0
    integrand=LFF*(box_xi_zeta-box_xi_one*theta)*jac_zeta + box_xi_one*theta*zeta_min*LFF2
    #print(box_xi_zeta,box_xi_one,integrand)

    return integrand

@jit(nopython=True)
def get_integrand1(X,kinematics,F,D,iphih,iphis,iupol,isiv,icol,qed_order=1,rotation=1):
    phis=-np.pi + 2*np.pi*X[2]
    phih=-np.pi + 2*np.pi*X[3]
    kinematics[6] = np.cos(phih)
    kinematics[7] = np.cos(phis)
    kinematics[8] = np.sin(phih)
    kinematics[9] = np.sin(phis)
    jac   = (2*np.pi)**2
    phase = np.sin(iphih*phih+iphis*phis)*jac
    if qed_order==0:
        int0=get_int0(1,1,kinematics,iupol,isiv,icol,rotation)
        out=int0*phase

    elif qed_order==1:
        integrand0=get_integrand0(X,kinematics,F,D,iupol,isiv,icol,rotation)
        out=integrand0*phase

    if np.isnan(out): out=0
    return out

def get_xsec(kinematics
             ,iphih,iphis
             ,iupol,isiv,icol
             ,qed_order
             ,dglap_order,rotation
             ,nbiter,nburn,neiter,neval,ftol=1e-8):
    Q=kinematics[1]
    mu2=Q**2
    F,D=dglap.evolve(mu2,dglap_order)
    trial=1
    
    checkdir('data')
    fname='data/'
    fname+='rs:%f'%rs
    fname+='-x:%f'%x
    fname+='-z:%f'%z
    fname+='-Q:%f'%Q
    fname+='-qT_Q:%f'%qT_Q
    fname+='-iphih:%d'%iphih
    fname+='-iphis:%d'%iphis
    fname+='-iupol:%d'%iupol
    fname+='-isiv:%d'%isiv
    fname+='-icol:%d'%icol
    fname+='-qed_order:%d'%qed_order
    fname+='-dglap_order:%d'%dglap_order
    fname+='-rotation:%d'%rotation
    fname+='.npy'

    #imin=5 ;imax=10
    imin=10;imax=100000

    results=[]

    while 1:    
        integ = vegas.Integrator([[0, 1] for _ in range(4)])
        func=lambda X:get_integrand1(X,kinematics,F,D,iphih,iphis,iupol,isiv,icol,qed_order,rotation)
        integ(func, nitn=nbiter, neval=nburn)        # burn
        result=integ(func, nitn=neiter, neval=neval) # final
        results.append(result.val)
        if len(results)<=imin:
            lprint('size=%d'%len(results))
            continue
        mean=np.mean(results)
        std =np.std(results)
        np.save(fname,[mean,std])

        lprint('trial=%d size=%d mean=%0.2e std=%0.2e rel=%0.2e'%(trial,len(results),mean,std,std/mean))
        if (std/mean)<ftol:break        
        if len(results)>imax:
            results=[]
            nburn *=10
            neiter*=10
            trial +=1
            print()
            

            
    return mean,std

if __name__=="__main__":
    
    rs    = float(sys.argv[1])
    x     = float(sys.argv[2])
    z     = float(sys.argv[3])
    Q     = float(sys.argv[4])
    qT_Q  = float(sys.argv[5])
    
    iphih  = int(sys.argv[6])
    iphis  = int(sys.argv[7])
    
    iupol  = int(sys.argv[8])
    isiv   = int(sys.argv[9])
    icol   = int(sys.argv[10])
    
    qed_order   = int(sys.argv[11])
    dglap_order = int(sys.argv[12])
    rotation    = int(sys.argv[13])
    
    
    
    y    = Q**2/rs**2/x
    qT   = qT_Q*Q
    PhT  = z*qT
    ST   = 1
    kinematics=np.zeros(11)
    kinematics[0] = rs**2
    kinematics[1] = Q
    kinematics[2] = x
    kinematics[3] = y
    kinematics[4] = z
    kinematics[5] = PhT
    kinematics[10]= ST


    
    nbiter,nburn = 1,1000
    neiter,neval = 1,10000
    
    ftol=1e-2
    get_xsec(kinematics
             ,iphih,iphis
             ,iupol,isiv,icol
             ,qed_order
             ,dglap_order,rotation
             ,nbiter,nburn,neiter,neval,ftol)

