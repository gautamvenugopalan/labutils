# Functions for lossy homodyne spectral density
# These are Equations 5.7--5.12 from Buonanno and Chen 2001, https://doi.org/10.1103/PhysRevD.64.042006
# IFO is a class object with the following parameters:
#	 IFO.Pin  ------ Input power to the back of PRM
#	 IFO.Litm ------ ITM loss
#	 IFO.Letm ------ ETM loss
#	 IFO.Lsrm ------ SRM loss
#	 IFO.T_I ------- ITM Power Transmissivity
#	 IFO.T_E ------- ETM Power Transmissivity
#	 IFO.rs -------- SRM Amplitude Reflectivity
#	 IFO.ts -------- SRM Amplitude Transmissivity

# Imports
import numpy as np
import scipy.constants as scc


def armRefl(IFO, printInfo=False):
    '''
    Function to calculate the reflectivity of a (lossy) two mirror
    Fabry Perot cavity. Takes as input the (power) transmissivity
    and loss of the input and output couplers.
    '''
    ti = IFO.ti
    te = IFO.te
    ri = IFO.ri
    re = IFO.re
    etai = IFO.Litm
    etae = IFO.Letm
    refl = (re*(ti**2 + ri**2) - ri)/(1 - ri*re)
    if printInfo:
        print('Ri={} %, Ti={} %, Re={} %, Te={} ppm, loss_i={} ppm, loss_e={} ppm'.format(round(100*(ri**2),3),round(100*Ti,3), 
            round(100*(re**2),3), round(1e6*Te,3), 1e6*etai, 1e6*etae))
        print('Amplitude reflectivity is {}'.format(refl))
    return refl
def PRG(IFO,  printInfo=False):
    '''
    Function to calculate the RPG for a given arm reflectivity
    and PRM (power) transmissivity.
    '''
    tp = IFO.tp
    rp = IFO.rp
    rarm = armRefl(IFO)
    gp = tp / (1-rp*rarm)
    if printInfo:
        print('PRG is {} for Tp= {} %, Rarm = {} %'.format(round(gp**2, 3), round(100*Tp,3), round(100*rarm**2,3)))
    return gp**2
def armPole(IFO):
    '''
    Computes the pole frequency [Hz] of a single arm cavity 
    '''
    return scc.c*np.log(1/(IFO.ri*IFO.re))/2/IFO.larm/2/np.pi
def hSQL(IFO,ff):
    '''
    Computes the SQL [1/rtHz] for GW detection at frequency ff [Hz]
    '''
    omega = 2*np.pi*ff
    return np.sqrt(8*scc.hbar/IFO.m_TM/omega**2/IFO.larm**2)
def I0SQL(IFO):
    '''
    Computes power required by CONVENTIONAL (i.e. non SR) IFO to reach SQL 
    at a GW frequency equal to the arm cavity pole. 
    '''
    gam = 2*np.pi*armPole(IFO)
    return IFO.m_TM * IFO.larm**2 * gam**4 / 4 / (2*np.pi*scc.c/1064e-9)

def kappa(IFO,ff):
    '''
    Calculates the optomechanical coupling 
    at a frequency ff [Hz]
    '''
    omega = 2*np.pi*ff
    wp = 2*np.pi*armPole(IFO)
    isq = I0SQL(IFO)
    ra = armRefl(IFO)
    gp = PRG(IFO)
    I0 = IFO.Pin * gp
    return 2*(I0/isq)*wp**4 / (omega**2 * (wp**2 + omega**2))
def beta(IFO,ff):
    '''
    Calculates the net phase gain for the sideband in the arm cavity
    '''
    omega = 2*np.pi*ff
    wp = 2*np.pi*armPole(IFO)
    return np.arctan(omega/wp)
def eps(IFO):
    '''
    Calculates the quantity epsilon given ITM/ETM losses and ITM power transmissivity.
    '''
    return 2*(IFO.Litm+IFO.Letm)/IFO.T_I
def Cmatrix_lossy(IFO,ph,ff,Lpd=0.1):
    '''
    Calculates the C matrix for a lossy SR interferometer, assuming readout chain losses of Lpd.
    Returns matrix elements row-wise.
    '''
    phi = np.deg2rad(ph) 
    bb = beta(IFO,ff)
    K = kappa(IFO,ff)
    Lsr = 2*IFO.Lsrm
    epsilon = eps(IFO)
    rs = IFO.rs
    ts = IFO.ts
    C11 = np.sqrt(1.-Lpd)*( (1+rs**2)*(np.cos(2*phi) + K*np.sin(2*phi)/2) - 2*rs*np.cos(2*bb) - (epsilon/4)*(-2*rs*(1+np.exp(2j*bb))**2 + 
        4*(1+rs**2)*np.cos(bb)**2 * np.cos(2*phi) + (3+np.exp(2j*bb))*K*(1+rs**2)*np.sin(2*phi)) + 
                           Lsr*(np.exp(2j*bb) - ((1+rs**2)/2)*(np.cos(2*phi) + K*np.sin(2*phi)/2)) )
    C22 = C11
    C12 = ts**2 * np.sqrt(1-Lpd)*( -(np.sin(2*phi) + K*np.sin(phi)**2) + 
                             (epsilon*np.sin(phi)/2)*(K*np.sin(phi)*(3+np.exp(2j*bb)) + 4*np.cos(phi)*np.cos(bb)**2) 
                            + (Lsr/2)*(np.sin(2*phi) + K*np.sin(phi)**2) )
    C21 = ts**2 * np.sqrt(1-Lpd)*( (np.sin(2*phi) - K*np.cos(phi)**2) + 
                             (epsilon*np.cos(phi)/2)*(K*np.cos(phi)*(3+np.exp(2j*bb)) - 4*np.sin(phi)*np.cos(bb)**2) 
                            + (Lsr/2)*(-np.sin(2*phi) + K*np.cos(phi)**2) )
    return C11,C12,C21,C22
def Dmatrix_lossy(IFO,ph,ff,Lpd=0.1):
    '''
    Calculates the D matrix for a lossy SR interferometer, assuming readout chain losses of Lpd.
    Returns matrix elements row-wise.
    '''
    phi = np.deg2rad(ph) 
    bb = beta(IFO,ff)
    K = kappa(IFO,ff)
    Lsr = 2*IFO.Lsrm
    epsilon = eps(IFO)
    rs = IFO.rs
    ts = IFO.ts
    D1 = np.sqrt(1.-Lpd)* ( -(1+rs*np.exp(2j*bb))*np.sin(phi) + (epsilon*np.sin(phi)/4)*(3 + rs + 2*rs*np.exp(4j*bb) + np.exp(2j*bb)*(1 + 5*rs)) 
                          + Lsr*np.exp(2j*bb)*rs*np.sin(phi)/2 )
    D2 = np.sqrt(1.-Lpd)* ( -(-1+rs*np.exp(2j*bb))*np.cos(phi) + (epsilon*np.cos(phi)/4)*(-3 + rs + 2*rs*np.exp(4j*bb) + np.exp(2j*bb)*(-1 + 5*rs)) 
                          + Lsr*np.exp(2j*bb)*rs*np.cos(phi)/2 )
    return D1,D2
def Pmatrix_lossy(IFO,ph,ff,Lpd=0.1):
    '''
    Calculates the P matrix for a lossy SR interferometer, assuming readout chain losses of Lpd.
    Returns matrix elements row-wise.
    '''
    phi = np.deg2rad(ph) 
    bb = beta(IFO,ff)
    K = kappa(IFO,ff)
    Lsr = 2*IFO.Lsrm
    epsilon = eps(IFO)
    rs = IFO.rs
    ts = IFO.ts
    P11 = 0.5*np.sqrt(1-Lpd)*np.sqrt(Lsr)*ts*( -2*rs*np.exp(2j*bb) + 2*np.cos(2*phi) + K*np.sin(2*phi) )
    P22 = P11
    P12 = -np.sqrt(1-Lpd)*np.sqrt(Lsr)*ts*np.sin(phi)*(2*np.cos(phi) + K*np.sin(phi))
    P21 = np.sqrt(1-Lpd)*np.sqrt(Lsr)*ts*np.cos(phi)*(2*np.sin(phi) - K*np.cos(phi))
    return P11,P12,P21,P22
def Qmatrix_lossy(IFO,ph,ff,Lpd=0.1):
    '''
    Calculates the Q matrix for a lossy SR interferometer, assuming readout chain losses of Lpd.
    Returns matrix elements row-wise.
    '''
    phi = np.deg2rad(ph) 
    bb = beta(IFO,ff)
    K = kappa(IFO,ff)
    Lsr = 2*IFO.Lsrm
    epsilon = eps(IFO)
    rs = IFO.rs
    ts = IFO.ts
    Q11 = np.sqrt(Lpd) * ( np.exp(-2j*bb) + rs**2*np.exp(2j*bb) - rs*(2*np.cos(2*phi) + K*np.sin(2*phi))
    +(epsilon*rs/2)*(np.exp(-2j*bb)*np.cos(2*phi) + np.exp(2j*bb)*(-2*rs -2*rs*np.cos(2*bb) +np.cos(2*phi) + K*np.sin(2*phi)) + 
    2*np.cos(2*phi) + 3*K*np.sin(2*phi)) - (Lsr*rs/2)*(2*rs*np.exp(2j*bb) -2*np.cos(2*phi) -K*np.sin(2*phi)) )
    Q22 = Q11
    Q12 = np.zeros(len(Q11))
    Q21 = np.zeros(len(Q11))
    return Q11, Q12, Q21, Q22
def Nmatrix_lossy(IFO,ph,ff,Lpd=0.1):
    '''
    Calculates the N matrix for a lossy SR interferometer, assuming readout chain losses of Lpd.
    Returns matrix elements row-wise.
    '''
    phi = np.deg2rad(ph) 
    bb = beta(IFO,ff)
    K = kappa(IFO,ff)
    Lsr = 2*IFO.Lsrm
    epsilon = eps(IFO)
    rs = IFO.rs
    ts = IFO.ts
    N11 = ts*np.sqrt((1.-Lpd)*epsilon/2)* ( K*np.sin(phi)*(1+rs*np.exp(2j*bb)) + 2*np.cos(bb)*( np.exp(-1j*bb)*np.cos(phi)
        -rs*np.exp(1j*bb)*(np.cos(phi) + K*np.sin(phi))))
    N22 = -ts*np.sqrt((1.-Lpd)*2*epsilon) * (-np.exp(-1j*bb) + rs*np.exp(1j*bb)) * np.cos(bb) * np.cos(phi)
    N12 = -ts*np.sqrt((1.-Lpd)*2*epsilon) * (np.exp(-1j*bb) + rs*np.exp(1j*bb)) * np.cos(bb) * np.sin(phi)
    N21 = ts*np.sqrt((1.-Lpd)*epsilon/2)* ( -K*np.cos(phi)*(1+rs) + 2*np.cos(bb*(np.exp(-1j*bb) + rs*np.exp(1j*bb)))*np.cos(bb)*np.sin(phi) )
    return N11,N12,N21,N22
def homodyneASD_lossy(IFO,ph,zz,ff,Lpd=0.1):
    '''
    Computes the ASD for homodyne readout at homodyne angle zeta (degrees)
    for an SRC tuning of ph [deg] for a LOSSY SR interferometer
    '''
    zeta = np.deg2rad(zz)
    K = kappa(IFO,ff)
    C11, C12, C21, C22 = Cmatrix_lossy(IFO,ph,ff,Lpd)
    D1, D2 = Dmatrix_lossy(IFO,ph,ff,Lpd)
    P11, P12, P21, P22 = Pmatrix_lossy(IFO,ph,ff,Lpd)
    Q11, Q12, Q21, Q22 = Qmatrix_lossy(IFO,ph,ff,Lpd)
    N11, N12, N21, N22 = Nmatrix_lossy(IFO,ph,ff,Lpd)
    hsql = hSQL(IFO,ff)
    num = np.abs(C11*np.sin(zeta) + C21*np.cos(zeta))**2 + np.abs(C12*np.sin(zeta) + C22*np.cos(zeta))**2
    num = num + np.abs(P11*np.sin(zeta) + P21*np.cos(zeta))**2 + np.abs(P12*np.sin(zeta) + P22*np.cos(zeta))**2
    num = num + np.abs(Q11*np.sin(zeta) + Q21*np.cos(zeta))**2 + np.abs(Q12*np.sin(zeta) + Q22*np.cos(zeta))**2
    num = num + np.abs(N11*np.sin(zeta) + N21*np.cos(zeta))**2 + np.abs(N12*np.sin(zeta) + N22*np.cos(zeta))**2
    num = num * hsql**2
    den = 2 * K * IFO.ts**2 * np.abs(D1*np.sin(zeta) + D2*np.cos(zeta))**2
    return np.sqrt(num/den)
def vacuumNoise(IFO,ph,zz,ff,Lpd=0.1):
    '''
    Computes the vacuum noise at homodyne angle zeta (degrees)
    for an SRC tuning of ph [deg] for a LOSSY SR interferometer
    '''
    zeta = np.deg2rad(zz)
    K = kappa(IFO,ff)
    C11, C12, C21, C22 = Cmatrix_lossy(IFO,ph,ff,Lpd)
    mm = M_lossy(IFO, ph, ff, Lpd)
    num = np.abs(C11*np.sin(zeta) + C21*np.cos(zeta))**2 + np.abs(C12*np.sin(zeta) + C22*np.cos(zeta))**2
    return num / (np.abs(mm**2))
def M_lossy(IFO,ph,ff,Lpd=0.1):
    '''
    Calculates M for a lossy SR interferometer, assuming readout chain losses of Lpd.
    '''
    phi = np.deg2rad(ph) 
    bb = beta(IFO,ff)
    K = kappa(IFO,ff)
    Lsr = 2*IFO.Lsrm
    epsilon = eps(IFO)
    rs = IFO.rs
    ts = IFO.ts
    M = (1 + rs**2*np.exp(4j*bb) - 2*rs*np.exp(2j*bb)*(np.cos(2*phi) + (K/2)*np.sin(2*phi)) +
         Lsr*rs*np.exp(2j*bb)*(-rs*np.exp(2j*bb) + np.cos(2*phi) + (K/2)*np.sin(2*phi)) +
         epsilon*rs*np.exp(2j*bb)*(2*np.cos(bb)**2 *(-rs*np.exp(2j*bb) + np.cos(2*phi))
                                      + (K/2)*(3+np.exp(2j*bb))*np.sin(2*phi)))
    return M
