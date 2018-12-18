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

class ifo:
    '''
    Class object for the IFO object, with which to calculate various quantum noises.
    All default values are taken from the GWINC Voyager parameter file.

    Attributes:
    ============
    T_I: float
        ITM power transmissivity. Defaults to 0.12436875%.
    T_E: float
        ETM power transmissivity. Defaults to 5ppm.
    T_P: float
        PRM power transmissivity. Defaults to 3%.
    T_S: float
        SRM power transmissivity. Defaults to 5.86%.
    Litm: float
        ITM loss. Defaults to 10 ppm.
    Letm: float
        ETM loss. Defaults to 10 ppm.
    Lprm: float
        Loss in the PRC. Defaults to 500 ppm.
    Lsrm: float
        SRM loss. Defaults to 500 ppm.
    larm: float
        Arm length in meters. Defaults to 3.995 km.
    Pin: float
        Input power on the back of PRM. Defaults to 144.6848 W.
    m_TM: float
        Mass of the test masses. Defaults to 200 kg.
    '''
    def __init__(self, T_I=0.12436875e-2, T_E=5e-6, T_P=3e-2, T_S=5.86e-2, Litm=10e-6, Letm=10e-6, Lprm=500e-6, Lsrm=500e-6, larm=3995, Pin=144.6848, m_TM=200):
        self.__dict__.update(locals()) # This hack assigns all the arguments passed to init.
        self.ti = np.sqrt(self.T_I)
        self.te = np.sqrt(self.T_E)
        self.ri = np.sqrt(1 - self.T_I - self.Litm)
        self.re = np.sqrt(1 - self.T_E - self.Letm)
        self.tp = np.sqrt(self.T_P)
        self.ts = np.sqrt(self.T_S)
        self.rp = np.sqrt(1 - self.T_P - self.Lprm)
        self.rs = np.sqrt(1 - self.T_S - self.Lsrm)
        return

# Class definition for the amplifier object
class amp:
    '''
    Class object for the triangular amplifier cavity geometry.
    
    Attributes:
    -----------
    Ti: float
      Input coupler (Power) transmissivity. Defaults to 0.35%.
    Te: float
      Power transmissivities of the other two mirrors in the ring. Defaults to 5 ppm.
    Lrt: float
      Round-trip loss for the amplifier cavity. Defaults to 30ppm.
    Ip: float
      Pump power. Defaults to 50 W.
    m: float
      Mass of the light mirror [kg]. Defaults to 100g.
    f0: float
      Pendulum resonant frequency [Hz]. Defaults to 0.6 Hz.
    thetas: array_like
      Angles of incidence, in degrees, on the cavity mirrors. Defaults to [45,0,45].
    phi: float
      Round-trip phase (radians) for the pump field. Defaults to 0 rad.
    lam: float
      Wavelength of light [m]. Defaults to 2 um.
    
    Methods:
    --------
    __init__(self, Ti=0.35e-2, Te=5e-6, Lrt=30e-6, Ip=50, m=100e-3, f0=1, thetas=[45,0,45], phi=0, lam=2e-6):
        Initializes a cavity. Also calculates amplitude transmissivities, reflectivities, cavity gain and sets these as attributes.
 
    cavGain(self, phi=0):
        Calculates the cavity field gain, for a round-trip phase of phi.
    ringdownTrans(self, phi=0, alpha=-1., t=1e-6*np.linspace(-50,200,500), noisy=False, sigma=5e-4):
        Calculates the ringdown in transmission monitored as power at MC2. Option to add some Gaussian noise to the time series.
    ringdownRefl(self, phi=0, alpha=-1., t=1e-6*np.linspace(-50,200,500), noisy=False, sigma=5e-4):
        Calculates the ringdown in transmission monitored as power at MC2. Option to add some Gaussian noise to the time series.     
    '''
    def __init__(self, Ti=0.35e-2, Te=5e-6, Lrt=30e-6, Ip=50, m=100e-3, f0=0.6, thetas=[45,0,45], phi=0, lam=2e-6):
        self.__dict__.update(locals()) # This hack assigns all the arguments passed to init.
        self.t1 = np.sqrt(self.Ti)
        self.t3 = np.sqrt(self.Te)
        self.t2 = np.sqrt(self.Te)
        self.r1 = np.sqrt(1 - self.t1**2 - self.Lrt/3)
        self.r2 = np.sqrt(1 - self.t2**2 - self.Lrt/3)
        self.r3 = np.sqrt(1 - self.t3**2 - self.Lrt/3)
        rho = 1. - self.t1**2 - self.t2**2 - self.t3**2 - self.Lrt
        self.finesse = np.pi/(2*np.arcsin((1-np.sqrt(rho))/(2*rho**0.25)))
        self.alphas = np.cos(np.deg2rad(self.thetas))**2
        num = self.t1
        den = 1. - np.exp(-1j*self.phi)*self.r1*self.r2*self.r3
        self.cavGain = num/den
        self.Icirc = self.Ip * np.abs(self.cavGain)**2
        self.omega0 = 2*np.pi*scc.c/self.lam
        return
    
    def chi(self, ff=np.logspace(1,4,505)):
      '''
      Compute the mechanical TF from force to displacement of the light amplifier mirror.
      
      Parameters:
      ===========
      ff: array_like
        Frequency vector on which to evaluate the transfer function
      '''
      Omega0 = 2*np.pi*self.f0
      Omega = 2*np.pi*ff
      self.chiTF = -1 / (self.m * (Omega**2 - Omega0**2))
      return(ff, self.chiTF)
    
    def Q(self, ff=np.logspace(1,4,505)):
      '''
      
      
      Parameters:
      ===========
      ff: array_like
        Frequency vector on which to evaluate the transfer function
      '''
      self.QTF = 4 * self.omega0 * np.sqrt(2*self.Icirc)/ scc.c / np.sqrt(self.Ti * scc.hbar * omega0)
      return(ff, QTF)
    
    def R(self, ff=np.logspace(1,4,505)):
      '''
      Calculate the amplifier optomechanical gain.
      
      Parameters:
      ===========
      ff: array_like
        Frequency vector on which to evaluate the transfer function
      '''
      _, _ = self.chi()
      self.RTF = np.abs((32 * self.omega0 * self.Icirc / (self.Ti *scc.c**2))*(self.chi(ff)[1]))
      return(ff, self.RTF)
    def freqDepRIN(self, ff=np.logspace(1,4,505), cornerFreq=100, cornerLevel=1e-9):
      '''
      Returns frequency dependent RIN ASD
      
      Parameters:
      ------------
      ff: array_like
        Frequency vector
      cornerFreq: float
        Above this, the RIN ASD is assumed flat. Rises as 1/f below.
      cornerLevel: float
        ASD at cornerFreq in 1/rtHz
    
      Returns:
      --------
      RIN: array_like
        Frequency dependent RIN
      '''
      ww, hh = sig.freqs_zpk([2*np.pi*cornerFreq],[0],cornerLevel, worN=2*np.pi*ff)
      self.RIN = np.abs(hh)
      return(ff, np.abs(hh))
    
    
    def printSummary(self):
      pprint.pprint(vars(self)) 
      return

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
