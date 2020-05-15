'''
Set of python utilities for the 40m interferometer
'''

import generalUtils as genUtils
import scipy.signal as sig
import numpy as np
import scipy.constants as scc

def AAfilt(ff, f0_z=1e3*np.array([13.5305051680, 15.5318357741, 23.1746351749]),
        Q_z=np.array([423.6130434049e6, 747.6895990654e3, 1.5412966100e6]),
        f0_p=1e3*np.array([5.2860544577, 5.9752193716, 8.9271953580, 8.2181747850, 182.1403534923]),
        Q_p=np.array([503.1473053928e-3, 1.0543411596, 3.5384364788, 3.4220607928, 1.1187869426]),
        k=989.1003181564e-3):
    '''
    Function that calculates the frequency response for the AA filters
    used at the 40m (D000076).
    
    Parameters:
    -----------
    ff: array_like
        Frequency vector on which to evaluate
    f0_z: array_like
        Array of poles [Hz]. Defaults to those measured in Elog 14939.
    Q_z: array_like
        Array of Q of poles. Defaults to those measured in Elog 14939.
    f0_p: array_like
        Array of zeros [Hz]. Defaults to those measured in Elog 14939.
    Q_p: array_like
        Array of Q of zeros. Defaults to those measured in Elog 14939.
    k: float
        Overall gain. Defaults to that measured in Elog 14939.

    Returns:
    ---------
    ww: array_like
        Vector of angular frequencies [rad/sec] on which response is evaluated.
    hh: array_like
        Response
    '''
    poles = genUtils.fq2reim(f0_p, Q_p).flatten()
    zeros = genUtils.fq2reim(f0_z, Q_z).flatten()
    k *= np.abs(np.prod(poles)/np.prod(zeros))
    ww, hh = sig.freqs_zpk(zeros, poles, k, 2*np.pi*ff)
    return(ww, hh)

def AIfilt(ff, f0_z=1e3*np.array([16.5683962778, 32.6135551736]),
        Q_z=np.array([1.0003811032e9, 85.8262574773]),
        f0_p=1e3*np.array([8.4550003709, 16.3154891526, 31.2574514935]),
        Q_p=np.array([1.5939230324, 2.3053074845, 2.2233195780]),
        k=933.1654699165e-3):
    '''
    Function that calculates the frequency response for the AI filters
    used at the 40m (D000186).
    
    Parameters:
    -----------
    ff: array_like
        Frequency vector on which to evaluate. Defaults to 801 logspaced 
        points between 0.1 Hz and 100 kHz.
    f0_z: array_like
        Array of poles [Hz]. Defaults to those measured in Elog 8591.
    Q_z: array_like
        Array of Q of poles. Defaults to those measured in Elog 8591.
    f0_p: array_like
        Array of zeros [Hz]. Defaults to those measured in Elog 8591.
    Q_p: array_like
        Array of Q of zeros. Defaults to those measured in Elog 8591.
    k: float
        Overall gain. Defaults to that measured in Elog 8591.

    Returns:
    ---------
    ww: array_like
        Vector of angular frequencies [rad/sec] on which response is evaluated.
    hh: array_like
        Response
    '''
    poles = genUtils.fq2reim(f0_p, Q_p).flatten()
    np.hstack((poles, 5.5846336114e3)) # Tack on the real pole
    zeros = genUtils.fq2reim(f0_z, Q_z).flatten()
    k *= np.abs(np.prod(poles)/np.prod(zeros))
    ww, hh = sig.freqs_zpk(zeros, poles, k, 2*np.pi*ff)
    return(ww, hh)

def DARMpole_RSE(Ts=9.903e-2, Ti=1.384e-2, Te=13.7e-6, Eta_arm=50e-6, Eta_SRM=1000e-6, Larm=37.795):
    '''
    Function that calculates the RSE pole for the 40m, under
    the single-pole approximation for cavity transfer functions.

    Parameters:
    ------------
    Ts: float
        SRM power transmissivity. Defaults to 9.903 %.
    Ti: float
        ITM power transmissivity, Defaults to 1.384 %.
    Te: float
        ETM power transmissivity. Defaults to 13.7 ppm.
    Eta_arm: float
        Arm cavity loss. Split equally between ITM and ETM.
        Defaults to 50 ppm, round trip.
    Eta_SRM: float
        SRC cavity scatter loss. Defaults to 1000ppm.
    Larm: float
        Length of arm cavity [m]. Defaults to 37.795 m.

    Returns:
    -----------
    fP: float
        DARM pole [Hz], assuming the SRC is configured for RSE.
    '''
    ts = np.sqrt(Ts)
    rs = np.sqrt(1 - Ts - Eta_SRM)
    ti = np.sqrt(Ti)
    ri = np.sqrt(1 - Ti - Eta_arm/2)
    te = np.sqrt(Te)
    re = np.sqrt(1 - Te - Eta_arm/2)
    num = 1 - ri*rs
    den = re*ri - re*rs*(ti**2 + ri**2)
    fP = scc.c / 2 / Larm 
    fP *= np.log(num/den)
    return(fP/2/np.pi)

def DARMpole_ESR(Ts=9.903e-2, Ti=1.384e-2, Te=13.7e-6, Eta_arm=50e-6, Eta_SRM=1000e-6, Larm=37.795):
    '''
    Function that calculates the ESR pole for the 40m, under
    the single-pole approximation for cavity transfer functions.

    Parameters:
    ------------
    Ts: float
        SRM power transmissivity. Defaults to 9.903 %.
    Ti: float
        ITM power transmissivity, Defaults to 1.384 %.
    Te: float
        ETM power transmissivity. Defaults to 13.7 ppm.
    Eta_arm: float
        Arm cavity loss. Split equally between ITM and ETM.
        Defaults to 50 ppm, round trip.
    Eta_SRM: float
        SRC cavity scatter loss. Defaults to 1000ppm.
    Larm: float
        Length of arm cavity [m]. Defaults to 37.795 m.

    Returns:
    -----------
    fP: float
        DARM pole [Hz], assuming the SRC is configured for ESR.
    '''
    ts = np.sqrt(Ts)
    rs = np.sqrt(1 - Ts - Eta_SRM)
    ti = np.sqrt(Ti)
    ri = np.sqrt(1 - Ti - Eta_arm/2)
    te = np.sqrt(Te)
    re = np.sqrt(1 - Te - Eta_arm/2)
    num = 1 + ri*rs
    den = re*ri + re*rs*(ti**2 + ri**2)
    fP = scc.c / 2 / Larm 
    fP *= np.log(num/den)
    return(fP/2/np.pi)


def DARMpole_PRFPMI(Ti=1.384e-2, Te=13.7e-6, Eta_arm=50e-6, Larm=37.795):
    '''
    Function that calculates the DARM pole for the 40m configured as a PRFPMI, 
    under the single-pole approximation for cavity transfer functions.

    Parameters:
    ------------
    Ti: float
        ITM power transmissivity, Defaults to 1.384 %.
    Te: float
        ETM power transmissivity. Defaults to 13.7 ppm.
    Eta_arm: float
        Arm cavity loss. Split equally between ITM and ETM.
        Defaults to 50 ppm, round trip.
    Larm: float
        Length of arm cavity [m]. Defaults to 37.795 m.

    Returns:
    -----------
    fP: float
        DARM pole [Hz], assuming a PRFPMI.
    '''
    ti = np.sqrt(Ti)
    ri = np.sqrt(1 - Ti - Eta_arm/2)
    te = np.sqrt(Te)
    re = np.sqrt(1 - Te - Eta_arm/2)
    num = 1
    den = re*ri
    fP = scc.c / 2 / Larm 
    fP *= np.log(num/den)
    return(fP/2/np.pi)
