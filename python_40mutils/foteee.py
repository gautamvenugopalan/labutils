'''
Set of python utilities for the 40m interferometer
'''

import generalUtils as genUtils
import scipy.signal as sig
import numpy as np
import scipy.constants as scc
import readFotonFilterFile

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

def CARMpole(Ti=1.384e-2, Te=13.7e-6, Tp=5.637e-2, Eta_arm=50e-6, Eta_PRC=0.5e-2, Larm=37.795):
    '''
    Function that calculates the CARM pole for the 40m.

    Parameters:
    ------------
    Ti: float
        ITM power transmissivity, Defaults to 1.384 %.
    Te: float
        ETM power transmissivity. Defaults to 13.7 ppm.
    Tp: float
        PRM power transmissivity. Defaults to 5.637 %.
    Eta_arm: float
        Arm cavity loss. Split equally between ITM and ETM.
        Defaults to 50 ppm, round trip.
    Eta_PRC: float
        Loss in the power recycling cavity, including mode
        mismatch. Defaults to 0.5%.
    Larm: float
        Length of arm cavity [m]. Defaults to 37.795 m.

    Returns:
    --------
    fP: float
        CARM pole [Hz].
    '''
    tp = np.sqrt(Tp)
    rp = np.sqrt(1 - Tp - Eta_PRC)
    ti = np.sqrt(Ti)
    ri = np.sqrt(1 - Ti - Eta_arm/2)
    te = np.sqrt(Te)
    re = np.sqrt(1 - Te - Eta_arm/2)
    num = 1 + ri*rp
    den = re*ri + re*rp*(ti**2 + ri**2)
    fP = np.log(num/den) * scc.c / 2 / Larm / 2 / np.pi
    return(fP)

def DAIfilt(ff, sosCoeffs=np.array([[1,-1.41346289716898, 0.99893884152400, 
    1, -1.71662585474518, 0.78495484219691],
    [1, 0.00000127375260, 0.99819981588176, 1, -1.68385964238855, 0.93734519457266]]), 
    overallGain=0.014805052402446, fs=2*np.pi*2**16):
    '''
    Calculate the frequency response of the CDS Digital anti-aliasing filter.

    Parameters:
    ------------
    ff: array_like
        Frequency vector on which to evaluate
    sosCoeffs: array_like
        Coefficients describing the digital AI filter.
        Defaults to what is found in 40m elog 3961.
    overallGain: float
        Overall gain of the digital AI filter.
        Defaults to what is found in 40m elog 3961.
    fs: int
        Sampling frequency [Hz]. Defaults to 2**16 Hz.

    Returns:
    --------
    ww: array_like
        Vector of angular frequencies [rad/sec] on which response is evaluated.
    hh: array_like
        Response
    '''
    ww, hh = sig.sosfreqz(sosCoeffs, worN=2*np.pi*ff, whole=True, fs=fs)
    return(ww,hh*overallGain)

def DACTF(ff, fs=2**16):
    '''
    DAC transfer function in going from discrete time to
    a voltage when sampled at fs.

    Parameters:
    -----------
    ff: array_like
        Frequency vector on which to evaluate the response
    fs: int
        Sampling frequency [Hz]. Defaults to 64 kHz.

    Returns:
    ---------
    D: array_like
        The transfer function from cts to V.
    '''
    x = np.pi*ff/fs
    D = np.exp(-1j*x)
    D *= np.sin(x)/x
    return(D)

def CDSfilts(ff):
    '''
    Convenience function to evaluate the combination of
    AA --> DAA --> DAI ---> AI --> DAC.
    This is mainly to avoid repeated evaluations of these 
    functions, which are typically fixed.

    Parameters:
    -----------
    ff: array_like
        Frequency vector on which to evaluate the TF
    
    Returns:
    ----------
    Hdigital: array_like
        The transfer function of the combo indicated.
    '''
    # Construct the overall loop
    _, Hdigital = AAfilt(ff)                   # ADC AA
    Hdigital *= DAIfilt(ff)[1]                 # Digital AA
    Hdigital *= DAIfilt(ff)[1]                 # Digital AI
    Hdigital *= DACTF(ff, fs=2**16)            # DAC transfer function
    Hdigital *= AIfilt(ff)[1]
    return(Hdigital)

def DARMloop(ff, filtFile, filts, optGain=1e8, actGainDC=1e-9, delay=200e-6, returnDict=True):
    '''
    Calculate the DARM loop.
    Primarily to invert the noise curve for noise budgeting.

    Parameters:
    ------------
    ff: array_like
        Frequency vector [Hz]
    filtFile: str
        Path to the foton filter file to extract the CDS 
        DARM loop.
    filts: array_like
        Indices of filters turned on in the DARM filter bank.
    optGain: float
        DARM optical gain [W/m]. Defaults to 0.1 GW/m.
    actGainDC: float
        DC scaling [nm/ct] for the DARM pendulum transfer function.
        Defaults to 1 nm/ct.
    delay: float
        Delay [s] in the DARM loop. Defaults to 200us.
    returnDict: bool
        Optionally, return a breakdown of the various pieces
        composing the DARM loop. This defaults to True.

    Return:
    --------
    H: array_like
        DARM transfer function
    '''
    # Pendulum TF
    actTF = actGainDC / ff**2
    # Load the digital filter 
    filtDict = readFotonFilterFile.readFilterFile(filtFile)
    digitalFilt = np.array([1,0,0,1,0,0])
    for ii in filts:
        digitalFilt = np.vstack((digitalFilt, filtDict['LSC_DARM'][ii]['sosCoeffs']))
    _, digitalFiltResp = sig.sosfreqz(digitalFilt, worN=2*np.pi*ff, whole=True, fs=2*np.pi*filtDict['fs'])
    # DARM pole
    darmPole = DARMpole_PRFPMI()
    darmPoleTF = sig.freqs_zpk([], [-2*np.pi*darmPole], 2*np.pi*darmPole, worN=2*np.pi*ff)[1]
    # Construct the overall loop
    H = CDSfilts(ff)                    # ADC AA + DAA + DAI + DAC + AI
    H *= digitalFiltResp                # Digital loop
    H *= actTF                          # Pendulum TF
    H *= optGain                        # Optical gain
    H *= darmPoleTF                     # DARM pole
    
    if returnDict:
        # Build up a dictionary of the breakdown
        darmDict = {}
        darmDict['AA'] = AAfilt(ff)[1]
        darmDict['DAA'] = DAIfilt(ff)[1]
        darmDict['CDSfilt'] = digitalFiltResp
        darmDict['DAI'] = DAIfilt(ff)[1]
        darmDict['DAC'] = DACTF(ff)
        darmDict['AI'] = AIfilt(ff)[1]
        darmDict['actuatorTF'] = actTF
        darmDict['optGain'] = np.ones_like(ff) * optGain * darmPoleTF
        darmDict['Total'] = H
        return(H, darmDict)
    else:
        return(H)

