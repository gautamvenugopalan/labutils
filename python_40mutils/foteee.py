'''
Set of python utilities for the 40m interferometer
'''

import generalUtils as genUtils
import scipy.signal as sig
import numpy as np

def AAfilt(f0_z=1e3*np.array([13.5305051680, 15.5318357741, 23.1746351749]),
        Q_z=np.array([423.6130434049e6, 747.6895990654e3, 1.5412966100e6]),
        f0_p=1e3*np.array([5.2860544577, 5.9752193716, 8.9271953580, 8.2181747850, 182.1403534923]),
        Q_p=np.array([503.1473053928e-3, 1.0543411596, 3.5384364788, 3.4220607928, 1.1187869426]),
        k=989.1003181564e-3, ff=np.logspace(1,5,801)):
    '''
    Function that calculates the frequency response for the AA filters
    used at the 40m (D000076).
    
    Parameters:
    -----------
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
    ff: array_like
        Frequency vector on which to evaluate. Defaults to 801 logspaced 
        points between 0.1 Hz and 100 kHz.

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

def AIfilt(f0_z=1e3*np.array([16.5683962778, 32.6135551736]),
        Q_z=np.array([1.0003811032e9, 85.8262574773]),
        f0_p=1e3*np.array([8.4550003709, 16.3154891526, 31.2574514935]),
        Q_p=np.array([1.5939230324, 2.3053074845, 2.2233195780]),
        k=933.1654699165e-3, ff=np.logspace(1,5,801)):
    '''
    Function that calculates the frequency response for the AI filters
    used at the 40m (D000186).
    
    Parameters:
    -----------
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
    ff: array_like
        Frequency vector on which to evaluate. Defaults to 801 logspaced 
        points between 0.1 Hz and 100 kHz.

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
