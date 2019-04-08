# Some functions to facilitate commonly used unit conversions that are not readily available in Numpy

def NF(vn, R=50, T=298):
    '''
    RF amplifier noise floor, only considering Johnson noise of 50 ohms.

    Parameters:
    -----------
    vn: float or array_like
        Voltage noise level
    R: float
        Resistance whose Johnson noise is the "Noise" in the SNR calculation. Defaults to 50 ohms.
    T: float
        Temperature at which Johnson noise is calculated. Defaults to 298 K.
    Returns:
    --------
    noiseFig: float or array_like
        Noise figure of amplifier
    '''
    num = 4 * scc.k * T * R/2 + vn**2
    den = 4 * scc.k * T * R/2
    return(10*np.log10(num/den))
