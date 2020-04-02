'''
Collection of useful Gaussian beam functions
'''
import numpy as np

def gaussBeamProfile(w0,z, lam=1064e-9):
    '''
    Gaussian beam width evolution along propagation direction.
    Assumes propagation in a vacuum.

    Parameters:
    ----------
    w0: float or array_like
        waist size [m]
    z: float or array_like
        Distance from waist position [m].
    lam: float
        Wavelength [m]

    Returns:
    --------
    w(z): float or array_like
        Gaussian beam width [m] as a function of the propagation
        coordinate, z.
    zR: float
        The Rayleigh range for the beam.
    '''
    zR = np.pi * w0**2 / lam
    w = w0*np.sqrt(1+(z/zR)**2)
    return(w, zR)

def gaussIntensity(P0, w0, r, z, lam=1064e-9):
    '''
    Gaussian beam intensity [W/m^2].
    
    Parameters:
    -----------
    P0: float
        Total power [W].
    w0: float
        Waist size (assumed at z=0).
    r: float or array_like
        Coordinate along beam cross section [m].
    z: float or array_like
        Coordinate along propagation direction [m].
    lam: float
        Wavelength [m].

    Returns:
    ---------
    I(r,z): float or array_like
        Intensity as a function of propagation coordinate and beam cross section.
    I0: float
        Peak intensity [W/m^2]
    '''
    # Get the peak intensity for the given total power
    I0 = 2*P0/np.pi/w0**2
    # First, calculate the beam width as a function of propagaiton coordinate.
    w, zR = gaussBeamProfile(w0, z, lam)
    I = I0 * (w0/w)**2 * np.exp(-2*r**2/w**2)
    return(I, I0)

def gaussPower(P0, w0, R, z, lam=1064e-9):
    '''
    Power through an aperture of radius R.
    It is assumed that the beam is centered on the aperture.

    Parameters:
    -----------
    P0: float
        Total power [W].
    w0: float
        Waist size (assumed at z=0).
    R: float or array_like
        Radius of aperture [m].
    z: float or array_like
        Coordinate along propagation direction [m].
    lam: float
        Wavelength [m].

    Returns:
    ---------
    P(R,z): float or array_like
        Power through aperture [W].
    '''
    w, zR = gaussBeamProfile(w0, z, lam)
    P = P0 * (1 - np.exp(-2*r**2/w**2))
    return(P)
