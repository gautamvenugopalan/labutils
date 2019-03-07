'''
Collection of generally useful functions for making ELOG plots etc
'''
import sys, os, subprocess
import numpy as np
import scipy.signal as sig
import nds2
from matplotlib.ticker import FormatStrFormatter

def johnson(R,T=298):
    '''
    Calculate the Johnson noise for a resistance R at a temperature T.
    Parameters:
    ------------
    R: float or array_like
        Resistance in Ohms
    T: float or array_like
        Temperature in Kelvin. Defaults to 298.
    Returns:
    --------
    vj: float or array_like
        Voltage Johnson noise
    ij: float or array_like
        Current Johnson noise
    '''
    vj = np.sqrt(4*scc.k*T*R)
    ij = vj / R
    return(vj,ij)


def fq2reim(f0,q):
    '''
    Function that converts f0/q pole/zero to ReIm representation.
    Example usage:
        poles_cplx = fq2reim(f0,Q)
    will return an array of complex pole pairs (s-plane) 
    corresponding to the inputs f0 (in Hz) and Q.
    '''
    b = -1j * f0 / q
    c = -1. * f0**2
    pp = -(b/2) + 0.5*np.sqrt(b**2 - 4*c)
    pn = -(b/2) - 0.5*np.sqrt(b**2 - 4*c)
    return np.squeeze((-2 * np.pi / 1j)*np.array([pp,pn]))

def tarballz(outFile, fileList):
    '''
    Function to make a tarball of files in fileList.
    Example usage:
        tarballz('MC2_radPress.tgz',['file1.txt', 'file2.txt'])
    will make a tarball called MC2_radPress.tgz from file1.txt and file2.txt
    '''
    tarCommand = 'tar -cvzf '+outFile+' '
    for ff in fileList:
        tarCommand = tarCommand + ff + ' '
    FNULL = open(os.devnull, 'w')
    subprocess.call(tarCommand, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    return

def makeLSCspecgram(channel, tStart, tStop, t_fft=4, win=('tukey', 0.25), medianNorm=False,  ndsServer='131.215.115.200', ndsPort=31200):
    '''
    Computes a spectrogram by downloading data from NDS2.
    Example usage:
        ff, tt, Sxx = makeLSCspecgram(['C1:LSC-XARM_IN1_DQ'],1208808916,1208808917, t_fft=4, win=('tukey', 0.25))
    will download the data, and make a spectrogram with 4 second FFT time.
    Option to specify windowing function used by scipy.signal.spectrogram.
    Note that Power Spectral Density is returned, use sqrt where appropriate!
    '''
    conn = nds2.connection(ndsServer, ndsPort)
    print('Getting NDS data')
    dat = conn.fetch(tStart, tStop, channel)
    print('Got NDS data...')
    Fs = dat[0].channel.sample_rate
    ff, tt, Sxx = sig.spectrogram(dat[0].data,fs=Fs,nperseg=int(t_fft*Fs),window=win)
    Sxx = np.sqrt(Sxx)
    if medianNorm:
        Sxx = Sxx / np.median(Sxx, axis=-1)[:,None]
    return ff, tt, Sxx

def makeASD(channel, tStart, tStop, t_fft=4, win=('tukey', 0.25)):
    '''
    Computes an ASD by downloading data from NDS2 using pwelch.
    Example usage:
        ff, Pxx = makeASD(['C1:LSC-XARM_IN1_DQ'],1208808916,1208808917, t_fft=4, win=('tukey', 0.25))
    will download the data, and make a ASD with 4 second FFT time.
    Option to specify windowing function used by scipy.signal.welch
    Note that Power Spectral Density is returned, use sqrt where appropriate!
    '''
    conn = nds2.connection('nds40.ligo.caltech.edu', 31200)
    print('Getting NDS data')
    dat = conn.fetch(tStart, tStop, channel)
    print('Got NDS data...')
    Fs = dat[0].channel.sample_rate
    ff, Pxx = sig.welch(dat[0].data,fs=Fs,nperseg=int(t_fft*Fs),window=win)
    return ff, Pxx

def cum_rms(ff, asd): 
    df = np.roll(ff, -1) - ff
    df = df[:-1]
    cum_ms = np.zeros(len(ff))
    for ii in range(len(cum_ms)):
        cum_ms[ii] = np.trapz(asd[ii:]**2, ff[ii:])
    return cum_ms**0.5

def plotSpec(ff,ax,col=[1],scale=1,doFormat=False,xlabel='Frequency [Hz]',ylabel='$\mathrm{V_{rms}}/\sqrt{\mathrm{Hz}}$',**kwargs):
    '''
    Function for plotting a spectrum from SR785 or AG4395.
    Expects column 0 to be frequency, column 1/2 to be Vrms/rtHz.
    Plotting arguments, e.g. color, are passed as **kwargs.
    First argument is either a file name (to a file downloaded from SR785/AG4395)
    or an array of 2 or 3 columns, with the first being frequency.
    Option to scale all y-axis values by some fixed number
    '''
    try:
        dat = np.loadtxt(ff)
    except:
        dat=ff
    
    for ii in col:
        ax.loglog(dat[:,0],dat[:,int(ii)]*scale,**kwargs)
    if doFormat:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid('on', which='both', linestyle='--',alpha=0.4)
    return

def plotTF(ff,axMag, axPh, cols='dbdeg',doFormat=False,magLabel='Magnitude [dB]',**kwargs):
    '''
    Function for plotting a transfer function from SR785 or AG4395 or modeling.
    Input has to be a 3 column file, with first column frequency.
    Plotting arguments, e.g. color, are passed as **kwargs
    First argument is either a file name (to a file downloaded from SR785/AG4395)
    or an array of 3 columns, [freq,mag,phase].
    '''
    try:
        dat = np.loadtxt(ff)
    except:
        dat = ff

    if cols=='dbdeg':
        axMag.semilogx(dat[:,0],dat[:,1],**kwargs)
        axPh.semilogx(dat[:,0], dat[:,2],**kwargs)
    elif cols=='absdeg':
        axMag.loglog(dat[:,0],dat[:,1],**kwargs)
    elif cols=='reim':
        cplx = dat[:,1] + 1j*dat[:,2]
        axMag.semilogx(dat[:,0],20*np.log10(np.abs(cplx)),**kwargs)
        axPh.semilogx(dat[:,0],np.rad2deg(np.angle(cplx)))
    else:
        print('Unrecognized file formatting. Please use either dbdeg, absdeg or reim')
        return
    if doFormat:
        axMag.set_ylabel(magLabel)
        axPh.set_ylabel('Phase [degrees]')
        axPh.set_xlabel('Frequency [Hz]')
        axPh.set_yticks(np.linspace(-180,180,9))
        axPh.yaxis.set_major_formatter(FormatStrFormatter("%2d"))
        axMag.grid(True,which='both',linestyle='--',alpha=0.4)
        axPh.grid(True,which='both',linestyle='--',alpha=0.4)
    return

def getTopNoisesLISO(lisoFile, ff=100, nNoise=3):
    '''
    Function to identify the top nNoise noise contributions from a given LISO output file.
    Returns the indices of these noise columns in the LISO file.
    Example usage:
        ind = getTopNoisesLISO('noise.out',100,3)
    returns a 3-tuple with the top 3 noises at 100Hz in noise.out
    '''
    dat = np.loadtxt(lisoFile)
    if nNoise > np.shape(dat)[1]:
        print('To few noises in this file to identify the number you requested.')
        print(getTopNoisesLISO.__doc__)
        return
    xInd = np.argmin(np.abs(dat[:,0]-ff))
    arr = dat[xInd,:]
    return np.flip(np.argsort(arr[1:-1])+1,axis=0)[:nNoise]
