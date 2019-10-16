### Set of functions useful in doing complex demodulation
# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.widgets as mw
from matplotlib.patches import Ellipse
import nds2
import nbutils as nbu
import yaml
import glob
import h5py
import scipy.signal as sig
import timeit
from astropy.time import Time
import uncertainties as uc
from uncertainties import umath
import pickle
import datetime, os
import socket
import tqdm

# Read in the YAML parameter file
def importParams(paramFile):
    '''
    Function to read in a parameter file with list
    of LSC photodiodes and actuator drives/frequencies.
    ------------
    Returns:
        params: A dictionary with the names and parameters
    '''
    with open(paramFile,'r') as f:
        params = yaml.safe_load(f)
    return params

def dlData(paramFile):
    '''
    Executes downloading of data from an nds server,
    for paramters specified in the paramFile.
    -------------
    Returns:
        Status of download
    '''
    par = importParams(paramFile)
    if par['filename']+'.hdf5' in glob.glob('*hdf5'):
        print('Data file already exists...')
    else:
        DoFdict = par['DoFs']
        PDdict = par['PDs']
        print('Trying to download {} seconds of data from {}'.format(par['duration']
            ,par['tStart']))
        tic = timeit.default_timer()
        try:
            conn = nds2.connection(par['ndsServer'],par['ndsPort'])
        except:
            print('Cant open NDS connection ')
            return
        with h5py.File(par['filename']+'.hdf5','w') as f:
            # Download the LSC photodiode data
            for PD, PDd in tqdm.tqdm(PDdict.items()):
                Ichan = 'C1:LSC-'+PD+'_I_ERR_DQ'
                Qchan = 'C1:LSC-'+PD+'_Q_ERR_DQ'
                dat = conn.fetch(par['tStart'],par['tStart']+par['duration'],
                                 [Ichan, Qchan])
                IQdata = np.vstack((dat[0].data,dat[1].data)).T # Column #1 is I, Column #2 is Q
                PDdata = f.create_dataset(PD, data=IQdata)
                for k, v in PDd.items():
                    PDdata.attrs[k] = v
            # Download the suspension output signal data
            for dof, dofd in DoFdict.items():
                optic = dofd['act']
                actChan = 'C1:SUS-'+optic+'_LSC_OUT_DQ'
                actData = conn.fetch(par['tStart'],par['tStart']+par['duration'],[actChan]);
                Adata = f.create_dataset(optic, data=actData[0].data)
                for k, v in dofd.items():
                    Adata.attrs[k] = v
            print('Download successful, closing file...')
            f.close()
            toc = timeit.default_timer() - tic
            print('Elapsed time is {} seconds'.format(round(toc,2)))
    return

def fft(dat, fDemod, fs=2**14, tFFT=5, win=('tukey',0.25),
        nOverlap=0, detrend='constant', median=False):
    '''
    Single frequency digital demodulation.

    Parameters:
    -----------
    dat: array_like
         Data to demodulate. May contain multiple time series
    fDemod: float
         Frequency at which to do the demodulation.
    fs: float, optional
        Sampling frequency of the time series. Defaults to 2**14=16384 Hz.
    tFFT: float, optional
        Segment length (seconds) to evaluate the FFT. Defaults to 5 s.
    win: tuple, optional
        Input to scipy.signal window function. Defaults to Tukey window with alpha=0.25
    nOverlap: int, optional
        Number of samples to overlap window. Defaults to 0.
    detrend: string, optional
        Input to scipy.signal detrend function. Defaults to 'constant'
    median: Bool, optional
        Median averaging of final result. Defaults to False.
        
    Returns:
    --------
    result: complex
        Result of the digital demodulation.
    TODO: error handling...
    '''
    
    dat = np.asarray(dat)
    if dat.size==0:
        return(np.empty(dat.shape[-1]))
    nperseg = int(np.round(tFFT*fs)) # Number of segments in a sample segment.
    nOverlap = int(nOverlap);
    # Make the LO time series
    tt = np.arange(len(dat))/fs
    LO = np.exp(-1j*2*np.pi*fDemod*tt)
    
    # Compute the step to take as we stride through the segments
    step = nperseg - nOverlap
    segShape = ((dat.shape[-1]-nOverlap)//step, nperseg)
    datStrides = (step*dat.strides[-1], dat.strides[-1])
    LOStrides = (step*LO.strides[-1], LO.strides[-1])
    dS = np.lib.stride_tricks.as_strided(dat, shape=segShape, strides=datStrides)
    LOS = np.lib.stride_tricks.as_strided(LO, shape=segShape, strides=LOStrides)
    
    # Detrend the data
    data = sig.detrend(dS, type=detrend)
    
    # Demodulate the (windowed) data
    wind = sig.get_window(win,nperseg)
    result = data * wind * LOS
    result = result.sum(axis=-1)
    scale = 2*np.sqrt(1/(wind.sum()**2))
    result *= scale
    return result

def demodStats(dat):
    '''
    Computes some uncertainties on the
    demodulated complex amplitude
    
    Parameters:
    -----------
        dat: array like
            Array of demodulated complex amplitudes
    Returns:
    ----------
        tot_unc: 2-tuple
            Mean and uncertainties in the magnitude and phase [rad]
    '''
    data = np.copy(dat)
    N = len(data)
    stats_real = uc.ufloat(np.mean(data.real), np.std(data.real)/np.sqrt(N))
    stats_imag = uc.ufloat(np.mean(data.imag), np.std(data.imag)/np.sqrt(N))
    tot_unc = [(stats_real**2 + stats_imag**2)**0.5, umath.atan2(stats_imag,stats_real)]
    return tot_unc
  
def demodData(paramFile):
    par = importParams(paramFile)
    datFile = par['filename']+'.hdf5'
    # Open the HDF5 file
    f = h5py.File(datFile,'r')
    DoFdict = par['DoFs']
    PDdict = par['PDs']
    PDresults = {}
    fs = 2**14
    tFFT = 5
    # First, look at the actuator signals
    for dof, dofd in DoFdict.items():
        dofd['actData'] = f[dofd['act']][:]
        demod = fft(dofd['actData'], dofd['freq'], fs=fs, tFFT=tFFT)
        # Do an SNR check
        dofd['demodOut'] = demodStats(demod)
        if dofd['demodOut'][0].nominal_value / dofd['demodOut'][0].std_dev < 10:
            print('Poor SNR in '+ dof + 'actuator signal. Check demod frequency maybe?')
    # Repeat for the Photodiodes...
    for PD, PDd in PDdict.items():
        PDresults[PD] = dict()
        dat = f[PD][:] # This is the stacked I and Q data
        print('Analyzing '+PD+'...')
        for dof, dofd in DoFdict.items():
            # Data quality check, using coherence b/w suspension and PD as a measure...
            _, coh = sig.coherence(dat.T, dofd['actData'],fs=fs, nperseg=fs)
            fInd = np.round(dofd['freq'])
            #if any(coh[:,fInd] > 0.5):
            if True:
                # Demodulate the PD outputs
                demod_I = fft(dat[:,0], dofd['freq'], fs=fs, tFFT=tFFT)
                demod_Q = fft(dat[:,1], dofd['freq'], fs=fs, tFFT=tFFT)
                # Convert the cts/rtHz peak height into cts/m
                demod_actuator = dofd['demodOut']
                demod_I /= demod_actuator[0].nominal_value
                demod_Q /= demod_actuator[0].nominal_value
                # Make it a complex number
                demod_I *= np.exp(-1j*demod_actuator[1].nominal_value)
                demod_Q *= np.exp(-1j*demod_actuator[1].nominal_value)
                mag_I = np.abs(demod_I)*np.sign(np.mean(demod_I.real))
                mag_Q = np.abs(demod_Q)*np.sign(np.mean(demod_Q.real))
                PDresults[PD][dof] = demodStats(mag_I + 1j*mag_Q)
            else:
                print('Poor SNR of '+dof+' in '+PD+'! So Im not demodulating.')
    print('Demodulation complete - saving data to file...')
    f.close()
    with open(par['filename']+'.p','wb') as ff:
        pickle.dump(PDresults,ff)
    print('Data saved... Time to plot...')
    return 

def plotData(paramFile, saveFig=False):
    par = importParams(paramFile)
    demodFile = par['filename']+'.p'
    with open(demodFile,'rb') as f:
        PDresults = pickle.load(f)
    DoFdict = par['DoFs']
    PDdict = par['PDs']
    #Pull the info from the demodulation and package it into arrays.
    n_PDs = len(PDdict.keys())
    physPDs = list(par['PDs'].keys())
    mag = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    phase = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    magU = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    phaseU = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    for kk, PD in enumerate(PDdict.keys()):
        result = PDresults[PD]
        for ii, dof in enumerate(DoFdict.keys()):
            conv = 10/2**15  # ADC Volts/count conversion
            conv /= 10**(PDdict[PD]['gain']/20.)  # PD whitening gain
            conv /= DoFdict[dof]['mconv']*DoFdict[dof]['freq']**-2  # Actuator meters/DAC count
            mag[kk,ii] = result[dof][0].nominal_value * conv
            phase[kk,ii] = result[dof][1].nominal_value - np.deg2rad(PDdict[PD]['angle'])
            magU[kk,ii] = result[dof][0].std_dev * conv
            phaseU[kk,ii] = result[dof][1].std_dev
    #Make the plot
    plt.style.use('default')
    rc('font', weight=900)
    thetaticks = np.arange(0,360,45)
    figs,ax = plt.subplots(2,3,subplot_kw=dict(projection='polar'),figsize=(16,12))

    for iii in range(n_PDs):
        if iii<3:
            ax[0,iii].plot([0,phase[iii,2]],[0,np.log10(mag[iii,2])], label='SRCL',marker='.', markersize=6)
            ax[0,iii].plot([0,phase[iii,0]],[0,np.log10(mag[iii,0])], label='MICH',marker='.', markersize=6)
            ax[0,iii].plot([0,phase[iii,1]],[0,np.log10(mag[iii,1])], label='PRCL',marker='.', markersize=6)
            #Add the uncertainty ellipses
            ax[0,iii].add_patch(Ellipse(xy=(phase[iii][2],np.log10(mag[iii][2])),
                    width=10*np.abs(phaseU[iii][2]),
                    height=10*np.abs(np.log10(1+magU[iii][2]/mag[iii][2])),
                    angle=0.,alpha=0.4, color='#1f77b4'))
            ax[0,iii].add_patch(Ellipse(xy=(phase[iii][0],np.log10(mag[iii][0])),
                    width=10*np.abs(phaseU[iii][0]),
                    height=10*np.abs(np.log10(1+magU[iii][0]/mag[iii][0])),
                    angle=0.,alpha=0.4, color='#ff7f0e'))
            ax[0,iii].add_patch(Ellipse(xy=(phase[iii][1],np.log10(mag[iii][1])),
                    width=10*np.abs(phaseU[iii][1]),
                    height=10*np.abs(np.log10(1+magU[iii][1]/mag[iii][1])),
                    angle=0.,alpha=0.4, color='#2ca02c'))
            ax[0,iii].set_title(physPDs[iii],fontsize=20,fontweight='bold',y=1.15)
            ax[0,iii].tick_params(labelsize=16)
            ax[0,iii].set_thetagrids(thetaticks)#,frac=1.2)
            ax[0,iii].grid(linestyle='--', linewidth=0.7, alpha = 0.9)
            #Add a line indicating the I and Q phases for this photodiode
            ax[0,iii].plot([0,-np.deg2rad(PDdict[list(PDdict.keys())[iii]]['angle'])],[0,9],linewidth=6,linestyle='--',alpha=0.4,color='grey',label='PD_I')
            ax[0,iii].plot([0,np.pi/2-np.deg2rad(PDdict[list(PDdict.keys())[iii]]['angle'])],[0,9],linewidth=6,linestyle=':',alpha=0.4,color='grey',label='PD_Q')
        else:
            ax[1,iii-3].plot([0,phase[iii,2]],[0,np.log10(mag[iii,2])],label='SRCL',marker='.', markersize=6)
            ax[1,iii-3].plot([0,phase[iii,0]],[0,np.log10(mag[iii,0])], label='MICH',marker='.', markersize=6)
            ax[1,iii-3].plot([0,phase[iii,1]],[0,np.log10(mag[iii,1])], label='PRCL',marker='.', markersize=6)
            ax[1,iii-3].set_title(physPDs[iii],fontsize=20, fontweight='bold',y=1.15)
            #ax[1,iii-3].set_yticklabels([])
            ax[1,iii-3].tick_params(labelsize=16)
            ax[1,iii-3].set_thetagrids(thetaticks)#,frac=1.2)
            #ax[1,iii-3].set_yticks(range(1,10))
            ax[1,iii-3].grid(linestyle='--', linewidth=0.7, alpha = 0.9)
            ax[1,iii-3].plot([0,-np.deg2rad(PDdict[list(PDdict.keys())[iii]]['angle'])],[0,9],linewidth=6,linestyle='--',alpha=0.4,color='grey',label='PD_I')
            ax[1,iii-3].plot([0,np.pi/2-np.deg2rad(PDdict[list(PDdict.keys())[iii]]['angle'])],[0,9],linewidth=6,linestyle=':',alpha=0.4,color='grey',label='PD_Q')
            #Add the uncertainty ellipses
            ax[1,iii-3].add_patch(Ellipse(xy=(phase[iii][2],np.log10(mag[iii][2])),
                    width=10*np.abs(phaseU[iii][2]),
                    height=10*np.abs(np.log10(1+magU[iii][2]/mag[iii][2])),
                    angle=0.,alpha=0.4, color='#1f77b4'))
            ax[1,iii-3].add_patch(Ellipse(xy=(phase[iii][0],np.log10(mag[iii][0])),
                    width=10*np.abs(phaseU[iii][0]),
                    height=10*np.abs(np.log10(1+magU[iii][0]/mag[iii][0])),
                    angle=0.,alpha=0.4, color='#ff7f0e'))
            ax[1,iii-3].add_patch(Ellipse(xy=(phase[iii][1],np.log10(mag[iii][1])),
                    width=10*np.abs(phaseU[iii][1]),
                    height=10*np.abs(np.log10(1+magU[iii][1]/mag[iii][1])),
                    angle=0.,alpha=0.4, color='#2ca02c'))

    ax[1,2].axis('off')
    figs.subplots_adjust(hspace=0.5, wspace=0.33)
    handles, labels = ax[1,1].get_legend_handles_labels()
    legend = ax[1,2].legend(handles,labels,loc='center',fontsize=20)
    ax[1,2].text(2,1.5,'Radial axes are\n$\mathrm{log}_{10}(\mathrm{mag}).$\nUnits are [V/m].\nUncertainties multiplied by 10.',fontsize=14, weight='extra bold')

    #Print the nominal values of the sensing matrix onto the plot.
    mich = mag[3][0] * np.abs(np.cos(phase[4][0] + np.deg2rad(90) + np.deg2rad(PDdict['AS55']['angle']))) #MICH in AS55_Q
    prcl = mag[0][1] * np.abs(np.cos(phase[0][1])) #PRCL in REFL11_I
    srcl = mag[4][2] * np.abs(np.cos(phase[1][2])) #SRCL in REFL55_I
    ax[1,2].text(2.4,0.9,'$\mathrm{MICH}_{\mathrm{AS55Q}} = %.3E \mathrm{ V/m} $\n$\mathrm{PRCL}_{\mathrm{REFL11I}} = %.3E \mathrm{V/m}$\n$\mathrm{SRCL}_{\mathrm{REFL55I}} = %.3E  \mathrm{V/m}$\n'% (mich,prcl,srcl),fontsize=12, weight='extra bold')
    if saveFig:
        figs.savefig(figDir+'sensMat.pdf', bbox_inches='tight')
        print('Sensing matrix pdf saved to {}'.format(figDir+'sensMat.pdf'))
    return
