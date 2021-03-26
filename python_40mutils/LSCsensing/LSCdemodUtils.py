### Set of functions useful in doing complex demodulation
# Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.widgets as mw
from matplotlib.patches import Ellipse
import nds2
#import nbutils as nbu
import yaml
#import glob
import h5py
import scipy.signal as sig
import scipy.stats as scst
import timeit
from astropy.time import Time
import uncertainties as uc
from uncertainties import umath
import pickle
import datetime, os
import socket
import tqdm
import logging
import pandas as pd

logging.basicConfig(
    level=os.getenv('LOG_LEVEL', 'INFO'),
    format="%(levelname)s \n%(message)s")

# Directory setup
globFigDir = 'Figures/'
globDataDir = 'Data/'
if globFigDir.strip('/') not in os.listdir():
    #print('Figures directory not found, making it...')
    logging.debug('Figures directory not found, making it...')
    os.mkdir(globFigDir)
if globDataDir.strip('/') not in os.listdir():
    #print('Data directory not found, making it...')
    logging.debug('Data directory not found, making it...')
    os.mkdir(globDataDir)



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
    #if par['filename']+'.hdf5' in glob.glob(dataDir+'*hdf5'):
    if par['filename']+'.hdf5' in os.listdir(dataDir):
        print('Data file already exists...')
    else:
        DoFdict = par['DoFs']
        PDdict = par['PDs']
        print('Trying to download {} seconds of data from {}'.format(par['duration']
            ,par['tStart']))
        tic = timeit.default_timer()
        if 'pianosa' in socket.gethostname() or 'rossa' in socket.gethostname():
            server, port = 'fb', 8088
            print('On martian network, so using FB')
        else:
            server, port = par['ndsServer'], par['ndsPort']
            print('Not On martian network, so using {}:{} for NDS'.format(server, port))
        try:
            conn = nds2.connection(server, port)
        except:
            print('Cant open NDS connection ')
            return
        with h5py.File(dataDir+par['filename']+'.hdf5','w') as f:
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
    datFile = dataDir+par['filename']+'.hdf5'
    demodFile = dataDir+par['filename']+'_demod'+'.hdf5'
    # Open the HDF5 file
    f = h5py.File(datFile,'r')
    # Delete the demod file to avoid any conflicts
    if os.path.exists(demodFile):
        print('Demod file exists, deleting it...')
        os.remove(demodFile)
    f2 = h5py.File(demodFile,'w')
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
        # Also save to the HDF5 file so that we can avoid pickling
        f2.create_dataset(dof+'_demodOut', data=demod)
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
                # Also save to the HDF5 file so that we can avoid pickling
                f2.create_dataset(PD+'_'+dof+'_demodOut_I', data=demod_I)
                f2.create_dataset(PD+'_'+dof+'_demodOut_Q', data=demod_Q)
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
                # Also save to the HDF5 file so that we can avoid pickling
                f2.create_dataset(PD+'_'+dof+'_mag_I', data=mag_I)
                f2.create_dataset(PD+'_'+dof+'_mag_Q', data=mag_Q)
            else:
                print('Poor SNR of '+dof+' in '+PD+'! So Im not demodulating.')
    print('Demodulation complete - saving data to file...')
   # with open(par['filename']+'.p','wb') as ff:
   #     pickle.dump(PDresults,ff)
    print('Data saved... Time to plot...')
    f.close()
    return

def plotData(paramFile, saveFig=False):
    par = importParams(paramFile)
    DoFdict = par['DoFs']
    PDdict = par['PDs']
    demodFile = dataDir+par['filename']+'_demod'+'.hdf5'
    actTbl = []
    PDtbl = []
    for kk, vv in par['PDs'].items():
        try:
            PDtbl.append([kk, vv['gain'], vv['angle'], vv['Z'], vv['conv gain']])
        except KeyError as e:
            print(e)
            print('The parameter file does not specify all parameters required to convert measured TF to W/m.')
    for kk, vv in par['DoFs']:
        try:
            actTbl.append([kk, vv['act'], vv['mconv']])
        except KeyError as e:
            print(e)
            print('The parameter file does not specify all parameters required to convert measured TF to W/m.')
    physPDs = list(par['PDs'].keys())
    mag = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    phase = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    magU = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    phaseU = np.zeros([len(PDdict.keys()), len(DoFdict.keys())])
    fil = h5py.File(demodFile,'r')
    # Get the stats
    for kk, PD in enumerate(PDdict.keys()):
        for ii, dof in enumerate(DoFdict.keys()):
            conv = 10/2**15  # ADC Volts/count conversion
            conv /= 10**(PDdict[PD]['gain']/20.)  # PD whitening gain
            conv /= DoFdict[dof]['mconv']*DoFdict[dof]['freq']**-2  # Actuator meters/DAC count
            result = demodStats(fil[PD +'_'+dof+'_mag_I'][:] + 1j*fil[PD +'_'+dof+'_mag_Q'][:])
            mag[kk,ii] = result[0].nominal_value * conv
            phase[kk,ii] = result[1].nominal_value - np.deg2rad(PDdict[PD]['angle'])
            magU[kk,ii] = result[0].std_dev * conv
            phaseU[kk,ii] = result[1].std_dev
    # Make the plot
    rc('font', weight=900)
    thetaticks = np.arange(0,360,45)
    if len(PDdict) < 3:
        figs,ax = plt.subplots(1,3,subplot_kw=dict(projection='polar'),figsize=(16,12))
    else:
        figs,ax = plt.subplots(2,3,subplot_kw=dict(projection='polar'),figsize=(16,12))
    
    for iii, aa in enumerate(ax.flatten()):
        try:
            aa.set_title(physPDs[iii],fontsize=20,fontweight='bold',y=1.15)
            aa.plot([0,-np.deg2rad(PDdict[list(PDdict.keys())[iii]]['angle'])],[0,9],linewidth=6,linestyle='--',alpha=0.4,color='grey',label='PD I')
            aa.plot([0,np.pi/2-np.deg2rad(PDdict[list(PDdict.keys())[iii]]['angle'])],[0,9],linewidth=6,linestyle=':',alpha=0.4,color='grey',label='PD Q')
            aa.grid(linestyle='--', linewidth=0.7, alpha=0.9)
            aa.set_thetagrids(thetaticks)#,frac=1.2)
            for kkk, dof in enumerate(DoFdict.keys()):
                print('Plotting {} in {}'.format(dof, physPDs[iii]))
                ll = aa.plot([0,phase[iii,kkk]], [0,np.log10(mag[iii,kkk])], label=dof, marker='.', markersize=6)
                aa.add_patch(Ellipse(xy=(phase[iii][kkk],np.log10(mag[iii][kkk])),
                    width=10*np.abs(phaseU[iii][kkk]),
                    height=10*np.abs(np.log10(1+magU[iii][kkk]/mag[iii][kkk])),
                    angle=0.,alpha=0.4, color=ll[0].get_color()))
            if iii==0:
                handles, labels = aa.get_legend_handles_labels()
        except:
            pass
        
    figs.subplots_adjust(hspace=0.5, wspace=0.33)
    aa.axis('off') # Now aa is the last axis.
    # Add the tables.
    tbl1 = plt.table(cellText=PDtbl, fontsize=14, colLabels=['PD', 'Wht Gain [dB]','$\phi^{\circ}$','Z$\Omega$', r'\frac{V_{\mathrm{IF}}}{V_{\mathrm{RF}}}'],
         loc='upper left', colWidths=[0.1,0.2,0.2, 0.2,0.2], rowLoc='center', colLoc='center')
    tbl2 = plt.table(cellText=actTbl, fontsize=14, colLabels=['DoF', 'Actuator','DC gain [m/ct]'],
         loc='upper left', colWidths=[0.1,0.15,0.15], rowLoc='center', colLoc='center')
    legend = aa.legend(handles,labels,loc='best',fontsize=14, ncol=2)
    aa.text(2,1.5,'Radial axes are\n$\mathrm{log}_{10}(\mathrm{mag}).$\nUnits are [V/m].\nUncertainties multiplied by 10.',fontsize=14, weight='extra bold')
    if saveFig:
        figs.savefig(figDir+par['filename']+'sensMat.pdf', bbox_inches='tight')
        print('Sensing matrix pdf saved to {}'.format(figDir+par['filename']+'sensMat.pdf'))
    return()

def printMatrixElements(paramFile, sensorList, dofList):
    '''
    Function to print out the demodulated sensing matrix elements

    Parameters:
    -----------
    paramFile: str
        Path to the parameter file.
    sensorList: array_like
        List of sensors
    dofList: array_like
        List of DoFs
    '''
    par = importParams(paramFile)
    DoFdict = par['DoFs']
    PDdict = par['PDs']
    demodFile = dataDir+par['filename']+'_demod'+'.hdf5'
    fil = h5py.File(demodFile, 'r')
    for ss, dd in zip(sensorList, dofList):
        PD = ss.split(sep='_')[0]
        quad = ss.split(sep='_')[1]
        resp = np.mean(np.abs(fil[PD+'_'+dd+'_demodOut_'+quad]))
        resp /= DoFdict[dd]['mconv']*DoFdict[dd]['freq']**-2
        outStr = 'Response of {} in {}_{} is {:.2E} cts/m'.format(dd, PD, quad, resp)
        print(outStr)
    return()

def printMatrix(paramFile):
    '''
    Function to print out the demodulated sensing matrix.

    Parameters:
    -----------
    paramFile: str
        Path to the parameter file.

    Returns:
    --------
    matrix: pandas dataframe
        A pandas dataframe that is the sensing matrix
    '''
    def _round_expr(expr, num_digits):
        return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(sympy.Number)})
    par = importParams(paramFile)
    DoFdict = par['DoFs']
    PDdict = par['PDs']
    demodFile = dataDir+par['filename']+'_demod'+'.hdf5'
    matrix = np.ones((2*len(PDdict), len(DoFdict)))
    fil = h5py.File(demodFile, 'r')
    for jj, ss in enumerate(PDdict):
        for kk, dd in enumerate(DoFdict):
            resp_I = np.mean(np.abs(fil[ss+'_'+dd+'_demodOut_I']))
            resp_I /= DoFdict[dd]['mconv']*DoFdict[dd]['freq']**-2
            resp_Q = np.mean(np.abs(fil[ss+'_'+dd+'_demodOut_Q']))
            resp_Q /= DoFdict[dd]['mconv']*DoFdict[dd]['freq']**-2
            matrix[2*jj, kk] = resp_I
            matrix[2*jj+1, kk] = resp_Q
    cols = [aa+'{}'.format(bb) for aa in PDdict for bb in ['_I', '_Q']]
    print('Sensing matrix in cts/m is')
    df = pd.DataFrame(matrix.T, columns=cols, index=[aa for aa in DoFdict])
    print(df)
    return(df)

def sensingHistograms(paramFile, sensorList, dofList, saveFig=False):
    '''
    Function to make histograms of the demodulated outputs,
    for the pairings given in sensorList and dofList
    Parameters:
    -----------
    paramFile: str
        Path to the parameter file.
    sensorList: array_like
        List of sensors
    dofList: array_like
        List of DoFs
    '''
    par = importParams(paramFile)
    DoFdict = par['DoFs']
    PDdict = par['PDs']
    demodFile = dataDir+par['filename']+'_demod'+'.hdf5'
    # Make a figure. nRows = number of sensors, nCols = 2
    if len(dofList) <= 4:
        fig, ax = plt.subplots(2,2, figsize=(16,9))
    else:
        fig, ax = plt.subplots(3,2, figsize=(16,9))
    fil = h5py.File(demodFile, 'r')
    for jj, (ss, dd) in enumerate(zip(sensorList, dofList)):
        resp_I = np.abs(fil[ss+'_'+dd+'_demodOut_I']) / (DoFdict[dd]['mconv']*DoFdict[dd]['freq']**-2) / 1e9
        resp_Q = np.abs(fil[ss+'_'+dd+'_demodOut_Q']) / (DoFdict[dd]['mconv']*DoFdict[dd]['freq']**-2) / 1e9
        # Use scipy stats to estimate some probability density function
        density_I = scst.gaussian_kde(resp_I)
        density_Q = scst.gaussian_kde(resp_Q)
        ax.flatten()[jj].hist(resp_I, bins=20, alpha=0.6, label='I', density=True)
        ax.flatten()[jj].plot(sorted(resp_I), density_I(sorted(resp_I)), 
                linestyle='--', color='xkcd:charcoal')
        ax.flatten()[jj].hist(resp_Q, bins=20, alpha=0.6, label='Q', density=True)
        ax.flatten()[jj].plot(sorted(resp_Q), density_Q(sorted(resp_Q)), 
                
                linestyle='--', color='xkcd:charcoal')
        ax.flatten()[jj].set_xlabel('{} in {} [cts/nm]'.format(dd, ss))
        ax.flatten()[jj].set_ylabel('Normalized density')
    for aa in ax.flatten():
        aa.grid(True, which='both', linestyle='--', alpha=0.4)
    ax[0,0].legend(loc='best')
    ax[2,1].axis('off')
    fig.suptitle('LSC sensing responses')
    fig.subplots_adjust(top=0.96, hspace=0.25, wspace=0.25)
    if saveFig:
        fig.savefig(figDir+par['filename']+'sensMatHistograms.pdf', bbox_inches='tight')
        print('Sensing matrix pdf saved to {}'.format(figDir+par['filename']+'sensMatHistograms.pdf'))
    return()

