#from gwpy.timeseries import TimeSeries, TimeSeriesDict
from matplotlib import pyplot as plt
import numpy as np
import scipy.signal as sig
import os
import nds2
import scipy.io as scio
import noisesub as n

nds_server = 'nds40.ligo.caltech.edu'
nds_port = 31200

def esplit(estring):
    '''
    Split an e-formatted string (e.g., 7.63e-08) and separately return the
    significand and the exponent.

    '''
    significand, exponent = estring.split('e')
    if exponent[0] == '+':
        sign = '+'
        exponent = exponent[1:]
    elif exponent[0] == '-':
        sign = '-'
        exponent = exponent[1:]
    while exponent[0] == '0':
        exponent = exponent[1:]
    if sign == '-':
        exponent = '-' + exponent
    return significand, exponent

def efmt(number):
    '''
    Take a significand and exponent (as strings) and return the properly
    LaTeX-formatted string for the overall number in scientific notation.
    '''
    significand, exponent = esplit('{:.3e}'.format(number))
    return '{}\times10^{{ {} }}'.format(significand, exponent)

def logbin(arr, N):
    idx = np.logspace(0, np.log10(len(arr)), N)
    idx = np.unique(np.ceil(idx))
    idx = np.concatenate([[0], idx])
    arr2 = np.zeros(len(idx)-1)
    for ii in range(len(idx)-1):
        arr2[ii] = np.sqrt(np.mean(arr[idx[ii]:idx[ii+1]]**2))
    return arr2

def set_grid(*axs):
    if len(axs) < 1:
        axs = [plt.gca()]
    for ax in axs:
        ax.grid('on', which='both', alpha=0.2)
        ax.grid(which='major', linestyle='solid', alpha=0.8)

def get_complex_interp(fname, ff, **kwargs):
    data = scio.loadmat(fname)
    f1 = data['f'][:,0]
    mag = data['mag'][:,0]
    phase = data['phase'][:,0]
    z = mag * np.exp(1j*np.deg2rad(phase))
    re = np.real(z)
    im = np.imag(z)
    #f1, re, im = np.loadtxt(fname, unpack=1)
    re2 = np.interp(ff, f1, re, **kwargs)
    im2 = np.interp(ff, f1, im, **kwargs)
    return re2+1j*im2

def bode_plot(ff, tfObj):
    h = plt.figure()
    axMag = h.add_subplot(211)
    axPha = h.add_subplot(212, sharex=axMag)
    if isinstance(tfObj, list):
        for tf in tfObj:
            axMag.loglog(ff, np.abs(tf))
            axPha.semilogx(ff, np.angle(tf, deg=True))
    else:
        axMag.loglog(ff, np.abs(tfObj))
        axPha.semilogx(ff, np.angle(tfObj, deg=True))
    myyticks = np.arange(-225, 226, 45)
    myyticklabels = ['$'+str(tick)+'^\circ$' for tick in myyticks]
    axPha.set_yticks(myyticks)
    axPha.set_yticklabels(myyticklabels)
    axPha.set_ylim(-190, 190)
    set_grid(axMag)
    set_grid(axPha)
    axMag.set_ylabel('Magnitude')
    axPha.set_ylabel('Phase')
    axPha.set_xlabel('Frequency [Hz]')
    return h, axMag, axPha

def get_med_asds(TSobj, fftlength=10, overlap=5, method='welch', window='hann',
               binNum=None, zpk=None, fullOutput=False, **kwargs):
    '''
    Takes a GWpy timeseries object and returns 10, 50, and 90 percentile ASDs.
    Optional logarithmic rebinning.

    '''
    myspec = (TSobj.spectrogram(fftlength, fftlength=fftlength, overlap=overlap, window=window))**0.5
    if zpk != None:
        myspec = myspec.zpk(*zpk)
    ff = np.arange(0, len(myspec.value[0])) * myspec.df.value + myspec.f0.value
    myspecdata = np.copy(myspec.value)
    myspecdata2 = myspecdata
    myASD10 = np.percentile(myspecdata2, 10, axis=0)/np.sqrt(np.log(2))
    myASD50 = np.percentile(myspecdata2, 50, axis=0)/np.sqrt(np.log(2))
    myASD90 = np.percentile(myspecdata2, 90, axis=0)/np.sqrt(np.log(2))
    myASDMean = np.sqrt(np.mean(myspecdata2**2, axis=0))
    if binNum != None:
        ff = logbin(ff, binNum)
        myASD10 = logbin(myASD10, binNum)
        myASD50 = logbin(myASD50, binNum)
        myASD90 = logbin(myASD90, binNum)
        myASDMean = logbin(myASDMean, binNum)
    if fullOutput==False:
        return ff, myASD10, myASD50, myASD90, myASDMean
    else:
        return ff, myASD10, myASD50, myASD90, myASDMean, myspecdata2

def aa_16k(ff):
    # coefficients for the 16 kHz IIR filtering
    g = 0.014805052402446
    a = [-1.71662585474518, 0.78495484219691, -1.68385964238855,    0.93734519457266]
    b = [-1.41346289716898,   0.99893884152400, 0.00000127375260,   0.99819981588176]
    # number of sample
    #ff = np.linspace(0, 7444, 7445)
    N = 7444 + 1
    f_end = 7444
    fs = 16384.0*4 # 64 kHz system
    z = np.exp(-1j*2*np.pi*ff/fs)
    # IOP decimation filter
    H = g * (1.0 + b[0]*z + b[1]*z**2)*(1.0 + b[2]*z + b[3]*z**2)\
        /(1.0 + a[0]*z + a[1]*z**2)/(1.0 + a[2]*z + a[3]*z**2)

    # ADC response
    if ff[0] == 0:
        digi = np.exp(-1j*np.pi*ff/fs)
        digi[1:] *= np.sin(np.pi*ff[1:]/fs)/(np.pi*ff[1:]/fs)
        digi[0] = digi[1]
    else:
        digi = np.exp(-1j*np.pi*ff/fs) * np.sin(np.pi*ff/fs)/(np.pi*ff/fs)

    aaTotal = digi*H
    aaTotal[0] = 1
    return aaTotal

def stitch_spec(flist):
    '''
    Stitch together multiple ASD traces.

    Must be given with frequency vectors in ascending order.

    '''
    ff_list = []
    asd_list = []
    ff_full = None
    asd_full = None
    for fname in flist:
        ff, asd = np.loadtxt(fname, unpack=1)
        if ff_full == None:
            ff_full = ff
            asd_full = asd
        else:
            mask = (ff>ff_full[-1])
            ff_full = np.concatenate([ff_full, ff[mask]])
            asd_full = np.concatenate([asd_full, asd[mask]])
    return np.array([ff_full, asd_full])

def freqresp(zpkTup, ff):
    '''
    Like scipy.signal.freqresp, but without having to pass around aggrivating factors
    of -2*np.pi.

    '''
    zArr, pArr, k = zpkTup
    zArr = np.array(zArr)
    pArr = np.array(pArr)
    _, tf = sig.freqresp((-2*np.pi*zArr, -2*np.pi*pArr,
        k*np.prod(-2*np.pi*pArr)/np.prod(-2*np.pi*zArr)), w=2*np.pi*ff)
    return ff, tf
'''
def get_ifo_data(channel, start, stop, savedir, full_output=False, **kwargs):
    savechan = channel.split(':')[1]
    savetime = start.replace(':', '').replace(' ', '_')
    savepath = savedir+savechan+'_'+savetime+'.txt'
    if os.path.isfile(savepath):
        print('Loading {}'.format(savepath))
        ff, asd_10, asd_50, asd_90, asd_mean = np.loadtxt(savepath, unpack=1)
    else:
        print('I did not find saved data at {}. Fetching data from NDS...'.format(savepath))
        data_ts = TimeSeries.fetch(channel, start, stop, verbose=True,
                host=nds_server, port=31200)
        print('Making ASDs...')
        ff, asd_10, asd_50, asd_90, asd_mean = get_med_asds(data_ts, **kwargs)
        header = channel+'\nData from {} to {}\n'.format(start, stop) \
                + 'Mean: {:.3g} \n'.format(data_ts.mean()) \
                + '\t'.join(['Frequency [Hz]',
                        '10th percentile ASD',
                        '50th percentile ASD',
                        '90th percentile ASD',
                        'Mean ASD',
                            ])
        # Undo the aa filter, if desired
        if 'aa' in kwargs.keys():
            print('Undoing AA filter...')
            if kwargs['aa'] == '16k':
                asd_10 /= np.abs(aa_16k(ff))
                asd_50 /= np.abs(aa_16k(ff))
                asd_90 /= np.abs(aa_16k(ff))
                asd_mean /= np.abs(aa_16k(ff))
            else:
                print('Got unknown value for "aa" keyword.')
        if 'calibration' in kwargs.keys():
            print('Applying calibration...')
            cal_interp = np.interp(ff, *kwargs['calibration'])
            asd_10 *= cal_interp
            asd_50 *= cal_interp
            asd_90 *= cal_interp
            asd_mean *= cal_interp
        np.savetxt(savepath,
                np.c_[ff,
                      asd_10,
                      asd_50,
                      asd_90,
                      asd_mean,
                      ],
                header=header)
    if full_output == True:
        return ff, asd_10, asd_50, asd_90, asd_mean
    else:
        return ff, asd_50
'''
def get_noise_coupling(stim_chan, resp_chan, inj_start, inj_stop, quiet_start, quiet_stop,
        proj_start, proj_stop, stim_name, resp_name, savedir, full_output=False, **kwargs):
    savepath = '_'.join([savedir, stim_name, resp_name, 'coupling']) + '.txt'
    if os.path.isfile(savepath):
        f_coup, coup_asd = np.loadtxt(savepath, unpack=1)
    else:
        f_proj, proj_asd = get_ifo_data(stim_chan, proj_start, proj_stop, savedir)
        f_tf, tf = get_whitenoise_tf(stim_chan, resp_chan, inj_start, inj_stop, quiet_start, quiet_stop,
                stim_name, resp_name, savedir)
        if not np.array_equal(f_proj, f_tf):
            tf = np.interp(f_proj, f_tf, tf, left=0, right=0)
        coup_asd = proj_asd * tf
        f_coup = f_proj
        if 'mask' in kwargs.keys():
            mask = kwargs['mask']
            mask = np.all([f_coup>mask[0], f_coup<mask[1]], axis=0)
            #f_coup = f_coup
            coup_asd = coup_asd*mask
    return f_coup, coup_asd

'''
# Generate paths for stimulus/response channel TF, and coupling ASD
savepath_base = '_'.join([stim_name, resp_name, 'coupling', inj_start.replace(' ', '_').replace(':', '')])
savepath_base_coup = '_'.join([stim_name, resp_name, 'coupling', proj_start.replace(' ', '_').replace(':', '')])
savepath_coup_asd = 'Data/' + savepath_base_coup + '_asd.txt'
savepath_tf = 'Data/' + savepath_base + '_tf.txt'
savepath_fig = 'Data/' + savepath_base + '_tf.txt'
if os.path.isfile(savepath_coup_asd) and os.path.isfile(savepath_tf):
    ff, coupling = np.loadtxt(savepath_asd, unpack=1)
    f_tf, tf = np.loadtxt(savepath_tf, unpack=1)
else:
    print('Trying to fetch {} and {} from {} to {}'.format(stim_chan, resp_chan, inj_start, inj_stop))
    inj_ts_dict = TimeSeriesDict.fetch([stim_chan, resp_chan], inj_start, inj_stop, verbose=True)
    print('Trying to fetch {} and {} from {} to {}'.format(stim_chan, resp_chan, quiet_start, quiet_stop))
    quiet_ts_dict = TimeSeriesDict.fetch([stim_chan, resp_chan], quiet_start, quiet_stop, verbose=True)
    proj_ts = TimeSeries.fetch(stim_chan, proj_start, proj_stop, verbose=True)
    f_stim, stim_inj_asd_10, stim_inj_asd_50, stim_inj_asd_90, stim_inj_asd_mean = \
            get_med_asds(inj_ts_dict[stim_chan], fftlength=10, overlap=5, method='welch',
                    window='hann', binNum=10000)
    f_resp, resp_inj_asd_10, resp_inj_asd_50, resp_inj_asd_90, resp_inj_asd_mean = \
            get_med_asds(inj_ts_dict[resp_chan], fftlength=10, overlap=5, method='welch',
                    window='hann', binNum=10000)
    _, stim_quiet_asd_10, stim_quiet_asd_50, stim_quiet_asd_90, stim_quiet_asd_mean = \
            get_med_asds(quiet_ts_dict[stim_chan], fftlength=10, overlap=5, method='welch',
                    window='hann', binNum=10000)
    _, resp_quiet_asd_10, resp_quiet_asd_50, resp_quiet_asd_90, resp_quiet_asd_mean = \
            get_med_asds(quiet_ts_dict[resp_chan], fftlength=10, overlap=5, method='welch',
                    window='hann', binNum=10000)
    print(np.any([f_resp != f_stim]))
    if np.any([f_resp != f_stim]):
        if f_resp[-1] > f_stim[-1]:
            print('Interpolating response ASD...')
            # Interpolate the response ASDs
            resp_inj_asd_50 = np.interp(f_stim, f_resp, resp_inj_asd_50)
            resp_quiet_asd_50 = np.interp(f_stim, f_resp, resp_quiet_asd_50)
            ff = f_stim
        else:
            # Interpolate the stimulus ASDs
            print('Interpolating stimulus ASD...')
            stim_inj_asd_50 = np.interp(f_resp, f_stim, stim_inj_asd_50)
            stim_quiet_asd_50 = np.interp(f_resp, f_stim, stim_quiet_asd_50)
            ff = f_resp
    else:
        ff = f_resp
    # Plot and save the white noise TF
    print('Plotting white noise TF from {} to {}'.format(stim_name, resp_name))
    h_whitenoise = plt.figure()
    ax_stim = h_whitenoise.add_subplot(211)
    ax_resp = h_whitenoise.add_subplot(212, sharex=ax_stim)
    ax_stim.loglog(ff, stim_quiet_asd_50, label='Quiet')
    ax_stim.loglog(ff, stim_inj_asd_50, label='Injection')
    ax_stim.legend()
    set_grid(ax_stim)
    ax_stim.set_ylabel('ASD')
    ax_stim.set_title(stim_name.replace('_', '\_'))
    ax_resp.loglog(ff, resp_quiet_asd_50, label='Quiet')
    ax_resp.loglog(ff, resp_inj_asd_50, label='Injection')
    ax_resp.set_xlabel('Frequency [Hz]')
    ax_resp.set_ylabel('ASD')
    ax_resp.set_xlim([7, 8e3])
    set_grid(ax_resp)
    ax_resp.set_title(resp_name.replace('_', '\_'))
    print('Saving white noise spectra plot...')
    h_whitenoise.tight_layout()
    h_whitenoise.savefig(savepath_fig)
    print('Saving white noise TF data...')
    # Get the white noise TF, save it, and save the projected ASD
    tf = np.sqrt(np.abs(resp_inj_asd_50**2-resp_quiet_asd_50**2) /
            np.abs(stim_inj_asd_50**2-stim_quiet_asd_50**2))
    coupling = stim_quiet_asd_50 * tf
    np.savetxt(savepath_tf, np.c_[ff, whitenoise_tf])
    np.savetxt(savepath_asd, np.c_[ff, coupling])
if full_output == True:
    return ff, coupling, f_tf, tf
else:
    return ff, coupling


def get_whitenoise_tf(stim_chan, resp_chan, inj_start, inj_stop, quiet_start, quiet_stop,
        stim_name, resp_name, savedir):
    savepath_base = savedir + '_'.join([stim_name, resp_name, 'whitenoise_tf', inj_start.replace(' ', '_').replace(':', '')])
    savepath_data = savepath_base + '.txt'
    savepath_fig = 'Figures/' + savepath_base + '.pdf'
    if os.path.isfile(savepath_data):
        print('Loading locally saved coupling TF...')
        return np.loadtxt(savepath_data, unpack=1)
    else:
        print('I did not find a saved TF at {}. Getting data from NDS instead.'.format(savepath_data))
        print('Trying to fetch {} and {} from {} to {}'.format(stim_chan, resp_chan, inj_start, inj_stop))
        inj_ts_dict = TimeSeriesDict.fetch([stim_chan, resp_chan], inj_start, inj_stop,
                host=nds_server, port=nds_port, verbose=True)
        print('Trying to fetch {} and {} from {} to {}'.format(stim_chan, resp_chan, quiet_start, quiet_stop))
        quiet_ts_dict = TimeSeriesDict.fetch([stim_chan, resp_chan], quiet_start, quiet_stop,
                host=nds_server, port=nds_port, verbose=True)
        f_stim, stim_inj_asd_10, stim_inj_asd_50, stim_inj_asd_90, stim_inj_asd_mean = \
                get_med_asds(inj_ts_dict[stim_chan], fftlength=10, overlap=5, method='welch',
                        window='hann', binNum=10000)
        f_resp, resp_inj_asd_10, resp_inj_asd_50, resp_inj_asd_90, resp_inj_asd_mean = \
                get_med_asds(inj_ts_dict[resp_chan], fftlength=10, overlap=5, method='welch',
                        window='hann', binNum=10000)
        _, stim_quiet_asd_10, stim_quiet_asd_50, stim_quiet_asd_90, stim_quiet_asd_mean = \
                get_med_asds(quiet_ts_dict[stim_chan], fftlength=10, overlap=5, method='welch',
                        window='hann', binNum=10000)
        _, resp_quiet_asd_10, resp_quiet_asd_50, resp_quiet_asd_90, resp_quiet_asd_mean = \
                get_med_asds(quiet_ts_dict[resp_chan], fftlength=10, overlap=5, method='welch',
                        window='hann', binNum=10000)
        print(np.any([f_resp != f_stim]))
        if np.any([f_resp != f_stim]):
            if f_resp[-1] > f_stim[-1]:
                print('Interpolating response ASD...')
                # Interpolate the response ASDs
                resp_inj_asd_50 = np.interp(f_stim, f_resp, resp_inj_asd_50)
                resp_quiet_asd_50 = np.interp(f_stim, f_resp, resp_quiet_asd_50)
                ff = f_stim
            else:
                # Interpolate the stimulus ASDs
                print('Interpolating stimulus ASD...')
                stim_inj_asd_50 = np.interp(f_resp, f_stim, stim_inj_asd_50)
                stim_quiet_asd_50 = np.interp(f_resp, f_stim, stim_quiet_asd_50)
                ff = f_resp
        else:
            ff = f_resp
        print('Plotting white noise TF from {} to {}'.format(stim_name, resp_name))
        h_whitenoise = plt.figure()
        ax_stim = h_whitenoise.add_subplot(211)
        ax_resp = h_whitenoise.add_subplot(212, sharex=ax_stim)
        ax_stim.loglog(ff, stim_quiet_asd_50, label='Quiet')
        ax_stim.loglog(ff, stim_inj_asd_50, label='Injection')
        ax_stim.legend()
        set_grid(ax_stim)
        ax_stim.set_ylabel('ASD')
        ax_stim.set_title(stim_name.replace('_', '\_'))
        ax_resp.loglog(ff, resp_quiet_asd_50, label='Quiet')
        ax_resp.loglog(ff, resp_inj_asd_50, label='Injection')
        ax_resp.set_xlabel('Frequency [Hz]')
        ax_resp.set_ylabel('ASD')
        ax_resp.set_xlim([7, 8e3])
        set_grid(ax_resp)
        ax_resp.set_title(resp_name.replace('_', '\_'))
        print('Saving white noise spectra plot...')
        h_whitenoise.tight_layout()
        h_whitenoise.savefig(savepath_fig)
        print('Saving white noise TF data...')
        whitenoise_tf = np.sqrt(np.abs(resp_inj_asd_50**2-resp_quiet_asd_50**2) /
                np.abs(stim_inj_asd_50**2-stim_quiet_asd_50**2))
        np.savetxt(savepath_data, np.c_[ff, whitenoise_tf])
        return ff, whitenoise_tf
'''
def cum_rms(ff, asd):
    df = np.roll(ff, -1) - ff
    df = df[:-1]
    cum_ms = np.zeros(len(ff))
    for ii in range(len(cum_ms)):
        cum_ms[ii] = np.trapz(asd[ii:]**2, ff[ii:])
    return cum_ms**0.5


def get_ifo_data_nds(channel, start, stop, **kwargs):
    '''
    For grabbing channel data using NDS rather than GWPy
    '''
    conn = nds2.connection('nds40.ligo.caltech.edu', 31200)
    data = conn.fetch(start, stop, [channel])
    fs = data[0].channel.sample_rate  #Hz, downsampled everything to this frequency...
    t_fft = 20
    nperseg = t_fft*fs
    window = sig.get_window('hanning', int(nperseg))
    ff1, psd = sig.welch(data[0].data, fs, window, int(nperseg))
    return ff1,np.sqrt(psd)
'''
def get_OL_coupling(start, stop, chans, cpl_mag, cpl_phase):
    '''
'''
    For estimating the Oplev coupling to a given channel in loop
'''
'''
    cpl_mag = 10.**(cpl_mag/20.0)
    cpl_phase = np.deg2rad(cpl_phase)
    cpl_cmplx = cpl_mag * np.exp(1j*cpl_phase)
    conn = nds2.connection('nds40.ligo.caltech.edu', 31200)#coupling gains to dB
    data = conn.fetch(start, stop, chans)
    fs = int(data[1].channel.sample_rate)  #Hz, downsampled everything to this frequency...
    t_fft = 20
    nperseg = t_fft*fs
    window = sig.get_window('hanning', nperseg)
    pMICH_OL_coh = 0
    for ii in range(0,np.shape(cpl_cmplx)[0]):
        ff_temp, psd_temp = sig.welch(data[ii].data,fs, window, int(nperseg))
        asd_temp = np.sqrt(psd_temp)
        pMICH_OL_coh += (asd_temp*cpl_cmplx[ii])**2
        #pMICH_OL_coh += (psd_temp*cpl_cmplx[ii])
    pMICH_OL_coh = np.sqrt(np.abs(pMICH_OL_coh))
    return ff_temp, pMICH_OL_coh
'''
