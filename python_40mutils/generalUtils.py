'''
Collection of generally useful functions for making ELOG plots etc
'''
import sys, os, subprocess
import numpy as np
import scipy.signal as sig
import nds2

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

def makeLSCspecgram(channel, tStart, tStop, t_fft=4, win=('tukey', 0.25)):
	'''
	Computes a spectrogram by downloading data from NDS2.
	Example usage:
		tt, ff, Sxx = makeLSCspecgram(['C1:LSC-XARM_IN1_DQ'],1208808916,1208808917, t_fft=4, win=('tukey', 0.25))
	will download the data, and make a spectrogram with 4 second FFT time.
	Option to specify windowing function used by scipy.signal.spectrogram.
	Note that Power Spectral Density is returned, use sqrt where appropriate!
	'''
	conn = nds2.connection('nds40.ligo.caltech.edu', 31200)
	print('Getting NDS data')
	dat = conn.fetch(tStart, tStop, channel)
	print('Got NDS data...')
	Fs = dat[0].channel.sample_rate
	ff, tt, Sxx = sig.spectrogram(dat[0].data,fs=Fs,nperseg=int(t_fft*Fs),window=win)
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
