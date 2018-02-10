'''
Python script that calculates the THD for some number of harmonics
and makes a nice plot of the power spectrum with harmonics marked.
Algorithm:
	1. Get data (either online or offline)
	2. Compute power spectrum (Kaiser window, beta=38)
	3. Use scipy.signal.find_peaks_cwt to find peaks in the power spectrum
	4. Return all prominent peaks that are not powerline harmonics etc.
Need to figure out error handline
'''

import numpy as np
import scipy.signal as sig
import matplotlib.pyplot as plt
import h5py

plt.style.use('gvELOG')

def thd(x,nAvg,fs=16384,N=5):
	#Some parameters for the peak finding function
	wid = np.arange(1,10)
	#wid = np.array([1,2,3,4])
	min_snr = 50
	noise_perc = 5 

	#Compute the power spectrum
	print('Computing power spectrum...')
	nfft = len(x)/nAvg
	w = sig.kaiser(nfft,38)
	ff, Pxx = sig.periodogram(x, fs, w, scaling='spectrum',nfft=nfft)

	#Find the peaks 
	print('Finding peaks....')
	peak_ind = sig.find_peaks_cwt(Pxx, wid, min_snr=min_snr, noise_perc=noise_perc)
	if len(peak_ind) < N:
		print('Not enough peaks of sufficient SNR found - check measurement!')
		return
	
	print('Found {} peaks with SNR > {} ...'.format(len(peak_ind),min_snr))
	
	#Eliminate powerline harmonics
	for ii, pp in enumerate(peak_ind):
		if np.mod(ff[pp],60)==0:
			peak_ind[ii] = 0
	peak_ind = peak_ind[peak_ind!=0]
    
	fund = ff[peak_ind[np.argmax(Pxx[peak_ind])]]
	print('Fundamental frequency is {} Hz'.format(np.round(fund,decimals=2)))
    
	#Eliminate peaks below the fundamental frequency
	for ii, pp in enumerate(peak_ind):
		if ff[pp]<fund:
			peak_ind[ii] = 0
	peak_ind = peak_ind[peak_ind!=0]
	
	#Eliminate non-harmonic peaks
	for ii, pp in enumerate(peak_ind):
		if 1 < np.mod(np.round(ff[pp],decimals=2),np.round(fund,decimals=2)) < fund-1:
			peak_ind[ii] = 0
	peak_ind = peak_ind[peak_ind!=0]
	
	#Get the N largest peaks
	#peak_ind = sorted(peak_ind, reverse=True)
	ff_pk = ff[peak_ind[np.argsort(-1*Pxx[peak_ind])]][0:N]
	Pxx_pk = sorted(Pxx[peak_ind], reverse=True)[0:N]

	#Calculate the THD
	print('List of {} highest peaks are:'.format(N))
	for ii, fff in enumerate(ff_pk):
		print('f = {} Hz, Vsq = {}'.format(np.round(ff_pk[ii],decimals=2),np.round(Pxx_pk[ii], decimals=2)))
	
	thd = np.sqrt(np.sum(Pxx_pk[1:]))/np.sqrt(Pxx_pk[0])
	print('Measured THD with {} harmonics is {}%'.format(N, 1e2*thd))

	fig, ax = plt.subplots(1,1,figsize=(12,8))
	#ax.loglog(ff, Pxx, rasterized=True, alpha=0.6, label='Measurement', color='xkcd:olive green')
	#ax.loglog(ff_pk[0],Pxx_pk[0],'o',markersize=20,alpha=0.6,fillstyle='none', label='Fundamental', color='xkcd:blood')
	#ax.loglog(ff_pk[1:],Pxx_pk[1:],'o',markersize=20,alpha=0.6,fillstyle='none',label='Harmonics', color='xkcd:sky blue')
	ax.semilogy(ff, Pxx, rasterized=True, alpha=0.6, label='Measurement', color='xkcd:olive green')
	ax.semilogy(ff_pk[0],Pxx_pk[0],'o',markersize=20,alpha=0.6,fillstyle='none', label='Fundamental', color='xkcd:blood',linewidth=3)
	ax.semilogy(ff_pk[1:],Pxx_pk[1:],'o',markersize=20,alpha=0.6,fillstyle='none',label='Harmonics', color='xkcd:sky blue', linewidth=3)
	ax.grid('on',which='both',linestyle='--',alpha=0.4)
	ax.legend(loc='best')
	ax.set_xlabel('Frequency [Hz]')
	ax.set_ylabel('Power Spectrum [$\mathrm{cts}^2$]')
	ax.set_title('THD measurement, $f_0 = ${} Hz, N = {}'.format(np.round(fund, decimals=2),N))
	ax.set_ylim([0.01*Pxx_pk[-1], 10*Pxx_pk[0]])
	return ff_pk, Pxx_pk
	
def main():
	f = h5py.File('/users/gautam/2018_02_D990694THD/ALS_Y_I_THD.hdf5','r')
	#fs = 16384
	#Downsample for speed
	fs = 16384/8
	x = sig.decimate(f['data'][:], 8, zero_phase=True)
	a,b = thd(x,10, N=8, fs=fs)
	return

if __name__ == "__main__":
	main()
