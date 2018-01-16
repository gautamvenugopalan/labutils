'''
Python script that fits a Gaussian to 
Oplev Error Signal vs Cavity Power timeseries 
Assumption is that scan is done slow enough
to be outside OL bandwidth.

Usage:
	python OLcalibFactor.py tStart tStop optic DoF

ToDo: Need to add more formulae for various optic tilts
GV 18Nov2017

Changelog:
	gautamv 15 Jan 2018: Added functions for ITMs, ETMs and BS
'''

import numpy as np
from scipy.optimize import curve_fit
import sys
import ezca
import cdsutils as cds
import matplotlib.pyplot as plt
from matplotlib import rc

z = ezca.Ezca()
if 'gvELOG' in plt.style.available:
	plt.style.use('gvELOG')
else:
	plt.style.use('bmh')

rc('text',usetex=False)

#Global IFO params
lam = 1064e-9  #Laser wavelength [m]
L_arm = 37.79  #Arm cavity length [m]
R_ETMX = 59.48 #ETMX RoC, ATF measurement [m]
R_ETMY = 60.26 #ETMY RoC, ATF measurement [m]

#Parse command line arguments
tStart = long(sys.argv[1])
tStop = long(sys.argv[2])
opt = sys.argv[3]
dof = sys.argv[4]

#In general, 
# a = overall scaling for the arm transmission
# b = shift of origin of angular misalignment
# c = scale factor to convert Oplev misalignment (error sig) to urad
# x = oplev error signal during misalignment

def tilt_BS(x,a,b,c):
	#Expected variation of arm power as function of BS tilt.
	#c is the scaling factor for converting current OL readout, x
	#to the physical tilt, in rad. 
	#Factor of 1e-6 is for enabling good fitting
	#Assumption is that translation of cavity axis has negligible effect.
	w0 = 3e-3      #Arm cavity waist [m]
	alpha0 = lam / np.pi / w0
	return a*np.exp(-4* (1e-6*c*(x-b)/alpha0)**2) #Factor of 4 because tilt of BS shows up as tilt of cavity axis by twice as much

#functions for other optics are defined in an analogous manner to the above...
def tilt_ITMX(x,a,b,c):
	w0 = np.sqrt(lam/np.pi) * ((L_arm*(R_ETMX-L_arm))**0.25)
	alpha0 = lam / np.pi / w0
	rot = 1e-6 * c * (x-b)         #Tilt in urad
	trans = (R_ETMX - L_arm) * rot #Shift of beam axis
	return a * np.exp(-1.*(rot/alpha0)**2) * np.exp(-1.*(trans/w0)**2)

def tilt_ITMY(x,a,b,c):
	w0 = np.sqrt(lam/np.pi) * ((L_arm*(R_ETMY-L_arm))**0.25)
	alpha0 = lam / np.pi / w0
	rot = 1e-6 * c * (x-b)         #Tilt in urad
	trans = (R_ETMY - L_arm) * rot #Shift of beam axis
	return a * np.exp(-1.*(rot/alpha0)**2) * np.exp(-1.*(trans/w0)**2)

def tilt_ETMX(x,a,b,c):
	w0 = np.sqrt(lam/np.pi) * ((L_arm*(R_ETMX-L_arm))**0.25)
	alpha0 = lam / np.pi / w0
	rot = 1e-6 * c * (x-b)         #Tilt in urad
	trans = R_ETMX * rot           #Shift of beam axis
	return a * np.exp(-1.*(trans/w0)**2)

def tilt_ETMY(x,a,b,c):
	w0 = np.sqrt(lam/np.pi) * ((L_arm*(R_ETMY-L_arm))**0.25)
	alpha0 = lam / np.pi / w0
	rot = 1e-6 * c * (x-b)         #Tilt in urad
	trans = R_ETMY * rot           #Shift of beam axis
	return a * np.exp(-1.*(trans/w0)**2)

def fitOLcalibFactor(tStart, tStop, opt, dof):
	if opt in ['ITMX','ETMX','BS']:
		chans = ['SUS-'+opt+'_OPLEV_'+dof, 'LSC-TRX_OUT_DQ']
	elif opt in ['ITMY','ETMY']:
		chans = ['SUS-'+opt+'_OPLEV_'+dof, 'LSC-TRY_OUT_DQ']
	tDur = int(tStop - tStart)
	print('Grabbing data for {} and {} from {} to {}'.format(chans[0],chans[1],tStart, tStop))
	data = cds.getdata(chans, tDur, long(tStart))
	x = data[0].data
	y = data[1].data[::8] #Re-sampling 16k data
	#Fit the data to the model
	if opt=='BS':
		popt,pcov = curve_fit(tilt_BS, x,y, p0=(1,1,10))
	elif opt=='ITMX':
		popt, pcov = curve_fit(tilt_ITMX, x,y, p0=(1,1,10))
	elif opt=='ITMY':
		popt, pcov = curve_fit(tilt_ITMY, x,y, p0=(1,1,10))
	elif opt=='ETMX':
		popt, pcov = curve_fit(tilt_ETMX, x,y, p0=(1,1,10))
	elif opt=='ETMY':
		popt, pcov = curve_fit(tilt_ETMY, x,y, p0=(1,1,10))
	#Make a plot
	fig,ax = plt.subplots(2,1,sharex=True,figsize=(12,8))
	ax[0].plot(x,y,'.',alpha=0.04, label='Data',rasterized=True, color='xkcd:wine red')
	if opt=='BS':
		ax[0].plot(x,tilt_BS(x,*popt),'-k', label='fit to model',linewidth=3)
		ax[1].plot(x,y-tilt_BS(x,*popt),'.', color='xkcd:blush', alpha=0.03, label='residuals', rasterized=True)
	elif opt=='ITMX':
		ax[0].plot(x,tilt_ITMX(x,*popt),'-k', label='fit to model',linewidth=3)
		ax[1].plot(x,y-tilt_ITMX(x,*popt),'.', color='xkcd:blush',  alpha=0.03, label='residuals', rasterized=True)
	elif opt=='ITMY':
		ax[0].plot(x,tilt_ITMY(x,*popt),'-k', label='fit to model',linewidth=3)
		ax[1].plot(x,y-tilt_ITMY(x,*popt),'.', color='xkcd:blush',  alpha=0.03, label='residuals', rasterized=True)
	elif opt=='ETMX':
		ax[0].plot(x,tilt_ETMX(x,*popt),'-k', label='fit to model',linewidth=3)
		ax[1].plot(x,y-tilt_ETMX(x,*popt),'.', color='xkcd:blush',  alpha=0.03, label='residuals', rasterized=True)
	elif opt=='ETMY':
		ax[0].plot(x,tilt_ETMY(x,*popt),'-k', label='fit to model',linewidth=3)
		ax[1].plot(x,y-tilt_ETMY(x,*popt),'.', color='xkcd:blush',  alpha=0.03, label='residuals', rasterized=True)
	
	ax[1].set_xlabel('Uncalibrated tilt [a.u.]')
	ax[0].set_ylabel('{} [a.u.]'.format(chans[1]))
	ax[0].set_title('{} vs {} from {} to {}'.format(chans[1], chans[0], tStart, tStop))
	ax[1].set_ylabel('Fit Residuals [a.u.]')
	ax[0].legend(loc='best')
	plt.savefig('Figures/OL_calib_'+opt+'_'+dof+'.pdf',dpi=80)
	#print the result
	if dof=='PERROR':
		currentValue = z.read('SUS-'+opt+'_OL_PIT_CALIB')
	elif dof=='YERROR':
		currentValue = z.read('SUS-'+opt+'_OL_YAW_CALIB')
	proposedValue = np.abs(currentValue*popt[2])
	print('Fit yielded the following values: a={}, b={}, c={}'.format(popt[0],popt[1],popt[2]))
	print('Current OL calibration is {}. Fit says it should be {}. Scaling factor is {}'.format(currentValue, proposedValue, np.abs(popt[2])))
	return

fitOLcalibFactor(tStart,tStop,opt,dof)


