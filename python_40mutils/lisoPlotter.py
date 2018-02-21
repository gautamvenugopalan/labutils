#!/usr/local/bin/python
'''
Function to plot liso .out files
Example usage:
	python lisoPlotter 'file1.out,file2.out,...,fineN.out'
'''

import numpy as np 
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import FormatStrFormatter

fileNames = sys.argv[1].split(',')
legends = sys.argv[2].split(',')

plt.style.use('gvELOG')
fig,ax = plt.subplots(2,1,sharex=True,figsize=(12,8))

for ii,f in enumerate(fileNames):
	dat = np.loadtxt(f)
	ax[0].semilogx(dat[:,0],dat[:,1],label=legends[ii])
	ax[1].semilogx(dat[:,0],dat[:,2],label=legends[ii])

#Beautification
ax[1].set_ylim([-180,180])
ax[1].set_yticks(np.linspace(-180,180,9))
ax[1].yaxis.set_major_formatter(FormatStrFormatter('%d'))
ax[1].set_xlabel('Frequency [Hz]')
ax[1].set_ylabel('Phase [deg]')
ax[0].set_ylabel('Magnitude [dB]')
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%d'))
ax[1].legend(loc='best')
plt.suptitle(sys.argv[3])

plt.savefig('./lisoPlot.pdf')
