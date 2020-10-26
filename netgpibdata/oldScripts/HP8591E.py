#!/usr/bin/env python

# this module allows setting of simple scan parameters over GPIB for HP8591E spectrum analyzer 
# Nichin Sreekantaswamy: July 2014

import netgpib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import argparse
import yaml
import time
import math

#==================================================================
# Read data from the .yml file 
#==================================================================
def readParams(paramFile):
    with open(paramFile,'r') as f:
        reader = yaml.load_all(f)
        params = reader.next()
        reader.close()
    return(params)

#==================================================================

#Usage text
usage = """usage: %prog [options]

This command will run a scan on the HP8591E spectrum analyzer and download the result(s) as a (series of) .dat files.  This can be used to increase the effective resolution of the NA and to perform a scan with variable resolution over different frequency ranges. Enter all your parameters in the .yml file... follow the param_NWAG4395A.yml and param_[PDname].yml
"""

class optionsclass:
	pass

#==================================================================
#Check if the numpy version is greater than 1.6.0
#==================================================================	
def version_cmp(version1, version2):
   	parts1 = [int(x) for x in version1.split('.')]
	parts2 = [int(x) for x in version2.split('.')]

	# fill up the shorter version with zeros ...
	lendiff = len(parts1) - len(parts2)
	if lendiff > 0:
	        parts2.extend([0] * lendiff)
	elif lendiff < 0:
	        parts1.extend([0] * (-lendiff))

	for i, p in enumerate(parts1):
	        ret = cmp(p, parts2[i])
        	if ret: return ret
    	return 0

#==================================================================

#==================================================================
# Plot data 
#==================================================================
def plotstuff(dbmdata,freq,outDir,peakF1,peakA1,peakF2,peakA2):
	matplotlib.rc('xtick',labelsize =17)
	matplotlib.rc('ytick',labelsize =17)
	a1=plt.plot(freq,dbmdata)
	plt.setp(a1, color = 'g')
        plt.title('HP8591E Spectrum Analyzer\n'+'Peak1: '+str(peakA1)+'dBm @ '+str(freq[peakF1])+'MHz\n'+'Peak2: '+str(peakA2)+'dBm @ '+str(freq[peakF2])+'MHz',fontsize=16)
	plt.ylabel ('Magnitude [dBm]',fontsize=20)
	plt.xlabel ('Frequency [MHz]',fontsize=20)
	plt.ylim((-70,-10))
	plt.grid()
        plt.plot(freq[peakF1],dbmdata[peakF1],'gD')
        plt.plot(freq[peakF2],dbmdata[peakF2],'yD')
	#plt.tight_layout(pad=0.4,h_pad=0.1,w_pad=0.1)
	
        plt.savefig(outDir+'/HP8591E_View.pdf', format='pdf', bbox_inches='tight')
	plt.clf()

#==================================================================
#Main function 
#==================================================================
def main(ymlFile):
   vercmp = version_cmp('1.6.0',np.__version__)
   if vercmp == 1:
        print'numpy version is not compatible with data plotting code'
   else:
	params=readParams(ymlFile)
	options = optionsclass()
	
	ipAddress = params['ipAddress']
	gpibAddress = int(params['gpibAddress'])


	print('Connecting to '+str(ipAddress))

	#Create a netGPIB object
	gpibObj = netgpib.netGPIB(ipAddress, gpibAddress, '\004',0)
	time.sleep(0.1)

	print('Connected!')

	startF =int(params['startFreq'])
	endF = int(params['stopFreq'])
	gpibObj.command('IP')
	time.sleep(0.1)
	
	freq = np.linspace(startF,endF,num=401) #create frequency points
	print('Start frequency: '+str(startF)+'MHz')
	gpibObj.command('FA '+str(startF)+'MZ')
	time.sleep(0.1)

	print('Stop frequency: '+str(endF)+'MHz')
	gpibObj.command('FB '+str(endF)+'MZ')
	time.sleep(0.1)
	
	gpibObj.command('SNGLS')
	time.sleep(0.1)
	
	try:
	   	while 1:
			gpibObj.command('TS')
			time.sleep(0.1)
			dbmdata = np.fromstring(gpibObj.query('TRA?'),sep=',') #get data points
                        peakA1 = dbmdata[1]
                        peakA2 = dbmdata[1]
                        peakF1 = 1
                        peakF2 = 1
                        for j in range(2,400):
                            if dbmdata[j] > dbmdata[j-1] and dbmdata[j] > dbmdata[j+1]:
                                if peakA1 < dbmdata[j]:
                                    peakA2 = peakA1
                                    peakF2 = peakF1
                                    peakA1 = dbmdata[j]
                                    peakF1 = j
                                elif peakA2 < dbmdata[j]: 
                                    peakA2 = dbmdata[j]
                                    peakF2 = j  

			plotstuff(dbmdata,freq,params['saveDir'],peakF1,peakA1,peakF2,peakA2)
	except KeyboardInterrupt:
		gpibObj.command('IP')
		time.sleep(0.1)
		gpibObj.command('FA '+str(startF)+'MZ')
		time.sleep(0.1)
		gpibObj.command('FB '+str(endF)+'MZ')
		time.sleep(0.1)
		gpibObj.command('CONTS')
		time.sleep(0.1)

if __name__ == "__main__":
    templateFile = '/opt/rtcds/caltech/c1/scripts/general/netgpibdata/HP8591E_param.yml' 
    parser = argparse.ArgumentParser() 
    parser.add_argument('paramFile', nargs='?',
                        help = 'The Parameter file for the measurement.' \
                               ' If not specified, uses the template parameter'\
                               ' file in scripts/general/netgpibdata/',
                        default=templateFile)
    parser.add_argument('--template', help='Copy template parameter file to'\
                                           ' current dir, no measurement is'\
                                           ' made.',
                        action='store_true')
    args = parser.parse_args()
    if args.template:
        import shutil
        print 'Copying ' +templateFile+ ' to ' + os.getcwd()
        shutil.copyfile(templateFile, os.getcwd()+'/param_HP8591E.yml')
        print 'Done!'
    else:    
        main(args.paramFile)
#==================================================================


