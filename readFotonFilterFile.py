'''
Python transcription of readFilterFile.m
Example usage:
	filtDict = readFilterFile(fileName)
returns a dictionary with the stuff in the 
Foton .txt file specified as fileName
'''

import numpy as np
import scipy.signal as sig

def readFilterFile(fileName):
	try:
		f = open(fileName,'r')
	except:
		print('No such filter file!!')
		return
	#Initialize the dictionary
	p = {}
	p['fileName'] = fileName
	#Loop through file
	s = f.readline()
	while 1:
		if not s: break
		arg = s.split()
		#Check type of line
		if len(arg)<2:
			pass
		elif len(arg)>2 and arg[0]=='#':
			if arg[1]=='MODULES':
				#initialization of all filters
				for ii in range(2,len(arg)):
					p[arg[ii]] = emptyFilter()
			elif arg[1]=='SAMPLING':
				#Sampling rate of this RT model
				p['fs'] = float(arg[3])
			elif arg[1]=='DESIGN':
				fname = arg[2]
				ind = int(arg[3])
				p[fname][ind]={}
				p[fname][ind]['design'] = ''.join(arg[4:len(arg)+1])
		elif len(arg)==12:
			#This is an actual filter SOS definition
			fname=arg[0]
			index = int(arg[1])
			p[fname][index]['name']=arg[6]
			p[fname][index]['gain']=float(arg[7])
			gain = float(arg[7])
			soscoef=np.array(arg[8:12],dtype='float')
			order = int(arg[3])
			for kk in range(0,order-1):
				s = f.readline()
				arg = s.split()
				temp = np.array(arg[0:4],dtype='float')
				soscoef = np.vstack((soscoef,temp))
			if order==1:
				soscoef = np.vstack((soscoef,soscoef)) #for indexing convenience
			#reshape the SOS coeffs
			coef = np.zeros((order,6))
			for jj in range(order):
				if jj==0:
					coef[jj][0] = gain
					coef[jj][1] = gain*soscoef[jj][2]
					coef[jj][2] = gain*soscoef[jj][3]
				else:
					coef[jj][0] = 1.
					coef[jj][1] = soscoef[jj][2]
					coef[jj][2] = soscoef[jj][3]
				coef[jj][3] = 1.
				coef[jj][4] = soscoef[jj][0]
				coef[jj][5] = soscoef[jj][1]
			p[fname][index]['sosCoeffs'] = np.squeeze(coef)
		s = f.readline()

	f.close()
	return p

def emptyFilter():
	'''
	Returns an "empty" filter bank structure for initialization
	'''
	fb = {}
	for ii in range(10):
		fb[ii] = {}
		fb[ii]['name']='<empty>'
		fb[ii]['sosCoeffs'] = np.array([1.,0.,0.,1.,0.,0.])
		fb[ii]['fs'] = 16384
		fb[ii]['design'] = 'none'
	return fb
	


