'''
Collection of generally useful functions for making ELOG plots etc
'''
import sys, os, subprocess
import numpy as np

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


