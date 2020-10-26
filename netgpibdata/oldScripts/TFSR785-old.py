#! /usr/bin/env python
#
# TFSR785.py [-f filename] [-i ip_address] [-a gpib_address] [{-s|--startfreq} start_freq]
#                 [{-e|--stopfreq} stop_freq] [-n|--numpoints num_of_points] 
#                 [{-x|--excamp} excitation_amplitude] [{-c|--settlecycle} settle_cycle]
#                 [{-t|--intcycle} integration_cycle] 
#                 
# Command a SR785 to execute a transfer function measurement.
#
# Yoichi Aso  Sep 22 2008

import re
import sys
import math
from optparse import OptionParser
from socket import *
import time
import pdb
import gpib

#Parse options
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="Output file name without an extension", default="data")
parser.add_option("-i", "--ip",
                  dest="ipAddress", default="gpib01",
                  help="IP address/Host name")
parser.add_option("-a", "--address",
                  dest="gpibAddress", type="int", default=10,
                  help="GPIB device address")
parser.add_option("-s", "--startfreq",
                  dest="startFreq", default="100kHz",
                  help="Start frequency")
parser.add_option("-e", "--stopfreq",
                  dest="stopFreq", default="10kHz",
                  help="Stop frequency")
parser.add_option("-n", "--numpoints",
                  dest="numOfPoints", type="int", default="50",
                  help="Number of frequency points")
parser.add_option("-x", "--excamp",
                  dest="excAmp", default="1mV",
                  help="Excitation amplitude")
parser.add_option("-c", "--settlecycle",
                  dest="settleCycles", type="int", default="10",
                  help="Settle cycles")
parser.add_option("-t", "--intcycle",
                  dest="intCycles", type="int",default="20",
                  help="Integration cycles")
(options, args) = parser.parse_args()

# Open socket

print('Connecting to '+str(options.ipAddress)+' ...'),

netAddr=(options.ipAddress, 1234)
netSock = socket(AF_INET, SOCK_STREAM)
netSock.connect(netAddr)
netSock.setblocking(0) # Set non-blocking socket
print('done.')

#open files
dataFileName=options.filename+'.dat'
paramFileName=options.filename+'.par'
dataFile = open(dataFileName,'w')
paramFile = open(paramFileName,'w')

print('Data will be written into '+dataFileName)
print('Parameters will be written into '+paramFileName+'\n')
print('Setting up parameters for the measurement') 
#pdb.set_trace()
#
#Prepare for the TF measurement
netSock.send("++addr "+str(options.gpibAddress)+"\n")
time.sleep(0.1)
netSock.send("++eos 3\n")
time.sleep(0.1)
netSock.send("++eot_enable 0\n")
time.sleep(0.1)
netSock.send("++mode 1\n")
time.sleep(0.1)
netSock.send("++clr\n")
time.sleep(0.1)
netSock.send("++auto 0\n")
time.sleep(0.1)
netSock.send("++ifc\n")
time.sleep(0.2)
netSock.send("++read_tmo_ms 4000\n")
time.sleep(0.1)
netSock.send("++addr "+str(options.gpibAddress)+"\n")
time.sleep(0.1)



#Set output to GPIB
netSock.send("OUTX0\n")
time.sleep(0.1)

# Set dump destination
netSock.send("POUT2\n") # ASCII DUMP Mode
time.sleep(0.1)
netSock.send("PDST3\n") # Destination = GPIB
time.sleep(0.1)
netSock.send("PCIC0\n") # GPIB Control = Host

#Set measurement parameters
netSock.send('DFMT1\n') # Dual display
time.sleep(0.1)
netSock.send('ACTD0\n') # Active display 0
time.sleep(0.1)

netSock.send('MGRP2,3\n') # Measurement Group = Swept Sine
time.sleep(0.1)
netSock.send('MEAS2,47\n') # Frequency Resp
time.sleep(0.1)
netSock.send('VIEW0,0\n') # Disp 0 = LogMag
time.sleep(0.1)
netSock.send('VIEW1,5\n') # Dsip 1 = Phase
time.sleep(0.1)
netSock.send('UNDB0,1\n') # dB ON
time.sleep(0.1)
netSock.send('UNPK0,0\n') # PK Unit Off
time.sleep(0.1)
netSock.send('UNDB1,0\n') # dB OFF
time.sleep(0.1)
netSock.send('UNPK1,0\n') # PK Unit Off
time.sleep(0.1)
netSock.send('UNPH1,0\n') # Phase Unit deg.
time.sleep(0.1)
netSock.send('DISP0,1\n') # Live display on
time.sleep(0.1)
netSock.send('DISP1,1\n') # Live display on
time.sleep(0.1)

netSock.send('SSCY2,'+str(options.settleCycles)+'\n') # Settle cycles
time.sleep(0.1)
netSock.send('SICY2,'+str(options.intCycles)+'\n') # Integration cycles
time.sleep(0.1)

#pdb.set_trace()
netSock.send('SSTR2,'+options.startFreq+'\n') #Start frequency
time.sleep(0.1)
netSock.send('SSTP2,'+options.stopFreq+'\n') #Stop frequency
time.sleep(0.1)
netSock.send('SNPS2,'+str(options.numOfPoints)+'\n') #Stop frequency
time.sleep(0.1)

netSock.send('SSAM'+options.excAmp+'\n') #Source Amplitude
time.sleep(0.1)

netSock.send('A1RG1\n') # Ch1 AutoRange On
time.sleep(0.1)
netSock.send('A2RG1\n') # Ch2 AutoRange On
time.sleep(0.1)

#Start measurement
print('Transfer function measurement started ... ')
netSock.send("++eot_enable 1\n")
netSock.send("++eot_char 4\n")
time.sleep(0.1)
netSock.send('STRT\n') #Source Amplitude
time.sleep(0.1)

#Wait for the measurement to end
measuring = True
while measuring:
    netSock.send('DSPS?4\n')  #Get status 
    netSock.send('++read eoi\n')
    measuring = not int(gpib.gpibGetData(netSock, 10, '\004'))
    time.sleep(1)
        
print('done')


#Download Data
#pdb.set_trace()
netSock.send('++clr\n')
netSock.send('++addr '+str(options.gpibAddress)+'\n')
netSock.send("++auto 0\n")  #Turn off automatic mode
netSock.send("++eot_char 4\n")
time.sleep(0.1)
netSock.send("++eot_enable 1\n")
time.sleep(3)


data=[]
for disp in range(2):
    print('Downloading data from display #'+str(disp))
    netSock.send("DSPN?"+str(disp)+"\n")
    time.sleep(0.5)
    netSock.send("++read eoi\n") #Change to listening mode
    numPoints=int(gpib.gpibGetData(netSock, 1024, '\004', 1))
    while 1:
        netSock.send("ACTD"+str(disp)+"\n") #Set active display
        time.sleep(1)
        netSock.send("DUMP\n") #Request data
        time.sleep(0.5)
        netSock.send("++read eoi\n") #Change to listening mode
        time.sleep(0.5)
        netSock.send("DFMT?\n") #Request data
        time.sleep(0.5)
        netSock.send("++read eoi\n") #Change to listening mode

        # Receive data
        dumpedData = gpib.gpibGetData(netSock, 1024, '\004', 1)        
        dumpedData = dumpedData[0:len(dumpedData)-2]
        

        #Parse data
        #Matching to the second column of dumpedData
        parsedData=re.findall(r'^[^\s]*\s*([^\s]*)', dumpedData, re.M)
        if len(parsedData) == numPoints :
            data.append(parsedData)
            break
        
        print('retry')
        time.sleep(1)

    
#Construct frequency list (first column of dumpedData)
freq=re.findall(r'^[^\s]*', dumpedData, re.M)
        
#Change the measurement to Norm Variance
print('Switching to normalized variance display')
netSock.send('MEAS0,44\n')  #Norm Variance Ch1
time.sleep(0.1)
netSock.send('MEAS1,45\n')  #Norm Variance Ch2
time.sleep(2)
# Download Normarized Variance
for disp in range(2):
    print('Downloading data from display #'+str(disp))

    netSock.send("DSPN?"+str(disp)+"\n")
    time.sleep(0.5)
    netSock.send("++read eoi\n") #Change to listening mode
    numPoints=int(gpib.gpibGetData(netSock, 1024, '\004', 1))
    while 1:

        netSock.send("ACTD"+str(disp)+"\n") #Set active display
        time.sleep(0.5)
        netSock.send("DUMP\n") #Request data
        time.sleep(0.5)
        netSock.send("++read eoi\n") #Change to listening mode
        netSock.send("DFMT?\n") #Request data
        time.sleep(0.5)
        netSock.send("++read eoi\n") #Change to listening mode

        # Receive data
        dumpedData = gpib.gpibGetData(netSock, 1024, '\004', 1)
        dumpedData = dumpedData[0:len(dumpedData)-2]


        #Parse data
        #Matching to the second column of dumpedData
        parsedData=re.findall(r'^[^\s]*\s*([^\s]*)', dumpedData, re.M)
        if len(parsedData) == numPoints :
            data.append(parsedData)
            break
        
        print('retry')
        time.sleep(1)

#Construct frequency list (first column of dumpedData)
freq=re.findall(r'^[^\s]*', dumpedData, re.M)

print(len(freq))
print(len(data[0]))
print(len(data[1]))
print(len(data[2]))
print(len(data[3]))
    
#Write to the data file
print('Writing data into the data file ...\n')
    
for i in range(len(freq)):
    dataFile.write(freq[i])
    for disp in range(4):
        dataFile.write(' '+data[disp][i])
    dataFile.write('\n')

#Write to parameter file
print('Writing measurement parameters into the parameter file ...\n')
paramFile.write('#SR785 Transfer function measurement\n')
paramFile.write('#Column format = [Freq  Mag(dB) Phase(deg) NormVar1  NormVar2]\n')
paramFile.write('Start frequency = '+options.startFreq+'\n')
paramFile.write('Stop frequency = '+options.stopFreq+'\n')
paramFile.write('Number of frequency points = '+str(options.numOfPoints)+'\n')
paramFile.write('Excitation amplitude = '+options.excAmp+'\n')
paramFile.write('Settling cycles = '+str(options.settleCycles)+'\n')
paramFile.write('Integration cycles = '+str(options.intCycles)+'\n')

netSock.send("++eot_enable 0\n")

dataFile.close()
paramFile.close()
netSock.close()


