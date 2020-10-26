# This module provides data access to SR785 analyzer

import re
import sys
import math
from optparse import OptionParser
from socket import *
import time
import pdb
import gpib

def getdata(netSock, gpibAddr, dataFile, paramFile):


    # Initialization of the GPIB-Ethernet converter box
    print('Initializing the GPIB-Ethernet converter\n')
    netSock.setblocking(0)
    netSock.send("++addr "+str(gpibAddr)+"\n")
    time.sleep(0.1)
    netSock.send("++eos 3\n")
    time.sleep(0.1)
    netSock.send("++mode 1\n")
    time.sleep(0.1)
    #netSock.send("++clr\n")
    time.sleep(0.1)
    netSock.send("++auto 0\n")
    time.sleep(0.1)
    netSock.send("++ifc\n")
    time.sleep(0.1)
    netSock.send("++read_tmo_ms 4000\n")
    time.sleep(0.1)
    netSock.send("++eot_char 4\n")
    netSock.send("++eot_enable 1\n")
    netSock.send("++addr "+str(gpibAddr)+"\n")
    #pdb.set_trace()
    
    # Get number of displays
    netSock.send("OUTX0\n")
    time.sleep(0.1)
    netSock.send("DFMT?\n")
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    
    numOfDisp=int(gpib.gpibGetData(netSock,100,'\004'))
    if numOfDisp == 3:
       numOfDisp=1        
    numOfDisp = numOfDisp + 1

    # Set dump destination
    netSock.send("++auto 0\n")
    time.sleep(0.1)
    netSock.send("POUT2\n") # ASCII DUMP Mode
    time.sleep(0.1)
    netSock.send("PDST3\n") # Destination = GPIB
    time.sleep(0.1)
    netSock.send("PCIC0\n") # GPIB Control = Host

    netSock.send("++eot_char 4\n")
    time.sleep(0.1)
    netSock.send("++eot_enable 1\n")
    time.sleep(0.1)
    #pdb.set_trace()
    # Download data
    
    data=[]
    freq=[]
    for disp in range(numOfDisp):
        print('Downloading data from display #'+str(disp))
        netSock.send("DSPN?"+str(disp)+"\n")
        time.sleep(0.5)
        netSock.send("++read eoi\n") #Change to listening mode
        numPoints=int(gpib.gpibGetData(netSock, 1024, '\004', 1))
        while 1:

            netSock.send("++auto 0\n")  #Turn off automatic mode
            netSock.send("ACTD"+str(disp)+"\n") #Set active display
            netSock.send("DUMP\n") #Request data
            netSock.send("++read eoi\n") #Change to listening mode
            time.sleep(0.5)
            netSock.send("DFMT?\n") #Request data
            time.sleep(0.5)
            netSock.send("++read eoi\n") #Change to listening mode

            # Receive data
            dumpedData = gpib.gpibGetData(netSock,8192,'\004', 1)
            dumpedData = dumpedData[0:len(dumpedData)-2]

            #Parse data
            #Matching to the second column of dumpedData
            parsedData=re.findall(r'^[^\s]*\s*([^\s]*)', dumpedData, re.M)

            print len(parsedData)
            print numPoints
            if len(parsedData) == numPoints :
            #if True:
                data.append(parsedData)
                #Construct frequency list (first column of dumpedData)
                freq.append(re.findall(r'^[^\s]*', dumpedData, re.M))
                break
        
            print('retry')
            time.sleep(1)
            

    #Check the measurement group
    netSock.send( 'MGRP?'+str(disp)+'\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    
    if i == 0:  #If FFT group
        for disp in range(numOfDisp):
            dataFile.write('Display #'+str(disp+1)+'\n')

            for j in range(len(freq[disp])):
                dataFile.write(freq[disp][j])
                dataFile.write(' '+data[disp][j]+'\n')
        
    else:  #Else
    #Write to the data file
        print('Writing data into the data file ...\n')
    
        for i in range(len(freq[0])):
            dataFile.write(freq[0][i])
        #        for disp in range(numOfDisp):
            for disp in range(numOfDisp):
                dataFile.write(' '+data[disp][i])
            dataFile.write('\n')


    netSock.send("++eot_enable 0\n")

def getparam(netSock, gpibAddr, filename, dataFile, paramFile):
    #Get measurement parameters
    
    print('Reading measurement parameters\n')
    

    # Initialization of the GPIB-Ethernet converter box
    netSock.send("++addr "+str(gpibAddr)+"\n")
    netSock.send("++eos 3\n")
    netSock.send("++eot_enable 1\n")
    netSock.send("++eot_char 4\n")
    netSock.send("++mode 1\n")
    netSock.send("++clr\n")
    netSock.send("++auto 0\n")
    netSock.send("++ifc\n")
    netSock.send("++read_tmo_ms 4000\n")
    netSock.send("++addr "+str(gpibAddr)+"\n")
    
    #Get the display format
    netSock.send( "DFMT?\n")
    netSock.send("++read eoi\n")
    numOfDisp=int(gpib.gpibGetData(netSock, 100, '\004'))
    if numOfDisp == 3:
        numOfDisp=1
    numOfDisp = numOfDisp + 1

    #Get display parameters for each display
    measGrp=[]
    measurement=[]
    view=[]
    unit=[]

    time.sleep(0.1)

    for disp in range(numOfDisp):
        netSock.send( 'MGRP?'+str(disp)+'\n')
        time.sleep(0.1)
        netSock.send("++read eoi\n")
        i=int(gpib.gpibGetData(netSock, 100, '\004'))
        measGrp.append({0: 'FFT' , 
                         1: 'Correlation', 
                         2: 'Octave', 
                         3: 'Swept Sine', 
                         4: 'Order', 
                         5: 'Time/Histogram'}[i])

    #Get measurement
        netSock.send( 'MEAS?'+str(disp)+'\n')
        time.sleep(0.1)
        netSock.send("++read eoi\n")
        i=int(gpib.gpibGetData(netSock, 100, '\004'))
        measurement.append(
        {0: 'FFT 1',
         1: 'FFT 2',
         2: 'Power Spectrum 1',
         3: 'Power Spectrum 2',
         4: 'Time 1',
         5: 'Time 2',
         6: 'Windowed Time 1',
         7: 'Windowed Time 2',
         8: 'Orbit',
         9: 'Coherence',
         10: 'Cross Spectrum',
         11: 'Frequency Response',
         12: 'Capture Buffer 1',
         13: 'Capture Buffer 2',
         14: 'FFT User Function 1',
         15: 'FFT User Function 2',
         16: 'FFT User Function 3',
         17: 'FFT User Function 4',
         18: 'FFT User Function 5',
         19: 'Auto Correlation 1',
         20: 'Auto Correlation 2',
         21: 'Cross Correlation',
         22: 'Time 1',
         23: 'Time 2',
         24: 'Windowed Time 1',
         25: 'Windowed Time 2',
         26: 'Capture Buffer 1',
         27: 'Capture Buffer 2',
         28: 'Correlation Function 1',
         29: 'Correlation Function 2',
         30: 'Correlation Function 3',
         31: 'Correlation Function 4',
         32: 'Correlation Function 5',
         33: 'Octave 1',
         34: 'Octave 2',
         35: 'Capture 1',
         36: 'Capture 2',
         37: 'Octave User Function 1',
         38: 'Octave User Function 2',
         39: 'Octave User Function 3',
         40: 'Octave User Function 4',
         41: 'Octave User Function 5',
         42: 'Spectrum 1',
         43: 'Spectrum 2',
         44: 'Normalized Variance 1',
         45: 'Normalized Variance 2',
         46: 'Cross Spectrum',
         47: 'Frequency Response',
         48: 'Swept Sine User Function 1',
         49: 'Swept Sine User Function 2',
         50: 'Swept Sine User Function 3',
         51: 'Swept Sine User Function 4',
         52: 'Swept Sine User Function 5',
         53: 'Linear Spectrum 1',
         54: 'Linear Spectrum 2',
         55: 'Power Spectrum 1',
         56: 'Power Spectrum 2',
         57: 'Time 1',
         58: 'Time 2',
         59: 'Windowed Time 1',
         60: 'Windowed Time 2',
         61: 'RPM Profile',
         62: 'Orbit',
         63: 'Track 1',
         64: 'Track 2',
         65: 'Capture Buffer 1',
         66: 'Capture Buffer 2',
         67: 'Order User Function 1',
         68: 'Order User Function 2',
         69: 'Order User Function 3',
         70: 'Order User Function 4',
         71: 'Order User Function 5',
         72: 'Histogram 1',
         73: 'Histogram 2',
         74: 'PDF 1',
         75: 'PDF 2',
         76: 'CDF 1',
         77: 'CDF 2',
         78: 'Time 1',
         79: 'Time 2',
         80: 'Capture Buffer 1',
         81: 'Capture Buffer 2',
         82: 'Histogram User Function 1',
         83: 'Histogram User Function 2',
         84: 'Histogram User Function 3',
         85: 'Histogram User Function 4',
         86: 'Histogram User Function 5'
         }[i])

        #View information
        netSock.send( 'VIEW?'+str(disp)+'\n')
        time.sleep(0.1)
        netSock.send("++read eoi\n")
        i=int(gpib.gpibGetData(netSock, 100, '\004'))
        view.append({0: 'Log Magnitude',
                     1: 'Linear Magnitude',
                     2: 'Magnitude Squared',
                     3: 'Real Part',
                     4: 'Imaginary Part',
                     5: 'Phase',
                     6: 'Unwrapped Phase',
                     7: 'Nyquist',
                     8: 'Nichols'}[i])

        #Units
        netSock.send( 'UNIT?'+str(disp)+'\n')
        time.sleep(0.1)
        netSock.send("++read eoi\n")
        result=gpib.gpibGetData(netSock, 512, '\004')
        result=result[0:len(result)-1]  # Chop a new line character
        unit.append(result)   


    #Input Source
    netSock.send( 'ISRC?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    time.sleep(0.1)
    inputSource={0: 'Analog',
                 1: 'Capture'}[i]

    #Input Mode
    netSock.send( 'I1MD?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH1inputMode={0: 'Single ended',
                 1: 'Differential'}[i]
         
    netSock.send( 'I2MD?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH2inputMode={0: 'Single ended',
                 1: 'Differential'}[i]

    #Grounding
    netSock.send( 'I1GD?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH1Grounding={0: 'Float',
                 1: 'Grounded'}[i]
    netSock.send( 'I2GD?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH2Grounding={0: 'Float',
                 1: 'Grounded'}[i]

    #Coupling
    netSock.send( 'I1CP?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH1Coupling={0: 'DC',
                 1: 'AC',
                  2:'ICP'}[i]
    netSock.send( 'I2CP?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH2Coupling={0: 'DC',
                 1: 'AC',
                  2:'ICP'}[i]

    #Input Range
    netSock.send( 'I1RG?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    result=gpib.gpibGetData(netSock, 100, '\004')
    match=re.search(r'^\s*([-+\d]*),.*',result)
    CH1Range=str(float(match.group(1)))
    match=re.search(r'\d,(\d)',result)
    i=int(match.group(1))
    CH1Range=CH1Range+{0: 'dBVpk', 1: 'dBVpp', 2: 'dBVrms', 3: 'Vpk', 4: 'Vpp', 5: 'Vrms', 
                       6: 'dBEUpk', 7: 'dBEUpp', 8: 'dBEUrms', 9: 'EUpk', 10: 'EUpp', 11: 'EUrms'}[i]

    netSock.send( 'I2RG?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    result=gpib.gpibGetData(netSock, 100, '\004')
    match=re.search(r'^\s*([-+\d]*),.*',result)
    CH2Range=str(float(match.group(1)))
    match=re.search(r'\d,(\d)',result)
    i=int(match.group(1))
    CH2Range=CH2Range+{0: 'dBVpk', 1: 'dBVpp', 2: 'dBVrms', 3: 'Vpk', 4: 'Vpp', 5: 'Vrms', 
                       6: 'dBEUpk', 7: 'dBEUpp', 8: 'dBEUrms', 9: 'EUpk', 10: 'EUpp', 11: 'EUrms'}[i]

    #Auto Range
    netSock.send( 'A1RG?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH1AutoRange={0: 'Off', 1: 'On'}[i]
    netSock.send( 'I1AR?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH1AutoRangeMode={0: 'Up Only', 1: 'Tracking'}[i]
    netSock.send( 'A2RG?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH2AutoRange={0: 'Off', 1: 'On'}[i]
    netSock.send( 'I2AR?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH2AutoRangeMode={0: 'Normal', 1: 'Tracking'}[i]

    #Anti-Aliasing Filter
    netSock.send( 'I1AF?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH1AAFilter={0: 'Off', 1: 'On'}[i]
    netSock.send( 'I2AF?\n')
    time.sleep(0.1)
    netSock.send("++read eoi\n")
    i=int(gpib.gpibGetData(netSock, 100, '\004'))
    CH2AAFilter={0: 'Off', 1: 'On'}[i]

    #Write to the parameter file
    #Header
    paramFile.write('#SR785 parameter file\n#This file contains measurement parameters for the data saved in '
                    +filename+'.dat\n')
    #For the ease of getting necessary information for plotting the data, several numbers and strings are put first
    #The format is the number of channels comes first, then the title of the channels and the units follow, 
    #one per line.
    paramFile.write('#The following lines are for a matlab plotting function\n')
    paramFile.write(str(numOfDisp)+'\n')
    for disp in range(numOfDisp):
        paramFile.write(view[disp]+'\n')
        paramFile.write('Frequency\n')
    for disp in range(numOfDisp):
        paramFile.write(unit[disp]+'\n')
        paramFile.write('Hz\n')

    #Human readable section begin
    paramFile.write('##################### Measurement Parameters #############################\n')
    paramFile.write('Measurement Group: ')
    for disp in range(numOfDisp):
        paramFile.write(' "'+measGrp[disp]+'"')
    paramFile.write('\n')
    
    paramFile.write('Measurements: ')
    for disp in range(numOfDisp):
        paramFile.write(' "'+measurement[disp]+'"')
    paramFile.write('\n')
    
    paramFile.write('View: ')
    for disp in range(numOfDisp):
        paramFile.write(' "'+view[disp]+'"')
    paramFile.write('\n')
        
    paramFile.write('Unit: ')
    for disp in range(numOfDisp):
        paramFile.write(' "'+unit[disp]+'"')
    paramFile.write('\n\n')
    
    paramFile.write('##################### Input Parameters #############################\n')

    paramFile.write('Input Source: ')
    paramFile.write(inputSource+'\n')

    paramFile.write('Input Mode: ')
    paramFile.write(CH1inputMode+', '+CH2inputMode+'\n')
    
    paramFile.write('Input Grounding: ')
    paramFile.write(CH1Grounding+', '+CH2Grounding+'\n')
    
    paramFile.write('Input Coupling: ')
    paramFile.write(CH1Coupling+', '+CH2Coupling+'\n')

    paramFile.write('Input Range: ')
    paramFile.write(CH1Range+', '+CH2Range+'\n')

    netSock.send("++eot_enable 0\n")
