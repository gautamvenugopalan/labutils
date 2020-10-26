#!/usr/bin/env python

# This creates a spectrum using an AG4395A analyzer.
# If you want to plot the data, tell it to with --plot. It opens a Matlab session to plot.
# Also make sure the AG4395A.py code prints the channel name in the data file as
# "dataFile.write('%Channel '+str(ch)+'\n')", otherwise Matlab can't read the data file.
# Plotting dual channels currently doesn't work since channel 2's data is in the
# same column as channel 1's data. I don't know python well enough to put them in
# separate columns :(
#
# SVN: $Id: SPAG4395A 42648 2013-05-22 00:32:14Z controls $
# Michael Rodruck
#
# modified by Rana   30-Oct-2013

import os
import re
import sys
import math
import optparse
import time
import pdb
import netgpib
import AG4395A
import termstatus
import numpy as np
import matplotlib.pyplot as plt


# Parse options
usage = """usage: %prog [options]

This program executes a spectral measurement using AG4395A.
Various measurement conditions can be set using options.
The measurement result will be saved in FILENAME.dat and the measurement parameters in FILENAME.par.
Optionally, it can plot the retrieved data.
"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-f", "--file", dest="filename",
                  help="Output file name without an extension", default="data")

parser.add_option("-i", "--ip",
                  dest="ipAddress", default="192.168.113.108",
                  help="IP address/Host name")

parser.add_option("-a", "--address",
                  dest="gpibAddress", type="int", default=10,
                  help="GPIB device address")

parser.add_option("-d", "--dualchannel",
                  dest="dual", default = False,
                  action="store_true",
                  help="Set to dual channel mode.")

#parser.add_option("--ch", "--channel",
#                  dest="chan", type="str", default="R",
#                  help="What channel are you looking at (R, A, B)")

parser.add_option("-R", "--channelR",
                  dest="chanR", default=False,
                  action="store_true",
                  help="What channel are you looking at (R, A, B)? Specify with '-R'")

parser.add_option("-A", "--channelA",
                  dest="chanA",default=False,
                  action="store_true",
                  help="What channel are you looking at (R, A, B)? Specify with '-A'")

parser.add_option("-B", "--channelB",
                  dest="chanB",default=False,
                  action="store_true",
                  help="What channel are you looking at (R, A, B)? Specify with '-B'")

parser.add_option("-v", "--averaging",
                   dest="numAvg", type="str",default="10",
                  help="Number of averages")

parser.add_option("--att", "--attenuation",
                  dest="dBatt", type="str",default="10",
                  help="dB of attenuation")

parser.add_option("--st", "--start",
                  dest="startF", type="str",default="1MHz",
                  help="Start frequency")

parser.add_option("--fin", "--end",
                  dest="endF", type="str",default="200MHz",
                  help="End frequency")

parser.add_option("--bw", "--Bandwidth",
                  dest="BW", type="str",default="30kHz",
                  help="Bandwidth. Changes number of points - max number of points is 801")

parser.add_option("--plot",
                  dest="plotData", default=False,
                  action="store_true",
                  help="Plot the downloaded data.")
parser.add_option("--lin",
                  dest="lin_p", default=None,
                  action="store_true",
                  help="Plot with linear axis")
parser.add_option("--semix",
                  dest="semix", default=None,
                  action="store_true",
                  help="Plot with logarithmic x axis")
parser.add_option("--semiy",
                  dest="semiy_p", default=None,
                  action="store_true",
                  help="Plot with logarithmic y axis")
parser.add_option("--loglog",
                  dest="loglog_p", default=None,
                  action="store_true",
                  help="Plot with both logarithmic axis")
parser.add_option("--title",
                  dest="title", type="string",default="",
                  help="Title of the measurement. The given string will be written into the parameter file.")



(options, args) = parser.parse_args()
# Create a netGPIB class object
print('Connecting to '+str(options.ipAddress)+' ...'),
gpibObj=netgpib.netGPIB(options.ipAddress, options.gpibAddress, '\004',0)
print('done.')

#File names
dataFileName  = options.filename + '.dat'
paramFileName = options.filename + '.par'

print('Data will be written into '       + dataFileName)
print('Parameters will be written into ' + paramFileName)
print('Setting up parameters for the measurement') 

from time import gmtime

gpibObj.command("PRES")
time.sleep(0.1)
gpibObj.command("SA")
time.sleep(0.1)

    
if options.dual:
    duac = "ON"
    if options.chanR:
        MS = "R"
        gpibObj.command("CHAN1")
        gpibObj.command("MEAS "+str(MS))
        if options.chanA:
            MS = "A"
            gpibObj.command("CHAN2")
            gpibObj.command("MEAS "+str(MS))
            if options.chanB:
                MS = "B"
                gpibObj.command("CHAN2")
                gpibObj.command("MEAS "+str(MS))
    elif options.chanB:
        MS = "B"
        gpibObj.command("CHAN1")
        gpibObj.command("MEAS "+str(MS))
        if options.chanA:
            MS = "A"
            gpibObj.command("CHAN2")
            gpibObj.command("MEAS "+str(MS))
    gpibObj.command("DUAC "+str(duac))
    time.sleep(0.1)
    gpibObj.command("CHAN1")
    time.sleep(0.1)
    gpibObj.command("AVERFACT "+str(options.numAvg))
    time.sleep(0.1)
    gpibObj.command("AVER OFF")
    gpibObj.command("STAR "+str(options.startF))
    time.sleep(0.1)
    gpibObj.command("STOP "+str(options.endF))
    time.sleep(0.1)
    gpibObj.command("BW "+str(options.BW))
    time.sleep(0.1)
    gpibObj.command("CHAN2")
    time.sleep(0.1)
    gpibObj.command("AVERFACT "+str(options.numAvg))
    time.sleep(0.1)
    gpibObj.command("AVER OFF")
    gpibObj.command("STAR "+str(options.startF))
    time.sleep(0.1)
    gpibObj.command("STOP " + str(options.endF))
    time.sleep(0.1)
    gpibObj.command("BW " + str(options.BW))
    time.sleep(1)
    print('Parameters set')
    print('try to open data file')
    dataFile = open(dataFileName,'w')
    print('1')
    paramFile = open(paramFileName,'w')
    print('2')
    gpibObj.command("AVER ON")
    print('3')
    time.sleep(0.1)
    gpibObj.command("CHAN2")
    time.sleep(0.1)
    gpibObj.command("AVER ON")
    tim = gpibObj.query("SWET?")
    tot_time = float(tim)*int(options.numAvg)+int(options.numAvg)
    print('Run time will be approximately '+str(tot_time)+'s')    
    gpibObj.command("AVERREST") #Start measurement
    gpibObj.command("CHAN1")
    gpibObj.command("AVERREST")
    print('Running...')
    time.sleep((tot_time))
    a = gmtime()
    date_time = '%02d/%02d, %02d:%02d:%02d' %(a[1],a[2],a[3],a[4],a[5])
    print('Done')
else:
    if options.chanR:
        MS = "R"
    elif options.chanA == 1:
        MS = "A"
    elif options.chanB == 1:
        MS = "B"
    gpibObj.command("MEAS "+str(MS))
    time.sleep(0.1)
    gpibObj.command("AVERFACT "+str(options.numAvg))
    time.sleep(0.1)
    gpibObj.command("AVER OFF")
    gpibObj.command("STAR "+str(options.startF))
    time.sleep(0.1)
    gpibObj.command("STOP "+str(options.endF))
    time.sleep(0.1)
    gpibObj.command("ATTAUTO OFF")
    time.sleep(0.1)
    gpibObj.command("ATT" +str(MS) +" " +str(options.dBatt) +"DB")
    #gpibObj.command("ATTB 10DB")
    time.sleep(0.11)
    gpibObj.command("BW "+str(options.BW))
    time.sleep(5)
    print('Parameters set')
    gpibObj.command("AVER ON")
    dataFile = open(dataFileName,'w')
    paramFile = open(paramFileName,'w')
    tim = gpibObj.query("SWET?")
    tot_time = float(tim)*int(options.numAvg)+int(options.numAvg)
    print('Run time is '+str(round(tot_time))+' s')
    gpibObj.command("AVERREST") # Start measurement
    print('Running...')
    time.sleep((tot_time))
    a = gmtime()
    date_time = '%02d/%02d, %02d:%02d:%02d' %(a[1],a[2],a[3],a[4],a[5])
    print('Done')


print('Getting data (if this takes more than ~15s it probably crashed, try again)')
AG4395A.getdata(gpibObj, dataFile, paramFile)
time.sleep(0.1)
#Deal with an empty title
if options.title == "":
    options.title = options.filename
time.sleep(0.1)
print('Print parameters')
#Parameters
paramFile.write('Title: '+options.title+'\n')
paramFile.write('############## Spectral Measurement Parameters #########################\n')
fSpan = gpibObj.query("SPAN?")

paramFile.write('Frequency Span: '+fSpan+'\n')

fStart = gpibObj.query("STAR?")
paramFile.write('Start Frequency: '+fStart+'\n')

fStop = gpibObj.query("STOP?")
paramFile.write('Stop Frequency: '+fStop+'\n')

fRes = int(gpibObj.query("POIN?"))
paramFile.write('Number of Points: '+str(fRes)+'\n')

nAvg = int(gpibObj.query("AVERFACT?"))
paramFile.write('Number of Averages: '+str(nAvg)+'\n')

AAtt = int(gpibObj.query("ATTAUTO?"))
paramFile.write('Auto Attenuation: '+str(AAtt)+'\n')


paramFile.write('####################################################################\n')

AG4395A.getparam(gpibObj, options.filename, dataFile, paramFile)

dataFile.close()
paramFile.close()
gpibObj.close()
print('Done!')
if options.lin_p:
    plot_type = "plot"
if options.semix:
    plot_type = "semilogx"
if options.semiy_p:
    plot_type = "semilogy"
if options.loglog_p:
    plot_type = "loglog"
    
#os.system('matlab -r "load '+dataFileName+', loglog('+options.filename+'(:,1),'+options.filename+'(:,2))"')
if options.plotData:
    #os.system('matlab -r "load '+dataFileName+', '+plot_type+'('+options.filename+'(:,1),'+options.filename+'(:,2)), xlabel(\'Frequency (Hz)\'), ylabel(\'Power\'), title({[\'\\bf{AG4395A Spectrum at '+date_time+' UTC}\'];[\'\\rm{Start Frequency: '+options.startF+', Stop Frequency: '+options.endF+'}\'];[\'Number of Averages: '+options.numAvg+', Bandwidth: '+options.BW+'\']})"')    
    #import gpibplot
    #gpibplot.plotAG4395A(options.filename, xlog=options.xlog, ylog=options.ylog)
    #raw_input('Press enter to quit:')
    data = np.loadtxt(dataFileName, delimiter=',', skiprows=1)
    f    = data[:,0]
    # convert from Watts to Volts
    v    = np.sqrt(data[:,1]*50)

    plt.figure(9)
    if   plot_type == "plot":
        plt.plot(f/1e6,     v, color='m', linewidth=3)
    elif plot_type == "semilogx":
        plt.semilogx(f/1e6, v, color='m', linewidth=3)
    elif plot_type == "semilogy":
        plt.semilogy(f/1e6, v, color='m', linewidth=3)
    elif plot_type == "loglog":
        plt.loglog(f/1e6,   v, color='m', linewidth=3)


    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Input Signal [V]")
    plt.grid()
    plt.savefig('out.pdf', transparent=True, bbox_inches='tight')
    #plt.show()
    #time.sleep(3)
    #plt.close()
    raw_input('Press enter to quit:')
