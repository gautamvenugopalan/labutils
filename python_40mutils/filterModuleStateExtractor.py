'''
Craig Cahillane, July 12, 2018
filterModuleStateExtractor.py takes in a gpstime and a list of filter modules and returns their full state.  Also delves into the foton filter archives, finds the correct foton file given your gpstime, and returns the second-order seconds of the foton file.

Example for knowing the state of filter module ASC-MICH_P and ASC-MICH_Y today, July 12, 2018 19:52:06 UTC:

python filterModuleStateExtractor.py 1215460344 ASC-MICH_P ASC-MICH_Y

{FilterModuleName}_SWSTAT contains the state of all the filter module buttons.  There are 16 bits (2 bytes) for every SWSTAT.
I call the furthest RIGHT bit, the least significant bit, bit 0, and the furthest left bit, the most significant bit, bit 15.

Bit 0 = FM1 State
Bit 1 = FM2 State
...
Bit 9 = FM10 State
Bit 10 = INPUT on/off
Bit 11 = OFFSET on/off
Bit 12 = OUTPUT on/off
Bit 13 = LIMIT on/off
Bit 14 = HOLD OUTPUT on/off
Bit 15 = DECIMATION on/off
'''
from __future__ import division
import numpy as np
import os, sys, time
import argparse
import nds2
import subprocess
import readFotonFilterFile
import pickle
from bisect import bisect
import scipy.signal as sig
import pylab

class bcolors: # For printing to console with color
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def parseArgs():
    '''
    Parses user shell arguments input to this script
    '''
    parser = argparse.ArgumentParser(description=__doc__)
    
    parser.add_argument('gpstime', help='GPS time to retrieve the filter module state.')
    parser.add_argument('filter_modules', nargs=argparse.REMAINDER, help='Filter Modules to check the state of')
    parser.add_argument('--debug', '-d', action='store_true', help='Flag. If set, prints verbose debug messages.')
    parser.add_argument('--hostServer', '-s', type=str, default='nds.ligo-wa.caltech.edu', help='String. Name of the server you want to try to access data from.  Default is nds.ligo-wa.caltech.edu.  Also automatically tries h1nds1 at 8088 if this server fails.')
    parser.add_argument('--portNumber', '-p', type=int, default=31200, help='Int. Port number of server selected.  Default is 31200.')
    parser.add_argument('--fotonDictLocation', type=str, default=os.path.expanduser('~'), help='String corresponding to the location of the foton lookup dictionary.  Default is your home directory.')
    parser.add_argument('--fotonArchiveLocation', type=str, default='/opt/rtcds/lho/h1/chans/filter_archive', help='String corresponding to the location of the foton lookup dictionary.  Default is your home directory.')
    parser.add_argument('--plotTF', action='store_true', help='Flag. If set, interactively plots your filter module TFs')
    parser.add_argument('--saveTF', action='store_true', help='Flag. If set, saves your filter module TF data')
    parser.add_argument('--savePlot', action='store_true', help='Flag. If set, saves your plots')
    parser.add_argument('--saveDir', type=str, default='{0}/{1}'.format(os.path.expanduser('~'), 'fotonTFs'), help='String.  Directory to save your data and plots.  Will create a (saveDir)/data and (saveDir)/plots directory for you if they do not already exist.  Default is {0}'.format('{0}/{1}'.format(os.path.expanduser('~'), 'fotonTFs')))
    args = parser.parse_args()

    if args.gpstime == 'now':
        args.gpstime = tconvertNow()
    else:
        args.gpstime = int(args.gpstime)

    for ii, FM in enumerate(args.filter_modules):
        if '--' in FM:
            print 'Deleting {0}'.format(FM)
            args.filter_modules = np.delete(args.filter_modules, ii)

    if args.debug:
        print 'Arguments:'
        print 'gpstime = {0}'.format(args.gpstime)
        print 'filter_modules = {0}'.format(args.filter_modules)
        print 'debug = {0}'.format(args.debug)
        print 'hostServer = {0}'.format(args.hostServer)
        print 'portNumber = {0}'.format(args.portNumber)
    return args

def tconvertNow():
    '''Calls tconvert from shell'''
    return int(subprocess.check_output(['tconvert', 'now']))

def extractState(args, filter_module):
    '''
    Recover the state of filter_module at gpstime
    '''
    gpstime = args.gpstime
    hostServer = args.hostServer
    portNumber = args.portNumber
    #hostServer = 'h1nds1'
    #portNumber = 8088
    FMname = '{0}{1}'.format('H1:', filter_module)    

    allChans = np.array(['{0}{1}'.format(FMname, '_SWSTAT'),
                         '{0}{1}'.format(FMname, '_GAIN'), 
                         '{0}{1}'.format(FMname, '_OFFSET'),
                         '{0}{1}'.format(FMname, '_TRAMP'),
                         '{0}{1}'.format(FMname, '_LIMIT')])
    
    chanData = np.array([])
    try:
        conn = nds2.connection(hostServer, portNumber) 
        chanData = conn.fetch(gpstime-1, gpstime, allChans)
        conn.close()
    except RuntimeError:
        print 
        print 'Failed to retrieve filter module named {0} on nds.ligo-wa.caltech.edu'.format(filter_module)
        print 'Trying h1nds1...'
        try:
            conn.close()
        except NameError:
            print 'Variable "conn" was not defined'
        try: 
            conn2 = nds2.connection('h1nds1', 8088)
            chanData = conn2.fetch(gpstime, gpstime+1, allChans)
            conn2.close()
        except RuntimeError:
            print 'Failed to retrieve {0} again'.format(filter_module)
            print 'Wait two seconds and try again...'
            time.sleep(2)
            try: 
                conn2 = nds2.connection('h1nds1', 8088)
                chanData = conn2.fetch(gpstime, gpstime+1, allChans)
                conn2.close()
            except RuntimeError:
                print 'Failed to retrieve {0} for a third time'.format(filter_module)
                print 'Should {0} have a dash instead?'.format(filter_module)
                print 'Exiting...'
                sys.exit()

    return chanData

def onOff(inp):
    '''
    Returns 'ON' if inp == 1 or inp == '1'
    Else if inp == 0, returns 'OFF'
    Else returns not sure
    '''
    inp = int(inp)
    if inp == 1:
        return '{0}ON{1}'.format(bcolors.GREEN, bcolors.ENDC)
    elif inp == 0:
        return '{0}OFF{1}'.format(bcolors.RED, bcolors.ENDC)
    else:
        return 'Not Sure'

def query_yes_no(question, default='no'):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None: 
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        print question,
        sys.stdout.write(prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid: 
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
    return


def storeState(FMchanData):
    '''Put the state into a convenient dict'''
    SWSTAT = FMchanData[0].data[-1]
    GAIN = FMchanData[1].data[-1]
    OFFSET = FMchanData[2].data[-1]
    TRAMP = FMchanData[3].data[-1]
    LIMIT = FMchanData[4].data[-1]
    SWSTAT = np.binary_repr(SWSTAT)

    FMdict = {'ONOFF': {'INPUT': SWSTAT[5],
                        'OUTPUT': SWSTAT[3],
                        'OFFSET': SWSTAT[4],
                        'HOLD OUTPUT': SWSTAT[1],
                        'DECIMATION': SWSTAT[0],
                        'LIMIT': SWSTAT[2]
                       },
              'VALUES': {'SWSTAT': SWSTAT,
                         'OFFSET': OFFSET,
                         'TRAMP': TRAMP,
                         'LIMIT': LIMIT,
                         'GAIN': GAIN
                        }
             }
    for ii in range(1,11):
        FMdict['ONOFF']['FM{0}'.format(ii)] = SWSTAT[-ii] 
    return FMdict

def createFotonLookupDict(args, fotonDictName='fotonLookupTable.dict'):
    ''' Creates a fast lookup dictionary for filter modules, with keys being filter module names and values being directories containing those FMs.  Should only need to be run once '''
    fullFotonDictName = '{0}/{1}'.format(args.fotonDictLocation,fotonDictName)
    if os.path.isfile(fullFotonDictName):
        print 'Foton dictionary lookup table already exists at {0}'.format(fullFotonDictName)
        return fullFotonDictName
        #answer = query_yes_no('Okay to overwrite?', default='yes')
        #if not answer:
        #    print 'NO'
        #    print 'Not overwriting {0}'.format(fullFotonDictName)
        #    return

    lookupDict = {}
    fotonList = os.listdir(args.fotonArchiveLocation)
    for fotonSubdir in fotonList:
        print fotonSubdir
        newDir = '{0}/{1}'.format(args.fotonArchiveLocation, fotonSubdir)
        curFotonFileList = os.listdir(newDir) # get foton files in specific dir for all time
        # Find the most recent foton file
        gpsFotonFiles = np.array([])
        for fotonFile in curFotonFileList:
            curSplit = fotonFile.split('_')
            if len(curSplit) == 2:  # if there's only one underscore
                try:
                    curGpstime = int(curSplit[1].split('.')[0]) # remove .txt
                except ValueError: # in case the second element is not an int
                    continue
                gpsFotonFiles = np.append(gpsFotonFiles, curGpstime)
        if len(gpsFotonFiles) == 0: # if no gpstimes in foton file names, continue
            continue
        sortedGpsFotonFiles = np.sort(gpsFotonFiles)
        latestFotonGpstime = sortedGpsFotonFiles[-1]
        latestFotonFile = '{0}_{1}{2}'.format(curSplit[0], int(latestFotonGpstime), '.txt')
        #print 'latestFotonFile = {0}'.format(latestFotonFile)

        # Now that latest file is found, read in this file using readFotonFilterFile.py
        curFotonDict = readFotonFilterFile.readFilterFile('{0}/{1}'.format(newDir, latestFotonFile))
        curModules = curFotonDict.keys() # filter modules names are the keys
        for curModule in curModules: # not very pythonic but it works
            lookupDict[curModule] = fotonSubdir

    # Dump the dictionary
    with open(fullFotonDictName, 'wb') as filename:
        pickle.dump(lookupDict, filename)
    return fullFotonDictName

def findFotonFilename(args, filter_module, fotonDictName='fotonLookupTable.dict'):
    ''' Finds the filterModule foton file you need based on gpstime and filter module name'''
    # read in lookup dictionary
    fullFotonDictName = '{0}/{1}'.format(args.fotonDictLocation, fotonDictName)
    with open(fullFotonDictName, 'rb') as filename:
        fotonDict = pickle.load(filename)
    subDir = fotonDict[filter_module]
    fullDir = '{0}/{1}'.format(args.fotonArchiveLocation, subDir)
    curFotonFileList = os.listdir(fullDir)
    # Sort foton files by gpstime
    gpsFotonFiles = np.array([])
    for fotonFile in curFotonFileList:
        curSplit = fotonFile.split('_')
        if len(curSplit) == 2:  # if there's only one underscore
            try:
                curGpstime = int(curSplit[1].split('.')[0]) # remove .txt
            except ValueError:
                continue
            gpsFotonFiles = np.append(gpsFotonFiles, curGpstime)
    sortedGpsFotonFiles = np.sort(gpsFotonFiles) 
    # Find foton file corresponding to user gpstime
    index = bisect(sortedGpsFotonFiles, args.gpstime) # bisect returns index between two values
    if index is not 0:
        index = index - 1 # if not first index (unlikely), to go past index
    correctGpstime = int(sortedGpsFotonFiles[index])
    correctFotonFile = '{fullDir}/{0}_{1}{2}'.format(curSplit[0], correctGpstime, '.txt', fullDir=fullDir)
    if os.path.isfile(correctFotonFile):
        return correctFotonFile
    else:
        print '{0} is not a file'.format(correctFotonFile)
        return

def readFotonFile(fotonFilename):
    '''Reads a fotonFile using readFotonFilterFile utils'''    
    fotonDict = readFotonFilterFile.readFilterFile(fotonFilename)
    return fotonDict

def createFilterModuleTF(args, FMdicts, requestedFotonDict, ff=np.logspace(-2, np.log10(2**13), 300), nyquist=2.0**14, saveTF=False, plotTF=True, savePlot=False):
    '''Now that we have the foton dictionary and the filter module state at the gpstime of interest, we create a transfer function for the whole filter module relevant at that time.'''
    for filter_module in args.filter_modules:
        FMdict = FMdicts[filter_module]
        fDict = requestedFotonDict[filter_module]
        GAIN = FMdict['VALUES']['GAIN']
        INPUT = FMdict['ONOFF']['INPUT']
        OUTPUT = FMdict['ONOFF']['OUTPUT']

        wwNorm = 2*np.pi*ff/nyquist
        totalTF = np.array([])
        totalFMs = np.array([])
        for FM in fDict.keys(): # recall that FM1 == 0, FM2 == 1... in the keys of this dict
            print 'Getting FM{0}'.format(FM+1)
            if int(FMdict['ONOFF']['FM{0}'.format(FM+1)]) == 0:
                print 'FM{0} is OFF, skipping...'.format(FM+1)
                continue
            else: 
                totalFMs = np.append(totalFMs, 'FM{0}'.format(FM+1))
            sos = fDict[FM]['sosCoeffs']
            ww, TF = sig.sosfreqz(sos, worN=wwNorm)
            if not np.count_nonzero > 0:
                print 'FM{0} is just a FM of gain 0'.format(FM+1)
                print 'Not including in total transfer function calculation'
                continue
            if len(totalTF) == 0:
                totalTF = np.copy(TF)
            else:
                totalTF *= TF
        if len(totalTF) == 0:
            print 'No FMs active, TF = 1 for all frequencies'
            totalTF = np.ones(len(ff))
        if not GAIN == 0:
            totalTF *= GAIN
        else:
            print '{cb}WARNING: {0} FM GAIN = 0{ce}'.format(filter_module, cb=bcolors.RED, ce=bcolors.ENDC)

        filename = '{0}_TF_GPStime_{1}'.format(filter_module.replace('-','_'), args.gpstime)
        if args.saveTF:
            saveData = np.vstack((ff, np.abs(totalTF), np.angle(totalTF))).T
            saveDirData = '{0}/{1}'.format(args.saveDir,'data')
            if not os.path.isdir(saveDirData):
                os.mkdir(saveDirData)
            np.savetxt(saveDirData+'/'+filename+'.txt', saveData, header='{0:18s}{1:18s}{2:18s}'.format('Frequency [Hz]', 'Magnitude', 'Phase [rads]'))

        pylab.ion()
        fig, axes = pylab.subplots(2, sharex=True)
        fig.suptitle('{0} Total TF at GPStime = {1}\nFMs Engaged: {2}, Gain = {3:.3f}, Input = {4}, Output = {5}'.format(filter_module, args.gpstime, ' '.join(totalFMs), GAIN, INPUT, OUTPUT))
        axes[0].loglog(ff, np.abs(totalTF))
        axes[1].semilogx(ff, 180/np.pi*np.angle(totalTF))

        axes[0].set_ylabel('Magnitude')
        axes[1].set_ylabel('Phase [degs]')
        axes[1].set_xlabel('Frequency [Hz]')

        axes[0].set_xlim([min(ff), max(ff)])
        axes[1].set_xlim([min(ff), max(ff)])
        minYAxis2 = 45.0*np.floor(min(180/np.pi*np.angle(totalTF))/45.0)
        maxYAxis2 = 45.0*np.ceil(max(180/np.pi*np.angle(totalTF))/45.0)
        axes[1].set_ylim([minYAxis2, maxYAxis2])
        axes[1].set_yticks(np.linspace(minYAxis2, maxYAxis2, np.round((maxYAxis2-minYAxis2)/45)+1))

        axes[0].grid()
        axes[0].grid(which='minor', ls='--')
        axes[1].grid()
        axes[1].grid(which='minor', ls='--')
        if args.savePlot:
            saveDirPlots = '{0}/{1}'.format(args.saveDir,'plots')
            if not os.path.isdir(saveDirPlots):
                os.mkdir(saveDirPlots)
            pylab.savefig(saveDirPlots+'/'+filename+'.pdf', bbox_inches='tight')
    pylab.show(block=True)
    return ff, totalTF

def printState(FMchanData):
    SWSTAT = FMchanData[0].data[-1]
    GAIN = FMchanData[1].data[-1]
    OFFSET = FMchanData[2].data[-1]
    TRAMP = FMchanData[3].data[-1]
    LIMIT = FMchanData[4].data[-1]
    SWSTAT = np.binary_repr(SWSTAT)
    
    print
    print bcolors.YELLOW+'=============================================='+bcolors.ENDC
    print bcolors.PURPLE+'  {0}   (at GPS Time = {1})'.format(FMchanData[0].name[3:-7], FMchanData[0].gps_seconds)+bcolors.PURPLE
    print bcolors.YELLOW+'=============================================='+bcolors.ENDC
    print 'SWSTAT = {0}'.format(SWSTAT)
    print
    print 'INPUT is {0}'.format(onOff(SWSTAT[5]))
    print 'OUTPUT is {0}'.format(onOff(SWSTAT[3]))
    print 'OFFSET is {0}, OFFSET VALUE = {1}'.format(onOff(SWSTAT[4]), OFFSET)
    print
    print 'HOLD OUTPUT is {0}'.format(onOff(SWSTAT[1]))
    print 'DECIMATION is {0}'.format(onOff(SWSTAT[0]))
    print 'LIMIT is {0}, LIMIT VALUE = {1}'.format(onOff(SWSTAT[2]), LIMIT) 
    print
    print '{1}GAIN = {0}{2}'.format(GAIN, bcolors.BOLD, bcolors.ENDC)
    print 'TRAMP = {0}'.format(TRAMP)
    print
    for ii in range(1,6):
        print 'FM{0}  '.format(ii),
    print
    for ii in range(1,6):
        print '{0:14}'.format(onOff(SWSTAT[-ii])),
    print
    print '----------------------------'
    for ii in range(6,11):
        print '{0:14}'.format(onOff(SWSTAT[-ii])),
    print
    for ii in range(6,11):
        print 'FM{0}  '.format(ii),
    print
    print 
    return

def getFMState(filter_modules, gpstime, debug=False, hostServer='nds.ligo-wa.caltech.edu', portNumber=31200, fotonDictLocation=os.path.expanduser('~'), fotonArchiveLocation='/opt/rtcds/lho/h1/chans/filter_archive', plotTF=True, saveTF=False, savePlot=False, saveDir='{0}/{1}'.format(os.path.expanduser('~'), 'fotonTFs')):
    '''
    Main function of filterModuleStateExtractor library
    Acquires the filter module at the gpstime from user inputs, 
    and returns two dictionaries and a namespace corresponding to the filter module information.
    OUTPUTS:
    The first is the info from the medm screen, like what buttons are on, etc.
    The second is the foton files containing the second order sections of the filters.    
    The third is the args namespace corresponding to your inputs.

    You have to create a lookup dictionary using createFotonLookupDict() in this library to get this function to work.

    Example Usages:
    getFMState(['IMC-MCL'], 'now')
    returns a current foton dict for the IMC-MCL filter module    

    getFMState(['ASC-INP1_P', 'ASC-INP1_Y'], 1210012345)
    returns the foton dicts for INP1 pitch and yaw at GPS = 1210012345

    Enter filter module names with a dash first, then underscores.  

    Requires tconvert on your machine if you want to use 'now' for gpstimes, 
    and a version of scipy that has signal.sosfreqz() if you want to make transfer functions.
    '''
    # Create an artificial args class
    class Test:
        def __init__(self, FMs, gps, db, hs, pn, fdl, fal, pTF, sTF, sp, sd):
            self.gpstime = gps
            self.filter_modules = FMs
            self.debug = db
            self.hostServer = hs
            self.portNumber = pn
            self.fotonDictLocation = fdl
            self.fotonArchiveLocation = fal
            self.plotTF = pTF
            self.saveTF = sTF
            self.savePlot = sp
            self.saveDir = sd
    
    args = Test(filter_modules, gpstime, debug, hostServer, portNumber, fotonDictLocation, fotonArchiveLocation, plotTF, saveTF, savePlot, saveDir)
    # Do same argument conditioning as in parseArgs()
    if args.gpstime == 'now':
        args.gpstime = tconvertNow()
    else:
        args.gpstime = int(args.gpstime)

    for ii, FM in enumerate(args.filter_modules):
        if '--' in FM: 
            print 'Deleting {0}'.format(FM)
            args.filter_modules = np.delete(args.filter_modules, ii) 

    if args.debug:
        print 'Arguments:'
        print 'gpstime = {0}'.format(args.gpstime)
        print 'filter_modules = {0}'.format(args.filter_modules)
        print 'debug = {0}'.format(args.debug)
        print 'hostServer = {0}'.format(args.hostServer)
        print 'portNumber = {0}'.format(args.portNumber)
        print 'fotonDictLocation = {0}'.format(args.fotonDictLocation)
        print 'fotonArchiveLocation = {0}'.format(args.fotonArchiveLocation)
        print 'plotTF = {0}'.format(args.plotTF)
        print 'saveTF = {0}'.format(args.saveTF)
        print 'savePlot = {0}'.format(args.savePlot)
        print 'saveDir = {0}'.format(args.saveDir)

    # If foton lookup dict does not exist, create it
    createFotonLookupDict(args, fotonDictName='fotonLookupTable.dict')

    FMdicts = {}
    requestedFotonDict = {} # make a smaller dict with only the requested FMs
    for filter_module in args.filter_modules:
        FMchanData = extractState(args, filter_module)
        if not FMchanData:
            print 'Channels for {0} not acquired'.format(filter_module)
            print 'Skipping {0}'.format(filter_module)
        else:
            FMdict = storeState(FMchanData)
            FMdicts[filter_module] = FMdict
            # Get foton file, read relevant sos, make TF
            # Foton names everything with underscores, no dashes :(
            filter_module_underscore = filter_module.replace('-', '_')
            fotonFilename = findFotonFilename(args, filter_module_underscore)
            fotonDict = readFotonFile(fotonFilename)
            curDict = fotonDict[filter_module_underscore] # switch name back to dashes
            requestedFotonDict[filter_module] = curDict

    if args.plotTF:
        createFilterModuleTF(args, FMdicts, requestedFotonDict)

    return FMdicts, requestedFotonDict, args
#########################################################################################
if __name__=='__main__':
    startTime = time.time()
    args = parseArgs()
    
    FMdicts, requestedFotonDict, args = getFMState(args.filter_modules, args.gpstime, args.debug, args.hostServer, args.portNumber, args.fotonDictLocation, args.fotonArchiveLocation)

    print
    print 'Done in {0} minutes'.format((time.time() - startTime)/60.0)
    print
