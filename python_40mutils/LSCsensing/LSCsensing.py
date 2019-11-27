import argparse
import LSCdemodUtils as utils
from LSCdemodUtils import *
from cycler import cycler

parser = argparse.ArgumentParser(description=
        '''Usage:
            python LSCsensing.py --paramFile <path_to_parameter_file>
           will use the parameter to do the following operations:
           1. Download the data for the times, and from the NDS server specified in the param file, downsample it as specified,and save to a .hdf5 file in the Data directory.
           2. Demodulate the LSC sensing photodiode signals at the known frequencies at which the DoFs were driven.
           3. Infer the physical motion of the drive using the known actuator calibration.
           4. Construct a sensing matrix, which is the map from the DoF motion to the response in some sensors.
        ''')
parser.add_argument('--paramFile', type=str, help='Path to the parameter file specifying which optic to analyze', nargs=1, required=True)
args = parser.parse_args()

# Global setup
customCycler = (cycler(color=['#30a2da','#fc4f30','#e5ae38','#6d904f','#8b8b8b']))
plt.rc('axes', prop_cycle=customCycler)
paramFile = args.paramFile[0]
par = importParams(paramFile)
if str(par['tStart']) not in os.listdir(globDataDir):
    logging.debug('Making subdirectory for {} in {}'.format(par['tStart'], globDataDir))
    os.mkdir(globDataDir+str(par['tStart']))
if str(par['tStart']) not in os.listdir(globFigDir):
    logging.debug('Making subdirectory for {} in {}'.format(par['tStart'], globFigDir))
    os.mkdir(globFigDir+str(par['tStart']))
# Define the directory Macros
utils.dataDir = globDataDir + str(par['tStart']) + '/'
utils.figDir = globFigDir + str(par['tStart']) + '/'
# Build the FFT dictionary
fftParams = {}
fftParams['window'] = ('tukey',0.25)
fftParams['tFFT'] = 1024


# Download the data
logging.debug('Downloading data..')
dlData(paramFile)
# Demodulate the lines
demodData(paramFile)
# Make the sensing matrix
plotData(paramFile, saveFig=True)

# Finally, print the matrix that should go into the EPICS screen
#logging.info('Computed matrix that will best diagonalize the sensor signals is')
