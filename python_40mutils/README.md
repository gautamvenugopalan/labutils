# Filter Module State Extractor (python3)
filterModuleStateExtractor3.py reads in foton files on the LHO CDS machines.

## Testing
```bash
$ bash test_script_fmse3.sh
```
should run a test call to the filterModuleStateExtractor3.py using the scipyenv anaconda environement on CDS machines.

If you are on LHO CDS machines, you can call

## Instructions:
```bash
$ source /opt/rtcds/userapps/release/cds/h1/scripts/setup_anaconda
$ conda activate scipyenv
$ cd /ligo/home/craig.cahillane/Git/labutils/python_40mutils
$ python filterModuleStateExtractor3.py --help
```
will print the help and arguments that filterModuleStateExtractor3.py accepts.




## Examples:
- plot the April 1 2021 ASC-MICH_P and ASC-MICH_Y filter module settings
```bash
$ python filterModuleStateExtractor3.py --plotTF --debug 1301354613 ASC-MICH_P ASC-MICH_Y
```

### WARNING: 
The filter module names (e.g. ASC-MICH_P) *must* be the final arguments.
`argparse` will read in all the final positional arguments as filter module names.
All optional flags must be before the positional.


## Help:
```bash
usage: filterModuleStateExtractor3.py [-h] [--debug] [--hostServer HOSTSERVER]
                                      [--portNumber PORTNUMBER]
                                      [--fotonDictLocation FOTONDICTLOCATION]
                                      [--fotonArchiveLocation FOTONARCHIVELOCATION]
                                      [--plotTF] [--saveTF] [--savePlot]
                                      [--saveDir SAVEDIR] [--fflow FFLOW]
                                      [--ffhigh FFHIGH]
                                      [--ffnumPoints FFNUMPOINTS]
                                      gpstime ...

Craig Cahillane, July 12, 2018 filterModuleStateExtractor.py takes in a
gpstime and a list of filter modules and returns their full state. Also delves
into the foton filter archives, finds the correct foton file given your
gpstime, and returns the second-order seconds of the foton file. Example for
knowing the state of filter module ASC-MICH_P and ASC-MICH_Y today, July 12,
2018 19:52:06 UTC: python filterModuleStateExtractor.py 1215460344 ASC-MICH_P
ASC-MICH_Y {FilterModuleName}_SWSTAT contains the state of all the filter
module buttons. There are 16 bits (2 bytes) for every SWSTAT. I call the
furthest RIGHT bit, the least significant bit, bit 0, and the furthest left
bit, the most significant bit, bit 15. Bit 0 = FM1 State Bit 1 = FM2 State ...
Bit 9 = FM10 State Bit 10 = INPUT on/off Bit 11 = OFFSET on/off Bit 12 =
OUTPUT on/off Bit 13 = LIMIT on/off Bit 14 = HOLD OUTPUT on/off Bit 15 =
DECIMATION on/off

positional arguments:
  gpstime               GPS time to retrieve the filter module state.
  filter_modules        Filter Modules to check the state of

optional arguments:
  -h, --help            show this help message and exit
  --debug, -d           Flag. If set, prints verbose debug messages.
  --hostServer HOSTSERVER, -s HOSTSERVER
                        String. Name of the server you want to try to access
                        data from. Default is nds.ligo-wa.caltech.edu. Also
                        automatically tries h1nds1 at 8088 if this server
                        fails.
  --portNumber PORTNUMBER, -p PORTNUMBER
                        Int. Port number of server selected. Default is 31200.
  --fotonDictLocation FOTONDICTLOCATION
                        String corresponding to the location of the foton
                        lookup dictionary. Default is your home directory.
  --fotonArchiveLocation FOTONARCHIVELOCATION
                        String corresponding to the location of the foton
                        lookup dictionary. Default is your home directory.
  --plotTF              Flag. If set, interactively plots your filter module
                        TFs. Need this set to save anything.
  --saveTF              Flag. If set, saves your filter module TF data
  --savePlot            Flag. If set, saves your plots
  --saveDir SAVEDIR     String. Directory to save your data and plots. Will
                        create a (saveDir)/data and (saveDir)/plots directory
                        for you if they do not already exist. Default is
                        /ligo/home/craig.cahillane/fotonTFs
  --fflow FFLOW, -ffl FFLOW
                        Int. Log 10 low frequency. Default is -2.
  --ffhigh FFHIGH, -ffh FFHIGH
                        Int. Log 10 high frequency. Default is 5.
  --ffnumPoints FFNUMPOINTS, -ffn FFNUMPOINTS
                        Int. Number of points in freq vec. Default is 1000.