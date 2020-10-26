#!/usr/bin/env python
"""
FSS_watcher.py: monitors the C1:PSL-FSS_FAST channel; if it crosses a certain THRESHOLD (in absolute value), it launches the getMCerr4395.sh script

Dependencies: -

Usage: "FSS_watcher.py THRESHOLD POLLING_TIME"

Diego Bersanetti; Jan 22 2015
"""


from __future__ import division
import subprocess
import sys
import time
import numpy as np

sys.path.append('/opt/rtcds/cdsutils/release/lib')
sys.path.append('/ligo/apps/ubuntu12/nds2-client/lib/python2.7/site-packages/')

import nds2
from ezca import Ezca


# Brief 'usage' function
def usage():
        msg = """\
Usage: FSS_watcher.py
    Example usage:
        FSS_watcher.py THRESHOLD POLLING_TIME
            monitors the C1:PSL-FSS_FAST channel; if it crosses a certain THRESHOLD (in absolute value), it launches the getMCerr4395.sh script
"""
        print >> sys.stderr, msg


# Check input arguments
if len(sys.argv) != 3 or sys.argv[1] == '-h':
	usage()
	sys.exit()
else:
    try:
        FSS_FAST_threshold = float(sys.argv[1])
        polling_time = float(sys.argv[2])
    except (ValueError, TypeError):
        print "ERROR: the input arguments must be numeric!"
        usage()
        sys.exit()


# Initialize ezca
ez = Ezca(ifo=None)

# Connect to fb
conn = nds2.connection('fb',8088)
    
    
while True:
    
    FSS_FAST_current = ez.read('C1:PSL-FSS_FAST')
    print "C1:PSL-FSS_FAST = " + str(FSS_FAST_current)
    
    if np.abs(FSS_FAST_current) >= np.abs(FSS_FAST_threshold):
        print "---> Channel value crossed threshold: executing getMCerr4395.sh script...\n"
        try:
            get_FSS_data = subprocess.Popen(["./getMCerr4395.sh"])
        except:
            print "ERROR: the getMCerr4395 script was not found or working properly... Are you in the right folder?"
        get_FSS_data.terminate()
        print "---> Execution completed! Waiting for " + str(polling_time) + " s ...\n"
    
    time.sleep(polling_time)
