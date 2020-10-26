#!/usr/bin/env bash

TNOW=`date +%y%m%d_%H%M%S`

#./SPAG4395A.py -f Data/MCerrorSpectra/${TNOW} -i vanna.martian -A -v 4 --att=20 --start=1kHz --end=3MHz --bw=1kHz

echo "Dowloading Agilent 4395 data without changing the params."
echo
./netgpibdata -f Data/MCOLG_${TNOW} -i vanna -d AG4395A -a 10
echo


