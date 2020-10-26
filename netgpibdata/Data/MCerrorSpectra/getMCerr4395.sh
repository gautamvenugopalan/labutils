#!/usr/bin/env bash

TNOW=`date +%g%m%d_%H%M%S`

echo "Getting Low Frequency Data"
echo
../../SPAG4395A.py -f ${TNOW}LF -i vanna.martian -B -v 5 --att=0 --start=10kHz --end=1MHz --bw=2kHz


echo
echo "Getting High Frequency Data"
echo
../../SPAG4395A.py -f ${TNOW}HF -i vanna.martian -B -v 5 --att=0 --start=1MHz --end=11MHz --bw=2kHz
echo

