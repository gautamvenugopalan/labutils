'''
A script to interact with a marconi via netgpib
'''
import numpy as np
import sys
sys.path.append('/users/gautam/GIT/labutils/netgpibdata/')
sys.path.reverse()
import netgpib
import cdsutils as cds
import scipy.signal as sig
import os
import ezca 
import argparse

class Marconi:
    '''
    Wrapper with some common commands to a Marconi via GPIB
    '''
    def __init__(self, ip='192.168.113.105', gpibAdd=10, **kwargs):
        print('Establishing connection with {}:{}'.format(ip, gpibAdd))
        try:
            self.gpibObj = netgpib.netGPIB(ip, gpibAdd)
            print('Established connection with {}:{}'.format(ip, gpibAdd))
        except Exception as e:
            print(e)
            print('Unable to connect to {}:{}.'.format(ip,gpibAdd))
    def RFset(self, RFlvl):
        self.gpibObj.command(f'RFLV:VALUE {RFlvl}')
        return()

    def CarrSet(self, CFRQ):
        self.gpibObj.command(f'CFRQ:VALUE {CFRQ}')
        return()

    def outputON(self):
        self.gpibObj.command('RFLV:ON')
        return()

    def outputOFF(self):
        self.gpibObj.command('RFLV:OFF')
        return()

    def CarrGet(self):
        queryCF = self.gpibObj.query('CFRQ?')
        CF = float(queryCF.split(';')[0].split(' ')[-1])
        print('Marconi frequency is {} MHz'.format(CF/1e6))
        return()

    def RFget(self):
        queryLvl = self.gpibObj.query('RFLV?')
        Lvl = float(queryLvl.split(';')[0].split(' ')[-1])
        print('Marconi RF level is {} dBm'.format(Lvl))
        return()