#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

dataFileName = "131030_0309.dat"

data = np.loadtxt(dataFileName,delimiter=',',skiprows=1)
f = data[:,0]
# convert from Watts to Volts
v = np.sqrt(data[:,1]*50)
l1 = plt.semilogy(f, v, color='m', linewidth=3)
lx = plt.xlabel("Frequency [MHz]")
ly = plt.ylabel("Input Signal [V]")
plt.show()

