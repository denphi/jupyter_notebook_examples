#!/usr/bin/env python
import sys
import numpy as np
from scipy import interpolate, optimize
from subprocess import call

def dump_hdf5(name, v, desc=''):
    np.set_printoptions(threshold=np.nan)
    print 'HDF5:%s:5FDH' % repr({'name': name, 'desc': desc, 'value':v})

# Calculates the Yield Strength from the
# stress-strain.dat file using the offset yield
# method.
# http://www.engineeringarchives.com/les_mom_offsetyieldmethod.html

call(['../slip-gb'] + sys.argv[1:])

# data = np.genfromtxt('stress-strain.dat')
with open('stress-strain.dat') as f:
    data = np.array([[float(a) for a in line.split()] for line in f])

x = data[:,0]
y = data[:,1]
numpts=500
f = interpolate.interp1d(x, y, kind='linear')
xx = np.linspace(np.min(x), np.max(x), numpts)
xline = np.array([np.min(x), np.max(x)])

# calculate slope using points with strain <= 0.02
sx = x[x<=.02]
sy = y[x<=.02]
A = np.array([sx, np.ones(len(sx))]).T
slope = np.linalg.lstsq(A, sy)[0][0]
offset = .002
linef = lambda x: (x-offset)*slope
yline = linef(xline)

yint = np.inf
g = np.array([(s, f(s)-linef(s)) for s in xx])
for i in range(0, len(g)):
    if i < numpts-2 and g[i][1]*g[i+1][1] <= 0:
        xint = optimize.brentq(lambda x: f(x)-linef(x), g[i][0], g[i+1][0])
        yint = float(f(xint))

dump_hdf5('stress', yint, "Yield Stress (Pa)")

