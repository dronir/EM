#!/usr/bin/python

import pylab as pl
import numpy as np

from numpy import abs, floor, sqrt, cos, arccos, sin, arcsin, pi
from numpy.random import rand

from math import radians as rad, degrees as deg

from gather import *


TWO_PI      = 2.0*pi

def randCartHemiSphUniform(d,n):

    d[:,2] = rand(n)

    phi  = rand(n) * 2.0 * pi
    r    = sqrt(1.0 - d[:,2]**2)

    d[:,0] = r * cos(phi)
    d[:,1] = r * sin(phi)


theta_in = rad(0.0)
resTheta = 45

nRays  = 500000
albedo = 1.0

I = 1.0/float(nRays)
L = gth_hemisphere(resTheta,1,1)

rD    = np.zeros([nRays,3])

randCartHemiSphUniform(rD,nRays)

for nI in range(nRays):
    L.addDataCar(rD[nI,:], I)

L.divideBySolidAngle()

La = L.toArray() * 2.0 * pi

print L.nCells
print La.min()
print La.max()

print 1.0 / sqrt(nRays/float(L.nCells))

pl.gray()
pl.imshow(La, vmin=0.5, vmax=1.5, interpolation='nearest')
pl.show()
