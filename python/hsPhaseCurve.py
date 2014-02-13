#!/usr/bin/python

import math
import sys

import pylab as pl
import numpy as np

from gather import *

from math import radians as rad
from math import degrees as deg

def normalize(v):
    return v / np.sqrt(np.dot(v,v))

if(len(sys.argv) < 2):
  print "Usage: hsPhaseCurve hemisphereFileName"
  sys.exit()

hs = gth_hemisphere()
hs.load(sys.argv[1])

nSamples = 10000
dw       = 2.0 * math.pi / float(nSamples)

## Generate nSamples random sample directions with uniform distribution over a sphere.
##
##
samples      = np.zeros([nSamples,3], dtype='float64')
samples[:,2] = np.random.rand(nSamples)

phi          = np.random.rand(nSamples) * 2.0 * math.pi
r            = np.sqrt(1.0 - samples[:,2]**2)

samples[:,0] = r * np.cos(phi)
samples[:,1] = r * np.sin(phi)


## Compute the cos(theta) between the lighsource and a sample point
##
##
I      = np.array([0.0, 0.0, 1.0],dtype='float64')
iCosTh = np.dot(samples, I)
iTheta = np.arccos(iCosTh)
lit    = iCosTh > 0.0

p      = np.zeros([3,180],dtype='float64')

Ip   = np.zeros([nSamples,3], dtype='float64')
Ep   = np.zeros([nSamples,3], dtype='float64')
ePhi = np.zeros([nSamples], dtype='float64')

for a in range(0,90):
    E = np.array([np.sin(rad(a)), 0.0, np.cos(rad(a))],dtype='float64')

    ## Compute the cos(theta) between the viewer and a sample point
    ##
    ##
    eCosTh = np.dot(samples, E)
    eTheta = np.arccos(eCosTh)
    vis    = eCosTh > 0.0

    ok = np.logical_and(vis, lit)

    ## Compute the azimuth angles.
    ##
    
    for i in range(nSamples):
      Ip[i,:] =  samples[i,:] * iCosTh[i]
      Ep[i,:] =  samples[i,:] * eCosTh[i]

      Ip[i] = normalize(I - Ip[i,:])
      Ep[i] = normalize(E - Ep[i,:])
      
      ePhi[i] = np.arccos(np.dot(Ip[i,:], Ep[i,:]) - 1e-8)

    iTh = iTheta[ok]
    eTh = eTheta[ok]
    ePh = ePhi[ok]

    for i in range(len(iTh)):
      p[1,a] += hs.eval(deg(iTh[i]), deg(eTh[i]), deg(ePh[i]))

    p[1,a] *= dw / 15.0

    p[0,a] = a
    p[2,a] = np.sum(iTh*eTh / (iTh+eTh)) * dw

    

#for i in range(nSamples):
#  pl.plot([samples[i,0]], [samples[i,1]],'o', c=str(ePhi[i] / math.pi))

pl.plot(p[0,:],p[1,:],'o')
pl.plot(p[0,:],p[2,:],'o')
pl.show()
