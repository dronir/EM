#!/usr/bin/python

import pylab as pl
import sys

from math   import radians as rad
from math   import degrees as deg
from math   import pi

from gather import *

test = gth_hemisphere() 

for fName in sys.argv[1:]:

  test.load(fName)

  print fName
  
  a = test.toArray()
  pl.plot(a.sum(1))
  pl.show()


  thSlice1 = test.thetaSlice(math.pi * 0.00, 0)
  thSlice2 = test.thetaSlice(math.pi * 0.00, 25)
  thSlice3 = test.thetaSlice(math.pi * 0.00, 45)
  thSlice4 = test.thetaSlice(math.pi * 0.00, 75)

  phSlice1 = test.phiSlice(rad(15.0), 15)
  phSlice2 = test.phiSlice(rad(25.0), 25)
  phSlice3 = test.phiSlice(rad(45.0), 45)

  slice1 = test.thetaSlice(math.pi * 0.00, 0)
  slice2 = test.thetaSlice(math.pi * 0.00, 1)
  slice3 = test.thetaSlice(math.pi * 0.00, 2)
  slice4 = test.thetaSlice(math.pi * 0.00, 3)

  #slice1 = test.phiSlice(rad(15.0), 15)
  #slice2 = test.phiSlice(rad(25.0), 25)
  #slice3 = test.phiSlice(rad(45.0), 45)

  pl.plot(slice1[:,0]*(180.0/pi),slice1[:,1], label='i = 15')
  pl.plot(slice2[:,0]*(180.0/pi),slice2[:,1], label='i = 25')
  pl.plot(slice3[:,0]*(180.0/pi),slice3[:,1], label='i = 45')
  pl.plot(slice4[:,0]*(180.0/pi),slice4[:,1], label='i = 75')

  #pl.plot([65.0],[test.eval(45, 45.0, 65.0)/ 15.0], 'o')

#pl.ylim([0.0, 0.3])
#pl.xlim([0.0, 360.0])
#pl.xlim([0.0, 0.5 * math.pi])

pl.title(sys.argv[1].split('/')[-1].split('.')[0])
pl.legend()

pl.show()
