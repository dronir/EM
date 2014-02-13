#!/usr/bin/python

import pylab as pl
import sys

from gather import *

test = gth_hemisphere() 

for fName in sys.argv[1:]:

  test.load(fName)

  print fName
  
#  slice1 = test.thetaSlice(math.pi * 0.0)
#  slice2 = test.thetaSlice(math.pi * 1.0)

#  pl.plot(slice1[:,0],slice1[:,1])
#  pl.plot(slice2[:,0],slice2[:,1])

#pl.ylim([0.0,1.1])
#pl.xlim([0.0, 0.5 * math.pi])

#pl.show()