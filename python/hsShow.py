#!/usr/bin/env python

import math
import pylab as pl
import sys

import Image
import ImageFilter

from numpy import arange

from gather import *

if(len(sys.argv) < 2):
  print "Usage: hsShow hemisphereFileName"
  sys.exit()

norm = 1.0 / (2.0 * 4.0 * math.pi)

dSet = 0
dLvl = 0

if(len(sys.argv) == 3):
  dSet = int(sys.argv[2])

if(len(sys.argv) == 4):
  dLvl = int(sys.argv[3])

hs = gth_hemisphere() 
hs.load(sys.argv[1])

#hs.divideBySolidAngle()

plots = []

pl.figure(1, [5.0 * len(sys.argv[2:]), 1.25 * hs.nLevels])

pl.gray()

v = arange(0.0, 0.55, 0.05) * norm #[0.0, 0.1*norm, 0.2*norm, 0.3*norm]

plotId = 1
for setId in sys.argv[2:]:

  hsArray = hs.toArray(set=int(setId), lvl=0)

  for levelId in range(1,hs.nLevels):
    hsArray += hs.toArray(set=int(setId), lvl=levelId)

  plots.append(pl.subplot(hs.nLevels, len(sys.argv[2:]), plotId))
  pl.imshow(hsArray, extent=[0,360,90,0], interpolation='nearest')#, vmin=0.0, vmax=norm*hs.w)
  plotId += 1

pl.subplots_adjust(left=0.05, right=0.95, hspace=0.0, wspace=0.05)
pl.show()

sys.exit()

plotId = 1
for levelId in range(hs.nLevels):
  for setId in sys.argv[2:]: #
  #for setId in range(hs.nLevels):
    hsArray = hs.toArray(set=int(setId), lvl=levelId)

    print hsArray.min(), hsArray.max()

    plots.append(pl.subplot(hs.nLevels, len(sys.argv[2:]), plotId))

    #pl.imshow(hsArray, \
    #          extent=[0,360,90,0], \
    #          interpolation='nearest', \
    #          norm = colors.Normalize(vmin = 0.0, vmax = 1.0/(8.0*math.pi), \
    #          clip = False))

    #pl.contourf(hsArray, v, extent=[0,360,90,0])

    pl.imshow(hsArray, extent=[0,360,90,0], interpolation='nearest', vmin=0.0, vmax=norm*hs.w)

    plotId += 1

#plots.append(pl.subplot(hs.nLevels + 1, len(sys.argv[2:]), plotId))

#pl.imshow(hs.toArray(set=int(setId), lvl=0) + hs.toArray(set=int(setId), lvl=1), \
#            extent=[0,360,90,0], \
#            interpolation='nearest', \
#            norm = colors.Normalize(vmin = 0.0, vmax = 1.0/(8.0*math.pi), \
#            clip = False))

pl.subplots_adjust(left=0.05, right=0.95, hspace=0.0, wspace=0.05)
pl.show()
