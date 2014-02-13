#!/usr/bin/python

import math
import pylab as pl
import sys

from numpy import arange

from scipy.ndimage import gaussian_filter as smooth

from gather import *

if(len(sys.argv) < 2):
  print "Usage: hsShow hemisphereFileName"
  sys.exit()

norm = 1.0 / (4.0 * math.pi)

dSet = 0
dLvl = 0

if(len(sys.argv) == 3):
  dSet = int(sys.argv[2])

if(len(sys.argv) == 4):
  dLvl = int(sys.argv[3])

hs = gth_hemisphere() 
hs.load(sys.argv[1])

hs.divideBySolidAngle()

lPlots = []
iPlots = []
cPlots = []

figLine = pl.figure(1, [5.0 * len(sys.argv[2:]), 1.25 * hs.nLevels])

## Plot the line plots
##
##

pl.figure(1)

thE = np.linspace(0,90,hs.resTheta)

plotId = 1
for levelId in range(hs.nLevels):

  for setId in sys.argv[2:]:

    aArr = hs.toArray(set=int(setId), lvl=2)
    hsArray = hs.toArray(set=int(setId), lvl=levelId)

    relI = hsArray.mean(1) / aArr.mean(1)

    lPlots.append(pl.subplot(hs.nLevels, len(sys.argv[2:]), plotId))
    pl.plot(thE, relI)
    #pl.ylim(0.02, 2.5)

    if levelId < hs.nLevels-1:
      pl.xticks([])
    else:
      pl.xticks([0, 45, 90])

    #pl.yticks(np.linspace(relI.min(), relI.max(), 4))

    if (levelId % 2 == 0):
      lI = "a"
    else:
      lI = "b"

    pl.text(0.95, 0.8, hs.Elements[levelId/2]+" K "+lI, va = "top", ha = "right", transform = lPlots[plotId-1].transAxes)

    plotId += 1

pl.subplots_adjust(left=0.05, right=0.95, hspace=0.05, wspace=0.05)


pl.show()
