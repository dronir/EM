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

hs = xrHemisphere() 
hs.load(sys.argv[1])

#hs.divideBySolidAngle()


lPlots = []
iPlots = []
cPlots = []

figLine = pl.figure(1, [5.0 * len(sys.argv[2:]), 1.25 * hs.nLevels])
figImg  = pl.figure(2, [5.0 * len(sys.argv[2:]), 1.25 * hs.nLevels])
figCnt  = pl.figure(3, [5.0 * len(sys.argv[2:]), 1.25 * hs.nLevels])

## Plot the line plots
##
##

pl.figure(0)

pl.subplot(5,1,1)
#pl.semilogy(hs.Spectrum)
pl.plot(hs.Spectrum)

#pl.imshow(hs.muAbsCDF.transpose(), aspect='auto', interpolation='nearest')

pl.subplot(5,1,2)
pl.plot(hs.SpectrumCdf)

pl.subplot(5,1,3)
pl.plot(hs.SpectrumCdfInv)

#pl.subplot(5,1,4)
#pl.hist(hs.SpectrumTest, 100)


pl.subplot(5,1,5)
pl.plot(hs.energy, hs.muAbs.sum(1))
for l in hs.muAbs.transpose():
  pl.plot(hs.energy, l,'--', c='b')

#pl.show()
#sys.exit()

pl.figure(1)

thE = np.linspace(0,90,hs.resTheta)

plotId = 1
for levelId in range(0,hs.nLevels,2):
  for setId in sys.argv[2:]:

    Ka = hs.toArray(set=int(setId), lvl=levelId)
    Kb = hs.toArray(set=int(setId), lvl=levelId+1)

    lPlots.append(pl.subplot(hs.nLevels/2, len(sys.argv[2:]), plotId))
    pl.plot(thE, Ka.mean(1),'--',c='b')
    pl.plot(thE, Kb.mean(1),'--',c='b')
    pl.plot(thE, Ka.mean(1)+Kb.mean(1),c='b')

    #pl.ylim(0.0, 40.0e-2)

    if levelId < hs.nLevels-1:
      pl.xticks([])
    else:
      pl.xticks([0, 45, 90])

    if setId == sys.argv[2]:
      pl.yticks(np.linspace(0, 10, 3))
    else:
      pl.yticks([])

    if (levelId % 2 == 0):
      lI = "a"
    else:
      lI = "b"

    pl.text(0.95, 0.8, hs.Elements[levelId/2]+" K "+lI, va = "top", ha = "right", transform = lPlots[plotId-1].transAxes)

    plotId += 1


pl.subplots_adjust(left=0.10, right=0.95, bottom=0.05, top=0.95, hspace=0.05, wspace=0.05)

## Plot the image plots
##
##

pl.figure(2)
pl.gray()

plotId = 1
for levelId in range(hs.nLevels):
  for setId in sys.argv[2:]: #

    hsArray = smooth(hs.toArray(set=int(setId), lvl=levelId), 1.0)

    iPlots.append(pl.subplot(hs.nLevels, len(sys.argv[2:]), plotId))

    pl.imshow(hsArray, \
              extent=[0,360,90,0], \
              interpolation='nearest', \
              #vmin=0.0, vmax=1.5e-2, \
              aspect = 'auto')

    if levelId < hs.nLevels-1:
      pl.xticks([])
    else:
      pl.xticks([0, 180, 360])

    if setId == sys.argv[2]:
      pl.yticks([15, 45, 75])
    else:
      pl.yticks([])

    plotId += 1

pl.subplots_adjust(left=0.05, right=0.95, hspace=0.05, wspace=0.05)

## Plot the contour plots
##
##

pl.figure(3)
pl.gray()

#v = arange(0.0, 0.55, 0.05) * norm #[0.0, 0.1*norm, 0.2*norm, 0.3*norm]
v = np.linspace(0.0, 2.5, 10)

plotId = 1
for levelId in range(hs.nLevels):
  for setId in sys.argv[2:]: 

    hsArray = smooth(hs.toArray(set=int(setId), lvl=levelId), 1.0)

    print hsArray.min(), hsArray.max()

    cPlots.append(pl.subplot(hs.nLevels, len(sys.argv[2:]), plotId))

    pl.contourf(hsArray, v, extent=[0,360,90,00])
    pl.contour(hsArray, v, extent=[0,360,90,00], colors='black')

    if levelId < hs.nLevels-1:
      pl.xticks([])
    else:
      pl.xticks([0, 180, 360])

    if setId == sys.argv[2]:
      pl.yticks([75, 45, 15],[15, 45, 75])
    else:
      pl.yticks([])

    plotId += 1

pl.subplots_adjust(left=0.05, right=0.95, hspace=0.05, wspace=0.05)
pl.show()
