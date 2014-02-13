#!/usr/bin/python

import Scientific.IO.NetCDF as nc
import numpy as np
import sys
import math
import pylab as pl
import matplotlib.colors as colors

from numpy import floor, sqrt, sin, cos, arccos, arctan2, pi

fName = sys.argv[1]

try:
  dFile = nc.NetCDFFile(fName, "r")
except IOError:
  print "Error reading file, exiting."
  sys.exit()

if "fluorescenceData" not in dFile.variables.keys():
  print "Error: not a proper fluorescent file."
  sys.exit()

if "Elements" in  dir(dFile):
  elements = str(dFile.Elements).split()

data           = np.array(dFile.variables['fluorescenceData'].getValue())
data_an        = np.array(dFile.variables['fluorescenceData_analytic'].getValue())

spectrum       = np.array(dFile.variables['Spectrum'].getValue())
#spectrumCDF    = np.array(dFile.variables['SpectrumCdf'].getValue())
#spectrumCDFInv = np.array(dFile.variables['SpectrumCdfInv'].getValue())

#spcE           = np.array(dFile.SpectrumEnergy)

#testSpc        = np.array(dFile.TestSpectrum)

e = dFile.Elements.split()

iSi = e.index('Si')*2
iCa = e.index('Ca')*2
iTi = e.index('Ti')*2

dFile.close()

#for i,l in enumerate(data):
  #print elements[i/2], l
#  pl.plot(l/data[8,:])

Na = data[2,:]+data[3,:]
Mg = data[4,:]+data[5,:]
Al = data[6,:]+data[7,:]
Si = data[iSi,:]+data[iSi+1,:]
K  = data[12,:]+data[13,:]
Ca = data[iCa,:]+data[iCa+1,:]
Ti = data[iTi,:]+data[iTi+1,:]
Fe = data[20,:]+data[21,:]

Ca_an = data_an[iCa,:]+data_an[iCa+1,:]
Si_an = data_an[iSi,:]+data_an[iSi+1,:]
Ti_an = data_an[iTi,:]+data_an[iTi+1,:]

#FeaCaa = data[18,:]/data[12,:]

#FeCa = Fe/Ca

CaSi = Ca/Si
TiSi = Ti/Si
TiCa = Ti/Ca
SiAl = Si/Al
SiMg = Si/Mg
SiNa = Si/Na

CaSi_an = Ca_an/Si_an
TiSi_an = Ti_an/Si_an
TiCa_an = Ti_an/Ca_an

pl.figure(1)
#pl.plot(FeaCaa / FeaCaa[0])

pl.plot(CaSi / CaSi[0], '-', label='Ca/Si')
pl.plot(TiSi / TiSi[0], '--', label='Ti/Si')
pl.plot(TiCa / TiCa[0], '-.', label='Ti/Ca')

#pl.plot(CaSi_an / CaSi_an[0], '-', label='Ca/Si an')
#pl.plot(TiSi_an / TiSi_an[0], '--', label='Ti/Si an')
#pl.plot(TiCa_an / TiCa_an[0], '-.', label='Ti/Ca an')

#pl.plot(SiNa / SiNa[0], '-', label='Si/Na')
#pl.plot(SiAl / SiAl[0], '--', label='Si/Al')
#pl.plot(SiMg / SiMg[0], '-.', label='Si/Mg')

pl.legend(loc=0)
pl.title('Fe55 source')

#pl.figure(2)
#pl.subplot(1,4,1)
#pl.plot(spcE, spectrum)

#pl.subplot(1,4,2)
#pl.plot(spcE, spectrumCDF)

#pl.subplot(1,4,3)
#pl.plot(spectrumCDFInv)

#pl.subplot(1,4,4)
#pl.plot(testSpc)

print TiSi
print TiCa

#print CaSi/CaSi[0] - TiSi/TiSi[0]

pl.show()

print
print "    %12s %12s %12s" %("Ca/Si", "Ti/Si", "Ti/Ca")
for i in range(Ca.size):
  print "%3i %12.7f %12.7f %12.7f" %(i*5, CaSi[i], TiSi[i], TiCa[i])
