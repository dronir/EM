#!/usr/bin/python
import sys
import numpy as np
import Scientific.IO.NetCDF as nc
import pylab as pl

if(len(sys.argv) != 2):
  print "Usage: plotRho mediumFileName"
  sys.exit()

fName = sys.argv[1]

try:
  dFile = nc.NetCDFFile(fName, "r")
except IOError:
  print "Error reading file, exiting."
  sys.exit()

for mName in dFile.variables.keys():
   if(mName.find("rho") > 0):
     print mName
     rho = dFile.variables[mName].getValue()[1,:]
     h   = dFile.variables[mName].getValue()[0,:]

     print np.mean(rho[50:150]), np.std(rho[50:150])

     pl.plot(h,rho)

pl.show()