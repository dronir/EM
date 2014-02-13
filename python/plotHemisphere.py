#!/usr/bin/python
import numpy as np
import Scientific.IO.NetCDF as nc
import pylab as pl

dFile = nc.NetCDFFile("vScatter.nc", "r")

for hName in dFile.variables.keys():
    print dFile.variables[hName].getValue()
#  rho = dFile.variables[mName].porosity	
#  h   = dFile.variables[mName].porosityHeight	

#  print np.mean(rho[50:150]), np.std(rho[50:150])

#  pl.plot(h,rho)

#pl.show()