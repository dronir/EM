#!/usr/bin/env python
# coding: utf8
#
# Tool for plotting the packing density of a medium generated
# by the medgen program.
#
# Uses either the Scientific package or the pynetcdf package
# for reading the NetCDF file and matplotlib for plotting.
# These two packages both implement the same NetCDFFile class.
#
# Olli Wilkman Jan 2012


try:
	from Scientific.IO.NetCDF import NetCDFFile
except ImportError:
	try:
		from pynetcdf import NetCDFFile
	except ImportError:
		print "ImportError: Neither Scientific.IO.NetCDF nor pynetcdf was found"

from matplotlib import pyplot as plot
from numpy import array
from sys import argv

infilename = argv[1]

try:
	infile = NetCDFFile(infilename, "r")
except IOError:
	print "Failed to read file %s!" % (infilename)
	exit()

data = infile.variables["medUniform_004_001_001_rho"]

plot.plot(data[1,:], data[0,:])
plot.show()