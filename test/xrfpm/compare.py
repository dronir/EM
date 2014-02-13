#!/usr/bin/python

import Scientific.IO.NetCDF as nc
import numpy                as np
import pylab                as pl
import sys

if len(sys.argv) != 2 or '-h' in sys.argv or '--help' in sys.argv:
    print "\nUsage:\tcompare.py <file.nc>\n"
    print "Author:\tHannu Parviainen"
    print "\thannu.p.parviainen@helsinki.fi\n"
    sys.exit(1)

filename = sys.argv[1]

f = nc.NetCDFFile(filename, 'r') 

if f.result_type == 'hemisphere':
    ds = np.array(f.variables['res__fld_hemisphere_simulated'].getValue(), np.float64)

elif f.result_type == 'tabulated':
    ds = np.array(f.variables['res__fld_table_simulated'].getValue(), np.float64)


print '\nSimulation type: %s' %f.simulation_type
print 'Result type:     %s' %f.result_type

print ds.mean(), ds.max(), ds.min()

pl.gray()

nLines  = ds.shape[0]/2 - 1
elNames = f.mat__element_names.split()

pIdx = 1
for ln in range(1,nLines + 1):
    sll = ds[0,:,:] / ds[ln*2,:,:]
    vmi = sll.min()
    vma = sll.max()

    pl.subplot(nLines,1,pIdx)
    pl.imshow(sll) #, vmin=vmi, vmax=vma)
    #pl.imshow(ds[ln*2, :, :]) #, vmin=vmi, vmax=vma)

    pl.text(0.5, 1.5, elNames[0] + '/' + elNames[ln])
    pIdx +=1

pl.show()
f.close
