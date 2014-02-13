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
    da = np.array(f.variables['res__fld_hemisphere_analytical'].getValue(), np.float64)

elif f.result_type == 'tabulated':

    ds = np.array(f.variables['res__fld_table_simulated'].getValue(), np.float64)
    da = np.array(f.variables['res__fld_table_analytical'].getValue(), np.float64)


de = np.abs((da - ds) / da) * 100.0

print '\nSimulation type: %s' %f.simulation_type
print 'Result type:     %s' %f.result_type

print '\nError:\tmean: \t%6.3f\n\tmax: \t%6.3f\n\tmin: \t%6.3f\n' %(de.mean(), de.max(), de.min())


pl.gray()

nLines  = ds.shape[0]/2 - 1
elNames = f.mat__element_names.split()

pIdx = 1
for ln in range(1,nLines + 1):
    sll = ds[0,:,:] / ds[ln*2,:,:]
    all = da[0,:,:] / da[ln*2,:,:]

    vmi = np.min(sll.min(), all.min())
    vma = np.max(sll.max(), all.max())

    pl.subplot(nLines,3,pIdx)
    if(ln==1): pl.title('Simulation')
    pl.imshow(sll, vmin=vmi, vmax=vma)
    pl.text(0.5, 1.5, elNames[0] + '/' + elNames[ln])
    pIdx +=1

    pl.subplot(nLines,3,pIdx)
    if(ln==1): pl.title('Analytical')
    pl.imshow(all, vmin=vmi, vmax=vma)
    pIdx +=1

    pl.subplot(nLines,3,pIdx)
    if(ln==1): pl.title('Abs(Sim - An)')
    pl.imshow(np.abs(sll - all), vmin=0.0, vmax=0.1*vma)
    pIdx +=1

pl.show()
f.close
