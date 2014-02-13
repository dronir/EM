#!/usr/bin/python2.4

import Scientific.IO.NetCDF as nc
import numpy                as np
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

f.close
