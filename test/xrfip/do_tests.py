#!/usr/bin/python
import numpy as np
import sys
import os

from string import Template

xrfBin          = '../../build_32/bin/xrfip'
nThreads        = 4

verbose         = '.True.'

parTemplateFile = open('parTemplate.txt')
parTemplate     = Template(parTemplateFile.read())

outputTypes     = ['hemisphere', 'tabulated']
methods         = ['analytic', 'first_order', 'mc_force', 'mc_peel']
sourceTypes     = ['spectrum','line']

outputLabels    = ['hs','tb']
methodLabels    = ['an','fo','mf','mp']
sourceLabels    = ['sp','ln']

spcFile         = "../../data/ti_tube_25keV.dat"

nSamples        = [100, 1000,    25000, 1500]
nSrcSamples     = [100, 10000, 25000, 1500]

nOrders         = 1
srcLineEnergy   = 9e3

nThetaIn        = 5
thetaIn         = str(np.linspace(0,85,nThetaIn)).strip('[]')

nEmAngles       = 10
hsThetaRes      = nEmAngles
thetaEm         = str(np.linspace(0,85,nEmAngles)).strip('[]')
phiEm           = str(np.zeros(nEmAngles)).strip('[]')

matName         = "Test Si/Al"
elements        = 'Si Al'
nElements       = len(elements.split())
elemFracs       = '0.5 0.5'

obsHalfAngle    = 2.5
energyMax       = 2.4e4

parFiles        = []
resFiles        = []


for outi, outType in enumerate(outputTypes):
    for srci, srcType in enumerate(sourceTypes):
        for simi, simMethod in enumerate(methods):

            nameBase = outputLabels[outi]+'_'+sourceLabels[srci]+'_'+methodLabels[simi]

            parFilename = nameBase+'.in'
            resFilename = nameBase+'.nc'

            parFiles.append(parFilename)
            resFiles.append(resFilename)

            parFile = open(parFilename, 'w')

            parFile.write(parTemplate.substitute(result_file         = resFilename,
                                                 src_spectrum_file   = spcFile,
                                                 source_type         = srcType,
                                                 source_line_energy  = srcLineEnergy,
                                                 n_source_samples    = nSrcSamples[simi],
                                                 n_incidence_angles  = nThetaIn,
                                                 incidence_angles    = thetaIn,
                                                 output_type         = outType,
                                                 hs_theta_resolution = hsThetaRes,
                                                 n_emergence_angles  = nEmAngles,
                                                 theta_em            = thetaEm,
                                                 phi_em              = phiEm,
                                                 material_name       = matName,
                                                 n_elements          = nElements,
                                                 elements            = elements,
                                                 element_fractions   = elemFracs,
                                                 n_samples           = nSamples[simi],
                                                 n_orders            = nOrders,
                                                 simulation_method   = simMethod,
                                                 obs_half_angle      = obsHalfAngle,
                                                 energy_max          = energyMax,
                                                 n_threads           = nThreads,
                                                 verbose             = verbose))

            parFile.close()

if '--run-tests' in sys.argv:
    for parFile in parFiles[0:4]:
        print 
        os.system('export OMP_NUM_THREADS=%i; echo \"Running with\" $OMP_NUM_THREADS \"threads\";%s %s'%(nThreads, xrfBin, parFile))
    for resFile in resFiles[0:4]:
        os.system('python compare.py %s' %resFile)
