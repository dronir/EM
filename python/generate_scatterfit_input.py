#!/usr/bin/env python

import sys
import os

from glob import glob
from os.path import join
from optparse import OptionParser
from string import Template

import templates

opt = OptionParser()
opt.add_option('','--n-threads', dest='n_threads', type='int', default=1)
opt.add_option('','--n-fields', dest='n_fields', type='int', default=1)
opt.add_option('','--n-samples', dest='n_samples', type='int', default=1000)
opt.add_option('','--medium-root', dest='medroot', type='str', default='')
opt.add_option('','--result-dir', dest='resdir', type='str', default='')

op, args = opt.parse_args()

nThreads = op.n_threads
nf = op.n_fields
ns = op.n_samples
medRoot = op.medroot
resDir  = op.resdir

parTemplate = Template(templates.t_dscatter)
subTemplate = Template(templates.t_submit)

sTypes = ['smooth','Gaussian','fBm']
HVals  = [0.4, 0.6, 0.8]
lVals  = [0.5, 1.0, 1.5]

SVals  = [0.03, 0.04, 0.05, 0.06]
rVals  = ['020','030','040','050','055']

appFields = '.true.'

for sType in sTypes:
    if sType == 'smooth':
        for rho in rVals:
            medFilename = join(medRoot,'medium_rho_'+rho+'.nc')
            rhoStr = 'R'+rho 
            subFilename = 'scatter_smooth_' + rhoStr + '.sub'
            parFilename = 'scatter_smooth_' + rhoStr + '.in'
            outFilename = 'scatter_smooth_' + rhoStr + '.scangles'
            outFilename = join(resDir, outFilename)

            f = open(parFilename, 'w')
            f.write(parTemplate.substitute(mFilename=medFilename,
                                         oFilename=outFilename,
                                         nSamples=ns*nf,
                                         nFields=nf,
                                         H=0.8,
                                         std=0.0,
                                         applyFields = appFields,
                                         sType = 'constant',
                                         nThreads = nThreads))
            f.close()
    else: 
        P1 = lVals if sType == 'Gaussian' else HVals
        for P in P1:
            for S in SVals:
                for rho in rVals:

                    medFilename = join(medRoot,'medium_rho_'+rho+'.nc')

                    rhoStr = 'R'+rho 
                    hStr = '__P' + ('%4.2f' %P)#.replace('.','_')#[-1]
                    sStr = '__S' + ('%4.2f' %S)#[-2:]

                    subFilename = 'scatter_' + sType + '__' + rhoStr + hStr + sStr + '.sub'
                    parFilename = 'scatter_' + sType + '__' + rhoStr + hStr + sStr + '.in'
                    outFilename = 'scatter_' + sType + '__' + rhoStr + hStr + sStr + '.scangles'
                    outFilename = join(resDir, outFilename)

                    f = open(parFilename, 'w')
                    f.write(parTemplate.substitute(mFilename=medFilename,
                                                   oFilename=outFilename,
                                                   nSamples=ns,
                                                   nFields=nf,
                                                   H=P,
                                                   std=S,
                                                   applyFields = appFields,
                                                   sType = sType,
                                                   nThreads = nThreads))
                    f.close()


for sType in sTypes:
    infiles = glob('*%s*.in'%sType)
    f = open('distribute_%s.condor'%sType,'w')
    f.write('executable = /home/hannu/soft/xr/bin/fit_2hg\n')
    f.write('getenv = True\n')
    f.write('universe   = vanilla\n\n')
    f.write('Requirements = Arch == "X86_64"\n')
    for infile in infiles:
        f.write('arguments = %s\n' %infile)
        f.write('output = %s.out\n'%infile[:-3])
        f.write('error =  %s.err\n'%infile[:-3])
        f.write('log =    %s.log\n'%infile[:-3])
        f.write('queue\n\n')
    f.close()
